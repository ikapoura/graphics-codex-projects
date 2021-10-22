/** \file App.cpp */
#include "App.h"

#include <numeric>

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[])
{
	initGLG3D(G3DSpecification());

	GApp::Settings settings(argc, argv);

	// Change the window and other startup parameters by modifying the
	// settings class.  For example:
	settings.window.caption = argv[0];

	// Set enable to catch more OpenGL errors
	// settings.window.debugContext     = true;

	settings.window.fullScreen = false;
	if (settings.window.fullScreen) {
		settings.window.width = 1920;
		settings.window.height = 1080;
	} else {
		settings.window.height = int(OSWindow::primaryDisplayWindowSize().y * 0.95f);
		// Constrain ultra widescreen aspect ratios
		settings.window.width = min(settings.window.height * 1920 / 1080, int(OSWindow::primaryDisplayWindowSize().x * 0.95f));

		// Make even
		settings.window.width -= settings.window.width & 1;
		settings.window.height -= settings.window.height & 1;
	}
	settings.window.resizable = !settings.window.fullScreen;
	settings.window.framed = !settings.window.fullScreen;
	settings.window.defaultIconFilename = "icon.png";

	// Set to true for a significant performance boost if your app can't
	// render at the display frequency, or if you *want* to render faster
	// than the display.
	settings.window.asynchronous = true;

	// Render slightly larger than the screen so that screen-space refraction, bloom,
	// screen-space AO, and screen-space reflection to look good at screen edges. Set to zero for
	// maximum performance and free memory. Increase the second argument to improve AO without affecting
	// color. The third argument is the resolution. Set to 0.5f to render at half-res and upscale,
	// 2.0f to supersample.
	settings.hdrFramebuffer.setGuardBandsAndSampleRate(64, 0, 1.0f);

	settings.renderer.deferredShading = true;
	settings.renderer.orderIndependentTransparency = true;

	settings.dataDir = FileSystem::currentDirectory();

	settings.screenCapture.outputDirectory = FilePath::concat(FileSystem::currentDirectory(), "../journal");
	if (!FileSystem::exists(settings.screenCapture.outputDirectory)) {
		settings.screenCapture.outputDirectory = "";
	}
	settings.screenCapture.includeAppRevision = false;
	settings.screenCapture.includeG3DRevision = false;
	settings.screenCapture.filenamePrefix = "_";

	return App(settings).run();
}



// Get a color from Radiance, saturation and gamma
Color3 PostProcess::pixelValue(const Radiance3& L, const float k, const float gamma)
{
	Radiance3 outL = L;

	// Adjust for constant sensitivity
	outL *= 1.0f / k;

	// Maximum radiance at any frequency
	float m = fmaxf(fmaxf(outL.r, outL.g), fmaxf(outL.b, 1.0f));

	// Normalize the input
	outL *= 1.0f / m;

	// Restore magnitude, but fade towards white when the maximum value approaches 1.0
	m = clamp((m - 1.0f) * 0.2f, 0.0f, 1.0f);
	outL = outL * (1.0f - m) + Radiance3(m, m, m);

	// Gamma encode 
	return Color3(pow(outL.r, 1.0f / gamma),
				  pow(outL.g, 1.0f / gamma),
				  pow(outL.b, 1.0f / gamma));
}

PinholeCamera::PinholeCamera(const CoordinateFrame& frame, const Projection& projection) :
	m_frame(frame), m_projection(projection)
{
}

Ray PinholeCamera::getPrimaryRay(float x, float y, int imageWidth, int imageHeight) const
{
	// Init the properties.
	const FOVDirection fovDir = m_projection.fieldOfViewDirection();
	const float fovAngle = m_projection.fieldOfViewAngle();
	debugAssertM(fovDir != FOVDirection::DIAGONAL, "Diagonal FOVDirection is not supported yet");

	const float zNear = m_projection.nearPlaneZ();

	// Compute the origin of the ray.
	const float side = -2.0f * tanf(fovAngle / 2.0f);

	Point3 P(zNear * (x / imageWidth - 0.5f) * side,
			 zNear * -(y / imageHeight - 0.5f) * side,
			 zNear);

	if (fovDir == FOVDirection::HORIZONTAL) {
		P.y *= (imageHeight / (float)imageWidth);
	} else if (fovDir == FOVDirection::VERTICAL) {
		P.x *= (imageWidth / (float)imageHeight);
	}

	// The incoming direction is simply that from the origin to P.
	Vector3 w = P.direction();

	// Transform based on the orientation.
	P = (m_frame.toMatrix4() * Vector4(P, 1.0f)).xyz();
	w = m_frame.rotation * w;

	return Ray(P, w);
}






RayTracer::RayTracer(const Settings& settings, const shared_ptr<Scene>& scene) :
	m_settings(settings),
	m_scene(scene)
{
	Array<shared_ptr<Surface>> sceneSurfaces;
	m_scene->onPose(sceneSurfaces);

	m_sceneTriTree = TriTreeBase::create(false);
	m_sceneTriTree->setContents(sceneSurfaces);
	m_lastTreeBuildTime = max(m_scene->lastEditingTime(), m_scene->lastStructuralChangeTime(), m_scene->lastVisibleChangeTime());
}

chrono::milliseconds RayTracer::traceImage(const shared_ptr<Camera>& activeCamera, shared_ptr<Image>& image)
{
	rebuildTreeStructureBasedOnLastChange();

	// Measure the time of the rendering and not the scene setup.
	chrono::milliseconds elapsedTime(0.0);
	Stopwatch stopwatch("Image trace");
	stopwatch.tick();

	// Initialize the image
	runConcurrently(Point2int32(0, 0), Point2int32(image->width(), image->height()),
		[image](Point2int32 point) -> void {
			image->set(point.x, point.y, Color3::black());
		},
		!m_settings.multithreading);


	const PinholeCamera camera(activeCamera->frame(), activeCamera->projection());

	Array<Radiance3> radianceBuffer;
	radianceBuffer.resize(m_settings.numLightTrasportPaths);

	Array<Radiance3> modulationBuffer;
	modulationBuffer.resize(m_settings.numLightTrasportPaths);

	Array<Ray> rayBuffer;
	rayBuffer.resize(m_settings.numLightTrasportPaths);

	Array<Ray> shadowRayBuffer;
	shadowRayBuffer.resize(m_settings.numLightTrasportPaths);

	Array<Biradiance3> biradianceBuffer;
	biradianceBuffer.resize(m_settings.numLightTrasportPaths);

	Array<bool> lightShadowedBuffer;
	lightShadowedBuffer.resize(m_settings.numLightTrasportPaths);

	Array<shared_ptr<Surfel>> surfelBuffer;
	surfelBuffer.resize(m_settings.numLightTrasportPaths);

	const float pixelOffsetMin = 0.0f;
	const float pixelOffsetMax = 1.0f - std::numeric_limits<float>::min();

	for (int x = 0, width = image->width(); x < width; x++) {
		for (int y = 0, height = image->height(); y < height; y++) {
			const int pixelIndex1D = x * height + y;

			int activeLightTransportPaths = m_settings.numLightTrasportPaths;
			// Generate the primary rays
			// Always generate a ray through the center of the pixel
			rayBuffer[0] = camera.getPrimaryRay(float(x) + 0.5f, float(y) + 0.5f, width, height);
			// Generate the rest
			runConcurrently(1, m_settings.numLightTrasportPaths,
				[&](int i) -> void {
					float offsetX = uniformRandom(pixelOffsetMin, pixelOffsetMax);
					float offsetY = uniformRandom(pixelOffsetMin, pixelOffsetMax);

					rayBuffer[i] = camera.getPrimaryRay(float(x) + offsetX, float(y) + offsetY, width, height);
				},
				!m_settings.multithreading
			);

			// Initialization
			radianceBuffer.setAll(Radiance3());
			modulationBuffer.setAll(Radiance3(1.0f / float(m_settings.numLightTrasportPaths)));

			// We know that the primary rays are close together.
			TriTree::IntersectRayOptions rayOptions = TriTree::COHERENT_RAY_HINT;

			for (int scatterIdx = 0; scatterIdx < m_settings.maxScatterEvents; scatterIdx++) {
				m_sceneTriTree->intersectRays(rayBuffer, surfelBuffer, rayOptions);

				addEmittedRadiance(rayBuffer, surfelBuffer, modulationBuffer, radianceBuffer);

				if (m_scene->lightingEnvironment().lightArray.size() > 0) {
					sampleDirectLights(surfelBuffer, shadowRayBuffer, biradianceBuffer);

					// If a lightShadowedBuffer value is true, that means that the corresponding surfel is hidden from the light because
					// there is an intersection in between.
					m_sceneTriTree->intersectRays(shadowRayBuffer, lightShadowedBuffer, TriTree::COHERENT_RAY_HINT | TriTree::DO_NOT_CULL_BACKFACES | TriTree::OCCLUSION_TEST_ONLY);

					addDirectIllumination(rayBuffer, surfelBuffer, biradianceBuffer, shadowRayBuffer, lightShadowedBuffer, modulationBuffer, radianceBuffer);
				}

				if (scatterIdx < m_settings.maxScatterEvents - 1) {
					scatterRays(surfelBuffer, rayBuffer, modulationBuffer);
				}

				// Reset the ray options after the first event because the rays won't be coherent anymore.
				rayOptions = 0;
			}

			// No need for a division because we already have accounted for that in the initialization of the
			// modulation buffer.
			Radiance3 sum = std::accumulate(radianceBuffer.begin(), radianceBuffer.end(), Color3::black());
			image->set(x, y, sum);
		}
	}

	stopwatch.tock();

	elapsedTime = stopwatch.elapsedDuration<chrono::milliseconds>();

	return elapsedTime;
}

void RayTracer::addEmittedRadiance(const Array<Ray>& rayBuffer, const Array<shared_ptr<Surfel>>& surfelBuffer, const Array<Radiance3>& modulationBuffer, Array<Radiance3>& radianceBuffer) const
{
	runConcurrently(0, surfelBuffer.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = surfelBuffer[i];
	
			const Vector3& wi = rayBuffer[i].direction();
	
			if (notNull(surfel)) {
				Random& random = Random::threadCommon();
	
				radianceBuffer[i] += surfel->emittedRadiance(-wi) * modulationBuffer[i];
			} else {
				// Rays that didn't intersect anything.
				radianceBuffer[i] += Radiance3(m_settings.environmentBrightness) * modulationBuffer[i];
				// radianceBuffer[i] += randomColorFromDirection(wi) * modulationBuffer[i];
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::sampleDirectLights(const Array<shared_ptr<Surfel>>& surfelBuffer, Array<Ray>& shadowRayBuffer, Array<Biradiance3>& biradianceBuffer) const
{
	const float eps = 1e-4f;

	runConcurrently(0, surfelBuffer.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = surfelBuffer[i];
	
			if (notNull(surfel)) {
				Random& random = Random::threadCommon();
	
				int lightIdx = random.integer(0, m_scene->lightingEnvironment().lightArray.size() - 1);
				const shared_ptr<Light>& light = m_scene->lightingEnvironment().lightArray[lightIdx];

				if (light->producesDirectIllumination()) {
					const Point3& surfelPos = surfel->position;
					const Point3& lightPos = light->position().xyz();

					Vector3 lightToSurfelDir = (surfelPos + surfel->geometricNormal * eps) - lightPos;
					const float lightToSurfelDist = lightToSurfelDir.length();
					lightToSurfelDir /= lightToSurfelDist;

					Point3 shadowRayPos = lightPos + eps * lightToSurfelDir;
					shadowRayBuffer[i] = Ray(shadowRayPos, lightToSurfelDir, 0.0f, lightToSurfelDist - eps);

					biradianceBuffer[i] = light->biradiance(surfelPos);
				}
			} else {
				shadowRayBuffer[i] = Ray(Point3(), Vector3(), 0.0f, 0.0f);
				biradianceBuffer[i] = Biradiance3();
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::addDirectIllumination(const Array<Ray>& rayBuffer, const Array<shared_ptr<Surfel>>& surfelBuffer, const Array<Biradiance3>& biradianceBuffer,
	const Array<Ray>& shadowRayBuffer, const Array<bool>& lightShadowedBuffer, const Array<Radiance3>& modulationBuffer, Array<Radiance3>& radianceBuffer) const
{
	runConcurrently(0, surfelBuffer.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = surfelBuffer[i];
	
			const bool isVisibleFromLight = !lightShadowedBuffer[i];

			if (notNull(surfel) && isVisibleFromLight) {
				const Vector3 surfelToLightDir = -shadowRayBuffer[i].direction();
				const Vector3& surfelNormal = surfel->shadingNormal;

				const Color3 f = surfel->finiteScatteringDensity(surfelToLightDir, -rayBuffer[i].direction());
				const Color3 cosFactor = Color3(fabs(surfelToLightDir.dot(surfelNormal)));
				radianceBuffer[i] += biradianceBuffer[i] * f * cosFactor * modulationBuffer[i];
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::scatterRays(const Array<shared_ptr<Surfel>>& surfelBuffer, Array<Ray>& rayBuffer, Array<Radiance3>& modulationBuffer) const
{
	const float eps = 1e-4f;

	runConcurrently(0, surfelBuffer.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = surfelBuffer[i];

			if (notNull(surfel)) {
				Random& random = Random::threadCommon();

				Color3 scatterWeight;
				Vector3 scatterDir;
				surfel->scatter(PathDirection::EYE_TO_SOURCE, rayBuffer[i].direction(), false, random, scatterWeight, scatterDir);

				const Point3 P = surfel->position + eps * surfel->geometricNormal * sign(surfel->geometricNormal.dot(scatterDir));

				rayBuffer[i] = Ray(P, scatterDir);
				modulationBuffer[i] *= scatterWeight;
			} else {
				modulationBuffer[i] = Radiance3(0.0f);
			}
		},
		!m_settings.multithreading
	);
}

Biradiance3 RayTracer::randomColorFromDirection(const Vector3& w) const
{
	// For rays that hit the sky, generate a random hue from the direction of the ray.
	const Vector3 wiNormal = w / 2.0f + Vector3(0.5f, 0.5f, 0.5f); // Move all components to [0, 1]

	const float someRandomHue = (wiNormal.x * 0.3f + wiNormal.y * 0.2f + wiNormal.z * 0.5f) + 0.2f;
	return Radiance3::rainbowColorMap(someRandomHue);
}

void RayTracer::rebuildTreeStructureBasedOnLastChange()
{
	RealTime lastSceneChangeTime = max(m_scene->lastEditingTime(), m_scene->lastStructuralChangeTime(), m_scene->lastVisibleChangeTime());
	if (lastSceneChangeTime > m_lastTreeBuildTime) {
		m_lastTreeBuildTime = lastSceneChangeTime;

		Array<shared_ptr<Surface>> sceneSurfaces;
		m_scene->onPose(sceneSurfaces);
		m_sceneTriTree->clear();
		m_sceneTriTree->setContents(sceneSurfaces);
	}
}

App::App(const GApp::Settings& settings) :
	GApp(settings)
{
}

// Called before the application loop begins.  Load data here and
// not in the constructor so that common exceptions will be
// automatically caught.
void App::onInit()
{
	GApp::onInit();

	setFrameDuration(1.0f / 240.0f);

	// Call setScene(shared_ptr<Scene>()) or setScene(MyScene::create()) to replace
	// the default scene here.

	showRenderingStats = false;

	loadScene(

#       ifndef G3D_DEBUG
		// "G3D Debug Teapot"
		"G3D Simple Cornell Box (Area Light)" // Load something simple
#       else
		"G3D Simple Cornell Box (Area Light)" // Load something simple
#       endif
	);

	// Make the GUI after the scene is loaded because loading/rendering/simulation initialize
	// some variables that advanced GUIs may wish to reference with pointers.
	makeGUI();
}




// This default implementation is a direct copy of GApp::onGraphics3D to make it easy
// for you to modify. If you aren't changing the hardware rendering strategy, you can
// delete this override entirely.
void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces)
{
	if (!scene()) {
		if ((submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) && (!rd->swapBuffersAutomatically())) {
			swapBuffers();
		}
		rd->clear();
		rd->pushState(); {
			rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());
			drawDebugShapes();
		} rd->popState();
		return;
	}

	BEGIN_PROFILER_EVENT("App::onGraphics3D");
	GBuffer::Specification gbufferSpec = m_gbufferSpecification;

	extendGBufferSpecification(gbufferSpec);
	m_gbuffer->setSpecification(gbufferSpec);
	m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
	m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.hdrFramebuffer.depthGuardBandThickness, m_settings.hdrFramebuffer.colorGuardBandThickness);
	debugAssertGLOk();

	m_renderer->render(rd,
		activeCamera(),
		m_framebuffer,
		scene()->lightingEnvironment().ambientOcclusionSettings.enabled ? m_depthPeelFramebuffer : nullptr,
		scene()->lightingEnvironment(), m_gbuffer,
		allSurfaces,
		[&]() -> decltype(auto) { return scene()->tritree(); }); // decltype(auto) for correct return type deduction in the lambda.

	// Debug visualizations and post-process effects
	rd->pushState(m_framebuffer); {
		// Call to make the App show the output of debugDraw(...)
		rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());
		drawDebugShapes();
		const shared_ptr<Entity>& selectedEntity = (notNull(developerWindow) && notNull(developerWindow->sceneEditorWindow)) ? developerWindow->sceneEditorWindow->selectedEntity() : nullptr;
		scene()->visualize(rd, selectedEntity, allSurfaces, sceneVisualizationSettings(), activeCamera());

		onPostProcessHDR3DEffects(rd);
	} rd->popState();

	// We're about to render to the actual back buffer, so swap the buffers now.
	// This call also allows the screenshot and video recording to capture the
	// previous frame just before it is displayed.
	if (submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) {
		swapBuffers();
	}

	// Clear the entire screen (needed even though we'll render over it, since
	// AFR uses clear() to detect that the buffer is not re-used.)
	rd->clear();

	// Perform gamma correction, bloom, and AA, and write to the native window frame buffer
	m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0),
		settings().hdrFramebuffer.trimBandThickness().x,
		settings().hdrFramebuffer.depthGuardBandThickness.x,
		Texture::opaqueBlackIfNull(notNull(m_gbuffer) ? m_gbuffer->texture(GBuffer::Field::SS_POSITION_CHANGE) : nullptr),
		activeCamera()->jitterMotion());
	END_PROFILER_EVENT();
}


void App::onAI()
{
	GApp::onAI();
	// Add non-simulation game logic and AI code here
}


void App::onNetwork()
{
	GApp::onNetwork();
	// Poll net messages here
}


void App::onSimulation(RealTime rdt, SimTime sdt, SimTime idt)
{
	GApp::onSimulation(rdt, sdt, idt);

	// Example GUI dynamic layout code.  Resize the debugWindow to fill
	// the screen horizontally.
	debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


bool App::onEvent(const GEvent& event)
{
	// Handle super-class events
	if (GApp::onEvent(event)) { return true; }

	// If you need to track individual UI events, manage them here.
	// Return true if you want to prevent other parts of the system
	// from observing this specific event.
	//
	// For example,
	// if ((event.type == GEventType::GUI_ACTION) && (event.gui.control == m_button)) { ... return true; }
	// if ((event.type == GEventType::KEY_DOWN) && (event.key.keysym.sym == GKey::TAB)) { ... return true; }
	// if ((event.type == GEventType::KEY_DOWN) && (event.key.keysym.sym == 'p')) { ... return true; }

	if ((event.type == GEventType::KEY_DOWN) && (event.key.keysym.sym == 'p')) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		r->setDeferredShading(!r->deferredShading());
		return true;
	}

	return false;
}


void App::onUserInput(UserInput* ui)
{
	GApp::onUserInput(ui);
	(void)ui;
	// Add key handling here based on the keys currently held or
	// ones that changed in the last frame.
}


void App::onPose(Array<shared_ptr<Surface> >& surface, Array<shared_ptr<Surface2D> >& surface2D)
{
	GApp::onPose(surface, surface2D);

	// Append any models to the arrays that you want to later be rendered by onGraphics()
}


void App::onGraphics2D(RenderDevice* rd, Array<shared_ptr<Surface2D> >& posed2D)
{
	// Render 2D objects like Widgets.  These do not receive tone mapping or gamma correction.
	Surface2D::sortAndRender(rd, posed2D);
}


void App::onCleanup()
{
	// Called after the application loop ends.  Place a majority of cleanup code
	// here instead of in the constructor so that exceptions can be caught.
}

void App::drawDebugShapes()
{
	GApp::drawDebugShapes();

	const shared_ptr<DefaultRenderer>& defaultRenderer = dynamic_pointer_cast<DefaultRenderer>(m_renderer);

	if (isNull(defaultRenderer)) {
		return;
	}

	const Array<shared_ptr<DDGIVolume>>& volumeArray = scene()->lightingEnvironment().ddgiVolumeArray;

	// Draw debug shapes for the DDGIVolume.
	for (int i = 0; i < volumeArray.size(); ++i) {
		if (defaultRenderer->m_showProbeLocations[i]) {
			// TODO: enable probe visualization radius.
			float radius = 0.1f;
			// TODO: enable depth visualization from second to last parameter.
			defaultRenderer->m_ddgiVolumeArray[i]->debugRenderProbeVisualization(renderDevice, activeCamera(), false, radius);
		}
	}
}

void App::makeGUI()
{
	debugWindow->setVisible(true);
	developerWindow->videoRecordDialog->setEnabled(true);
	GuiPane* infoPane = debugPane->addPane("Info", GuiTheme::ORNATE_PANE_STYLE);

	// Example of how to add debugging controls
	infoPane->addLabel("You can add GUI controls");
	infoPane->addLabel("in App::onInit().");
	infoPane->addButton("Exit", [this]() { m_endProgram = true; });
	infoPane->pack();

	GuiPane* rendererPane = debugPane->addPane("DefaultRenderer", GuiTheme::ORNATE_PANE_STYLE);

	// showInTextureBrowser("G3D::GBuffer/CS_NORMAL");

	GuiCheckBox* deferredBox = rendererPane->addCheckBox("Deferred Shading",
		Pointer<bool>([&]() {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		return r && r->deferredShading();
	},
			[&](bool b) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (r) { r->setDeferredShading(b); }
	}));
	rendererPane->addCheckBox("Order-Independent Transparency",
		Pointer<bool>([&]() {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		return r && r->orderIndependentTransparency();
	},
			[&](bool b) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (r) { r->setOrderIndependentTransparency(b); }
	}));

	GuiPane* giPane = rendererPane->addPane("Ray Tracing", GuiTheme::SIMPLE_PANE_STYLE);
	giPane->addCheckBox("Diffuse",
		Pointer<bool>([&]() {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		return r && r->enableDiffuseGI();
	},
			[&](bool b) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (r) { r->setEnableDiffuseGI(b); }
	}));
	giPane->addCheckBox("Glossy",
		Pointer<bool>([&]() {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		return r && r->enableGlossyGI();
	},
			[&](bool b) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (r) { r->setEnableGlossyGI(b); }
	}));
	giPane->addCheckBox("Show Probes",
		Pointer<bool>([&]() {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (notNull(r)) {
			bool allEnabled = r->m_ddgiVolumeArray.size() > 0;
			for (int i = 0; i < r->m_ddgiVolumeArray.size(); ++i) {
				allEnabled = allEnabled && r->m_showProbeLocations[i];
			}
			return allEnabled;
		}
		return false;
	},
			[&](bool b) {
		const shared_ptr<DefaultRenderer>& r = dynamic_pointer_cast<DefaultRenderer>(m_renderer);
		if (notNull(r)) {
			for (int i = 0; i < r->m_ddgiVolumeArray.size(); ++i) {
				r->m_showProbeLocations[i] = b;
			}
		}
	}), GuiTheme::TOOL_CHECK_BOX_STYLE);
	giPane->pack();
	giPane->moveRightOf(deferredBox, 100);

	rendererPane->pack();
	rendererPane->moveRightOf(infoPane, 10);

	GuiPane* raytracePane = debugPane->addPane("Offline Ray Trace", GuiTheme::ORNATE_PANE_STYLE);
	m_rayTraceSettings.resolutionList = raytracePane->addDropDownList("Resolution", Array<String>({ "1 x 1", "20 x 20", "320 x 200", "640 x 400", "1280 x 720", "1920 x 1080" }));
#ifndef G3D_DEBUG
	m_rayTraceSettings.resolutionList->setSelectedIndex(3);
#else
	m_rayTraceSettings.resolutionList->setSelectedIndex(1);
#endif // !G3D_DEBUG
	raytracePane->addCheckBox("Multithreading", &m_rayTraceSettings.multithreading);
	GuiNumberBox<int>* lightTransportPathsSlider = raytracePane->addNumberBox<int>("Light transport paths per pixel", &m_rayTraceSettings.numLightTrasportPaths, "", GuiTheme::LINEAR_SLIDER, 1, 4096);
	lightTransportPathsSlider->setWidth(320.0F);
	lightTransportPathsSlider->setCaptionWidth(180.0F);
	GuiNumberBox<int>* scatterEventsSlider = raytracePane->addNumberBox<int>("Maximum scatter events per path", &m_rayTraceSettings.maxScatterEvents, "", GuiTheme::LINEAR_SLIDER, 1, 36);
	scatterEventsSlider->setWidth(340.0F);
	scatterEventsSlider->setCaptionWidth(200.0F);
	GuiNumberBox<float>* envBrightnessSlider = raytracePane->addNumberBox<float>("Environment brightness", &m_rayTraceSettings.environmentBrightness, "", GuiTheme::LINEAR_SLIDER, 0.0f, 1.0f);
	envBrightnessSlider->setWidth(340.0F);
	envBrightnessSlider->setCaptionWidth(200.0F);
	raytracePane->addButton("Render", this, &App::render);
	raytracePane->pack();
	raytracePane->moveRightOf(rendererPane, 10);

	// More examples of debugging GUI controls:
	// debugPane->addCheckBox("Use explicit checking", &explicitCheck);
	// debugPane->addTextBox("Name", &myName);
	// debugPane->addNumberBox("height", &height, "m", GuiTheme::LINEAR_SLIDER, 1.0f, 2.5f);
	// button = debugPane->addButton("Run Simulator");
	// debugPane->addButton("Generate Heightfield", [this](){ generateHeightfield(); });
	// debugPane->addButton("Generate Heightfield", [this](){ makeHeightfield(imageName, scale, "model/heightfield.off"); });

	debugWindow->pack();
	debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}

Vector2int32 App::resolution() const
{
	return Vector2int32::parseResolution(m_rayTraceSettings.resolutionList->selectedValue());
}

String App::durationToString(chrono::milliseconds duration) const
{
	const int h = int(std::chrono::duration_cast<chrono::hours>(duration).count());
	duration -= std::chrono::duration_cast<chrono::milliseconds>(chrono::hours(h));

	const int m = int(std::chrono::duration_cast<chrono::minutes>(duration).count());
	duration -= std::chrono::duration_cast<chrono::milliseconds>(chrono::minutes(m));

	const int s = int(std::chrono::duration_cast<chrono::seconds>(duration).count());
	duration -= std::chrono::duration_cast<chrono::milliseconds>(chrono::seconds(s));

	const int ms = int(duration.count());

	const size_t BUFFER_SIZE = 128;
	char buffer[BUFFER_SIZE];
	snprintf(buffer, BUFFER_SIZE, "%dh %dm %ds %dms", h, m, s, ms);

	return String(buffer);
}

void App::render()
{
	drawMessage("Raytracing current scene. Please wait.");

	const Vector2int32 res = resolution();
	shared_ptr<Image> image = Image::create(res.x, res.y, ImageFormat::RGB32F());

	RayTracer rayTracer(m_rayTraceSettings, scene());

	activeCamera()->filmSettings().setBloomStrength(0.0f);

	const chrono::milliseconds durationMs = rayTracer.traceImage(activeCamera(), image);

	// Convert to texture to post process.
	shared_ptr<Texture> src = Texture::fromImage("Render result", image);
	if (m_result) {
		m_result->resize(res.x, res.y);
	}

	m_film->exposeAndRender(renderDevice, activeCamera()->filmSettings(), src,
		settings().hdrFramebuffer.trimBandThickness().x,
		settings().hdrFramebuffer.depthGuardBandThickness.x, m_result);

	const String durationPrintOutput = String("Render duration: ") + durationToString(durationMs);

	consolePrint(durationPrintOutput);

	show(m_result, durationPrintOutput);
}

