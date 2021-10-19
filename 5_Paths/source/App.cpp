/** \file App.cpp */
#include "App.h"

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

void PinholeCamera::getPrimaryRay(float x, float y, int width, int height, Point3& P, Vector3& w) const
{
	// Init the properties.
	const FOVDirection fovDir = m_projection.fieldOfViewDirection();
	const float fovAngle = m_projection.fieldOfViewAngle();
	debugAssertM(fovDir != FOVDirection::DIAGONAL, "Diagonal FOVDirection is not supported yet");

	const float zNear = m_projection.nearPlaneZ();

	// Compute the origin of the ray.
	const float side = -2.0f * tanf(fovAngle / 2.0f);

	P = Point3(zNear * (x / width - 0.5f) * side,
			   zNear * -(y / height - 0.5f) * side,
			   zNear);

	if (fovDir == FOVDirection::HORIZONTAL) {
		P.y *= (height / (float)width);
	} else if (fovDir == FOVDirection::VERTICAL) {
		P.x *= (width / (float)height);
	}

	// The incoming direction is simply that from the origin to P.
	w = P.direction();

	// Transform based on the orientation.
	P = (m_frame.toMatrix4() * Vector4(P, 1.0f)).xyz();
	w = m_frame.rotation * w;
}






RayTracer::RayTracer(const Settings& settings, const shared_ptr<Scene>& scene) :
	m_settings(settings),
	m_scene(scene)
{
	m_scene->onPose(m_sceneSurfaces);

	m_sceneTriTree = TriTreeBase::create(false);
	m_sceneTriTree->setContents(m_sceneSurfaces);
}

chrono::milliseconds RayTracer::traceImage(const shared_ptr<Camera>& activeCamera, shared_ptr<Image>& image)
{
	chrono::milliseconds elapsedTime(0.0);

	const PinholeCamera camera(activeCamera->frame(), activeCamera->projection());

	// Measure the time of the rendering and not the scene setup.
	Stopwatch stopwatch("Image trace");
	stopwatch.tick();

	// Main loop over the pixels.
	runConcurrently(Point2int32(0, 0), Point2int32(image->width(), image->height()),
		[this, camera, image](Point2int32 point) -> void {
			Random& treadRandom = Random::threadCommon();

			// Always cast a ray through the center of the pixel.
			Radiance3 meanRadiance = getLightTransportPathRadiance(float(point.x) + 0.5f, float(point.y) + 0.5f, camera, image);

			// The pixel numbering starts from (index, index) and spans 1 in each direction. We constrain the max to avoid
			// sampling the same image point from different pixels, because 2 + 1 == 3 and 3 + 0 == 3.
			const float offsetMin = 0.0f;
			const float offsetMax = 1.0f - std::numeric_limits<float>::min();
			// Process the rest of the transport paths.
			for (int i = 1; i < m_settings.numLightTrasportPaths; ++i) {
				float offsetX = treadRandom.uniform(offsetMin, offsetMax);
				float offsetY = treadRandom.uniform(offsetMin, offsetMax);

				meanRadiance += getLightTransportPathRadiance(float(point.x) + offsetX, float(point.y) + offsetY, camera, image);
			}

			meanRadiance /= float(m_settings.numLightTrasportPaths);

			image->set(point.x, point.y, meanRadiance);
		},
		!m_settings.multithreading);

	stopwatch.tock();

	elapsedTime = stopwatch.elapsedDuration<chrono::milliseconds>();

	return elapsedTime;
}

Radiance3 RayTracer::getLightTransportPathRadiance(float imageX, float imageY, const PinholeCamera& camera, const shared_ptr<Image>& image) const
{
	Point3 P;
	Vector3 w;

	camera.getPrimaryRay(imageX, imageY, image->width(), image->height(), P, w);

	const shared_ptr<UniversalSurfel> mainSurfel = findIntersection(P, w, finf(), IntersectionMode::Nearest);

	return L_i(mainSurfel, w, m_settings.maxScatterEvents);
}

shared_ptr<UniversalSurfel> RayTracer::findIntersection(const Point3& X, const Vector3& wi, const float maxDistance, const IntersectionMode mode) const
{
	shared_ptr<UniversalSurfel> result;
	float t = maxDistance;

	intersectTriangulatedSurfaces(X, wi, mode, result, t);

	return result;
}

void RayTracer::intersectTriangulatedSurfaces(const Point3& X, const Vector3& wi, const IntersectionMode mode, shared_ptr<UniversalSurfel>& result, float& t) const
{
	TriTree::Hit hit;
	if (m_sceneTriTree->intersectRay(Ray::fromOriginAndDirection(X, wi, 0.0f, t), hit)) {
		// Sample the triangle.
		shared_ptr<UniversalSurfel> universalSurfel = std::make_shared<UniversalSurfel>();
		shared_ptr<Surfel> surfel = universalSurfel;
		m_sceneTriTree->sample(hit, surfel);

		// Alpha testing.
		if (universalSurfel->coverage == 1.0f) {
			result = universalSurfel;
			t = hit.distance;

			// Return immediately if we find an intersection.
			if (mode == IntersectionMode::First) {
				return;
			}
		}
	}
}

bool RayTracer::visible(const shared_ptr<UniversalSurfel>& s, const Point3& from) const
{
	const float eps = 1e-4f;

	Vector3 dir = (s->position + s->geometricNormal * eps) - from;
	const float distance = dir.length();
	dir /= distance;

	return isNull(findIntersection(from + eps * dir, dir, distance - eps, IntersectionMode::First));
}

Radiance3 RayTracer::L_i(const shared_ptr<UniversalSurfel>& s, const Vector3& wi, const int depth) const
{
	Radiance3 toReturn;

	if (depth >= 0) {
		if (notNull(s)) {
			toReturn = L_o(s, -wi, depth);
		} else {
			toReturn = randomColorFromDirection(wi);
		}
	} else {
		toReturn = Radiance3::black();
	}

	return toReturn;
}

Radiance3 RayTracer::L_o(const shared_ptr<UniversalSurfel>& s, const Vector3& wo, const int depth) const
{
	return s->emittedRadiance(wo) + L_direct(s, wo) + L_indirect(s, wo, depth);
}

Radiance3 RayTracer::L_direct(const shared_ptr<UniversalSurfel>& s, const Vector3& wo) const
{
	const Point3& surfelPos = s->position;
	const Vector3& surfelNormal = s->shadingNormal;

	Radiance3 direct(0.0f, 0.0f, 0.0f);

	for (const shared_ptr<Light>& light : m_scene->lightingEnvironment().lightArray) {
		if (light->producesDirectIllumination()) {
			const Point3& lightPos = light->position().xyz();

			const bool isSurfelVisible = (light->castsShadows()) ? visible(s, lightPos) : true;

			if (isSurfelVisible) {
				const Biradiance3& biradiance = light->biradiance(surfelPos);

				const Vector3& surfelToLightDir = (lightPos - surfelPos).direction();
				const Color3& f = s->finiteScatteringDensity(surfelToLightDir, wo);

				const float cosFactor = fabs(surfelToLightDir.dot(surfelNormal));

				direct += biradiance * f * cosFactor;
			}
		}
	}

	return direct;
}

Radiance3 RayTracer::L_indirect(const shared_ptr<UniversalSurfel>& s, const Vector3& wo, const int depth) const
{
	const float eps = 1e-5f;
	const Vector3& surfelNormal = s->shadingNormal;

	Radiance3 indirect(0.0f, 0.0f, 0.0f);
	if (notNull(s)) {
		Color3 scatterWeight;
		Vector3 scatterDir;
		s->scatter(PathDirection::EYE_TO_SOURCE, -wo, true, Random::threadCommon(), scatterWeight, scatterDir);

		if (scatterWeight != Color3::zero()) {
			const Point3 P = s->position + eps * s->geometricNormal * sign(s->geometricNormal.dot(scatterDir));

			shared_ptr<UniversalSurfel> scatterSurfel = findIntersection(P, scatterDir, finf(), IntersectionMode::Nearest);

			if (notNull(scatterSurfel)) {
				const Radiance3 scatterRadiance = L_i(scatterSurfel, scatterDir, depth - 1);

				indirect += scatterRadiance * scatterWeight;
			}
		}
	}

	return indirect;
}

Biradiance3 RayTracer::randomColorFromDirection(const Vector3& w) const
{
	// For rays that hit the sky, generate a random hue from the direction of the ray.
	const Vector3 wiNormal = w / 2.0f + Vector3(0.5f, 0.5f, 0.5f); // Move all components to [0, 1]

	const float someRandomHue = (wiNormal.x * 0.3f + wiNormal.y * 0.2f + wiNormal.z * 0.5f) + 0.2f;
	return Radiance3::rainbowColorMap(someRandomHue);
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
	m_rayTraceSettings.resolutionList = raytracePane->addDropDownList("Resolution", Array<String>({ "1 x 1", "20 x 20", "320 x 200", "640 x 400", "1920 x 1080" }));
	m_rayTraceSettings.resolutionList->setSelectedIndex(2);
	raytracePane->addCheckBox("Multithreading", &m_rayTraceSettings.multithreading);
	GuiNumberBox<int>* lightTransportPathsSlider = raytracePane->addNumberBox<int>("Light transport paths per pixel", &m_rayTraceSettings.numLightTrasportPaths, "", GuiTheme::LINEAR_SLIDER, 1, 256);
	lightTransportPathsSlider->setWidth(320.0F);
	lightTransportPathsSlider->setCaptionWidth(180.0F);
	GuiNumberBox<int>* scatterEventsSlider = raytracePane->addNumberBox<int>("Maximum scatter events per path", &m_rayTraceSettings.maxScatterEvents, "", GuiTheme::LINEAR_SLIDER, 0, 512);
	scatterEventsSlider->setWidth(340.0F);
	scatterEventsSlider->setCaptionWidth(200.0F);
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

