/** \file App.cpp */
#include "App.h"

#include <cmath>
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


int PathBuffers::size() const
{
	return imageCoordinates.size();
}

bool PathBuffers::empty() const
{
	return size() == 0;
}

void PathBuffers::resize(size_t size)
{
	imageCoordinates.resize(size);
	modulation.resize(size);
	rays.resize(size);
	surfels.resize(size);
	surfelShadowed.resize(size);
}

void PathBuffers::resizeShadowRays(size_t size)
{
	shadowRays.resize(size);
	shadowRayOriginalIndex.resize(size);
	shadowRayResult.resize(size);
	biradiance.resize(size);
}

void PathBuffers::clearRedundantSurfels(float minModulationSum)
{
	for (int i = 0; i < surfels.size(); i++) {
		if (isNull(surfels[i]) || (modulation[i].sum() < minModulationSum)) {
			fastRemove(i);
			i--;
		}
	}
}

void PathBuffers::clearRedundantShadowRays()
{
	for (int i = 0; i < shadowRays.size(); i++) {
		if (shadowRayOriginalIndex[i] == -1) {
			fastRemoveShadows(i);
			i--;
		}
	}
}

void PathBuffers::postProcessShadowResults()
{
	for (int i = 0; i < shadowRays.size(); i++) {
		const int originalIdx = shadowRayOriginalIndex[i];
		surfelShadowed[originalIdx] = shadowRayResult[i];
	}
}

void PathBuffers::fastRemove(int index)
{
	imageCoordinates.fastRemove(index);
	modulation.fastRemove(index);
	rays.fastRemove(index);
	surfels.fastRemove(index);
	surfelShadowed.fastRemove(index);
}

void PathBuffers::fastRemoveShadows(int index)
{
	shadowRays.fastRemove(index);
	shadowRayOriginalIndex.fastRemove(index);
	shadowRayResult.fastRemove(index);
	biradiance.fastRemove(index);
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

chrono::milliseconds RayTracer::traceImage(const shared_ptr<Camera>& activeCamera, shared_ptr<Image>& radianceImage)
{
	rebuildTreeStructureBasedOnLastChange();

	// Measure the time of the rendering and not the scene setup.
	chrono::milliseconds elapsedTime(0.0);
	Stopwatch stopwatch("Image trace");
	stopwatch.tick();

	// Initialize the image
	radianceImage->setAll(Color3::black());

	// Cache the direct illumination lights.
	Array<shared_ptr<Light>> directLightsBuffer;
	for (shared_ptr<Light> light : m_scene->lightingEnvironment().lightArray) {
		if (light->enabled() && light->visible() && light->type() != Light::Type::AREA) {
			directLightsBuffer.push_back(light);
		}
	}

	const PinholeCamera camera(activeCamera->frame(), activeCamera->projection());

	const Point2int32 imageSize(radianceImage->width(), radianceImage->height());

	shared_ptr<Image> pathWeightImage = Image::create(imageSize.x, imageSize.y, ImageFormat::R32F());

	PathBuffers buffers;

	// We fully trace one path for all pixels and then move to the next.
	for (int pathIdx = 0; pathIdx < m_settings.numLightTrasportPaths; pathIdx++) {
		buffers.resize(imageSize.x * imageSize.y);

		initializeTransportPaths(camera, imageSize, pathIdx, buffers, pathWeightImage);

		// The first rays we intersect are the primary rays and we know that they are close together.
		TriTree::IntersectRayOptions rayOptions = TriTree::COHERENT_RAY_HINT;

		for (int scatterIdx = 0; scatterIdx < m_settings.maxScatterEvents && !buffers.empty(); scatterIdx++) {
			m_sceneTriTree->intersectRays(buffers.rays, buffers.surfels, rayOptions);

			addEmittedRadiance(buffers, radianceImage);

			// We cull after the emitted radiance to have some default contribution from the rays that didn't hit anything.
			buffers.clearRedundantSurfels(0.03f);

			if (directLightsBuffer.size() > 0) {
				buffers.resizeShadowRays(buffers.surfels.size());

				sampleDirectLights(buffers, directLightsBuffer, imageSize, pathIdx, scatterIdx);

				buffers.clearRedundantShadowRays();

				// If a lightShadowedBuffer value is true, that means that the corresponding surfel is hidden from the light because
				// there is an intersection in between.
				rayOptions = TriTree::COHERENT_RAY_HINT | TriTree::DO_NOT_CULL_BACKFACES | TriTree::OCCLUSION_TEST_ONLY;
				m_sceneTriTree->intersectRays(buffers.shadowRays, buffers.shadowRayResult, rayOptions);

				buffers.postProcessShadowResults();

				addDirectIllumination(buffers, radianceImage);
			}

			if (scatterIdx < m_settings.maxScatterEvents - 1) {
			 	scatterRays(buffers);
			}

			// Reset the ray options after the first event because the rays won't be coherent anymore.
			rayOptions = 0;
		}

		debugPrintf("Raytracing scene. Current transport path: %d of %d total.\n", pathIdx+1, m_settings.numLightTrasportPaths);
	}

	// The final weighting of the radiance.
	runConcurrently(Point2int32(0, 0), imageSize, [&](Point2int32 pixel) {
		radianceImage->set(pixel, radianceImage->get<Color3>(pixel) / max(0.0001f, pathWeightImage->get<Color1>(pixel).value));
	});

	stopwatch.tock();

	elapsedTime = stopwatch.elapsedDuration<chrono::milliseconds>();

	return elapsedTime;
}

void RayTracer::initializeTransportPaths(const PinholeCamera& camera, const Point2int32& imageSize, int pathIdx,
	PathBuffers& buffers, shared_ptr<Image>& pathWeightImage) const
{
	const bool isFirstPath = (pathIdx == 0);

	constexpr float pixelOffsetMin = 0.0f;
	constexpr float pixelOffsetMax = 1.0f - std::numeric_limits<float>::min();

	runConcurrently(0, buffers.rays.size(),
		[&](int i) -> void {
			Random& random = Random::threadCommon();

			const int x = i / imageSize.y;
			const int y = i - x * imageSize.y;

			float offsetX = 0.5f;
			float offsetY = 0.5f;

			// If it is the first path, we send a ray through the center of the pixel.
			if (!isFirstPath) {
				float offsetX = random.uniform(pixelOffsetMin, pixelOffsetMax);
				float offsetY = random.uniform(pixelOffsetMin, pixelOffsetMax);
			}

			Point2 imgCoords(float(x) + offsetX, float(y) + offsetY);

			buffers.rays[i] = camera.getPrimaryRay(imgCoords.x, imgCoords.y, imageSize.x, imageSize.y);

			// Initialization
			buffers.modulation[i] = Radiance3(1.0f);

			buffers.imageCoordinates[i] = imgCoords;
		},
		!m_settings.multithreading
	);

	// We need the total weights that have been applied to all pixels to normalize the final value.
	const Color1 one(1.0f);
	for (const Point2& imgCoords : buffers.imageCoordinates) {
		pathWeightImage->bilinearIncrement(imgCoords, one);
	}
}

void RayTracer::addEmittedRadiance(PathBuffers& buffers, shared_ptr<Image>& image) const
{
	runConcurrently(0, buffers.surfels.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = buffers.surfels[i];
	
			if (notNull(surfel)) {
				Random& random = Random::threadCommon();

				const Vector3& wi = buffers.rays[i].direction();
	
				image->bilinearIncrement(buffers.imageCoordinates[i], surfel->emittedRadiance(-wi) * buffers.modulation[i]);
			} else {
				// Rays that didn't intersect anything.
				image->bilinearIncrement(buffers.imageCoordinates[i], Radiance3(m_settings.environmentBrightness) * buffers.modulation[i]);
				// const Vector3& wi = buffers.rays[i].direction();
				// image->bilinearIncrement(buffers.imageCoordinates[i], randomColorFromDirection(wi) * buffers.modulation[i]);
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::sampleDirectLights(PathBuffers& buffers, const Array<shared_ptr<Light>>& directLightsBuffer, const Point2int32& imageSize, int pathIdx, int scatterIdx) const
{
	runConcurrently(0, buffers.surfels.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = buffers.surfels[i];

			const int serializedPixelId = int(floor(buffers.imageCoordinates[i].x)) * imageSize.x + int(floor(buffers.imageCoordinates[i].y));
			int pixelIdx = serializedPixelId * m_settings.maxScatterEvents + scatterIdx;

			int lightIdx = 0;
			float lightWeight = 1.0f;
			float lightAreaTimesPdfValue = 1.0f;

			Point3 lightPos = lightImportanceSampling(directLightsBuffer, surfel, pixelIdx, pathIdx, Random::threadCommon(), lightIdx, lightWeight, lightAreaTimesPdfValue);

			// Sample the light
			const shared_ptr<Light>& selectedLight = directLightsBuffer[lightIdx];

			Biradiance3 biradiance = selectedLight->biradiance(surfel->position, lightPos) * (lightWeight / lightAreaTimesPdfValue);

			if (biradiance.sum() > 0.0f) {
				buffers.biradiance[i] = biradiance;
				buffers.shadowRays[i] = generateShadowRay(surfel->position, surfel->geometricNormal, lightPos);
				buffers.shadowRayOriginalIndex[i] = i;
				buffers.surfelShadowed[i] = false; // Simple initialization.
			} else {
				buffers.shadowRayOriginalIndex[i] = -1;
				buffers.surfelShadowed[i] = true; // Mark as shadowed to avoid redundant computations.
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::addDirectIllumination(PathBuffers& buffers, shared_ptr<Image>& image) const
{
	runConcurrently(0, buffers.surfels.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = buffers.surfels[i];
	
			const bool isVisibleFromLight = !buffers.surfelShadowed[i];

			if (isVisibleFromLight) {
				const Vector3 surfelToLightDir = -buffers.shadowRays[i].direction();
				const Vector3& surfelNormal = surfel->shadingNormal;

				const Color3 f = surfel->finiteScatteringDensity(surfelToLightDir, -buffers.rays[i].direction());
				const Color3 cosFactor = Color3(fabs(surfelToLightDir.dot(surfelNormal)));

				image->bilinearIncrement(buffers.imageCoordinates[i], buffers.biradiance[i] * f * cosFactor * buffers.modulation[i]);
			}
		},
		!m_settings.multithreading
	);
}

void RayTracer::scatterRays(PathBuffers& buffers) const
{
	const float eps = 1e-4f;

	runConcurrently(0, buffers.surfels.size(),
		[&](int i) -> void {
			const shared_ptr<Surfel>& surfel = buffers.surfels[i];

			Random& random = Random::threadCommon();

			Color3 scatterWeight;
			Vector3 scatterDir;
			surfel->scatter(PathDirection::EYE_TO_SOURCE, -buffers.rays[i].direction(), false, random, scatterWeight, scatterDir);

			const Point3 P = surfel->position + eps * surfel->geometricNormal * sign(surfel->geometricNormal.dot(scatterDir));

			buffers.rays[i] = Ray(P, scatterDir);
			buffers.modulation[i] *= scatterWeight;

		},
		!m_settings.multithreading
	);
}

Point3 RayTracer::lightImportanceSampling(const Array<shared_ptr<Light>>& directLightsBuffer, const shared_ptr<Surfel>& surfel, int pixelIdx, int pathIdx, Random& random,
	int& selectedLightIdx, float& selectedLightWeight, float& selectedLightAreaTimesPdfValue) const
{
	const int numLights = directLightsBuffer.size();

	const Point3& surfelPos = surfel->position;
	// Sum the biradiance terms of each light to get an estimate of how much they contribute
	// to the surfel's outgoing radiance.
	thread_local Array<Radiance3> biradianceTimesAreaTimesPdfValue;
	biradianceTimesAreaTimesPdfValue.resize(numLights);

	thread_local Array<float> areaTimesPdfValue;
	areaTimesPdfValue.resize(numLights);

	thread_local Array<float> biradianceSums;
	biradianceSums.resize(numLights);

	thread_local Array<Point3> newLightPositions;
	newLightPositions.resize(numLights);

	Point3 resultPosition;

	if (numLights > 1) {
		for (int l = 0; l < numLights; l++) {
			newLightPositions[l] = directLightsBuffer[l]->lowDiscrepancySolidAnglePosition(pixelIdx, l, pathIdx, m_settings.numLightTrasportPaths, surfel->position, areaTimesPdfValue[l]);

			biradianceTimesAreaTimesPdfValue[l] = directLightsBuffer[l]->biradiance(surfelPos, newLightPositions[l]) / areaTimesPdfValue[l];

			biradianceSums[l] = biradianceTimesAreaTimesPdfValue[l].sum();
		}

		float accumulatedBiradiance = std::accumulate(biradianceSums.begin(), biradianceSums.end(), 0.0f);

		// Select a light with probability relative to the previous sums.
		float randomBiradiance = random.uniform(0.0f, accumulatedBiradiance);

		for (int l = 0; l < numLights; l++) {
			randomBiradiance -= biradianceSums[l];

			if (randomBiradiance < 0) {
				selectedLightIdx = l;
				// It's the inverse probability of selecting it. In more verbose form:
				// 1.0f / (biradianceSum[i] / accumulatedBiradiance)
				selectedLightWeight = accumulatedBiradiance / max(biradianceSums[l], 0.0001f);
				selectedLightAreaTimesPdfValue = areaTimesPdfValue[l];
				resultPosition = newLightPositions[l];
				break;
			}
		}
	} else {
		selectedLightIdx = 0;
		selectedLightWeight = 1.0f;
		resultPosition = directLightsBuffer[0]->lowDiscrepancySolidAnglePosition(pixelIdx, selectedLightIdx, pathIdx, m_settings.numLightTrasportPaths, surfel->position, selectedLightAreaTimesPdfValue);
	}

	return resultPosition;
}

Ray RayTracer::generateShadowRay(const Point3& surfelPos, const Vector3& surfelGNormal, const Point3& lightPos) const
{
	const float eps = 1e-4f;

	Vector3 lightToSurfelDir = (surfelPos + surfelGNormal * eps) - lightPos;
	const float lightToSurfelDist = lightToSurfelDir.length();
	lightToSurfelDir /= lightToSurfelDist;

	Point3 shadowRayPos = lightPos + eps * lightToSurfelDir;
	return Ray(shadowRayPos, lightToSurfelDir, 0.0f, lightToSurfelDist - eps);
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

	setFrameDuration(1.0f / 60.0f);

	// Call setScene(shared_ptr<Scene>()) or setScene(MyScene::create()) to replace
	// the default scene here.

	showRenderingStats = false;

	loadScene(

#       ifndef G3D_DEBUG
		// "G3D Debug Teapot"
		"G3D Simple Cornell Box (Spheres)" // Load something simple
#       else
		"G3D Simple Cornell Box (Spheres)" // Load something simple
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

	GuiPane* raytracePane = debugPane->addPane("Offline Path Trace", GuiTheme::ORNATE_PANE_STYLE);
	m_rayTraceSettings.resolutionList = raytracePane->addDropDownList("Resolution", Array<String>({ "1 x 1", "20 x 20", "320 x 200", "640 x 400", "1280 x 720", "1920 x 1080" }));
#ifndef G3D_DEBUG
	m_rayTraceSettings.resolutionList->setSelectedIndex(3);
#else
	m_rayTraceSettings.resolutionList->setSelectedIndex(2);
#endif // !G3D_DEBUG
	raytracePane->addCheckBox("Multithreading", &m_rayTraceSettings.multithreading);
	GuiNumberBox<int>* lightTransportPathsSlider = raytracePane->addNumberBox<int>("Light transport paths per pixel", &m_rayTraceSettings.numLightTrasportPaths, "", GuiTheme::LINEAR_SLIDER, 1, 4096);
	lightTransportPathsSlider->setWidth(320.0F);
	lightTransportPathsSlider->setCaptionWidth(180.0F);
	GuiNumberBox<int>* scatterEventsSlider = raytracePane->addNumberBox<int>("Maximum scatter events per path", &m_rayTraceSettings.maxScatterEvents, "", GuiTheme::LINEAR_SLIDER, 1, 36);
	scatterEventsSlider->setWidth(340.0F);
	scatterEventsSlider->setCaptionWidth(200.0F);
	GuiNumberBox<float>* envBrightnessSlider = raytracePane->addNumberBox<float>("Environment brightness", &m_rayTraceSettings.environmentBrightness, "", GuiTheme::LINEAR_SLIDER, 0.0f, 1000.0f);
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

	float prevSensitivity = activeCamera()->filmSettings().sensitivity();
	activeCamera()->filmSettings().setSensitivity(prevSensitivity * 0.5f);

	const chrono::milliseconds durationMs = rayTracer.traceImage(activeCamera(), image);

	// Convert to texture to post process.
	shared_ptr<Texture> src = Texture::fromImage("Render result", image);
	if (m_result) {
		m_result->resize(res.x, res.y);
	}

	m_film->exposeAndRender(renderDevice, activeCamera()->filmSettings(), src,
		settings().hdrFramebuffer.trimBandThickness().x,
		settings().hdrFramebuffer.depthGuardBandThickness.x, m_result);

	activeCamera()->filmSettings().setSensitivity(prevSensitivity);

	const String durationPrintOutput = String("Render duration: ") + durationToString(durationMs);

	consolePrint(durationPrintOutput);

	show(m_result, durationPrintOutput);
}

