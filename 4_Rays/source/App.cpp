/** \file App.cpp */
#include "App.h"

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[]) {
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
Color3 pixelValue(const Radiance3& L, const float k, const float gamma) {
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



// If ray P + tw hits triangle V[0], V[1], V[2], then the function return true, stores the
// barycentric coordinates in b[] and stores the distance to the intersection in t. Otherwise,
// returns false and the other output parameters are undefined.
bool rayTriangleIntersect(const Point3& P, const Vector3& w, const Point3 V[3], float b[3], float& t) {
	// Precision threshold depends on scene scale. Too small leaves holes at edges, too 
	// large expands triangles.
	const float eps = 1e-6;
	
	// Edge vector
	const Vector3& e_1 = V[1] - V[0];
	const Vector3& e_2 = V[2] - V[0];

	// Face normal 
	const Vector3& n = e_1.cross(e_2).direction();

	// Nothing that has any physical meaning. Just a very optimized
	// way to compute 3 cross products.
	const Vector3& q = w.cross(e_2);
	const float a = e_1.dot(q);

	// Backfacing / nearly parallel or close to the limit of precision/
	if ((n.dot(w) >= 0) || (fabsf(a) <= eps)) {
		return false;
	}

	const Vector3& s = (P - V[0]) / a;
	const Vector3& r = s.cross(e_1);

	// Barycentrics
	b[0] = s.dot(q);
	b[1] = r.dot(w);
	b[2] = 1.0f - b[0] - b[1];

	// Intersected outside triangle?
	if ((b[0] < 0.0f) || (b[1] < 0.0f) || (b[2] < 0.0f)) {
		return false;
	}

	t = e_2.dot(r);
	return (t >= 0.0f);
}


PinholeCamera::PinholeCamera(float z_near, float verticalFieldOfView) :
	m_zNear(z_near), m_verticalFieldOfView(verticalFieldOfView)
{
}


void PinholeCamera::getPrimaryRay(float x, float y, int width, int height, Point3& P, Vector3& w) const {
	// Compute the side of a square at z = -1 on our vertical top-to-bottom field of view; the result 
	// is negative because of our convention to have the image plane on the negative axis
	const float side = -2.0f * tanf(m_verticalFieldOfView / 2.0f);
	
	// Invert the y-axis because we're moving from the 2D=down to the 3D y=up coordinate system
	P = Point3(m_zNear * (x / width - 0.5f) * side * width / height,
			   m_zNear * -(y / height - 0.5f) * side,
			   m_zNear);

	// The incoming direction is simply that from the origin to P
	w = P.direction();
}




App::App(const GApp::Settings& settings) : GApp(settings) {
}

// Called before the application loop begins.  Load data here and
// not in the constructor so that common exceptions will be
// automatically caught.
void App::onInit() {
	GApp::onInit();

	setFrameDuration(1.0f / 240.0f);

	// Call setScene(shared_ptr<Scene>()) or setScene(MyScene::create()) to replace
	// the default scene here.

	showRenderingStats = false;

	loadScene(

#       ifndef G3D_DEBUG
		"G3D Sponza"
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
void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces) {
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


void App::onAI() {
	GApp::onAI();
	// Add non-simulation game logic and AI code here
}


void App::onNetwork() {
	GApp::onNetwork();
	// Poll net messages here
}


void App::onSimulation(RealTime rdt, SimTime sdt, SimTime idt) {
	GApp::onSimulation(rdt, sdt, idt);

	// Example GUI dynamic layout code.  Resize the debugWindow to fill
	// the screen horizontally.
	debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


bool App::onEvent(const GEvent& event) {
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


void App::onUserInput(UserInput* ui) {
	GApp::onUserInput(ui);
	(void)ui;
	// Add key handling here based on the keys currently held or
	// ones that changed in the last frame.
}


void App::onPose(Array<shared_ptr<Surface> >& surface, Array<shared_ptr<Surface2D> >& surface2D) {
	GApp::onPose(surface, surface2D);

	// Append any models to the arrays that you want to later be rendered by onGraphics()
}


void App::onGraphics2D(RenderDevice* rd, Array<shared_ptr<Surface2D> >& posed2D) {
	// Render 2D objects like Widgets.  These do not receive tone mapping or gamma correction.
	Surface2D::sortAndRender(rd, posed2D);
}


void App::onCleanup() {
	// Called after the application loop ends.  Place a majority of cleanup code
	// here instead of in the constructor so that exceptions can be caught.
}

void App::drawDebugShapes() {
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

void App::makeGUI() {
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
	m_raytraceSettings.resolutionList = raytracePane->addDropDownList("Resolution", Array<String>({ "1 x 1", "320 x 200", "640 x 400", "1920 x 1080" }));
	m_raytraceSettings.resolutionList->setSelectedIndex(1);
	raytracePane->addCheckBox("Add fixed primitives", &m_raytraceSettings.addFixedPrimitives);
	raytracePane->addCheckBox("Multithreading", &m_raytraceSettings.multithreading);
	GuiNumberBox<int>* raysSlider = raytracePane->addNumberBox<int>("Indirect rays per pixel", &m_raytraceSettings.indirectRaysPerPixel, "", GuiTheme::LINEAR_SLIDER, 0, 2048);
	raysSlider->setWidth(290.0F);
	raysSlider->setCaptionWidth(140.0F);
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

void App::render() {
	drawMessage("Raytracing current scene. Please wait.");

	PinholeCamera camera(-0.001f, 45.0f);

	const Vector2int32 resolution = Vector2int32::parseResolution(m_raytraceSettings.resolutionList->selectedValue());

	shared_ptr<Image> image = Image::create(resolution.x, resolution.y, ImageFormat::RGB32F());
	render(camera, image);

	show(image);
}

void App::render(const PinholeCamera& camera, shared_ptr<Image>& image) const {
	const int width = image->width();
	const int height = image->height();

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			Point3 P;
			Vector3 w;

			// Find the ray through (x,y) and the center of projection
			camera.getPrimaryRay(float(x) + 0.5f, float(y) + 0.5f, width, height, P, w);

			image->set(x, y, L_i(P, w));
		}
	}
}

Radiance3 App::L_i(const Point3& X, const Vector3& wi) const {
	// Find the first intersection with the scene
	const shared_ptr<UniversalSurfel>& s = findFirstIntersection(X, wi);

	if (notNull(s)) {
		return Radiance3::one();
	} else {
		return Radiance3::zero();
	}

}

const shared_ptr<UniversalSurfel>& App::findFirstIntersection(const Point3& X, const Vector3& wi) const {
	// THIS IS A BUG
	return shared_ptr<UniversalSurfel>(new UniversalSurfel());
}
