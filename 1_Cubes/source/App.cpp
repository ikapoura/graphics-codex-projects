/** \file App.cpp */
#include "App.h"

#include <iostream>
#include <fstream>
#include <string>

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

void generateCubeOff();
void generateStaircaseFile();
void pixelsToCubes(const char* imageName, const char* sceneName);

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
#       // ifndef G3D_DEBUG
		// #endif
		"Heightfield"
	);

	// Make the GUI after the scene is loaded because loading/rendering/simulation initialize
	// some variables that advanced GUIs may wish to reference with pointers.
	makeGUI();
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

	giPane->moveRightOf(deferredBox);
	giPane->moveBy(100, 0);

	rendererPane->pack();

	// More examples of debugging GUI controls:
	// debugPane->addCheckBox("Use explicit checking", &explicitCheck);
	// debugPane->addTextBox("Name", &myName);
	// debugPane->addNumberBox("height", &height, "m", GuiTheme::LINEAR_SLIDER, 1.0f, 2.5f);
	// button = debugPane->addButton("Run Simulator");
	// debugPane->addButton("Generate Heightfield", [this](){ generateHeightfield(); });
	// debugPane->addButton("Generate Heightfield", [this](){ makeHeightfield(imageName, scale, "model/heightfield.off"); });

	GuiPane* meshesPane = debugPane->addPane("Meshes", GuiTheme::ORNATE_PANE_STYLE);

	GuiPane* cylinderPane = meshesPane->addPane("Cylinder", GuiTheme::ORNATE_PANE_STYLE);
	cylinderPane->addNumberBox("Radius", &m_cylinderSettings.radius, "m");
	cylinderPane->addNumberBox("Height", &m_cylinderSettings.height, "m");
	cylinderPane->addNumberBox("Steps", &m_cylinderSettings.steps);
	cylinderPane->addButton("Generate", [this]() { generateCylinderOff(m_cylinderSettings.height, m_cylinderSettings.radius, m_cylinderSettings.steps); });
	cylinderPane->pack();

	GuiPane* heightfieldPane = meshesPane->addPane("Heightfield", GuiTheme::ORNATE_PANE_STYLE);
	GuiPane* heightfieldScalePane = heightfieldPane->addPane("Size", GuiTheme::SIMPLE_PANE_STYLE);
	GuiButton* fileDialogButton = heightfieldScalePane->addButton("Select image",
		[this]() {
			m_heightfieldSettings.validFilename = FileDialog::getFilename(m_heightfieldSettings.imageFilename, "jpg", false);
			if (m_heightfieldSettings.validFilename) {
				String& file = m_heightfieldSettings.imageFilename;
				const size_t nameIndex = file.find_last_of("\\") + 1;
				const size_t nameLength = file.npos - nameIndex;
				m_heightfieldSettings.filenameLabel->setCaption(file.substr(nameIndex, nameLength));
			}
	});
	m_heightfieldSettings.filenameLabel = heightfieldScalePane->addLabel("No image selected.");
	m_heightfieldSettings.filenameLabel->moveRightOf(fileDialogButton, 5.0F);
	heightfieldScalePane->addNumberBox("Scale X", &m_heightfieldSettings.scaleX);
	heightfieldScalePane->addNumberBox("Scale Z", &m_heightfieldSettings.scaleZ);
	heightfieldScalePane->addNumberBox("Height (Y=1)", &m_heightfieldSettings.scaleY);
	heightfieldScalePane->pack();

	GuiPane* heightfieldMeshPane = heightfieldPane->addPane("Subdivisions", GuiTheme::SIMPLE_PANE_STYLE);
	heightfieldMeshPane->addNumberBox("X", &m_heightfieldSettings.subdivisionsX, "", GuiTheme::LINEAR_SLIDER, 10, 500);
	heightfieldMeshPane->addNumberBox("Z", &m_heightfieldSettings.subdivisionsZ, "", GuiTheme::LINEAR_SLIDER, 10, 500);
	heightfieldMeshPane->pack();
	heightfieldMeshPane->moveRightOf(heightfieldScalePane, 5.0f);

	heightfieldPane->addButton("Generate",
		[this]() {
			generateHeightfieldOff(m_heightfieldSettings.scaleX, m_heightfieldSettings.scaleY, m_heightfieldSettings.scaleZ);
	});
	heightfieldPane->pack();
	heightfieldPane->moveRightOf(cylinderPane, 5.0F);

	GuiPane* glassPane = meshesPane->addPane("Glass", GuiTheme::ORNATE_PANE_STYLE);
	glassPane->addNumberBox("Steps", &m_glassSettings.steps);
	glassPane->addNumberBox("Quality", &m_glassSettings.quality, "%", GuiTheme::LINEAR_SLIDER, 20, 100);
	glassPane->addButton("Generate", [this]() { generateGlassOff(m_glassSettings.steps, m_glassSettings.quality / 100.0F); });
	glassPane->pack();
	glassPane->moveRightOf(heightfieldPane, 5.0F);

	meshesPane->pack();

	// Final layout
	rendererPane->moveRightOf(infoPane, 10.0F);

	meshesPane->moveRightOf(rendererPane, 10.0F);

	debugWindow->pack();
	debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
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

void App::generateCylinderOff(float height, float radius, int steps)
{
	if (steps < 3) {
		steps = 3;
	}

	drawMessage("Generating cylinder...");

	// Generate the geometry
	const float halfHeight = height / 2.0F;

	Array<Point3> vertices;
	Array<FaceIndices> faces;
	computeCylindricalShape(vertices, faces, steps, { {0.0F, -halfHeight}, {1.0F, -halfHeight}, {1.0F, halfHeight}, {0.0F, halfHeight} });

	// Output to file
	const std::string dir = FileSystem::currentDirectory().c_str();

	std::ofstream ofstream(dir + "\\model\\cylinder.off", std::ios::out);

	ofstream << "OFF\n";
	// No. vertices No. faces No. edges(or 0)
	ofstream << vertices.size() << " " << faces.size() << " 0" << '\n';

	for (const Point3& v : vertices) {
		ofstream << v.x << " " << v.y << " " << v.z << '\n';
	}

	for (const FaceIndices& t : faces) {
		ofstream << t.num;
		for (int i = 0; i < t.num; i++) {
			ofstream << " " << t.indices[i];
		}
		ofstream << '\n';
	}

	ofstream.close();

	// Cleanup
	ArticulatedModel::clearCache();
	
	// Reload the current scene
	loadScene(scene()->name());
}

void App::computeCylindricalShape(Array<Point3>& vertices, Array<FaceIndices>& faces, const int steps, const std::vector<Point2>& rings) const
{
	debugAssertM(rings.size() > 2, "Need at least three rings (top, bottom and in-between) to generate a cylindrical shape.");

	const float deltaAngleDegrees = 360.0F / steps;
	const float deltaAngleRad = deltaAngleDegrees * (float)G3D::pi() / 180.0F;

	// Create a map for the sin and cos to not recalculate them for multiple rings.
	std::vector<std::pair<float, float>> sin_cos;

	for (int i = 0; i < steps; ++i) {
		const float currentAngleRad = deltaAngleRad * i;
		sin_cos.push_back(std::pair<float, float>(std::cosf(currentAngleRad), std::sinf(-currentAngleRad)));
	}

	// Vertices
	for (int r = 1; r < rings.size() - 1; r++) {
		for (int i = 0; i < steps; i++) {
			vertices.append({
				sin_cos[i].first * rings[r].x,
				rings[r].y,
				sin_cos[i].second * rings[r].x
			});
		}
	}

	// Special case for those two because they don't define a ring
	// -- first vertex
	vertices.append({
		0.0F,
		rings[0].y,
		0.0F,
	});
	// -- last vertex
	vertices.append({
		0.0F,
		rings[rings.size() - 1].y,
		0.0F,
	});

	// Faces
	const int centerTopIndex = steps * 2;
	const int centerBottomIndex = steps * 2 + 1;

	for (int r = 1; r < rings.size() - 2; r++) {
		for (int i = 0; i < steps - 1; i++) {
			const int index = (r - 1) * steps + i;

			const int currFirst = index;
			const int currSecond = currFirst + 1;
			const int nextFirst = currFirst + steps;
			const int nextSecond = nextFirst + 1;

			faces.append({ 4, {currFirst, currSecond, nextSecond, nextFirst } });
		}

		// Last face needs to close the last and first edges.
		const int index = (r - 1) * steps + steps - 1;

		const int currFirst = index;
		const int currSecond = (r - 1) * steps;
		const int nextFirst = currFirst + steps;
		const int nextSecond = currSecond + steps;
		faces.append({ 4, {currFirst, currSecond, nextSecond, nextFirst } });
	}

	// Opening and closing triangle fans
	for (int i = 0; i < steps - 1; i++) {
		const int openIndex = vertices.size() - 2;
		const int openFirst = i; // On the first ring
		const int openSecond = i + 1;  // On the first ring
		faces.append({ 3, openIndex, openSecond, openFirst});

		// Last 
		const int closeIndex = vertices.size() - 1;
		const int closeFirst = (rings.size() - 3) * steps + i; // On the last ring
		const int closeSecond = closeFirst + 1; // On the last ring
		faces.append({ 3, closeIndex, closeFirst, closeSecond});
	}

	const int openIndex = vertices.size() - 2;
	const int openFirst = steps - 1; // On the first ring
	const int openSecond = 0;  // On the first ring
	faces.append({ 3, openIndex, openSecond, openFirst});

	// Last 
	const int closeIndex = vertices.size() - 1;
	const int closeFirst = (rings.size() - 3) * steps + steps - 1; // On the last ring
	const int closeSecond = (rings.size() - 3) * steps; // On the last ring
	faces.append({ 3, closeIndex, closeFirst, closeSecond});
}

void App::generateHeightfieldOff(float scaleX, float scaleY, float scaleZ)
{
	if (m_heightfieldSettings.imageFilename.empty()) {
		drawMessage("Please select an image");
		return;
	}

	// Can't let them go negative
	const float minValue = 0.001F;
	scaleX = fmaxf(scaleX, minValue);
	scaleY = fmaxf(scaleY, minValue);
	scaleZ = fmaxf(scaleZ, minValue);

	drawMessage("Generating heightfield...");

	shared_ptr<const Image> imageInput = Image::fromFile(m_heightfieldSettings.imageFilename);
	const int inputWidth = imageInput->width();
	const int inputHeight = imageInput->height();

	const int widthSteps = m_heightfieldSettings.subdivisionsX;
	const float widthOffset = 1.0F / widthSteps;

	const int heightSteps = m_heightfieldSettings.subdivisionsZ;
	const float heightOffset = 1.0F / heightSteps;

	const int numOfVertices = (widthSteps + 1) * (heightSteps + 1);
	const int numOfFaces = widthSteps * heightSteps * 2; // Triangles
	const int numOfEdges = widthSteps * heightSteps * 3 + widthSteps + heightSteps;

	const std::string dir = FileSystem::currentDirectory().c_str();
	std::ofstream ofstream(dir + "\\model\\heightfield.off", std::ios::out);

	ofstream << "OFF\n";
	// No. vertices No. faces No. edges(or 0)
	ofstream << numOfVertices << " " << numOfFaces << " " << numOfEdges << "\n";

	// Write the vertices
	for (int i = 0; i < widthSteps + 1; i++) {
		for (int j = 0; j < heightSteps + 1; j++) {
			const float mapU = i * widthOffset;
			const float mapV = j * heightOffset;

			Point3 position;
			position.x = (mapU - 0.5F) * m_heightfieldSettings.scaleX;
			position.z = (mapV - 0.5F) * m_heightfieldSettings.scaleZ;

			const float imageU = mapU * (inputWidth - 1);
			const float imageV = mapV * (inputHeight - 1);

			const float fracU = imageU - (int)imageU;
			const float fracV = imageV - (int)imageV;

			const int pixelX = (int)imageU % inputWidth;
			const int pixelY = (int)imageV % inputHeight;

			// The 4 boundary vertices are standing in the middle of the corresponding pixel.
			// The rest are interpolated in between using the fractional part.
			Color3& color(imageInput->get<Color3>(pixelX, pixelY));
			if (fracU > 0.0F) {
				color = color * (1.0F - fracU) + fracU * imageInput->get<Color3>(pixelX + 1, pixelY);
			}
			if (fracV > 0.0F) {
				color = color * (1.0F - fracV) + fracV * imageInput->get<Color3>(pixelX, pixelY + 1);
			}
			const float luminance = (color.r * 0.3F + color.g * 0.59F + color.b * 0.11F);  // "luminance"
			position.y = luminance * m_heightfieldSettings.scaleY;

			ofstream << position.x << "  " << position.y << " " << position.z << "\n";
		}
	}

	// Write the faces (triangles)
	for (int i = 0; i < widthSteps; i++) {
		for (int j = 0; j < heightSteps; j++) {
			const int baseIndex = j + i * (heightSteps + 1);
			const int i0 = baseIndex;
			const int i1 = i0 + 1;
			const int i2 = baseIndex + (heightSteps + 1);
			const int i3 = i2 + 1;

			// Triangle 1
			ofstream << "3 " << i0 << " " << i1 << " " << i3 << "\n";

			// Triangle 2
			ofstream << "3 " << i0 << " " << i3 << " " << i2 << "\n";
		}
	}

	ofstream.close();

	// Cleanup
	ArticulatedModel::clearCache();
	
	// Reload the current scene
	loadScene(scene()->name());
}

void App::generateGlassOff(int steps, float qualityPercent)
{
	drawMessage("Generating glass...");

	shared_ptr<const Image> imageInput = Image::fromFile("glass_contour.png");

	const int imageWidth = imageInput->width();
	const int imageHeight = imageInput->height();

	const int topFirstBlackLine = getTopFirstBlackLine(imageInput);
	const int bottomFirstBlackLine = getBottomFirstBlackLine(imageInput);

	if (topFirstBlackLine != -1 && bottomFirstBlackLine != -1) {
		std::vector<Point2int32> contour;
		contour.push_back(Point2int32(0, bottomFirstBlackLine));

		const Color3 white = Color3::white();
		const Color3 black = Color3::black();

		// The segment between the two points we found which is the base and the stem.
		for (int h = bottomFirstBlackLine; h >= topFirstBlackLine; h--) {
			const int firstWhite = getFirstInLine(imageInput, white, black, h, 0);
			contour.push_back(Point2int32(firstWhite, h));
		}

		// The top part of the glass (with the liquid)
		std::vector<Point2int32> contourTail;

		int height = 0;

		for (int h = topFirstBlackLine - 1; h >= 0; h--) {
			const int firstBlack = getFirstInLine(imageInput, black, white, h, 0);

			if (firstBlack != -1) {
				contourTail.push_back(Point2int32(firstBlack, h));

				const int firstWhiteAfterBlack = getFirstInLine(imageInput, white, black, h, firstBlack);
				const int lastBlack = firstWhiteAfterBlack;
				contour.push_back(Point2int32(lastBlack, h));
			} else {
				// Encountered a line without a black pixel and we terminate.
				height = bottomFirstBlackLine - (h - 1);
				break;
			}
		}

		// Finish the contour by copying the tail into the contour inverted.
		for (auto iter = contourTail.rbegin(); iter != contourTail.rend(); iter++) {
			contour.push_back(*iter);
		}

		contour.push_back(Point2int32(0, topFirstBlackLine));

		// Max width to normalize the points
		int maxWidth = 0;
		for (const Point2int32& p : contour) {
			if (p.x > maxWidth) {
				maxWidth = p.x;
			}
		}

		const float widthHeightRatio = (float)maxWidth / (float)height;

		const float cosLimit = 0.98F;
		const float distanceLimit = 0.4F;

		// Compute the cylinder from the points
		std::vector<Point2> glassRings;
		float accumulatedDiffSq = std::numeric_limits<float>::max();

		glassRings.push_back({ contour[0].x / (float)maxWidth * widthHeightRatio, (bottomFirstBlackLine - contour[0].y) / (float)height }); // Insert the first.

		for (int i = 1; i < contour.size() - 1; i++) {
			const Point2int32& p = contour[i];
			const Point2 ring = { p.x / (float)maxWidth * widthHeightRatio, (bottomFirstBlackLine - p.y) / (float)height };

			const int currentRings = glassRings.size();
			if (currentRings > 3) {
				accumulatedDiffSq += (ring - glassRings[currentRings - 1]).length();
				const bool distanceCriterion = accumulatedDiffSq > (distanceLimit * (1.0F - qualityPercent));

				const Vector2 dirPrevious = (glassRings[currentRings - 1] - glassRings[currentRings - 3]).direction();
				const Vector2 dirNext = (ring - glassRings[currentRings - 2]).direction();
				const float cos = dirNext.dot(dirPrevious);
				const bool angleCriterion = (cos < (cosLimit + qualityPercent * (1.0F - cosLimit)));

				if (distanceCriterion || angleCriterion) {
					glassRings.push_back(ring);

					accumulatedDiffSq = 0.0F;
				}
			} else {
				glassRings.push_back(ring);
			}
		}

		glassRings.push_back({ contour[contour.size()-1].x / (float)maxWidth * widthHeightRatio, (bottomFirstBlackLine - contour[contour.size()-1].y) / (float)height }); // Insert the last.

		Array<Point3> vertices;
		Array<FaceIndices> faces;
		computeCylindricalShape(vertices, faces, steps, glassRings);

		// Output to file
		const std::string dir = FileSystem::currentDirectory().c_str();

		std::ofstream ofstream(dir + "\\model\\glass.off", std::ios::out);

		ofstream << "OFF\n";
		// No. vertices No. faces No. edges(or 0)
		ofstream << vertices.size() << " " << faces.size() << " 0" << '\n';

		for (const Point3& v : vertices) {
			ofstream << v.x << " " << v.y << " " << v.z << '\n';
		}

		for (const FaceIndices& t : faces) {
			ofstream << t.num;
			for (int i = 0; i < t.num; i++) {
				ofstream << " " << t.indices[i];
			}
			ofstream << '\n';
		}

		ofstream.close();
	}

	// Cleanup
	ArticulatedModel::clearCache();
	
	// Reload the current scene
	loadScene(scene()->name());
}

int App::getTopFirstBlackLine(shared_ptr<const Image> image) const
{
	const int imageHeight = image->height();

	// Find the first black line from the top
	int h = 0;

	while (h < imageHeight) {
		const Color3& color(image->get<Color3>(0, h));

		if (color == Color3::black()) {
			return h;
		}

		h++;
	}

	return -1;
}

int App::getBottomFirstBlackLine(shared_ptr<const Image> image) const
{
	const int imageHeight = image->height();

	// Find the first black line from the bottom
	int h = imageHeight - 1;

	while (h >= 0) {
		const Color3& color(image->get<Color3>(0, h));

		if (color == Color3::black()) {
			return h;
		}

		h--;
	}

	return -1;
}

int App::getFirstInLine(shared_ptr<const Image> image, const Color3& color, const Color3& colorPrv, int line, int startingOffset) const
{
	debugAssertM(startingOffset >= 0 && startingOffset < image->width(), "Starting offset is not within image bounds.");

	int currentColorIdx = startingOffset;

	int w = startingOffset + 1;
	while (w < image->width()) {
		const Color3& pixelColor(image->get<Color3>(w, line));

		if (pixelColor == colorPrv) {
			currentColorIdx = w;
		}

		if (pixelColor == color) {
			return w * 0.5F + currentColorIdx * 0.5F;
		}

		w++;
	}

	return -1;
}

void generateCubeOff()
{
	const std::string dir = FileSystem::currentDirectory().c_str();

	std::ofstream cubeOfstream(dir + "\\model\\cube.off", std::ios::out);

	cubeOfstream << "OFF\n";
	cubeOfstream << "8 6 0\n";
	cubeOfstream << "-0.5 -0.5 -0.5\n";
	cubeOfstream << " 0.5 -0.5 -0.5\n";
	cubeOfstream << " 0.5  0.5 -0.5\n";
	cubeOfstream << "-0.5  0.5 -0.5\n";
	cubeOfstream << "-0.5 -0.5  0.5\n";
	cubeOfstream << " 0.5 -0.5  0.5\n";
	cubeOfstream << " 0.5  0.5  0.5\n";
	cubeOfstream << "-0.5  0.5  0.5\n";
	cubeOfstream << "4 0 1 2 3\n";
	cubeOfstream << "4 1 5 6 2\n";
	cubeOfstream << "4 5 4 7 6\n";
	cubeOfstream << "4 4 0 3 7\n";
	cubeOfstream << "4 1 0 4 5\n";
	cubeOfstream << "4 2 6 7 3\n";

	cubeOfstream.close();
}

void generateStaircaseFile()
{
	const float stepX = 0.3F;
	const float stepY = 0.1F;
	const float stepZ = 2.0F;

	const std::string dir = FileSystem::currentDirectory().c_str();

	std::ofstream staircaseOfstream(dir + "\\scene\\staircase.Scene.Any", std::ios::out);

	// Write the start of the file with all the configuration pretyped.
	staircaseOfstream << "// -*- c++ -*-\n";
	staircaseOfstream << "{\n";
	staircaseOfstream << "name = \"Staircase\";\n";
	staircaseOfstream << "models = {\n";
	staircaseOfstream << "stepModel = ArticulatedModel::Specification {\n";
	staircaseOfstream << "filename = \"model/crate/crate.obj\";\n";
	staircaseOfstream << "preprocess = {\n";
	staircaseOfstream << "setMaterial(\n";
	staircaseOfstream << "all(),\n";
	staircaseOfstream << "UniversalMaterial::Specification {\n";
	staircaseOfstream << "lambertian = \"marble_texture.jpg\";\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << ");\n";
	staircaseOfstream << "transformGeometry(all(), Matrix4::scale(" << stepX << ", " << stepY << ", " << stepZ << "));\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << "entities = {\n";
	staircaseOfstream << "skybox = Skybox {\n";
	staircaseOfstream << "texture = \"cubemap/whiteroom/whiteroom-*.png\";\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << "sun = Light {\n";
	staircaseOfstream << "attenuation = (0.0001, 0, 1);\n";
	staircaseOfstream << "biradiance = Biradiance3(0.5, 0.6, 0.5);\n";
	staircaseOfstream << "castsShadows = false;\n";
	staircaseOfstream << "enabled = true; \n";
	staircaseOfstream << "frame = CFrame::fromXYZYPRDegrees(1, 1, 1, -30, 45, 0);\n";
	staircaseOfstream << "producesDirectIllumination = true;\n";
	staircaseOfstream << "producesIndirectIllumination = true;\n";
	staircaseOfstream << "spotHalfAngleDegrees = 180;\n";
	staircaseOfstream << "type = \"DIRECTIONAL\";\n";
	staircaseOfstream << "};\n";
	staircaseOfstream << "camera = Camera {\n";
	staircaseOfstream << "frame = CFrame::fromXYZYPRDegrees(0, 0, 5);\n";
	staircaseOfstream << "};\n";

	// Write the staircase.
	const float yawDeltaDegrees = 17.0F;
	const int steps = 50;

	for (int i = 0; i < steps; ++i) {
		const float currentY = i * stepY;
		const float currentYawDegrees = i * yawDeltaDegrees;

		staircaseOfstream << "step" << i << " = VisibleEntity { \n";
		staircaseOfstream << "model = \"stepModel\";\n";
		staircaseOfstream << "frame = CFrame::fromXYZYPRDegrees(0, " << currentY << ", 0, " << currentYawDegrees << ", 0, 0);\n";
		staircaseOfstream << "};\n";
	}

	// Write the closing brackets.
	staircaseOfstream << "};\n";
	staircaseOfstream << "};\n";

	staircaseOfstream.close();
}

void pixelsToCubes(const char* imageName, const char* sceneName)
{
	shared_ptr<const Image> imageInput = Image::fromFile(imageName);

	const int inputWidth = imageInput->width();
	const int inputHeight = imageInput->height();

	std::string outDir = FileSystem::currentDirectory().c_str();
	outDir += "\\scene\\";
	outDir += sceneName;
	outDir += ".Scene.Any";

	std::ofstream sceneOfstream(outDir, std::ios::out);

	// Write the start of the file with all the configuration pretyped.
	sceneOfstream << "// -*- c++ -*-\n";
	sceneOfstream << "{\n";
	sceneOfstream << "name = \"" << sceneName << "\";\n";
	sceneOfstream << "models = {\n";
	sceneOfstream << "cubeModel = ArticulatedModel::Specification {\n";
	sceneOfstream << "filename = \"model/crate/crate.obj\";\n";
	sceneOfstream << "preprocess = {\n";
	sceneOfstream << "setMaterial(\n";
	sceneOfstream << "all(),\n";
	sceneOfstream << "UniversalMaterial::Specification {\n";
	sceneOfstream << "lambertian = Color3(0.5, 0.5, 0.5);\n";
	sceneOfstream << "};\n";
	sceneOfstream << ");\n";
	sceneOfstream << "};\n";
	sceneOfstream << "};\n";
	sceneOfstream << "};\n";
	sceneOfstream << "entities = {\n";
	sceneOfstream << "skybox = Skybox {\n";
	sceneOfstream << "texture = \"cubemap/whiteroom/whiteroom-*.png\";\n";
	sceneOfstream << "};\n";
	sceneOfstream << "sun = Light {\n";
	sceneOfstream << "attenuation = (0.0001, 0, 1);\n";
	sceneOfstream << "biradiance = Biradiance3(0.5, 0.6, 0.5);\n";
	sceneOfstream << "castsShadows = false;\n";
	sceneOfstream << "enabled = true; \n";
	sceneOfstream << "frame = CFrame::fromXYZYPRDegrees(1, 1, 1, -30, -45, 0);\n";
	sceneOfstream << "producesDirectIllumination = true;\n";
	sceneOfstream << "producesIndirectIllumination = true;\n";
	sceneOfstream << "spotHalfAngleDegrees = 180;\n";
	sceneOfstream << "type = \"DIRECTIONAL\";\n";
	sceneOfstream << "};\n";
	sceneOfstream << "camera = Camera {\n";
	sceneOfstream << "frame = CFrame::fromXYZYPRDegrees(0, 0, " << sqrt(inputWidth + inputHeight) << ");\n"; // Move the camera farther to get a better view.
	sceneOfstream << "};\n";

	// Offset to center the cubes to origin.
	const int widthOffset = -inputWidth / 2;
	const int heightOffset = -inputHeight / 2;

	// The model that is used has a scale of 1.
	const float sizeX = 1.0F;
	const float sizeY = 1.0F;
	const float sizeZ = 1.0F;

	for (int i = 0; i < inputWidth; i++) {
		for (int j = 0; j < inputHeight; j++) {
			const Color3& color(imageInput->get<Color3>(i, j));

			const float currentY = (j + widthOffset) * sizeY;
			const float currentX = (i + heightOffset) * sizeX;
			const float currentZ = (color.r * 0.3F + color.g * 0.59F + color.b * 0.11F) * sizeZ; // "luminance"

			sceneOfstream << "cuby" << i << "_" << j << " = VisibleEntity { \n";
			sceneOfstream << "model = \"cubeModel\";\n";
			sceneOfstream << "frame = CFrame::fromXYZYPRDegrees(" << currentX << ", " << currentY << ", " << currentZ << ", 0, 0, 0);\n";
			sceneOfstream << "};\n";
		}
	}

	// Write the closing brackets.
	sceneOfstream << "};\n";
	sceneOfstream << "};\n";

	sceneOfstream.close();
}

