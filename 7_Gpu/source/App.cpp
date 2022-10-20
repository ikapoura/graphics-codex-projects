/**
 C++ framework for GPU ray marching.
 Copyright 2014-2016 Morgan McGuire
 mcguire@cs.williams.edu
*/
#include "App.h"

//////////////////////////////////////////////////////////////////////////////

App::App(const GApp::Settings& settings) :
	GApp(settings),
	m_dirCurrent(FileSystem::currentDirectory() + String("\\")),
	m_dirTextures(m_dirCurrent + String("textures\\")),
	m_dirShaders(m_dirCurrent + String("..\\source\\shaders\\"))
{
}


void App::onInit()
{
	GApp::onInit();

	setFrameDuration(1.0f / 60.0f);
	setLowerFrameRateInBackground(false);
	showRenderingStats = false;

	developerWindow->sceneEditorWindow->setVisible(false);
	developerWindow->cameraControlWindow->setVisible(false);
	developerWindow->setVisible(false);
	developerWindow->videoRecordDialog->setEnabled(true);

	Texture::Specification e;
	e.filename = System::findDataFile(m_dirTextures + String("cubemap\\cornellbox\\empty-RG-*.png"));
	e.encoding.format = ImageFormat::SRGB8();
	e.dimension = Texture::DIM_CUBE_MAP;

	m_environmentMap = Texture::create(e);

	for (int i = 0; i < 4; ++i) {
		m_veniceTexture.append(Texture::fromFile(m_dirTextures + format("iChannel%d.jpg", i)));
	}

	debugCamera()->setFieldOfViewAngle(40 * units::degrees());

	// Disable needless post-processing.
	debugCamera()->filmSettings().setAntialiasingEnabled(false);
	debugCamera()->filmSettings().setVignetteSizeFraction(0.0f);
	debugCamera()->filmSettings().setBloomStrength(0.0f);
}


void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces)
{
	// Bind the main framebuffer
	rd->push2D(m_framebuffer); {
		rd->clear();

		Args args;

		// Prepare the arguments for the shader function invoked below
		args.setUniform("cameraToWorldMatrix", activeCamera()->frame());

		m_environmentMap->setShaderArgs(args, "environmentMap.", Sampler::cubeMap());
		args.setUniform("environmentMap_MIPConstant", std::log2(float(m_environmentMap->width() * sqrt(3.0f))));

		args.setUniform("tanHalfFieldOfViewY", float(tan(activeCamera()->projection().fieldOfViewAngle() / 2.0f)));

		// Projection matrix, for writing to the depth buffer. This
		// creates the input that allows us to use the depth of field
		// effect below.
		Matrix4 projectionMatrix;
		activeCamera()->getProjectUnitMatrix(rd->viewport(), projectionMatrix);
		args.setUniform("projectionMatrix22", projectionMatrix[2][2]);
		args.setUniform("projectionMatrix23", projectionMatrix[2][3]);

		// Textures for the Venice example
		for (int i = 0; i < m_veniceTexture.size(); ++i) {
			args.setUniform(format("iChannel%d", i), m_veniceTexture[i], Sampler::defaults());
		}

		// Set the domain of the shader to the viewport rectangle
		args.setRect(rd->viewport());

		// Call the program in trace.pix for every pixel within the
		// domain, using many threads on the GPU. This call returns
		// immediately on the CPU and the code executes asynchronously
		// on the GPU.
		const String shaderFile =
			m_dirShaders + String(
				// "trace-minimal.pix"
				// "trace-analytic.pix"
				// "trace-raymarch.pix"
				"trace-venice.pix"
			);

		// TODO: This reloads the shader in EVERY FRAME to enable automatic shader reloading.
		// Replace LAUNCH_SHADER_PTR with LAUNCH_SHADER to measure performance and actually ship.
		const shared_ptr<G3D::Shader> shader = G3D::Shader::getShaderFromPattern(shaderFile);
		LAUNCH_SHADER_PTR_WITH_HINT(shader, args, "");
		// LAUNCH_SHADER(shaderFile, args);

		// Post-process special effects
		m_depthOfField->apply(rd, m_framebuffer->texture(0),
			m_framebuffer->texture(Framebuffer::DEPTH),
			activeCamera(), Vector2int16());

	} rd->pop2D();

	swapBuffers();

	rd->clear();

	// Perform gamma correction, bloom, and SSAA, and write to the native window frame buffer
	m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0), 0, 0);
}


//////////////////////////////////////////////////////////////////////


G3D_START_AT_MAIN();

int main(int argc, const char* argv[])
{
	GApp::Settings settings(argc, argv);

	settings.window.caption = "Simple GPU Ray Marcher";
	settings.window.width = 1440;
	settings.window.height = 1080;
	settings.hdrFramebuffer.depthGuardBandThickness = Vector2int16(0, 0);
	settings.hdrFramebuffer.colorGuardBandThickness = Vector2int16(0, 0);

	return App(settings).run();
}

