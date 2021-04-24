/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3D.h>

class CylinderSettings
{
public:
	float height{ 1.0F };
	float radius{ 1.0F };

	const int defaultSteps{ 6 };
	int steps{ defaultSteps };
};

class HeightfieldSettings
{
public:
	GuiLabel* filenameLabel{ nullptr };
	bool validFilename{ false };
	String imageFilename{ "" };

	float scaleX{ 1.0F };
	float scaleY{ 1.0F };
	float scaleZ{ 1.0F };

	const int subdivisionsMin{ 10 };
	const int subdivisionsMax{ 500 };
	int subdivisionsX{ 30 };
	int subdivisionsZ{ 40 };
};

class GlassSettings
{
public:
	int steps{ 4 };
	int quality{ 50 };
};

/** \brief Application framework. */
class App : public GApp {
protected:

	/** Called from onInit */
	void makeGUI();

public:

	App(const GApp::Settings& settings = GApp::Settings());

	void onInit() override;
	void onAI() override;
	void onNetwork() override;
	void onSimulation(RealTime rdt, SimTime sdt, SimTime idt) override;
	void onPose(Array<shared_ptr<Surface> >& posed3D, Array<shared_ptr<Surface2D> >& posed2D) override;

	// You can override onGraphics if you want more control over the rendering loop.
	// void onGraphics(RenderDevice* rd, Array<shared_ptr<Surface> >& surface, Array<shared_ptr<Surface2D> >& surface2D) override;

	void onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& surface3D) override;
	void onGraphics2D(RenderDevice* rd, Array<shared_ptr<Surface2D> >& surface2D) override;

	bool onEvent(const GEvent& e) override;
	void onUserInput(UserInput* ui) override;
	void onCleanup() override;

	void drawDebugShapes() override;

private:
	void generateCylinderOff(float height, float radius, int steps);

	struct FaceIndices {
		int num{ 0 };
		int indices[4];
	};
	void computeCylindricalShape(Array<Point3>& vertices, Array<FaceIndices>& faces, const int steps, const std::vector<Point2>& rings) const;

	void generateHeightfieldOff(float scaleX, float scaleY, float scaleZ);

	void generateGlassOff(int steps, float qualityPercent);
	int getTopFirstBlackLine(shared_ptr<const Image> image) const;
	int getBottomFirstBlackLine(shared_ptr<const Image> image) const;
	int getFirstInLine(shared_ptr<const Image> image, const Color3& color, const Color3& colorPrv, int line, int startingOffset) const;

private:
	CylinderSettings m_cylinderSettings;
	HeightfieldSettings m_heightfieldSettings;
	GlassSettings m_glassSettings;
};
