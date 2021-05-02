/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3D.h>

class PostProcess {
public:
	Color3 pixelValue(const Radiance3& L, const float k, const float gamma);
};

class Intersector {
public:
	bool rayTriangleIntersect(const Point3& P, const Vector3& w, const Point3 V[3], float b[3], float& t);
};

class PinholeCamera {
public:
	PinholeCamera(float z_near, float verticalFieldOfView);

	// x, y, width and height in pixels; P in meters
	void getPrimaryRay(float x, float y, int width, int height, Point3& P, Vector3& w) const;

protected:
	// Negative
	float m_zNear;

	float m_verticalFieldOfView;
};

class BRDF {
public:
	BRDF() = default;
	virtual ~BRDF() = default;

	Radiance3 L_i(const Point3& X, const Vector3& wi, const shared_ptr<UniversalSurfel>& s) const;
};

class RayTraceSettings {
public:
	bool addFixedPrimitives{ false };
	int  indirectRaysPerPixel{ 1 };
	bool multithreading{ false };
	GuiDropDownList* resolutionList{ nullptr };
};

class RayTracer {
public:
	RayTracer() = delete;
	RayTracer(shared_ptr<BRDF> brdf);
	virtual ~RayTracer() = default;

	void render(const PinholeCamera& camera, shared_ptr<Image>& image) const;

	const shared_ptr<UniversalSurfel>& findFirstIntersection(const Point3& X, const Vector3& wi) const;

	RayTraceSettings& settings();
	const RayTraceSettings& settings() const;

	Vector2int32 resolution() const;

private:
	RayTraceSettings m_settings;
	shared_ptr<BRDF> m_brdf;
};


 /** \brief Application framework. */
class App : public GApp {
public:

	App(const GApp::Settings& settings = GApp::Settings());

	virtual void onInit() override;
	virtual void onAI() override;
	virtual void onNetwork() override;
	virtual void onSimulation(RealTime rdt, SimTime sdt, SimTime idt) override;
	virtual void onPose(Array<shared_ptr<Surface> >& posed3D, Array<shared_ptr<Surface2D> >& posed2D) override;

	// You can override onGraphics if you want more control over the rendering loop.
	// virtual void onGraphics(RenderDevice* rd, Array<shared_ptr<Surface> >& surface, Array<shared_ptr<Surface2D> >& surface2D) override;

	virtual void onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& surface3D) override;
	virtual void onGraphics2D(RenderDevice* rd, Array<shared_ptr<Surface2D> >& surface2D) override;

	virtual bool onEvent(const GEvent& e) override;
	virtual void onUserInput(UserInput* ui) override;
	virtual void onCleanup() override;

	virtual void drawDebugShapes() override;

protected:
	/** Called from onInit */
	void makeGUI();

private:
	void render();

private:
	unique_ptr<RayTracer> m_rayTracer;
};
