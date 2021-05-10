/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3D.h>

class PostProcess
{
public:
	Color3 pixelValue(const Radiance3& L, const float k, const float gamma);
};

class Intersector
{
public:
	static bool raySphereIntersect(const Point3& P, const Vector3& w, const Sphere& s, float& t);
	static bool rayTriangleIntersect(const Point3& P, const Vector3& w, const Point3 V[3], float b[3], float& t);
};

class PinholeCamera
{
public:
	PinholeCamera(const CoordinateFrame& frame, const Projection& projection);

	// x, y, width and height in pixels; P in meters
	void getPrimaryRay(float x, float y, int width, int height, Point3& P, Vector3& w) const;

protected:
	CoordinateFrame m_frame;
	Projection m_projection;
};

class SpherePrimitive
{
public:
	Sphere sphere;
	Color3 color;
};

class RayTracer
{
public:
	class Settings {
	public:
		bool addFixedPrimitives{ false };
		int  indirectRaysPerPixel{ 1 };
		bool multithreading{ true };
		GuiDropDownList* resolutionList{ nullptr };
	};

	RayTracer() = delete;
	RayTracer(const Settings& settings, const shared_ptr<Scene>& scene);
	virtual ~RayTracer() = default;

	void addFixedSphere(const Point3& center, float radius, const Color3& color);

	chrono::milliseconds traceImage(const shared_ptr<Camera>& camera, shared_ptr<Image>& image);

private:
	enum class IntersectionMode {
		Nearest,
		First
	};
	shared_ptr<UniversalSurfel> findIntersection(const Point3& X, const Vector3& wi, const float maxDistance, const IntersectionMode mode) const;
	void intersectFixedPrimitives(const Point3& X, const Vector3& wi, const IntersectionMode mode, shared_ptr<UniversalSurfel>& result, float& t) const;
	void intersectTriangulatedSurfaces(const Point3& X, const Vector3& wi, const IntersectionMode mode, shared_ptr<UniversalSurfel>& result, float& t) const;

	bool visibleFromLight(const shared_ptr<UniversalSurfel>& s, const Point3& lightPosition) const;

	Radiance3 L_i(const shared_ptr<UniversalSurfel>& s, const Vector3& wi) const;
	Radiance3 L_o(const shared_ptr<UniversalSurfel>& s, const Vector3& wo) const;
	Radiance3 randomColorFromDirection(const Vector3& w) const;

private:
	const Settings m_settings;

	const shared_ptr<Scene>& m_scene;
	Array<shared_ptr<Surface>> m_sceneSurfaces;
	shared_ptr<TriTree> m_sceneTriTree;

	Array<SpherePrimitive> m_fixedSpheres;
};


 /** \brief Application framework. */
class App : public GApp
{
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

	Vector2int32 resolution() const;

	String durationToString(chrono::milliseconds duration) const;

private:
	void render();

private:
	RayTracer::Settings m_rayTraceSettings;
};

