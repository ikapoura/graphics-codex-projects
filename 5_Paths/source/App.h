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

class PinholeCamera
{
public:
	PinholeCamera(const CoordinateFrame& frame, const Projection& projection);

	// x, y, width and height in pixels; P in meters
	Ray getPrimaryRay(float x, float y, int width, int height) const;

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
		int  numLightTrasportPaths{ 1 };
		int  maxScatterEvents{ 2 };
#ifndef G3D_DEBUG
		bool multithreading{ true };
#else
		bool multithreading{ false };
#endif
		float environmentBrightness{ 0.0f };
		GuiDropDownList* resolutionList{ nullptr };
	};

	RayTracer() = delete;
	RayTracer(const Settings& settings, const shared_ptr<Scene>& scene);
	virtual ~RayTracer() = default;

	chrono::milliseconds traceImage(const shared_ptr<Camera>& camera, shared_ptr<Image>& image);

private:
	enum class IntersectionMode {
		Nearest,
		First
	};
	Radiance3 randomColorFromDirection(const Vector3& w) const;

	void rebuildTreeStructureBasedOnLastChange();

	void addEmittedRadiance(const Array<Ray>& rayBuffer, const Array<shared_ptr<Surfel>>& surfelBuffer, const Array<Radiance3>& modulationBuffer, Array<Radiance3>& radianceBuffer) const;
	void sampleDirectLights(const Array<shared_ptr<Surfel>>& surfelBuffer, Array<Ray>& shadowRayBuffer, Array<Biradiance3>& biradianceBuffer) const;
	void addDirectIllumination(const Array<Ray>& rayBuffer, const Array<shared_ptr<Surfel>>& surfelBuffer, const Array<Biradiance3>& biradianceBuffer, const Array<Ray>& shadowRayBuffer,
							   const Array<bool>& lightShadowedBuffer, const Array<Radiance3>& modulationBuffer, Array<Radiance3>& radianceBuffer) const;
	void scatterRays(const Array<shared_ptr<Surfel>>& surfelBuffer, Array<Ray>& rayBuffer, Array<Radiance3>& modulationBuffer) const;

private:
	const Settings m_settings;

	const shared_ptr<Scene>& m_scene;
	shared_ptr<TriTree> m_sceneTriTree;
	RealTime m_lastTreeBuildTime;
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
	/** Allocated by expose and render */
	shared_ptr<Texture> m_result;

	RayTracer::Settings m_rayTraceSettings;
};

