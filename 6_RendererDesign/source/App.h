/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3D.h>

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

class PathBuffers
{
public:
	PathBuffers() = default;
	~PathBuffers() = default;

	int size() const;
	bool empty() const;

	void resize(size_t size);

	void resizeShadowRays(size_t size);

	void clearRedundantSurfels(float minModulationSum);
	void clearRedundantShadowRays();

	void postProcessShadowResults();

public:
	// x, y coordinates of the remaining surfels
	Array<Point2> imageCoordinates;
	Array<Radiance3> modulation;
	Array<Ray> rays;
	Array<shared_ptr<Surfel>> surfels;

	// Whether a surfel is in shadow from the sampled light.
	Array<bool> surfelShadowed;

	Array<Ray> shadowRays;
	Array<int> shadowRayOriginalIndex; // The index for each shadow ray before culling.
	Array<bool> shadowRayResult;
	Array<Biradiance3> biradiance;

private:
	void fastRemove(int index);
	void fastRemoveShadows(int index);
};

class RayTracer
{
public:
	class Settings {
	public:
		int  numLightTrasportPaths{ 1 };
#ifndef G3D_DEBUG
		int  maxScatterEvents{ 1 };
#else
		int  maxScatterEvents{ 4 };
#endif
#ifndef G3D_DEBUG
		bool multithreading{ true };
#else
		bool multithreading{ false };
#endif
		float participatingMediumDensity{ 0.01f };
		float environmentBrightness{ 0.0f };
		GuiDropDownList* resolutionList{ nullptr };
	};

	RayTracer() = delete;
	RayTracer(const Settings& settings, const shared_ptr<Scene>& scene);
	virtual ~RayTracer() = default;

	chrono::milliseconds traceImage(const shared_ptr<Camera>& camera, shared_ptr<Image>& image);

private:
	Radiance3 randomColorFromDirection(const Vector3& w) const;

	void rebuildTreeStructureBasedOnLastChange();

	void initializeTransportPaths(const PinholeCamera& camera, const Point2int32& imageSize, int pathIdx,
								  PathBuffers& buffers, shared_ptr<Image>& pathWeightImage) const;

	void processParticipatingMedium(Array<Ray>& rays, Array<shared_ptr<Surfel>>& surfels, float participatingMediumDensity) const;

	void addEmittedRadiance(PathBuffers& buffers, shared_ptr<Image>& image) const;
	void sampleDirectLights(PathBuffers& buffers, const Array<shared_ptr<Light>>& directLightsBuffer, const Point2int32& imageSize, int pathIdx, int scatterIdx) const;
	void addDirectIllumination(PathBuffers& buffers, shared_ptr<Image>& image) const;
	void scatterRays(PathBuffers& buffers) const;

	Point3 lightImportanceSampling(const Array<shared_ptr<Light>>& directLightsBuffer, const shared_ptr<Surfel>& surfel, int pixelIdx, int pathIdx, Random& random,
		int& selectedLightIdx, float& selectedLightWeight, float& selectedLightAreaTimesPdfValue) const;

	Ray generateShadowRay(const Point3& surfelPos, const Vector3& surfelGNormal, const Point3& lightPos) const;

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

