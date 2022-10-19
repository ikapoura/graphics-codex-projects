#ifndef DEBUG_RAYMARCHSPHERE_GLSL
#define DEBUG_RAYMARCHSPHERE_GLSL

struct Sphere
{
	vec3  center;
	float radius;
};

float distanceSphere(Sphere sphere, vec3 point)
{
	return length(point - sphere.center) - sphere.radius;
}


bool raymarchSphere(Sphere sphere, Ray ray, inout float t)
{
	const float maxDistance   = 1e10f;
	const int   maxIterations = 100;
	const float closeEnough   = 1e-2f; // Lowering this value gets us closer to the analytic solution.

	t = 0.0f;
	for (int i = 0; i < maxIterations; ++i) {
		float dt = distanceSphere(sphere, ray.origin + ray.direction * t);
		t += dt;

		if (dt < closeEnough) {
			return true;
		}
	}

	return false;
}

vec3 dbgRaymarchSpherePosition(Ray ray)
{
	Sphere sphere;
	sphere.center = vec3(0.5f, 0.5f, 0.0f);
	sphere.radius = 0.5f;

	float t = 0.0f;
	if (raymarchSphere(sphere, ray, t)) {
		vec3 hitPoint = ray.origin + ray.direction * vec3(t);
		return hitPoint;
	}

	return vec3(0.0f);
}

vec3 dbgRaymarchSphereNormal(Ray ray)
{
	Sphere sphere;
	sphere.center = vec3(0.5f, 0.5f, 0.0f);
	sphere.radius = 0.5f;

	float t = 0.0f;
	if (raymarchSphere(sphere, ray, t)) {
		vec3 hitPoint = ray.origin + ray.direction * vec3(t);
		vec3 normal = normalize(hitPoint - sphere.center);
		return normal;
	}

	return vec3(0.0f);
}

#endif // DEBUG_RAYMARCHSPHERE_GLSL
