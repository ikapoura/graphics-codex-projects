#ifndef DEBUG_ANALYTICSPHERE_GLSL
#define DEBUG_ANALYTICSPHERE_GLSL

struct Sphere
{
	vec3  center;
	float radius;
};

bool intersectSphere(Sphere sphere, Ray ray, inout float dist)
{
	vec3  v = ray.origin - sphere.center;
	float b = 2.0f * dot(ray.direction, v);
	float c = dot(v, v) - square(sphere.radius);
	float d = square(b) - 4.0f * c;
	if (d < 0.0f) { return false; }

	float dsqrt = sqrt(d);

	// Choose the first positive intersection
	float t = min(infIfNegative((-b - dsqrt) / 2.0f), infIfNegative((-b + dsqrt) / 2.0f));

	if (t > 0.0f) {
		dist = t;
		return true;
	} else {
		return false;
	}
}

vec3 dbgAnalyticSpherePosition(Ray ray)
{
	Sphere sphere;
	sphere.center = vec3(0.5f, 0.5f, 0.0f);
	sphere.radius = 0.5f;

	float t = 0.0f;
	if (intersectSphere(sphere, ray, t)) {
		vec3 hitPoint = ray.origin + ray.direction * vec3(t);
		return hitPoint;
	}

	return vec3(0.0f);
}

vec3 dbgAnalyticSphereNormal(Ray ray)
{
	Sphere sphere;
	sphere.center = vec3(0.5f, 0.5f, 0.0f);
	sphere.radius = 0.5f;

	float t = 0.0f;
	if (intersectSphere(sphere, ray, t)) {
		vec3 hitPoint = ray.origin + ray.direction * vec3(t);
		vec3 normal = normalize(hitPoint - sphere.center);
		return normal;
	}

	return vec3(0.0f);
}


#endif // DEBUG_ANALYTICSPHERE_GLSL
