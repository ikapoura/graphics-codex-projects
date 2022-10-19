#ifndef DEBUG_EYEDIRECTIONS_GLSL
#define DEBUG_EYEDIRECTIONS_GLSL

#include <g3dmath.glsl>

vec3 dbgEyeDirections(mat4x3 cameraToWorldMatrix, float tanHalfFieldOfViewY)
{
	// Generate an eye ray in camera space
	vec3 direction;
	direction.xy = (gl_FragCoord.xy - g3d_FragCoordExtent / 2.0) * Vector2(1, -1);
	direction.z  = g3d_FragCoordExtent.y / (-2.0 * tanHalfFieldOfViewY);
	direction = normalize(direction);

	// Transform the ray to world space
	direction = (cameraToWorldMatrix * vec4(direction, 0.0f)).xyz;

	return (direction + vec3(1.0f)) / vec3(2.0f);
}

#endif // DEBUG_EYEDIRECTIONS_GLSL
