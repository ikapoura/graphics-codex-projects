#ifndef DEBUG_CHECKERBOARD_GLSL
#define DEBUG_CHECKERBOARD_GLSL

vec3 dbgCheckerboardPattern()
{
	ivec2 pixelCoordinates = ivec2(gl_FragCoord.xy);

	const int TILE_SIZE = 16;
	const int MAX_SIZE = 2 * TILE_SIZE;

	vec3 pixelColor;

	if ((pixelCoordinates.x % MAX_SIZE < TILE_SIZE)) {
		if ((pixelCoordinates.y % MAX_SIZE < TILE_SIZE)) {
			pixelColor = vec3(0.0);
		} else {
			pixelColor = vec3(1.0);
		}
	} else {
		if ((pixelCoordinates.y % MAX_SIZE < TILE_SIZE)) {
			pixelColor = vec3(1.0);
		} else {
			pixelColor = vec3(0.0);
		}
	}

	return pixelColor;
}

#endif // DEBUG_CHECKERBOARD_GLSL
