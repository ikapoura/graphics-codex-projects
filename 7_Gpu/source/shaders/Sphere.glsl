#ifndef Sphere_glsl // -*- c++ -*-
#define Sphere_glsl

#include "Surfel.glsl"

struct Sphere {
    Point3      center;
    float       radius;
    Material    material;
};

/* Updates the distance to the first surfel, if closer than distance,
  and returns true on intersection or false if no intersection.
  This uses the implicit equation for the surface of the sphere.
  See Graphics Codex [raySphr] */
bool intersect(Sphere sphere, Point3 P, Vector3 w, inout float distance, inout Surfel surfel) {
     Point3  C = sphere.center;
     float   r = sphere.radius;
     
     Vector3 v = P - C;
     float b = 2.0 * dot(w, v);
     float c = dot(v, v) - square(r);
     float d = square(b) - 4.0 * c;
     if (d < 0.0) { return false; }
     
     float dsqrt = sqrt(d);
     
     // Choose the first positive intersection
     float t = min(infIfNegative((-b - dsqrt) / 2.0),   
                   infIfNegative((-b + dsqrt) / 2.0));

     if (t < distance) { 
         surfel.position = P + w * t;
         surfel.normal   = normalize(surfel.position - C);
         surfel.material = sphere.material;
         
         distance = t;
         return true;
     } else {
         return false;
     }
}


/**
 Signed distance function for a sphere.  Negative inside of the sphere.
*/
float distanceEstimate(Sphere sphere, Point3 P, inout Material material) {
    P -= sphere.center;
    return length(P) - sphere.radius;
}

#endif
