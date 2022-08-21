#ifndef Surfel_glsl // -*- c++ -*-
#define Surfel_glsl

/** Blinn-Phong material. */
struct Material {
    /** Base color, applied to Lambertian term for non-metals and glossy term for
        metals.*/
    Color3      color;

    /** 0 = dielectric, 1 = conductor.  Intermediate models interpolate between the models. 
        In this application, currently drives the Lambertian magnitude but may later
        be adjusted to affect Fresnel or other factors. */
    float       metal;

    float       glossyCoefficient;

    /** Higher = smoother. 1e7 = almost mirror-like, 50 = dull plastic */
    float       glossyExponent;
};


struct Surfel {
    Point3      position;
    Vector3     normal;
    Material    material;
};


const Material defaultMaterial = Material(Color3(0.9), 0.5, 1.0, 500);

#endif
