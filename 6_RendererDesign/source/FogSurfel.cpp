#include "FogSurfel.h"

FogSurfel::FogSurfel(const Color3& _color) :
	color(_color)
{
}

G3D::Color3 FogSurfel::finiteScatteringDensity(const Vector3& wi, const Vector3& wo, const ExpressiveParameters& expressiveParameters /*= ExpressiveParameters()*/) const
{
	return Color3(0.5f / G3D::pif());
}

void FogSurfel::getImpulses(PathDirection direction, const Vector3& w, ImpulseArray& impulseArray, const ExpressiveParameters& expressiveParameters /*= ExpressiveParameters()*/) const
{
}

bool FogSurfel::scatter(PathDirection pathDirection, const Vector3& w_before, bool russianRoulette, Random& rng, Color3& weight,
	Vector3& w_after, bool& impulseScattered /*= ignoreBool*/, float& probabilityHint /*= ignore*/, const ExpressiveParameters& expressiveParameters /*= ExpressiveParameters()*/) const
{
	rng.sphere(w_after.x, w_after.y, w_after.z);
	weight = color;
	return true;
}
