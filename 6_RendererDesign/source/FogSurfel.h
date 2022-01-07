#pragma once

#include <G3D/G3D.h>

class FogSurfel : public G3D::Surfel
{
public:
	FogSurfel() = default;
	FogSurfel(const Color3& _color);
	~FogSurfel() = default;
	FogSurfel(const FogSurfel& other) = default;
	FogSurfel(FogSurfel&& other) = default;

	G3D::Color3 finiteScatteringDensity(const Vector3& wi, const Vector3& wo, const ExpressiveParameters& expressiveParameters = ExpressiveParameters()) const override;

	void getImpulses(PathDirection direction, const Vector3& w, ImpulseArray& impulseArray, const ExpressiveParameters& expressiveParameters = ExpressiveParameters()) const override;

	bool scatter(PathDirection pathDirection, const Vector3& w_before, bool russianRoulette, Random& rng, Color3& weight,
		Vector3& w_after, bool& impulseScattered = ignoreBool, float& probabilityHint = ignore,
		const ExpressiveParameters& expressiveParameters = ExpressiveParameters()) const override;

public:
	Color3 color = Color3(1.0);
};
