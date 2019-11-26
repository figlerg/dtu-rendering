// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "HitInfo.h"
#include "int_pow.h"
#include "GlossyVolume.h"

using namespace optix;

#ifndef M_1_PIf
#define M_1_PIf 0.31830988618379067154
#endif

float3 GlossyVolume::shade(const Ray& r, HitInfo& hit, bool emit) const
{
  // Compute the specular part of the glossy shader and attenuate it
  // by the transmittance of the material if the ray is inside (as in
  // the volume shader).

	if (hit.trace_depth >= max_depth)
		return make_float3(0.0f);

	float R;
	Ray reflected, refracted;
	HitInfo hit_reflected, hit_refracted;
	tracer->trace_reflected(r, hit, reflected, hit_reflected);
	tracer->trace_refracted(r, hit, refracted, hit_refracted, R);


	return shade_new_ray(refracted, hit_refracted) + Phong::shade(r, hit, emit);
	return R * shade_new_ray(reflected, hit_reflected) + (1.0f - R) * shade_new_ray(refracted, hit_refracted) + Phong::shade(r, hit, emit);


  return Volume::shade(r, hit, emit)+Phong::shade(r,hit,emit);

}
