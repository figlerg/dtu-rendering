// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "mt_random.h"
#include "sampler.h"
#include "HitInfo.h"
#include "MCGlossy.h"

using namespace optix;

#ifndef M_1_PIf
#define M_1_PIf 0.31830988618379067154
#endif

float3 MCGlossy::shade(const Ray& r, HitInfo& hit, bool emit) const
{
  if(hit.trace_depth >= max_depth)
    return make_float3(0.0f);

  float3 rho_d = get_diffuse(hit);
  float3 result = make_float3(0.0f,0.0f,0.0f);



  // Implement a path tracing shader here.
  //
  // Input:  r          (the ray that hit the material)
  //         hit        (info about the ray-surface intersection)
  //         emit       (passed on to Emission::shade)
  //
  // Return: radiance reflected to where the ray was coming from
  //
  // Relevant data fields that are available (see Mirror.h and HitInfo.h):
  // max_depth          (maximum trace depth)
  // tracer             (pointer to ray tracer)
  // hit.trace_depth    (number of surface interactions previously suffered by the ray)
  //
  // Hint: Use the function shade_new_ray(...) to pass a newly traced ray to
  //       the shader for the surface it hit.

  float R;
  Ray reflected, refracted;
  HitInfo hit_reflected, hit_refracted;
  tracer->trace_reflected(r, hit, reflected, hit_reflected);
  tracer->trace_refracted(r, hit, refracted, hit_refracted, R);

  float3 new_dir = sample_cosine_weighted(hit.shading_normal);
  Ray new_ray = Ray(hit.position, new_dir, 0, 1e-4, r.tmax);
  HitInfo new_hit = HitInfo();
  new_hit.trace_depth = hit.trace_depth + 1;
  new_hit.ray_ior = hit.ray_ior;
  
  bool has_hit = tracer->trace_to_closest(new_ray, new_hit);

  //float3 test_dir = sample_cosine_weighted(hit.shading_normal);

  // return R * shade_new_ray(reflected, hit_reflected) + (1.0f - R) * shade_new_ray(refracted, hit_refracted) + Phong::shade(r, hit, emit);

  if (has_hit) {
	  throw;
	  return result + Phong::shade(r, hit, emit) + shade_new_ray(new_ray, new_hit, false);
  }
  else return result + Phong::shade(r, hit, emit);
}
