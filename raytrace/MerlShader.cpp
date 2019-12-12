// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "mt_random.h"
#include "sampler.h"
#include "HitInfo.h"
#include "ObjMaterial.h"
#include "MerlTexture.h"
#include "MerlShader.h"

using namespace optix;

#ifndef M_1_PIf
#define M_1_PIf 0.31830988618379067154
#endif

float3 MerlShader::shade(const Ray& r, HitInfo& hit, bool emit) const
{
  if(hit.trace_depth >= max_depth)
    return make_float3(0.0f);

  const ObjMaterial* m = hit.material;
  MerlTexture* tex = brdfs && m && m->has_texture ? (*brdfs)[m->tex_name] : 0;
  float3 rho_d = get_diffuse(hit);
  float3 result = make_float3(0.0f);

  // Implement a path tracing shader here.
  //
  // Input:  r          (the ray that hit the material)
  //         hit        (info about the ray-surface intersection)
  //         emit       (passed on to Emission::shade)
  //
  // Return: radiance reflected to where the ray was coming from
  //
  // Relevant data fields that are available:
  // tracer             (pointer to ray tracer)
  //
  // Hint: Use the function tex->brdf_lookup(...) to retrieve the value of
  //       the measured BRDF for a given light-view configuration. Ensure
  //       that tex is non-zero and has texture before using it.

  if (!tex) return result; // if it's either 0 or has no texture this is true

  float P_r = (rho_d.x + rho_d.y + rho_d.z) / 3; // probability of reflection? 
  float xi = mt_random();

  
  if (xi <= P_r) {
	  //float3 dir = sample_cosine_weighted(hit.shading_normal);
	  float3 dir = tex->importance_sampler(-normalize(r.direction), normalize(hit.shading_normal));
	  if (dir.x == -1000.0f) {
		  return result;
	  }

	  Ray r_out = Ray(hit.position, dir, 0, 1e-04, RT_DEFAULT_MAX);
	  HitInfo hit_out = HitInfo();

	  hit_out.trace_depth = hit.trace_depth + 1;
	  float3 reflection = shade_new_ray(r_out, hit_out, emit);

	  float3 brdf_val = tex->brdf_lookup(hit.shading_normal, -r.direction, dir);

	  result = M_PIf * reflection * brdf_val / P_r;
  }
  // else it stays zero because of absorption
  
  //result = tex->brdf_lookup(hit.shading_normal,r.direction,)




  return result;
}
