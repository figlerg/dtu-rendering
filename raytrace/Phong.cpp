// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "HitInfo.h"
#include "Phong.h"

using namespace optix;

#ifndef M_1_PIf
#define M_1_PIf 0.31830988618379067154
#endif

float3 Phong::shade(const Ray& r, HitInfo& hit, bool emit) const
{
  float3 rho_d = get_diffuse(hit);
  float3 rho_s = get_specular(hit);
  float s = get_shininess(hit);

  // Implement Phong reflection here.
  //
  // Input:  r          (the ray that hit the material)
  //         hit        (info about the ray-surface intersection)
  //         emit       (passed on to Emission::shade)
  //
  // Return: radiance reflected to where the ray was coming from
  //
  // Relevant data fields that are available (see Lambertian.h, Ray.h, and above):
  // lights             (vector of pointers to the lights in the scene)
  // hit.position       (position where the ray hit the material)
  // hit.shading_normal (surface normal where the ray hit the material)
  // rho_d              (difuse reflectance of the material)
  // rho_s              (specular reflectance of the material)
  // s                  (shininess or Phong exponent of the material)
  //
  // Hint: Call the sample function associated with each light in the scene.

  float3 final = make_float3(0, 0, 0);

  float3 result;

  for (int i = 0; i < lights.size(); i++) {
	  float3 accumulator = make_float3(0, 0, 0);
	  float3 direction = make_float3(0, 0, 0);
	  float3 intensity = make_float3(0, 0, 0);

	  bool not_occluded = lights[i]->sample(hit.position, direction, intensity);

	  if (not_occluded) {
		  float sthg = dot(normalize(- hit.position), reflect(r.direction, hit.shading_normal));
		  float3 new_term = rho_s * (s + 2) * 0.5 * M_1_PIf * pow(sthg, s);
		  float3 diffuse = (rho_d * M_1_PIf + new_term) * intensity * dot(hit.shading_normal, direction);

		  float cosine = dot(hit.shading_normal, direction);
		  if (cosine > 0) {
			  accumulator += diffuse;
		  }
	  }

	  final += accumulator;
  }
  result = final;



  return result + Emission::shade(r, hit, emit);


  return Lambertian::shade(r, hit, emit);
}
