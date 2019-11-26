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

	//return Volume::shade(r, hit, emit);

	//float3 rho_d = get_diffuse(hit);
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

	float3 result = Volume::shade(r, hit, emit);
	for (auto light:lights) {
		float3 light_intensity;
		float3 light_direction;
		if (light->sample(hit.position, light_direction, light_intensity)) {
			float3 direction_to_observer = normalize(reflect(-light_direction, hit.shading_normal));
			result += rho_s * ((s + 2) / (2 * M_PIf)) * powf(fmax(dot(direction_to_observer, -r.direction), 0), (int)s) * light_intensity *
				fmax(dot(light_direction, hit.shading_normal), 0);
		}
	}




	//for (int i = 0; i < lights.size(); i++) {
	//	float3 accumulator = make_float3(0, 0, 0);
	//	float3 direction = make_float3(0, 0, 0);
	//	float3 intensity = make_float3(0, 0, 0);

	//	bool not_occluded = lights[i]->sample(hit.position, direction, intensity);

	//	if (not_occluded) {
	//		float sthg = fmax(dot(normalize(-r.direction), reflect(-direction, hit.shading_normal)), 0);


	//		float3 new_term = rho_s * (s + 2) * 0.5 * M_1_PIf * pow(sthg, s);
	//		float3 diffuse = (rho_d * M_1_PIf + new_term) * intensity * dot(hit.shading_normal, direction);
	//		

	//		float cosine = dot(hit.shading_normal, direction);
	//		if (cosine > 0) {
	//			float3 direction_to_observer = reflect(-direction, hit.shading_normal);
	//			
	//			float3 transmittance = get_transmittance(hit);
	//			intensity *= transmittance;

	//			//accumulator += rho_s * ((s + 2) / (2 * M_PIf)) * powf(fmax(dot(direction_to_observer, -r.direction), 0), (int)s) * intensity *
	//			//fmax(dot(direction, hit.shading_normal), 0);



	//			//accumulator += diffuse;
	//		}
	//	}

	//	final += accumulator;
	//}
	//result = final;



	return result + Emission::shade(r, hit, emit);


  return Volume::shade(r, hit, emit)+Phong::shade(r,hit,emit);

}



