// 02576 Rendering Framework
// Written by Jeppe Revall Frisvad, 2010
// Copyright (c) DTU Informatics 2010

#include <optix_world.h>
#include "HitInfo.h"
#include "Volume.h"

using namespace optix;

float3 Volume::shade(const Ray& r, HitInfo& hit, bool emit) const
{
  // If inside the volume, Find the direct transmission through the volume by using
  // the transmittance to modify the result from the Transparent shader.

	//dot(normalize(r.direction), normalize(hit.shading_normal)) > 0 // check this?
	
	//float s = hit.dist; // this makes no sense I think? This is the hit distance from the first ray hitting the volume?
	//float s = 0.0f;

	float3 no_absorption = Transparent::shade(r, hit, emit);

	float3 transmittance = make_float3(1.0f);
	if (dot(normalize(r.direction), normalize(hit.shading_normal)) > 0) {
		transmittance = get_transmittance(hit);
	}

	float3 result = no_absorption * transmittance; // apparently this is component wise according to web. check!


  return result;
}

float3 Volume::get_transmittance(const HitInfo& hit) const
{
  if(hit.material)
  {

    // Compute and return the transmittance using the diffuse reflectance of the material.
    // Diffuse reflectance rho_d does not make sense for a specular material, so we can use 
    // this material property as an absorption coefficient. Since absorption has an effect
    // opposite that of reflection, using 1/rho_d-1 makes it more intuitive for the user.
    float3 rho_d = make_float3(hit.material->diffuse[0], hit.material->diffuse[1], hit.material->diffuse[2]);
	
	// TODO not sure because of dim = 3?
	// What to do for the corner case rho_d == 0 ? some - INF value?
	//float3 transmittance = make_float3(0.0f);
	float x, y, z;
	if (rho_d.x == 0.0f) x = 0.0f;
	else x = 1 / rho_d.x - 1;
	if (rho_d.y == 0.0f) y = 0.0f;
	else y = 1 / rho_d.y - 1;
	if (rho_d.z == 0.0f) z = 0.0f;
	else  z = 1 / rho_d.z - 1;



	float3 transmittance = expf(-hit.dist * make_float3(x, y, z));
	return transmittance;

  }
  return make_float3(1.0f);
}
