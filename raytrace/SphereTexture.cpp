// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "SphereTexture.h"

using namespace std;
using namespace optix;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float4 SphereTexture::sample_nearest(const float3& d) const
{
  float u, v;
  project_direction(d, u, v);
  return Texture::sample_nearest(make_float3(u, v, 0.0f));
}

float4 SphereTexture::sample_linear(const float3& d) const
{
  float u, v;
  project_direction(d, u, v);
  return Texture::sample_linear(make_float3(u, v, 0.0f));
}

void SphereTexture::project_direction(const float3& d, float& u, float& v) const
{
  // Implement the angular map from direction to texture uv-coordinates.
  // Remember to handle the singularity.
	// from slides 10-7: project lookup direction to the texture coordinates

	float r = (d.x==0.0f && d.y == 0.0f) ? 0.0f : acosf(-d.z) / (2 * M_PI * sqrtf(powf(d.x, 2) + powf(d.y, 2)));
	u = 0.5 + r * d.x;
	v = 0.5 + r * d.y;

}
