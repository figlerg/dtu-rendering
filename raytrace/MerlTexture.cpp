// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <iostream>
#include <optix_world.h>
#include "BRDF.h"
#include "MerlTexture.h"
#include "mt_random.h"

using namespace std;
using namespace optix;

void MerlTexture::load(const char* filename)
{
  const unsigned int MERL_SIZE = BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2*3;
  float3 rho_d;
  brdf.resize(MERL_SIZE);
  if(!read_brdf(filename, &brdf[0], MERL_SIZE, rho_d))
    cerr << "Error reading file " << filename << endl;
  printf("test");
  M_matrix = initialize_M(&brdf[0], nr_bins);
  marginal_pdf_matrix = initialize_marginal_densities(&brdf[0],nr_bins);
  max_ratios = find_ratio_bound(marginal_pdf_matrix, nr_bins);
  float theta_r = 0.21 * M_PIf;
  float theta_i = 0.35 * M_PIf;
  float phi_diff = 0.74 * M_PIf;
  printf("\n Szenario1 (theta_r,theta_i,phi_diff) = (%f,%f,%f) => mdf = %f", theta_r, theta_i, phi_diff, read_matrix(marginal_pdf_matrix, theta_r, theta_i,nr_bins));
  width = height = 1;
  channels = 3;
  data = (unsigned char *)malloc(width*height*channels);
  data[0] = static_cast<unsigned char>(rho_d.x*255 + 0.5);
  data[1] = static_cast<unsigned char>(rho_d.y*255 + 0.5);
  data[2] = static_cast<unsigned char>(rho_d.z*255 + 0.5);

  delete [] fdata;
  fdata = new float4[1];
  fdata[0] = make_float4(rho_d);
  tex_handle = SOIL_create_OGL_texture(data, width, height, channels, tex_handle, SOIL_FLAG_INVERT_Y | SOIL_FLAG_TEXTURE_REPEATS);
  tex_target = GL_TEXTURE_2D;
}

float4 MerlTexture::sample_nearest(const float3& texcoord) const
{
  if(!fdata)
    return make_float4(0.0f);
  else
    return fdata[0];
}

float4 MerlTexture::sample_linear(const float3& texcoord) const
{
  return sample_nearest(texcoord);
}

float3 MerlTexture::brdf_lookup(const float3& n, const float3& normalized_wi, const float3& normalized_wo) const
{
  return lookup_brdf_val(&brdf[0], n, normalized_wi, normalized_wo);
}

float3 MerlTexture::importance_sampler(const float3 dir, const float3 n) {

	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees

	Matrix3x3 basis_change = get_trafo(n);

	float3 from_camera = basis_change * dir; // conversion from world to tangent space

	float tmp = acosf(from_camera.z);
	float theta_r = tmp;
	if (theta_r > M_PI_2f) {
		return make_float3(-1000.0f);
	}

	float theta_i = sample_theta_i(theta_r);

	float phi_start = atan2f(from_camera.y, from_camera.x); 
	// we need the original phi from the input dir, 
	// since we assume in tangent space that the input dir has 0 as phi
	float phi_diff = phi_start + sample_phi_diff(theta_r, theta_i);

	float x = sinf(theta_r) * cosf(phi_diff);
	float y = sinf(theta_i) * sinf(phi_diff);
	float z = cos(theta_i);

	float3 o = make_float3(x, y, z);

	// convert with transposed matrix 
	float3 result = basis_change.transpose()* o;

	return result;

}

float MerlTexture::sample_theta_i(const float theta_r) {
	// sample a theta_i from just the ray that comes from camera (r in shader) and normal

	float full_theta_r = M_PI_2f; // 90 degrees

	uint j = (int)(theta_r * nr_bins / full_theta_r) % nr_bins; // idx for theta_r lookup

	float theta_i = -1.0f;
	float f_theta_i = -1.0f;
	float u = -1.0f;
	do {
		u = mt_random();
		theta_i = mt_random() * M_PI_2f;
		f_theta_i = read_matrix(marginal_pdf_matrix, theta_i, theta_r, nr_bins);
	} while (u > f_theta_i / (max_ratios[j] * 2 / M_PIf));

	

	return theta_i;

}

float MerlTexture::sample_phi_diff(float theta_r, float theta_i) {
	
	float phi_diff;
	float u;
	float M;
	float f_sum;

	do {
		u = mt_random();
		phi_diff = mt_random() * M_PI_2f;
		M = read_matrix(M_matrix, theta_i, theta_r, nr_bins);
		float3 f = lookup_brdf_val_2(&brdf[0], theta_i, -1.0f, theta_r, phi_diff);
		f_sum = f.x + f.y + f.z;


	} while (u > f_sum / M);

	return phi_diff;

}

