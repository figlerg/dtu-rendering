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
  printf("\ntesting:\n");
  float theta_r = 1.1;
  float theta_i = 0.35 * M_PIf;
  float phi_diff = 0.74 * M_PIf;
  M_matrix = initialize_M(&brdf[0], nr_bins);
  marginal_pdf_matrix = initialize_marginal_densities(&brdf[0],nr_bins);
  max_ratios = find_ratio_bound(marginal_pdf_matrix, nr_bins);
  //float expected = 0.0f;
  //for (int i = 0; i < 10000; i++) {
	 // float a = sample_theta_i(theta_r);
	 // float b = sample_phi_diff(theta_r, a);
	 // expected += b;
  //}
  //printf("\nexpected = %f\n", expected / (10000));

  pdf_matrix = initialize_pdf(&brdf[0], nr_bins);
  printf("\n");
	for (int j = 0; j < 20; j++) {
		for (int k = 0;k< 20;k++){
			printf("%f", pdf_matrix[3][j][k]);
		}
		printf("\n");
	}
	  //printf("\n");
  //printf("%.3f,", max_ratios[3]);
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

float3 MerlTexture::importance_sampler(const float3 dir, const float3 n, float &pdf_val) {

	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees

	Matrix3x3 basis_change = get_trafo(n);

	float3 from_camera = basis_change * dir; // conversion from world to tangent space

	float tmp = safe_acosf(from_camera.z);
	//printf("\nn = (%f,%f,%f), tmp = %f", n.x, n.y, n.z, tmp);
	float theta_r = tmp;
	if (theta_r > M_PI_2f || theta_r < 0) {
		//printf("\ncaught invalid theta_r in importance sampler!\n");
		return make_float3(-1000.0f);
	}

	float theta_i = sample_theta_i(theta_r);

	float phi_start = atan2f(from_camera.y, from_camera.x); 
	// we need the original phi from the input dir, 
	// since we assume in tangent space that the input dir has 0 as phi
	float phi_diff = sample_phi_diff(theta_r, theta_i);

	float phi_in = phi_start + phi_diff;

	float x = sinf(theta_i) * cosf(phi_in);
	float y = sinf(theta_i) * sinf(phi_in);
	float z = cos(theta_i);

	float3 o = make_float3(x, y, z);

	// convert with transposed matrix 
	float3 result = basis_change.transpose()* o;
	pdf_val = read_matrix_3d(pdf_matrix, theta_i, theta_r, phi_diff, nr_bins);

	return result;

}

float MerlTexture::sample_theta_i(const float theta_r) {
	// sample a theta_i from just the ray that comes from camera (r in shader) and normal

	float full_theta_r = M_PI_2f; // 90 degrees

	uint j = (int)(theta_r * nr_bins / full_theta_r); // idx for theta_r lookup

	//printf("j = %i", j);

	//if (j >= nr_bins) {
	//	printf("j = %i", j);
	//	throw invalid_argument("\nj too high in sample_theta_i");
	//}
	//if (j < 0) {
	//	printf("j = %i", j);
	//	throw invalid_argument("\nj too low in sample_theta_i");
	//}
	

	float theta_i = -1.0f;
	float f_theta_i = -1.0f;
	float u = -1.0f;
	do {
		u = mt_random();
		theta_i = mt_random() * M_PI_2f;
		//u = 0.0f;
		//theta_i = 0.7 * M_PI_2f;
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
		phi_diff = mt_random() * 2 * M_PIf;
		//u = 0.0f;
		//phi_diff = 0.4f * 2* M_PIf;
		M = read_matrix(M_matrix, theta_i, theta_r, nr_bins);
		float3 f = lookup_brdf_val_2(&brdf[0], theta_i, -1.0f, theta_r, phi_diff);
		f_sum = f.x + f.y + f.z;


	} while (u > f_sum / M);

	return phi_diff;

}


float MerlTexture::safe_acosf(float x) const {
	if (x < -1.0) x = -1.0;
	else if (x > 1.0) x = 1.0;
	return acosf(x);
}
