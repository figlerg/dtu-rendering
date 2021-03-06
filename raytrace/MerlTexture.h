// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2019
// Copyright (c) DTU Informatics 2019

#ifndef MERLTEXTURE_H
#define MERLTEXTURE_H

#include <vector>
#include <optix_world.h>
#include "Texture.h"

class MerlTexture : public Texture
{
public:
  MerlTexture() : Texture() { }
  ~MerlTexture() { free(static_cast<void*>(data)); delete [] fdata; }

  // Load texture from file
  void load(const char* filename);

  // Look up the texel using texture space coordinates
  virtual optix::float4 sample_nearest(const optix::float3& texcoord) const;
  virtual optix::float4 sample_linear(const optix::float3& texcoord) const;

  optix::float3 brdf_lookup(const optix::float3& n, const optix::float3& normalized_wi, const optix::float3& normalized_wo) const;

  // my functions for project

  optix::float3 importance_sampler(const optix::float3 dir, const optix::float3 n, float &pdf_val);
  float sample_theta_i(const float theta_r);
  float sample_phi_diff(float theta_r, float theta_i);
  float safe_acosf(float x) const;

  // MINE
  float*** pdf_matrix; // saved values of conditional probabilities for triple of all three angles



protected:
  // Pointers to image data
  std::vector<float> brdf;
  optix::float3 rho_d;

  // new members for project!
  int nr_bins = 100;
  float** M_matrix;
  float** marginal_pdf_matrix;
  float* max_ratios;
};

#endif // MERLTEXTURE_H
