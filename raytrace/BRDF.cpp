#include <fstream>
#include <optix_world.h>
#include "sampler.h"
#include "BRDF.h"

using namespace std;
using namespace optix;

// Rotate vector around an axis
float3 rotate_vector(const float3& vector, const float3& axis, float angle)
{
  // Implement quaternion-based rotation by an angle of a vector around an axis
  // Hint: RixShadingUtils.h has many useful helper functions (such as RixSinCos, Cross, Dot, and Normalize).

	return cosf(angle) * vector + sinf(angle) * cross(axis, vector) + (1 - cosf(angle)) * dot(axis, vector) * axis;
	
	//return make_float3(0.0f);
}

// Convert vectors in tangent space to half/difference coordinates
void vectors_to_half_diff_coords(const float3& in, const float3& out,
  float& theta_half, float& phi_half, float& theta_diff, float& phi_diff)
{
  // compute halfway vector
	float3 h = (in + out) / length(in + out);

  // compute theta_half, phi_half
	theta_half = acosf(h.z);
	phi_half = atan2(h.y, h.x);


  // compute diff vector
	float3 e3 = make_float3(0, 0, 1);
	float3 inner_rotated = rotate_vector(in, e3, -phi_half);

	float3 e2 = make_float3(0, 1, 0);
	float3 d = rotate_vector(inner_rotated, e2, -theta_half);

  // compute theta_diff, phi_diff	
	theta_diff = acosf(d.z);
	phi_diff = atan2(d.y, d.x);
}

float3 lookup_brdf_val(const float* brdf, const float3& n, const float3& normalized_wi, const float3& normalized_wo)
{
  // Transform vectors to tangent space and use the vectors_to_half_diff_coords function
  // to get input for the other version of the lookup_brdf_val function.
  // Hint: There is a built-in member function for 3-vectors called CreateOrthonormalBasis
  //       and a function called RixChangeBasisTo for change of basis.
	// NEUE ERSATZ ONB FUNKTION
	
	float3 t = make_float3(0.0f);
	float3 b = make_float3(0.0f);

	onb(n, t, b);
	
	Matrix3x3 basis_change = Matrix3x3();
	basis_change.setRow(0,t);
	basis_change.setRow(1,b);
	basis_change.setRow(2,n);

	float3 i = basis_change * normalized_wi;
	float3 o = basis_change * normalized_wo;

	float theta_half, phi_half, theta_diff, phi_diff;
	vectors_to_half_diff_coords(i, o, theta_half, phi_half, theta_diff, phi_diff);

	return lookup_brdf_val(brdf, theta_half,phi_half,theta_diff,phi_diff);

	//return make_float3(0.0f); // BY ME BECAUSE OF DEBUGGING PROBLEM
}

float3 integrate_brdf(float* brdf, unsigned int N)
{
  float3 sum = make_float3(0.0f);
  float3 n = make_float3(0.0f, 0.0f, 1.0f);

  for (int i = 0; i < N; i++) {
	  float3 wi = sample_cosine_weighted(n);
	  float3 wo = sample_cosine_weighted(n);
	  float3 f_rk = lookup_brdf_val(brdf, n, wi, wo);

	  sum += f_rk * M_1_PIf/N;	 // reingezogen, siehe formel für rho merlbrdf.pdf s4
  }


  // Implement Monte Carlo integration to estimate the bihemisphercial diffuse reflectance (rho).
  // Use N as the number of samples.
  // Hint: Pass the random number generator seed t to the sample_cosine_hemisphere function.

  return sum;
}

// Read BRDF data
bool read_brdf(const char *filename, float* brdf, unsigned int size, float3& rho_d)
{
  FILE* f = 0;
  if(fopen_s(&f, filename, "rb"))
    return false;

  int dims[3];
  fread(dims, sizeof(int), 3, f);
  unsigned int n = 3*dims[0]*dims[1]*dims[2];
  if(size != n)
    return false;

  double* brdf_d = new double[n];
  fread(&brdf_d[0], sizeof(double), n, f);
  for(unsigned int i = 0; i < size; ++i)
    brdf[i] = static_cast<float>(brdf_d[i]);
  fclose(f);
  delete[] brdf_d;

  rho_d = integrate_brdf(brdf, 100000);
  return true;
}
