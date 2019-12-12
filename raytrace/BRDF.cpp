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

	return cosf(angle) * vector + sinf(angle) * cross(axis, vector) + (1.0f - cosf(angle)) * dot(axis, vector) * axis;
	
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
	phi_half = atan2f(h.y, h.x);


  // compute diff vector
	float3 e3 = make_float3(0, 0, 1);
	float3 inner_rotated = rotate_vector(in, e3, -phi_half);

	float3 e2 = make_float3(0, 1, 0);
	float3 d = rotate_vector(inner_rotated, e2, -theta_half);

  // compute theta_diff, phi_diff	
	theta_diff = acosf(d.z);
	phi_diff = atan2f(d.y, d.x);
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

Matrix3x3 get_trafo(float3 n) {
	float3 t = make_float3(0.0f);
	float3 b = make_float3(0.0f);

	onb(n, t, b);

	Matrix3x3 basis_change = Matrix3x3();
	basis_change.setRow(0, t);
	basis_change.setRow(1, b);
	basis_change.setRow(2, n);
	
	return basis_change;
}

float3 lookup_brdf_val_2(const float* brdf, const float theta_i, const float phi_i, const float theta_o, const float phi_o)
{
	// to be conform with notation in paper for brdf project, we define a new lookup function that works just with the 
	// angles and acts as an interface between our parameterisation and the paper's.
	// we simply get both vectors from the angles and use vectors_to_half_diff from there (isotropic! rot around z doesnt matter)

	// phi_i not used

	float3 i = make_float3(sinf(theta_i), 0, cosf(theta_i));
	float3 o = make_float3(sinf(theta_o) * cosf(phi_o), sinf(theta_o) * sinf(phi_o), cosf(theta_o));

	float theta_half, phi_half, theta_diff, phi_diff;
	vectors_to_half_diff_coords(i, o, theta_half, phi_half, theta_diff, phi_diff);

	return lookup_brdf_val(brdf, theta_half, phi_half, theta_diff, phi_diff);
}



float3 integrate_brdf(float* brdf, unsigned int N)
{
  float3 sum = make_float3(0.0f);
  float3 n = make_float3(0.0f, 0.0f, 1.0f);

  for (int i = 0; i < N; i++) {
	  float3 wi = sample_cosine_weighted(n);
	  float3 wo = sample_cosine_weighted(n);
	  float3 f_rk = lookup_brdf_val(brdf, n, wi, wo);

	  sum += f_rk;	 // formula for rho, merlbrdf.pdf s4
  }

  sum *= M_PIf / (float) N;

  printf("\n rho_d = ( %f, %f, %f )\n", sum.x, sum.y, sum.z);

  float3 sample1 = normalize(make_float3(0.0f, 0.3f, 1.0f));
  float3 sample2 = normalize(make_float3(0.0f, -0.3f, 1.0f));
  float3 test_lookup = lookup_brdf_val(brdf, n, sample1, sample2);
  printf("\n lookup : %f, %f, %f", test_lookup.x, test_lookup.y, test_lookup.z);


  // Implement Monte Carlo integration to estimate the bihemisphercial diffuse reflectance (rho).
  // Use N as the number of samples.
  // Hint: Pass the random number generator seed t to the sample_cosine_hemisphere function.

  return sum;
}



// Read BRDF data
bool read_brdf(const char *filename, float* brdf, unsigned int size, float3& rho_d)
{
  FILE* f = 0;
  //if(fopen_s(&f, filename, "rb"))
  f = fopen(filename, "rb");
  if(!f)
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

  rho_d = integrate_brdf(brdf, 1000000);
  return true;
}


float** initialize_M(const float* brdf, int n_bins) {
	float** result = new float* [n_bins];
	for (int i = 0; i < n_bins; i++) {
		result[i] = new float[n_bins];
	}

	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees
	float full_phi_diff = M_PIf; // 180 degrees because the values are the same for the rest  of the circle (hopefully)

	// this order is kept in loops: from outer to inner

	float theta_step = full_theta_i / n_bins; // same for both thetas
	float phi_diff_step = full_phi_diff / n_bins;

#pragma omp parallel for

	for (int i = 0; i < n_bins; i++) {

		for (int j = 0; j < n_bins; j++) {
			float max = 0.0f;

			for (int k = 0; k < n_bins; k++) {
				float theta_i = i * theta_step;
				float theta_r = j * theta_step;
				float phi_diff = k * phi_diff_step;

				float3 f_rgb = lookup_brdf_val_2(brdf, theta_i, -1.0f, theta_r, phi_diff);
				float f_sum = f_rgb.x + f_rgb.y + f_rgb.z;
				f_sum *= cosf(theta_i);
				max = (max < f_sum) ? f_sum : max;
			}

			result[i][j] = max;
		}
		printf("%i", i);

	}

	printf("M-Matrix for importance sampling initialized");
	return result;
}

float** initialize_marginal_densities(const float* brdf, int n_bins) {
	// same bins and structure as above, but now we save marginal densities for given theta_i and theta_r
	float** result = new float* [n_bins];
	for (int i = 0; i < n_bins; i++) {
		result[i] = new float[n_bins];
	}

	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees
	float full_phi_diff = M_PIf; // 180 degrees because the values are the same for the rest  of the circle (hopefully)

	// this order is kept in loops: from outer to inner

	float theta_step = full_theta_i / n_bins; // same for both thetas
	float phi_diff_step = full_phi_diff / n_bins;

#pragma omp parallel for

	for (int i = 0; i < n_bins; i++) {

		for (int j = 0; j < n_bins; j++) {
			float max = 0.0f;

			float theta_i = i * theta_step;
			float theta_r = j * theta_step;

			result[i][j] = get_marginal_density(brdf,theta_r,theta_i,n_bins);
		}
		printf("%i", i);
	}
	printf("^\nMatrix with marginal densities for importance sampling initialized");
	return result;
}


float get_marginal_density(const float* brdf, const float theta_r, const float theta_i, int n_bins) {

	float** result = new float* [n_bins];
	for (int i = 0; i < n_bins; i++) {
		result[i] = new float[n_bins];
	}

	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees
	float full_phi_diff = M_PIf; // 180 degrees because the values are the same for the rest  of the circle (hopefully)

	// this order is kept in loops: from outer to inner

	float theta_step = full_theta_i / n_bins; // same for both thetas
	float phi_diff_step = full_phi_diff / n_bins;

	float enumerator = 0.0f;
	float denominator = 0.0f;

			
	for (int k = 0; k < n_bins; k++) {
		float phi_diff = k * phi_diff_step;

		float3 f_rgb = lookup_brdf_val_2(brdf, theta_i, -1.0f, theta_r, phi_diff);
		float f_sum = f_rgb.x + f_rgb.y + f_rgb.z;
		f_sum *= cosf(theta_i);

		enumerator += f_sum;

	}
	enumerator *= phi_diff_step; // my definition of this discrete integral

	for (int j = 0; j < n_bins; j++) {

		for (int k = 0; k < n_bins; k++) {
			float phi_diff = k * phi_diff_step;
			float theta_i_inner = j * theta_step;

			float3 f_rgb = lookup_brdf_val_2(brdf, theta_i_inner, -1.0f, theta_r, phi_diff);
			float f_sum = f_rgb.x + f_rgb.y + f_rgb.z;
			f_sum *= cosf(theta_i_inner);

			denominator += f_sum;
		}

	}

	denominator *= theta_step * phi_diff_step;

	return enumerator / denominator;

}

float read_matrix(float** matrix, const float theta_i, const float theta_r, const int nr_bins) {
	float full_theta_i = M_PI_2f; // 90 degrees
	float full_theta_r = M_PI_2f; // 90 degrees	

	int i = (int)(theta_i * nr_bins / full_theta_i) % nr_bins;
	int j = (int)(theta_r * nr_bins / full_theta_r) % nr_bins;

	return matrix[i][j];

}

float* find_ratio_bound(float ** marginal_densities, int nr_bins) {
	float full_theta_i = M_PI_2f; // 90 degrees always i index!
	float full_theta_r = M_PI_2f; // 90 degrees always j index!

	float* max_stored = new float[nr_bins];

	// looking for upper bound of ration f(theta_i)/g(theta_i)

#pragma omp parallel for
	
	for (int j = 0; j < nr_bins; j++) {
		float max = 1.0f; // apparently closer to one reduces samples (wiki: rejection sampling)
		for (int i = 0; i < nr_bins; i++) {
			float tmp = marginal_densities[i][j] *M_PI_2f;
			max = (tmp > max) ? tmp : max;
		}
		max_stored[j] = max;

	}
	printf("\nArray with ratio upper bounds for 1st rejection sampling initialized");

	return max_stored; 
}

