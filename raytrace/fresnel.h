// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#ifndef FRESNEL_H
#define FRESNEL_H

inline float fresnel_r_s(float cos_theta1, float cos_theta2, float ior1, float ior2)
{
  // Compute the perpendicularly polarized component of the Fresnel reflectance

	float ni_cosi = ior1 * cos_theta1;
	float nt_cost = ior2 * cos_theta2;


  return (ni_cosi-nt_cost)/(ni_cosi + nt_cost);
}

inline float fresnel_r_p(float cos_theta1, float cos_theta2, float ior1, float ior2)
{
  // Compute the parallelly polarized component of the Fresnel reflectance
	float nt_cosi = ior2 * cos_theta1;
	float ni_cost = ior1 * cos_theta2;


	return (nt_cosi - ni_cost) / (nt_cosi + ni_cost);
	
	//return 0.0;
}

inline float fresnel_R(float cos_theta1, float cos_theta2, float ior1, float ior2)
{
  // Compute the Fresnel reflectance using fresnel_r_s(...) and fresnel_r_p(...)
	float r1 = fresnel_r_s(cos_theta1, cos_theta2, ior1, ior2);
	float r2 = fresnel_r_p(cos_theta1, cos_theta2, ior1, ior2);


	return 0.5*(powf(r1,2.0f) + powf(r2, 2.0f));
}

#endif // FRESNEL_H
