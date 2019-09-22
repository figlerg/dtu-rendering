// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "HitInfo.h"
#include "Triangle.h"

using namespace optix;

bool intersect_triangle(const Ray& ray, 
                        const float3& v0, 
                        const float3& v1, 
                        const float3& v2, 
                        float3& n,
                        float& t,
                        float& v,
                        float& w)
{
  // Implement ray-triangle intersection here (see Listing 1 in the lecture note).
  // Note that OptiX also has an implementation, so you can get away
  // with not implementing this function. However, I recommend that
  // you implement it for completeness.
	//if(abs(dot(ray.direction, n) <= 1e-07)) 
	//	return false;

	//float t_prime = dot(v0 - ray.origin, n) / dot(ray.direction, n);
	//float3 e0 = v1 - v2;
	//float3 e1 = v0 - v2;

	//// TODO ask why we already have v, w?

	//float v_new = dot(cross((v0 - ray.origin), ray.direction), make_float3(1, 0, 0)) / dot(ray.direction, n);	
	//float w_new = dot(cross((v0 - ray.origin), ray.direction), make_float3(1, 0, 0)) / dot(ray.direction, n);

	//float u = 1 - v - w;

	//if (v >= 0 && w >= 0 && u <= 0)
	//	return true;

  return false;
}


bool Triangle::intersect(const Ray& r, HitInfo& hit, unsigned int prim_idx) const
{
  // Implement ray-triangle intersection here.
  //
  // Input:  r                    (the ray to be checked for intersection)
  //         prim_idx             (index of the primitive element in a collection, not used here)
  //
  // Output: hit.has_hit          (set true if the ray intersects the triangle)
  //         hit.dist             (distance from the ray origin to the intersection point)
  //         hit.position         (coordinates of the intersection point)
  //         hit.geometric_normal (the normalized normal of the triangle)
  //         hit.shading_normal   (the normalized normal of the triangle)
  //         hit.material         (pointer to the material of the triangle)
  //        (hit.texcoord)        (texture coordinates of intersection point, not needed for Week 1)
  //
  // Return: True if the ray intersects the triangle, false otherwise
  //
  // Relevant data fields that are available (see Triangle.h)
  // r                            (the ray)
  // v0, v1, v2                   (triangle vertices)
  // (t0, t1, t2)                 (texture coordinates for each vertex, not needed for Week 1)
  // material                     (material of the triangle)
  //
  // Hint: Use the function intersect_triangle(...) to get the hit info.
  //       Note that you need to do scope resolution (optix:: or just :: in front
  //       of the function name) to choose between the OptiX implementation and
  //       the function just above this one.
	
	float3 n;
	float t, v, w;

	bool is_hit = optix::intersect_triangle(r, v0 , v1, v2, n, t, v, w);

	if (is_hit && (t >= r.tmin & t <= r.tmax)) {
		hit.dist = t;
		hit.position = r.origin + r.direction * t;
		hit.geometric_normal = n;
		hit.shading_normal = n;
		hit.material = &material;
		hit.has_hit = true; // otherwise it overwrites hitinfo wrongly

		// TODO looks wrong!

		return true;
		}

	return false;
}

void Triangle::transform(const Matrix4x4& m) 
{ 
  v0 = make_float3(m*make_float4(v0, 1.0f)); 
  v1 = make_float3(m*make_float4(v1, 1.0f)); 
  v2 = make_float3(m*make_float4(v2, 1.0f)); 
}

Aabb Triangle::compute_bbox() const
{
  Aabb bbox;
  bbox.include(v0);
  bbox.include(v1);
  bbox.include(v2);
  return bbox;
}
