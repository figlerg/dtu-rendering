// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <optix_world.h>
#include "IndexedFaceSet.h"
#include "ObjMaterial.h"
#include "mt_random.h"
#include "cdf_bsearch.h"
#include "HitInfo.h"
#include "AreaLight.h"
#include "int_pow.h"


using namespace optix;

bool AreaLight::sample(const float3& pos, float3& dir, float3& L) const

{
  const IndexedFaceSet& normals = mesh->normals;
  L = make_float3(0.0f);


  //float3 light_pos = mesh->compute_bbox().center(); // original line
  // if I understand correctly, here the sampling is not random but from the center?
  // instead, choose a random triangle and then a random point on it

  // choose one of the triangles in the mesh:
  int n = mesh->geometry.no_faces();
  float xi = mt_random()*n; // uniform distribution of [0,n]
  int chosen = (int)xi % n; // this is now the index of the randomly chosen triangle
   
  // sample a random point on the triangle using barycentric coordinates
  float xi1 = mt_random();
  float xi2 = mt_random();

  // these are the barycentric coordinates
  float u = 1 - sqrtf(xi1);
  float v = (1 - xi2) * sqrtf(xi1);
  float w = xi2 * sqrtf(xi1);

  // use them to get actual point in triangle
  const uint3& geo_face = mesh->geometry.face(chosen); // like in trimesh.cpp
  float3 light_pos = u * mesh->geometry.vertex(geo_face.x)+ v * mesh->geometry.vertex(geo_face.y)+ w * mesh->geometry.vertex(geo_face.z);


  //const uint3& norm_face = mesh->normals.face(chosen);
  //float

  //float3 light_pos = normalize()


  float dist = length(light_pos - pos);

  dir = normalize(light_pos - pos);
  float cutoff = length(pos - light_pos) - 0.0001;

  Ray shadow_ray = Ray(pos, dir, 0, 0.001, cutoff);
  HitInfo info = HitInfo();

  bool shadowed = tracer->trace_to_any(shadow_ray, info);


  float distance = optix::length(light_pos - pos);
  //L = intensity / (pow(distance, 2));


  for(int i = 0; i < mesh->geometry.no_faces(); i++) {
	  float3 emission = get_emission(i);
	  uint3 face = normals.face(i);
	  float3 normal = normalize(normals.vertex(face.x) + normals.vertex(face.y) + normals.vertex(face.z));
		
	  L += dot(-dir, normal) * emission * mesh->face_areas[i];
  }
  L /= int_pow(dist, 2);


  return !shadowed;



  // Compute output and return value given the following information.
  //
  // Input:  pos  (the position of the geometry in the scene)
  //
  // Output: dir  (the direction toward the light)
  //         L    (the radiance received from the direction dir)
  //
  // Return: true if not in shadow
  //
  // Relevant data fields that are available (see Light.h and above):
  // shadows             (on/off flag for shadows)
  // tracer              (pointer to ray tracer)
  // normals             (indexed face set of vertex normals)
  // mesh->face_areas    (array of face areas in the light source)
  //
  // Hints: (a) Use mesh->compute_bbox().center() to get the center of 
  //        the light source bounding box.
  //        (b) Use the function get_emission(...) to get the radiance
  //        emitted by a triangle in the mesh.

  return false;  
}

bool AreaLight::emit(Ray& r, HitInfo& hit, float3& Phi) const
{
  // Generate and trace a ray carrying radiance emitted from this area light.
  //
  // Output: r    (the new ray)
  //         hit  (info about the ray-surface intersection)
  //         Phi  (the flux carried by the emitted ray)
  //
  // Return: true if the ray hits anything when traced
  //
  // Relevant data fields that are available (see Light.h and Ray.h):
  // tracer              (pointer to ray tracer)
  // geometry            (indexed face set of triangle vertices)
  // normals             (indexed face set of vertex normals)
  // no_of_faces         (number of faces in triangle mesh for light source)
  // mesh->surface_area  (total surface area of the light source)
  // r.origin            (starting position of ray)
  // r.direction         (direction of ray)

  // Get geometry info
  const IndexedFaceSet& geometry = mesh->geometry;
	const IndexedFaceSet& normals = mesh->normals;
	const float no_of_faces = static_cast<float>(geometry.no_faces());

  // Sample ray origin and direction
 
  // Trace ray
  
  // If a surface was hit, compute Phi and return true

  return false;
}

float3 AreaLight::get_emission(unsigned int triangle_id) const
{
  const ObjMaterial& mat = mesh->materials[mesh->mat_idx[triangle_id]];
  return make_float3(mat.ambient[0], mat.ambient[1], mat.ambient[2]);
}
