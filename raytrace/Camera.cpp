#include <iostream>
#include <optix_world.h>
#include "my_glut.h"
#include "Camera.h"

using namespace optix;

void Camera::set(const float3& eye_point, const float3& view_point, const float3& up_vector, float camera_constant)
{
  eye = eye_point; // e
  lookat = view_point; // p
  up = up_vector;
  cam_const = camera_constant;

  ip_normal = normalize(lookat - eye);
  ip_xaxis = normalize(cross(ip_normal, up));
  ip_yaxis = normalize(cross(ip_xaxis, ip_normal));
   

  // Compute camera coordinate frame (image plane normal and basis).
  //
  // Relevant data fields that are available (see Camera.h)
  // ip_normal  (the viewing direction and the 3rd basis vector, v)
  // ip_xaxis   (the 1st basis vector, b_1)
  // ip_yaxis   (the 2nd basis vector, b_2)
  //
  // Hint: The OptiX math library has a function normalize(v) which returns
  //       the vector v normalized and a function cross(v, w) which returns
  //       the cross product of the vectors v and w.

  fov = 2 * atan(0.5 / cam_const)* 180 * M_1_PIf;
  // vertical field of view with basic tangens calc -> to degrees

  // Assume that the height and the width of the film is 1.
  // With this information, use the pinhole camera model to compute
  // the field of view (fov) in degrees.
  //
  // Relevant data fields that are available (see Camera.h)
  // camera_constant  (the camera constant, d)
  // fov              (the vertical field of view of the camera, phi)
  //
  // Hints: (a) The inverse tangent function is called atan(x).
  //        (b) The OptiX math library has a constant for 1/pi called M_1_PIf.
  //        (c) The following constant is the field of view that you should 
  //            get with the default scene (if you set a breakpoint and run
  //            in Debug mode, you can use it to check your code).
  
}

/// Get direction of viewing ray from image coords.
float3 Camera::get_ray_dir(const float2& coords) const
{
	float3 q_ip = make_float3(coords.x, coords.y, cam_const);
	
	Matrix3x3 transformation_matrix;
	transformation_matrix.setCol(0, ip_xaxis);
	transformation_matrix.setCol(1, ip_yaxis);
	transformation_matrix.setCol(2, ip_normal);

	float3 dir_vec = normalize(transformation_matrix * q_ip);
	   	 


  // Given the image plane coordinates, compute the normalized ray
  // direction by a change of basis.
  // return make_float3(0.0f);
	return dir_vec;
}

/// Return the ray corresponding to a set of image coords
Ray Camera::get_ray(const float2& coords) const
{
	float x = coords.x;
	float y = coords.y;
	float3 dir_vector = get_ray_dir(coords);

	Ray output_ray = Ray(eye, dir_vector, 0, 0);

  // Use the function get_ray_dir(...) to construct a ray starting at
  // the eye and going through the position in the image plane given
  // by the image plane coordinates received as argument.
  //
  // Hint: You can set the ray type to 0 as the framework currently
  //       does not use different ray types.
  return output_ray; 
}

// OpenGL

void Camera::glSetPerspective(unsigned int width, unsigned int height) const
{
  GLdouble aspect = width/static_cast<float>(height);    

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fov, aspect, cam_const*NEAR_PLANE, FAR_PLANE);

  glMatrixMode(GL_MODELVIEW);
}

void Camera::glSetCamera() const
{
  gluLookAt(eye.x,   eye.y,   eye.z, 
	          lookat.x, lookat.y, lookat.z, 
	          up.x,    up.y,    up.z);
}

