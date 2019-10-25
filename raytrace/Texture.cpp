// 02562 Rendering Framework
// Written by Jeppe Revall Frisvad, 2011
// Copyright (c) DTU Informatics 2011

#include <iostream>
#include <optix_world.h>
#include "my_glut.h"
#include "../SOIL/SOIL.h"
#include "Texture.h"

using namespace std;
using namespace optix;

void Texture::load(const char* filename)
{
  SOIL_free_image_data(data);
  data = SOIL_load_image(filename, &width, &height, &channels, SOIL_LOAD_AUTO);
  if(!data)
  {
    cerr << "Error: Could not load texture image file." << endl;
    return;
  }
  int img_size = width*height;
  delete [] fdata;
  fdata = new float4[img_size];
  for(int i = 0; i < img_size; ++i)
    fdata[i] = look_up(i);
  tex_handle = SOIL_create_OGL_texture(data, width, height, channels, tex_handle, SOIL_FLAG_INVERT_Y | SOIL_FLAG_TEXTURE_REPEATS);
  tex_target = GL_TEXTURE_2D;
}

void Texture::load(GLenum target, GLuint texture)
{
  glBindTexture(target, texture);    
  glGetTexLevelParameteriv(target, 0, GL_TEXTURE_WIDTH, &width);
  glGetTexLevelParameteriv(target, 0, GL_TEXTURE_HEIGHT, &height);
  delete [] fdata;
  fdata = new float4[width*height];
  glGetTexImage(target, 0, GL_RGBA, GL_FLOAT, &fdata[0].x);
  tex_handle = texture;
  tex_target = target;
}

float4 Texture::sample_nearest(const float3& texcoord) const
{
  if(!fdata)
    return make_float4(0.0f);

  float u = texcoord.x;
  float v = texcoord.y;

  float s = u - floor(u);
  float t = v - floor(v);

  int a = int(s * width);
  int b = height- int(t*height) -1; // my original one
  int idx = (int)b * width + a;


  //float a = s * width;
  //float b = t * height; // slides
  //int U = (int)(a + 0.5) % width;
  //int V = (int)(b + 0.5) % height;
  //int idx = U + V * width;







  //int a = int(s * width + 0.5);
  //int b = int((1 - t) * height + 0.5) - 1;
  //or use (int) (s*width + 0.5) % width (actually closest)
  // otherwise you try to access index 512 of array 0-511. like this it starts again at 0


  //int idx = (int) b * width + a % (width*height);

  //if (idx >= (width) * (height)) return make_float4(1, 1, 1, 1);
  //if ( idx < 0) return make_float4(0, 0, 0, 0);

  return fdata[idx];


	//float s = texcoord.x - floor(texcoord.x);
	//float t = texcoord.y - floor(texcoord.y);
	//s = s < 0 ? s + 1 : s;
	//t = t < 0 ? t + 1 : t;

	//int a = roundf((1 - s) * height);
	//int b = roundf(t * width) * height;
	//float4 texel_color = fdata[a + b];
	//return texel_color;
  
  // Implement texture look-up of nearest texel.
  //
  // Input:  texcoord      (texture coordinates: u = texcoord.x, v = texcoord.y)
  //
  // Return: texel color found at the given texture coordinates
  //
  // Relevant data fields that are available (see Texture.h)
  // fdata                 (flat array of texture data: texel colors in float4 format)
  // width, height         (texture resolution)
  //
  // Hint: Remember to revert the vertical axis when finding the index
  //       into fdata.

  //return make_float4(0.0f);
}

float4 Texture::sample_linear(const float3& texcoord) const
{
  if(!fdata)
    return make_float4(0.0f);

  //float4 texel_color = fdata[idx];

  float u = texcoord.x;
  float v = texcoord.y;

  float s = u - floor(u);
  float t = v - floor(v);

  int a = int(s * width);
  //int b = int((1 - t) * height);
  int b = height - int(t * height) - 1;

  float4 t_00 = fdata[(int)b * width + a];
  float4 t_01 = fdata[(int)b * width + (a + 1) % width];
  float4 t_10 = fdata[(int)(b + 1) % height * width + a];
  float4 t_11 = fdata[(int)(b + 1) % height * width + (a + 1) % width];

  float c1 = s * width - floor(s * width);
  float c2 = t * height - floor(t * height);
   
  float4 texel_color = bilerp(t_00, t_01, t_10, t_11, c1, c2);

  //return sample_nearest(texcoord);
  return texel_color;


  // Implement texture look-up which returns the bilinear interpolation of
  // the four nearest texel neighbors.
  //
  // Input:  texcoord      (texture coordinates: u = texcoord.x, v = texcoord.y)
  //
  // Return: texel color found at the given texture coordinates after
  //         bilinear interpolation
  //
  // Relevant data fields that are available (see Texture.h)
  // fdata                 (flat array of texture data: texel colors in float4 format)
  // width, height         (texture resolution)
  //
  // Hint: Use three lerp operations (or one bilerp) to perform the
  //       bilinear interpolation.

}

float4 Texture::look_up(unsigned int idx) const
{
  idx *= channels;
  switch(channels)
  {
  case 1: 
  {
    float v = convert(data[idx]);
    return make_float4(v, v, v, 1.0f);
  }
  case 2: 
    return make_float4(convert(data[idx]), convert(data[idx]), convert(data[idx]), convert(data[idx + 1]));
  case 3: 
    return make_float4(convert(data[idx]), convert(data[idx + 1]), convert(data[idx + 2]), 1.0f);
  case 4: 
    return make_float4(convert(data[idx]), convert(data[idx + 1]), convert(data[idx + 2]), convert(data[idx + 3]));
  }
  return make_float4(0.0f);
}

float Texture::convert(unsigned char c) const
{
  return (c + 0.5f)/256.0f;
}
