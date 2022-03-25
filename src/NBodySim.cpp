/*
C++ implementation of gravitational N-Body Problem
- https://medium.com/swlh/create-your-own-n-body-simulation-with-matlab-22344954228e
- we are creating a simulation of dynamical system of particles interacting gravitationally
*/

// include Eigen --> future
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <time.h>
#include <vector>

#define UNUSED(x) (void)(x)

struct vec3
{
  float x = 0.0f;
  float y = 0.0f;
  float z = 0.0f;

  static vec3 rand()
  {
    vec3 result = { .x = float(::rand()) / RAND_MAX * 2.0f - 1.0f,
                    .y = float(::rand()) / RAND_MAX * 2.0f - 1.0f,
                    .z = float(::rand()) / RAND_MAX * 2.0f - 1.0f };
    return result;
  }
};

vec3 operator*(float f, const vec3& v)
{
  vec3 result = {};
  result.x = f * v.x;
  result.y = f * v.y;
  result.z = f * v.z;

  return result;
}

vec3 operator*(const vec3& v, float f)
{
  return f * v;
}

vec3 operator/(const vec3& v, float f)
{
  return (1.0f / f) * v;
}

vec3 operator-(const vec3& v1, const vec3& v2)
{
  vec3 result = {};
  result.x = v1.x - v2.x;
  result.y = v1.y - v2.y;
  result.z = v1.z - v2.z;
  return result;
}

// computes the mean of a vector of floats
float MeanOfVecf32(const std::vector<float>& v)
{
  float result = std::accumulate(v.begin(), v.end(), 0.0f);
  return result / float(v.size());
}

// computing the product of a vector of floats and vector of vec3
std::vector<vec3> VecF32TimeVecVec3(const std::vector<float>& f, const std::vector<vec3>& v)
{
  std::vector<vec3> resultVec(f.size());
  for (size_t i = 0; i < f.size(); i++) {
    resultVec[i].x = f[i] * v[i].x;
    resultVec[i].y = f[i] * v[i].y;
    resultVec[i].z = f[i] * v[i].z;
  }
  return resultVec;
}

// computes the mean of a vector of vec3 and returns a vec3
vec3 MeanOfVecVec3(const std::vector<vec3>& v)
{
  vec3 resultVec3 = {};
  for (const auto& vec : v) {
    resultVec3.x += vec.x;
    resultVec3.y += vec.y;
    resultVec3.z += vec.z;
  }
  resultVec3.x /= float(v.size());
  resultVec3.y /= float(v.size());
  resultVec3.z /= float(v.size());

  return resultVec3;
}

// Force Calculation
// - Each particle reacts to other particles graviational attraction -> Newton's law of universal gravitation
std::vector<vec3> getAcc(std::vector<vec3>& pos, std::vector<float> mass, float G, float softening)
{
  auto N = pos.size();
  auto a = std::vector<vec3>(N);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      auto dx = pos[j].x - pos[i].x;
      auto dy = pos[j].y - pos[i].y;
      auto dz = pos[j].z - pos[i].z;
      auto inv_r3p1 = (exp2(dx) + exp2(dy) + exp2(dz) + exp2(softening));
      auto inv_r3 = pow(inv_r3p1, -1.5f);
      a[i].x += G * (dx * inv_r3) * mass[j];
      a[i].y += G * (dx * inv_r3) * mass[j];
      a[i].z += G * (dx * inv_r3) * mass[j];
    }
  }
  return a;
}

// Energy Calculation

int main(int /*argc*/, char** /*argv*/)
{
  srand((unsigned int)time(NULL));
  std::cout << "Hello World" << std::endl;

  // Simulation Parameters
  int N = 100;
  // Number of particles
  int t = 0;                 // Current Time of the Simulation
  float tEnd = 10.0f;        // Time at which simulation ends
  float dt = 0.01f;          // Timestep
  float softening = 0.1f;    // Softening Length
  float G = 1.0f;            // Newton's Gravitational Constant
  bool plotRealTime = false; // I will not be using this

  // unused section
  UNUSED(t);
  UNUSED(tEnd);
  UNUSED(dt);
  UNUSED(plotRealTime);

  // total mass of particles is 20?
  auto mass = std::vector<float>(N);
  // fills pos with random numbers between -1 and 1
  auto pos = std::vector<vec3>(N);
  auto vel = std::vector<vec3>(N);
  for (int i = 0; i < N; i++) {
    mass[i] = 20.0f * 1.0f / float(N);
    pos[i] = vec3::rand();
    vel[i] = vec3::rand();
  }

  // convert to center-of-mass frame - the velocity of the center of mass (duh)
  // vel -= np.mean(mass, vel, 0) / np.mean(mass)
  // this is computing the mean of mass * vel
  vec3 meanOfMTV = MeanOfVecVec3(VecF32TimeVecVec3(mass, vel));
  float meanOfMass = MeanOfVecf32(mass);
  auto quotientOfMeanMVM = meanOfMTV / meanOfMass;
  for (size_t i = 0; i < vel.size(); i++) {
    vel[i] = vel[i] - quotientOfMeanMVM;
  }

  auto acc = getAcc(pos, mass, G, softening);

  return 0;
}
