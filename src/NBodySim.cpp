/*
C++ implementation of gravitational N-Body Problem
- https://medium.com/swlh/create-your-own-n-body-simulation-with-matlab-22344954228e
- we are creating a simulation of dynamical system of particles interacting gravitationally
*/

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <time.h>
#include <vector>

#define UNUSED(x) (void)(x)

// Force Calculation
// - Each particle reacts to other particles graviational attraction -> Newton's law of universal gravitation
std::vector<Eigen::Vector3f> getAcc(const Eigen::MatrixX3f& pos,
                                    const Eigen::VectorXf& mass,
                                    float G,
                                    float softening)
{
  UNUSED(mass);
  UNUSED(G);
  UNUSED(softening);

  // auto N = pos.size();
  /*
  Eigen::MatrixXf zeros = Eigen::MatrixXf::Zero(pos.rows(), 1);

  std::cout << "xStuff:\n"
            << pos(Eigen::all, 0) << "\nyStuff:\n"
            << pos(Eigen::all, 1) << "\nzStuff:\n"
            << pos(Eigen::all, 2) << std::endl;

  Eigen::MatrixX3f x(pos.rows(), 3);
  x << pos(Eigen::all, 0), zeros, zeros;
  Eigen::MatrixX3f y(pos.rows(), 3);
  y << zeros, pos(Eigen::all, 1), zeros;
  Eigen::MatrixX3f z(pos.rows(), 3);
  z << zeros, zeros, pos(Eigen::all, 2);

  std::cout << "X:\n" << x << "\nY:\n" << y << "\nZ:\n" << z << std::endl;

  std::cout << "HERE" << std::endl;
  std::cout << "x:\n" << x << "\ny:\n" << y << "\nz:\n" << z << std::endl;
  std::cout << "xT:\n"
            << x.transpose() << "\nyT:\n"
            << y.transpose() << "\nzT:\n"
            << z.transpose() << std::endl;

  auto dX = x.transpose() - x;
  auto dY = y.transpose() - y;
  auto dZ = z.transpose() - z;
  */

  Eigen::MatrixXf x = pos.col(0);
  Eigen::MatrixXf y = pos.col(1);
  Eigen::MatrixXf z = pos.col(2);

  std::cout << "X:\n" << x << "\nY:\n" << y << "\nZ:\n" << z << std::endl;

  Eigen::MatrixXf xT = x.transpose();
  Eigen::MatrixXf yT = y.transpose();
  Eigen::MatrixXf zT = z.transpose();

  std::cout << "xT:\n" << xT << "\nyT:\n" << yT << "\nzT:\n" << zT << std::endl;

  Eigen::MatrixXf xTR = xT.replicate(pos.rows(), 1);
  Eigen::MatrixXf yTR = yT.replicate(pos.rows(), 1);
  Eigen::MatrixXf zTR = zT.replicate(pos.rows(), 1);

  std::cout << "xTR:\n" << xTR << "\nyTR:\n" << yTR << "\nzTR:\n" << zTR << std::endl;

  Eigen::MatrixXf dX = xTR - x.colwise();
  Eigen::MatrixXf dY = yTR - y.colwise();
  Eigen::MatrixXf dZ = zTR - z.colwise();

  std::cout << "dX:\n" << dX << "\ndY:\n" << dY << "\ndZ:\n" << dZ << std::endl;

  /*
  auto dX = x.transpose().replicate(pos.rows(), 1) - x.colwise();
  auto dY = y.transpose().replicate(pos.rows(), 1) - y.colwise();
  auto dZ = z.transpose().replicate(pos.rows(), 1) - z.colwise();
  */

  // std::cout << "dX:\n" << dX << "\ndY:\n" << dY << "\ndZ:\n" << dZ << std::endl;

  /*

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Eigen::Vector3f dVec = pos[j] - pos[i];
      auto inv_r3 = powf((exp2f(dVec(0)) + exp2f(dVec(1)) + exp2f(dVec(2)) + exp2f(softening)), -1.5f);
      a[i] += G * (dVec * inv_r3) * mass[j];
    }
  }
  */
  return {};
}

// Energy Calculation

int main(int /*argc*/, char** /*argv*/)
{
  srand((unsigned int)time(NULL));
  std::cout << "Hello World" << std::endl;

  // Simulation Parameters
  constexpr int N = 10;
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
  UNUSED(N);
  UNUSED(softening);
  UNUSED(G);

  // total mass of particles is 20?
  constexpr const float initialMass = 20.0f * 1.0f / float(N);
  Eigen::VectorXf mass = initialMass * Eigen::VectorXf::Ones(N);

  // fills pos with random numbers between -1 and 1
  Eigen::MatrixX3f pos = Eigen::MatrixX3f::Random(N, 3);
  Eigen::MatrixX3f vel = Eigen::MatrixX3f::Random(N, 3);

  // clang-format off
  /*
  std::cout << "Mass:\n" << mass
            << "\nVel:\n" << vel
            << "\nVel * mass:\n" << vel.cwiseProduct(mass.replicate(1, 3))
            << "\nMean(Vel * mass):\n" << vel.cwiseProduct(mass.replicate(1, 3)).colwise().mean()
            << "\nMean(mass):\n" << mass.mean()
            << "\nMean(Vel * mass) / Mean(mass):\n" << vel.cwiseProduct(mass.replicate(1, 3)).colwise().mean() / mass.mean()
            << "new_vel:\n" << vel - (vel.cwiseProduct(mass.replicate(1, 3)).colwise().mean() / mass.mean()).replicate(N, 1)
            << std::endl;
  */
  // clang-format on

  vel = vel - (vel.cwiseProduct(mass.replicate(1, 3)).colwise().mean() / mass.mean()).replicate(N, 1);

  // Computing initial gravitational accelerations
  auto acc = getAcc(pos, mass, G, softening);

  /*
  // convert to center-of-mass frame - the velocity of the center of mass (duh)
  // vel -= np.mean(mass, vel, 0) / np.mean(mass)
  // this is computing the mean of mass * vel
  Eigen::Vector3f meanOfMTV = MeanOfVecVec3(VecF32TimeVecVec3(mass, vel));
  float meanOfMass = MeanOfVecf32(mass);
  auto quotientOfMeanMVM = meanOfMTV / meanOfMass;
  for (size_t i = 0; i < vel.size(); i++) {
    vel[i] = vel[i] - quotientOfMeanMVM;
  }


*/
  return 0;
}
