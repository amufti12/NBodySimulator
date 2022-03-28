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

constexpr const int N = 1000;

// Force Calculation
// - Each particle reacts to other particles graviational attraction -> Newton's law of universal gravitation
Eigen::MatrixXf getAcc(const Eigen::MatrixX3f& pos, const Eigen::VectorXf& mass, float G, float softening)
{
  Eigen::MatrixXf x = pos.col(0);
  Eigen::MatrixXf y = pos.col(1);
  Eigen::MatrixXf z = pos.col(2);

  Eigen::MatrixXf xT = x.transpose();
  Eigen::MatrixXf yT = y.transpose();
  Eigen::MatrixXf zT = z.transpose();

  Eigen::MatrixXf xTR = xT.replicate(x.rows(), 1);
  Eigen::MatrixXf yTR = yT.replicate(y.rows(), 1);
  Eigen::MatrixXf zTR = zT.replicate(z.rows(), 1);

  Eigen::MatrixXf dX = xTR - x.replicate(1, xTR.cols());
  Eigen::MatrixXf dY = yTR - y.replicate(1, yTR.cols());
  Eigen::MatrixXf dZ = zTR - z.replicate(1, zTR.cols());

  Eigen::MatrixXf dX2 = dX.array().pow(2);
  Eigen::MatrixXf dY2 = dY.array().pow(2);
  Eigen::MatrixXf dZ2 = dZ.array().pow(2);
  Eigen::MatrixXf softening2 = (Eigen::MatrixXf::Ones(dX.rows(), dX.cols()) * softening).array().pow(2);

  Eigen::MatrixXf inv_r3 = (dX2 + dY2 + dZ2 + softening2);

  inv_r3 = (inv_r3.array() > 0.0f).select(inv_r3.array().pow(-1.5), inv_r3);

  Eigen::MatrixXf aX = (G * (dX * inv_r3)) * mass;
  Eigen::MatrixXf aY = (G * (dY * inv_r3)) * mass;
  Eigen::MatrixXf aZ = (G * (dZ * inv_r3)) * mass;

  Eigen::MatrixXf resultVector(aX.rows(), aX.cols() + aY.cols() + aZ.cols());
  resultVector << aX, aY, aZ;

  return resultVector;
}

// Energy Calculation

int main(int /*argc*/, char** /*argv*/)
{
  srand((unsigned int)time(NULL));
  std::cout << "Hello World" << std::endl;

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
