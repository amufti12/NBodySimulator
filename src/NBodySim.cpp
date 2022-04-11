/*
C++ implementation of gravitational N-Body Problem
- https://medium.com/swlh/create-your-own-n-body-simulation-with-matlab-22344954228e
- we are creating a simulation of dynamical system of particles interacting gravitationally
*/

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

#include <imgui.h>

#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include <implot.h>

#include <glad/glad.h>

#include <GLFW/glfw3.h>

#define UNUSED(x) (void)(x)

constexpr const int N = 256;
uint32_t width = 1280;
uint32_t height = 720;

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

  Eigen::MatrixXf aX = (G * (dX.cwiseProduct(inv_r3))) * mass;

  Eigen::MatrixXf aY = (G * (dY.cwiseProduct(inv_r3))) * mass;
  Eigen::MatrixXf aZ = (G * (dZ.cwiseProduct(inv_r3))) * mass;

  Eigen::MatrixXf resultVector(aX.rows(), aX.cols() + aY.cols() + aZ.cols());
  resultVector << aX, aY, aZ;

  return resultVector;
}

// Energy Calculation
std::tuple<float, float> getEnergy(const Eigen::MatrixX3f& pos,
                                   const Eigen::MatrixX3f& vel,
                                   const Eigen::VectorXf& mass,
                                   float G)
{
  // Kinetic Energy computation
  Eigen::MatrixXf velsq2 = vel.array().pow(2);
  auto KE = 0.5f * ((velsq2.cwiseProduct(mass.replicate(1, 3))).sum());

  // initializing the r positions of particles
  Eigen::MatrixXf x = pos.col(0);
  Eigen::MatrixXf y = pos.col(1);
  Eigen::MatrixXf z = pos.col(2);

  // creating  the matrix that stores all pairwise particle separations: r_j - r_i
  // transposing each position
  Eigen::MatrixXf xT = x.transpose();
  Eigen::MatrixXf yT = y.transpose();
  Eigen::MatrixXf zT = z.transpose();

  // resizing to perform subtraction computation
  Eigen::MatrixXf xTR = xT.replicate(x.rows(), 1);
  Eigen::MatrixXf yTR = yT.replicate(y.rows(), 1);
  Eigen::MatrixXf zTR = zT.replicate(z.rows(), 1);

  // performing r_j - r_i
  Eigen::MatrixXf dX = xTR - x.replicate(1, xTR.cols());
  Eigen::MatrixXf dY = yTR - y.replicate(1, yTR.cols());
  Eigen::MatrixXf dZ = zTR - z.replicate(1, zTR.cols());

  // initializing matrix for 1/r for particle pairwise particle separation
  Eigen::MatrixXf dX2 = dX.array().pow(2);
  Eigen::MatrixXf dY2 = dY.array().pow(2);
  Eigen::MatrixXf dZ2 = dZ.array().pow(2);

  // performing 1/r
  Eigen::MatrixXf inv_r = (dX2 + dY2 + dZ2).cwiseSqrt();
  inv_r = (inv_r.array() > 0.0f).select(1 / inv_r.array(), inv_r);

  // PE computation -> "sum over upper triangle to count each interaction only once"
  Eigen::MatrixXf massT = mass.transpose();
  Eigen::MatrixXf MTMt = -(mass * mass.transpose());
  Eigen::MatrixXf MTmtTinv = MTMt.cwiseProduct(inv_r);

  // loop to create an upper triangle starting at 1
  int startingpoint = 1;
  for (int i = 0; i < MTmtTinv.rows(); i++) {
    for (int j = 0; j < MTmtTinv.cols(); j++) {
      if (j < startingpoint) {
        MTmtTinv(i, j) = 0.0f;
      }
    }
    startingpoint++;
  }

  auto PE = G * MTmtTinv.sum();
  return std::make_tuple(KE, PE);
}

void recordPosAtTime(Eigen::MatrixXf& pos_save, const Eigen::MatrixXf& pos, int timeIndex)
{
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp(pos);
  Eigen::Map<Eigen::VectorXf> tempMap(temp.data(), temp.size());
  pos_save.block<N * 3, 1>(0, timeIndex) = tempMap;
}

void drawScatter(Eigen::MatrixX3f& pos)
{
  // separate x's and y's into arrays of floats
  float* xdata = pos.col(0).data();
  float* ydata = pos.col(1).data();

  float aspect_ratio = (float)height / (float)width;

  // draw
  if (ImPlot::BeginPlot("Scatter Plot",
                        ImVec2((float)width - 20, (float)height - 40),
                        ImPlotFlags_CanvasOnly | ImPlotFlags_NoInputs)) {
    ImPlot::SetupAxisLimits(ImAxis_X1, -5, 5);
    ImPlot::SetupAxisLimits(ImAxis_Y1, -5 * aspect_ratio, 5 * aspect_ratio);
    ImPlot::PlotScatter("Positions", xdata, ydata, N);
    ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
    ImPlot::SetNextMarkerStyle(
      ImPlotMarker_Circle, 6, ImPlot::GetColormapColor(1), IMPLOT_AUTO, ImPlot::GetColormapColor(1));
    ImPlot::PopStyleVar();
    ImPlot::EndPlot();
  }
}

int main(int /*argc*/, char** /*argv*/)
{
  srand((unsigned int)time(NULL));

  // Number of particles
  float t = 0.0f;            // Current Time of the Simulation
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

  ///////////////////////////////////////////
  // Creating GLFW Window Section
  if (glfwInit() != GLFW_TRUE) {
    std::exit(-1);
  }

  // setting window hints
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, false);
  glfwWindowHint(GLFW_SAMPLES, 4); // Anti-aliasing

  // glfwcreateWindow(width, height, title, monitor, share)
  GLFWwindow* window = glfwCreateWindow(width, height, "NBodySimulator", nullptr, nullptr);
  glfwMakeContextCurrent(window);
  // Disable vsync
  glfwSwapInterval(0);

  // Loading all function pointers for OpenGL
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::exit(-1);
  }

  glEnable(GL_MULTISAMPLE); // Turning anti-aliasing on

  ///////////////////////////////////////////
  //// IMGUI SECTION
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  io.ConfigFlags |= ImGuiConfigFlags_NoMouse;
  io.ConfigFlags |= ImGuiConfigFlags_NavNoCaptureKeyboard;
  UNUSED(io);
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 410");

  //////////////////////////////////////////

  // total mass of particles is 20?
  constexpr const float initialMass = 10.0f * 1.0f / float(N);
  Eigen::VectorXf mass = initialMass * Eigen::VectorXf::Ones(N);

  // fills pos with random numbers between -1 and 1
  Eigen::MatrixX3f pos = Eigen::MatrixX3f::Random(N, 3);
  Eigen::MatrixX3f vel = Eigen::MatrixX3f::Random(N, 3);

  // center of mass frame conversion
  vel = vel - (vel.cwiseProduct(mass.replicate(1, 3)).colwise().mean() / mass.mean()).replicate(N, 1);

  // Computing initial gravitational accelerations
  auto acc = getAcc(pos, mass, G, softening);
  auto [KE, PE] = getEnergy(pos, vel, mass, G);

  // number of timesteps
  int Nt = int((tEnd / dt) + 0.5f);

  // save energies for trails
  Eigen::MatrixXf pos_save = Eigen::MatrixXf::Zero(N * 3, Nt + 1);
  std::vector<std::pair<float, float>> KEPE_save(Nt + 1);
  KEPE_save[0] = std::make_pair(KE, PE);

  // set initial positions in pos_save
  recordPosAtTime(pos_save, pos, 0);

  // simulation loop
  for (int i = 0; i < Nt; i++) {
    glfwPollEvents();

    glClear(GL_COLOR_BUFFER_BIT);

    // Displaying the UI Window
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();

    ImGui::NewFrame();
    {
      ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(),
                                   ImGuiDockNodeFlags_NoSplit | ImGuiDockNodeFlags_NoResize);

      ImGui::Begin("NBody Simulation Visualization");
      drawScatter(pos);
      ImGui::End();
    }
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);

    // 1/2 kick (idk what this means)
    vel += acc * dt / 2.0;

    // drift
    pos += vel * dt;

    // update accelerations
    acc = getAcc(pos, mass, G, softening);

    // 1/2 kick again not sure what this means
    vel += acc * dt / 2.0;

    // update time
    t += dt;

    // get energy of system
    auto [sKE, sPE] = getEnergy(pos, vel, mass, G);
    UNUSED(sKE);
    UNUSED(sPE);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();

  glfwTerminate();

  return 0;
}
