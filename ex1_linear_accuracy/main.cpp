/**
 * @file main.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-04-16
 */

#include <iostream>
#include <string>
#include <sstream>
#include "DGFEMSpace1D.h"

double I0(const double mu, const double x, const double t) {
  return sin(2*M_PI*x);
}

double BL(const double mu, const double x, const double t) {
  return sin(2*M_PI*(-1*t));
}

double BR(const double mu, const double x, const double t) {
  return sin(2*M_PI*(1-1*t));
}

double sigma_t(const double x) {
  return 0;
}

double sigma_s(const double x) {
  return 0;
}

double q(const double mu, const double x, const double t) {
  return 0;
}

int main(int argc, char *argv[]) {
  if(argc != 4) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
    abort();
  }

  clock_t t1, t2;
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  std::cout << "Set up problem ..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info ..." << std::endl;
  Problem.BuildQuad(K);
  std::cout << "Build universal matrix ..." << std::endl;
  Problem.BuildUniMAT();
  std::cout << "Initialize ..." << std::endl;
  Problem.init(I0);
  std::cout << "Start to solve ..." << std::endl;
  t1 = clock();
  Problem.run_unsteady(sigma_t, sigma_s, q, BL, BR, 3e-3);
  std::cout << "Finished ..." << std::endl;
  t2 = clock();
  std::stringstream s;
  s << "ex1_Nx" << Nx << "_K" << K << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  Problem.print_solution_integral(out);
  out.close();
  std::cout << "Time consumed: "
    << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;

  return 0;
}

