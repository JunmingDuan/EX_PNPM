#include "DGFEMSpace1D.h"

void DGFEMSpace1D::scaling_limiter(const u_int cell, EVEC& u) {
  int flag(0);
  double min(1);
  for(u_int g = 0; g < K; ++g) {
    min = std::min(min, u[g]);
  }
  if(min < 0) flag = 1;
  VEC<double> w = QUADINFO[cell].weight();
  double av(0);
  if(flag == 1) {
    for(u_int g = 0; g < K; ++g) {
      av += u[g]*w[g];
    }
    av /= 2.0;
    if(av < 0) {
      std::cout << "av < 0 !!!" << std::endl;
      abort();
    }
    else {
      double k = fabs(av / (av - min));
      for(u_int g = 0; g < K; ++g) {
        u[g] = k*u[g] + (1-k)*av;
      }
      //abort();
    }
  }
}

