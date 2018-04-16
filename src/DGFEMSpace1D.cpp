/**
 * @file DGFEMSpace1D.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0 for 1D RTE
 * @
 * @date 2018-04-16
 */

#include <stdio.h>
#include <stdlib.h>
#include "DGFEMSpace1D.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr)
  : Nx(Nx), xl(xl), xr(xr) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  mu.resize(M);
  wgt.resize(M);
  I.resize(Nx); I1.resize(Nx); In.resize(Nx);
  W.resize(K*K);
  W_uw.resize(K*K);
  for(u_int i = 0; i < Nx; ++i) {
    I[i].resize(M); I1[i].resize(M); In[i].resize(M);
    for(u_int m = 0; m < M; ++m) {
      I[i][m].resize(K); I1[i][m].resize(K); In[i][m].resize(K);
    }
  }
  M1.resize(K, K);
  Mx.resize(K, K*K);
  Fx_plus.resize(K, K*K);
  Fx_minus.resize(K, K*K);
  K1.resize(K*K, K*K);
  Kx.resize(K*K, K*K);
  F0.resize(K*K, K);
  A.resize(K, K);
  rhs.resize(K);
}

void DGFEMSpace1D::BuildQuad(u_int np) {
  VEC<VEC<double> > pw;
  //Legendre Gauss points
  pw = QUADINFO[0].LG(np);
  TemQuad_Gauss.set_np(np);
  TemQuad_Gauss.set_points(pw[0]);
  TemQuad_Gauss.set_weight(pw[1]);
  QUADINFO.resize(Nx);
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    QUADINFO[i].set_np(np);
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    QUADINFO[i].set_weight(pw[1]);
    VEC<double> x(np);
    for(u_int j = 0; j < np; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[j]);
    }
    QUADINFO[i].set_points(x);
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv),
        global_to_local_jacobian(lv, gv) );
  }
  //direction
  sum_wm = 2;
  if(M > 1) {
    pw = TemQuad_Gauss.LG(M);
    for(u_int m = 0; m < M; ++m) {
      mu[m] = pw[0][m];
      wgt[m] = pw[1][m];
    }
  }
  else {
    mu[0] = 1;
    wgt[0] = 2;
  }

}

void DGFEMSpace1D::st2coe(const u_int q, u_int& qx, u_int& qt) {
  qx = q%K;
  qt = (q-qx)/K;
}

void DGFEMSpace1D::BuildUniMAT() {
  u_int qx, qt, px, pt;
  VEC<double> x = TemQuad_Gauss.points();
  VEC<double> w = TemQuad_Gauss.weight();
  EVEC pl = Lagrange_Poly(-1);
  EVEC pr = Lagrange_Poly(1);
  EVEC prime;
  M1.setZero();
  for(u_int j = 0; j < K; ++j) {
    M1(j,j) = 0.5*w[j];
  }
  Mx.setZero();
  for(u_int l = 0; l < K; ++l) {
    for(u_int q = 0; q < K*K; ++q) {
      st2coe(q, qx, qt);
      prime = Lagrange_PolyG(x[qx]);
      Mx(l,q) = 0.5*w[qx]*w[qt]*prime[l];
    }
  }
  Fx_plus.setZero();
  Fx_minus.setZero();
  for(u_int l = 0; l < K; ++l) {
    for(u_int q = 0; q < K*K; ++q) {
      st2coe(q, qx, qt);
      Fx_plus(l,q) = 0.5*w[qt]*pr[qx]*pr[l];
      Fx_minus(l,q) = 0.5*w[qt]*pl[qx]*pl[l];
    }
  }
  K1.setZero();
  for(u_int p = 0; p < K*K; ++p) {
    for(u_int q = 0; q < K*K; ++q) {
      st2coe(p, px, pt);
      st2coe(q, qx, qt);
      if(px == qx) {
        prime = Lagrange_PolyG(x[qt]);
        K1(p,q) = 0.5*w[px]*pr[pt]*pr[qt]
        - 0.5*w[qx]*w[qt]*prime[pt];
      }
    }
  }
  Kx.setZero();
  for(u_int p = 0; p < K*K; ++p) {
    for(u_int q = 0; q < K*K; ++q) {
      st2coe(p, px, pt);
      st2coe(q, qx, qt);
      if(pt == qt) {
        prime = Lagrange_PolyG(x[px]);
        Kx(p,q) = 0.5*w[px]*w[pt]*prime[qx];
      }
    }
  }
  F0.setZero();
  for(u_int p = 0; p < K*K; ++p) {
    for(u_int l = 0; l < K; ++l) {
      st2coe(p, px, pt);
      if(l == px) {
        F0(p,l) = 0.5*w[l]*pl[pt];
      }
    }
  }

  //std::cout << "M1:" << std::endl;
  //std::cout << M1 << std::endl;
  //std::cout << "Mx:" << std::endl;
  //std::cout << Mx << std::endl;
  //std::cout << "Fx_plus:" << std::endl;
  //std::cout << Fx_plus << std::endl;
  //std::cout << "Fx_minus:" << std::endl;
  //std::cout << Fx_minus << std::endl;
  //std::cout << "K1:" << std::endl;
  //std::cout << K1 << std::endl;
  //std::cout << "Kx:" << std::endl;
  //std::cout << Kx << std::endl;
  //std::cout << "F0:" << std::endl;
  //std::cout << F0 << std::endl;
  //abort();
}

void DGFEMSpace1D::Projection(u_int cell, func I0, double t, bU& u) {
  VEC<double> p = QUADINFO[cell].points();
  for(u_int m = 0; m < M; ++m) {
    for(u_int g = 0; g < K; ++g) {
      u[m][g] = I0(mu[m], p[g], t);
    }
  }
}

EVEC DGFEMSpace1D::Composition(const SOL& I, u_int cell, double x) {
  EVEC u(M);
  u.setZero();
  VEC<double> lv(2), gv(2);
  EVEC V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];

  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Lagrange_Poly(lp);
  for(u_int m = 0; m < M; ++m) {
    for(u_int k = 0; k < K; ++k) {
      u[m] += I[cell][m][k]*V[k];
    }
  }

  return u;
}

void DGFEMSpace1D::init(func I0) {
  for(u_int i = 0; i < Nx; ++i) {
    Projection(i, I0, 0, I[i]);
  }
}

double DGFEMSpace1D::cal_dt(const SOL& I, const double t) {
  return 0.5*(mesh[1]-mesh[0])/mu[0]/(2*K-1);
}

void DGFEMSpace1D::print_DG_coe(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  VEC<double> x = TemQuad_Gauss.points();
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int m = 0; m < M; ++m) {
      for(u_int k = 0; k < K; ++k) {
        os << I[i][m][k] << " ";
      }
      os << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_integral(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  VEC<double> x = TemQuad_Gauss.points();
  for(u_int i = 0; i < Nx; ++i) {
    VEC<double> w = QUADINFO[i].weight();
    VEC<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < x.size(); ++g) {
      os << p[g] << " "  << w[g];
      for(u_int m = 0; m < M; ++m) {
        os  << " " << I[i][m][g];
      }
      os << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_average(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  for(u_int i = 0; i < Nx; ++i) {
      os << 0.5*(mesh[i]+mesh[i+1]) << " " << I[i][0] << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

