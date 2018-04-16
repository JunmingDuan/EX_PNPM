#include "DGFEMSpace1D.h"

double DGFEMSpace1D::cal_norm_I(const SOL& s1, const SOL& s2, int n) {
  double RES(0);
  VEC<double> norm(M,0);
  EVEC tmp1(M), tmp2(M);
  if(n == 10) {//max norm
    //for(u_int i = 0; i < Nx; ++i) {
      //VEC<double> p = QUADINFO[i].points();
      //VEC<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g]);
        //tmp2 = Composition(s2,i,p[g]);
        //for(u_int d = 0; d < M; ++d) {
          //norm[d] = std::max(fabs(tmp1[d]-tmp2[d]), norm[d]);
        //}
      //}
    //}
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int d = 0; d < M; ++d) {
        norm[d] = std::max(fabs(s1[i][d][0]-s2[i][d][0]), norm[d]);
        RES = std::max(norm[d], RES);
        //norm[d] = std::max((s1[i][d]-s2[i][d]).norm(), norm[d]);
      }
    }
    return RES;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return RES;
  }
}

double exact(const double mu, const double x, const double t) {
  return pow(mu,2)*pow(cos(M_PI*(x+t)), 4);
}

double DGFEMSpace1D::cal_err(const SOL& s1, int n, double t_end) {
  EVEC tmp1;
  double tmp2, NORM(0);
  double pnt;
  VEC<double> norm(M,0);
  u_int Ne(8);
  VEC<VEC<double>> pw = TemQuad_Gauss.LG(Ne);
  VEC<VEC<double>> x(Nx);
  VEC<double> jab(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    x[i].resize(Ne);
  }
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    for(u_int j = 0; j < Ne; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[i][j]);
    }
    jab[i] = local_to_global_jacobian(lv, gv);
  }

  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] += pow(tmp1[m]-tmp2, 2)*pw[1][g];
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, sqrt(norm[m]/Nx));
    }
    return NORM;
  }
  else if(n == 1) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] += fabs(tmp1[m]-tmp2)*pw[1][g];
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, norm[m]/Nx);
    }
    return NORM;
  }
  else if(n == 0) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] = std::max(norm[m], tmp1[m]-tmp2);
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, norm[m]);
    }
    return NORM;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return NORM;
  }
}


