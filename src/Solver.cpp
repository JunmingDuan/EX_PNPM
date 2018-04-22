#include "DGFEMSpace1D.h"

void DGFEMSpace1D::solve_leqn(const EMAT& A, const EVEC& rhs, EVEC& u) {
  solver.compute(A);
  u = solver.solve(rhs);
}

void DGFEMSpace1D::STDG_predictor(const SOL& I, const double dt, const u_int cell, const u_int m, ST_ele& W) {
  double hi = mesh[cell+1] - mesh[cell];
  double para = -1e1;
  //A = hi*K1 + (dt*mu[m])*Kx - para*(hi*dt)*MM;//all
  A = hi*K1 + (dt*mu[m])*Kx;//without source
  //A = hi*K1 - (hi*dt)*MM;//without advection
  rhs = hi*F0*I[cell][m];
  solve_leqn(A, rhs, W);
}

double DGFEMSpace1D::RAD_BE_unsteady(SOL& I,
    func_para sigma_t, func_para sigma_s, src q,
    const double t, const double dt, func BL, func BR) {//period
  VEC<ST_ele> W_list;
  W_list.resize(Nx);
  Reconstruct(I, I1);
  //I1 = I;
  for(u_int m = 0; m < M; ++m) {//for each light direction
    Boundary(BL, BR, BD_L, BD_R, mu[m], t, dt);
    if(mu[m] >= 0) {
      for(u_int i = 0; i < Nx; ++i) {
        STDG_predictor(I1, dt, i, m, W_list[i]);
      }
      for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}^{l}
        double hi = mesh[i+1] - mesh[i];

        A = hi*M1;
        if(i == 0) {//!!! notice that we use I1 below
          rhs = hi*M1*I[i][m] + (mu[m]*dt)*Mx*W_list[i] - dt*Fx_plus*W_list[i] + dt*Fx_minus*W_list[Nx-1];
        }
        else {
          rhs = hi*M1*I[i][m] + (mu[m]*dt)*Mx*W_list[i] - dt*Fx_plus*W_list[i] + dt*Fx_minus*W_list[i-1];
        }
        double para = -1e1;
        //rhs += (para*hi*dt)*Mt*W_list[i];
        solve_leqn(A, rhs, I[i][m]);

      }//end i
      //abort();
    }//end ->
    else {
      u_int i = Nx-1;
      while(i >= 0) {
        if(i-- == 0) break;
      }//end i
    }//end <-
  }//end M

  return 1;
}

//double DGFEMSpace1D::RAD_BE_unsteady(const SOL& In, SOL& I,
//func_para sigma_t, func_para sigma_s, func q,
//const double t, const double dt, SOL& I_new, func BL, func BR) {
//for(u_int m = 0; m < M; ++m) {//for each light direction
//Boundary(BL, BR, BD_L, BD_R, mu[m], t, dt);
//if(mu[m] >= 0) {
//for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}^{l}
//std::cout << "########i: " << i << " ########" << std::endl;
//double hi = mesh[i+1] - mesh[i];
//STDG_predictor(In, dt, i, m, W);

//A = hi*M1;
//if(i == 0) {
//FLUX = BD_L*Lagrange_Poly(-1);
//rhs = hi*M1*In[i][m] + (mu[m]*dt)*Mx*W - dt*Fx_plus*W + FLUX;
//}
//else {
//rhs = hi*M1*In[i][m] + (mu[m]*dt)*Mx*W - dt*Fx_plus*W + dt*Fx_minus*W_uw;
//}
//solve_leqn(A, rhs, I_new[i][m]);

//W_uw = W;//update the upwindign W

//}//end i
//}//end ->
//else {
//u_int i = Nx-1;
//while(i >= 0) {
//if(i-- == 0) break;
//}//end i
//}//end <-
//}//end M

//return 1;
//}

void DGFEMSpace1D::run_unsteady(func_para sigma_t, func_para sigma_s, src q, func BL, func BR, double t_end) {
  std::cout.precision(16);
  std::cout << std::showpos;
  std::cout.setf(std::ios::scientific);
  double t(0), dt(0);

  while (t < t_end) {
    dt = cal_dt(I, t);
    dt = std::min(dt, t_end-t);
    RAD_BE_unsteady(I, sigma_t, sigma_s, q, t, dt, BL, BR);
    t += dt;
    std::cout << "t: " << t << ", dt: " << dt << std::endl;
  }
}

