#include "DGFEMSpace1D.h"
#include "conWENO3.h"
#include "conWENO5.h"
#include "conWENO7.h"

void DGFEMSpace1D::Reconstruct(const SOL& I0, SOL& I) {
  if(K == 1) {
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        I[i][m][0] = I0[i][m][0];
      }
    }
  }
  else if(K == 2) {//WENO3
    u_int n1, n2, n3;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-1+Nx)%Nx;
        n2 = (i+Nx)%Nx;
        n3 = (i+1+Nx)%Nx;
        I[i][m][0] = WENONI21(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
        I[i][m][1] = WENONI22(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
      }
    }
  }
  else if(K == 3) {//WENO5
    u_int n1, n2, n3, n4, n5;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-2+Nx)%Nx;
        n2 = (i-1+Nx)%Nx;
        n3 = (i+Nx)%Nx;
        n4 = (i+1+Nx)%Nx;
        n5 = (i+2+Nx)%Nx;
        I[i][m][0] = WENONI31(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[i][m][1] = WENONI32(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[i][m][2] = WENONI33(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
      }
    }
  }
  else if(K == 4) {//WENO7
    u_int n1, n2, n3, n4, n5, n6, n7;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-3+Nx)%Nx;
        n2 = (i-2+Nx)%Nx;
        n3 = (i-1+Nx)%Nx;
        n4 = (i+Nx)%Nx;
        n5 = (i+1+Nx)%Nx;
        n6 = (i+2+Nx)%Nx;
        n7 = (i+3+Nx)%Nx;
        I[i][m][0] = WENONI41(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][1] = WENONI42(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][2] = WENONI43(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][3] = WENONI44(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
      }
    }
  }
}


