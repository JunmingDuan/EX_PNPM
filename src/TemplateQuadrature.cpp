/**
 * @file Quadrature.cpp
 * @brief Quadrature points and weights on [-1,1]
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "TemplateQuadrature.h"

VEC<VEC<double> > TemplateQuadrature::LGL(u_int np) {
  VEC<VEC<double> > pw;
  pw.resize(2);
  pw[0].resize(np);
  pw[1].resize(np);
  switch (np) { //from wiki
    case 2: {
              pw[0][0] = -1, pw[0][1] = 1;
              pw[1][0] = 1, pw[1][1] = 1;
              break;
            }
    case 3: {
              pw[0][0] = -1, pw[0][1] = 0, pw[0][2] = 1;
              pw[1][0] = 1./3, pw[1][1] = 4./3, pw[1][2] = 1./3;
              break;
            }
    case 4: {
              pw[0][0] = -1, pw[0][1] = -1/sqrt(5), pw[0][2] = 1/sqrt(5), pw[0][3] = 1;
              pw[1][0] = 1./6,  pw[1][1] = 5./6, pw[1][2] = 5./6, pw[1][3] = 1./6;
              break;
            }
    case 5: {
              pw[0][0] = -1, pw[0][1] = -sqrt(3./7), pw[0][2] = 0, pw[0][3] = sqrt(3./7), pw[0][4] = 1;
              pw[1][0] = 1./10, pw[1][1] = 49./90, pw[1][2] = 32./45, pw[1][3] = 49./90, pw[1][4] = 1./10;
              break;
            }
    case 6: {
              pw[0][0] = -1, pw[0][1] = -sqrt(1./3+2*sqrt(7)/21), pw[0][2] = -sqrt(1./3-2*sqrt(7)/21),
              pw[0][3] = sqrt(1./3-2*sqrt(7)/21), pw[0][4] = sqrt(1./3+2*sqrt(7)/21), pw[0][5] = 1;
              pw[1][0] = 1./15, pw[1][1] = (14-sqrt(7))/30, pw[1][2] = (14+sqrt(7))/30,
              pw[1][3] = (14+sqrt(7))/30, pw[1][4] = (14-sqrt(7))/30, pw[1][5] = 1./15;
              break;
            }
    case 7: {
              pw[0][0] = -1, pw[0][1] = -sqrt(5./11+2./11*sqrt(5./3)), pw[0][2] = -sqrt(5./11-2./11*sqrt(5./3)),
              pw[0][3] = 0,
              pw[0][4] = sqrt(5./11-2./11*sqrt(5./3)), pw[0][5] = sqrt(5./11+2./11*sqrt(5./3)), pw[0][6] = 1;
              pw[1][0] = 1./21, pw[1][1] = (124-7*sqrt(15))/350, pw[1][2] = (124+7*sqrt(15))/350,
              pw[1][3] = 256./525,
              pw[1][4] = (124+7*sqrt(15))/350, pw[1][5] = (124-7*sqrt(15))/350, pw[1][6] = 1./21;
              break;
            }
    default:{
               std::cout << "Wrong Legendre Gauss Lobatto choice!" << std::endl;
               break;
            }
  }
  return pw;
}

VEC<VEC<double> > TemplateQuadrature::LG(u_int np) {
  VEC<VEC<double> > pw;
  pw.resize(2);
  pw[0].resize(np);
  pw[1].resize(np);
  switch (np) { //from wiki
    case 1: {
              pw[0][0] = 0;
              pw[1][0] = 2;
              break;
            }
    case 2: {
              pw[0][0] = -sqrt(3)/3, pw[0][1] = sqrt(3)/3;
              pw[1][0] = 1, pw[1][1] = 1;
              break;
            }
    case 3: {
              pw[0][0] = -sqrt(15)/5, pw[0][1] = 0, pw[0][2] = sqrt(15)/5;
              pw[1][0] = 5./9, pw[1][1] = 8./9, pw[1][2] = 5./9;
              break;
            }
    case 4: {
              pw[0][0] = -0.8611363115940520, pw[0][1] = -0.3399810435848560, pw[0][2] = 0.3399810435848560, pw[0][3] = 0.8611363115940520;
              pw[1][0] = 0.3478548451374530,  pw[1][1] = 0.6521451548625460, pw[1][2] = 0.6521451548625460, pw[1][3] = 0.3478548451374530;
              break;
            }
    case 5: {
              pw[1][0] = 0.5688888888888889; pw[0][0] = 0.0000000000000000;
              pw[1][1] = 0.4786286704993665; pw[0][1] = -0.5384693101056831;
              pw[1][2] = 0.4786286704993665; pw[0][2] = 0.5384693101056831;
              pw[1][3] = 0.2369268850561891; pw[0][3] = -0.9061798459386640;
              pw[1][4] = 0.2369268850561891; pw[0][4] = 0.9061798459386640;
              break;
            }
    case 6: {
              pw[0][0] = -0.932469514203152, pw[0][1] = -0.661209386466264, pw[0][2] = -0.238619186083197,
                pw[0][3] =  0.238619186083197, pw[0][4] =  0.661209386466264, pw[0][5] =  0.932469514203152;
              pw[1][0] = 0.171324492379170, pw[1][1] = 0.360761573048138, pw[1][2] = 0.467913934572691,
                pw[1][3] = 0.467913934572691, pw[1][4] = 0.360761573048138, pw[1][5] = 0.171324492379170;
              break;
            }
    case 8: {
              pw[0][0] = -0.960289856497536, pw[0][1] = -0.796666477413627, pw[0][2] = -0.525532409916329, pw[0][3] = -0.183434642495650,
                pw[0][4] =  0.183434642495650, pw[0][5] =  0.525532409916329, pw[0][6] =  0.796666477413627, pw[0][7] =  0.960289856497536;
              pw[1][0] =  0.101228536290376, pw[1][1] =  0.222381034453374, pw[1][2] =  0.313706645877887, pw[1][3] =  0.362683783378362,
                pw[1][4] =  0.362683783378362, pw[1][5] =  0.313706645877888, pw[1][6] =  0.222381034453374, pw[1][7] =  0.101228536290376;
              break;
            }
    default:{
              std::cout << "Wrong Legendre Gauss choice!" << std::endl;
              break;
            }
  }
  return pw;
}


void TemplateQuadrature::set_np(u_int n) {
  np = n;
}

void TemplateQuadrature::set_weight(const VEC<double>& w) {
  wei = w;
}

void TemplateQuadrature::set_points(const VEC<double>& p) {
  pnt = p;
}

VEC<double> TemplateQuadrature::points() {
  return pnt;
}

VEC<double> TemplateQuadrature::weight() {
  return wei;
}

void TemplateQuadrature::print(std::ostream& os) {
  os << "Number of quadrature points: " << np << std::endl;
  for(u_int i = 0; i < np; ++i) {
    os << pnt[i] << "\t";
  }
  os << "\n";
  for(u_int i = 0; i < np; ++i) {
    os << wei[i] << "\t";
  }
  os << std::endl;
}

