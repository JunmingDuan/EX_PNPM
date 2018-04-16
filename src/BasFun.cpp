#include "BasFun.h"

EVEC Legendre_Poly(double x) {
  EVEC poly(K);
  switch (K) {
    case 1:
      poly[0] = 1;
      return poly;
    case 2:
      poly[0] = 1;
      poly[1] = x;
      return poly;
    case 3:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      return poly;
    case 4:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      poly[3] = (5*x*x*x-3*x)/2;
      return poly;
    case 5:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      poly[3] = (5*x*x*x-3*x)/2;
      poly[4] = (35*x*x*x*x-30*x*x+3)/8;
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

/**
 * @brief PG Gradient of polynomial P.
 *
 * @param x
 *
 * @return A vector of length DIM.
 */
EVEC Legendre_PolyG(double x) {
  EVEC poly(K);
  switch (K) {
    case 1:
      poly[0] = 0;
      return poly;
    case 2:
      poly[0] = 0;
      poly[1] = 1;
      return poly;
    case 3:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      return poly;
    case 4:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      poly[3] = (15*x*x-3)/2;
      return poly;
    case 5:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      poly[3] = (15*x*x-3)/2;
      poly[4] = (140*x*x*x-60*x)/8;
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

EVEC Lagrange_Poly(double x) {
  EVEC poly(K);
  double x0, x1, x2, x3, x4;
  switch (K) {
    case 1:
      poly[0] = 1;
      return poly;
    case 2:
      x0 = -sqrt(3)/3;
      x1 = sqrt(3)/3;
      poly[0] = (x-x1)/(x0-x1);
      poly[1] = (x-x0)/(x1-x0);
      return poly;
    case 3:
      x0 = -sqrt(15)/5;
      x1 = 0;
      x2 = sqrt(15)/5;
      poly[0] = (x-x1)*(x-x2)/((x0-x1)*(x0-x2));
      poly[1] = (x-x0)*(x-x2)/((x1-x0)*(x1-x2));
      poly[2] = (x-x0)*(x-x1)/((x2-x0)*(x2-x1));
      return poly;
    case 4:
      x0 = -0.8611363115940520;
      x1 = -0.3399810435848560;
      x2 = 0.3399810435848560;
      x3 = 0.8611363115940520;
      poly[0] = (x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3));
      poly[1] = (x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3));
      poly[2] = (x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3));
      poly[3] = (x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2));
      return poly;
    case 5:
      x0 = 0.0000000000000000;
      x1 = -0.5384693101056831;
      x2 = 0.5384693101056831;
      x3 = -0.9061798459386640;
      x4 = 0.9061798459386640;
      poly[0] = (x-x1)*(x-x2)*(x-x3)*(x-x4)/((x0-x1)*(x0-x2)*(x0-x3)*(x0-x4));
      poly[1] = (x-x0)*(x-x2)*(x-x3)*(x-x4)/((x1-x0)*(x1-x2)*(x1-x3)*(x1-x4));
      poly[2] = (x-x0)*(x-x1)*(x-x3)*(x-x4)/((x2-x0)*(x2-x1)*(x2-x3)*(x2-x4));
      poly[3] = (x-x0)*(x-x1)*(x-x2)*(x-x4)/((x3-x0)*(x3-x1)*(x3-x2)*(x3-x4));
      poly[4] = (x-x0)*(x-x1)*(x-x2)*(x-x3)/((x4-x0)*(x4-x1)*(x4-x2)*(x4-x3));
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

/**
 * @brief PG Gradient of polynomial P.
 *
 * @param x
 *
 * @return A vector of length DIM.
 */
EVEC Lagrange_PolyG(double x) {
  EVEC poly(K);
  double x0, x1, x2, x3, x4;
  switch (K) {
    case 1:
      poly[0] = 0;
      return poly;
    case 2:
      x0 = -sqrt(3)/3;
      x1 = sqrt(3)/3;
      poly[0] = 1./(x0-x1);
      poly[1] = 1./(x1-x0);
      return poly;
    case 3:
      x0 = -sqrt(15)/5;
      x1 = 0;
      x2 = sqrt(15)/5;
      poly[0] = (2.*x-x1-x2)/((x0-x1)*(x0-x2));
      poly[1] = (2.*x-x0-x2)/((x1-x0)*(x1-x2));
      poly[2] = (2.*x-x0-x1)/((x2-x0)*(x2-x1));
      return poly;
    case 4:
      x0 = -0.8611363115940520;
      x1 = -0.3399810435848560;
      x2 = 0.3399810435848560;
      x3 = 0.8611363115940520;
      poly[0] = (3.*x*x - 2.*x*(x1+x2+x3) + x1*x2+x2*x3+x3*x1)/((x0-x1)*(x0-x2)*(x0-x3));
      poly[1] = (3.*x*x - 2.*x*(x0+x2+x3) + x0*x2+x2*x3+x3*x0)/((x1-x0)*(x1-x2)*(x1-x3));
      poly[2] = (3.*x*x - 2.*x*(x0+x1+x3) + x0*x1+x1*x3+x3*x0)/((x2-x0)*(x2-x1)*(x2-x3));
      poly[3] = (3.*x*x - 2.*x*(x0+x1+x2) + x0*x1+x1*x2+x2*x0)/((x3-x0)*(x3-x1)*(x3-x2));
      return poly;
    case 5:
      x0 = 0.0000000000000000;
      x1 = -0.5384693101056831;
      x2 = 0.5384693101056831;
      x3 = -0.9061798459386640;
      x4 = 0.9061798459386640;
      poly[0] = ((x-x2)*(x-x3)*(x-x4) + (x-x1)*(x-x3)*(x-x4) + (x-x1)*(x-x2)*(x-x4) + (x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3)*(x0-x4));
      poly[1] = ((x-x2)*(x-x3)*(x-x4) + (x-x0)*(x-x3)*(x-x4) + (x-x0)*(x-x2)*(x-x4) + (x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3)*(x1-x4));
      poly[2] = ((x-x1)*(x-x3)*(x-x4) + (x-x0)*(x-x3)*(x-x4) + (x-x0)*(x-x1)*(x-x4) + (x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3)*(x2-x4));
      poly[3] = ((x-x1)*(x-x2)*(x-x4) + (x-x0)*(x-x2)*(x-x4) + (x-x0)*(x-x1)*(x-x4) + (x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2)*(x3-x4));
      poly[4] = ((x-x1)*(x-x2)*(x-x3) + (x-x0)*(x-x2)*(x-x3) + (x-x0)*(x-x1)*(x-x3) + (x-x0)*(x-x1)*(x-x2))/((x4-x0)*(x4-x1)*(x4-x2)*(x4-x3));
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

