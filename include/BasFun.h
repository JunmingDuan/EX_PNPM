/**
 * @file BasFun.h
 * @brief Basis function on [-1,1], i.e. Legendre polynomial not normalized.
 * or Lagrange_intepolated polynomial based on Legendre-Gauss points
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-04
 */

#ifndef BASFUN_H
#define BASFUN_H

#include <vector>
#include <iostream>
#include "para.h"

EVEC Legendre_Poly(double x);
EVEC Legendre_PolyG(double x);

EVEC Lagrange_Poly(double x);
EVEC Lagrange_PolyG(double x);

#endif //BASFUN_H

