/**
 * @file TemplateQuadrature.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef TEMPLATEQUADRTURE_H
#define TEMPLATEQUADRTURE_H

#include <iostream>
#include <cmath>
#include "VEC.h"

class TemplateQuadrature {
  protected:
    u_int np;
    VEC<double> pnt;
    VEC<double> wei;
  public:
    VEC<VEC<double> > LGL(u_int np);//legendre-Gauss-Lobatto
    VEC<VEC<double> > LG(u_int np);//legendre-Gauss
    void set_np(u_int np);
    void set_points(const VEC<double>& p);
    void set_weight(const VEC<double>& w);
    VEC<double> points();
    VEC<double> weight();
    void print(std::ostream&);

};

#endif

