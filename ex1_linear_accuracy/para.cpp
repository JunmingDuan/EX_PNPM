#include "para.h"

u_int M = 1;
u_int K = 1;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
u_int BDL = 2; u_int BDR = 2;
double EPS = 0;

double c = 1;

double tol_I = 1e-14;
double TOL = 1e-14;
u_int MAXITE = 1e2;

