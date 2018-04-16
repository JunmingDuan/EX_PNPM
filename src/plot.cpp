#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "DGFEMSpace1D.h"
#define  BUFSIZE 256

int DGFEMSpace1D::plot_I(engine *ep, const VEC<double>& mesh, const SOL& In)
{
  u_int N = In.size();
  u_int Np = In[0][0].size();
  VEC<double> x0(N*Np), I0(N*Np);
  for(u_int i = 0; i < N; ++i) {
    x0[2*i] = mesh[i];
    x0[2*i+1] = mesh[i+1];
    I0[2*i] = In[i][0][0];
    I0[2*i+1] = In[i][0][1];
  }
	mxArray *I = NULL, *x = NULL;

	x = mxCreateDoubleMatrix(1, N*Np, mxREAL);
	I = mxCreateDoubleMatrix(1, N*Np, mxREAL);
  memcpy((void *)mxGetPr(x), (void *)(&x0[0]), x0.size()*sizeof(double));
  memcpy((void *)mxGetPr(I), (void *)(&I0[0]), I0.size()*sizeof(double));
	/*
	 * Place the variable T into the MATLAB workspace
	 */
	engPutVariable(ep, "I", I);
	engPutVariable(ep, "x", x);

	engEvalString(ep, "plot(x,I,'-b');");
	engEvalString(ep, "title('Marshak wave');");
	engEvalString(ep, "xlabel('x');");
	engEvalString(ep, "ylabel('I');");
	engEvalString(ep, "hold on;");
	mxDestroyArray(I);
	mxDestroyArray(x);

	return EXIT_SUCCESS;
}

