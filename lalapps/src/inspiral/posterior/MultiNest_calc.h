#include "nest_calc.h"

extern int multinest_seg;

extern void __nested_MOD_nestrun(int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims, int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, char *root, int *seed, 
int *pWrap, int *fb, int *resume, void (*Loglike)(double *Cube, int *n_dim, int *n_par, double *lnew), int *context);

void MultiNestRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes, int updInt, double Ztol, char
root[], int seed, int *pWrap, int fb, int resume, void (*LogLike)(double *Cube, int *n_dim, int *n_par, double *lnew), int context);

void LogLike(double *Cube, int *ndim, int *npars, double *lnew);

void MultiNestZ(UINT4 Nlive, LALMCMCInput *MCMCinput);
