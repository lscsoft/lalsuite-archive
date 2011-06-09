#include "nest_calc.h"

extern int MNSeg;

extern void __nested_MOD_nestrun(int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, 
char *, int *, int *, int *, int *, int *, int *, double *, void (*Loglike)(double *, int *, int *, double *), void (*dumper)(int *, int *, 
int *, double **, double **, double *, double *, double *), int *context);

void MultiNestRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes,
int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, int initMPI, double logZero, 
void (*LogLike)(double *, int *, int *, double *), void (*dumper)(int *, int *, int *, double **, double **, double *, 
double *, double *, double *), int context);

void LogLike(double *Cube, int *ndim, int *npars, double *lnew);

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr);

void MultiNestZ(UINT4 Nlive, LALMCMCInput *MCMCinput);
