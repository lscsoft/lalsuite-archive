#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

extern void __nested_MOD_nestrun(int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
    int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, char *root,
    int *seed, int *pWrap, int *fb, int *resume, void (*Loglike)(double *Cube, int *n_dim, int *n_par,
    double *lnew), int *context);

void run(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes,
    int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume,
    void (*LogLike)(double *Cube, int *n_dim, int *n_par, double *lnew), int context)
    {
	int i;
	for (i = strlen(root); i < 100; i++)
		root[i] = ' ';
	
        __nested_MOD_nestrun(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, LogLike, &context);
    }   

// Now an example, sample an egg box likelihood

void LogLike(double *Cube, int *ndim, int *npars, double *lnew)
{
    double chi = 1.0;
    int i;
    for(i = 0; i < *ndim; i++)
    {
        double x = Cube[i]*10.0*M_PI;
        chi *= cos(x/2.0);
            Cube[i] = x;
    }
    *lnew = powf(chi + 2.0, 5.0);
}

int main(int argc, char *argv[])
{
    int mmodal = 1;
    int ceff = 0;
    int nlive = 1000;
    double efr = 0.8;
    double tol = 0.5;
    int ndims = 2;
    int nPar = 2;
    int nClsPar = 2;
    int updInt = 100;
    double Ztol = -1.e90;
    int maxModes = 100;
    int pWrap[] = {0,0};
    char root[100] = "chains/1-";
    int seed = -1;
    int fb = 1;
    int resume = 1;
    int context = 0;

    run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb,
            resume, LogLike, context);
}
