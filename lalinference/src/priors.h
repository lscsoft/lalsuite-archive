#include "math.h"

double skillingScalePrior(double r);
double logPrior(double r, double x1, double x2);
double flatPrior(double r, double x1, double x2);
double powerPrior(double power, double r, double x1, double x2);
double gaussianPrior(double r, double mu, double sigma);
double sinPrior(double r, double x1, double x2);
double cauchyPrior(double r, double x0, double gamma);
double logNormalPrior(double r, double a, double sigma);
double dierfc(double y);

