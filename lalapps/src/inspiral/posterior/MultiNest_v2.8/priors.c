#include "priors.h"

double skillingScalePrior(double r)
{
	return ( ( 1.0 - r ) / r );
}

double logPrior(double r, double x1, double x2)
{
	double lx1, lx2;
	lx1 = log( x1 );
	lx2 = log( x2 );
	return exp( lx1 + r * ( lx2 - lx1 ) );
}


double flatPrior(double r, double x1, double x2)
{
	return x1 + r * ( x2 - x1 );
}


double gaussianPrior(double r, double x1, double x2)
{
	if( x2 <= 0.0 )
	{
		printf("sigma <= 0 in routine gaussianPrior\n");
		exit(1);
	}
	
	if( r <= 1.E-16 || ( 1.0 - r ) <= 1.E-16 )
       		return ( -1.0E32 );
      	else
       		return ( x1 + x2 * sqrt( 2.0 ) * dierfc( 2.0 * ( 1.0 - r ) ) );
}


double sinPrior(double r, double x1, double x2)
{
	double cx1, cx2;
      	double deg2rad = 0.017453292;

      	cx1 = cos( x1 * deg2rad );
      	cx2 = cos( x2 * deg2rad );
      	return 1.0 * acos( cx1 + r * ( cx2 - cx1 ) );
}


/* Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma] */
double cauchyPrior(double r, double x0, double gamma)
{
	double Pi = 4.0 * atan( 1.0 );
	return x0 + gamma * tan( Pi * ( r - 0.5 ) );
}
		

/* Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma] */
double logNormalPrior(double r, double a, double sigma)
{
	double bracket = sigma * sigma + sigma * sqrt( 2.0 ) * dierfc( 2.0 * r );
      	return a * exp( bracket );
}

/* Return inverse of the complimentary error function */
double dierfc(double y)
{
	double infinity = 5.0;
	double qa = 9.16461398268964E-01,
      	qb = 2.31729200323405E-01, 
      	qc = 4.88826640273108E-01, 
      	qd = 1.24610454613712E-01, 
      	q0 = 4.99999303439796E-01, 
      	q1 = 1.16065025341614E-01, 
      	q2 = 1.50689047360223E-01, 
      	q3 = 2.69999308670029E-01, 
      	q4 = -7.28846765585675E-02,
      	pa = 3.97886080735226000E+00, 
      	pb = 1.20782237635245222E-01, 
      	p0 = 2.44044510593190935E-01, 
      	p1 = 4.34397492331430115E-01, 
      	p2 = 6.86265948274097816E-01, 
      	p3 = 9.56464974744799006E-01, 
      	p4 = 1.16374581931560831E+00, 
      	p5 = 1.21448730779995237E+00, 
      	p6 = 1.05375024970847138E+00, 
      	p7 = 7.13657635868730364E-01, 
      	p8 = 3.16847638520135944E-01, 
      	p9 = 1.47297938331485121E-02, 
      	p10 = -1.05872177941595488E-01, 
      	p11 = -7.43424357241784861E-02,
      	p12 = 2.20995927012179067E-03, 
      	p13 = 3.46494207789099922E-02, 
      	p14 = 1.42961988697898018E-02, 
      	p15 = -1.18598117047771104E-02, 
      	p16 = -1.12749169332504870E-02, 
      	p17 = 3.39721910367775861E-03, 
      	p18 = 6.85649426074558612E-03, 
      	p19 = -7.71708358954120939E-04, 
      	p20 = -3.51287146129100025E-03, 
      	p21 = 1.05739299623423047E-04, 
      	p22 = 1.12648096188977922E-03;
      
	if( y == 0.0 ) return infinity;
	double z, w, u, s, t, x;
	z = y;
	if( y > 1.0 ) z = 2.0 - y;
	w = qa - log( z );
	u = sqrt( w );
	s = ( qc + log( u ) ) / w;
	t = 1.0 / ( u + qb );
	x = u * ( 1.0 - s * ( 0.5 + s * qd ) ) - ( ( ( ( q4 * t + q3 ) * t + q2 ) * t + q1 ) * t + q0 ) * t;
	t = pa / ( pa + x );
	u = t - 0.5;
	s = ( ( ( ( ( ( ( ( ( p22 * u + p21 ) * u + p20 ) * u + p19 ) * u + p18 ) * u + p17 ) * u + p16 ) * 
		u + p15 ) * u + p14 ) * u + p13 ) * u + p12;
	s = ( ( ( ( ( ( ( ( ( ( ( ( s * u + p11 ) * u + p10 ) * u + p9 ) * u + p8 ) * u + p7 ) * u + p6 ) * 
		u + p5 ) * u + p4 ) * u + p3 ) * u + p2 ) * u + p1 ) * u + p0 ) * t - z * exp( x * x - pb );
      	x = x + s * ( 1.0 + x * s );
      	if( y > 1.0 ) x = - x;
	return x;
}
