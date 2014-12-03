#include <math.h>

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>

//Definitions in SI Units
/*
#define G (6.6726e-11)
#define c (299792458.)
#define Msun (1.9889e30)
#define Msun_sec (G*Msun/(c*c*c))
#define GPC_sec (3.08567818585e25/c)
#define ER_sec (2.12514687e-2)		//Earth radius in seconds
#define YEAR (31557600.)
#define DAY (86400.)
#define pi (M_PI)
#define gamma (M_EULER)*/

#define c ((REAL8)(LAL_C_SI))
#define G ((REAL8)(LAL_G_SI))
#define pi ((REAL8)(LAL_PI))
#define pc ((REAL8)(LAL_PC_SI))
#define gamma ((REAL8)(LAL_GAMMA))
#define Msun ((REAL8)(LAL_MSUN_SI))
#define Msun_sec ((REAL8)(G*Msun/(c*c*c)))
#define GPC_sec ((REAL8)(pc*pow(10.,9.)/c))

//Function prototypes
INT4 XLALSearchArr(REAL8 *arr, REAL8 time);
REAL8 XLALInitialP(REAL8 p0, void *params);
INT4 XLALComputeEnergyAng(REAL8 p, REAL8 e, REAL8 a, INT4 dir, REAL8 *Eng, REAL8 *Ang);
REAL8 XLALDFDr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
REAL8 XLALDVphiVtDr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
REAL8 XLALD2FDrDE(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
REAL8 XLALD2FDrDL(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
REAL8 XLALEdot3PN_Effective(REAL8 a, REAL8 p, REAL8 e, REAL8 nu);
REAL8 XLALLdot3PN_Effective(REAL8 a, REAL8 p, REAL8 e, REAL8 nu, INT4 dir);
REAL8 XLALEnergyAngDot(REAL8 a, REAL8 p, REAL8 e, REAL8 nu, INT4 dir, REAL8* dEdt, REAL8* dLdt);
INT4 XLALPEdot(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 nu, REAL8 dir, REAL8 *pdot,
	REAL8 *edot);
REAL8 XLALDphiDt_geo(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
REAL8 XLALDphiDt_sf(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi, REAL8 nu);
REAL8 XLALDpsiDt_geo(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi);
REAL8 XLALDpsiDt_sf(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi, REAL8 nu);
REAL8 XLALPlso(REAL8 a, REAL8 e, INT4 dir);
INT4 XLALPlusEquations(const gsl_vector * x, void *params, gsl_vector * f);
INT4 XLALCrossEquations(const gsl_vector * z, void *params, gsl_vector * g);
INT4 XLALIMRIGenerator(
	REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        REAL8 phi0,                     //Initial phi
        REAL8 dt,                       //Time step (in seconds)
        REAL8 m1,                       //BH mass (passed in kg)
        REAL8 m2,                       //CO mass (passed in kg)
        REAL8 f_min,                    //Initial frequency (passed in Hz)
        REAL8 r,                        //Distance to system (passed in meters)
        REAL8 inc,                      //Inclination angle between line of sight \hat n and BH spin axis
        REAL8 s1z,                      //BH reduced spin
	LALSimInspiralTestGRParam *testGR
        );
INT4 XLALIMRIstart(REAL8 m, REAL8 M, REAL8 a, REAL8 Dist, REAL8 Sdotn, REAL8 phi0, REAL8 psi0, REAL8 p0, REAL8 e0,
	REAL8Sequence *hplus, REAL8Sequence *hcross, REAL8 dt, INT4 dir, INT4 Npts);

//Assign an integer to each parameter label, to be used as indexes for params[]
enum {COMass, SMBHMass, SMBHSpin, DIR, DIST, COS_INC, THETA_NOUGHT, PHI_NOUGHT, PHI, ECC_ANOMALY, SEMILATUS, ECCENTRICITY, NParams};

//Data type to hold current simulation regime
enum stage {INSPIRAL, TRANSITION, PLUNGE};

INT4 XLALSearchArr(REAL8 *arr, REAL8 time) {

	INT4 i = 0;
	INT4 i_match;

	while (arr[i] < time) {
		i++;
		}

	if (fabs(arr[i-1] - time) < fabs(arr[i] - time))
		i_match = i-1;
	else
		i_match = i;

	return i_match;

	}

//Structure template used to match plunging waveform onto QNM ringdown
struct rparams {
		REAL8 a;
		REAL8 M;
		REAL8 Sdotn;
		REAL8 final_mass;
		REAL8 final_a;
		REAL8 w0;
		REAL8 w1;
		REAL8 w2;
		REAL8 wi0;
		REAL8 wi1;
		REAL8 wi2;
		REAL8 ab_factor;
		REAL8 hp_1;
		REAL8 dhpdi_1;
		REAL8 hp_2;
		REAL8 dhpdi_2;
		REAL8 hp_3;
		REAL8 dhpdi_3;
		REAL8 hx_1;
		REAL8 dhxdi_1;
		REAL8 hx_2;
		REAL8 dhxdi_2;
		REAL8 hx_3;
		REAL8 dhxdi_3;
		REAL8 dt;
		};

struct pParams {
	REAL8 m;
	REAL8 M;
	REAL8 e0;
	REAL8 a;
	REAL8 f0;
	INT4 dir; 
	};

REAL8 XLALInitialP(REAL8 p0, void *params) {

	struct pParams *p = (struct pParams *) params;
        REAL8 m = p->m;
        REAL8 M = p->M;
        REAL8 e0 = p->e0;
        REAL8 a = p->a;
	REAL8 f0 = p->f0;
        INT4 dir = p->dir;

	REAL8 E, L, nu;
	XLALComputeEnergyAng(p0,e0,a,dir,&E,&L);
        nu = (m*M)/pow(m+M,2.);

        return XLALDphiDt_sf(E,L,a,p0,e0,0.,nu) - pi*f0*(m+M)*Msun_sec;

        }

//============================================
// Energy and Ang. Momentum
// Following (gr-qc/0509101) A3-A5
//============================================

INT4 XLALComputeEnergyAng(REAL8 p, REAL8 e, REAL8 a, INT4 dir, REAL8 *Eng, REAL8 *Ang) {

	//=====================================================
	// Returns energy as defined in (gr-qc/0509101) Eq. A3
	//=====================================================

	//Equatorial Plane
	REAL8 theta_min = pi/2.;
	REAL8 zminus = pow(cos(theta_min),2.);

	//Get periastron and apastron radii
	REAL8 r1 = p/(1.+e);
	REAL8 r2 = p/(1.-e);

	//Determinant helper functions (evaluated at r1)
	REAL8 d1 = (r1*r1-2.*r1+a*a)*(r1*r1 + zminus*a*a);
	REAL8 f1 = pow(r1,4.) + a*a*(r1*r1 + 2.*r1 + zminus*(r1*r1-2.*r1+a*a));
	REAL8 g1 = 2.*a*r1;
	REAL8 h1 = r1*r1 - 2.*r1 + zminus*(r1*r1-2.*r1+a*a)/(1.-zminus);

	//Equations A12-A16
	REAL8 kappa;
	REAL8 epsilon;
	REAL8 rho;
	REAL8 eta;
	REAL8 sigma;
	REAL8 iota;

	if (e==0.) {

		//Helper function derivatives
		REAL8 dprime = 2.*(2.*r1-3.)*pow(r1,2.) + 2.*pow(a,2.)*((1.+zminus)*r1-zminus);
		REAL8 fprime = 4.*pow(r1,3.) + 2.*pow(a,2.)*((1.+zminus)*r1+(1.-zminus));
		REAL8 gprime = 2.*a;
		REAL8 hprime = 2.*(r1-1.)/(1.-zminus);

		//In the e=0 case, evaluation at r1 is the same as evaluating at p
		kappa	= d1*hprime - dprime*h1;
		epsilon = d1*gprime - dprime*g1;
		rho		= f1*hprime - fprime*h1;
		eta		= f1*gprime - fprime*g1;
		sigma	= g1*hprime - gprime*h1;
		iota	= f1*dprime - fprime*d1;
		
		}

	else {
	
		//Helper functions (evaluated at r2)
		REAL8 d2 = (r2*r2-2.*r2+a*a)*(r2*r2 + zminus*a*a);
		REAL8 f2 = pow(r2,4.) + a*a*(r2*r2 + 2.*r2 + zminus*(r2*r2-2.*r2+a*a));
		REAL8 g2 = 2.*a*r2;
		REAL8 h2 = r2*r2 - 2.*r2 + zminus*(r2*r2-2.*r2+a*a)/(1.-zminus);

		kappa 	= d1*h2 - d2*h1;
		epsilon = d1*g2 - d2*g1;
		rho 	= f1*h2 - f2*h1;
		eta 	= f1*g2 - f2*g1;
		sigma 	= g1*h2 - g2*h1;	
		iota	= f1*d2 - f2*d1;
		
		}

	//Define components of A3
	REAL8 Energy_numerator = kappa*rho + 2.*epsilon*sigma - 2.*dir*sqrt(sigma*(sigma*pow(epsilon,2.) + rho*epsilon*kappa - eta*pow(kappa,2.)));
	REAL8 Energy_denom = pow(rho,2.) + 4.*eta*sigma;
	REAL8 Ang_numerator = -(rho*iota-2.*eta*epsilon) - 2.*dir*sqrt(eta*(eta*pow(epsilon,2.) - rho*epsilon*iota - sigma*pow(iota,2.)));
	REAL8 Ang_denom = pow(rho,2.) + 4.*eta*sigma;
	
	*Eng = sqrt(Energy_numerator/Energy_denom);
	*Ang = sqrt(Ang_numerator/Ang_denom);

	return 0;

	}

//===========================
// Potential Derivatives
//===========================

REAL8 XLALDFDr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//========================================
	// Dimensionless dF/dr, where F=R/(V_t)^2
	//========================================

	REAL8 R = pow(E*(a*a+r*r)-a*Lz,2.) - (r*r-2.*r+a*a)*((Lz-a*E)*(Lz-a*E) + r*r);
	REAL8 Vt = a*(Lz-a*E) + ((r*r+a*a)/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);
	REAL8 dRdr = 4.*E*r*(E*(a*a+r*r)-a*Lz) - 2.*(r-1.)*(pow(Lz-a*E,2.)+r*r) - 2.*r*(r*r-2.*r+a*a);
	REAL8 dVtdr = 2.*r*(E*(r*r+a*a)-a*Lz)/(r*r-2.*r+a*a) + (r*r+a*a)*(2.*E*r)/(r*r-2.*r+a*a) - (r*r+a*a)*(E*(r*r+a*a)-a*Lz)*(2.*r-2.)/pow(r*r-2.*r+a*a,2.);

	return ( dRdr/pow(Vt,2.) - 2.*R*dVtdr/pow(Vt,3.) );
	}

REAL8 XLALDVphiVtDr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//=========================
	// Dimensionless d(Vphi/Vt)/dr
	//=========================

	REAL8 num = -2.*a*pow(Lz-E*a,2.)+6.*E*(Lz-E*a)*pow(r,2.)-2.*E*Lz*pow(r,3.);
	REAL8 denom = pow(2.*Lz*a - E*(pow(r,3.)+pow(a,2.)*(r+2.)),2.);

	return num/denom;
	}

REAL8 XLALD2FDrDE(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//======================================
	// Dimensionless d^2F/drdE, as given by Mathematica
	//======================================

	//Coefficients for various powers of r
	REAL8 R0;
	REAL8 R1;
	REAL8 R2;
	REAL8 R3;
	REAL8 R4;
	REAL8 R5;
	REAL8 R6;

	R0 = (-12.*Lz+pow(Lz,3.))*pow(a,3.) + (12.*E-E*Lz*Lz)*pow(a,4.) + Lz*pow(a,5.) + 2.*E*pow(a,6.);
	R1 = -8.*pow(Lz,3.)*a + 8.*E*Lz*Lz*a*a + 4.*Lz*pow(a,3.) + (-4.*E+E*Lz*Lz)*pow(a,4.) + E*pow(a,6.);
	R2 = 7.*pow(Lz,3.)*a - 6.*E*Lz*Lz*a*a + 10.*Lz*pow(a,3.) - 3.*E*pow(a,4.);
	R3 = -12.*Lz*a + 2.*E*Lz*Lz*a*a + 2.*E*pow(a,4.);
	R4 = -5.*E*Lz*Lz + 9.*Lz*a - 8.*E*a*a;
	R5 = E*Lz*Lz + E*a*a;
	R6 = -3.*E;

	return -4.*pow(r*r - 2.*r + a*a,2.)*(R0 + R1*r + R2*pow(r,2.) + R3*pow(r,3.) + R4*pow(r,4.) + R5*pow(r,5.) + R6*pow(r,6.))
			/pow(-2.*Lz*a + 2.*E*a*a + E*a*a*r + E*r*r*r,4.);
	}

REAL8 XLALD2FDrDL(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//=========================================
	// Dimensionless d^2F/drdL, as given by Mathematica
	//=========================================

	//Coefficients for various powers of r
	REAL8 R0;
	REAL8 R1;
	REAL8 R2;
	REAL8 R3;
	REAL8 R4;
	REAL8 R5;

	R0 = -12.*Lz*pow(a,2.) + (12.*E+E*Lz*Lz)*pow(a,3.) - E*E*Lz*pow(a,4.) + 3.*E*pow(a,5.);
	R1 = -8.*E*Lz*Lz*a + (12.*Lz+8.*E*E*Lz)*pow(a,2.) - 12.*E*pow(a,3.) + E*E*Lz*pow(a,4.);
	R2 = 7.*E*Lz*Lz*a - 6.*E*E*Lz*pow(a,2.) + 6.*E*pow(a,3.);
	R3 = -12.*E*a + 2.*E*E*Lz*pow(a,2.);
	R4 = -5.*E*E*Lz + 3.*E*a;
	R5 = E*E*Lz;

	return 4.*pow(r*r-2.*r+a*a,2.)*(R0 + R1*r + R2*pow(r,2.) + R3*pow(r,3.) + R4*pow(r,4.) + R5*pow(r,5.))
			/pow(-2.*Lz*a + 2.*E*a*a + E*a*a*r + E*r*r*r,4.);
	}

//=====================================
// Energy and Angular Momentum Fluxes
//=====================================

REAL8 XLALEdot3PN_Effective(REAL8 a, REAL8 p, REAL8 e, REAL8 nu) {

	//=================================================
	// dE/dt for equatorial, elliptic orbits.
	//=================================================

	//Define enhancement coefficients
	REAL8 phi2,phi4,phi6,phi8,phi10;
	REAL8 psi2,psi4,psi6,psi8,psi10;
	REAL8 zet2,zet4,zet6,zet8,zet10;
	REAL8 kap2,kap4,kap6,kap8,kap10;

	phi2 = 991./192.;
	phi4 = -2099./256.;
	phi6 = -2.82574;
	phi8 = 4.61120;
	phi10 = 0.251162;

	psi2 = -80325./8191.;
	psi4 = -15201315./524224.;
	psi6 = 24.6550;
	psi8 = 12.5662;
	psi10 = 0.615425;

	zet2 = 668761./48972.;
	zet4 = 3244111./261184.;
	zet6 = -19.1172;
	zet8 = -7.08972;
	zet10 = -0.938933;

	kap2 = (76595216.-73821440.*log(2.)+73712835.*log(3.))/5604528.;
	kap4 = 7.*(63167701.+916168240.*log(2.)-494929035.*log(3.))/22418112.;
	kap6 = -6.44950;
	kap8 = 72.9269;
	kap10 = -87.8989;

	//Flux coefficients
	REAL8 g3A,g4A,g5A,g6A,g7A,g8A;			//Higher order e-dependence of O(\nu) radiative components
	REAL8 g1b,g2b,g3b,g4b,g5b,g6b,g7b,g8b;	//Second order O(\nu^2) radiative components
	REAL8 gEA,gEb;							//3PN first and second-order radiative terms
	REAL8 gElogA,gElogb;						//Logarithmic dependence at first and second order
	REAL8 g1,g2,g3,g4,g5,g6,g7,g8,gE,gElog;	//Final combined PN coefficients

	//Define higher eccentricity corrections
	g3A = -(809./128.)*pow(e,4.) -(8609./5376)*pow(e,6.);

	g4A = (1./(48.*pow(1.-e*e,2.)))*pow(e,2.)*(
		- 991.+192.*phi2
		+ pow(e,2.)*(2558.+192.*phi4)
		+ pow(e,4.)*(-1375.+192.*phi6)
		+ pow(e,6.)*192.*phi8
		+ pow(e,8.)*192.*phi10);

	g5A = (1./193536.)*(-1451520.
		- 2237760.*pow(e,2.)
		+ 27633560.*pow(e,4.)
		- 1221612.*pow(e,6.)
		- 1114881.*pow(e,8.)
		+ sqrt(1.-pow(e,2.))*(1451520.
		+ 2963520.*pow(e,2.)
		- 3855600.*pow(e,4.)
		- 559440.*pow(e,6.)));

	g6A = pow(e,4.)*(1465./128.) + pow(e,6.)*(883./768.);

	g7A = (pow(e,2.)/(672.*pow(1.-e*e,2.)))*(
		- 2919. + 16128.*phi2 + 8191.*psi2
		+ pow(e,2.)*(169933. + 61824.*phi2 + 32256.*phi4 - 8191.*psi2 + 8191.*psi4)
		+ pow(e,4.)*(-89062. + 45696.*phi4 + 48384.*phi6 - 8191.*psi4 + 8191.*psi6)
		+ pow(e,6.)*(29568.*phi6 + 64512.*phi8 - 8191.*(psi6 - psi8))
		+ pow(e,8.)*(80640.*phi10 + 13440.*phi8 + 8191.*(psi10 - psi8))
		- pow(e,10.)*(2688.*phi10 + 8191.*psi10));

	g8A = 0;

	//Define higher mass-ratio corrections
	g1b = 0.;

	g2b = (1./64.)*(-3432.
		- pow(e,2.)*9544.
		- pow(e,4.)*3411.
		- pow(e,6.)*98.);

	g3b = (1./768.)*(664.
		+ pow(e,2.)*12008.
		+ pow(e,4.)*9387.
		+ pow(e,6.)*481.);

	g4b = 0.;

	g5b = (1./129024.)*(2161088.
		- pow(e,2.)*4462736.
		- pow(e,4.)*5449464.
		+ pow(e,6.)*4260315.
		+ pow(e,8.)*358689.
		+ sqrt(1.-e*e)*(
			- 387072.
			- pow(e,2.)*790272.
			+ pow(e,4.)*1028160.
			+ pow(e,6.)*149184.));
			
	g6b = (1./384.)*(-10448.
		- pow(e,2.)*37960.
		- pow(e,4.)*18330.
		- pow(e,6.)*735.);

	g7b = (1./(24.*pow(-1 + e*e,2)))*(363.
		- 1707.*pow(e,2.)
		+ 583.*pow(e,10.)*zet10
		- 583.*pow(e,12.)*zet10
		+ 583.*pow(e,2.)*zet2
		- 583.*pow(e,4.)*zet2
		+ 583.*pow(e,4.)*zet4
		- 583.*pow(e,6.)*zet4
		+ 583.*pow(e,6.)*zet6
		- 583.*pow(e,8.)*zet6
		+ 583.*pow(e,8.)*zet8
		- 583.*pow(e,10.)*zet8
		- 120.*pow(e,8.)*phi10
		- 940.*pow(e,10.)*phi10
		- 284.*pow(e,12.)*phi10
		- 24.*phi2
		- 364.*pow(e,2.)*phi2
		- 956.*pow(e,4.)*phi2
		- 48.*pow(e,2.)*phi4
		- 508.*pow(e,4.)*phi4
		- 788.*pow(e,6.)*phi4
		- 72.*pow(e,4.)*phi6
		- 652.*pow(e,6.)*phi6
		- 620.*pow(e,8.)*phi6
		- 96.*pow(e,6.)*phi8
		- 796.*pow(e,8.)*phi8
		- 452.*pow(e,10.)*phi8);

	g8b = 0.;

	//Define 3PN terms
	gEA = (1./(17882726400.*sqrt(1.-e*e)))*(
		- 256.*(4424511234. - 11068250753.*sqrt(1.-e*e) + 10644480.*sqrt(1-e*e)*(107.*gamma-35.*pow(pi,2.)))
		- 128.*pow(e,2.)*(42116770542. - 319720104403.*sqrt(1.-e*e) + 301593600.*sqrt(1.-e*e)*(107.*gamma - 35.*pow(pi,2)) + 4438786176.*kap2)
		- 64.*pow(e,4.)*(57526950789. - 1178152279409.*sqrt(1.-e*e) + 1146720960.*sqrt(1.-e*e)*(107.*gamma - 35.*pow(pi,2.)) + 8877572352.*kap4)
		- 96.*pow(e,6.)*(-7.*(15484930967. + 23518217253.*sqrt(1.-e*e)) + 258867840.*sqrt(1.-e*e)*(107.*gamma - 35.*pow(pi,2.)) + 5918381568.*kap6)
		- 6.*pow(e,8.)*(-113485018521.*sqrt(1.-e*e) + 131725440.*sqrt(1.-e*e)*(107.*gamma - 35.*pow(pi,2.)) + 45056.*(1475341. + 2101698.*kap8))
		+ 9.*pow(e,10.)*(37078693925.*sqrt(1.-e*e) - 1056.*(38915625. + 59781632.*kap10))
		- 47456640.*sqrt(1.-e*e)*(3072. + 43520.*pow(e,2.) + 82736.*pow(e,4.) + 28016.*pow(e,6.) + 891.*pow(e,8.))*(
			log(64.) + 3.*log(1.-e*e) - 2.*log(1.+sqrt(1.-e*e))));

	gEb = (1./(4644864.*sqrt(1. - pow(e,2.))))*(231451776.
		+ 717700032.*pow(e,2.)
		- 1346174496.*pow(e,4.)
		- 320607360.*pow(e,6.)
		+ 651679200.*pow(e,8.)
		+ 65950848.*pow(e,10.)
		- 2521083152.*sqrt(1. - pow(e,2.))
		- 11007222480.*pow(e,2.)*sqrt(1. - pow(e,2.))
		- 2502051186.*pow(e,4.)*sqrt(1. - pow(e,2.))
		+ 3018130094.*pow(e,6.)*sqrt(1. - pow(e,2.))
		- 282120903.*pow(e,8.)*sqrt(1. - pow(e,2.))
		- 43482825.*pow(e,10.)*sqrt(1. - pow(e,2.))
		- 2975616.*pow(pi,2.)
		- 3099600.*pow(e,2.)*pow(pi,2.)
		+ 13979196.*pow(e,4.)*pow(pi,2.)
		- 6757128.*pow(e,6.)*pow(pi,2.)
		- 1146852.*pow(e,8.)*pow(pi,2.)
		+ 95467680.*sqrt(1. - pow(e,2.))*pow(pi,2.)
		+ 377283312.*pow(e,2.)*sqrt(1. - pow(e,2.))*pow(pi,2.)
		+ 85208004.*pow(e,4.)*sqrt(1. - pow(e,2.))*pow(pi,2.)
		- 79768206.*pow(e,6.)*sqrt(1. - pow(e,2.))*pow(pi,2.)
		- 4602906.*pow(e,8.)*sqrt(1. - pow(e,2.))*pow(pi,2.));

	//Define log(v) terms
	gElogA = (-107./20160.)*(3072.
		+ pow(e,2.)*43520.
		+ pow(e,4.)*82736.
		+ pow(e,6.)*28016.
		+ pow(e,8.)*891.);

	gElogb = 0.;

	//Define coefficients
	g1 = 1. + (73./24.)*pow(e,2.) + (37./96.)*pow(e,4.) + nu*g1b;
	g2 = (73./12.) + (823./24.)*pow(e,2.) + (949./32.)*pow(e,4.) + (491./192.)*pow(e,6.) + nu*g2b;
	g3 = (1247./336.) + (9181./672.)*pow(e,2.) + g3A + nu*g3b;
	g4 = 4. + (1375./48.)*pow(e,2.) + g4A + nu*g4b;
	g5 = (44711./9072.) + (172157./2592.)*pow(e,2.) + g5A + nu*g5b;
	g6 = (33./16.) + (359./32.)*pow(e,2.) + g6A + nu*g6b;
	g7 = (8191./672.) + (44531./336.)*pow(e,2.) + g7A + nu*g7b;
	g8 = (3749./336.) - (5143./168.)*pow(e,2.) + g8A + nu*g8b;
	gE = gEA + nu*gEb;
	gElog = gElogA + nu*gElogb;
	
	//Analogous flux coefficients evaluated at e=0. All g3A_0 = 0.
	REAL8 g1b_0,g2b_0,g3b_0,g4b_0,g5b_0,g6b_0,g7b_0,g8b_0;
	REAL8 gEA_0,gEb_0;
	REAL8 gElogA_0,gElogb_0;
	REAL8 g1_0,g2_0,g3_0,g4_0,g5_0,g6_0,g7_0,g8_0,gE_0,gElog_0;
	
	g1b_0 = 0.;
	g2b_0 = (1./64.)*(-3432.);
	g3b_0 = (1./768.)*(664.);	
	g4b_0 = 0.;	
	g5b_0 = (1./129024.)*(2161088. + sqrt(1.)*(- 387072.));	
	g6b_0 = (1./384.)*(-10448.);	
	g7b_0 = (1./24.)*(363. - 24.*phi2);	
	g8b_0 = 0.;
	gEA_0 = (6643739519./69854400.)-(1712./105.)*gamma+(16./3.)*pow(pi,2.)-(1712./105.)*log(4.);
	gEb_0 = (-143101961.+5780754*pow(pi,2.))/290304.;
	gElogA_0 = (-1712./105.);
	gElogb_0 = 0.;
	
	
	//JUST GO AHEAD AND COMBINE THESE (with the above) PROBABLY....
	g1_0 = 1. + nu*g1b_0;
	g2_0 = (73./12.) + nu*g2b_0;
	g3_0 = (1247./336.) + nu*g3b_0;
	g4_0 = 4. + nu*g4b_0;
	g5_0 = (44711./9072.) + nu*g5b_0;
	g6_0 = (33./16.) + nu*g6b_0;
	g7_0 = (8191./672.) + nu*g7b_0;
	g8_0 = (3749./336.) + nu*g8b_0;
	gE_0 = gEA_0 + nu*gEb_0;
	gElog_0 = gElogA_0 + nu*gElogb_0;

	REAL8 eccentricTerms;
	REAL8 circularTerms;

	eccentricTerms =  g1 - g3/p + (pi*g4-a*g2)/pow(p,3./2.) + (a*a*g6-g5)/pow(p,2.) + (a*g8-pi*g7)/pow(p,5./2.) 
		+ gE/pow(p,3.) + gElog*log(1./sqrt(p))/pow(p,3.);

	circularTerms =  g1_0 - g3_0/p + (pi*g4_0-a*g2_0)/pow(p,3./2.) + (a*a*g6_0-g5_0)/pow(p,2.) + (a*g8_0-pi*g7_0)/pow(p,5./2.) 		+ gE_0/pow(p,3.) + gElog_0*log(1./sqrt(p))/pow(p,3.);

	return -(32./5.)*nu*pow(p,-5.)*(eccentricTerms-circularTerms);

	}

REAL8 XLALLdot3PN_Effective(REAL8 a, REAL8 p, REAL8 e, REAL8 nu, INT4 dir) {

	//=================================================
	// dL/dt for equatorial, elliptic orbits.
	//=================================================

	//Define enhancement coefficients
	REAL8 phiSg2,phiSg4,phiSg6,phiSg8,phiSg10;
	REAL8 psiSg2,psiSg4,psiSg6,psiSg8,psiSg10;
	REAL8 zetSg2,zetSg4,zetSg6,zetSg8,zetSg10;
	REAL8 kapSg2,kapSg4,kapSg6,kapSg8,kapSg10;

	phiSg2 = -15./32.;
	phiSg4 = -749./128.;
	phiSg6 = 9.85010;
	phiSg8 = -5.15135;
	phiSg10 = 0.610049;

	psiSg2 = -74753./8191.;
	psiSg4 = 4611875./524224.;
	psiSg6 = 11.3765;
	psiSg8 = -16.7590;
	psiSg10 = 4.58258;

	zetSg2 = 45237./8162.;
	zetSg4 = -3195515./261184.;
	zetSg6 = -0.296213;
	zetSg8 = 9.39768;
	zetSg10 = -3.34773;

	kapSg2 = -(15./3736352.)*(-1284371. + 2193072.*log(2.) - 2184084.*log(3.));
	kapSg4 = (35./7472704.)*(-1751415. + 36083824.*log(2.) - 21996846.*log(3.));
	kapSg6 = -12.8413;
	kapSg8 = 1.40457;
	kapSg10 = 5.82671;

	//Low order expansion of fitted exponential

	REAL8 g11A,g12A,g13A,g14A,g15A,g16A;
	REAL8 g9b,g10b,g11b,g12b,g13b,g14b,g15b,g16b;
	REAL8 gLA,gLb;
	REAL8 gLlogA,gLlogb;
	REAL8 g9,g10,g11,g12,g13,g14,g15,g16,gL,gLlog;

	//Define higher eccentricity corrections (cannot get g14A and higher from matching)
	g11A = -(10751./2688.)*pow(e,4.);

	g12A = (1./(8.*pow(1.-e*e,4.)))*(-32.
		+ pow(e,2.)*31.
		+ pow(e,4.)*196.
		- pow(e,6.)*454.
		+ pow(e,8.)*356.
		- pow(e,10.)*97.
		+ sqrt(1.-pow(e,2.))*(32.
		+ pow(e,2.)*32.*phiSg2
		+ pow(e,4.)*32.*phiSg4
		+ pow(e,6.)*32.*phiSg6
		+ pow(e,8.)*32.*phiSg8
		+ pow(e,10.)*32.*phiSg10));

	g13A = (1./48384.)*(-362880.
		+ pow(e,2.)*226800.
		+ pow(e,4.)*1743550.
		- pow(e,6.)*629733.
		+ sqrt(1.-pow(e,2.))*(362880.
		- pow(e,2.)*45360.
		- pow(e,4.)*317520.));

	g14A = (311./128.)*pow(e,4.);

	g15A = (-1./(1344.*pow(1.-pow(e,2.),7./2.)))*(-16382.
		+ pow(e,2.)*(-155650. - 32256.*phiSg2 - 16382.*psiSg2)
		+ pow(e,4.)*(-139776.*phiSg2 - 64512.*phiSg4 + 16382.*psiSg2 - 16382.*psiSg4)
		+ pow(e,6.)*(-107520.*phiSg4 - 96768.*phiSg6 + 16382.*psiSg4 - 16382.*psiSg6)
		+ pow(e,8.)*(-75264.*phiSg6 - 129024.*phiSg8 + 16382.*psiSg6 - 16382.*psiSg8)
		+ pow(e,10.)*(-161280.*phiSg10 - 43008.*phiSg8 - 16382.*psiSg10 + 16382.*psiSg8)
		+ pow(e,12.)*(-10752.*phiSg10 + 16382.*psiSg10)
		+ pow(1.-e*e,7./2.)*(16382. + pow(e,2.)*48361.));

	g16A = 0.;

	//Define higher mass ratio corrections (missing spin-dependent g10b, cannot get g14b and higher)
	g9b = 0.;

	g10b = (1./96.)*(-2908. - 3716.*pow(e,2.) - 401.*pow(e,4.));

	g11b = (1./192.)*(366. + 1679*pow(e,2.) + 305.*pow(e,4.));

	g12b = 0.;

	g13b = (1./64512.)*(129208.
		- 2517684.*pow(e,2.)
		+ 1092843.*pow(e,4.)
		+ 476577.*pow(e,6.)
		- 193536.*sqrt(1. - pow(e,2.))
		+ 24192.*pow(e,2.)*sqrt(1. - pow(e,2.))
		+ 169344.*pow(e,4.)*sqrt(1. - pow(e,2.)));

	g14b = (-15./64.)*(72. + 128.*pow(e,2.) + 17.*pow(e,4.));

	g15b = (-1./(24.*pow(1.-e*e,7./2.)))*(-375.
		+ 1719.*pow(e,2.)
		- 583.*pow(e,10.)*zetSg10
		+ 583.*pow(e,12.)*zetSg10
		- 583.*pow(e,2.)*zetSg2
		+ 583.*pow(e,4.)*zetSg2
		- 583.*pow(e,4.)*zetSg4
		+ 583.*pow(e,6.)*zetSg4
		- 583.*pow(e,6.)*zetSg6
		+ 583.*pow(e,8.)*zetSg6
		- 583.*pow(e,8.)*zetSg8
		+ 583.*pow(e,10.)*zetSg8
		+ 120.*pow(e,8.)*phiSg10
		+ 928.*pow(e,10.)*phiSg10
		+ 296.*pow(e,12.)*phiSg10
		+ 24.*phiSg2
		+ 352.*pow(e,2.)*phiSg2
		+ 968.*pow(e,4.)*phiSg2
		+ 48.*pow(e,2.)*phiSg4
		+ 496.*pow(e,4.)*phiSg4
		+ 800.*pow(e,6.)*phiSg4
		+ 72.*pow(e,4.)*phiSg6
		+ 640.*pow(e,6.)*phiSg6
		+ 632.*pow(e,8.)*phiSg6
		+ 96.*pow(e,6.)*phiSg8
		+ 784.*pow(e,8.)*phiSg8
		+ 464*pow(e,10)*phiSg8);

	g16b = 0.;

	//Define 3PN terms
	gLA = (1./(12169195315200.*pow(1.-e*e,5./2.)))*(
		- 174208.*(1867780530. - 8511520049.*sqrt(1.-e*e) + 10644480.*sqrt(1.-e*e)*(107.*gamma - 35.*pow(pi,2.)))
		- 174208.*pow(e,2.)*(6850686150. - 58670843875.*sqrt(1.-e*e) + 9504.*sqrt(1.-e*e)*(5775.*(107.*gamma - 35.*pow(pi,2.)) + 233522.*kapSg2))
		+ 32.*pow(e,4.)*(148043381495625. - 590675694967855.*sqrt(1.-e*e) + 25869888.*sqrt(1.-e*e)*(18375.*(107.*gamma - 35.*pow(pi,2.)) - 467044.*kapSg4))
		+ 16.*pow(e,6.)*(-175.*(1140375922488. + 622597756831.*sqrt(1.-e*e)) + 12934944.*sqrt(1.-e*e)*(25025.*(107.*gamma - 35.*pow(pi,2.)) - 1868176.*kapSg6))
		- 9.*pow(e,8.)*(175.*(1288884769216. - 6632090872287.*sqrt(1.-e*e)) + 45990912.*sqrt(1.-e*e)*(20475.*(107.*gamma - 35.*pow(pi,2.)) + 934088.*kapSg8))
		- 8166.*pow(e,10.)*(7.*(-45644089920. + 41564211853.*sqrt(1.-e*e)) + 25344.*sqrt(1.-e*e)*(258405.*gamma - 84525.*pow(pi,2.) + 1868176.*kapSg10))
		+ 2143575.*pow(e,12.)*(-282062880. + 226213231.*sqrt(1.-e*e))
		- 387530922240.*pow(1.-e*e,5./2.)*(256. + 1832.*pow(e,2.) + 1308.*pow(e,4.) + 69.*pow(e,6.))*(log(64.) + 3.*log(1.-e*e) - 2.*log(1. + sqrt(1.-e*e))));

	gLb = (1./1161216.)*(-415871908.
		- 836752274.*pow(e,2.)
		+ 353749248.*pow(e,4.)
		+ 20730867.*pow(e,6.)
		- 27872577.*pow(e,8.)
		+ 66934944.*sqrt(1. - pow(e,2.))
		+ 74324736.*pow(e,2.)*sqrt(1. - pow(e,2.))
		- 101108304.*pow(e,4.)*sqrt(1. - pow(e,2.))
		- 40151376.*pow(e,6.)*sqrt(1. - pow(e,2.))
		+ 14692104.*pow(pi,2.)
		+ 28826280.*pow(e,2.)*pow(pi,2.)
		- 348705.*pow(e,4.)*pow(pi,2.)
		- 697410.*pow(e,6.)*pow(pi,2.)
		- 743904.*sqrt(1. - pow(e,2.))*pow(pi,2.)
		+ 92988.*pow(e,2.)*sqrt(1. - pow(e,2.))*pow(pi,2.)
		+ 650916.*pow(e,4.)*sqrt(1. - pow(e,2.))*pow(pi,2.));

	//Define log(v) terms
	gLlogA = (-107./1680.)*(256.
		+ pow(e,2.)*1832.
		+ pow(e,4.)*1308.
		+ pow(e,6.)*69.);

	gLlogb = 0.;

	//Define coefficients
	g9 = 1. + (7./8.)*pow(e,2.) + nu*g9b;
	g10 = (61./12.) + (119./8.)*pow(e,2.) + (183./32.)*pow(e,4.) + nu*g10b;
	g11 = (1247./336.) + (425./336.)*pow(e,2.) + g11A + nu*g11b;
	g12 = 4. + (97./8.)*pow(e,2.) + g12A + nu*g12b;
	g13 = (44711./9072.) + (302893./6048.)*pow(e,2.) + g13A + nu*g13b;
	g14 = (33./16.) + (95./16.)*pow(e,2.) + g14A + nu*g14b;
	g15 = (8191./672.) + (48361./1344.)*pow(e,2.) + g15A + nu*g15b;
	g16 = (417./56.) - (37241./672.)*pow(e,2.) + g16A + nu*g16b;
	gL = gLA + nu*gLb;
	gLlog = gLlogA + nu*gLlogb;

	//Evaluate the above coefficients (truncated to 2PN) at e=0.
	REAL8 g9b_0,g10b_0,g11b_0,g12b_0,g13b_0,g14b_0;
	REAL8 g9_0,g10_0,g11_0,g12_0,g13_0,g14_0;
	REAL8 gLA_0, gLb_0, gLlogA_0, gLlogb_0, gL_0, gLlog_0;

	g9b_0 = 0.;
	g10b_0 = (1./96.)*(-2908.);
	g11b_0 = (1./192.)*(366.);
	g12b_0 = 0.;
	g13b_0 = (1./64512.)*(-64328.);
	g14b_0 = (-15./64.)*(72.);

	gLA_0 = (6643739519./69854400.) - (1712./105.)*gamma + (16./3.)*pow(pi,2.) - (1712./105.)*log(4.);
	gLb_0 = -(87234241./290304.) + (3075./256.)*pow(pi,2.);
	gLlogA_0 = (-1712./105.);
	gLlogb_0 = 0.;

	g9_0 = 1. + nu*g9b_0;
	g10_0 = (61./12.) + nu*g10b_0;
	g11_0 = (1247./336.) + nu*g11b_0;
	g12_0 = 4. + nu*g12b_0;
	g13_0 = (44711./9072.) + nu*g13b_0;
	g14_0 = (33./16.) + nu*g14b_0;

	gL_0 = gLA_0 + nu*gLb_0;
	gLlog_0 = gLlogA_0 + nu*gLlogb_0;
	

	REAL8 eccentricTerms;
	REAL8 circularTerms;

	eccentricTerms = g9 - g11/p + (pi*g12-a*g10)/pow(p,3./2.) + (a*a*g14-g13)/pow(p,2.) + (a*g16-pi*g15)/pow(p,5./2.)
		+ (gL)/pow(p,3.) + (gLlog)*log(1./sqrt(p))/pow(p,3.);

	circularTerms = g9_0 - g11_0/p + (pi*g12_0-a*g10_0)/pow(p,3./2.) + (a*a*g14_0-g13_0)/pow(p,2.);

	return -dir*(32./5.)*nu*pow(p,-7./2.)*(eccentricTerms-circularTerms);

	}

REAL8 XLALEnergyAngDot(REAL8 a, REAL8 p, REAL8 e, REAL8 nu, INT4 dir, REAL8* dEdt, REAL8* dLdt) {

	//Equatorial
	REAL8 iota = 0;

	//Define Teukolsky fit corrections. Corrections taken from hLISAstart, which give more decimal places than Gair&Glamp
	REAL8 da1, db1, dc1;
	REAL8 da2, db2, dc2;
	REAL8 ca1, cb1, cc1;
	REAL8 ca2, cb2, cc2;
	REAL8 ca3, cb3, cc3;
	REAL8 ca4, cb4, cc4;
	REAL8 ca5, cb5, cc5;
	REAL8 ca6, cb6, cc6;
	REAL8 ca7, cb7, cc7;
	REAL8 ca8, cb8, cc8;
	REAL8 ca9, cb9, cc9;
	REAL8 fa1, fb1;
	REAL8 fa2, fb2;
	REAL8 fa3, fb3;
	REAL8 fa4, fb4;
	REAL8 fa5, fb5;
	REAL8 fa6, fb6;

	da1 = -10.741956;
	db1 = 28.5942157;
	dc1 = -9.077378144;
	da2=-1.428362761;
	db2=10.70029768;
	dc2=-33.70903016;

	ca1=-28.15174147;
	cb1=60.9607071973;
	cc1=40.99984205;
	ca2=-0.348161211;
	cb2=2.37258476;
	cc2=-66.65840948;
	ca3=-0.715392387;
	cb3=3.21592568;
	cc3=5.28887649;
	ca4=-7.6103411; 
	cb4=128.87778309;
	cc4=-475.4650442;
	ca5=12.290783385;
	cb5=-113.1250548;
	cc5=306.11883292;
	ca6=40.9258725;
	cb6=-347.2713496;       
	cc6=886.50332051;
	ca7=-25.48313727;
	cb7=224.22721861;
	cc7=-490.98212316;
	ca8=-9.006337706;
	cb8=91.17666278;
	cc8=-297.001939215;
	ca9=-0.64500047;
	cb9=-5.13591989;
	cc9=47.19818628;

	fb1 = 736.2086781;
	fa1 = -283.9553066;
	fb2 = -1325.1852209;
	fa2 = 483.266206498;
	fb3 = 634.49936445;
	fa3 = -219.223848944;
	fb4 = 82.07804475;
	fa4 = -25.82025864;
	fb5 = -904.16109275;
	fa5 = 301.477789146;
	fb6 = 827.31891826;
	fa6 = -271.9659423;

	//Teukolsky fit Ldot
	REAL8 Ldot_fit;
	REAL8 c1, c2, c3, c4, c5;
	REAL8 c6, c7, c8, c9, c10;
	REAL8 c11, c12, c13, c14, c15;
	REAL8 c16, c17;

	c1 = da1 + db1*pow(p,-1./2.) + dc1*pow(p,-1.);
	c2 = da2 + db2*pow(p,-1./2.) + dc2*pow(p,-1.);
	c3 = ca1 + cb1*pow(p,-1./2.) + cc1*pow(p,-1.);
	c4 = ca2 + cb2*pow(p,-1./2.) + cc2*pow(p,-1.);
	c5 = ca3 + cb3*pow(p,-1./2.) + cc3*pow(p,-1.);
	c6 = ca4 + cb4*pow(p,-1./2.) + cc4*pow(p,-1.);
	c7 = ca5 + cb5*pow(p,-1./2.) + cc5*pow(p,-1.);
	c8 = ca6 + cb6*pow(p,-1./2.) + cc6*pow(p,-1.);
	c9 = ca7 + cb7*pow(p,-1./2.) + cc7*pow(p,-1.);
	c10 = ca8 + cb8*pow(p,-1./2.) + cc8*pow(p,-1.);
	c11 = ca9 + cb9*pow(p,-1./2.) + cc9*pow(p,-1.);
	c12 = fa1 + fb1*pow(p,-1./2.);
	c13 = fa2 + fb2*pow(p,-1./2.);
	c14 = fa3 + fb3*pow(p,-1./2.);
	c15 = fa4 + fb4*pow(p,-1./2.);
	c16 = fa5 + fb5*pow(p,-1./2.);
	c17 = fa6 + fb6*pow(p,-1./2.);

	Ldot_fit = -dir*(32./5.)*nu*pow(p,-7./2.)*(cos(iota) + a*pow(p,-3./2.)*(61./24.-(61./8.)*pow(cos(iota),2.))
		+ (-1247./336.)*pow(p,-1.)*cos(iota) + 4.*pi*pow(p,-3./2.)*cos(iota)
		+ (-44711./9072.)*pow(p,-2.)*cos(iota) + a*a*pow(p,-2.)*cos(iota)*(33./16.-(45./8.)*pow(sin(iota),2.))
		+ pow(p,-5./2.)*(
			c1*a +
			c2*pow(a,3.) + 
			c3*cos(iota) +
			c4*pow(a,2.)*cos(iota) +
			c5*pow(a,4.)*cos(iota) +
			c6*a*pow(cos(iota),2.) +
			c7*pow(a,3.)*pow(cos(iota),2.) +
			c8*pow(a,2.)*pow(cos(iota),3.) +
			c9*pow(a,4.)*pow(cos(iota),3.) +
			c10*pow(a,3.)*pow(cos(iota),4.) +
			c11*pow(a,4.)*pow(cos(iota),5.))
		+ pow(p,-7./2.)*a*cos(iota)*(
			c12 +
			c13*a +
			c14*pow(a,2.) +
			c15*pow(cos(iota),2.) +
			c16*a*pow(cos(iota),2.) +
			c17*pow(a,2.)*pow(cos(iota),2.)));

	//Replace the 2PN circular (e=0) component of Ldot3PN with Ldot_fit
	REAL8 dLdt_0 = (XLALLdot3PN_Effective(a,p,0.,nu,dir) + Ldot_fit);
	*dLdt = pow(1.-e*e,3./2.)*(XLALLdot3PN_Effective(a,p,e,nu,dir) + Ldot_fit);

	//Evaluate energy and angular momentum for e=0
	//REAL8 E0 	= Energy(p,0.,a,dir);
	//REAL8 L0 	= Ang(p,0.,a,dir);
	REAL8 E0, L0;
	
	XLALComputeEnergyAng(p,0.,a,dir,&E0,&L0);

	REAL8 N1 = E0*pow(p,4.) + pow(a,2.)*E0*pow(p,2.) - 2.*a*(L0-a*E0)*p;
	REAL8 N2 = -L0*pow(p,2.) + 2.*(L0-a*E0)*p;

	*dEdt =  pow(1.-e*e,3./2.)*(XLALEdot3PN_Effective(a,p,e,nu) - (N2/N1)*dLdt_0);

	return 0;

	}
	
INT4 XLALPEdot(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 nu, REAL8 dir, REAL8 *pdot, REAL8 *edot) {
	
	//==================================================================
	// Computes time derivatives of semilatus p and eccentricity e.
	// Resulting values are stored in the input variables pdot and edot.
	//==================================================================
	
	//Compute angular momentum derivative
	//REAL8 dLdt = Ldot(a,p,e,nu,dir);
	REAL8 dLdt, dEdt;
	XLALEnergyAngDot(a,p,e,nu,dir,&dEdt,&dLdt);	

	//Circular case
	if (e==0.) {
			
		REAL8 denom = pow(p,3.) - 3.*pow(p,2.) + 2.*dir*a*pow(p,3./2.);
		REAL8 dLdp = (2.*p - dir*a/sqrt(p))/sqrt(denom) + (-1./2.)*(pow(p,2.)-2.*dir*a*sqrt(p)+pow(a,2.))*(3.*pow(p,2.)-6.*p+3.*dir*a*sqrt(p))/pow(denom,3./2.);
		
		*pdot = dir*(1./dLdp)*dLdt;
		*edot = 0.;
		
		}
	
	//Eccentric case	
	else {
	
		//REAL8 dEdt = Edot(a,p,e,nu,dir);
		REAL8 N1_a = E*pow(p/(1.-e),4.) + pow(a,2.)*E*pow(p/(1.-e),2.) - 2.*a*(Lz-a*E)*p/(1.-e);
		REAL8 N1_p = E*pow(p/(1.+e),4.) + pow(a,2.)*E*pow(p/(1.+e),2.) - 2.*a*(Lz-a*E)*p/(1.+e);
		REAL8 N2_a = -Lz*pow(p/(1.-e),2.) + 2.*(Lz-a*E)*p/(1.-e);
		REAL8 N2_p = -Lz*pow(p/(1.+e),2.) + 2.*(Lz-a*E)*p/(1.+e);
		REAL8 D_p = -2.*(1.-E*E)*pow(p/(1.+e),3.) + 3.*pow(p/(1.+e),2.) - (a*a*(1.-E*E)+pow(Lz,2.))*p/(1.+e) + pow(Lz-a*E,2.);
		REAL8 D_a = -2.*(1.-E*E)*pow(p/(1.-e),3.) + 3.*pow(p/(1.-e),2.) - (a*a*(1.-E*E)+pow(Lz,2.))*p/(1.-e) + pow(Lz-a*E,2.);

		//Of the form dp/dt = (dp/dE)*(dE/dt) + (dp/dL)*(dL/dt)
		*pdot = (-(1./2.)*(pow(1.+e,2.)*N1_p/D_p + pow(1.-e,2.)*N1_a/D_a))*dEdt + (-(1./2.)*(pow(1.+e,2.)*N2_p/D_p + pow(1.-e,2.)*N2_a/D_a))*dLdt;
		*edot = ((1.-e*e)*((1.+e)*N1_p/D_p - (1.-e)*N1_a/D_a)/(2.*p))*dEdt + ((1.-e*e)*((1.+e)*N2_p/D_p - (1.-e)*N2_a/D_a)/(2.*p))*dLdt;
		
		}
		
	return 0;
	
	}

REAL8 XLALDphiDt_geo(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//==================
	// Geodesic dphi/dt
	//==================

	REAL8 Vphi = (Lz-a*E) + (a/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);
	REAL8 Vt = a*(Lz-a*E) + ((r*r+a*a)/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);
	return Vphi/Vt;

	}

REAL8 XLALDphiDt_sf(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi, REAL8 nu) {

	//==========================================
	// dphi/dt with 3PN conservative corrections
	//==========================================

	//Conservative corrections
	REAL8 c0;
	REAL8 c1;
	REAL8 c15;
	REAL8 c2;
	REAL8 c25;
	REAL8 c3;

	c0 	= 0.;
	c1 	= (1./8.)*(1.-pow(e,2.));
	c15 = 5.*a;
	c2 	= -(3./8.)*(15. + 6.*pow(e,2.) - pow(e,4.) + 8.*pow(a,2.) + 8.*sqrt(1.-pow(e,2.))*(-1.+pow(e,2.)));
	c25 = (1./4.)*(267. + 56.*pow(e,2.) + pow(e,4.) + 88.*sqrt(1-pow(e,2.))*(-1.+pow(e,2.)))*a;
	c3 	= (1./384.)*(-57399.+1476.*pow(pi,2.)-52920.*pow(a,2.)
			+ pow(e,2.)*(-29043.+369.*pow(pi,2.)-7320.*pow(a,2.))
			+ pow(e,4.)*(-2613.-48.*pow(a,2.))
			+ pow(e,6.)*(735.)
			+ sqrt(1.-e*e)*(20920.-246.*pow(pi,2.)+11904.*pow(a,2.)
				+ pow(e,2.)*(-16240.+246.*pow(pi,2.)-11904.*pow(a,2.))
				+ pow(e,4.)*(-4680.)));

	REAL8 r = p/(1.+e*cos(psi));
	REAL8 Vphi = (Lz-a*E) + (a/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);
	REAL8 Vt = a*(Lz-a*E) + ((r*r+a*a)/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);

	return (Vphi/Vt)*(1. + nu*(c0 + c1/p + c15/pow(p,1.5) + c2/pow(p,2.) + c25/pow(p,2.5) + c3/pow(p,3.)));
	}

REAL8 XLALDpsiDt_geo(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi) {

	//=========================================
	// dpsi/dt with 3PN conservative corrections
	//=========================================
	
	//Instantaneous radius
	REAL8	r = p/(1.+e*cos(psi));

	//Expression of the form (Geodesic Value)*(1 + nu*(Conservative corrections))
	return  (sqrt(1.-E*E)*sqrt(p*(1.+e))*sqrt(fabs((p - (2./(1.-E*E)-(p/(1.+e)+p/(1.-e)))*(1.+e*cos(psi)))*(1.-e)))
					/((1.-e*e)*(E*(pow((r*r+a*a),2.)/(r*r-2.*r+a*a) - a*a) - 2.*r*a*Lz/(r*r-2.*r+a*a))));
	}

REAL8 XLALDpsiDt_sf(REAL8 E, REAL8 Lz, REAL8 a, REAL8 p, REAL8 e, REAL8 psi, REAL8 nu) {

	//=========================================
	// dpsi/dt with 3PN conservative corrections
	//=========================================
	
	//Instantaneous radius
	REAL8	r = p/(1.+e*cos(psi));

	//Conservative corrections
	REAL8 c0	= 0.;
	REAL8 c1	= (1./8.)*(1.-pow(e,2.));
	REAL8 c15	= 0;
	REAL8 c2	= (3./8.)*pow(1.-pow(e,2.),3./2.)*(8.+sqrt(1.-pow(e,2.)));
	REAL8 c25	= (1./4.)*(1.-pow(e,2.))*(1.-pow(e,2.)-88.*sqrt(1.-pow(e,2.)))*a;
	REAL8 c3	= (1./384.)*pow(1.-e*e,3./2.)*(20920.-246.*pow(pi,2.)+11904.*pow(a,2.)
						+ pow(e,2.)*4680.
						+ sqrt(1.-pow(e,2.))*(-927.-48.*pow(a,2.)
							+ pow(e,2.)*735.));
				
	//Expression of the form (Geodesic Value)*(1 + nu*(Conservative corrections))
	return  (sqrt(1.-E*E)*sqrt(p*(1.+e))*sqrt(fabs((p - (2./(1.-E*E)-(p/(1.+e)+p/(1.-e)))*(1.+e*cos(psi)))*(1.-e)))
					/((1.-e*e)*(E*(pow((r*r+a*a),2.)/(r*r-2.*r+a*a) - a*a) - 2.*r*a*Lz/(r*r-2.*r+a*a))))
							*(1. + nu*(c0 + c1/p + c15/pow(p,1.5) + c2/pow(p,2.)
								+ c25/pow(p,2.5) + c3/pow(p,3.)));

	}


REAL8 XLALPlso(REAL8 a, REAL8 e, INT4 dir) {

	//===================================================
	// LSO semilatus for eccentric Kerr orbits.
	// Follows Komorowsky et. al. (0903.3684) Equations 49-50
	//===================================================

	//Define component terms.
	REAL8 Zii;
	REAL8 ZiA, ZiB;	//Write Zi = ZiA + ZiB*sqrt(Zii). The first term is necessarily real; the second may be imaginary
	REAL8 phi;			//Polar angle made by Zi in the imaginary plane.
	REAL8 AbsZi;		//Abs(Zi)
	REAL8 ReZiThird;	//Real part of Zi^(1/3)
	REAL8 Z0;			//Real part of Zo

	ZiA = pow(3.+e,6.) + pow(a,2.)*(1.+e)*( pow(a,2.)*(1.+e)*(pow(a,2.)*(1.+e)*pow(3.-e,3.)+3.*pow(e,4.)+18.*e*e+459.) - 3.*(e*e+15.)*pow(3.+e,3.) );
	ZiB = 24.*sqrt(3.);
	Zii = pow(1.+e,4.)*pow(a,6.)*(1.-a*a)*((1.-e)*pow(e+3.,3.)-a*a*(1.+e)*pow(3.-e,3.));

	//Check for imaginary Zi (negative Zii)
	if (Zii < 0.0) {
		AbsZi = sqrt(pow(ZiA,2.)+pow(ZiB*sqrt(-1.*Zii),2.));
		phi = atan2(ZiB*sqrt(-1.*Zii),ZiA);
		}
	else {
		AbsZi = sqrt(pow(ZiA + ZiB*sqrt(Zii),2.));
		phi = atan2(0.,ZiA+ZiB*sqrt(Zii));
		}

	ReZiThird = pow(AbsZi,1./3.)*cos(phi/3.);

	Z0 = (1./3.)*pow(a,2.)*(1.+e)*(3.-e) + (1./3.)*pow(3.+e,2.)
		+ (  ((1./3.)*(pow(a,4.)*pow(1.+e,2.)*pow(3.-e,2.) - 2.*pow(a,2.)*(3.+e)*(1.+e)*(e*e+15.) + pow(3.+e,4.)))/pow(AbsZi,2./3.) + (1./3.)   )
		* ReZiThird;

	return (3.+e) + sqrt(Z0) - dir*sqrt(16.*a*a*(1.+e)/sqrt(Z0)-Z0+pow(3.+e,2.)+a*a*(1.+e)*(3.-e));

	}

INT4 XLALPlusEquations(const gsl_vector * x, void *params, gsl_vector * f) {

	//See definition rparams above for structure information
	REAL8 dt		= ((struct rparams *) params)->dt;
	REAL8 wi0 		= ((struct rparams *) params)->wi0;
	REAL8 w0 		= ((struct rparams *) params)->w0;
	REAL8 wi1 		= ((struct rparams *) params)->wi1;
	REAL8 w1		= ((struct rparams *) params)->w1;
	REAL8 wi2 		= ((struct rparams *) params)->wi2;
	REAL8 w2 		= ((struct rparams *) params)->w2;
	REAL8 final_mass 	= ((struct rparams *) params)->final_mass;
	REAL8 final_a 		= ((struct rparams *) params)->final_a;
	REAL8 Sdotn 		= ((struct rparams *) params)->Sdotn;
	REAL8 M 		= ((struct rparams *) params)->M;
	REAL8 ab_factor 	= ((struct rparams *) params)->ab_factor;

    	REAL8 hp_1 	= ((struct rparams *) params)->hp_1;
    	REAL8 dhpdi_1 	= ((struct rparams *) params)->dhpdi_1;
    	REAL8 hp_2 	= ((struct rparams *) params)->hp_2;
    	REAL8 dhpdi_2 	= ((struct rparams *) params)->dhpdi_2;
    	REAL8 hp_3 	= ((struct rparams *) params)->hp_3;
    	REAL8 dhpdi_3	= ((struct rparams *) params)->dhpdi_3;
        
    	//Leading amplitude factors
	REAL8 aone	= 1.+ 2.*Sdotn + Sdotn*Sdotn;
	REAL8 atwo	= 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;

	//Unknown ringdown amplitudes
    	const REAL8 a0n = gsl_vector_get (x, 0);
    	const REAL8 a0p = gsl_vector_get (x, 1);
	const REAL8 a1n = gsl_vector_get (x, 2);
    	const REAL8 a1p = gsl_vector_get (x, 3);
	const REAL8 a2n = gsl_vector_get (x, 4);
    	const REAL8 a2p = gsl_vector_get (x, 5);


	//Functions of the form: (Ringdown Quantity) - (Plunge Quantity). We seek the roots of these functions
	
	//h_plus at time step 1
	const REAL8 yy0 = ab_factor*a0n*(aone-(2./9.)*atwo*w0*final_a) - hp_1;

	//Derivative of h_plus at time step 1
    	const REAL8 yy1 = ab_factor*(aone-(2./9.)*atwo*w0*final_a)*(-a0p*w0 - a0n*wi0) - dhpdi_1*final_mass/(M*dt);

	//h_plus at time step 2
    	const REAL8 yy2 = (((aone - (2./9.)*atwo*w0*final_a)*(a0n*cos((M*dt*w0)/final_mass) - 1.*a0p*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1n*cos((M*dt*w1)/final_mass) - 1.*a1p*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass))
			- hp_2/ab_factor;

	//Derivative of h_plus at time step 2
    	const REAL8 yy3 = (((-1.*(aone - (2./9.)*atwo*w0*final_a)*wi0*(a0n*cos((M*dt*w0)/final_mass) - 1.*a0p*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_a)*(-1.*a0p*w0*cos((M*dt*w0)/final_mass) - 1.*a0n*w0*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_a)*wi1*(a1n*cos((M*dt*w1)/final_mass) - 1.*a1p*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(-1.*a1p*w1*cos((M*dt*w1)/final_mass) - 1.*a1n*w1*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)))
			- dhpdi_2*final_mass/(M*dt*ab_factor);

	//h_plus at time step 3
    	const REAL8 yy4 = (((aone - (2./9.)*atwo*w0*final_a)*(a0n*cos((2.*M*dt*w0)/final_mass) - 1.*a0p*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1n*cos((2.*M*dt*w1)/final_mass) - 1.*a1p*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_a)*(a2n*cos((2.*M*dt*w2)/final_mass) - 1.*a2p*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass))
			- hp_3/ab_factor;

	//Derivative of h_plus at time step 3
	const REAL8 yy5 = (((-1.*(aone - (2./9.)*atwo*w0*final_a)*wi0*(a0n*cos((2.*M*dt*w0)/final_mass) - 1.*a0p*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_a)*(-1.*a0p*w0*cos((2.*M*dt*w0)/final_mass) - 1.*a0n*w0*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_a)*wi1*(a1n*cos((2.*M*dt*w1)/final_mass) - 1.*a1p*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(-1.*a1p*w1*cos((2.*M*dt*w1)/final_mass) - 1.*a1n*w1*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w2*final_a)*wi2*(a2n*cos((2.*M*dt*w2)/final_mass) - 1.*a2p*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_a)*(-1.*a2p*w2*cos((2.*M*dt*w2)/final_mass) - 1.*a2n*w2*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)))
			- dhpdi_3*final_mass/(M*dt*ab_factor);

	gsl_vector_set (f, 0, yy0);
	gsl_vector_set (f, 1, yy1);
	gsl_vector_set (f, 2, yy2);
	gsl_vector_set (f, 3, yy3);
	gsl_vector_set (f, 4, yy4);  
	gsl_vector_set (f, 5, yy5);

	return GSL_SUCCESS;

	}

INT4 XLALCrossEquations(const gsl_vector * z, void *params, gsl_vector * g) {

	//See definition rparams above for structure information
	REAL8 dt		= ((struct rparams *) params)->dt;
	REAL8 wi0 		= ((struct rparams *) params)->wi0;
	REAL8 w0 		= ((struct rparams *) params)->w0;
	REAL8 wi1 		= ((struct rparams *) params)->wi1;
	REAL8 w1		= ((struct rparams *) params)->w1;
	REAL8 wi2 		= ((struct rparams *) params)->wi2;
	REAL8 w2 		= ((struct rparams *) params)->w2;
	REAL8 final_mass 	= ((struct rparams *) params)->final_mass;
	REAL8 final_a 		= ((struct rparams *) params)->final_a;
	REAL8 Sdotn 		= ((struct rparams *) params)->Sdotn;
	REAL8 M 		= ((struct rparams *) params)->M;
	REAL8 ab_factor 	= ((struct rparams *) params)->ab_factor;

	REAL8 hx_1 		= ((struct rparams *) params)->hx_1;
	REAL8 dhxdi_1 		= ((struct rparams *) params)->dhxdi_1;
	REAL8 hx_2 		= ((struct rparams *) params)->hx_2;
	REAL8 dhxdi_2 		= ((struct rparams *) params)->dhxdi_2;
	REAL8 hx_3 		= ((struct rparams *) params)->hx_3;
	REAL8 dhxdi_3 		= ((struct rparams *) params)->dhxdi_3;

	//Leading amplitude factors
	REAL8 aone		= 1.+ 2.*Sdotn + Sdotn*Sdotn;
	REAL8 atwo		= 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;

	//Unknown ringdown amplitudes
	const REAL8 a0c 	= gsl_vector_get (z, 0);
	const REAL8 a0cp 	= gsl_vector_get (z, 1);
	const REAL8 a1c 	= gsl_vector_get (z, 2);
	const REAL8 a1cp 	= gsl_vector_get (z, 3);
	const REAL8 a2c 	= gsl_vector_get (z, 4);
	const REAL8 a2cp 	= gsl_vector_get (z, 5);


	//Functions of the form: (Ringdown Quantity) - (Plunge Quantity). We seek the roots of these functions

	//h_cross at time step 1
	const REAL8 yy0 = ab_factor*a0c*(aone-(2./9.)*atwo*w0*final_a) + hx_1;

	//Derivative of h_cross at time step 1
    	const REAL8 yy1 = ab_factor*(aone-(2./9.)*atwo*w0*final_a)*(a0cp*w0 - a0c*wi0) + dhxdi_1*final_mass/(M*dt);

	//h_cross at time step 2
    	const REAL8 yy2 = (((aone - (2./9.)*atwo*w0*final_a)*(a0c*cos((M*dt*w0)/final_mass) + a0cp*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1c*cos((M*dt*w1)/final_mass) + a1cp*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass))
			+ hx_2/ab_factor;   

	//Derivative of h_cross at time step 2
    	const REAL8 yy3 = (((-1.*(aone - (2./9.)*atwo*w0*final_a)*wi0*(a0c*cos((M*dt*w0)/final_mass) + a0cp*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_a)*(a0cp*w0*cos((M*dt*w0)/final_mass) - 1.*a0c*w0*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_a)*wi1*(a1c*cos((M*dt*w1)/final_mass) + a1cp*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1cp*w1*cos((M*dt*w1)/final_mass) - 1.*a1c*w1*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)))
			+ dhxdi_2*final_mass/(M*dt*ab_factor);      

	//h_cross at time step 3
    	const REAL8 yy4 = (((aone - (2./9.)*atwo*w0*final_a)*(a0c*cos((2.*M*dt*w0)/final_mass) + a0cp*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1c*cos((2.*M*dt*w1)/final_mass) + a1cp*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_a)*(a2c*cos((2.*M*dt*w2)/final_mass) + a2cp*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass))
			+ hx_3/ab_factor;    

	//Derivative of h_cross at time step 3
    	const REAL8 yy5 = (((-1.*(aone - (2./9.)*atwo*w0*final_a)*wi0*(a0c*cos((2.*M*dt*w0)/final_mass) + a0cp*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_a)*(a0cp*w0*cos((2.*M*dt*w0)/final_mass) - 1.*a0c*w0*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_a)*wi1*(a1c*cos((2.*M*dt*w1)/final_mass) + a1cp*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_a)*(a1cp*w1*cos((2.*M*dt*w1)/final_mass) - 1.*a1c*w1*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w2*final_a)*wi2*(a2c*cos((2.*M*dt*w2)/final_mass) + a2cp*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_a)*(a2cp*w2*cos((2.*M*dt*w2)/final_mass) - 1.*a2c*w2*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)))
			+ dhxdi_3*final_mass/(M*dt*ab_factor);

	gsl_vector_set (g, 0, yy0);
	gsl_vector_set (g, 1, yy1);
	gsl_vector_set (g, 2, yy2);
	gsl_vector_set (g, 3, yy3);
	gsl_vector_set (g, 4, yy4);
	gsl_vector_set (g, 5, yy5);

    	return GSL_SUCCESS;

     	}


//=========================================
//=========================================
//
//	CENTRAL FUNCTIONS
//
//=========================================
//=========================================

//INT4 XLALIMRIstart(REAL8 *params, REAL8 t_i, REAL8 dt1, REAL8 dt2, INT4 Npts,
//	REAL8 *t, REAL8 *h_plus, REAL8 *h_cross) {

INT4 XLALIMRIstart(REAL8 m, REAL8 M, REAL8 a, REAL8 Dist, REAL8 Sdotn, REAL8 phi0, REAL8 psi0,
	REAL8 p0, REAL8 e0, REAL8Sequence *hplus, REAL8Sequence *hcross, REAL8 dt_i, INT4 dir, INT4 Npts) {

	/******************************************
	Main Function. Generates inspiral of eccentric, equatorial binary.

	INPUTS
	******************************************/

	FILE *logF;
	logF = fopen("logfile3.txt","w");


	//Get extrinsic source parameters and define trigonometric functions
	REAL8 cos2theta = 2.*Sdotn*Sdotn-1.;

	//Get intrinsic source parameters
	REAL8 	mu;			//Reduced mass in seconds
	REAL8 	nu;			//Reduced mass ratio (mu/(m+M))
	REAL8 	final_mass;	//Final post-merger mass in units of seconds, from Huerta & Gair (1009.1985v3) Eq.40
	REAL8	final_a;	//Final post-merger spin, from hLISAstart()

	mu	= ((m*M)/(m+M));
	nu	= mu/(m+M);
	final_mass 	= (M+m)*(1. + nu*(sqrt(8./9.)-1.) - 0.498*nu*nu);
	final_a	= a - 0.129*a*a*nu - 0.384*a*nu*nu - 2.686*a*nu + 2.*sqrt(3.)*nu - 3.454*nu*nu + 2.353*nu*nu*nu;

	//Get final light ring radius
	REAL8 radiusLR = 2.*(1.+cos((2./3.)*acos(-a)));

	/*******************************************
	Iterate through trajectory.
	Stage 1: Adiabatic Inspiral. Indicated by merging=False.
	Stage 2: Transition Phase. Initiated by merging=True.
	Stage 3: Plunge.
	*******************************************/

	REAL8 p, e, psi, E, L, dpdt, dedt, dpsidt;
	REAL8 p_K2, e_K2, psi_K2, E_K2, L_K2, dpdt_K2, dedt_K2, dphidt_K2, dpsidt_K2, r_K2, drdt_K2, d2rdt2_K2;
	REAL8 p_K3, e_K3, psi_K3, E_K3, L_K3, dpdt_K3, dedt_K3, dphidt_K3, dpsidt_K3, r_K3, drdt_K3, d2rdt2_K3;
	REAL8 p_K4, e_K4, psi_K4, E_K4, L_K4, dpdt_K4, dedt_K4, dphidt_K4, dpsidt_K4, r_K4, drdt_K4, d2rdt2_K4;
	REAL8 drdt=0.;
	REAL8 d2rdt2=0.;
	REAL8 dphidt_old=0.;
	REAL8 dp, de, dphi, dpsi;
	REAL8 t_lso = 0., p_lso = 0.;
	REAL8 E_lso, L_lso, Edot_lso, Ldot_lso;
	REAL8 a_lso, b_lso, k_lso, tau0, R0 = 0.;
	REAL8 dtdtau_lso;
	INT4 i_match, i_lightring = 0, i_plungeCounter = 0;
	REAL8 t_match;
	REAL8 E_plunge = 0., L_plunge = 0.;
	REAL8 d2phidt2 = 0.;

	fprintf(logF,"defining arrays\n");

	REAL8* r_Arr = malloc(Npts*sizeof(REAL8));
	REAL8* phi_Arr = malloc(Npts*sizeof(REAL8));
	REAL8* t = malloc(Npts*sizeof(REAL8));

	//REAL8 r_Arr[Npts];
	//REAL8 phi_Arr[Npts];
	//REAL8 t[Npts];
	REAL8 r=0., phi=0., dphidt=0.;
	REAL8 dt = dt_i;

	fprintf(logF,"dt:	%f\n",dt);
	fprintf(logF,"dir:	%d\n",dir);

//	fprintf(logF,"arrays defined\n");

	//Read off initial values
	p = p0;
	phi_Arr[0] = phi0;
	psi = psi0;
	e = e0;
	t[0] = 0.;

	//fprintf(logF,"starting Inspiral\n");
	//fclose(logF);

	enum stage currentStage = INSPIRAL;
	for (INT4 i=0; i<Npts; i++) {

		//==========================================
		// Direct system to the appropriate iterator
		//==========================================

		switch (currentStage) {

			//===========================
			// Stage 1: Adiabatic Inspiral
			//===========================
			case INSPIRAL:

				//Compute radius, energy, angular momentum
				r_Arr[i] = p/(1.+e*cos(psi));
				XLALComputeEnergyAng(p,e,a,dir,&E,&L);
				
				//Get derivatives
				//dphidt[i] = XLALDphiDt_sf(E[i],L[i],a,p[i],e[i],psi[i],nu);
				//dpsidt[i] = XLALDpsiDt_sf(E[i],L[i],a,p[i],e[i],psi[i],nu);
				dphidt = XLALDphiDt_geo(E,L,a,r_Arr[i]);
				dpsidt = XLALDpsiDt_geo(E,L,a,p,e,psi);
				XLALPEdot(E,L,a,p,e,nu,dir,&dpdt,&dedt);

				//Velocities & accelerations needed for GWs
				drdt = (r_Arr[i]-r_Arr[i-1])/(t[i]-t[i-1]);
				d2rdt2 = (drdt-(r_Arr[i-1]-r_Arr[i-2])/(t[i]-t[i-1]))/(t[i]-t[i-1]);
				d2phidt2 = (dphidt-dphidt_old)/(t[i]-t[i-1]);

				//Get second order RK4 terms
				p_K2 	= p+(dt/2.)*dpdt;
				e_K2 	= e+(dt/2.)*dedt;
				psi_K2	= psi+(dt/2.)*dpsidt;
				XLALComputeEnergyAng(p_K2,e_K2,a,dir,&E_K2,&L_K2);
				XLALPEdot(E_K2,L_K2,a,p_K2,e_K2,nu,dir,&dpdt_K2,&dedt_K2);
				//dphidt_K2 = XLALDphiDt_sf(E_K2,L_K2,a,p_K2,e_K2,psi_K2,nu);
				//dpsidt_K2 = XLALDpsiDt_sf(E_K2,L_K2,a,p_K2,e_K2,psi_K2,nu);
				dphidt_K2 = XLALDphiDt_geo(E_K2,L_K2,a,p_K2/(1.+e_K2*cos(psi_K2)));
				dpsidt_K2 = XLALDpsiDt_geo(E_K2,L_K2,a,p_K2,e_K2,psi_K2);

				//Get third order RK4 terms
				p_K3	= p+(dt/2.)*dpdt_K2;
				e_K3	= e+(dt/2.)*dedt_K2;
				psi_K3	= psi+(dt/2.)*dpsidt_K2;
				XLALComputeEnergyAng(p_K3,e_K3,a,dir,&E_K3,&L_K3);
				XLALPEdot(E_K3,L_K3,a,p_K3,e_K3,nu,dir,&dpdt_K3,&dedt_K3);
				//dphidt_K3 = XLALDpsiDt_sf(E_K3,L_K3,a,p_K3,e_K3,psi_K3,nu);
				//dpsidt_K3 = XLALDpsiDt_sf(E_K3,L_K3,a,p_K3,e_K3,psi_K3,nu);
				dphidt_K3 = XLALDphiDt_geo(E_K3,L_K3,a,p_K3/(1.+e_K3*cos(psi_K3)));
				dpsidt_K3 = XLALDpsiDt_geo(E_K3,L_K3,a,p_K3,e_K3,psi_K3);

				//Get fourth order RK4 terms
				p_K4	= p+dt*dpdt_K3;
				e_K4	= e+dt*dedt_K3;
				psi_K4	= psi+(dt)*dpsidt_K3;
				XLALComputeEnergyAng(p_K4,e_K4,a,dir,&E_K4,&L_K4);
				XLALPEdot(E_K4,L_K4,a,p_K4,e_K4,nu,dir,&dpdt_K4,&dedt_K4);
				//dphidt_K4 = XLALDpsiDt_sf(E_K4,L_K4,a,p_K4,e_K4,psi_K4,nu);
				//dpsidt_K4 = XLALDpsiDt_sf(E_K4,L_K4,a,p_K4,e_K4,psi_K4,nu);
				dphidt_K4 = XLALDphiDt_geo(E_K4,L_K4,a,p_K4/(1.+e_K4*cos(psi_K4)));
				dpsidt_K4 = XLALDpsiDt_geo(E_K4,L_K4,a,p_K4,e_K4,psi_K4);
				
				//GW amplitudes. Following Huerta & Gair (1009.1985v3) Eqs. 37 and 38.
				hplus->data[i] = (mu/(2.*Dist))*((1. - 2.*cos2theta*pow(cos(phi_Arr[i]),2.) - 3.*cos(2.*phi_Arr[i]))*drdt*drdt
						+ (3. + cos2theta)*(2.*cos(2.*phi_Arr[i])*dphidt*dphidt + sin(2.*phi_Arr[i])*d2phidt2)*r_Arr[i]*r_Arr[i]
						+ (4.*(3.+cos2theta)*sin(2.*phi_Arr[i])*dphidt*drdt
							+ (1.-2.*cos2theta*pow(cos(phi_Arr[i]),2.)-3.*cos(2.*phi_Arr[i]))*d2rdt2)*r_Arr[i]);
				hcross->data[i] = (-2.*mu/Dist)*Sdotn*(sin(2.*(phi_Arr[i]))*drdt*drdt
						+ r_Arr[i]*r_Arr[i]*(-2.*sin(2.*(phi_Arr[i]))*dphidt*dphidt + cos(2.*( phi_Arr[i]))*d2phidt2)
						+ r_Arr[i]*(4.*cos(2.*(phi_Arr[i]))*dphidt*drdt + sin(2.*(phi_Arr[i]))*d2rdt2));

				//Trial step. Update LSO location
				dp = (1./6.)*dt*(dpdt+2.*dpdt_K2+2.*dpdt_K3+dpdt_K4);
				de = (1./6.)*dt*(dedt+2.*dedt_K2+2.*dedt_K3+dedt_K4);
				dphi = (1./6.)*dt*(dphidt+2.*dphidt_K2+2.*dphidt_K3+dphidt_K4);
				dpsi = (1./6.)*dt*(dpsidt+2.*dpsidt_K2+2.*dpsidt_K3+dpsidt_K4);
				p_lso = XLALPlso(a,e+de,dir);

				//If we've skipped too far past the LSO, back up and take a smaller step
				if ((p_lso > p+dp) || (p+dp != p+dp)) {

					i = i-1;
					dt = dt/10.;

					break;

					}

				//Update Coordinates
				t[i+1] = t[i] + dt;
				phi_Arr[i+1] = phi_Arr[i] + dphi; 
				p += dp;
				e += de;
				psi += dpsi;
				dphidt_old = dphidt;

				//If this condition is met, then the current p[i] must very well approximate p_lso. Initiate transition.
				if (fabs(p-p_lso) <= 0.01) {

					fprintf(logF,"transitioning\n");
					fclose(logF);

					//Record LSO coordinates
					t_lso	= t[i+1];

					//LSO energy, ang, and fluxes
					XLALComputeEnergyAng(p_lso,e,a,dir,&E_lso,&L_lso);
					XLALEnergyAngDot(a,p_lso,e,nu,dir,&Edot_lso,&Ldot_lso);

					//Various helper constants from Sundararajan
					dtdtau_lso = pow(p_lso,-2.0)*( a*(L_lso-a*E_lso)
						+ ((p_lso*p_lso+a*a)/(p_lso*p_lso-2.*p_lso+a*a))*(E_lso*(p_lso*p_lso+a*a)-L_lso*a) ); //(Vt/Sigma)
					a_lso = (3./pow(p_lso,6.)) * (p_lso*p_lso + 2.*(a*a*(E_lso*E_lso-1.)-L_lso*L_lso)*p_lso
						+ 10.*pow(L_lso-a*E_lso,2.));
					b_lso = (2./pow(p_lso,4.)) * ( (L_lso - a*a*E_lso*(Edot_lso/Ldot_lso))*p_lso
						- 3.*(L_lso - a*E_lso)*(1. - a*(Edot_lso/Ldot_lso)) );
					k_lso = -(1./nu)*Ldot_lso*dtdtau_lso;
					tau0 = pow(a_lso*b_lso*k_lso,-1./5.);
					R0 = pow(b_lso*k_lso,2./5.)*pow(a_lso,-3./5.);

					//Find the time at T = -1 and the corresponding index i_match. Transition will begin at this index.
					i_match = XLALSearchArr(t,t_lso + (-1.)*(dtdtau_lso*tau0*pow(nu,-1./5.)));
					t_match = t[i_match];
					r = r_Arr[i_match+1];
					phi = phi_Arr[i_match+1];

					//Find initial velocity, acceleration
					drdt = (r_Arr[i_match]-r_Arr[i_match-1])/(t[i_match]-t[i_match-1]);
					dphidt_old = XLALDphiDt_geo(E_lso+(t_match-t_lso)*Edot_lso, L_lso+(t_match-t_lso)*Ldot_lso,a,r);

					//Jump back to one timestep before transition match. Switch to finer time increments
					i = i_match;
					dt = dt_i;
					currentStage = TRANSITION;

					free(r_Arr);
					free(phi_Arr);
					
					}

				break;

			/******************************
			Stage 2: Transition
			******************************/
			case TRANSITION:

				//Compute transition energy, angular momentum			
				L = L_lso + (t[i]-t_lso)*Ldot_lso;
				E = E_lso + (t[i]-t_lso)*Edot_lso;
				
				//Define dimensionless COORDINATE acceleration and velocities
				d2rdt2 = (1./2.)*(XLALD2FDrDE(E_lso,L_lso,a,r)*(E - E_lso) + XLALD2FDrDL(E_lso,L_lso,a,r)*(L - L_lso) + XLALDFDr(E_lso,L_lso,a,r));
				dphidt = XLALDphiDt_geo(E,L,a,r);
				d2phidt2 = (dphidt-dphidt_old)/(t[i]-t[i-1]);

				//Get second and third order RK4 terms. These are equal, as they are each evaluated at (t[i] + 1/2 dt) and d2rdt2 depends only on time
				L_K2 = L_lso + ((t[i]+dt/2.)-t_lso)*Ldot_lso;
				E_K2 = E_lso + ((t[i]+dt/2.)-t_lso)*Edot_lso;
				r_K2 = r + (dt/2.)*drdt;
				drdt_K2 = drdt + (dt/2.)*d2rdt2;
				dphidt_K2 = XLALDphiDt_geo(E_K2,L_K2,a,r_K2);
				d2rdt2_K2 = (1./2.)*(XLALD2FDrDE(E_lso,L_lso,a,r_K2)*(E_K2 - E_lso) + XLALD2FDrDL(E_lso,L_lso,a,r_K2)*(L_K2 - L_lso) + XLALDFDr(E_lso,L_lso,a,r_K2));

				//Get third order RK4 terms. For E,L,d2rdt2 these equal the second order terms as they are each evaluated at (t[i] + 1/2 dt)
				L_K3 = L_K2;
				E_K3 = E_K2;
				r_K3 = r + (dt/2.)*drdt_K2;
				drdt_K3 = drdt + (dt/2.)*d2rdt2_K2;
				dphidt_K3 = XLALDphiDt_geo(E_K3,L_K3,a,r_K3);
				d2rdt2_K3 = (1./2.)*(XLALD2FDrDE(E_lso,L_lso,a,r_K3)*(E_K3 - E_lso) + XLALD2FDrDL(E_lso,L_lso,a,r_K3)*(L_K3 - L_lso) + XLALDFDr(E_lso,L_lso,a,r_K3));

				//Get fourth order RK4 terms (evaluated at t = t[i] + dt)
				L_K4 = L_lso + ((t[i]+dt)-t_lso)*Ldot_lso;
				E_K4 = E_lso + ((t[i]+dt)-t_lso)*Edot_lso;
				r_K4 = r + dt*drdt_K3;
				drdt_K4 = drdt + dt*d2rdt2_K3;
				dphidt_K4 = XLALDphiDt_geo(E_K4,L_K4,a,r_K4);
				d2rdt2_K4 = (1./2.)*(XLALD2FDrDE(E_lso,L_lso,a,r_K4)*(E_K4 - E_lso) + XLALD2FDrDL(E_lso,L_lso,a,r_K4)*(L_K4 - L_lso) + XLALDFDr(E_lso,L_lso,a,r_K4));

				//GW amplitudes
				hplus->data[i] = (mu/(2.*Dist))*((1. - 2.*cos2theta*pow(cos(phi),2.) - 3.*cos(2.*phi))*drdt*drdt
								+ (3. + cos2theta)*(2.*cos(2.*phi)*dphidt*dphidt + sin(2.*phi)*d2phidt2)*r*r
								+ (4.*(3.+cos2theta)*sin(2.*phi)*dphidt*drdt
									+ (1.-2.*cos2theta*pow(cos(phi),2.)-3.*cos(2.*phi))*d2rdt2)*r);
				hcross->data[i] = (-2.*mu/Dist)*Sdotn*(sin(2.*(phi))*drdt*drdt
								+ r*r*(-2.*sin(2.*(phi))*dphidt*dphidt + cos(2.*( phi))*d2phidt2)
								+ r*(4.*cos(2.*(phi))*dphidt*drdt + sin(2.*(phi))*d2rdt2));

				//Update coordinates
				dphidt_old = dphidt;
				r += (1./6.)*dt*(drdt + 2.*drdt_K2 + 2.*drdt_K3 + drdt_K4);
				phi += (1./6.)*dt*(dphidt + 2.*dphidt_K2 + 2.*dphidt_K3 + dphidt_K4);
				drdt += (1./6.)*dt*(d2rdt2 + 2.*d2rdt2_K2 + 2.*d2rdt2_K3 + d2rdt2_K4);
				t[i+1] = t[i]+dt;

				//CHECK THE NU^(-2/5) (COPIED FROM MY FINAL RUNCODE ECCENTRIC IMRI VERSION)
				if (pow(nu,-2./5.)*(r-p_lso)/R0 < -1.0) {

					currentStage = PLUNGE;
					E_plunge = E_lso + (t[i+1]-t_lso)*Edot_lso;
					L_plunge = L_lso + (t[i+1]-t_lso)*Ldot_lso;

					}
				
				break;

			/******************************
			Stage 3: Plunge
			******************************/
			case PLUNGE:

				//Get RK4 time derivatives
				dphidt = XLALDphiDt_geo(E_plunge,L_plunge,a,r);
				d2rdt2 = (1./2.)*XLALDFDr(E_plunge,L_plunge,a,r);
				d2phidt2 = XLALDVphiVtDr(E_plunge,L_plunge,a,r)*drdt;

				r_K2 = r + (dt/2.)*drdt;
				drdt_K2 = drdt + (dt/2.)*d2rdt2;
				dphidt_K2 = XLALDphiDt_geo(E_plunge,L_plunge,a,r_K2);
				d2rdt2_K2 = (1./2.)*XLALDFDr(E_plunge,L_plunge,a,r_K2);

				r_K3 = r + (dt/2.)*drdt_K2;
				drdt_K3 = drdt + (dt/2.)*d2rdt2_K2;
				dphidt_K3 = XLALDphiDt_geo(E_plunge,L_plunge,a,r_K3);
				d2rdt2_K3 = (1./2.)*XLALDFDr(E_plunge,L_plunge,a,r_K3);

				r_K4 = r + (dt)*drdt_K3;
				drdt_K4 = drdt + (dt)*d2rdt2_K2;
				//drdt_K4 = drdt + (dt)*d2rdt2_K3;
				dphidt_K4 = XLALDphiDt_geo(E_plunge,L_plunge,a,r_K4);
				d2rdt2_K4 = (1./2.)*XLALDFDr(E_plunge,L_plunge,a,r_K4);

				//GW amplitudes
				hplus->data[i] = (mu/(2.*Dist))*((1. - 2.*cos2theta*pow(cos(phi),2.) - 3.*cos(2.*phi))*drdt*drdt
								+ (3. + cos2theta)*(2.*cos(2.*phi)*dphidt*dphidt + sin(2.*phi)*d2phidt2)*r*r
								+ (4.*(3.+cos2theta)*sin(2.*phi)*dphidt*drdt
									+ (1.-2.*cos2theta*pow(cos(phi),2.)-3.*cos(2.*phi))*d2rdt2)*r);
				hcross->data[i] = (-2.*mu/Dist)*Sdotn*(sin(2.*(phi))*drdt*drdt
								+ r*r*(-2.*sin(2.*(phi))*dphidt*dphidt + cos(2.*( phi))*d2phidt2)
								+ r*(4.*cos(2.*(phi))*dphidt*drdt + sin(2.*(phi))*d2rdt2));

				//Update coordinates
				r += (1./6.)*dt*(drdt + 2.*drdt_K2 + 2.*drdt_K3 + drdt_K4);
				phi += (1./6.)*dt*(dphidt + 2.*dphidt_K2 + 2.*dphidt_K3 + dphidt_K4);
				drdt += (1./6.)*dt*(d2rdt2 + 2.*d2rdt2_K2 + 2.*d2rdt2_K3 + d2rdt2_K4);
				t[i+1] = t[i]+dt;

				if (r < radiusLR) {

					//Record light ring
					if (i_lightring == 0) i_lightring = i+1;

					//After 10 iterations inside light ring, stop the loop
					i_plungeCounter += 1;

					if (i_plungeCounter > 10) i = Npts;

					}

				break;

			}

		}

	if (i_lightring==0) {
		return Npts;
		}

	/**********************************
	Get dominant ringdown frequencies
	**********************************/

	INT4 N = 11;

	//Final q values. Formatted as REAL8 for compatibility with gsl_spline_...
	REAL8 fq[11] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=0). From Berti et. al. Frequencies nondimensionlialized using final BH mass.
	REAL8 w0_arr[11] 	= {0.3737, 0.3870, 0.4021, 0.4195, 0.4398, 0.4641, 0.4940, 0.5326, 0.5860, 0.6716, 0.8254};
	REAL8 wi0_arr[11] = {0.0890, 0.0887, 0.0883, 0.0877, 0.0869, 0.0856, 0.0838, 0.0808, 0.0756, 0.0649, 0.0386};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=1). From Berti et. al.
	REAL8 w1_arr[11] 	= {0.3467, 0.3619, 0.3790, 0.3984, 0.4208, 0.4474, 0.4798, 0.5212, 0.5779, 0.6677, 0.8249};
	REAL8 wi1_arr[11]	= {0.2739, 0.2725, 0.2705, 0.2680, 0.2647, 0.2602, 0.2538, 0.2442, 0.2281, 0.1953, 0.1159};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=2). From Berti et. al.
	REAL8 w2_arr[11] 	= {0.3011, 0.3192, 0.3393, 0.3619, 0.3878, 0.4179, 0.4542, 0.4999, 0.5622, 0.6598, 0.8238};
	REAL8 wi2_arr[11]	= {0.4783, 0.4735, 0.4679, 0.4613, 0.4533, 0.4433, 0.4303, 0.4123, 0.3839, 0.3275, 0.1933};

	//Set up interpolator
	const gsl_interp_type *bn = gsl_interp_cspline;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();

	//Define a bunch of interpolation splines
	gsl_spline *spline0 = gsl_spline_alloc(bn,N);
	gsl_spline *spline1 = gsl_spline_alloc(bn,N);
	gsl_spline *spline2 = gsl_spline_alloc(bn,N);
	gsl_spline *spline3 = gsl_spline_alloc(bn,N);
	gsl_spline *spline4 = gsl_spline_alloc(bn,N);
	gsl_spline *spline5 = gsl_spline_alloc(bn,N);

	//Initialize splines
	gsl_spline_init(spline0, fq, w0_arr, N);
	gsl_spline_init(spline1, fq, wi0_arr, N);
	gsl_spline_init(spline2, fq, w1_arr, N);
	gsl_spline_init(spline3, fq, wi1_arr, N);
	gsl_spline_init(spline4, fq, w2_arr, N);
	gsl_spline_init(spline5, fq, wi2_arr, N);

	//Get the real and damping frequencies for our final BH
	REAL8 w0, wi0;
	REAL8 w1, wi1;
	REAL8 w2, wi2;

	w0 	= gsl_spline_eval(spline0, final_a, acc);
	wi0 	= gsl_spline_eval(spline1, final_a, acc);
	w1	= gsl_spline_eval(spline2, final_a, acc);
	wi1	= gsl_spline_eval(spline3, final_a, acc);
	w2	= gsl_spline_eval(spline4, final_a, acc);
	wi2	= gsl_spline_eval(spline5, final_a, acc);

	//Free memory
	gsl_spline_free(spline0);
	gsl_spline_free(spline1);
	gsl_spline_free(spline2);
	gsl_spline_free(spline3);
	gsl_spline_free(spline4);
	gsl_spline_free(spline5);

	/**************************************
	Find interpolation function for h across final plunge
	**************************************/

	//Define arrays to hold the final 20 values of hp, hx, and their derivatives during plunge. Declared as REAL8s for compatibility with gsl
	REAL8 i_interp[20];
	REAL8 hp_interp[20];
	REAL8 hx_interp[20];

	//Fit values for ringdown amplitude equation
	REAL8 a0n, a0p;
	REAL8 a1n, a1p;
	REAL8 a2n, a2p;
	REAL8 a0c, a0cp;
	REAL8 a1c, a1cp;
	REAL8 a2c, a2cp;

	//Read out the final twenty values of hplus and hcross.
	//Light ring bisects r[i_lightring] and r[i_lightring-1].
	//Note: i + i_lightring - 10 ranges from (i_lightring - 10) to (i_lightring + 9)
	for (INT4 i=0; i<20; i++) {

		i_interp[i]	= i+(i_lightring-10);
		hp_interp[i] 	= hplus->data[i+(i_lightring-10)]*Dist;
		hx_interp[i] 	= hcross->data[i+(i_lightring-10)]*Dist;

		}

	//Define hp, hx, and their first derivatives with respect to i at three steps across the light radius.
	//All have units of seconds

	REAL8 hp_1, dhpdi_1;
	REAL8 hx_1, dhxdi_1;
	REAL8 hp_2, dhpdi_2;
	REAL8 hx_2, dhxdi_2;
	REAL8 hp_3, dhpdi_3;
	REAL8 hx_3, dhxdi_3;

	gsl_spline *spline_hp = gsl_spline_alloc(gsl_interp_cspline, 20);	//To interpolate across hp_interp[]
	gsl_spline *spline_hx = gsl_spline_alloc(gsl_interp_cspline, 20);	//To interpolate across hx_interp[]
	gsl_spline_init(spline_hp, i_interp, hp_interp, 20);
	gsl_spline_init(spline_hx, i_interp, hx_interp, 20);

	//Evaluate first point at i_interp[9] = i_lightradius-1
	hp_1		= gsl_spline_eval(spline_hp, i_interp[9], acc);
	dhpdi_1		= gsl_spline_eval_deriv(spline_hp, i_interp[9], acc);
	hx_1		= gsl_spline_eval(spline_hx, i_interp[9], acc);
	dhxdi_1		= gsl_spline_eval_deriv(spline_hx, i_interp[9], acc);

	//Evaluate second point at i_interp[10] = i_lightradius
	hp_2		= gsl_spline_eval(spline_hp, i_interp[10], acc);
	dhpdi_2		= gsl_spline_eval_deriv(spline_hp, i_interp[10], acc);
	hx_2		= gsl_spline_eval(spline_hx, i_interp[10], acc);
	dhxdi_2		= gsl_spline_eval_deriv(spline_hx, i_interp[10], acc);

	//Evaluate third point at i_interp[11] = i_lightradius+1
	hp_3		= gsl_spline_eval(spline_hp, i_interp[11], acc);
	dhpdi_3		= gsl_spline_eval_deriv(spline_hp, i_interp[11], acc);
	hx_3		= gsl_spline_eval(spline_hx, i_interp[11], acc);
	dhxdi_3		= gsl_spline_eval_deriv(spline_hx, i_interp[11], acc);

	//Free memory
	gsl_spline_free(spline_hp);
	gsl_spline_free(spline_hx);

	//Create an instance of rparams, populate with literally every variable
	struct rparams par;
	par.dt		= dt;
	par.a 		= a;
	par.M		= M+m;
	par.Sdotn	= Sdotn;
	par.final_mass	= final_mass;
	par.final_a	= final_a;
	par.w0		= w0;
	par.w1		= w1;
	par.w2		= w2;
	par.wi0		= wi0;
	par.wi1		= wi1;
	par.wi2		= wi2;
	par.ab_factor	= (1./8.)*sqrt(5./pi)*(final_mass);
	par.hp_1	= hp_1;
	par.hp_2	= hp_2;
	par.hp_3	= hp_3;
	par.dhpdi_1	= dhpdi_1;
	par.dhpdi_2	= dhpdi_2;
	par.dhpdi_3	= dhpdi_3;
	par.hx_1	= hx_1;
	par.hx_2	= hx_2;
	par.hx_3	= hx_3;
	par.dhxdi_1	= dhxdi_1;
	par.dhxdi_2	= dhxdi_2;
	par.dhxdi_3	= dhxdi_3;

	//Initialize root solvers
	const gsl_multiroot_fsolver_type *T;
	const gsl_multiroot_fsolver_type *Q;
	gsl_multiroot_fsolver *plus_solver;
	gsl_multiroot_fsolver *cross_solver;
	INT4 statusp, statusc;
	size_t iter = 0;
	const size_t nmatch = 6;

	//Define functions
	gsl_multiroot_function f_plus = {&XLALPlusEquations, nmatch, &par};
	gsl_multiroot_function f_cross = {&XLALCrossEquations, nmatch, &par};

	//Define two 6D vectors, and feed initial guesses
	REAL8 x_init[6] = {-0.392254,  4.575194, -0.870431,  5.673678,  0.979356, -3.637467};

	REAL8 z_init[6] = {-0.392254,  4.575194, -0.870431,  5.673678,  0.979356, -3.637467};
	gsl_vector *xmr = gsl_vector_alloc(nmatch);
	gsl_vector *z 	= gsl_vector_alloc(nmatch);

	for (size_t i=0; i<nmatch; i++) {
		gsl_vector_set(xmr, i, x_init[i]);
		gsl_vector_set(z, i, z_init[i]);
		}

	//Set up solvers
	T = gsl_multiroot_fsolver_hybrids;
	Q = gsl_multiroot_fsolver_hybrids;
	plus_solver = gsl_multiroot_fsolver_alloc(T, 6);
	cross_solver = gsl_multiroot_fsolver_alloc(Q, 6);
	gsl_multiroot_fsolver_set(plus_solver, &f_plus, xmr);
	gsl_multiroot_fsolver_set(cross_solver, &f_cross, z);

	REAL8 fminval = 1.e-5;
	size_t MaxITS = 200000;

	//Find root of PlusEquations
	do {
		iter++;
		statusp = gsl_multiroot_fsolver_iterate(plus_solver);
		if (statusp) break;
		statusp = gsl_multiroot_test_residual(plus_solver->f, fminval);
		}
	while (statusp == GSL_CONTINUE && iter < MaxITS);

	//Find root of CrossEquations
	iter = 0;
	do {
		iter++;
		statusc = gsl_multiroot_fsolver_iterate(cross_solver);
		if (statusc) break;
		statusc = gsl_multiroot_test_residual(cross_solver->f, fminval);
		}
	while (statusc == GSL_CONTINUE && iter < MaxITS);

	//Read out fit parameters
	a0n	= gsl_vector_get(plus_solver->x,0);
	a0p	= gsl_vector_get(plus_solver->x,1);
	a1n	= gsl_vector_get(plus_solver->x,2);
	a1p	= gsl_vector_get(plus_solver->x,3);
	a2n	= gsl_vector_get(plus_solver->x,4);
	a2p	= gsl_vector_get(plus_solver->x,5);
	a0c	= gsl_vector_get(cross_solver->x,0);
	a0cp	= gsl_vector_get(cross_solver->x,1);
	a1c		= gsl_vector_get(cross_solver->x,2);
	a1cp	= gsl_vector_get(cross_solver->x,3);
	a2c		= gsl_vector_get(cross_solver->x,4);
	a2cp	= gsl_vector_get(cross_solver->x,5);

	//Free memory
	gsl_multiroot_fsolver_free(plus_solver);
	gsl_multiroot_fsolver_free(cross_solver);
	gsl_vector_free(xmr);
	gsl_vector_free(z);

	/*************************************************************************
	Recompute wave amplitudes and detector responses from (i_lightring-1) to i_max
	*************************************************************************/

	REAL8 aone, atwo;
	REAL8 hp0, hp1, hp2;
	REAL8 hc0, hc1, hc2;
	REAL8 ringdown_amp;
	REAL8 hp_ringdown, hc_ringdown;

	aone = 1.+ 2.*Sdotn + Sdotn*Sdotn;
	atwo = 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;
	ringdown_amp = (1./8.)*sqrt(5./pi)*(final_mass/Dist);

	for (INT4 i=i_lightring-1; i<Npts; i++) {

		//Update time
		t[i] 	 = t[i-1] + dt;

		//Compute coefficients. Factor of (M/final_mass) rescales the dimensionless frequencies w to match the scaling of dt (or vice-versa).
		//Note: final_q and the frequencies wi are both nondimensionalized by final_mass, so their direct product is scaled properly
		//In HG, the unknown constants are complex amplitude and phase. These unknowns correspond to two unknown amplitudes a0n and a0p of cos() and sin() contributions to Fourier superposition
		hp0 = exp(-(i-(i_lightring-1))*dt*(wi0*(m+M)/final_mass))*( aone-(2./9.)*final_a*w0*atwo )*( a0n*cos(w0*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) - a0p*sin(w0*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));
                hp1 = exp(-(i-(i_lightring-1))*dt*(wi1*(m+M)/final_mass))*( aone-(2./9.)*final_a*w1*atwo )*( a1n*cos(w1*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) - a1p*sin(w1*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));
                hp2 = exp(-(i-(i_lightring-1))*dt*(wi2*(m+M)/final_mass))*( aone-(2./9.)*final_a*w2*atwo )*( a2n*cos(w2*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) - a2p*sin(w2*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));
                hc0 = exp(-(i-(i_lightring-1))*dt*(wi0*(m+M)/final_mass))*( aone-(2./9.)*final_a*w0*atwo )*( a0c*cos(w0*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) + a0cp*sin(w0*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));
                hc1 = exp(-(i-(i_lightring-1))*dt*(wi1*(m+M)/final_mass))*( aone-(2./9.)*final_a*w1*atwo )*( a1c*cos(w1*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) + a1cp*sin(w1*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));
                hc2 = exp(-(i-(i_lightring-1))*dt*(wi2*(m+M)/final_mass))*( aone-(2./9.)*final_a*w2*atwo )*( a2c*cos(w2*(i-(i_lightring-1))*(dt*(m+M)/final_mass)) + a2cp*sin(w2*(i-(i_lightring-1))*(dt*(m+M)/final_mass)));	

		//Get hplus and hcross
		hp_ringdown	= (hp0+hp1+hp2);
		hc_ringdown	= -(hc0+hc1+hc2);

		//Record wave amplitudes
		hplus->data[i] 	= ringdown_amp*hp_ringdown;
		hcross->data[i] 	= ringdown_amp*hc_ringdown;

		if (fabs(hplus->data[i]) < 1.e-30) {
			free(t);
			hplus->length=i;
			hcross->length=i;
			return i;
			}

		}

	return Npts;

	}

//=========================================
//=========================================
//
//	MAIN
//
//=========================================
//=========================================

INT4 XLALIMRIGenerator(
	REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        REAL8 phi0,                     //Initial phi
        REAL8 dt,                       //Time step (in seconds)
        REAL8 m1,                       //BH mass (passed in kg)
        REAL8 m2,                       //CO mass (passed in kg)
        REAL8 f_min,                    //Initial frequency (passed in Hz)
        REAL8 r,                        //Distance to system (passed in meters)
        REAL8 inc,                      //Inclination angle between line of sight \hat n and BH spin axis
        REAL8 s1z,                       //BH reduced spin
	LALSimInspiralTestGRParam *testGR
        ) {

	FILE *logF;
	logF=fopen("logfile2.txt","w");
	fprintf(logF,"Starting generator\n");

        REAL8 a = s1z;
        REAL8 Sdotn = cos(inc);
        REAL8 M = m1/LAL_MSUN_SI; 
	REAL8 m = m2/LAL_MSUN_SI;
	REAL8 dist = r/(LAL_PC_SI*1.e9);

	fprintf(logF,"defining shit 1\n");

	REAL8 e0=0., psi0=0.;
	if (XLALSimInspiralTestGRParamExists(testGR,"ecc") == true)
		e0 = XLALSimInspiralGetTestGRParam(testGR,"ecc");
	else printf("No eccentricity...");

	if (XLALSimInspiralTestGRParamExists(testGR,"psi0") == true)
		psi0 = XLALSimInspiralGetTestGRParam(testGR,"psi0");
	else printf("No eccentricity...");
	

	//REAL8 e0 = 0.1;
	//REAL8 psi0=0.0;


	INT4 dir;
	if (a>0.0) dir = 1;
	else dir = -1;

	INT4 status;
	size_t iter=0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	struct pParams p0Params;
	p0Params.m = m;
	p0Params.M = M;
	p0Params.a = a;
	p0Params.e0 = e0;
	p0Params.dir = dir;
	p0Params.f0 = f_min;

	gsl_function F;
	F.function = &XLALInitialP;
	F.params = &p0Params;

	REAL8 p0 = 0.;
	REAL8 p_high = 1000.;
	REAL8 p_low = XLALPlso(a,e0,dir);

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s,&F,p_low,p_high);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		p0 = gsl_root_fsolver_root(s);
                p_low = gsl_root_fsolver_x_lower(s);
                p_high = gsl_root_fsolver_x_upper(s);
                status = gsl_root_test_interval (p_low, p_high, 0, 0.001);
                }
        while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fsolver_free(s);

        if (status != GSL_SUCCESS) {
                XLALPrintError("XLAL Error: Rootfinding failed. Could not find inital \
			radius 6M < r < 500M satisfying specified initial frequency.\n");
                XLAL_ERROR(XLAL_EFAILED); //"generic failure"

                }

	fprintf(logF,"found initial p\n");
	fprintf(logF,"%f",p0);

	REAL8 fDyne = 0.0;
	size_t length = 500/dt;
	static LIGOTimeGPS epoch;
	*hplus = XLALCreateREAL8TimeSeries("h_plus",&epoch,fDyne,dt,&lalStrainUnit,length);
	*hcross = XLALCreateREAL8TimeSeries("h_cross",&epoch,fDyne,dt,&lalStrainUnit,length);

	if (*hplus == NULL || *hcross == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	fprintf(logF,"\nStarting sim\n");
	fprintf(logF,"%f %f %f %f %f",m,M,a,e0,dt);

	fclose(logF);

	//IMRIstart....
	XLALIMRIstart(m*Msun_sec, M*Msun_sec, a, dist*GPC_sec, Sdotn, phi0, psi0, p0, e0, (*hplus)->data, (*hcross)->data,
		dt/((m+M)*Msun_sec), dir, length);
	

	return 0;

	}

/*
INT4 main(void) {

	printf("Running main\n");

	//1 for Yes, 0 for No
	//INT4 recordAscii = 0;

	INT4 direction = 1;
	REAL8 m = 10.;
	REAL8 M = 100.;
	REAL8 a = 0.3;
	REAL8 Dist = 0.1;
	REAL8 THETA_S = 0.0;
	REAL8 PHI_S = 0.0;
	REAL8 THETA_K = 0.0;
	REAL8 PHI_K	= 0.0;
	REAL8 THETA0 = 0.0;
	REAL8 PHI0 = 0.0;
	REAL8 phi0 = 0.0;
	REAL8 psi0 = 0.0;
	REAL8 e0= 0.3;

	REAL8 f0 = 20.;
	REAL8 p0 = pow(1./(pi*(f0*(m+M)*Msun_sec))-a,2./3.);
	
	REAL8 Sdotn = cos(THETA_S)*cos(THETA_K) + sin(THETA_S)*sin(THETA_K)*cos(PHI_S-PHI_K);
	printf("%f\n",Sdotn);
	//Step size
	REAL8 dt1 = 0.0001;
	//REAL8 dt2 = 0.0001;

	//Initialize data arrays
	printf("initializing...\n");
	size_t maxPts = 1000/dt1;
	REAL8 t[maxPts];
	REAL8 h_plus[maxPts];
	REAL8 h_cross[maxPts];

	t[0]=2.;
	h_plus[0]=2.;
	h_cross[0]=2.;

	//Run
	//REAL8 params[15] = {m,M,a,direction,Dist,Sdotn,THETA0,PHI0,phi0,psi0,p0,e0};
	printf("Starting...\n");
	printf("%d	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f",
		direction,m,M,a,Dist,THETA_S,PHI_S,THETA_K,PHI_K,THETA0,PHI0,phi0,psi0,e0,f0,p0);
	//INT4 finalInt = XLALIMRIstart(params,0.0,dt1/((m+M)*Msun_sec),dt2/((m+M)*Msun_sec),maxPts-1,t,h_plus,h_cross);

	//Write to ascii file
	if (recordAscii) {

		//Initialize files to hold output
		FILE *ascii_file;
		ascii_file = fopen("./asciiData/dat","w");

		//Record
		for (INT4 i=2; i<finalInt; i+=1) {
			fprintf(ascii_file,"%i	%f	%f	%f",i,(m+M)*Msun_sec*t[i],h_plus[i],h_cross[i]);
			}

		fclose(ascii_file);

		}

	}*/

