/*
*  Copyright (C) 2014 Chinmay Kalaghatgi, P. Ajith
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/


#include <lal/AVFactories.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <lal/LALStdlib.h>

typedef struct tagConstantDefinitions{
	REAL8 ChiPow217;	/* Powers of chi are defined so that they are computed only once during the computation of metric */
	REAL8 ChiPow26;
	REAL8 ChiPow783;
	REAL8 ChiPow152;
	REAL8 ChiPow3;
	REAL8 ChiPow45;
	REAL8 ChiPow74;
	REAL8 ChiPow75;
	REAL8 ChiPow7;
}
ConstantDefinitions;

typedef struct tagNormalizationParameters{
	REAL8 norm_merg;
	REAL8 norm_merg_mass;
	REAL8 norm_merg_eta;
	REAL8 norm_merg_chi;
	REAL8 norm_ring;
	REAL8 norm_ring_mass;
	REAL8 norm_ring_eta;
	REAL8 norm_ring_chi;
}
NormalizationParameters;

typedef struct tagAmplitudeParameters{
	REAL8 amp_coef_a;
	REAL8 amp_coef_b;
	REAL8 amp_coef_c;
	REAL8 amp_coef_a_mass;
	REAL8 amp_coef_b_mass;
	REAL8 amp_coef_c_mass;
	REAL8 amp_coef_a_eta;
	REAL8 amp_coef_b_eta;
	REAL8 amp_coef_c_eta;
	REAL8 amp_coef_a_chi;
	REAL8 amp_coef_b_chi;
	REAL8 amp_coef_c_chi;
}
AmplitudeParameters;

typedef struct tagPhaseParameters{
	REAL8 phase_coef_1_mass;
	REAL8 phase_coef_2_mass;
	REAL8 phase_coef_3_mass;
	REAL8 phase_coef_4_mass;
	REAL8 phase_coef_5_mass;
	REAL8 phase_coef_6_mass;
	REAL8 phase_coef_1_eta;
	REAL8 phase_coef_2_eta;
	REAL8 phase_coef_3_eta;
	REAL8 phase_coef_4_eta;
	REAL8 phase_coef_5_eta;
	REAL8 phase_coef_6_eta;
	REAL8 phase_coef_1_chi;
	REAL8 phase_coef_2_chi;
	REAL8 phase_coef_3_chi;
	REAL8 phase_coef_4_chi;
	REAL8 phase_coef_5_chi;
	REAL8 phase_coef_6_chi;
}
PhaseParameters;


static ConstantDefinitions * ChiPowList(
										const REAL8 chi
										){
	ConstantDefinitions *ChiPow	= (ConstantDefinitions *) malloc(sizeof(ConstantDefinitions));
	
	ChiPow->ChiPow217		=	pow(1. - chi,0.217);
	ChiPow->ChiPow26		=	pow(1. - chi,0.26);
	ChiPow->ChiPow783		=	pow(1. - chi,0.783);
	ChiPow->ChiPow152		=	pow(1. - chi,1.5230000000000001);
	ChiPow->ChiPow3			=	pow(1. - chi,0.3);
	ChiPow->ChiPow45		=	pow(1. - chi,0.45);
	ChiPow->ChiPow74		=	pow(1. - chi,0.74);
	ChiPow->ChiPow75		=	pow(1. - chi,0.75);
	ChiPow->ChiPow7			=	pow(1. - chi,0.7);
	
	return ChiPow;
	
}

/**
 * Compute the dimensionless Newtonian chirp time (Eq. XX of the paper LIGO-Doc-Number)
 */
static REAL8 ChirpTime_theta0(
							  const REAL8 mass,	/**< Total Mass of system */
							  const REAL8 eta,	/**< Symmetric mass ratio of system */
							  const REAL8 flow	/**< Lower Frequency Cut-off */
){
	REAL8	v0 =	cbrt(LAL_PI*mass*LAL_MTSUN_SI*flow);
	return 5.0/(128.0*eta*v0*v0*v0*v0*v0);
}

static REAL8 ChirpTime_theta3(
							  const REAL8 mass,	/**< Total Mass of system */
							  const REAL8 eta,	/**< Symmetric mass ratio of system */
							  const REAL8 flow	/**< Lower Frequency Cut-off */
){
	
	REAL8	v0	=	cbrt(mass*LAL_MTSUN_SI*flow);
	return (cbrt(LAL_PI)/(4.0*eta*v0*v0));
}


static REAL8 ChirpTime_theta3S(
							   const REAL8 mass,	/**< Total Mass of system */
							   const REAL8 chi, /**< Reduced spin parameter value	*/
							   const REAL8 flow	/**< Lower Frequency Cut-off */
){
	
	REAL8	v0	=	cbrt(mass*LAL_MTSUN_SI*flow);
	return (cbrt(LAL_PI)/(4.0*v0*v0))*(17022. - 9565.9*chi) ;
}

static REAL8 TransitionFrequencies_fmerg(
										 const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										 const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
										 const REAL8 chi,
										 const ConstantDefinitions *ChiPow
										 ){
	REAL8	M = mass*LAL_MTSUN_SI;
	const REAL8 chi_pow_217	=	ChiPow->ChiPow217;
	const REAL8 chi_pow_26	=	ChiPow->ChiPow26;
	
	return  ((0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*(chi*chi))*eta + (-0.01853136495384793 - 1.252422078178743*chi)*(eta*eta) - 2.25732638886097*(eta*eta*eta))/M);
	
	
}

static REAL8 TransitionFrequencies_fring(
										 const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										 const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
										 const REAL8 chi,
										 const ConstantDefinitions *ChiPow
										 ){
	REAL8	M = mass*LAL_MTSUN_SI;
	const REAL8	chi_pow_3 = ChiPow->ChiPow3	;
	
	return ((0.15915494309189532 - 0.10026761414789406*chi_pow_3 + (0.046759722280398854 - 0.039091637122231335*chi - 0.008305023240421283*(chi*chi))*eta + (-0.007925916165976387 + 0.054154060936448305*chi)*(eta*eta) + 0.7401341473545502*(eta*eta*eta))/M);
	
}

static REAL8 TransitionFrequencies_fcut(
										const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
										const REAL8 chi
										){
	REAL8	M = mass*LAL_MTSUN_SI;
	
	return ((0.10300507916907467 + (chi*chi)*(0.004285405997692373 + 0.046193130682991704*eta) - 0.042376595147648057*eta - 0.08638930311028079*(eta*eta) + 1.5667212597966178*(eta*eta*eta) + chi*(0.015576494280403797 - 0.026012283898939375*eta + 0.04070546824518315*(eta*eta)))/M);
}

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/

static REAL8 * MassPartialDerivativesWRTTheta0(
											   const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
											   const REAL8 theta3, /**< Theta3 component of Chirp-Time Co-ordinate system*/
											   const REAL8 theta3S,
											   const REAL8 flow	/**< Lower Frequency Cut-off */
){
	static  REAL8	coord_derivative_list[3];
	REAL8	mass, eta, mass_theta0, eta_theta0, chi_theta0;
	REAL8	theta0_pow_2_3 = cbrt(theta0*theta0);
	REAL8	theta3_pow_2_3 = cbrt(theta3*theta3);
	
	
	mass	=	0.0158314349441152767881061661265*(1./flow)*(theta3/theta0);
	eta		=	5.80732920322284746065959952653*(theta0_pow_2_3/(theta3_pow_2_3*theta3));
	
	mass_theta0		=	-mass/theta0;
	eta_theta0		=	(0.666666666666666666666666666667)*(eta/theta0);
	chi_theta0		=	0.000012000696668619612*(theta3S*theta3_pow_2_3/(theta0_pow_2_3*theta0));
	
	coord_derivative_list[0]		=		mass_theta0;
	coord_derivative_list[1]		=		eta_theta0;
	coord_derivative_list[2]		=		chi_theta0;
	
	return coord_derivative_list;
}

static REAL8 * MassPartialDerivativesWRTTheta3(
											   const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
											   const REAL8 theta3, /**< Theta3 component of Chirp-Time Co-ordinate system*/
											   const REAL8 theta3S,
											   const REAL8 flow	/**< Lower Frequency Cut-off */

){
	static  REAL8	coord_derivative_list[3];
	REAL8	eta, mass_theta3, eta_theta3, chi_theta3;
	REAL8	theta0_pow_2_3 = cbrt(theta0*theta0);
	REAL8	theta3_pow_2_3 = cbrt(theta3*theta3);
	
	eta		=	5.80732920322284746065959952653*(theta0_pow_2_3/(theta3_pow_2_3*theta3));
	mass_theta3		=	0.0158314349441152767881061661265*(1./flow)*(1/theta0);
	eta_theta3		=	-(5.0/3.0)*(eta/theta3);
	chi_theta3		=	-(0.666666666666666666666666666667)*0.00001800104500292942*(theta3S/(theta0_pow_2_3*cbrt(theta3)));
	
	coord_derivative_list[0]		=		mass_theta3;
	coord_derivative_list[1]		=		eta_theta3;
	coord_derivative_list[2]		=		chi_theta3;
	
	return coord_derivative_list;
}


static REAL8 * MassPartialDerivativesWRTTheta3S(
												const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
												const REAL8 theta3	/**< Theta3 component of Chirp-Time Co-ordinate system*/

){	static  REAL8	coord_derivative_list[3];
	
	REAL8	chi_theta3S;
	
	REAL8	theta0_pow_2_3 = cbrt(theta0*theta0);
	REAL8	theta3_pow_2_3 = cbrt(theta3*theta3);
	
	
	chi_theta3S		=	-0.00001800104500292942*(theta3_pow_2_3/(theta0_pow_2_3));
	
	coord_derivative_list[0]		=		0.0;
	coord_derivative_list[1]		=		0.0;
	coord_derivative_list[2]		=		chi_theta3S;
	
	return coord_derivative_list;
}

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/
static REAL8 CalculateDerivatives(
								  const REAL8 expr_mass,
								  const REAL8 expr_eta,
								  const REAL8 expr_chi,
								  const REAL8 *partial_der_list
								  ){
	REAL8	chain_rule;
	
	chain_rule	=	expr_mass*partial_der_list[0] + expr_eta*partial_der_list[1] + expr_chi*partial_der_list[2];
	return chain_rule;
}

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/
static NormalizationParameters * XLALSimIMRPhenomBNormalization(
																const REAL8 mass,
																const REAL8 eta,
																const REAL8 chi,
																const ConstantDefinitions *ChiPow
																){
	NormalizationParameters *Normalization 	= (NormalizationParameters *) malloc(sizeof(NormalizationParameters));
	
	REAL8	norm_merg, norm_merg_mass, norm_merg_eta, norm_merg_chi,pre_norm_ring, pre_norm_ring_mass, pre_norm_ring_eta, pre_norm_ring_chi;
	REAL8	M = mass*LAL_MTSUN_SI;
	
	
	const	REAL8	chi_pow_217		=	ChiPow->ChiPow217;
	const	REAL8	chi_pow_26		=	ChiPow->ChiPow26;
	const	REAL8	chi_pow_2		=	chi*chi;
	const	REAL8	chi_pow_783		=	ChiPow->ChiPow783;
	const	REAL8	chi_pow_152		=	ChiPow->ChiPow152;
	const	REAL8	chi_pow_3		=	ChiPow->ChiPow3;
	const	REAL8	chi_pow_45		=	ChiPow->ChiPow45;
	const	REAL8	chi_pow_74		=	ChiPow->ChiPow74;
	const	REAL8	chi_pow_75		=	ChiPow->ChiPow75;
	const	REAL8	chi_pow_7		=	ChiPow->ChiPow7;
	const	REAL8	eta_pow_2		=	eta*eta;
	const	REAL8	mass_pow_third	=	cbrt(M);
	
	const	REAL8	expr1			=	0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi -												0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta ;
	const	REAL8	expr1_pow_third	=	cbrt(expr1);
	const	REAL8	expr2			=	cbrt((1. - 4.455*chi_pow_217 + 3.521*chi_pow_26 + (0.64365 + 0.82696*chi - 0.27063*chi_pow_2)*eta + (-0.058218 - 3.9346*chi)*eta_pow_2 - 7.0916*									eta_pow_2*eta)*(1. - 4.455*chi_pow_217 + 3.521*chi_pow_26 + (0.64365 + 0.82696*chi - 0.27063*chi_pow_2)*eta + (-0.058218 - 3.9346*chi)*eta_pow_2 -									7.0916*eta_pow_2*eta));
	const	REAL8	expr3			=	1 + 2.1305418188357477*(-1.2990307279851514 + 1.*chi)*expr1_pow_third - 3.8938718645756443*(-0.9120806478268055 + 1.*chi)*expr1_pow_third*											expr1_pow_third;
	const	REAL8	expr4			=	cbrt((1 - 0.63*chi_pow_3)/(2.*M*LAL_PI) + (0.1469*eta - 0.12281*chi*eta - 0.026091*chi_pow_2*eta - 0.0249*eta_pow_2 + 0.17013*chi*eta_pow_2 + 2.3252								*eta_pow_2*eta)/(M*LAL_PI));
	const	REAL8	expr5			=	cbrt((0.15915494309189532 - 0.10026761414789406*chi_pow_3 + (0.046759722280398854 - 0.039091637122231335*chi - 0.008305023240421283*chi_pow_2)*eta + 								(-0.007925916165976387 + 0.054154060936448305*chi)*eta_pow_2 + 0.7401341473545502*eta_pow_2*eta)/M);
	const	REAL8	expr6			=	cbrt((1 - 4.455*chi_pow_217 + 3.521*chi_pow_26)/(M*LAL_PI) + (0.64365*eta + 0.82696*chi*eta - 0.27063*chi_pow_2*eta - 0.058218*eta_pow_2 - 3.9346*									chi*eta_pow_2 - 7.0916*eta_pow_2*eta)/(M*LAL_PI));
	
	norm_merg		=	(1 + ((-969 + 1804*eta)*expr2)/672. + 0.4961549999999999*chi*(-1.8409090909090908 + 1.*eta)*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta))/(expr3);
	
	norm_merg_mass	=	0.0;
	
	norm_merg_eta	=	(-((0.64365 - 0.27063*chi_pow_2 + chi*(0.82696 - 7.8692*eta) - 0.116436*eta - 21.2748*eta_pow_2)*(-1.8897 + 1.4547*chi - 5.317347306980865*(-0.9120806478268055 + 1.*chi)*expr1_pow_third)*(1 + ((-969 + 1804*eta)*expr2)/672. + 0.4961549999999999*chi*(-1.8409090909090908 + 1.*eta)*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta)))/(3.*expr2) + (expr3)*((451*pow(1. - 4.455*chi_pow_217 + 3.5209999999999995*chi_pow_26 + (0.6436500000000002 + 0.8269599999999999*chi - 0.27063*chi_pow_2)*eta + (-0.058218000000000006 - 3.9346*chi)*eta_pow_2 - 7.091600000000001*eta_pow_2*eta,0.6666666666666666))/168. - (25.997097497673163*(-0.5371396895787139 + 1.*eta)*(-0.030254103446330876 + 0.012720683625698009*chi_pow_2 + chi*(-0.03887040066181586 + 0.3698836181773742*eta) + 0.005472953917310622*eta + 1.*eta_pow_2))/expr1_pow_third + 0.49615499999999996*chi*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta) + 0.4961549999999999*chi*(-1.8409090909090908 + 1.*eta)*(-2.3783394302183796 + 1.*chi_pow_2 + 0.4302405498281787*eta + 78.61212725861878*eta_pow_2 + chi*(-3.055684883420168 + 29.0773380630381*eta))))/(expr3*expr3);
	
	norm_merg_chi	=	((expr3)*(chi*(3.375 - (11*eta)/6.)*(0.966735/chi_pow_783 - 0.91546/chi_pow_74 + (0.82696 - 0.54126*chi)*eta - 3.9346*eta_pow_2) + 13.001266666666666*(-1.8409090909090908 + 1.*eta)*(-0.14101190140447856 + 0.6282080207569519*chi_pow_217 - 0.49650290484516896*chi_pow_26 + (-0.09076231033899262 - 0.11661120198544758*chi + 0.038162050877094025*chi_pow_2)*eta + (0.008209430875965932 + 0.5548254272660613*chi)*eta_pow_2 + 1.*eta_pow_2*eta) - (4.807950242274655*(-0.5371396895787139 + 1.*eta)*(-0.24570096070756875*chi_pow_74 + 0.23266914044629697*chi_pow_783 + chi_pow_152*(-0.2101763838763788 + 0.13756417424897066*chi)*eta + 1.*chi_pow_152*eta_pow_2))/(chi_pow_152*(cbrt(-1.4180705429487876*chi_pow_217 + 1.120769109253127*chi_pow_26 + (0.2048801582421969 + 0.26322954347854755*chi - 0.08614420449791926*chi_pow_2)*eta + (-0.018531364953847926 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta + 1/LAL_PI)))) - ((1 + ((-969 + 1804*eta)*expr2)/672. + 0.4961549999999999*chi*(-1.8409090909090908 + 1.*eta)*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta))*(2.03451757159024 - 9.063775781434519*chi_pow_217 + 7.163536369569235*chi_pow_26 + (1.3095172349540583 + 1.6824646510022647*chi - 0.5506014903994666*chi_pow_2)*eta + (-0.11844554398284059 - 8.005012837178958*chi)*eta_pow_2 - 14.427984810689344*eta_pow_2*eta + (9.753542185587131*(-0.9120806478268055 + 1.*chi)*(-0.24570096070756875*chi_pow_74 + 0.23266914044629697*chi_pow_783 + chi_pow_152*(-0.2101763838763788 + 0.13756417424897066*chi)*eta + 1.*chi_pow_152*eta_pow_2)*expr1_pow_third)/chi_pow_152 - 11.681615593726933*cbrt(expr1)*expr1 + ((-1.8897 + 1.4547*chi)*(0.966735/chi_pow_783 - 0.91546/chi_pow_74 + (0.82696 - 0.54126*chi)*eta - 3.9346*eta_pow_2))/cbrt(LAL_PI*LAL_PI)))/(3.*expr1_pow_third*expr1_pow_third))/(expr3*expr3);
	
	
	
	pre_norm_ring	=	(5.5873167384795925*(-0.08711408460519897*chi_pow_45 + 0.05488187330127535*chi_pow_75 + (0.14279740748484215 + 0.012274722977210955*chi - 0.035131368039584644*chi_pow_2)*eta + (-0.6371872604362673 + 0.007028015889609032*chi)*eta_pow_2 + 1.*eta_pow_2*eta)*(-0.2568137922301612 - 0.5471525240001535*(-1.2990307279851518 + 1.*chi)*expr5*mass_pow_third + 1.*(-0.9120806478268055 + 1.*chi)*expr5*expr5*(mass_pow_third*mass_pow_third))*(cbrt((0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)/M)*cbrt((0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)/M)))/(expr5*expr5*M);
	
	pre_norm_ring_mass	=	((-2.415236754972846*(-1.5873015873015872 + 1.*chi_pow_3 + (-0.46634920634920635 + 0.38987301587301587*chi + 0.08282857142857143*chi_pow_2)*eta + (0.07904761904761905 - 0.5400952380952381*chi)*eta_pow_2 - 7.381587301587302*eta_pow_2*eta)*(-0.08711408460519895*chi_pow_45 + 0.05488187330127534*chi_pow_75 + (0.14279740748484215 + 0.012274722977210955*chi - 0.03513136803958464*chi_pow_2)*eta + (-0.6371872604362673 + 0.00702801588960903*chi)*eta_pow_2 + 1.*eta_pow_2*eta)*(-0.14101190140447856 + 0.6282080207569519*chi_pow_217 - 0.49650290484516896*chi_pow_26 + (-0.09076231033899262 - 0.11661120198544758*chi + 0.038162050877094025*chi_pow_2)*eta + (0.008209430875965932 + 0.5548254272660613*chi)*eta_pow_2 + 1.*eta_pow_2*eta)*(-0.2568137922301612 - 0.5471525240001535*(-1.2990307279851518 + 1.*chi)*expr5*mass_pow_third + 1.*(-0.9120806478268055 + 1.*chi)*expr5*expr5*(mass_pow_third*mass_pow_third)) + 3*(0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)*(0.07957747154594766*chi_pow_45 - 0.050133807073947025*chi_pow_75 + (-0.13044339135811742 - 0.011212784050710209*chi + 0.03209200272504978*chi_pow_2)*eta + (0.5820614578756795 - 0.006419992094440873*chi)*eta_pow_2 - 0.9134857113702426*eta_pow_2*eta)*(-7.829682347308064e-18*chi_pow_2*eta_pow_2 + chi*(-1.5659364694616128e-17 + eta*(3.914841173654032e-18 - 2.293457849610211e-17*expr5*mass_pow_third) + eta_pow_2*(5.8722617604810484e-18 + 1.7200933872076582e-17*expr5*mass_pow_third)) + eta_pow_2*(-9.78710293413508e-19 - 2.8668223120127636e-18*expr5*mass_pow_third) + chi_pow_2*chi*eta*(9.78710293413508e-19 + 2.8668223120127636e-18*expr5*mass_pow_third) + (4.586915699220422e-17 - 4.586915699220422e-17*chi_pow_3)*expr5*mass_pow_third)*expr5*mass_pow_third)*LAL_PI)/(6.*cbrt((0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)/M)*(expr5*expr5*(0.15915494309189532 - 0.10026761414789406*chi_pow_3 + (0.046759722280398854 - 0.039091637122231335*chi - 0.008305023240421283*chi_pow_2)*eta + (-0.007925916165976387 + 0.054154060936448305*chi)*eta_pow_2 + 0.7401341473545502*eta_pow_2*eta)/M)*M*M*M*M);
	
	pre_norm_ring_eta	=	((-0.0008425904274917489*(-4.064669708391192 + 1.*chi_pow_2 + chi*(-0.34939496131719894 - 0.4000991866693116*eta) + 36.27454870065463*eta - 85.39377107716723*eta_pow_2)*(-19.163696293741136 + 12.073128665056917*chi_pow_3 + (-5.630293971101146 + 4.7069870836686984*chi + 1.*chi_pow_2)*eta + (0.9543520754283086 - 6.52063930090836*chi)*eta_pow_2 - 89.11885324441378*eta_pow_2*eta)*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta)*(-0.2568137922301612 - 0.5471525240001535*(-1.2990307279851518 + 1.*chi)*expr5*mass_pow_third + 1.*(-0.9120806478268055 + 1.*chi)*expr5*expr5*(mass_pow_third*mass_pow_third)))/(M*M*M) + (0.0005617269516611657*(-5.630293971101146 + 1.*chi_pow_2 + chi*(4.7069870836686984 - 13.04127860181672*eta) + 1.9087041508566172*eta - 267.35655973324134*eta_pow_2)*(2.4796667327911126*chi_pow_45 - 1.562190041658401*chi_pow_75 + (-4.064669708391192 - 0.34939496131719894*chi + 1.*chi_pow_2)*eta + (18.137274350327313 - 0.2000495933346558*chi)*eta_pow_2 - 28.46459035905574*eta_pow_2*eta)*(-3.6950818460628905 + 16.461589624210177*chi_pow_217 - 13.010383179987437*chi_pow_26 + (-2.3783394302183796 - 3.055684883420168*chi + 1.*chi_pow_2)*eta + (0.21512027491408936 + 14.53866903151905*chi)*eta_pow_2 + 26.204042419539594*eta_pow_2*eta)*(-0.2568137922301612 - 0.5471525240001535*(-1.2990307279851518 + 1.*chi)*expr5*mass_pow_third + 1.*(-0.9120806478268055 + 1.*chi)*expr5*expr5*(mass_pow_third*mass_pow_third)))/(M*M*M) - (0.0005617269516611657*(-19.163696293741136 + 12.073128665056917*chi_pow_3 + (-5.630293971101146 + 4.7069870836686984*chi + 1.*chi_pow_2)*eta + (0.9543520754283086 - 6.52063930090836*chi)*eta_pow_2 - 89.11885324441378*eta_pow_2*eta)*(2.4796667327911126*chi_pow_45 - 1.562190041658401*chi_pow_75 + (-4.064669708391192 - 0.34939496131719894*chi + 1.*chi_pow_2)*eta + (18.137274350327313 - 0.2000495933346558*chi)*eta_pow_2 - 28.46459035905574*eta_pow_2*eta)*(-2.3783394302183796 + 1.*chi_pow_2 + 0.4302405498281787*eta + 78.61212725861878*eta_pow_2 + chi*(-3.055684883420168 + 29.0773380630381*eta))*(-0.2568137922301612 - 0.5471525240001535*(-1.2990307279851518 + 1.*chi)*expr5*mass_pow_third + 1.*(-0.9120806478268055 + 1.*chi)*expr5*expr5*(mass_pow_third*mass_pow_third)))/(M*M*M) + ((0.1469 - 0.026091*chi_pow_2 + chi*(-0.12281 + 0.34026*eta) - 0.0498*eta + 6.9756*eta_pow_2)*(0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)*(0.07957747154594766*chi_pow_45 - 0.050133807073947025*chi_pow_75 + (-0.13044339135811742 - 0.011212784050710209*chi + 0.03209200272504978*chi_pow_2)*eta + (0.5820614578756795 - 0.006419992094440873*chi)*eta_pow_2 - 0.9134857113702426*eta_pow_2*eta)*(0.15915494309189532 - 0.10026761414789406*chi_pow_3 + (0.046759722280398854 - 0.039091637122231335*chi - 0.008305023240421283*chi_pow_2)*eta + (-0.007925916165976387 + 0.054154060936448305*chi)*eta_pow_2 + 0.7401341473545502*eta_pow_2*eta)*(-1.8897 + 1.4547*chi - 5.317347306980865*(-0.9120806478268055 + 1.*chi)*expr5*mass_pow_third)*LAL_PI)/(pow((0.5 - 0.315*chi_pow_3 + (0.1469 - 0.12281*chi - 0.026091*chi_pow_2)*eta + (-0.0249 + 0.17013*chi)*eta_pow_2 + 2.3252*eta_pow_2*eta)/M,0.6666666666666666)*pow(M,3.6666666666666665)))/(6.*cbrt((0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_2*eta)/M)*(expr5*expr5*(0.15915494309189532 - 0.10026761414789406*chi_pow_3 + (0.046759722280398854 - 0.039091637122231335*chi - 0.008305023240421283*chi_pow_2)*eta + (-0.007925916165976387 + 0.054154060936448305*chi)*eta_pow_2 + 0.7401341473545502*eta_pow_2*eta)/M));
	
	pre_norm_ring_chi	=	((expr6*expr6)*(((1 - 0.63*chi_pow_3)*chi_pow_45)/(4.*M*LAL_PI) + (-0.4098*eta - 0.035226*chi*eta + 0.10082*chi_pow_2*eta + 1.8286*eta_pow_2 - 0.020169*chi*eta_pow_2 - 2.8698*eta_pow_2*eta)/(M*LAL_PI))*(2.1305418188357477*mass_pow_third*expr4 - 3.8938718645756443*(mass_pow_third*mass_pow_third)*expr4*expr4 + ((-1.8897 + 1.4547*chi)*mass_pow_third*(0.03008028424436822/(chi_pow_7*M) + (-0.12281*eta - 0.052182*chi*eta + 0.17013*eta_pow_2)/(M*LAL_PI))*cbrt(LAL_PI))/(3.*expr4*expr4) + (2*(1.6557 - 1.8153*chi)*(mass_pow_third*mass_pow_third)*(0.03008028424436822/(chi_pow_7*M) + (-0.12281*eta - 0.052182*chi*eta + 0.17013*eta_pow_2)/(M*LAL_PI))*cbrt(LAL_PI*LAL_PI))/(3.*expr4))*LAL_PI)/(2.*expr4*expr4) - ((0.03008028424436822/(chi_pow_7*M) + (-0.12281*eta - 0.052182*chi*eta + 0.17013*eta_pow_2)/(M*LAL_PI))*(expr6*expr6)*(((1 - 0.63*chi_pow_3)*chi_pow_45)/(4.*M*LAL_PI) + (-0.4098*eta - 0.035226*chi*eta + 0.10082*chi_pow_2*eta + 1.8286*eta_pow_2 - 0.020169*chi*eta_pow_2 - 2.8698*eta_pow_2*eta)/(M*LAL_PI))*(1 + (-1.8897 + 1.4547*chi)*mass_pow_third*expr4*cbrt(LAL_PI) + (1.6557 - 1.8153*chi)*(mass_pow_third*mass_pow_third)*expr4*expr4*cbrt(LAL_PI*LAL_PI))*LAL_PI)/(3.*pow((1 - 0.63*chi_pow_3)/(2.*M*LAL_PI) + (0.1469*eta - 0.12281*chi*eta - 0.026091*chi_pow_2*eta - 0.0249*eta_pow_2 + 0.17013*chi*eta_pow_2 + 2.3252*eta_pow_2*eta)/(M*LAL_PI),1.6666666666666667)) + (((-0.03580986219567645*(1 - 0.63*chi_pow_3))/(pow(1 - chi,0.55)*M) + 0.01504014212218411/(pow(1 - chi,0.24999999999999994)*M) + (-0.035226*eta + 0.20164*chi*eta - 0.020169*eta_pow_2)/(M*LAL_PI))*(expr6*expr6)*(1 + (-1.8897 + 1.4547*chi)*mass_pow_third*expr4*cbrt(LAL_PI) + (1.6557 - 1.8153*chi)*(mass_pow_third*mass_pow_third)*expr4*expr4*cbrt(LAL_PI*LAL_PI))*LAL_PI)/(2.*expr4*expr4) + (((0.966735/chi_pow_783 - 0.91546/chi_pow_74)/(M*LAL_PI) + (0.82696*eta - 0.54126*chi*eta - 3.9346*eta_pow_2)/(M*LAL_PI))*(((1 - 0.63*chi_pow_3)*chi_pow_45)/(4.*M*LAL_PI) + (-0.4098*eta - 0.035226*chi*eta + 0.10082*chi_pow_2*eta + 1.8286*eta_pow_2 - 0.020169*chi*eta_pow_2 - 2.8698*eta_pow_2*eta)/(M*LAL_PI))*(1 + (-1.8897 + 1.4547*chi)*mass_pow_third*expr4*cbrt(LAL_PI) + (1.6557 - 1.8153*chi)*(mass_pow_third*mass_pow_third)*expr4*expr4*cbrt(LAL_PI*LAL_PI))*LAL_PI)/(3.*expr6*expr4*expr4);
	
	Normalization->norm_merg				=	norm_merg;
	Normalization->norm_merg_mass			=	norm_merg_mass;
	Normalization->norm_merg_eta			=	norm_merg_eta;
	Normalization->norm_merg_chi			=	norm_merg_chi;
	Normalization->norm_ring				=	norm_merg*pre_norm_ring;
	Normalization->norm_ring_mass			=	norm_merg*pre_norm_ring_mass + norm_merg_mass*pre_norm_ring;
	Normalization->norm_ring_eta			=	norm_merg*pre_norm_ring_eta + norm_merg_eta*pre_norm_ring;
	Normalization->norm_ring_chi			=	norm_merg*pre_norm_ring_chi + norm_merg_chi*pre_norm_ring;
	
	return Normalization;
}


/**
 * Compute the co-efficients of inspiral phase of IMRPhenomB wavefrom amplitude in mass Co-ordinates.
 */


static AmplitudeParameters *XLALSimIMRPhenomBAmplitude_Inspiral(
																const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
																const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
																const REAL8 chi
																){
	
	AmplitudeParameters *Amplitude  = (AmplitudeParameters *) malloc(sizeof(AmplitudeParameters));
	REAL8 M = mass*LAL_MTSUN_SI;
	REAL8 eta_pow_half = sqrt(eta);
	REAL8 mass_pow_half	= sqrt(M);
	REAL8 mass_pow_5_6 = M/sqrt(cbrt(M));
	
	/* Inspiral Phase of waveform */
	Amplitude->amp_coef_a = (0.212787510139663399157865621657*eta_pow_half*mass_pow_5_6);
	Amplitude->amp_coef_b = eta_pow_half*(-0.65816363866878219734598675537 + 1.2253118721965769700847885518*eta)*M*sqrt(M);
	Amplitude->amp_coef_c = 0.212787510139663399157865621657*eta_pow_half*mass_pow_5_6*(10.6028752058655521798114214186*chi*M - 5.75958653158128760384817953601*chi*eta*M) ;
	
	/* Derivative of inspiral phase w.r.t theta0 co-ordinate */
	Amplitude->amp_coef_a_mass	= (0.177322925116386165964888018048*eta_pow_half)/sqrt(cbrt(M));
	Amplitude->amp_coef_b_mass	= -0.987245458003173296018980133062*eta_pow_half*mass_pow_half + 1.8379678082948654551271828277*eta_pow_half*eta*mass_pow_half ;
	Amplitude->amp_coef_c_mass	= chi*(4.13629226152578669789413343455 - 2.2468748087300569716955786558*eta)*eta_pow_half*mass_pow_5_6;
	
	/* Derivative of inspiral phase w.r.t theta3 co-ordinate */
	Amplitude->amp_coef_a_eta = (0.106393755069831699578932810829*mass_pow_5_6)/eta_pow_half;
	Amplitude->amp_coef_b_eta = ((-0.3290818193343910986729933777 + 1.837967808294865455127182828*eta)*M*sqrt(M))/eta_pow_half;
	Amplitude->amp_coef_c_eta = (chi*(1.128079707688850917607490937 - 1.838352116233682976841837082*eta)*(mass_pow_5_6*M))/eta_pow_half;
	
	/* Derivative of inspiral phase w.r.t theta3S co-ordinate */
	Amplitude->amp_coef_a_chi = 0.;
	Amplitude->amp_coef_b_chi = 0.;
	Amplitude->amp_coef_c_chi = -1.2255680774891219845612247213*eta_pow_half*(-1.8409090909090909090909090909 + 1.*eta)*(mass_pow_5_6*M);
	
	return Amplitude;
}

/**
 * Compute the co-efficients of merger phase of IMRPhenomB wavefrom amplitude in mass Co-ordinates.
 */

static AmplitudeParameters *XLALSimIMRPhenomBAmplitude_Merger(
															  const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
															  const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
															  const REAL8 chi,
															  const NormalizationParameters *Normalization,
															  const ConstantDefinitions *ChiPow
															  ){
	AmplitudeParameters *Amplitude  = (AmplitudeParameters *) malloc(sizeof(AmplitudeParameters));
	REAL8	M = mass*LAL_MTSUN_SI;
	const	REAL8	chi_pow_217		=	ChiPow->ChiPow217;
	const	REAL8	chi_pow_26		=	ChiPow->ChiPow26;
	const	REAL8	chi_pow_2		=	chi*chi;
	const	REAL8	chi_pow_783		=	ChiPow->ChiPow783;
	const	REAL8	chi_pow_152		=	ChiPow->ChiPow152;
	const	REAL8	chi_pow_74		=	ChiPow->ChiPow74;
	const 	REAL8	eta_pow_half	=	sqrt(eta);
	const 	REAL8	eta_pow_2		=	eta*eta;
	const 	REAL8	eta_pow_3		=	eta_pow_2*eta;
	
	const	REAL8	expr1			=	sqrt(0.31830988618379075 - 1.4180705429487876*chi_pow_217 + 1.1207691092531271*chi_pow_26 + (0.2048801582421969 + 0.2632295434785476*chi - 0.08614420449791928*chi_pow_2)*eta + (-0.01853136495384793 - 1.252422078178743*chi)*eta_pow_2 - 2.25732638886097*eta_pow_3);
	const	REAL8	expr2		=	(-0.22446689113355783 + 1.*chi_pow_217 - 0.7903479236812571*chi_pow_26 + (-0.14447811447811448 - 0.185625140291807*chi + 0.060747474747474745*chi_pow_2)*eta + (0.013068013468013468 + 0.8831874298540965*chi)*eta_pow_2 + 1.5918294051627384*eta_pow_3);
	
	REAL8 merg_coef1, merg_coef2, merg_coef3;
	REAL8 merg_coef1_mass, merg_coef2_mass, merg_coef3_mass;	/*  */
	REAL8 merg_coef1_eta, merg_coef2_eta, merg_coef3_eta;	/*  */
	REAL8 merg_coef1_chi, merg_coef2_chi, merg_coef3_chi;	/*  */
	
	/******************** Merger phase of waveform	******************************/
	
	merg_coef1	=	(0.21278751013966343*eta_pow_half*(cbrt(M)*M))/expr1;
	merg_coef2	=	merg_coef1*(2.1305418188357477*(-1.2990307279851514 + 1.*chi)*cbrt(M));
	merg_coef3	=	merg_coef1*(-3.8938718645756443*(-0.9120806478268055 + 1.*chi)*cbrt(M*M));
	
	Amplitude->amp_coef_a =	Normalization->norm_merg*merg_coef1;
	Amplitude->amp_coef_b =	Normalization->norm_merg*merg_coef2;
	Amplitude->amp_coef_c =	Normalization->norm_merg*merg_coef3;
	
	/******************** Derivative of merger phase w.r.t theta0 co-ordinate	******************************/
	merg_coef1_mass	=	merg_coef1*(4.0/(3.0*M));
	merg_coef2_mass	=	merg_coef1_mass*(2.6631772735446857*(-1.2990307279851514 + 1.*chi)*cbrt(M));
	merg_coef3_mass	=	merg_coef1_mass*(-5.840807796863467*(-0.9120806478268056 + 1.*chi)*cbrt(M*M));
	
	Amplitude->amp_coef_a_mass	=	Normalization->norm_merg*merg_coef1_mass + Normalization->norm_merg_mass*merg_coef1;
	Amplitude->amp_coef_b_mass	=	Normalization->norm_merg*merg_coef2_mass + Normalization->norm_merg_mass*merg_coef2;
	Amplitude->amp_coef_c_mass	=	Normalization->norm_merg*merg_coef3_mass + Normalization->norm_merg_mass*merg_coef3;
	
	/******************** Derivative of merger phase w.r.t theta3 co-ordinate	******************************/
	merg_coef1_eta	=	((-0.02388187543655033 + 0.10639375506983173*chi_pow_217 - 0.08408808341209371*chi_pow_26 + (2.1565808261191165e-18 + 4.313161652238233e-18*chi - 3.234871239178675e-18*chi_pow_2)*eta + (-0.0013903550241650875 - 0.09396562709265095*chi)*eta_pow_2 - 0.3387214156916807*eta_pow_3)*(cbrt(M)*M))/(eta_pow_half*expr1*expr2);
	merg_coef2_eta	=	((0.06609641677529171 - 0.29445953673392455*chi_pow_217 + 0.23272548346580207*chi_pow_26 + chi_pow_2*(8.873370157277154e-18 - 0.20019769805401816*eta)*eta + 0.0038480011918239324*eta_pow_2 + 0.9374586984073172*eta_pow_3 + chi*(-0.050881334329796725 + 0.22667634443924436*chi_pow_217 - 0.17915317817521426*chi_pow_26 - 4.436685078638577e-18*eta + 0.2571007519220506*eta_pow_2 - 0.7216601410663728*eta_pow_3))*pow(M,1.6666666666666667))/(eta_pow_half*expr1*expr2);
	merg_coef3_eta	=	((-0.08481708178650428 + 0.37786009935887654*chi_pow_217 - 0.29864094497028154*chi_pow_26 - 0.004937880867446704*eta_pow_2 + 0.3658901115732804*chi_pow_2*eta_pow_2 - 1.202977634394347*eta_pow_3 + chi*(0.09299296283568353 - 0.4142836494329701*chi_pow_217 + 0.32742822214444167*chi_pow_26 - 8.108680640276699e-18*eta - 0.3283074256868117*eta_pow_2 + 1.3189377904910664*eta_pow_3))*(M*M))/(eta_pow_half*expr1*expr2);
	
	Amplitude->amp_coef_a_eta	=	Normalization->norm_merg*merg_coef1_eta + Normalization->norm_merg_eta*merg_coef1;
	Amplitude->amp_coef_b_eta	=	Normalization->norm_merg*merg_coef2_eta + Normalization->norm_merg_eta*merg_coef2;
	Amplitude->amp_coef_c_eta	=	Normalization->norm_merg*merg_coef3_eta + Normalization->norm_merg_eta*merg_coef3;
	
	/******************** Derivative of merger phase w.r.t theta3S co-ordinate	******************************/
	merg_coef1_chi	=	(eta_pow_half*(0.023087444850153484*chi_pow_74 - 0.021862901687144366*chi_pow_783 + chi_pow_152*(0.019749355711009664 - 0.012926303898787232*chi)*eta - 0.09396562709265092*chi_pow_152*eta_pow_2)*(cbrt(M)*M))/(chi_pow_152*expr1*expr2);
	merg_coef2_chi	=	(eta_pow_half*(-0.06389771947126163*chi_pow_74 + 0.06050862570110856*chi_pow_783 - 0.10176266865959344*chi_pow_152 + 0.4533526888784887*pow(1. - 1.*chi,1.74) - 0.35830635635042846*pow(1. - 1.*chi,1.783) - 0.12015863449924255*chi_pow_152*eta - 4.436685078638576e-18*chi_pow_152*chi_pow_2*eta + 0.2659873804880869*chi_pow_152*eta_pow_2 + 0.7216601410663727*chi_pow_152*eta_pow_3 + chi*(0.049188766743316036*chi_pow_74 - 0.0465798263255557*chi_pow_783 - 0.006301481693574297*chi_pow_152*eta + 0.2001976980540181*chi_pow_152*eta_pow_2))*pow(M,1.6666666666666667))/(chi_pow_152*expr1*expr2);
	merg_coef3_chi	=	(eta_pow_half*(0.08199564156087623*chi_pow_74 - 0.07764664569227323*chi_pow_783 + 0.1859859256713671*chi_pow_152 - 0.8285672988659402*pow(1. - 1.*chi,1.74) + 0.6548564442888833*pow(1. - 1.*chi,1.783) + 0.189850175012543*chi_pow_152*eta - 0.34454901861791526*chi_pow_152*eta_pow_2 - 1.3189377904910666*chi_pow_152*eta_pow_3 + chi*(-0.08989955192695454*chi_pow_74 + 0.08513133775755487*chi_pow_783 + 0.03099336685883355*chi_pow_152*eta - 0.3658901115732804*chi_pow_152*eta_pow_2))*(M*M))/(chi_pow_152*expr1*expr2);
	
	Amplitude->amp_coef_a_chi	=	Normalization->norm_merg*merg_coef1_chi + Normalization->norm_merg_chi*merg_coef1;
	Amplitude->amp_coef_b_chi	=	Normalization->norm_merg*merg_coef2_chi + Normalization->norm_merg_chi*merg_coef2;
	Amplitude->amp_coef_c_chi	=	Normalization->norm_merg*merg_coef3_chi + Normalization->norm_merg_chi*merg_coef3;
	
	return Amplitude;
}


/**
 * Compute the co-efficients of ringdown phase of IMRPhenomB wavefrom amplitude in mass Co-ordinates.
 */

static AmplitudeParameters *XLALSimIMRPhenomBAmplitude_Ringdown(
																const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
																const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
																const REAL8 chi,
																const REAL8 fRing,
																const ConstantDefinitions *ChiPow
																){
	AmplitudeParameters *Amplitude  = (AmplitudeParameters *) malloc(sizeof(AmplitudeParameters));
	REAL8	M = mass*LAL_MTSUN_SI;
	const	REAL8	chi_pow_217		=	ChiPow->ChiPow217;
	const	REAL8	chi_pow_26		=	ChiPow->ChiPow26;
	const	REAL8	chi_pow_2		=	chi*chi;
	const	REAL8	chi_pow_783		=	ChiPow->ChiPow783;
	const	REAL8	chi_pow_75		=	ChiPow->ChiPow75;
	const	REAL8	chi_pow_7		=	ChiPow->ChiPow7;
	REAL8	f_ring, f_ring_mass, f_ring_eta, f_ring_chi;
	REAL8	sigma, sigma_mass, sigma_eta, sigma_chi;
	REAL8	amp_const, amp_const_mass, amp_const_eta, amp_const_chi;
	
	const REAL8	expr1	=	0.3183098861837907 - 1.4180705429487876*chi_pow_217 + 1.120769109253127*chi_pow_26 + eta*(0.2048801582421969 - 0.08614420449791926*chi_pow_2 + chi*(0.26322954347854755 - 1.252422078178743*eta) + (-0.018531364953847926 - 2.25732638886097*eta)*eta);
	
	/*******************************	Define the terms required for the lorentzian and it's derivatives	*************************************/
	f_ring			=	fRing;
	f_ring_mass		=	-f_ring/M;
	f_ring_eta		=	(0.04675972228039886 - 0.008305023240421283*chi_pow_2 + chi*(-0.03909163712223134 + 0.10830812187289664*eta) + eta*(-0.015851832331952777 + 2.2204024420636506*eta))/M;
	f_ring_chi		=	(0.03008028424436822/chi_pow_7 + (-0.039091637122231335 - 0.016610046480842567*chi + 0.05415406093644831*eta)*eta)/M;
	
	sigma			=	(0.07957747154594766*pow(1. - chi,0.45) - 0.050133807073947025*chi_pow_75 + (-0.13044339135811742 - 0.011212784050710209*chi + 0.03209200272504978*chi_pow_2)*eta + 				(0.5820614578756795 - 0.006419992094440873*chi)*(eta*eta) - 0.9134857113702426*(eta*eta*eta))/M;
	sigma_mass		=	-sigma/M;
	sigma_eta		=	(-0.13044339135811742 + 0.03209200272504978*chi_pow_2 + chi*(-0.01121278405071021 - 0.012839984188881748*eta) + 1.1641229157513593*eta - 2.7404571341107276*(eta*eta))/				M;
	sigma_chi		=	(-0.03580986219567645/pow(1 - chi,0.55) + 0.03760035530546027/sqrt(sqrt(1. - chi)) + (-0.01121278405071021 + 0.06418400545009956*chi)*eta - 0.006419992094440875*(eta*				eta))/M;
	
	amp_const		=	(0.21278751013966343*sqrt(eta)*(M*M))/(sqrt(cbrt(expr1))*expr1);
	amp_const_mass	=	2*amp_const/M;
	amp_const_eta	=	((0.016841105370463082 - 0.07502712442541304*chi_pow_217 + 0.05929753200940051*chi_pow_26 + eta*(-0.014453036628931423 + 0.0060769444618778996*chi_pow_2 + chi*(-0.0185692273295442 + 0.24296438169895485*eta) + (0.0035950033990112726 + 0.7165822970710561*eta)*eta))*(M*M))/(sqrt(eta)*sqrt(cbrt(expr1))*(-0.22446689113355778 + 1.*chi_pow_217 - 0.7903479236812568*chi_pow_26 + eta*(-0.14447811447811448 + 0.06074747474747474*chi_pow_2 + chi*(-0.18562514029180693 + 0.8831874298540965*eta) + eta*(0.013068013468013466 + 1.5918294051627384*eta)))*(-0.22446689113355778 + 1.*chi_pow_217 - 0.7903479236812568*chi_pow_26 + eta*(-0.14447811447811448 + 0.06074747474747474*chi_pow_2 + chi*(-0.18562514029180693 + 0.8831874298540965*eta) + eta*(0.013068013468013466 + 1.5918294051627384*eta))));
	amp_const_chi	=	(-0.2482520951629407*sqrt(eta)*(0.3077213078198869/chi_pow_783 - 0.29139996840581306/pow(1. - chi,0.74) + (0.26322954347854755 - 0.17228840899583853*chi - 1.252422078178743*eta)*eta)*(M*M))/(sqrt(cbrt(expr1))*expr1*expr1);
	
	Amplitude->amp_coef_a		=	f_ring;
	Amplitude->amp_coef_a_mass	=	f_ring_mass;
	Amplitude->amp_coef_a_eta	=	f_ring_eta;
	Amplitude->amp_coef_a_chi	=	f_ring_chi;
	Amplitude->amp_coef_b		=	sigma;
	Amplitude->amp_coef_b_mass	=	sigma_mass;
	Amplitude->amp_coef_b_eta	=	sigma_eta;
	Amplitude->amp_coef_b_chi	=	sigma_chi;
	Amplitude->amp_coef_c		=	amp_const;
	Amplitude->amp_coef_c_mass	=	amp_const_mass;
	Amplitude->amp_coef_c_eta	=	amp_const_eta;
	Amplitude->amp_coef_c_chi	=	amp_const_chi;
	
	return Amplitude;
}

/**
 * Compute the co-efficients of phase of IMRPhenomB wavefrom in mass Co-ordinates.
 */


static PhaseParameters *XLALSimIMRPhenomBPhase(
											   const REAL8 mass,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
											   const REAL8 eta,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
											   const REAL8 chi
											   ){
	PhaseParameters *Phase   = (PhaseParameters *) malloc(sizeof(PhaseParameters));
	REAL8	M = mass*LAL_MTSUN_SI;
	const	REAL8	 mass_pow_third		=	cbrt(M);
	const	REAL8	 mass_pow_2_third	=	mass_pow_third*mass_pow_third;
	const	REAL8	 eta_pow_2			=	eta*eta;
	const	REAL8	 eta_pow_3			=	eta_pow_2*eta;
	const	REAL8	 chi_pow_2			=	chi*chi;
	
	/******************************		Derivative of Phase w.r.t mass co-ordinate		*************************************/
	
	Phase->phase_coef_1_mass		=	-0.005796647796902313/(eta*pow(M,2.6666666666666665));
	Phase->phase_coef_2_mass		=	(-0.007460387957432594*(4.914021164021164 - 920.91*eta + 492.13*chi*eta + 135.03*chi_pow_2*eta + 6741.9*eta_pow_2 - 1053.4*chi*eta_pow_2 - 											13396.999999999998*eta_pow_3))/(eta*M*M);
	Phase->phase_coef_3_mass		=	(-0.007284282453678307*(-50.26548245743669 + 37.666666666666664*chi + 17022.*eta - 9565.9*chi*eta - 2182.1*chi_pow_2*eta - 121370.*eta_pow_2 + 20752.								*chi*eta_pow_2 + 238590.*eta_pow_3))/(eta*pow(M,1.6666666666666667));
	Phase->phase_coef_4_mass		=	(-0.005334250494181998*(30.103152950995213 - 50.625*chi_pow_2 - 125440.*eta + 75066.*chi*eta + 13382.*chi_pow_2*eta + 873540.*eta_pow_2 - 165730.*chi								*eta_pow_2 - 1.6936e6*eta_pow_3))/(eta*(mass_pow_third*M));
	Phase->phase_coef_5_mass		=	(0.0114421241215744*(-889770.*eta + 631020.*chi*eta + 50676.*chi_pow_2*eta + 5.9808e6*eta_pow_2 - 1.4148e6*chi*eta_pow_2 - 1.1279999999999998e7*									eta_pow_3))/(eta*mass_pow_2_third);
	Phase->phase_coef_6_mass		=	(0.03351608432985977*(869600.*eta - 670980.*chi*eta - 30082.*chi_pow_2*eta - 5.8379e6*eta_pow_2 + 1.5145e6*chi*eta_pow_2 + 1.0891e7*eta_pow_3))/(eta*								mass_pow_third);
	
	/******************************		Derivative of Phase w.r.t eta co-ordinate		*************************************/
	
	Phase->phase_coef_1_eta			=	-0.0034779886781413872/(eta_pow_2*pow(M,1.6666666666666667)) ;
	Phase->phase_coef_2_eta			=	50.2971895702148/M - (7.8587726743594954*chi)/M - 0.03666050431463239/(eta_pow_2*M) - (199.89363493144887*eta)/M;
	Phase->phase_coef_3_eta			=	(0.549221957835571 - 1326.1400421044043*eta_pow_2 + 5213.870851869322*eta_pow_3 + chi*(-0.41156195863282435 + 226.74514421809837*eta_pow_2))/										(eta_pow_2*mass_pow_2_third);
	Phase->phase_coef_4_eta			=	(-0.4817332755158474 + 0.8101392938038909*chi_pow_2 + 13979.043530063223*eta_pow_2 - 2652.136003202347*chi*eta_pow_2 - 54204.51982167979*eta_pow_3)/								(eta_pow_2*mass_pow_third);
	Phase->phase_coef_5_eta			=	(205299.1678389365 - 48564.95162161038*chi - 774402.960548155*eta)*mass_pow_third;
	Phase->phase_coef_6_eta			=	(-293495.32306393253 + 76140.16457635893*chi + 1.0950710233095083e6*eta)*mass_pow_2_third;
	
	/******************************		Derivative of Phase w.r.t chi co-ordinate		*************************************/
	
	Phase->phase_coef_1_chi			=	0. ;
	Phase->phase_coef_2_chi			=	(3.6714807254913024 + 2.014752371784246*chi - 7.8587726743594954*eta)/M;
	Phase->phase_coef_3_chi			=	(0.41156195863282435 + (-104.52107628546197 - 47.6850982265143*chi)*eta + 226.74514421809835*eta_pow_2)/(eta*mass_pow_2_third);
	Phase->phase_coef_4_chi			=	((1201.2625427887974 - 2652.136003202347*eta)*eta + chi*(-1.6202785876077814 + 428.2976406788609*eta))/(eta*mass_pow_third);
	Phase->phase_coef_5_chi			=	3479.0464919094256*(6.226024153445418 + 1.*chi - 13.959270660667771*eta)*mass_pow_third;
	Phase->phase_coef_6_chi			=	-3024.692546432525*(11.152516455022937 + 1.*chi - 25.17286084701815*eta)*mass_pow_2_third;
	
	return Phase;
}

/**
 * Function to compute the metric elements using waveform derivatives
 */

static REAL8 MetricCoeffs(REAL8Vector *Amp, REAL8Vector *dPsii, REAL8Vector *dPsij,
						  REAL8Vector *dAi, REAL8Vector*dAj, REAL8Vector *Sh, REAL8 hSqr, REAL8 df) {
	size_t k = Amp->length;
	REAL8 gij   = 0.;
	for (;k--;) {
		gij += df*(Amp->data[k]*Amp->data[k]*dPsii->data[k]*dPsij->data[k]
				   + dAi->data[k]*dAj->data[k])/(2.0*Sh->data[k]*hSqr);
	}
	return gij;
}


/**
 * Compute the three dimensional template-space metric of IMRPhenomB wavefrom in New Co-ordinates.
 */

int XLALSimIMRPhenomBMetricTheta0Theta3Theta3S(
											   REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
											   REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
											   REAL8 *gamma02,  /**< template metric coeff. 01/10 PN Chirp Time */
											   REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
											   REAL8 *gamma12,  /**< template metric coeff. 01/10 PN Chirp Time */
											   REAL8 *gamma22,  /**< template metric coeff. 01/10 PN Chirp Time */
											   const REAL8 Mass,     /**< Total Mass of the system */
											   const REAL8 eta,    /**< Symmetric mass ratio */
											   const REAL8 chi,	   /** Reduced Spin Parameter of the system **/
											   const REAL8 flow,   /**< low-frequency cutoff (Hz) */
											   const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
) {
	REAL8Vector *Amp=NULL, *dATheta0=NULL, *dATheta3=NULL, *dATheta3S=NULL;
	REAL8Vector *dAT0=NULL, *dAPhi=NULL, *dPhaseTheta0=NULL;
	REAL8Vector *dPhaseTheta3=NULL, *dPhaseTheta3S=NULL, *dPhaseT0=NULL, *dPhasePhi=NULL;
	REAL8 *der_list_theta0, *der_list_theta3, *der_list_theta3S;
	
	ConstantDefinitions *ChiPowers =	ChiPowList(chi);
	
	/* compute the chirp-time co-ordinates */
	const REAL8 theta0	=	ChirpTime_theta0(Mass,eta, flow);
	const REAL8 theta3	=	ChirpTime_theta3(Mass,eta, flow);
	const REAL8	theta3S	=	ChirpTime_theta3S(Mass, chi, flow);
	
	/* Compute the transition frequencies */
	
	const REAL8 fMerg  = TransitionFrequencies_fmerg(Mass,eta, chi, ChiPowers);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
	const REAL8 fRing  = TransitionFrequencies_fring(Mass,eta, chi, ChiPowers);	/**Frequency at which merger part transitions to ringdown part of the waveform*/
	const REAL8 fCut   = TransitionFrequencies_fcut(Mass,eta, chi);	/**Frequency at which ringdown part of the waveform is terminated*/
	
	/* Compute the list of partial derivatives */
	der_list_theta0		=	MassPartialDerivativesWRTTheta0(theta0, theta3, theta3S, flow);
	der_list_theta3		=	MassPartialDerivativesWRTTheta3(theta0, theta3, theta3S, flow);
	der_list_theta3S	=	MassPartialDerivativesWRTTheta3S(theta0, theta3);
	
	/*Compute the normalizations and their derivatives*/
	NormalizationParameters *Normalization		=		XLALSimIMRPhenomBNormalization(Mass,eta, chi, ChiPowers);
	const REAL8 norm_ring			=		Normalization->norm_ring;
	const REAL8	norm_ring_mass		=		Normalization->norm_ring_mass;
	const REAL8	norm_ring_eta		=		Normalization->norm_ring_eta;
	const REAL8	norm_ring_chi		=		Normalization->norm_ring_chi;
	
	AmplitudeParameters *Inspiral	=	XLALSimIMRPhenomBAmplitude_Inspiral(Mass,eta, chi);
	AmplitudeParameters *Merger		=	XLALSimIMRPhenomBAmplitude_Merger(Mass,eta, chi, Normalization, ChiPowers);
	AmplitudeParameters *Ringdown	=	XLALSimIMRPhenomBAmplitude_Ringdown(Mass,eta, chi, fRing, ChiPowers);
	PhaseParameters		*Phase		=	XLALSimIMRPhenomBPhase(Mass,eta, chi);
	
	
	/* make sure that the flow is lower than the fCut */
	if (fCut < flow) {
		XLALPrintError("IMRPhenomB fCut is less than the flow chosen");
		XLAL_ERROR(XLAL_EDOM);
	}
	
	REAL8 df = Sh->deltaF;
	REAL8 hSqr = 0.;
	int s = 0;
	
	
	/* create a view of the PSD between flow and fCut */
	size_t nBins = (fCut - flow) / df;
	size_t k = nBins;
	REAL8Vector Shdata = {nBins, Sh->data->data + (size_t) (flow / df)}; /* copy the Vector, including its pointer to the actual data */
	/* drop low-frequency samples */
	Shdata.length = nBins;  /* drop high-frequency samples */
	
	/* allocate memory for various vectors */
	Amp = XLALCreateREAL8Vector(nBins);
	dATheta0 = XLALCreateREAL8Vector(nBins);
	dATheta3 = XLALCreateREAL8Vector(nBins);
	dATheta3S = XLALCreateREAL8Vector(nBins);
	dAT0 = XLALCreateREAL8Vector(nBins);
	dAPhi = XLALCreateREAL8Vector(nBins);
	dPhaseTheta0 = XLALCreateREAL8Vector(nBins);
	dPhaseTheta3 = XLALCreateREAL8Vector(nBins);
	dPhaseTheta3S = XLALCreateREAL8Vector(nBins);
	dPhaseT0 = XLALCreateREAL8Vector(nBins);
	dPhasePhi = XLALCreateREAL8Vector(nBins);
	
	
	/* derivative of the ampl w.r.t t0 and phi0 are zero. Fill these vectors with zero */
	memset(dAT0->data, 0, nBins*sizeof(REAL8));
	memset(dAPhi->data, 0, nBins*sizeof(REAL8));
	
	/* compute derivatives of the amplitude and phase of the waveform */
	for (;k--;) {
		
		const REAL8 f = flow + k * df;
		REAL8	freq_coef_1		=	cbrt(1./(f*f*f*f*f));
		REAL8	freq_coef_2		=	1./f;
		REAL8	freq_coef_3		=	cbrt(1./(f*f));
		REAL8	freq_coef_4		=	cbrt(1./(f));
		REAL8	freq_coef_5		=	cbrt(f);
		REAL8	freq_coef_6		=	cbrt(f*f);
		
		REAL8	phase_mass		=	Phase->phase_coef_1_mass*freq_coef_1 + Phase->phase_coef_2_mass*freq_coef_2 + Phase->phase_coef_3_mass*freq_coef_3 + Phase->phase_coef_4_mass*										freq_coef_4 + Phase->phase_coef_5_mass*freq_coef_5 + Phase->phase_coef_6_mass*freq_coef_6;
		REAL8	phase_eta		=	Phase->phase_coef_1_eta*freq_coef_1 + Phase->phase_coef_2_eta*freq_coef_2 + Phase->phase_coef_3_eta*freq_coef_3 + Phase->phase_coef_4_eta*											freq_coef_4 + Phase->phase_coef_5_eta*freq_coef_5 + Phase->phase_coef_6_eta*freq_coef_6;
		REAL8	phase_chi		=	Phase->phase_coef_1_chi*freq_coef_1 + Phase->phase_coef_2_chi*freq_coef_2 + Phase->phase_coef_3_chi*freq_coef_3 + Phase->phase_coef_4_chi*											freq_coef_4 + Phase->phase_coef_5_chi*freq_coef_5 + Phase->phase_coef_6_chi*freq_coef_6;
		
		dPhaseTheta0->data[k]	=	CalculateDerivatives(phase_mass, phase_eta, phase_chi, der_list_theta0);
		dPhaseTheta3->data[k]	=	CalculateDerivatives(phase_mass, phase_eta, phase_chi, der_list_theta3);
		dPhaseTheta3S->data[k]	=	CalculateDerivatives(phase_mass, phase_eta, phase_chi, der_list_theta3S);
		dPhaseT0->data[k]		=	LAL_TWOPI * f;
		dPhasePhi->data[k]		=	1.;
		
		if (f <= fMerg){
			REAL8	freq_pow_7_6	=	cbrt(sqrt(1./(f*f*f*f*f*f*f)));
			REAL8	freq_pow_3_6	=	sqrt(1./(f));
			REAL8	freq_pow_1_6	=	cbrt(sqrt(1./(f)));
			
			REAL8	amp_insp_mass	=	Inspiral->amp_coef_a_mass*freq_pow_7_6	+	Inspiral->amp_coef_b_mass*freq_pow_3_6	+	Inspiral->amp_coef_c_mass*freq_pow_1_6;
			REAL8	amp_insp_eta	=	Inspiral->amp_coef_a_eta*freq_pow_7_6	+	Inspiral->amp_coef_b_eta*freq_pow_3_6	+	Inspiral->amp_coef_c_eta*freq_pow_1_6;
			REAL8	amp_insp_chi	=	Inspiral->amp_coef_a_chi*freq_pow_7_6	+	Inspiral->amp_coef_b_chi*freq_pow_3_6	+	Inspiral->amp_coef_c_chi*freq_pow_1_6;
			
			/* inspiral amplitude of the waveform */
			Amp->data[k]		=	Inspiral->amp_coef_a*freq_pow_7_6	+	Inspiral->amp_coef_b*freq_pow_3_6	+	Inspiral->amp_coef_c*freq_pow_1_6;
			
			/* inspiral waveform deratives with respect to theta0 */
			dATheta0->data[k]	=	CalculateDerivatives(amp_insp_mass, amp_insp_eta, amp_insp_chi, der_list_theta0);
			
			/* inspiral waveform deratives with respect to theta3 */
			dATheta3->data[k]	=	CalculateDerivatives(amp_insp_mass, amp_insp_eta, amp_insp_chi, der_list_theta3);
			
			
			/* inspiral waveform deratives with respect to theta3S */
			dATheta3S->data[k]	=	CalculateDerivatives(amp_insp_mass, amp_insp_eta, amp_insp_chi, der_list_theta3S);
			
			
		}
		else if ((fMerg<f) && (f<=fRing)){
			REAL8	freq_pow_2_3  =		cbrt(1./(f*f));
			REAL8	freq_pow_1_3  =		cbrt(1./(f));
			
			/* merger amplitude of the frequency-domain waveform */
			REAL8 amp_merg			=	Merger->amp_coef_a*freq_pow_2_3 + Merger->amp_coef_b*freq_pow_1_3 + Merger->amp_coef_c;
			REAL8 amp_merg_mass		=	Merger->amp_coef_a_mass*freq_pow_2_3 + Merger->amp_coef_b_mass*freq_pow_1_3 + Merger->amp_coef_c_mass;
			REAL8 amp_merg_eta		=	Merger->amp_coef_a_eta*freq_pow_2_3  + Merger->amp_coef_b_eta*freq_pow_1_3 	+ Merger->amp_coef_c_eta;
			REAL8 amp_merg_chi		=	Merger->amp_coef_a_chi*freq_pow_2_3  + Merger->amp_coef_b_chi*freq_pow_1_3	+ Merger->amp_coef_c_chi;
			
			Amp->data[k]	=	amp_merg;
			
			/* merger waveform deratives with respect to theta0 */
			
			dATheta0->data[k]			=	CalculateDerivatives(amp_merg_mass, amp_merg_eta, amp_merg_chi, der_list_theta0);
			
			/* merger waveform deratives with respect to theta3 */
			
			dATheta3->data[k]			=	CalculateDerivatives(amp_merg_mass, amp_merg_eta, amp_merg_chi, der_list_theta3);
			
			/* merger waveform deratives with respect to theta3S */
			
			dATheta3S->data[k]			=	CalculateDerivatives(amp_merg_mass, amp_merg_eta, amp_merg_chi, der_list_theta3S);
			
			
			
		}
		
		else{
			REAL8	f_merg			=	Ringdown->amp_coef_a;
			REAL8	f_merg_mass		=	Ringdown->amp_coef_a_mass;
			REAL8	f_merg_eta		=	Ringdown->amp_coef_a_eta;
			REAL8	f_merg_chi		=	Ringdown->amp_coef_a_chi;
			REAL8	sigma			=	Ringdown->amp_coef_b;
			REAL8	sigma_mass		=	Ringdown->amp_coef_b_mass;
			REAL8	sigma_eta		=	Ringdown->amp_coef_b_eta;
			REAL8	sigma_chi		=	Ringdown->amp_coef_b_chi;
			REAL8	amp_const		=	Ringdown->amp_coef_c;
			REAL8	amp_const_mass	=	Ringdown->amp_coef_c_mass;
			REAL8	amp_const_eta	=	Ringdown->amp_coef_c_eta;
			REAL8	amp_const_chi	=	Ringdown->amp_coef_c_chi;
			
			REAL8	lorentzian		=	1./((f - f_merg)*(f - f_merg) + sigma*sigma*0.25);
			REAL8	amp_ring		=	norm_ring*amp_const*(1./(2*LAL_PI))*(sigma*lorentzian);
			REAL8	amp_ring_mass	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_mass*lorentzian + (2*sigma*(f_merg_mass*(f - f_merg) - sigma*sigma_mass*0.25))*												(lorentzian*lorentzian) ) +	 norm_ring_mass*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_mass*sigma);
			REAL8	amp_ring_eta	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_eta*lorentzian + (2*sigma*(f_merg_eta*(f - f_merg) - sigma*sigma_eta*0.25))*													(lorentzian*lorentzian) ) +	 norm_ring_eta*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_eta*sigma);
			REAL8	amp_ring_chi	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_chi*lorentzian + (2*sigma*(f_merg_chi*(f - f_merg) - sigma*sigma_chi*0.25))*													(lorentzian*lorentzian) ) +	 norm_ring_chi*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_chi*sigma);
			
			Amp->data[k]		=	amp_ring;
			dATheta0->data[k]	=	CalculateDerivatives(amp_ring_mass, amp_ring_eta, amp_ring_chi, der_list_theta0);
			dATheta3->data[k]	=	CalculateDerivatives(amp_ring_mass, amp_ring_eta, amp_ring_chi, der_list_theta3);
			dATheta3S->data[k]	=	CalculateDerivatives(amp_ring_mass, amp_ring_eta, amp_ring_chi, der_list_theta3S);
			
		}
		
		hSqr += Amp->data[k] * Amp->data[k] / Shdata.data[k];
		
	}
	hSqr = hSqr*df;
	/* allocate memory, and initialize the Fisher matrix */
	gsl_matrix * g = gsl_matrix_calloc (5, 5);
	
	/* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
	gsl_matrix_set (g, 0,0, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta0, dATheta0, dATheta0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,1, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta3, dATheta0, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,2, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta3S, dATheta0, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,3, MetricCoeffs(Amp, dPhaseTheta0, dPhaseT0, dATheta0, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,4, MetricCoeffs(Amp, dPhaseTheta0, dPhasePhi, dATheta0, dAPhi, &Shdata, hSqr, df));
	
	
	gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
	gsl_matrix_set (g, 1,1, MetricCoeffs(Amp, dPhaseTheta3, dPhaseTheta3, dATheta3, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,2, MetricCoeffs(Amp, dPhaseTheta3, dPhaseTheta3S, dATheta3, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,3, MetricCoeffs(Amp, dPhaseTheta3, dPhaseT0, dATheta3, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,4, MetricCoeffs(Amp, dPhaseTheta3, dPhasePhi, dATheta3, dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
	gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
	gsl_matrix_set (g, 2,2, MetricCoeffs(Amp, dPhaseTheta3S, dPhaseTheta3S, dATheta3S, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 2,3, MetricCoeffs(Amp, dPhaseTheta3S, dPhaseT0, dATheta3S, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 2,4, MetricCoeffs(Amp, dPhaseTheta3S, dPhasePhi, dATheta3S, dAPhi, &Shdata, hSqr, df));
	
	
	gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
	gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
	gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
	gsl_matrix_set (g, 3,3, MetricCoeffs(Amp, dPhaseT0, dPhaseT0, dAT0,dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 3,4, MetricCoeffs(Amp, dPhaseT0, dPhasePhi, dAT0,dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
	gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
	gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
	gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
	gsl_matrix_set (g, 4,4, MetricCoeffs(Amp, dPhasePhi, dPhasePhi, dAPhi,dAPhi, &Shdata, hSqr, df));
	
	
	
	/* free the memory */
	XLALDestroyREAL8Vector(Amp);
	XLALDestroyREAL8Vector(dATheta0);
	XLALDestroyREAL8Vector(dATheta3);
	XLALDestroyREAL8Vector(dAT0);
	XLALDestroyREAL8Vector(dAPhi);
	XLALDestroyREAL8Vector(dPhaseTheta0);
	XLALDestroyREAL8Vector(dPhaseTheta3);
	XLALDestroyREAL8Vector(dPhaseT0);
	XLALDestroyREAL8Vector(dPhasePhi);
	
	
	{
		/* Form submatrices g1, g2, g3, g4, defined as:
		 *              g = [ g1 g2
		 *                    g4 g3 ]                           */
		gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
		gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
		gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
		gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);
		
		
		
		gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
		gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
		gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
		gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
		gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);
		
		/* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
		gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
		gsl_permutation * p = gsl_permutation_calloc (2);
		
		gsl_matrix_memcpy (g1, &g1v.matrix);
		gsl_matrix_memcpy (g2, &g2v.matrix);
		gsl_matrix_memcpy (g3, &g3v.matrix);
		gsl_matrix_memcpy (g4, &g4v.matrix);
		gsl_matrix_free (g);
		
		gsl_linalg_LU_decomp (g3, p, &s);
		gsl_linalg_LU_invert (g3, p, g3inv);
		gsl_permutation_free (p);
		gsl_matrix_free (g3);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
		gsl_matrix_free (g4);
		gsl_matrix_free (g3inv);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
		gsl_matrix_free (g2);
		gsl_matrix_free (g3invg4);
		
		gsl_matrix_sub (g1, g2g3invg4);
		gsl_matrix_free (g2g3invg4);
		
		*gamma00 = gsl_matrix_get(g1, 0, 0);
		*gamma01 = gsl_matrix_get(g1, 0, 1);
		*gamma02 = gsl_matrix_get(g1, 0, 2);
		*gamma11 = gsl_matrix_get(g1, 1, 1);
		*gamma12 = gsl_matrix_get(g1, 1, 2);
		*gamma22 = gsl_matrix_get(g1, 2, 2);
		
		gsl_matrix_free (g1);
	}
	
	return XLAL_SUCCESS;
}

/**
 * Compute the two dimensional template-space metric of IMRPhenomB wavefrom in New Co-ordinates.
 */

int XLALSimIMRPhenomBMetricTheta0Theta3(
										REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
										REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
										REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
										const REAL8 Mass,     /**< Total Mass of the system */
										const REAL8 eta,    /**< Symmetric mass ratio */
										const REAL8 flow,   /**< low-frequency cutoff (Hz) */
										const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
){
	REAL8Vector *Amp=NULL, *dATheta0=NULL, *dATheta3=NULL;
	REAL8Vector *dAT0=NULL, *dAPhi=NULL, *dPhaseTheta0=NULL;
	REAL8Vector *dPhaseTheta3=NULL, *dPhaseT0=NULL, *dPhasePhi=NULL;
	REAL8 *der_list_theta0, *der_list_theta3, *der_list_theta3S;
	
	ConstantDefinitions *ChiPowers =	ChiPowList(0.);
	
	/* compute the chirp-time co-ordinates */
	const REAL8 theta0	=	ChirpTime_theta0(Mass,eta, flow);
	const REAL8 theta3	=	ChirpTime_theta3(Mass,eta, flow);
	const REAL8	theta3S	=	ChirpTime_theta3S(Mass, 0., flow);
	
	/* Compute the transition frequencies */
	
	const REAL8 fMerg  = TransitionFrequencies_fmerg(Mass,eta, 0., ChiPowers);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
	const REAL8 fRing  = TransitionFrequencies_fring(Mass,eta, 0., ChiPowers);	/**Frequency at which merger part transitions to ringdown part of the waveform*/
	const REAL8 fCut   = TransitionFrequencies_fcut(Mass,eta, 0.);	/**Frequency at which ringdown part of the waveform is terminated*/
	
	/* Compute the list of partial derivatives */
	der_list_theta0		=	MassPartialDerivativesWRTTheta0(theta0, theta3, theta3S, flow);
	der_list_theta3		=	MassPartialDerivativesWRTTheta3(theta0, theta3, theta3S, flow);
	der_list_theta3S	=	MassPartialDerivativesWRTTheta3S(theta0, theta3);
	
	/*Compute the normalizations and their derivatives*/
	NormalizationParameters *Normalization		=		XLALSimIMRPhenomBNormalization(Mass,eta, 0., ChiPowers);
	const REAL8 norm_ring			=		Normalization->norm_ring;
	const REAL8	norm_ring_mass		=		Normalization->norm_ring_mass;
	const REAL8	norm_ring_eta		=		Normalization->norm_ring_eta;
	
	AmplitudeParameters *Inspiral	=	XLALSimIMRPhenomBAmplitude_Inspiral(Mass,eta, 0.);
	AmplitudeParameters *Merger		=	XLALSimIMRPhenomBAmplitude_Merger(Mass,eta, 0., Normalization, ChiPowers);
	AmplitudeParameters *Ringdown	=	XLALSimIMRPhenomBAmplitude_Ringdown(Mass,eta, 0., fRing, ChiPowers);
	PhaseParameters		*Phase		=	XLALSimIMRPhenomBPhase(Mass,eta, 0.);
	
	/* make sure that the flow is lower than the fCut */
	if (fCut < flow) {
		XLALPrintError("IMRPhenomB fCut is less than the flow chosen");
		XLAL_ERROR(XLAL_EDOM);
	}
	
	REAL8 df = Sh->deltaF;
	REAL8 hSqr = 0.;
	int s = 0;
	
	
	/* create a view of the PSD between flow and fCut */
	size_t nBins = (fCut - flow) / df;
	size_t k = nBins;
	REAL8Vector Shdata = {nBins, Sh->data->data}; /* copy the Vector, including its pointer to the actual data */
	/* drop low-frequency samples */
	Shdata.length = nBins;  /* drop high-frequency samples */
	
	/* allocate memory for various vectors */
	Amp = XLALCreateREAL8Vector(nBins);
	dATheta0 = XLALCreateREAL8Vector(nBins);
	dATheta3 = XLALCreateREAL8Vector(nBins);
	dAT0 = XLALCreateREAL8Vector(nBins);
	dAPhi = XLALCreateREAL8Vector(nBins);
	dPhaseTheta0 = XLALCreateREAL8Vector(nBins);
	dPhaseTheta3 = XLALCreateREAL8Vector(nBins);
	dPhaseT0 = XLALCreateREAL8Vector(nBins);
	dPhasePhi = XLALCreateREAL8Vector(nBins);
	
	
	/* derivative of the ampl w.r.t t0 and phi0 are zero. Fill these vectors with zero */
	memset(dAT0->data, 0, nBins*sizeof(REAL8));
	memset(dAPhi->data, 0, nBins*sizeof(REAL8));
	
	/* compute derivatives of the amplitude and phase of the waveform */
	for (;k--;) {
		
		const REAL8 f = flow + k * df;
		REAL8	freq_coef_1		=	cbrt(1./(f*f*f*f*f));
		REAL8	freq_coef_2		=	1./f;
		REAL8	freq_coef_3		=	cbrt(1./(f*f));
		REAL8	freq_coef_4		=	cbrt(1./(f));
		REAL8	freq_coef_5		=	cbrt(f);
		REAL8	freq_coef_6		=	cbrt(f*f);
		
		REAL8	phase_mass		=	Phase->phase_coef_1_mass*freq_coef_1 + Phase->phase_coef_2_mass*freq_coef_2 + Phase->phase_coef_3_mass*freq_coef_3 + Phase->phase_coef_4_mass*										freq_coef_4 + Phase->phase_coef_5_mass*freq_coef_5 + Phase->phase_coef_6_mass*freq_coef_6;
		REAL8	phase_eta		=	Phase->phase_coef_1_eta*freq_coef_1 + Phase->phase_coef_2_eta*freq_coef_2 + Phase->phase_coef_3_eta*freq_coef_3 + Phase->phase_coef_4_eta*											freq_coef_4 + Phase->phase_coef_5_eta*freq_coef_5 + Phase->phase_coef_6_eta*freq_coef_6;
		
		dPhaseTheta0->data[k]	=	CalculateDerivatives(phase_mass, phase_eta, 0., der_list_theta0);
		dPhaseTheta3->data[k]	=	CalculateDerivatives(phase_mass, phase_eta, 0., der_list_theta3);
		dPhaseT0->data[k]		=	LAL_TWOPI * f;
		dPhasePhi->data[k]		=	1.;
		
		if (f <= fMerg){
			REAL8	freq_pow_7_6	=	cbrt(sqrt(1./(f*f*f*f*f*f*f)));
			REAL8	freq_pow_3_6	=	sqrt(1./(f));
			REAL8	freq_pow_1_6	=	cbrt(sqrt(1./(f)));
			
			REAL8	amp_insp_mass	=	Inspiral->amp_coef_a_mass*freq_pow_7_6	+	Inspiral->amp_coef_b_mass*freq_pow_3_6	+	Inspiral->amp_coef_c_mass*freq_pow_1_6;
			REAL8	amp_insp_eta	=	Inspiral->amp_coef_a_eta*freq_pow_7_6	+	Inspiral->amp_coef_b_eta*freq_pow_3_6	+	Inspiral->amp_coef_c_eta*freq_pow_1_6;
			
			/* inspiral amplitude of the waveform */
			Amp->data[k]		=	Inspiral->amp_coef_a*freq_pow_7_6	+	Inspiral->amp_coef_b*freq_pow_3_6	+	Inspiral->amp_coef_c*freq_pow_1_6;
			
			/* inspiral waveform deratives with respect to theta0 */
			dATheta0->data[k]	=	CalculateDerivatives(amp_insp_mass, amp_insp_eta, 0., der_list_theta0);
			
			/* inspiral waveform deratives with respect to theta3 */
			dATheta3->data[k]	=	CalculateDerivatives(amp_insp_mass, amp_insp_eta, 0., der_list_theta3);
		}
		else if ((fMerg<f) && (f<=fRing)){
			REAL8	freq_pow_2_3  =		cbrt(1./(f*f));
			REAL8	freq_pow_1_3  =		cbrt(1./(f));
			
			/* merger amplitude of the frequency-domain waveform */
			REAL8 amp_merg			=	Merger->amp_coef_a*freq_pow_2_3 + Merger->amp_coef_b*freq_pow_1_3 + Merger->amp_coef_c;
			REAL8 amp_merg_mass		=	Merger->amp_coef_a_mass*freq_pow_2_3 + Merger->amp_coef_b_mass*freq_pow_1_3 + Merger->amp_coef_c_mass;
			REAL8 amp_merg_eta		=	Merger->amp_coef_a_eta*freq_pow_2_3  + Merger->amp_coef_b_eta*freq_pow_1_3 	+ Merger->amp_coef_c_eta;
			
			Amp->data[k]	=	amp_merg;
			dATheta0->data[k]			=	CalculateDerivatives(amp_merg_mass, amp_merg_eta, 0., der_list_theta0);
			dATheta3->data[k]			=	CalculateDerivatives(amp_merg_mass, amp_merg_eta, 0., der_list_theta3);
		}
		
		else{
			REAL8	f_merg			=	Ringdown->amp_coef_a;
			REAL8	f_merg_mass		=	Ringdown->amp_coef_a_mass;
			REAL8	f_merg_eta		=	Ringdown->amp_coef_a_eta;
			REAL8	sigma			=	Ringdown->amp_coef_b;
			REAL8	sigma_mass		=	Ringdown->amp_coef_b_mass;
			REAL8	sigma_eta		=	Ringdown->amp_coef_b_eta;
			REAL8	amp_const		=	Ringdown->amp_coef_c;
			REAL8	amp_const_mass	=	Ringdown->amp_coef_c_mass;
			REAL8	amp_const_eta	=	Ringdown->amp_coef_c_eta;
			
			REAL8	lorentzian		=	1./((f - f_merg)*(f - f_merg) + sigma*sigma*0.25);
			REAL8	amp_ring		=	norm_ring*amp_const*(1./(2*LAL_PI))*(sigma*lorentzian);
			REAL8	amp_ring_mass	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_mass*lorentzian + (2*sigma*(f_merg_mass*(f - f_merg) - sigma*sigma_mass*0.25))*												(lorentzian*lorentzian) ) +	 norm_ring_mass*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_mass*sigma);
			REAL8	amp_ring_eta	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_eta*lorentzian + (2*sigma*(f_merg_eta*(f - f_merg) - sigma*sigma_eta*0.25))*													(lorentzian*lorentzian) ) +	 norm_ring_eta*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_eta*sigma);
			
			Amp->data[k]		=	amp_ring;
			dATheta0->data[k]	=	CalculateDerivatives(amp_ring_mass, amp_ring_eta, 0., der_list_theta0);
			dATheta3->data[k]	=	CalculateDerivatives(amp_ring_mass, amp_ring_eta, 0., der_list_theta3);
			
		}
		
		hSqr += Amp->data[k] * Amp->data[k] / Shdata.data[k];
	}
	hSqr = hSqr*df;
	
	/* allocate memory, and initialize the Fisher matrix */
	gsl_matrix * g = gsl_matrix_calloc (4, 4);
	
	/* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
	gsl_matrix_set (g, 0,0, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta0, dATheta0, dATheta0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,1, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta3, dATheta0, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,2, MetricCoeffs(Amp, dPhaseTheta0, dPhaseT0, dATheta0, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,3, MetricCoeffs(Amp, dPhaseTheta0, dPhasePhi, dATheta0, dAPhi, &Shdata, hSqr, df));
	
	
	gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
	gsl_matrix_set (g, 1,1, MetricCoeffs(Amp, dPhaseTheta3, dPhaseTheta3, dATheta3, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,2, MetricCoeffs(Amp, dPhaseTheta3, dPhaseT0, dATheta3, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,3, MetricCoeffs(Amp, dPhaseTheta3, dPhasePhi, dATheta3, dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
	gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
	gsl_matrix_set (g, 2,2, MetricCoeffs(Amp, dPhaseT0, dPhaseT0, dAT0, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 2,3, MetricCoeffs(Amp, dPhaseT0, dPhasePhi, dAT0, dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
	gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
	gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
	gsl_matrix_set (g, 3,3, MetricCoeffs(Amp, dPhasePhi, dPhasePhi, dAPhi,dAPhi, &Shdata, hSqr, df));
	
	/* free the memory */
	XLALDestroyREAL8Vector(Amp);
	XLALDestroyREAL8Vector(dATheta0);
	XLALDestroyREAL8Vector(dATheta3);
	XLALDestroyREAL8Vector(dAT0);
	XLALDestroyREAL8Vector(dAPhi);
	XLALDestroyREAL8Vector(dPhaseTheta0);
	XLALDestroyREAL8Vector(dPhaseTheta3);
	XLALDestroyREAL8Vector(dPhaseT0);
	XLALDestroyREAL8Vector(dPhasePhi);
	
	{
		/* Form submatrices g1, g2, g3, g4, defined as:
		 *              g = [ g1 g2
		 *                    g4 g3 ]                           */
		gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 2, 2);
		gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 2, 2, 2);
		gsl_matrix_view g3v = gsl_matrix_submatrix (g, 2, 2, 2, 2);
		gsl_matrix_view g4v = gsl_matrix_submatrix (g, 2, 0, 2, 2);
		
		gsl_matrix * g1 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g2 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g4 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g2g3invg4 = gsl_matrix_calloc (2, 2);
		
		/* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
		gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
		gsl_permutation * p = gsl_permutation_calloc (2);
		
		gsl_matrix_memcpy (g1, &g1v.matrix);
		gsl_matrix_memcpy (g2, &g2v.matrix);
		gsl_matrix_memcpy (g3, &g3v.matrix);
		gsl_matrix_memcpy (g4, &g4v.matrix);
		gsl_matrix_free (g);
		
		gsl_linalg_LU_decomp (g3, p, &s);
		gsl_linalg_LU_invert (g3, p, g3inv);
		gsl_permutation_free (p);
		gsl_matrix_free (g3);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
		gsl_matrix_free (g4);
		gsl_matrix_free (g3inv);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
		gsl_matrix_free (g2);
		gsl_matrix_free (g3invg4);
		
		gsl_matrix_sub (g1, g2g3invg4);
		gsl_matrix_free (g2g3invg4);
		
		*gamma00 = gsl_matrix_get(g1, 0, 0);
		*gamma01 = gsl_matrix_get(g1, 0, 1);
		*gamma11 = gsl_matrix_get(g1, 1, 1);
		gsl_matrix_free (g1);
	}
	
	return XLAL_SUCCESS;
}


/**
 * Compute the two dimensional template-space metric of IMRPhenomB wavefrom in Mass Co-ordinates.
 */

int XLALSimIMRPhenomBMetricMassEtaChi(
									  REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
									  REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
									  REAL8 *gamma02,  /**< template metric coeff. 01/10 PN Chirp Time */
									  REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
									  REAL8 *gamma12,  /**< template metric coeff. 01/10 PN Chirp Time */
									  REAL8 *gamma22,  /**< template metric coeff. 01/10 PN Chirp Time */
									  const REAL8 Mass,     /**< Total Mass of the system */
									  const REAL8 eta,    /**< Symmetric mass ratio */
									  const REAL8 chi,	   /** Reduced Spin Parameter of the system **/
									  const REAL8 flow,   /**< low-frequency cutoff (Hz) */
									  const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
){
	REAL8Vector *Amp=NULL, *dATheta0=NULL, *dATheta3=NULL, *dATheta3S=NULL;
	REAL8Vector *dAT0=NULL, *dAPhi=NULL, *dPhaseTheta0=NULL;
	REAL8Vector *dPhaseTheta3=NULL, *dPhaseTheta3S=NULL, *dPhaseT0=NULL, *dPhasePhi=NULL;
	
	ConstantDefinitions *ChiPowers =	ChiPowList(chi);
	
	
	/* Compute the transition frequencies */
	
	const REAL8 fMerg  = TransitionFrequencies_fmerg(Mass,eta, chi, ChiPowers);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
	const REAL8 fRing  = TransitionFrequencies_fring(Mass,eta, chi, ChiPowers);	/**Frequency at which merger part transitions to ringdown part of the waveform*/
	const REAL8 fCut   = TransitionFrequencies_fcut(Mass,eta, chi);	/**Frequency at which ringdown part of the waveform is terminated*/
	
	
	/*Compute the normalizations and their derivatives*/
	NormalizationParameters *Normalization		=		XLALSimIMRPhenomBNormalization(Mass,eta, chi, ChiPowers);
	const REAL8 norm_ring			=		Normalization->norm_ring;
	const REAL8	norm_ring_mass		=		Normalization->norm_ring_mass;
	const REAL8	norm_ring_eta		=		Normalization->norm_ring_eta;
	const REAL8	norm_ring_chi		=		Normalization->norm_ring_chi;
	
	AmplitudeParameters *Inspiral	=	XLALSimIMRPhenomBAmplitude_Inspiral(Mass,eta, chi);
	AmplitudeParameters *Merger		=	XLALSimIMRPhenomBAmplitude_Merger(Mass,eta, chi, Normalization, ChiPowers);
	AmplitudeParameters *Ringdown	=	XLALSimIMRPhenomBAmplitude_Ringdown(Mass,eta, chi, fRing, ChiPowers);
	PhaseParameters		*Phase		=	XLALSimIMRPhenomBPhase(Mass,eta, chi);
	
	
	/* make sure that the flow is lower than the fCut */
	if (fCut < flow) {
		XLALPrintError("IMRPhenomB fCut is less than the flow chosen");
		XLAL_ERROR(XLAL_EDOM);
	}
	
	REAL8 df = Sh->deltaF;
	REAL8 hSqr = 0.;
	int s = 0;
	
	
	/* create a view of the PSD between flow and fCut */
	size_t nBins = (fCut - flow) / df;
	size_t k = nBins;
	REAL8Vector Shdata = {nBins, Sh->data->data + (size_t) (flow / df)}; /* copy the Vector, including its pointer to the actual data */
	/* drop low-frequency samples */
	Shdata.length = nBins;  /* drop high-frequency samples */
	
	/* allocate memory for various vectors */
	Amp = XLALCreateREAL8Vector(nBins);
	dATheta0 = XLALCreateREAL8Vector(nBins);
	dATheta3 = XLALCreateREAL8Vector(nBins);
	dATheta3S = XLALCreateREAL8Vector(nBins);
	dAT0 = XLALCreateREAL8Vector(nBins);
	dAPhi = XLALCreateREAL8Vector(nBins);
	dPhaseTheta0 = XLALCreateREAL8Vector(nBins);
	dPhaseTheta3 = XLALCreateREAL8Vector(nBins);
	dPhaseTheta3S = XLALCreateREAL8Vector(nBins);
	dPhaseT0 = XLALCreateREAL8Vector(nBins);
	dPhasePhi = XLALCreateREAL8Vector(nBins);
	
	
	/* derivative of the ampl w.r.t t0 and phi0 are zero. Fill these vectors with zero */
	memset(dAT0->data, 0, nBins*sizeof(REAL8));
	memset(dAPhi->data, 0, nBins*sizeof(REAL8));
	
	/* compute derivatives of the amplitude and phase of the waveform */
	for (;k--;) {
		
		const REAL8 f = flow + k * df;
		REAL8	freq_coef_1		=	cbrt(1./(f*f*f*f*f));
		REAL8	freq_coef_2		=	1./f;
		REAL8	freq_coef_3		=	cbrt(1./(f*f));
		REAL8	freq_coef_4		=	cbrt(1./(f));
		REAL8	freq_coef_5		=	cbrt(f);
		REAL8	freq_coef_6		=	cbrt(f*f);
		
		REAL8	phase_mass		=	Phase->phase_coef_1_mass*freq_coef_1 + Phase->phase_coef_2_mass*freq_coef_2 + Phase->phase_coef_3_mass*freq_coef_3 + Phase->phase_coef_4_mass*										freq_coef_4 + Phase->phase_coef_5_mass*freq_coef_5 + Phase->phase_coef_6_mass*freq_coef_6;
		REAL8	phase_eta		=	Phase->phase_coef_1_eta*freq_coef_1 + Phase->phase_coef_2_eta*freq_coef_2 + Phase->phase_coef_3_eta*freq_coef_3 + Phase->phase_coef_4_eta*											freq_coef_4 + Phase->phase_coef_5_eta*freq_coef_5 + Phase->phase_coef_6_eta*freq_coef_6;
		REAL8	phase_chi		=	Phase->phase_coef_1_chi*freq_coef_1 + Phase->phase_coef_2_chi*freq_coef_2 + Phase->phase_coef_3_chi*freq_coef_3 + Phase->phase_coef_4_chi*											freq_coef_4 + Phase->phase_coef_5_chi*freq_coef_5 + Phase->phase_coef_6_chi*freq_coef_6;
		
		dPhaseTheta0->data[k]	=	phase_mass;
		dPhaseTheta3->data[k]	=	phase_eta;
		dPhaseTheta3S->data[k]	=	phase_chi;
		dPhaseT0->data[k]		=	LAL_TWOPI * f;
		dPhasePhi->data[k]		=	1.;
		
		if (f <= fMerg){
			REAL8	freq_pow_7_6	=	cbrt(sqrt(1./(f*f*f*f*f*f*f)));
			REAL8	freq_pow_3_6	=	sqrt(1./(f));
			REAL8	freq_pow_1_6	=	cbrt(sqrt(1./(f)));
			
			REAL8	amp_insp_mass	=	Inspiral->amp_coef_a_mass*freq_pow_7_6	+	Inspiral->amp_coef_b_mass*freq_pow_3_6	+	Inspiral->amp_coef_c_mass*freq_pow_1_6;
			REAL8	amp_insp_eta	=	Inspiral->amp_coef_a_eta*freq_pow_7_6	+	Inspiral->amp_coef_b_eta*freq_pow_3_6	+	Inspiral->amp_coef_c_eta*freq_pow_1_6;
			REAL8	amp_insp_chi	=	Inspiral->amp_coef_a_chi*freq_pow_7_6	+	Inspiral->amp_coef_b_chi*freq_pow_3_6	+	Inspiral->amp_coef_c_chi*freq_pow_1_6;
			
			/* inspiral amplitude of the waveform */
			Amp->data[k]		=	Inspiral->amp_coef_a*freq_pow_7_6	+	Inspiral->amp_coef_b*freq_pow_3_6	+	Inspiral->amp_coef_c*freq_pow_1_6;
			
			/* inspiral waveform deratives with respect to theta0 */
			dATheta0->data[k]	=	amp_insp_mass;
			
			/* inspiral waveform deratives with respect to theta3 */
			dATheta3->data[k]	=	amp_insp_eta;
			
			
			/* inspiral waveform deratives with respect to theta3S */
			dATheta3S->data[k]	=	amp_insp_chi;
			
			
		}
		else if ((fMerg<f) && (f<=fRing)){
			REAL8	freq_pow_2_3  =		cbrt(1./(f*f));
			REAL8	freq_pow_1_3  =		cbrt(1./(f));
			
			/* merger amplitude of the frequency-domain waveform */
			REAL8 amp_merg			=	Merger->amp_coef_a*freq_pow_2_3 + Merger->amp_coef_b*freq_pow_1_3 + Merger->amp_coef_c;
			REAL8 amp_merg_mass		=	Merger->amp_coef_a_mass*freq_pow_2_3 + Merger->amp_coef_b_mass*freq_pow_1_3 + Merger->amp_coef_c_mass;
			REAL8 amp_merg_eta		=	Merger->amp_coef_a_eta*freq_pow_2_3  + Merger->amp_coef_b_eta*freq_pow_1_3 	+ Merger->amp_coef_c_eta;
			REAL8 amp_merg_chi		=	Merger->amp_coef_a_chi*freq_pow_2_3  + Merger->amp_coef_b_chi*freq_pow_1_3	+ Merger->amp_coef_c_chi;
			
			Amp->data[k]	=	amp_merg;
			
			/* merger waveform deratives with respect to theta0 */
			
			dATheta0->data[k]			=	amp_merg_mass;
			
			/* merger waveform deratives with respect to theta3 */
			
			dATheta3->data[k]			=	amp_merg_eta;
			
			/* merger waveform deratives with respect to theta3S */
			
			dATheta3S->data[k]			=	amp_merg_chi;
			
			
			
		}
		
		else{
			REAL8	f_merg			=	Ringdown->amp_coef_a;
			REAL8	f_merg_mass		=	Ringdown->amp_coef_a_mass;
			REAL8	f_merg_eta		=	Ringdown->amp_coef_a_eta;
			REAL8	f_merg_chi		=	Ringdown->amp_coef_a_chi;
			REAL8	sigma			=	Ringdown->amp_coef_b;
			REAL8	sigma_mass		=	Ringdown->amp_coef_b_mass;
			REAL8	sigma_eta		=	Ringdown->amp_coef_b_eta;
			REAL8	sigma_chi		=	Ringdown->amp_coef_b_chi;
			REAL8	amp_const		=	Ringdown->amp_coef_c;
			REAL8	amp_const_mass	=	Ringdown->amp_coef_c_mass;
			REAL8	amp_const_eta	=	Ringdown->amp_coef_c_eta;
			REAL8	amp_const_chi	=	Ringdown->amp_coef_c_chi;
			
			REAL8	lorentzian		=	1./((f - f_merg)*(f - f_merg) + sigma*sigma*0.25);
			REAL8	amp_ring		=	norm_ring*amp_const*(1./(2*LAL_PI))*(sigma*lorentzian);
			REAL8	amp_ring_mass	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_mass*lorentzian + (2*sigma*(f_merg_mass*(f - f_merg) - sigma*sigma_mass*0.25))*												(lorentzian*lorentzian) ) +	 norm_ring_mass*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_mass*sigma);
			REAL8	amp_ring_eta	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_eta*lorentzian + (2*sigma*(f_merg_eta*(f - f_merg) - sigma*sigma_eta*0.25))*													(lorentzian*lorentzian) ) +	 norm_ring_eta*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_eta*sigma);
			REAL8	amp_ring_chi	=	(1./(2*LAL_PI))*(norm_ring*amp_const*( sigma_chi*lorentzian + (2*sigma*(f_merg_chi*(f - f_merg) - sigma*sigma_chi*0.25))*													(lorentzian*lorentzian) ) +	 norm_ring_chi*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_chi*sigma);
			
			Amp->data[k]		=	amp_ring;
			dATheta0->data[k]	=	amp_ring_mass;
			dATheta3->data[k]	=	amp_ring_eta;
			dATheta3S->data[k]	=	amp_ring_chi;
			
		}
		
		hSqr += Amp->data[k] * Amp->data[k] / Shdata.data[k];
		
	}
	hSqr = hSqr*df;
	/* allocate memory, and initialize the Fisher matrix */
	gsl_matrix * g = gsl_matrix_calloc (5, 5);
	
	/* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
	gsl_matrix_set (g, 0,0, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta0, dATheta0, dATheta0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,1, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta3, dATheta0, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,2, MetricCoeffs(Amp, dPhaseTheta0, dPhaseTheta3S, dATheta0, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,3, MetricCoeffs(Amp, dPhaseTheta0, dPhaseT0, dATheta0, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 0,4, MetricCoeffs(Amp, dPhaseTheta0, dPhasePhi, dATheta0, dAPhi, &Shdata, hSqr, df));
	
	
	gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
	gsl_matrix_set (g, 1,1, MetricCoeffs(Amp, dPhaseTheta3, dPhaseTheta3, dATheta3, dATheta3, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,2, MetricCoeffs(Amp, dPhaseTheta3, dPhaseTheta3S, dATheta3, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,3, MetricCoeffs(Amp, dPhaseTheta3, dPhaseT0, dATheta3, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 1,4, MetricCoeffs(Amp, dPhaseTheta3, dPhasePhi, dATheta3, dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
	gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
	gsl_matrix_set (g, 2,2, MetricCoeffs(Amp, dPhaseTheta3S, dPhaseTheta3S, dATheta3S, dATheta3S, &Shdata, hSqr, df));
	gsl_matrix_set (g, 2,3, MetricCoeffs(Amp, dPhaseTheta3S, dPhaseT0, dATheta3S, dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 2,4, MetricCoeffs(Amp, dPhaseTheta3S, dPhasePhi, dATheta3S, dAPhi, &Shdata, hSqr, df));
	
	
	gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
	gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
	gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
	gsl_matrix_set (g, 3,3, MetricCoeffs(Amp, dPhaseT0, dPhaseT0, dAT0,dAT0, &Shdata, hSqr, df));
	gsl_matrix_set (g, 3,4, MetricCoeffs(Amp, dPhaseT0, dPhasePhi, dAT0,dAPhi, &Shdata, hSqr, df));
	
	gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
	gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
	gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
	gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
	gsl_matrix_set (g, 4,4, MetricCoeffs(Amp, dPhasePhi, dPhasePhi, dAPhi,dAPhi, &Shdata, hSqr, df));
	
	
	
	/* free the memory */
	XLALDestroyREAL8Vector(Amp);
	XLALDestroyREAL8Vector(dATheta0);
	XLALDestroyREAL8Vector(dATheta3);
	XLALDestroyREAL8Vector(dAT0);
	XLALDestroyREAL8Vector(dAPhi);
	XLALDestroyREAL8Vector(dPhaseTheta0);
	XLALDestroyREAL8Vector(dPhaseTheta3);
	XLALDestroyREAL8Vector(dPhaseT0);
	XLALDestroyREAL8Vector(dPhasePhi);
	
	
	{
		/* Form submatrices g1, g2, g3, g4, defined as:
		 *              g = [ g1 g2
		 *                    g4 g3 ]                           */
		gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
		gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
		gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
		gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);
		
		
		
		gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
		gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
		gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
		gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
		gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
		gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);
		
		/* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
		gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
		gsl_permutation * p = gsl_permutation_calloc (2);
		
		gsl_matrix_memcpy (g1, &g1v.matrix);
		gsl_matrix_memcpy (g2, &g2v.matrix);
		gsl_matrix_memcpy (g3, &g3v.matrix);
		gsl_matrix_memcpy (g4, &g4v.matrix);
		gsl_matrix_free (g);
		
		gsl_linalg_LU_decomp (g3, p, &s);
		gsl_linalg_LU_invert (g3, p, g3inv);
		gsl_permutation_free (p);
		gsl_matrix_free (g3);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
		gsl_matrix_free (g4);
		gsl_matrix_free (g3inv);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
		gsl_matrix_free (g2);
		gsl_matrix_free (g3invg4);
		
		gsl_matrix_sub (g1, g2g3invg4);
		gsl_matrix_free (g2g3invg4);
		
		*gamma00 = gsl_matrix_get(g1, 0, 0);
		*gamma01 = gsl_matrix_get(g1, 0, 1);
		*gamma02 = gsl_matrix_get(g1, 0, 2);
		*gamma11 = gsl_matrix_get(g1, 1, 1);
		*gamma12 = gsl_matrix_get(g1, 1, 2);
		*gamma22 = gsl_matrix_get(g1, 2, 2);
		
		gsl_matrix_free (g1);
	}
	
	return XLAL_SUCCESS;
}

