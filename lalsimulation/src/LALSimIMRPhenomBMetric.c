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


#define PI 3.1415926535897932384626433832795029

/**
*******************************************************************************************************************************
*/
/**
* Define the Chirp time co-ordinates
*/


static REAL8 ChirpTime_theta0(
	const REAL8 mass,	/**< Total Mass of system */
	const REAL8 eta,	/**< Symmetric mass ratio of system */
	const REAL8 f_low	/**< Lower Frequency Cut-off */
	){
	REAL8 pi;
	pi = PI;
	return 5.0/(128.0*eta*pow(pi*mass*LAL_MTSUN_SI*f_low,5.0/3.0));
	}

static REAL8 ChirpTime_theta3(
	const REAL8 mass,	/**< Total Mass of system */
	const REAL8 eta,	/**< Symmetric mass ratio of system */
	const REAL8 f_low	/**< Lower Frequency Cut-off */
	){
	REAL8 pi;
	pi = PI;
	return pow(pi,1.0/3.0)/(4.0*eta*pow(mass*LAL_MTSUN_SI*f_low,2.0/3.0));
	}
/**
*******************************************************************************************************************************
*/
/**
* Function to calculate the transition frequencies of IMRPhenomB waveform. F_merg, F_ring and F_cut.

	F_merg = Frequency at which waveform transitions from inspiral to merger mode.
	F_ring = Frequency at which waveform transitions from merger to ringdown mode.
	F_cut = Frequency at which IMRPhenomB waveform is terminated.
*/

static REAL8 TransitionFrequencies_fmerg(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 f_low	/**< Lower Frequency Cut-off */
	)
	{
	return (1.3270087368763253*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. - (1388.908285838913*pow(theta0,2))/pow(theta3,5) - (1.9634062693265484*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) + (3.7378874416543857*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3;


}

static REAL8 TransitionFrequencies_fring(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 f_low	/**< Lower Frequency Cut-off */
	)
	{
	return (3.719645701850315*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. + (455.39646148015123*pow(theta0,2))/pow(theta3,5) - (0.8397543046176621*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) + (0.8530966599534362*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3;

}

static REAL8 TransitionFrequencies_fcut(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 f_low	/**< Lower Frequency Cut-off */
	)
	{
	return (6.506364049290605*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. + (963.986488648419*pow(theta0,2))/pow(theta3,5) - (9.15298466960777*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) - (0.7731297368250576*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3;
}

/**
*******************************************************************************************************************************
*/
/**
*Define the Normalization constants and it's derivatives required for Merger and Ringdown parts of the IMRPhenomB amplitude and the amplitude derivatives
*/

static REAL8 XLALSimIMRPhenomBNormalization_Merger(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);

	return ((4.465858060541566*theta3_pow_166 + 9.415904762816128*
      pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
        (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
        (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
      (1.*theta0_pow_one_third - 0.09249341147745196*theta3_pow_166))/
   ((4.465858060541567 - 3.1035196288177636*cbrt(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*theta3_pow_166)*flow/flow);

	}


static REAL8 XLALSimIMRPhenomBNormalization_Ringdown(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = pow(theta3,6); 

	/*Recurring Expression*/
	REAL8	recurring_expr_1 = pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333);


	return ((-37427.24274521482*flow*theta0*(1.*theta0_pow_2 - 0.10972122263753406*theta0_pow_133*theta3_pow_166 + 
       0.004234161619441117*theta0_pow_one_third*theta3_pow_33 - 0.00016457382536589995*theta3_pow_5)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.6666666666666666)*(0.4742887882827215*theta3_pow_166 + 
       1.*recurring_expr_1*recurring_expr_1*
        (1.*theta0_pow_one_third - 0.09249341147745196*theta3_pow_166))*
     (4.465858060541567 - 3.1035196288177636*pow(theta3/(flow*theta0),0.3333333333333333)*
        pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
              17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
         0.3333333333333333) + 1.*pow(theta3/(flow*theta0),0.6666666666666666)*
        pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
              17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
         0.6666666666666666)))/((4.465858060541567 - 3.1035196288177636* recurring_expr_1 + 
       1.*recurring_expr_1*recurring_expr_1)*theta3_pow_6*theta3*cbrt(theta3)*cbrt(theta3)*
     pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
           17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
      0.6666666666666666)));
	}


static REAL8 XLALSimIMRPhenomBNormalization_Merger_Der_theta0(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_333 = pow(theta0,0.3333333333333333);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);

	/*Some Recurring Expressions*/
	REAL8 recurring_expr_1 = pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
           (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
           (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333);

	return ((4.895281031102956e9*theta0_pow_2*theta0_pow_2 + 1.384025938152937e7*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_166 - 
     2.6338983161110472e7*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33 + 
     theta0_pow_2*(-502488.4903868959 - 2.348553100590565e6*recurring_expr_1 + 
        1.088074683997717e6*recurring_expr_1*recurring_expr_1)*theta3_pow_5 + 
     theta0_pow_133*(-223233.04749560714 + 308517.60984538106*
         recurring_expr_1 - 
        49038.08804007761*recurring_expr_1*recurring_expr_1)*theta3_pow_33*theta3_pow_33 + 
     theta0_pow_one_third*(1008.9019847050173 + 3804.589089140962*
         recurring_expr_1 - 
        1999.603615336666*recurring_expr_1*recurring_expr_1)*pow(theta3,8.333333333333334) + 
     (242.52837069765866 - 241.8810122753211*recurring_expr_1 + 
        19.288678319713693*recurring_expr_1*recurring_expr_1)*theta3_pow_5*theta3_pow_5)/
   (theta3_pow_333*pow(4.465858060541567 - 3.1035196288177636*
        pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666),2)*
     pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
       (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
       (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*pow(theta3,11.666666666666666))*flow/flow);

	}


static REAL8 XLALSimIMRPhenomBNormalization_Merger_Der_theta3(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);

	/*Some Recurring Expressions*/
	REAL8	recurring_expr_1 = pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
             (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
             (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333);


	return ((theta0_pow_one_third*(theta0_pow_2*theta0_pow_2*(5.7243785765757225e10 - 
          9.222397892093138e9*recurring_expr_1) + 
       theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*(-2.49898008123127e9 - 2.60741677822487e7*
           recurring_expr_1)*theta3_pow_166 + 
       theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*(-2.6290334909565797e8 + 4.962096787529522e7*
           recurring_expr_1)*theta3_pow_33 + 
       theta0_pow_2*(5.089018756458765e6 + 946656.3339471949*recurring_expr_1 + 
          4.424524562888101e6*recurring_expr_1*recurring_expr_1)*theta3_pow_5 + 
       theta0_pow_133*(407146.9153272246 + 420556.8533427301*
           recurring_expr_1 - 
          581227.5406935291*recurring_expr_1*recurring_expr_1)*theta3_pow_33*theta3_pow_33 + 
       theta0_pow_one_third*(2267.989778592834 - 1900.7071254856517*
           recurring_expr_1 - 
          7167.603692829992*recurring_expr_1*recurring_expr_1)*pow(theta3,8.333333333333334) + 
       (-48.22169579928426 - 456.90801416377934*recurring_expr_1 + 
          455.68843209859017*recurring_expr_1*recurring_expr_1)*theta3_pow_5*theta3_pow_5))/
   (pow(4.465858060541567 - 3.1035196288177636*recurring_expr_1 + 
       1.*recurring_expr_1*recurring_expr_1,2)*theta3_pow_5*theta3*theta3*cbrt(theta3)*cbrt(theta3)*
     (-21044.06493695328*theta0_pow_2 - 29.748579838281113*theta0_pow_133*theta3_pow_166 + 
       56.63465820688478*theta0_pow_one_third*theta3_pow_33 + 1.*theta3_pow_5))*flow/flow);

	}

static REAL8 XLALSimIMRPhenomBNormalization_Ringdown_Der_theta0(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_merg,	/**< Normalization for merger waveform */
	const REAL8 norm_ring, /**< Normalization for ringdown waveform */
	const REAL8 norm_merg_theta0 /**< Derviative of norm_merg w.r.t theta0*/
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = pow(theta3,6);
	REAL8	theta0theta3_comb = pow(theta3/(flow*theta0),0.3333333333333333); 

	REAL8	norm_ring_theta0;

	/*Some Recurring Expressions*/

	REAL8 recurring_expr_1 = pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926418*theta0_pow_133*theta3_pow_166 + 
                 17.152526078154914*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
            0.3333333333333333);

	norm_ring_theta0 = ((4.81006881704554e8*flow*flow*theta0*(theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_166*
        (-0.26666856187261073 + 0.23761962032750067*theta0theta3_comb*
           recurring_expr_1 - 0.09341639683655707*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33*theta3_pow_33*
        (-0.0001659902885655082 + 0.00011805488681116238*theta0theta3_comb*
           recurring_expr_1 - 0.00003890933319150247*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_133*pow(theta3,11.666666666666666)*
        (-4.783524166754698e-9 + 3.2568071954548985e-9*theta0theta3_comb*
           recurring_expr_1 - 1.0276508948991564e-9*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + pow(theta3,15)*(3.274130801749477e-12 - 
          2.2753363570435764e-12*theta0theta3_comb*
           recurring_expr_1 + 7.331470811126562e-13*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*pow(theta3,13.333333333333334)*
        (1.3583519401264653e-10 - 9.672943317645832e-11*theta0theta3_comb*
           recurring_expr_1 + 3.1918941998055616e-11*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_2*theta3_pow_5*theta3_pow_5*
        (1.0819756796833929e-7 - 7.68710402064183e-8*theta0theta3_comb*
           recurring_expr_1 + 2.5310253760475244e-8*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*pow(theta3,8.333333333333334)*
        (1.7420618880578452e-6 - 1.3563925798001168e-6*theta0theta3_comb*
           recurring_expr_1 + 4.840150926665589e-7*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + pow(theta0,4)*theta3_pow_5*
        (0.0013337043441791188 - 0.0009363229675156456*theta0theta3_comb*
           recurring_expr_1 + 0.00030474967286614414*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33*
        (0.009194574209833916 - 0.007414061611724472*theta0theta3_comb*
           recurring_expr_1 + 0.002718981485432563*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + pow(theta0,6)*(3.091747888067239 - 
          2.6260550705381083*theta0theta3_comb*
           recurring_expr_1 + 1.*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1)))/
   (pow(theta3,12)*(1.*theta0_pow_2 - 0.0018440070919485247*theta0_pow_133*theta3_pow_166 + 
       0.0018733054209087634*theta0_pow_one_third*theta3_pow_33 + 0.00040623943233705464*theta3_pow_5)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.3333333333333333)*pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 
           16.884262106926414*theta0_pow_133*theta3_pow_166 + 
           17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
      0.6666666666666666)));

	return (norm_merg*norm_ring_theta0 + (norm_ring/norm_merg)*norm_merg_theta0);

	}

static REAL8 XLALSimIMRPhenomBNormalization_Ringdown_Der_theta3(
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_merg,	/**< Normalization for merger waveform */
	const REAL8 norm_ring, /**< Normalization for ringdown waveform */
	const REAL8 norm_merg_theta3 /**< Derviative of norm_merg w.r.t theta0*/
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = pow(theta3,6);
	REAL8	theta0theta3_comb = pow(theta3/(flow*theta0),0.3333333333333333); 

	REAL8	norm_ring_theta3;

	/*Some Recurring Expressions*/

	REAL8 recurring_expr_1 = pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
                 17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
            0.3333333333333333);


	norm_ring_theta3 = ((flow*flow*theta0_pow_2*(theta0_pow_2*theta0_pow_2*theta0_pow_2*(-2.974304021311726e9 + 
          2.6411329904792056e9*theta0theta3_comb*
           recurring_expr_1 - 1.0360148221328855e9*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33*
        (-8.483236005024768e6 + 7.12717031472055e6*theta0theta3_comb*
           recurring_expr_1 - 2.69338389090502e6*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_2*theta0_pow_2*theta3_pow_5*
        (-1.388190823915925e6 + 976106.5716062305*theta0theta3_comb*
           recurring_expr_1 - 318186.7317364536*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*pow(theta3,8.333333333333334)*
        (-1362.9020402713631 + 1122.416311616068*theta0theta3_comb*
           recurring_expr_1 - 418.13575617334664*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_2*pow(theta3,10)*
        (-106.62180006288368 + 76.11616463511724*theta0theta3_comb*
           recurring_expr_1 - 25.17663731273577*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*pow(theta3,13.333333333333334)*
        (-0.07943941508018493 + 0.05800974078945667*theta0theta3_comb*
           recurring_expr_1 - 0.01959503179167246*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + pow(theta3,15)*(-0.0015748794472423468 + 
          0.0010944524459305302*theta0theta3_comb*
           recurring_expr_1 - 0.0003526487913167943*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_133*pow(theta3,11.666666666666666)*
        (4.140774893012413 - 2.7964683432051944*theta0theta3_comb*
           recurring_expr_1 + 0.8749200456549843*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33*theta3_pow_33*
        (163931.3511918715 - 117171.0388392613*theta0theta3_comb*
           recurring_expr_1 + 38800.79496271701*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1) + theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_166*
        (2.387674490408703e8 - 2.2882156955475348e8*theta0theta3_comb*
           recurring_expr_1 + 9.399432106728905e7*theta0theta3_comb*theta0theta3_comb*
           recurring_expr_1*recurring_expr_1)))/
   (theta3_pow_6*theta3_pow_5*theta3*theta3*(1.*theta0_pow_2 - 0.0018440070919485247*theta0_pow_133*theta3_pow_166 + 
       0.0018733054209087634*theta0_pow_one_third*theta3_pow_33 + 0.00040623943233705464*theta3_pow_5)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.3333333333333333)*pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 
           16.884262106926414*theta0_pow_133*theta3_pow_166 + 
           17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
      0.6666666666666666)));

	return (norm_merg*norm_ring_theta3 + (norm_ring/norm_merg)*norm_merg_theta3);

	}


/**
*******************************************************************************************************************************
*/

/**
 * Frequency domain amplitude of the IMRPhenomB waveforms - The inspiral, merger and ringdown parts
 */

static REAL8 XLALSimIMRPhenomBAmplitude_Inspiral(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta3_pow_333 = pow(theta0,0.3333333333333333);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);


	REAL8 coef1, coef2, amp;

	coef1 = 0.016200730203430484*sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334);
	coef2 = (0.0341579447030346*sqrt(theta0_pow_one_third/theta3_pow_166)*sqrt(theta3/(flow*theta0)))/
    (flow*theta3_pow_333*cbrt(theta3)*cbrt(theta3)) - 
   (0.00315938483464183*sqrt(theta0_pow_one_third/theta3_pow_166)*theta3*sqrt(theta3/(flow*theta0)))/(flow*theta0);

	amp = coef1*pow(f,-7.0/6.0) + coef2*pow(f,-1.0/2.0);
	return amp;

	}

static REAL8 XLALSimIMRPhenomBAmplitude_Merger(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);

	REAL8 coef1, coef2, coef3, amp, norm_merger;

	coef1 = (0.016200730203430484*sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334))/
   sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
        ((-1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
          (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3);
	coef2 = coef1*(-0.694943633842542*theta3)/(flow*theta0*cbrt(theta3/(flow*theta0))*cbrt(theta3/(flow*theta0)));
	coef3 = coef1*(0.2239211337314047*theta3)/(flow*theta0*cbrt(theta3/(flow*theta0)));

	norm_merger = XLALSimIMRPhenomBNormalization_Merger(theta0,theta3,flow);

	amp = norm_merger*(coef1*pow(f,-2.0/3.0) + coef2*pow(f,-1.0/3.0) + coef3);
	return amp;
	}

static REAL8 XLALSimIMRPhenomBAmplitude_Ringdown(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = pow(theta3,6); 

	REAL8 amp, waveform, norm_ringdown;

	norm_ringdown = XLALSimIMRPhenomBNormalization_Ringdown(theta0,theta3,flow);

	waveform = ((0.002578426293574129*sqrt(theta0_pow_one_third/theta3_pow_166)*theta3_pow_6*theta3*
     (-11300.842322830984*theta0_pow_2 + 1239.942236495006*theta0_pow_133*theta3_pow_166 - 
       47.849592830686746*theta0_pow_one_third*theta3_pow_33 + 1.8598228509251575*theta3_pow_5))/
   (pow(theta3/(flow*theta0),0.16666666666666666)*pow((flow*theta0*
         (-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      1.1666666666666667)*(3.1927259301371995e7*pow(1.*flow*pow(theta0,3) - 
          0.10972122263753407*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
          0.004234161619441117*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 - 0.00016457382536589995*flow*theta0*theta3_pow_5,2) + 
       pow(f*theta3_pow_6 + flow*(-9156.289138283713*pow(theta0,3) + 
            16.884262106926414*pow(theta0,2.3333333333333335)*theta3_pow_166 - 
            17.15252607815491*pow(theta0,1.6666666666666667)*theta3_pow_33 - 3.719645701850315*theta0*theta3_pow_5),2))));

	amp = norm_ringdown*waveform;
	return amp;

	}
/**
*******************************************************************************************************************************
*/





/**
*******************************************************************************************************************************
*/
/**
 * Derivative of the amplitude with respect to \f$\theta0\f$
 */

/**
*Derivative of inspiral part of IMRPhenomB amplitude w.r.t theta0
*/
static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta0_Inspiral(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);


	REAL8 coef1, coef2, amp;

	coef1 = (-0.008100365101715244*sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334))/theta0;
	coef2 = (-0.017078972351517306*pow(theta0_pow_one_third/theta3_pow_166,1.5)*pow(theta3/(flow*theta0),1.5)*
     (1.*theta0_pow_one_third - 0.21581796011405455*theta3_pow_166))/pow(theta0,1.6666666666666667);
	amp = coef1*pow(f,-7.0/6.0) + coef2*pow(f,-1.0/2.0);
	return amp;

	}

/**
*Derivative of merger part of IMRPhenomB amplitude w.r.t theta0
*/
static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta0_Merger(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_merg,
	const REAL8 norm_merg_theta0,
	const REAL8 amp_merger
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = pow(theta3,6); 

	REAL8	coef1, coef2, coef3, amp_theta0, amp; 

	coef1 = (sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
     (-0.03240146040686097*theta0_pow_2 - 0.000038169805226880874*theta0_pow_133*theta3_pow_166 + 
       0.0000581334338539683*theta0_pow_one_third*theta3_pow_33 + 7.698479477214536e-7*theta3_pow_5))/
   (sqrt((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6)*
     (1.*pow(theta0,3) + 0.001413632771396877*pow(theta0,2.3333333333333335)*theta3_pow_166 - 
       0.0026912413726415843*pow(theta0,1.6666666666666667)*theta3_pow_33 - 0.00004751933635426132*theta0*theta3_pow_5));

	coef2 = coef1*((-0.810767572816299*pow(theta3/(flow*theta0),0.3333333333333333)*
     (1.*theta0_pow_2 + 0.0012116852326258948*theta0_pow_133*theta3_pow_166 - 
       0.0019223152661725599*theta0_pow_one_third*theta3_pow_33 - 0.000027153906488149293*theta3_pow_5))/
   (1.*theta0_pow_2 + 0.0011780273094973974*theta0_pow_133*theta3_pow_166 - 
     0.0017941609150943893*theta0_pow_one_third*theta3_pow_33 - 0.00002375966817713066*theta3_pow_5));

	coef3 = coef1*((0.29856151164187295*pow(theta3/(flow*theta0),0.6666666666666666)*
     (1.*theta0_pow_2 + 0.0012369286749722673*theta0_pow_133*theta3_pow_166 - 
       0.002018431029481188*theta0_pow_one_third*theta3_pow_33 - 0.000029699585221413194*theta3_pow_5))/
   (1.*theta0_pow_2 + 0.0011780273094973974*theta0_pow_133*theta3_pow_166 - 
     0.0017941609150943893*theta0_pow_one_third*theta3_pow_33 - 0.00002375966817713066*theta3_pow_5));


	amp_theta0 = norm_merg*(coef3 + coef1*pow(f,-2.0/3.0) + coef2*pow(f,-1.0/3.0));
	amp = norm_merg_theta0*amp_merger/norm_merg;

	return amp+amp_theta0 ;

	}


/**
*Derivative of ringdown part of IMRPhenomB amplitude w.r.t theta0
*/


static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta0_Ringdown(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_ringdown,
	const REAL8 norm_ringdown_theta0,
	const REAL8 amp_ringdown
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta0_pow_333 = theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third;
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta0_pow_266 = theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third;
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = theta3_pow_5*theta3;
	REAL8	theta3_pow_666 = theta3_pow_33*theta3_pow_33;
	REAL8	theta3_pow_833 = pow(theta3,8.333333333333334);
	REAL8	theta3_pow_116 = pow(theta3,11.666666666666666);
	REAL8	theta3_pow_133 = pow(theta3,13.333333333333334);
 
	REAL8	amp_theta0;


	amp_theta0 = ((-6.309330788947024e-11*flow*sqrt(theta0_pow_one_third/theta3_pow_166)*theta3_pow_6*theta3_pow_5*
     pow(theta3/(flow*theta0),1.8333333333333333)*(f*f*theta3_pow_6*theta3_pow_6*
        (1.2340282338796662e-9*theta0_pow_2*theta0_pow_2 - 2.2527748592813595e-10*theta0_pow_333*theta3_pow_166 + 
          1.3866741910127907e-11*theta0_pow_266*theta3_pow_33 - 4.791014881221891e-13*theta0_pow_2*theta3_pow_5 - 
          1.586443026896594e-14*theta0_pow_133*theta3_pow_666 + 
          7.894767665826844e-16*theta0_pow_one_third*theta3_pow_833 + 6.433761654177981e-18*theta3_pow_5*theta3_pow_5) + 
       f*flow*theta0*theta3_pow_6*(-0.00009039295451366248*theta0_pow_2*theta0_pow_2*theta0_pow_2 + 
          0.000011606998252519058*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_166 - 
          4.790497433533833e-7*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33 + 3.402203750822557e-10*theta0_pow_2*theta0_pow_2*theta3_pow_5 + 
          2.814424003062018e-9*theta0_pow_333*theta3_pow_666 - 
          1.3775867112252942e-10*theta0_pow_266*theta3_pow_833 + 3.1444151526879723e-12*theta0_pow_2*theta3_pow_5*theta3_pow_5 + 
          1.3374286585876486e-13*theta0_pow_133*theta3_pow_116 - 
          8.86454266853857e-15*theta0_pow_one_third*theta3_pow_133 - 1.196565694184624e-16*theta3_pow_5*theta3_pow_5*theta3_pow_5) + 
       pow(flow,2)*theta0_pow_2*(1.*pow(theta0,8) - 0.1760878215836856*pow(theta0,7.333333333333333)*theta3_pow_166 + 
          0.016434354889891443*pow(theta0,6.666666666666667)*theta3_pow_33 - 0.0008012371590087298*theta0_pow_2*theta0_pow_2*theta0_pow_2*theta3_pow_5 + 
          2.396325173672129e-6*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_666 + 
          4.292406610607157e-7*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_833 - 1.9043032274546434e-8*theta0_pow_2*theta0_pow_2*theta3_pow_5*theta3_pow_5 - 
          5.531895676049336e-9*theta0_pow_333*theta3_pow_116 + 
          2.628835629389233e-10*theta0_pow_266*theta3_pow_133 - 2.269291512881063e-12*theta0_pow_2*theta3_pow_5*theta3_pow_5*theta3_pow_5 - 
          2.610954820376939e-13*theta0_pow_133*pow(theta3,16.666666666666668) + 
          2.3048547012029015e-14*theta0_pow_one_third*pow(theta3,18.333333333333332) + 3.78318037515207e-16*theta3_pow_5*theta3_pow_5*theta3_pow_5*theta3_pow_5)))/
   (pow(1.*theta0_pow_2 + 0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5,2)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.16666666666666666)*pow(8.638197637157656e-9*f*f*theta3_pow_6*theta3_pow_6 + 
       f*flow*theta0*theta3_pow_6*(-0.00015818767039890933*theta0_pow_2 + 
          2.916991860744045e-7*theta0_pow_133*theta3_pow_166 - 
          2.963338204792056e-7*theta0_pow_one_third*theta3_pow_33 - 6.426206942557403e-8*theta3_pow_5) + 
       pow(flow,2)*theta0_pow_2*(1.*theta0_pow_2*theta0_pow_2 - 0.06319178654351876*theta0_pow_333*theta3_pow_166 + 
          0.00837150705535333*theta0_pow_266*theta3_pow_33 + 0.00023636648033447988*theta0_pow_2*theta3_pow_5 + 
          0.000016361044697260936*theta0_pow_133*theta3_pow_666 + 
          7.178925895994797e-7*theta0_pow_one_third*theta3_pow_833 + 1.2698581923826034e-7*theta3_pow_5*theta3_pow_5),2)));
	

	return (norm_ringdown*amp_theta0 + (norm_ringdown_theta0*amp_ringdown/norm_ringdown));
	}


/**
*******************************************************************************************************************************
*/





/**
*******************************************************************************************************************************
*/


/**
 * Derivative of the amplitude with respect to \f$\theta3\f$ 
 */

static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta3_Inspiral(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);

	REAL8 amp, coef1;

	coef1 = (-0.03415794470303461*pow(theta0_pow_one_third/theta3_pow_166,1.5)*sqrt(theta3/(flow*theta0))*
     (1.*theta0_pow_one_third + 0.06166227431830129*theta3_pow_166))/(flow*pow(theta0,1.6666666666666667));

	amp = coef1*pow(f,-1.0/2.0);
	return amp;
	}



static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta3_Merger(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_merger,
	const REAL8 norm_merger_theta3,
	const REAL8 amp_merger
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_433 = pow(theta3,4.333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_6 = theta3_pow_5*theta3;

	REAL8 coef1, coef2, coef3, amp, amp_theta3; 

	coef1 = (sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
     (0.04860219061029145*theta0_pow_2 + 0.000049620746794945154*theta0_pow_133*theta3_pow_166 - 
       0.00005813343385396831*theta0_pow_one_third*theta3_pow_33 - 3.8492397386072687e-7*theta3_pow_5))/
   (sqrt((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6)*
     (1.*theta0_pow_2*theta3 + 0.001413632771396877*theta0_pow_133*pow(theta3,2.6666666666666665) - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_433 - 0.00004751933635426132*theta3_pow_6));

	coef2 = coef1*(pow(theta3/(flow*theta0),0.3333333333333333)*(-0.77215959315838*theta0_pow_2 - 
       0.0008186625792278745*theta0_pow_133*theta3_pow_166 + 
       0.0010390339216949634*theta0_pow_one_third*theta3_pow_33 + 9.173127856615659e-6*theta3_pow_5))/
   (1.*theta0_pow_2 + 0.0010209570015644114*theta0_pow_133*theta3_pow_166 - 
     0.001196107276729593*theta0_pow_one_third*theta3_pow_33 - 7.919889392376888e-6*theta3_pow_5);

	coef3 = coef1*(pow(theta3/(flow*theta0),0.6666666666666666)*(0.27368138567171685*theta0_pow_2 + 
       0.0002989565721371088*theta0_pow_133*theta3_pow_166 - 
       0.0004017505462045103*theta0_pow_one_third*theta3_pow_33 - 4.138004760792837e-6*theta3_pow_5))/
   (1.*theta0_pow_2 + 0.0010209570015644114*theta0_pow_133*theta3_pow_166 - 
     0.001196107276729593*theta0_pow_one_third*theta3_pow_33 - 7.919889392376888e-6*theta3_pow_5);

	amp_theta3 = norm_merger*(coef1*pow(f,-2.0/3.0) + coef2*pow(f,-1.0/3.0) + coef3);
	amp = amp_theta3 + (norm_merger_theta3*amp_merger)/norm_merger;
	
	return amp;

	}



static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta3_Ringdown(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow,	/**< Lower Frequency Cut-off */
	const REAL8 norm_ringdown,
	const REAL8 norm_ringdown_theta3,
	const REAL8 amp_ringdown
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta0_pow_333 = theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third;
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta0_pow_266 = theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third;
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_6 = theta3_pow_5*theta3;
	REAL8	theta3_pow_666 = theta3_pow_33*theta3_pow_33;
	REAL8	theta3_pow_833 = pow(theta3,8.333333333333334);
	REAL8	theta3_pow_116 = pow(theta3,11.666666666666666);
	REAL8	theta3_pow_133 = pow(theta3,13.333333333333334);

	REAL8	amp, amp_theta3;


	amp_theta3 = ((1.0148182623362965e-28*theta3_pow_5*theta3*theta3*theta3*(f*f*theta3_pow_6*theta3_pow_6*
        (7.672203604734464e8*theta0_pow_2*theta0_pow_2 - 2.255052633647948e8*theta0_pow_333*theta3_pow_166 + 
          1.9955946557482686e7*theta0_pow_266*theta3_pow_33 - 847297.3107922212*theta0_pow_2*theta3_pow_5 - 
          17276.761525970855*theta0_pow_133*theta3_pow_666 + 
          948.924678555601*theta0_pow_one_third*theta3_pow_833 + 1.*theta3_pow_5*theta3_pow_5) + 
       f*flow*theta0*theta3_pow_6*(-9.834848034582369e13*theta0_pow_2*theta0_pow_2*theta0_pow_2 + 
          1.3397946998688654e13*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_166 - 
          5.98715176465292e11*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_33 + 1.275067347466667e10*theta0_pow_2*theta0_pow_2*theta3_pow_5 + 
          2.3522458859953704e9*theta0_pow_333*theta3_pow_666 - 
          1.726276444423315e8*theta0_pow_266*theta3_pow_833 + 6.218755128271468e6*theta0_pow_2*theta3_pow_5*theta3_pow_5 + 
          116529.57602908649*theta0_pow_133*theta3_pow_116 - 
          9022.053218358858*theta0_pow_one_third*theta3_pow_133 - 52.075039825904355*theta3_pow_5*theta3_pow_5*theta3_pow_5) + 
       pow(flow,2)*theta0_pow_2*(1.154623349117609e18*theta0_pow_2*theta0_pow_2*theta0_pow_2*theta0_pow_2 - 
          2.0514900773449933e17*pow(theta0,7.333333333333333)*theta3_pow_166 + 
          1.8907897656948656e16*pow(theta0,6.666666666666667)*theta3_pow_33 - 9.780388431773455e14*theta0_pow_2*theta0_pow_2*theta0_pow_2*theta3_pow_5 + 
          7.589655976348427e12*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_666 + 
          4.0039626591207336e11*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta0_pow_one_third*theta3_pow_833 - 4.233637953066569e10*theta0_pow_2*theta0_pow_2*theta3_pow_5*theta3_pow_5 - 
          4.087882103027116e9*theta0_pow_333*theta3_pow_116 + 
          3.0551750956245e8*theta0_pow_266*theta3_pow_133 - 8.532293257209271e6*theta0_pow_2*theta3_pow_5*theta3_pow_5*theta3_pow_5 - 
          200782.79369336233*theta0_pow_133*pow(theta3,16.666666666666668) + 
          21313.235117101256*theta0_pow_one_third*pow(theta3,18.333333333333332) + 191.10649228449162*theta3_pow_5*theta3_pow_5*theta3_pow_5*theta3_pow_5)))/
   (pow(flow,2)*pow(theta0_pow_one_third/theta3_pow_166,2.5)*pow(theta3/(flow*theta0),1.1666666666666667)*
     pow(1.*theta0_pow_2 + 0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5,2)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.16666666666666666)*pow(8.638197637157656e-9*f*f*theta3_pow_6*theta3_pow_6 + 
       f*flow*theta0*theta3_pow_6*(-0.00015818767039890933*theta0_pow_2 + 
          2.916991860744045e-7*theta0_pow_133*theta3_pow_166 - 
          2.963338204792056e-7*theta0_pow_one_third*theta3_pow_33 - 6.426206942557403e-8*theta3_pow_5) + 
       pow(flow,2)*theta0_pow_2*(1.*theta0_pow_2*theta0_pow_2 - 0.06319178654351876*theta0_pow_333*theta3_pow_166 + 
          0.00837150705535333*theta0_pow_266*theta3_pow_33 + 0.00023636648033447988*theta0_pow_2*theta3_pow_5 + 
          0.000016361044697260936*theta0_pow_133*theta3_pow_666 + 
          7.178925895994797e-7*theta0_pow_one_third*theta3_pow_833 + 1.2698581923826034e-7*theta3_pow_5*theta3_pow_5),2)));

	amp = norm_ringdown*amp_theta3 + (norm_ringdown_theta3*amp_ringdown/norm_ringdown);

	return amp;
	
	}


/**
*******************************************************************************************************************************
*/



/**
*******************************************************************************************************************************
*/

/**
 * Phase of IMRPhenomB and it's derivatives
 */

/**
* Derivative of IMRPhenomB Phase w.r.t theta0
*/

static REAL8 XLALSimIMRPhenomBPhase_Der_theta0(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_433 = pow(theta3,4.333333333333333);

	REAL8 coef_k_0, coef_k_1, coef_k_2, coef_k_3, coef_k_4, coef_k_5, coef_k_6, coef_k_7, phase; 

	coef_k_0 = (0.6000000000000003*theta3_pow_166)/(pow(theta0,1.6666666666666667)*pow(theta3/(flow*theta0),1.6666666666666667));
	coef_k_1 = 0.0;
	coef_k_2 = (-496796.3153004055*flow*theta0_pow_133)/theta3_pow_433 + 
   (30750.248778482535*flow*theta0_pow_one_third)/pow(theta3,2.6666666666666665) - (433.96861359260646*flow)/theta3 + 
   (0.13291697653291495*flow*pow(theta3,0.6666666666666666))/theta0_pow_one_third;
	coef_k_3 = (flow*pow(theta3/(flow*theta0),0.3333333333333333)*(2.788900203576315e6*theta0_pow_133 - 
       162863.7325707082*theta0_pow_one_third*theta3_pow_166 + 1966.6051704758727*theta3_pow_33))/
   theta3_pow_433;
	coef_k_4 = (flow*pow(theta3/(flow*theta0),0.6666666666666666)*(-6.066903031167109e6*theta0_pow_2 + 
       323306.0335697849*theta0_pow_133*theta3_pow_166 - 
       2664.8290736333292*theta0_pow_one_third*theta3_pow_33 - 0.1101206717478259*theta3_pow_5))/
   (theta0_pow_one_third*theta3_pow_433);
	coef_k_5 = 0.0;
	coef_k_6 = (flow*pow(theta3/(flow*theta0),1.3333333333333333)*(-3.2789137589417044e6*theta0_pow_133 + 
       99788.95476309463*theta0_pow_one_third*theta3_pow_166 + 2556.3746770489306*theta3_pow_33))/
   theta3_pow_433;
	coef_k_7 = (776165.1005429976*flow*pow(theta3/(flow*theta0),1.6666666666666667)*
     (1.*theta0_pow_133 - 0.0023675484834784735*theta3_pow_33))/theta3_pow_433;

	phase = coef_k_0*pow(f,-5.0/3.0) + coef_k_2*pow(f,-1.0) + coef_k_3*pow(f,-2.0/3.0) + coef_k_4*pow(f,-1.0/3.0) + coef_k_6*pow(f,1.0/3.0) + coef_k_7*pow(f,2.0/3.0) + coef_k_5 + coef_k_1;
	return phase; 

	}


/**
* Derivative of IMRPhenomB Phase w.r.t theta3
*/


static REAL8 XLALSimIMRPhenomBPhase_Der_theta3(
	const REAL8 f,		/**<Fourier Frequency*/
	const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
	const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
	const REAL8 flow	/**< Lower Frequency Cut-off */
	){
	REAL8	theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8	theta0_pow_2 = pow(theta0,2);
	REAL8	theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8	theta0_pow_333 = pow(theta0,0.3333333333333333);
	REAL8	theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8	theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8	theta3_pow_5 = pow(theta3,5);
	REAL8	theta3_pow_433 = pow(theta3,4.333333333333333);
	REAL8	theta3_pow_533 = pow(theta3,5.333333333333333);

	REAL8 coef_k_0, coef_k_1, coef_k_2, coef_k_3, coef_k_4, coef_k_5, coef_k_6, coef_k_7, phase; 

	coef_k_0 = 0;
	coef_k_1 = 0;
	coef_k_2 = (922621.7284150383*flow*pow(theta0,2.3333333333333335))/theta3_pow_533 - 
   (49200.39804557207*flow*pow(theta0,1.6666666666666667))/pow(theta3,3.6666666666666665) + (433.9686135926065*flow*theta0)/pow(theta3,2) + 
   (0.2658339530658298*flow*theta0_pow_333)/pow(theta3,0.3333333333333333);
	coef_k_3 = (flow*theta0_pow_333*pow(theta3/(flow*theta0),0.3333333333333333)*
     (-5.577800407152631e6*theta0_pow_2 + 285011.53199873934*theta0_pow_133*theta3_pow_166 - 
       1966.6051704758727*theta0_pow_one_third*theta3_pow_33 - 1.4999999999999998*theta3_pow_5))/theta3_pow_533;
	coef_k_4 = (flow*theta0_pow_333*pow(theta3/(flow*theta0),0.6666666666666666)*
     (1.3347186668567639e7*theta0_pow_2 - 646612.0671395698*theta0_pow_133*theta3_pow_166 + 
       2664.8290736333274*theta0_pow_one_third*theta3_pow_33 + 0.4404826869913036*theta3_pow_5))/theta3_pow_533;
	coef_k_5 = 0.0;
	coef_k_6 = (pow(theta3/(flow*theta0),0.3333333333333333)*(9.836741276825111e6*theta0_pow_133 - 
       399155.8190523788*theta0_pow_one_third*theta3_pow_166 - 2556.3746770489283*theta3_pow_33))/
   theta3_pow_433;
	coef_k_7 = (pow(theta3/(flow*theta0),0.6666666666666666)*(-3.1046604021719894e6*theta0_pow_133 + 
       107462.71729917097*theta0_pow_one_third*theta3_pow_166 + 1837.6085067194913*theta3_pow_33))/
   theta3_pow_433;

	phase = coef_k_0*pow(f,-5.0/3.0) + coef_k_2*pow(f,-1.0) + coef_k_3*pow(f,-2.0/3.0) + coef_k_4*pow(f,-1.0/3.0) + coef_k_6*pow(f,1.0/3.0) + coef_k_7*pow(f,2.0/3.0) + coef_k_5 + coef_k_1;
	return phase; 
	}


/**
*******************************************************************************************************************************
*/











/**
*******************************************************************************************************************************
*/

/**
 * Function to compute the metric elements using waveform derivatives
 */
static REAL8 MetricCoeffs(REAL8Vector *Amp, REAL8Vector *dPsii, REAL8Vector *dPsij,
        REAL8Vector *dAi, REAL8Vector*dAj, REAL8Vector *Sh, REAL8 hSqr, REAL8 df) {
    size_t k = Amp->length;
    REAL8 gij_1 = 0.;
    REAL8 gij_2 = 0.;
    REAL8 gij_3 = 0.;
    REAL8 gij   = 0.;
    for (;k--;) {
        gij_1 += (Amp->data[k]*Amp->data[k]*dPsii->data[k]*dPsij->data[k]
                + dAi->data[k]*dAj->data[k])/(2.0*Sh->data[k]);
        gij_2 += Amp->data[k]*dAi->data[k]/Sh->data[k];
        gij_3 += Amp->data[k]*dAj->data[k]/Sh->data[k];
    }
	
    gij =  df*gij_1/hSqr - df*df*(gij_2*gij_3)/(2.0*hSqr*hSqr) ;
    return gij;

	}
/**
*******************************************************************************************************************************
*/







/**
*******************************************************************************************************************************
*/



/**
 * Compute the template-space metric of IMRPhenomB wavefrom in PN Chirp Time Co-ordinates.
 */
int XLALSimIMRPhenomBMetricTheta0Theta3(
    REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
    REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
    const REAL8 Mass,     /**< Total Mass of the system */
    const REAL8 eta,    /**< Symmetric mass ratio */
    const REAL8 flow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
) {
    REAL8Vector *Amp=NULL, *dATheta0=NULL, *dATheta3=NULL;
    REAL8Vector *dAT0=NULL, *dAPhi=NULL, *dPhaseTheta0=NULL;
    REAL8Vector *dPhaseTheta3=NULL, *dPhaseT0=NULL, *dPhasePhi=NULL;

    /* compute the chirp-time co-ordinates */
    const REAL8 theta0 = ChirpTime_theta0(Mass,eta,flow);
    const REAL8 theta3 = ChirpTime_theta3(Mass,eta,flow);

    /* Compute the transition frequencies */

    const REAL8 fMerg  = TransitionFrequencies_fmerg(theta0,theta3,flow);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
    const REAL8 fRing  = TransitionFrequencies_fring(theta0,theta3,flow);	/**Frequency at which merger part transitions to ringdown part of the waveform*/ 
    const REAL8 fCut   = TransitionFrequencies_fcut(theta0,theta3,flow);	/**Frequency at which ringdown part of the waveform is terminated*/ 

    /*Compute the normalizations and their derivatives*/
    const REAL8	norm_merg		= XLALSimIMRPhenomBNormalization_Merger(theta0,theta3,flow);
    const REAL8	norm_ringdown		= XLALSimIMRPhenomBNormalization_Ringdown(theta0,theta3,flow);
    const REAL8	norm_merg_theta0	= XLALSimIMRPhenomBNormalization_Merger_Der_theta0(theta0,theta3,flow);	    
    const REAL8	norm_merg_theta3	= XLALSimIMRPhenomBNormalization_Merger_Der_theta3(theta0,theta3,flow);
    const REAL8	norm_ringdown_theta0	= XLALSimIMRPhenomBNormalization_Ringdown_Der_theta0(theta0,theta3,flow, norm_merg, norm_ringdown, norm_merg_theta0);
    const REAL8	norm_ringdown_theta3	= XLALSimIMRPhenomBNormalization_Ringdown_Der_theta3(theta0,theta3,flow, norm_merg, norm_ringdown, norm_merg_theta3);

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
	
        dPhaseTheta0->data[k] = XLALSimIMRPhenomBPhase_Der_theta0(f,theta0,theta3,flow);
        dPhaseTheta3->data[k] = XLALSimIMRPhenomBPhase_Der_theta3(f,theta0,theta3,flow);
        dPhaseT0->data[k] = LAL_TWOPI * f;
        dPhasePhi->data[k] = 1.;

		if (f <= fMerg){

        	/* inspiral amplitude of the waveform */
        	Amp->data[k] = XLALSimIMRPhenomBAmplitude_Inspiral(f,theta0,theta3,flow);

        	/* inspiral waveform deratives with respect to the parameters */
        	dATheta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Inspiral(f,theta0,theta3,flow);
        	dATheta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Inspiral(f,theta0,theta3,flow);

		}
		else if ((fMerg<f) && (f<=fRing)){

        	/* merger amplitude of the frequency-domain waveform */
        	Amp->data[k] = XLALSimIMRPhenomBAmplitude_Merger(f,theta0,theta3,flow);

        	/* merger waveform deratives with respect to the parameters */
        	dATheta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Merger(f,theta0,theta3,flow, norm_merg, norm_merg_theta0, Amp->data[k]);
        	dATheta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Merger(f,theta0,theta3,flow, norm_merg, norm_merg_theta3, Amp->data[k]);
		}

		else{
        	/* ringdown amplitude of the frequency-domain waveform */
        	Amp->data[k] = XLALSimIMRPhenomBAmplitude_Ringdown(f,theta0,theta3,flow);

        	/* ringdown waveform deratives with respect to the parameters */
        	dATheta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Ringdown(f,theta0,theta3,flow, norm_ringdown, norm_ringdown_theta0, Amp->data[k]);
        	dATheta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Ringdown(f,theta0,theta3,flow, norm_ringdown, norm_ringdown_theta3, Amp->data[k]);
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


