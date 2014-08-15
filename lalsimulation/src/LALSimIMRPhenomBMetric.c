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

/* Computed in Mathematica with N[\[Pi]^2, 35] */
#define Pi_p2 9.8696044010893586188344909998761511
#define Pi_p3 31.006276680299820175476315067101395
#define Pi_p4 97.409091034002437236440332688705111
#define Pi_p5 306.01968478528145326274131004343561
#define Pi_p4by3 4.6011511144704899609836928301589280
#define Pi_p5by3 6.7388085956981412852178979102191822
#define Pi_p10by3 45.411541289455155012188835024560824
#define Pi_p7by3 14.454942539276981063196177096150014
#define Pi_p1by3 1.4645918875615232630201425272637904
#define Pi_p2by3 2.1450293971110256000774441009412356
#define Pi_11by3 66.509374974201175524807943232081902
#define tenbyPi_p2by3 2.1638812222639821850478371484217674
#define twoPi_p2by3 3.4050219214767547368607032470366207

#define PI 3.1415926535897932384626433832795029

/**
*******************************************************************************************************************************
*/
/**
* Define the Chirp time co-ordinates
*/


static REAL8 ChirpTime_theta0(
	const REAL8 mass,
	const REAL8 eta,
	const REAL8 f_low){
	REAL8 pi;
	pi = PI;
	return (5.0/(128.0*eta*pow(pi*mass*LAL_MTSUN_SI*f_low,5.0/3.0)));
	}

static REAL8 ChirpTime_theta3(
	const REAL8 mass,
	const REAL8 eta,
	const REAL8 f_low){
	REAL8 pi;
	pi = PI;
	return (pow(pi,1.0/3.0)/(4.0*eta*pow(mass*LAL_MTSUN_SI*f_low,2.0/3.0)));
	}
/**
*******************************************************************************************************************************
*/
/**
* Function to calculate the transition frequencies of IMRPhenomB waveform. F_merg, F_ring and F_cut
*/

static REAL8 TransitionFrequencies_fmerg(
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 f_low)
	{
	return ((1.3270087368763253*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. - (1388.908285838913*pow(theta0,2))/pow(theta3,5) - (1.9634062693265484*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) + (3.7378874416543857*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3);


}

static REAL8 TransitionFrequencies_fring(
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 f_low)
	{
	return ((3.719645701850315*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. + (455.39646148015123*pow(theta0,2))/pow(theta3,5) - (0.8397543046176621*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) + (0.8530966599534362*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3);

}

static REAL8 TransitionFrequencies_fcut(
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 f_low)
	{
	return ((6.506364049290605*f_low*theta0)/theta3 + (20.106192982974676*f_low*theta0*(0. + (963.986488648419*pow(theta0,2))/pow(theta3,5) - (9.15298466960777*pow(theta0,1.3333333333333333))/pow(theta3,3.3333333333333335) - (0.7731297368250576*pow(theta0,0.6666666666666666))/pow(theta3,1.6666666666666667)))/theta3);
}

/**
*******************************************************************************************************************************
*/





/**
*******************************************************************************************************************************
*/

/**
 * Frequency domain amplitude of the IMRPhenomB waveforms - The inspiral, merger and ringdown parts
 */

static REAL8 XLALSimIMRPhenomBAmplitude_Inspiral(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta3_pow_six = pow(theta3,6);
	return (0.03580722744181748*sqrt(theta0_pow_one_third/pow(theta3,1.6666666666666667))*pow(theta3/(flow*theta0),0.8333333333333334)*
     (1.*theta0_pow_one_third*pow((f*theta3)/(flow*theta0),0.6666666666666666) + 
       pow(theta3,1.6666666666666667)*(0.4742887882827214 - 0.09249341147745196*pow((f*theta3)/(flow*theta0),0.6666666666666666))))/
   (f*pow(theta3,1.6666666666666667)*pow((flow*theta0*(-27925.658030729737*pow(theta0,2) - 
           39.476625355061934*pow(theta0,1.3333333333333333)*pow(theta3,1.6666666666666667) + 
           75.15468625054058*theta0_pow_one_third*pow(theta3,3.3333333333333335) + 1.3270087368763253*pow(theta3,5)))/theta3_pow_six,
      0.16666666666666666)*pow((f*theta3_pow_six)/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*pow(theta3,1.6666666666666667) + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*pow(theta3,3.3333333333333335) + 1.*flow*theta0*pow(theta3,5)),0.16666666666666666));

	}

static REAL8 XLALSimIMRPhenomBAmplitude_Merger(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);

	return ((-1.4770805194721893e-6*sqrt(theta0_pow_one_third/theta3_pow_166)*theta3_pow_33*
     pow(theta3/(flow*theta0),1.8333333333333333)*(1.*theta0_pow_one_third*
        pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
       (0.4742887882827214 - 0.09249341147745196*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
             (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
             (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*theta3_pow_166)*
     (4.465858060541567 - 3.1035196288177636*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 1.*pow((f*theta3)/(flow*theta0),0.6666666666666666)))/
   ((4.465858060541567 - 3.1035196288177636*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*
     (1.*theta0_pow_2 + 0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/pow(theta3,6),
      0.16666666666666666)*pow((f*pow(theta3,6))/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 + 1.*flow*theta0*theta3_pow_5),0.6666666666666666)));

	}

static REAL8 XLALSimIMRPhenomBAmplitude_Ringdown(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);

	return ((0.009420555376788691*flow*pow(theta0,0.3333333333333333)*pow(theta0_pow_one_third/theta3_pow_166,1.5)*theta3*
     pow(1.*theta0_pow_2 - 0.10972122263753406*theta0_pow_133*theta3_pow_166 + 
       0.004234161619441117*theta0_pow_one_third*theta3_pow_33 - 0.00016457382536589995*theta3_pow_5,2)*
     (0.4742887882827215*theta3_pow_166 + 1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
        (1.*theta0_pow_one_third - 0.09249341147745196*theta3_pow_166))*
     (4.465858060541567 - 3.1035196288177636*pow(theta3/(flow*theta0),0.3333333333333333)*
        pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
              17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
         0.3333333333333333) + 1.*pow(theta3/(flow*theta0),0.6666666666666666)*
        pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
              17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
         0.6666666666666666)))/((4.465858060541567 - 3.1035196288177636*
        pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*
     pow(theta3/(flow*theta0),0.16666666666666666)*sqrt((flow*theta0*
         (-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6)*
     pow((flow*theta0*(9156.289138283713*theta0_pow_2 - 16.884262106926414*theta0_pow_133*theta3_pow_166 + 
           17.15252607815491*theta0_pow_one_third*theta3_pow_33 + 3.719645701850315*theta3_pow_5))/theta3_pow_6,
      0.6666666666666666)*(8.638197637157656e-9*pow(f,2)*pow(theta3,12) + 
       f*flow*theta0*theta3_pow_6*(-0.00015818767039890933*theta0_pow_2 + 
          2.916991860744045e-7*theta0_pow_133*theta3_pow_166 - 
          2.963338204792056e-7*theta0_pow_one_third*theta3_pow_33 - 6.426206942557403e-8*theta3_pow_5) + 
       pow(flow,2)*theta0_pow_2*(1.*pow(theta0,4) - 0.06319178654351876*pow(theta0,3.3333333333333335)*theta3_pow_166 + 
          0.00837150705535333*pow(theta0,2.6666666666666665)*theta3_pow_33 + 0.00023636648033447988*theta0_pow_2*theta3_pow_5 + 
          0.000016361044697260936*theta0_pow_133*pow(theta3,6.666666666666667) + 
          7.178925895994797e-7*theta0_pow_one_third*pow(theta3,8.333333333333334) + 1.2698581923826034e-7*pow(theta3,10)))));

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
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);

	return ((0.7899601461412372*pow(theta0_pow_one_third/theta3_pow_166,1.5)*pow(theta3/(flow*theta0),0.8333333333333334)*
     pow((flow*theta0*(-27925.658030729737*pow(theta0,2) - 39.476625355061934*pow(theta0,1.3333333333333333)*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.8333333333333334)*pow((f*theta3_pow_6)/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 + 1.*flow*theta0*theta3_pow_5),0.8333333333333334)*
     (-0.017078972351517306*theta0_pow_one_third*pow((f*theta3)/(flow*theta0),0.6666666666666666) + 
       theta3_pow_166*(-0.008100365101715244 + 0.003685948973748802*pow((f*theta3)/(flow*theta0),0.6666666666666666))))/
   (pow(f,2)*pow(theta0,1.6666666666666667)));

	}

/**
*Derivative of merger part of IMRPhenomB amplitude w.r.t theta0
*/
static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta0_Merger(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);
	REAL8 theta3_pow_433 = pow(theta3,4.333333333333333);

	return ((1.2075798420906887*sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
     (flow*theta0*pow((f*theta3)/(flow*theta0),0.6666666666666666)*
        (5679.875272238851*pow(theta0,6.666666666666667) - 765.9461020324092*pow(theta0,6)*theta3_pow_166 - 
          41.34268940762763*pow(theta0,5.333333333333333)*theta3_pow_33 + 
          pow(theta0,4.666666666666667)*(4.946675698819638 - 1.2531591850457166e-15*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.3156166343015503*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*theta3_pow_5 + 
          pow(theta0,4)*(0.7880669667388006 - 0.19926726282829244*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.05756694741917788*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,6.666666666666667) + 
          pow(theta0,3.3333333333333335)*(-0.008900892505028683 - 
             0.00006636530752624308*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.0015733337431065527*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,8.333333333333334) + 
          pow(theta0,2.6666666666666665)*(-0.003037565997277338 + 
             0.0007427911992329659*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.0002636629301773663*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,10) + 
          theta0_pow_2*(-0.00003856500755870379 + 9.624879739594921e-6*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.000029880034039211766*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,11.666666666666666) + 
          theta0_pow_133*(2.940967874617337e-6 - 
             5.645669708682708e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             2.715041311884738e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,13.333333333333334) + 
          theta0_pow_one_third*(8.121326836363954e-8 - 
             1.0993408890779027e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             5.298340829647092e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,15) + 
          (5.49309992153288e-10 - 2.1089956641718695e-11*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             6.486870852376095e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,16.666666666666668)) + 
       f*(pow(theta0,6)*pow(theta3,2.6666666666666665)*(619.3767965838991 - 227.63323909215356*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,5.333333333333333)*theta3_pow_433*
           (37.079136752763475 - 14.637392421800635*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,6.666666666666667)*theta3*(-4933.991451827581 + 1907.766165619042*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,4)*pow(theta3,7.666666666666667)*(-0.6357989885982566 + 
             0.14108921353624862*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.04983928024275761*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             0.2332628570793198*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.046301917990925284*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.01922745269684038*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,3.3333333333333335)*pow(theta3,9.333333333333334)*
           (0.007196281693249607 - 0.0007968481598400191*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.0017121231425072405*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             0.0026443991938723582*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.0005283731581482105*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.0007510402043826716*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_133*pow(theta3,14.333333333333334)*
           (-2.488922199889174e-6 + 4.480844874143279e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             2.289460224383871e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             9.453903467320014e-7*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             1.6234044336877867e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             8.674408721170095e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_one_third*pow(theta3,16)*(-7.130054830506004e-8 + 
             8.54951535252515e-9*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             4.638038101847089e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             2.7762818696940254e-8*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             3.047904416115698e-9*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             1.80247868331797e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta3,17.666666666666668)*(-5.089859760706991e-10 + 
             1.954177480823683e-11*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             6.010679469890016e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             2.0500352702158878e-10*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             7.870811669266348e-12*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             2.4209124593887647e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_2*pow(theta3,12.666666666666666)*(0.00003208729522008504 - 
             5.789060657922041e-6*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.000025265930118383232*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             0.000012042483508440517*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             1.5754282231413067e-6*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             9.591342850168162e-6*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,2.6666666666666665)*pow(theta3,11)*(0.0024750907705611983 - 
             0.0005440374976117498*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.0002265509868677819*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             0.0009148466027326221*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.00018426723476536554*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.00008695646595002697*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,4.666666666666667)*theta3_pow_6*(-4.010822872230533 + 
             0.15780831715077542*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.3290036561140074*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             1.4770278609866396*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.101696355122387*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.14134646915458632*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)))))/
   (pow(flow,2)*pow(theta0,3)*pow(4.465858060541567 - 3.1035196288177636*
        pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666),2)*
     pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
       (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
       (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*pow(theta3,5.666666666666667)*
     pow((f*theta3)/(flow*theta0),0.6666666666666666)*pow(1.*theta0_pow_2 + 
       0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5,2)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.16666666666666666)*pow((f*theta3_pow_6)/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 + 1.*flow*theta0*theta3_pow_5),0.6666666666666666)));

	}


/**
*Derivative of ringdown part of IMRPhenomB amplitude w.r.t theta0
*/


static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta0_Ringdown(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666);
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta0_pow_333 = pow(theta0,0.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);


	return ((0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      ((0.23164787794751399*pow((3.719645701850315*flow*theta0)/theta3 + 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*theta3)/
         (flow*theta0_pow_2*pow(theta3/(flow*theta0),0.6666666666666666)) - 
        (0.14928075582093644*pow((3.719645701850315*flow*theta0)/theta3 + 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*theta3)/
         (flow*theta0_pow_2*pow(theta3/(flow*theta0),0.3333333333333333)) - 
        (0.23164787794751399*((3.719645701850315*flow)/theta3 + (20.106192982974676*flow*theta0*
                (0. + (910.7929229603025*theta0)/theta3_pow_5 - (1.1196724061568828*theta0_pow_333)/theta3_pow_33 + 
                  0.5687311066356241/(theta0_pow_333*theta3_pow_166)))/theta3 + 
             (20.106192982974676*flow*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3)*pow(theta3/(flow*theta0),0.3333333333333333))
          /pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666) + 
        (0.14928075582093644*((3.719645701850315*flow)/theta3 + (20.106192982974676*flow*theta0*
                (0. + (910.7929229603025*theta0)/theta3_pow_5 - (1.1196724061568828*theta0_pow_333)/theta3_pow_33 + 
                  0.5687311066356241/(theta0_pow_333*theta3_pow_166)))/theta3 + 
             (20.106192982974676*flow*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3)*pow(theta3/(flow*theta0),0.6666666666666666))
          /pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333))*
      (1. + 0.1352425763914989*(-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.0033751521257146845*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*theta3*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    (flow*theta0_pow_2*(0.25*pow((1.8598228509251575*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*pow(theta3/(flow*theta0),0.16666666666666666)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.0020250912754288105*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      ((1.3270087368763253*flow)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (2777.816571677826*theta0)/theta3_pow_5 - (2.617875025768731*theta0_pow_333)/theta3_pow_33 + 
             2.4919249611029235/(theta0_pow_333*theta3_pow_166)))/theta3 + 
        (20.106192982974676*flow*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
             (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      pow((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3,1.5)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.002700121700571747*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      ((3.719645701850315*flow)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (910.7929229603025*theta0)/theta3_pow_5 - (1.1196724061568828*theta0_pow_333)/theta3_pow_33 + 
             0.5687311066356241/(theta0_pow_333*theta3_pow_166)))/theta3 + 
        (20.106192982974676*flow*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
             (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,1.6666666666666667)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    + (0.008100365101715242*((1.8598228509251575*flow)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1124.1155729879047*theta0)/theta3_pow_5 + (82.22622336941672*theta0_pow_333)/theta3_pow_33 - 
             1.5865623383204817/(theta0_pow_333*theta3_pow_166)))/theta3 + 
        (20.106192982974676*flow*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      ((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.004050182550857621*(0.5*((1.8598228509251575*flow)/theta3 + 
           (20.106192982974676*flow*theta0*(0. - (1124.1155729879047*theta0)/theta3_pow_5 + 
                (82.22622336941672*theta0_pow_333)/theta3_pow_33 - 
                1.5865623383204817/(theta0_pow_333*theta3_pow_166)))/theta3 + 
           (20.106192982974676*flow*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3)*
         ((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3) + 
        2.*((-3.719645701850315*flow)/theta3 - (20.106192982974676*flow*theta0*
              (0. + (910.7929229603025*theta0)/theta3_pow_5 - (1.1196724061568828*theta0_pow_333)/theta3_pow_33 + 
                0.5687311066356241/(theta0_pow_333*theta3_pow_166)))/theta3 - 
           (20.106192982974676*flow*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3)*
         (f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3))*
      pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    (pow(0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2),2)*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    + (0.0013500608502858735*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    (theta0_pow_333*(0.25*pow((1.8598228509251575*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*theta3_pow_166*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))*
      ((-0.23164787794751399*((((1.3270087368763253*flow)/theta3 + 
                  (20.106192982974676*flow*theta0*(0. - (2777.816571677826*theta0)/theta3_pow_5 - 
                       (2.617875025768731*theta0_pow_333)/theta3_pow_33 + 
                       2.4919249611029235/(theta0_pow_333*theta3_pow_166)))/theta3 + 
                  (20.106192982974676*flow*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0) - 
             (1.*((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                     (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0_pow_2)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666) + 
        (0.14928075582093644*((((1.3270087368763253*flow)/theta3 + 
                  (20.106192982974676*flow*theta0*(0. - (2777.816571677826*theta0)/theta3_pow_5 - 
                       (2.617875025768731*theta0_pow_333)/theta3_pow_33 + 
                       2.4919249611029235/(theta0_pow_333*theta3_pow_166)))/theta3 + 
                  (20.106192982974676*flow*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0) - 
             (1.*((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                     (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0_pow_2)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      pow(1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666),2))
     + (0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(0. + (1.4056133797311472*
           pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                    (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                      (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                      (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))/
         (theta0_pow_333*theta3_pow_166) + 
        (0.09016171759433259*(-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
           ((((1.3270087368763253*flow)/theta3 + (20.106192982974676*flow*theta0*
                     (0. - (2777.816571677826*theta0)/theta3_pow_5 - (2.617875025768731*theta0_pow_333)/theta3_pow_33 + 
                       2.4919249611029235/(theta0_pow_333*theta3_pow_166)))/theta3 + 
                  (20.106192982974676*flow*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0) - 
             (1.*((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                     (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0_pow_2)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))));

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
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666) ;
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);

	return ((1.0482840157135205*pow(theta0_pow_one_third/theta3_pow_166,1.5)*pow(theta3/(flow*theta0),0.8333333333333334)*
     (-0.03415794470303461*pow(theta0,2.6666666666666665) - 0.0021545433464636586*theta0_pow_2*theta3_pow_166 + 
       0.00008894980069607421*theta0_pow_133*theta3_pow_33 + 
       7.291607649570083e-6*theta0_pow_one_third*theta3_pow_5 + 1.0008791375326495e-7*pow(theta3,6.666666666666667)))/
   (flow*pow(theta0,1.6666666666666667)*pow((f*theta3)/(flow*theta0),0.3333333333333333)*
     (1.*theta0_pow_2 + 0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.16666666666666666)*pow((f*theta3_pow_6)/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 + 1.*flow*theta0*theta3_pow_5),0.16666666666666666)));

	}



static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta3_Merger(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666) ;
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);
	REAL8 theta3_pow_433 = pow(theta3,4.333333333333333);

	return ((1.2075798420906887*pow(theta0_pow_one_third/theta3_pow_166,1.5)*pow(theta3/(flow*theta0),0.8333333333333334)*
     (flow*theta0*pow((f*theta3)/(flow*theta0),0.6666666666666666)*
        (-5679.875272238851*pow(theta0,6.666666666666667) + 1162.9703542376096*pow(theta0,6)*theta3_pow_166 + 
          31.27924094532669*pow(theta0,5.333333333333333)*theta3_pow_33 + 
          pow(theta0,4.666666666666667)*(-7.418126007613035 - 1.3624844617530238*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.15780831715077337*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*theta3_pow_5 + 
          pow(theta0,4)*(-1.2092118855254956 + 0.4756365645206837*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.059016986975879586*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,6.666666666666667) + 
          pow(theta0,3.3333333333333335)*(0.013526374411924266 + 
             0.0074439277732643114*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.0014087812645133934*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,8.333333333333334) + 
          pow(theta0,2.6666666666666665)*(0.004449887760107708 - 
             0.001616617648769473*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.00028514041283810563*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,10) + 
          theta0_pow_2*(0.0000507674754375469 - 0.000031829922080411194*
              pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.00003583945987914098*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,11.666666666666666) + 
          theta0_pow_133*(-3.5093860094614663e-6 + 
             9.301502070644329e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             3.311123783005197e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,13.333333333333334) + 
          theta0_pow_one_third*(-7.471855301676302e-8 + 
             1.9629223385873027e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             4.9920326925092795e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,15) + 
          (-2.74654996076644e-10 + 1.0544978320859357e-11*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             3.2434354261880464e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666))*pow(theta3,16.666666666666668)) + 
       f*(pow(theta0,6.666666666666667)*theta3*(4933.99145182758 - 1907.766165619043*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,5.333333333333333)*theta3_pow_433*
           (-30.08560730938491 + 12.38397363287465*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,6)*pow(theta3,2.6666666666666665)*(-895.2862731349995 + 316.5353597648054*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,4.666666666666667)*theta3_pow_6*(5.728341530684689 + 
             0.7890415857538726*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             2.784798188990482e-16*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             2.0304378160924874*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.20339271024477312*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.03533661728864687*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,2.6666666666666665)*pow(theta3,11)*(-0.00345657478857744 + 
             0.0011512976258003765*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.00024147662671382614*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             0.0012310952930591349*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.00037993544403006817*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.00009176572821711631*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_2*pow(theta3,12.666666666666666)*(-0.00004056732258965515 + 
             0.0000212203134718764*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.0000294073951671988*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             0.000014774873950192115*pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             6.54760647865073e-6*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) - 
             0.000010925784240633347*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta3,17.666666666666668)*(3.181162350441869e-10 - 
             1.2213609255148012e-11*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             3.7566746686812603e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             1.4350246891511216e-10*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             5.509568168486442e-12*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             1.694638721572135e-10*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_one_third*pow(theta3,16)*(6.678708722112878e-8 - 
             1.4550919658935874e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             4.4251712119489984e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             2.6308514673204913e-8*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             4.981645788551246e-9*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             1.733889817978917e-8*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          theta0_pow_133*pow(theta3,14.333333333333334)*
           (2.883940764059753e-6 - 7.021442300484052e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             2.7037039429341645e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             1.0726711799197896e-6*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             2.442022560910192e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             1.0009163348477206e-7*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,3.3333333333333335)*pow(theta3,9.333333333333334)*
           (-0.010410730897900558 - 0.00433014190900034*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
             0.00036027869742454396*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) + 
             0.003680142346518517*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.0011236189933541154*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.00008328163095897728*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
          pow(theta0,4)*pow(theta3,7.666666666666667)*(0.9284709688342044 - 
             0.33315030033688475*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) - 
             0.050846976001507276*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666) - 
             0.327566104759257*pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.10818684535444212*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333) + 
             0.019552147198332412*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
                (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
                (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*
              pow((f*theta3)/(flow*theta0),0.3333333333333333)))))/
   (pow(flow,2)*pow(theta0,2.6666666666666665)*pow(4.465858060541567 - 
       3.1035196288177636*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.3333333333333333) + 
       1.*pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
          (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
          (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666),2)*
     pow(1.3270087368763253 - (27925.658030729737*theta0_pow_2)/theta3_pow_5 - 
       (39.476625355061934*theta0_pow_133)/theta3_pow_33 + 
       (75.15468625054058*theta0_pow_one_third)/theta3_pow_166,0.6666666666666666)*theta3_pow_5*
     pow((f*theta3)/(flow*theta0),0.6666666666666666)*pow(1.*theta0_pow_2 + 
       0.001413632771396877*theta0_pow_133*theta3_pow_166 - 
       0.0026912413726415843*theta0_pow_one_third*theta3_pow_33 - 0.00004751933635426132*theta3_pow_5,2)*
     pow((flow*theta0*(-27925.658030729737*theta0_pow_2 - 39.476625355061934*theta0_pow_133*theta3_pow_166 + 
           75.15468625054058*theta0_pow_one_third*theta3_pow_33 + 1.3270087368763253*theta3_pow_5))/theta3_pow_6,
      0.16666666666666666)*pow((f*theta3_pow_6)/
       (-21044.06493695328*flow*pow(theta0,3) - 29.748579838281113*flow*pow(theta0,2.3333333333333335)*theta3_pow_166 + 
         56.63465820688478*flow*pow(theta0,1.6666666666666667)*theta3_pow_33 + 1.*flow*theta0*theta3_pow_5),0.6666666666666666)));

	}



static REAL8 XLALSimIMRPhenomBAmplitude_Der_theta3_Ringdown(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666) ;
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_6 = pow(theta3,6);
	REAL8 theta3_pow_433 = pow(theta3,4.333333333333333);

	return ((0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      ((-0.23164787794751399*pow((3.719645701850315*flow*theta0)/theta3 + 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333))/
         (flow*theta0*pow(theta3/(flow*theta0),0.6666666666666666)) + 
        (0.14928075582093644*pow((3.719645701850315*flow*theta0)/theta3 + 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666))/
         (flow*theta0*pow(theta3/(flow*theta0),0.3333333333333333)) - 
        (0.23164787794751399*((-3.719645701850315*flow*theta0)/pow(theta3,2) - 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
             (20.106192982974676*flow*theta0*(0. - (2276.9823074007563*theta0_pow_2)/theta3_pow_6 + 
                  (2.799181015392207*theta0_pow_133)/theta3_pow_433 - 
                  (1.4218277665890604*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*pow(theta3/(flow*theta0),0.3333333333333333))
          /pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666) + 
        (0.14928075582093644*((-3.719645701850315*flow*theta0)/pow(theta3,2) - 
             (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                  (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                  (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
             (20.106192982974676*flow*theta0*(0. - (2276.9823074007563*theta0_pow_2)/theta3_pow_6 + 
                  (2.799181015392207*theta0_pow_133)/theta3_pow_433 - 
                  (1.4218277665890604*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*pow(theta3/(flow*theta0),0.6666666666666666))
          /pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333))*
      (1. + 0.1352425763914989*(-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    + (0.0033751521257146845*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    (flow*theta0*(0.25*pow((1.8598228509251575*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*pow(theta3/(flow*theta0),0.16666666666666666)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.0020250912754288105*((-1.3270087368763253*flow*theta0)/pow(theta3,2) - 
        (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
             (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
        (20.106192982974676*flow*theta0*(0. + (6944.541429194564*theta0_pow_2)/theta3_pow_6 + 
             (6.544687564421828*theta0_pow_133)/theta3_pow_433 - 
             (6.22981240275731*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*
      pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      pow((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3,1.5)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.002700121700571747*((-3.719645701850315*flow*theta0)/pow(theta3,2) - 
        (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
             (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
        (20.106192982974676*flow*theta0*(0. - (2276.9823074007563*theta0_pow_2)/theta3_pow_6 + 
             (2.799181015392207*theta0_pow_133)/theta3_pow_433 - 
             (1.4218277665890604*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*
      pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,1.6666666666666667)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    + (0.008100365101715242*((-1.8598228509251575*flow*theta0)/pow(theta3,2) - 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
        (20.106192982974676*flow*theta0*(0. + (2810.2889324697617*theta0_pow_2)/theta3_pow_6 - 
             (205.5655584235418*theta0_pow_133)/theta3_pow_433 + 
             (3.966405845801204*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*
      ((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.004050182550857621*(0.5*((-1.8598228509251575*flow*theta0)/pow(theta3,2) - 
           (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
           (20.106192982974676*flow*theta0*(0. + (2810.2889324697617*theta0_pow_2)/theta3_pow_6 - 
                (205.5655584235418*theta0_pow_133)/theta3_pow_433 + 
                (3.966405845801204*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*
         ((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3) + 
        2.*((3.719645701850315*flow*theta0)/pow(theta3,2) + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) - 
           (20.106192982974676*flow*theta0*(0. - (2276.9823074007563*theta0_pow_2)/theta3_pow_6 + 
                (2.799181015392207*theta0_pow_133)/theta3_pow_433 - 
                (1.4218277665890604*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*
         (f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3))*
      pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    (pow(0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2),2)*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.0033751521257146845*theta0_pow_one_third*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3,2.6666666666666665)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666)))\
    - (0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(1. + 0.1352425763914989*
         (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))*
      ((-0.23164787794751399*(((1.3270087368763253*flow*theta0)/theta3 + 
                (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                     (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                     (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)/(flow*theta0) + 
             (((-1.3270087368763253*flow*theta0)/pow(theta3,2) - 
                  (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
                  (20.106192982974676*flow*theta0*(0. + (6944.541429194564*theta0_pow_2)/theta3_pow_6 + 
                       (6.544687564421828*theta0_pow_133)/theta3_pow_433 - 
                       (6.22981240275731*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*theta3)/(flow*theta0)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666) + 
        (0.14928075582093644*(((1.3270087368763253*flow*theta0)/theta3 + 
                (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                     (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                     (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)/(flow*theta0) + 
             (((-1.3270087368763253*flow*theta0)/pow(theta3,2) - 
                  (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
                  (20.106192982974676*flow*theta0*(0. + (6944.541429194564*theta0_pow_2)/theta3_pow_6 + 
                       (6.544687564421828*theta0_pow_133)/theta3_pow_433 - 
                       (6.22981240275731*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*theta3)/(flow*theta0)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      pow(1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666),2))
     + (0.004050182550857621*pow((1.8598228509251575*flow*theta0)/theta3 + 
        (20.106192982974676*flow*theta0*(0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
             (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
             (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2)*
      sqrt(theta0_pow_one_third/theta3_pow_166)*pow(theta3/(flow*theta0),0.8333333333333334)*
      (1. - 0.694943633842542*pow((3.719645701850315*flow*theta0)/theta3 + 
           (20.106192982974676*flow*theta0*(0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.3333333333333333)*
         pow(theta3/(flow*theta0),0.3333333333333333) + 0.2239211337314047*
         pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
                (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
                (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
         pow(theta3/(flow*theta0),0.6666666666666666))*(0. - (3.5140334493278687*theta0_pow_one_third*
           pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                    (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                      (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                      (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))/
         pow(theta3,2.6666666666666665) + (0.09016171759433259*
           (-1.4419642857142858 + (15.589913515794665*theta0_pow_one_third)/theta3_pow_166)*
           (((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                   (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                     (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                     (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)/(flow*theta0) + 
             (((-1.3270087368763253*flow*theta0)/pow(theta3,2) - 
                  (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                       (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                       (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/pow(theta3,2) + 
                  (20.106192982974676*flow*theta0*(0. + (6944.541429194564*theta0_pow_2)/theta3_pow_6 + 
                       (6.544687564421828*theta0_pow_133)/theta3_pow_433 - 
                       (6.22981240275731*theta0_pow_one_third)/pow(theta3,2.6666666666666665)))/theta3)*theta3)/(flow*theta0)))/
         pow((((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
                  (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333)))/
    ((0.25*pow((1.8598228509251575*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
              (0. - (562.0577864939523*theta0_pow_2)/theta3_pow_5 + 
                (61.669667527062536*theta0_pow_133)/theta3_pow_33 - 
                (2.3798435074807225*theta0_pow_one_third)/theta3_pow_166))/theta3,2) + 
        pow(f - (3.719645701850315*flow*theta0)/theta3 - (20.106192982974676*flow*theta0*
             (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - 
               (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
               (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,2))*
      pow((3.719645701850315*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. + (455.39646148015123*theta0_pow_2)/theta3_pow_5 - (0.8397543046176621*theta0_pow_133)/theta3_pow_33 + 
             (0.8530966599534362*theta0_pow_one_third)/theta3_pow_166))/theta3,0.6666666666666666)*
      sqrt((1.3270087368763253*flow*theta0)/theta3 + (20.106192982974676*flow*theta0*
           (0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
             (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*
      (1. - 0.694943633842542*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.3333333333333333) + 
        0.2239211337314047*pow((((1.3270087368763253*flow*theta0)/theta3 + 
               (20.106192982974676*flow*theta0*(0. - (1388.908285838913*theta0_pow_2)/theta3_pow_5 - 
                    (1.9634062693265484*theta0_pow_133)/theta3_pow_33 + 
                    (3.7378874416543857*theta0_pow_one_third)/theta3_pow_166))/theta3)*theta3)/(flow*theta0),0.6666666666666666))));

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
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666) ;
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta0_pow_333 = pow(theta0,0.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);

	
	return ((f*flow*theta0_pow_333*theta3*(theta0_pow_2*(2.788900203576317e6 - 6.066903031167111e6*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_one_third*theta3_pow_33*
         (1966.605170475875 - 2664.8290736333274*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta3_pow_5*(-5.163501364588971e-16 - 0.1101206717478258*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_133*theta3_pow_166*
         (-162863.73257070832 + 323306.03356978495*pow((f*theta3)/(flow*theta0),0.3333333333333333))) + 
     pow(f,2)*pow(theta3,2)*(99788.95476309484*theta0_pow_one_third*theta3_pow_166 + 
        theta3_pow_33*(2556.3746770489315 - 1837.6085067194897*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_133*(-3.2789137589417063e6 + 776165.1005429972*pow((f*theta3)/(flow*theta0),0.3333333333333333))) + 
     pow(flow,2)*theta0_pow_133*(-496796.31530040567*theta0_pow_2*pow((f*theta3)/(flow*theta0),0.6666666666666666) + 
        30750.248778482557*theta0_pow_133*theta3_pow_166*pow((f*theta3)/(flow*theta0),0.6666666666666666) - 
        433.96861359260686*theta0_pow_one_third*theta3_pow_33*pow((f*theta3)/(flow*theta0),0.6666666666666666) + 
        theta3_pow_5*(0.6000000000000003 + 0.13291697653291518*pow((f*theta3)/(flow*theta0),0.6666666666666666))))/
   (pow(flow,2)*pow(theta0,3)*theta3_pow_33*pow((f*theta3)/(flow*theta0),1.6666666666666667)));

	}


/**
* Derivative of IMRPhenomB Phase w.r.t theta3
*/


static REAL8 XLALSimIMRPhenomBPhase_Der_theta3(
	const REAL8 f,
	const REAL8 theta0,
	const REAL8 theta3,
	const REAL8 flow){
	REAL8 theta0_pow_one_third = pow(theta0,0.6666666666666666) ;
	REAL8 theta0_pow_2 = pow(theta0,2);
	REAL8 theta0_pow_133 = pow(theta0,1.3333333333333333);
	REAL8 theta0_pow_333 = pow(theta0,0.3333333333333333);
	REAL8 theta3_pow_33 = pow(theta3,3.3333333333333335);
	REAL8 theta3_pow_166 = pow(theta3,1.6666666666666667);
	REAL8 theta3_pow_5 = pow(theta3,5);
	REAL8 theta3_pow_433 = pow(theta3,4.333333333333333);

	
	return ((pow(flow,2)*theta0_pow_133*pow((f*theta3)/(flow*theta0),0.6666666666666666)*
      (922621.7284150389*theta0_pow_2 - 49200.398045572096*theta0_pow_133*theta3_pow_166 + 
        433.9686135926066*theta0_pow_one_third*theta3_pow_33 + 0.26583395306582985*theta3_pow_5) + 
     pow(f,2)*pow(theta3,2)*(theta0_pow_133*
         (9.836741276825117e6 - 3.1046604021719885e6*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta3_pow_33*(-2556.374677048931 + 1837.6085067194917*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_one_third*theta3_pow_166*
         (-399155.8190523791 + 107462.71729917094*pow((f*theta3)/(flow*theta0),0.3333333333333333))) + 
     f*flow*theta0_pow_333*theta3*(theta0_pow_133*theta3_pow_166*
         (285011.5319987395 - 646612.06713957*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta3_pow_5*(-1.5000000000000007 + 0.4404826869913037*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_one_third*theta3_pow_33*
         (-1966.605170475874 + 2664.829073633325*pow((f*theta3)/(flow*theta0),0.3333333333333333)) + 
        theta0_pow_2*(-5.577800407152632e6 + 1.3347186668567646e7*pow((f*theta3)/(flow*theta0),0.3333333333333333))))/
   (pow(flow,2)*theta0_pow_2*theta3_pow_433*pow((f*theta3)/(flow*theta0),1.6666666666666667)));

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
    for (;k--;) {
        gij_1 += df*(Amp->data[k]*Amp->data[k]*dPsii->data[k]*dPsij->data[k]
                + dAi->data[k]*dAj->data[k])/(2.0*Sh->data[k]*hSqr);
    }

    for (;k--;) {
        gij_2 += df*(Amp->data[k]*dAi->data[k])/(Sh->data[k]*hSqr);
    }


    for (;k--;) {
        gij_3 += df*(Amp->data[k]*dAj->data[k])/(Sh->data[k]*hSqr);
    }

    return (gij_1 - (gij_2*gij_3)/2.0 );

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
    printf("starting calculating the metric");
    REAL8Vector *Amp=NULL, *dA_theta0=NULL, *dA_theta3=NULL;
    REAL8Vector *dA_t0=NULL, *dA_phi=NULL, *dPhase_theta0=NULL;
    REAL8Vector *dPhase_theta3=NULL, *dPhase_t0=NULL, *dPhase_phi=NULL;


    REAL8 theta0, theta3;
    REAL8 fMerg, fRing, fCut;

    const REAL8 df = Sh->deltaF;
    REAL8 hSqr = 0.;
    int s = 0;

    /* compute the chirp-time co-ordinates */
    theta0 = ChirpTime_theta0(Mass,eta,flow);
    theta3 = ChirpTime_theta3(Mass,eta,flow);

    /* Compute the transition frequencies */
    fMerg = TransitionFrequencies_fmerg(theta0,theta3,flow);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
    fRing = TransitionFrequencies_fring(theta0,theta3,flow);	/**Frequency at which merger part transitions to ringdown part of the waveform*/ 
    fCut = TransitionFrequencies_fcut(theta0,theta3,flow);	/**Frequency at which ringdown part of the waveform is terminated*/ 

    printf("The Transition frequencies are fMerg = %f, fRing = %f , fCut = %f",fMerg,fRing,fCut);

    /* create a view of the PSD between flow and fCut */
    size_t nBins = (fCut - flow) / df;
    size_t k = nBins;
    REAL8Vector Shdata = {nBins, Sh->data->data + (size_t) (flow / df)}; /* copy the Vector, including its pointer to the actual data */
    /* drop low-frequency samples */
    Shdata.length = nBins;  /* drop high-frequency samples */

    /* allocate memory for various vectors */
    Amp = XLALCreateREAL8Vector(nBins);
    dA_theta0 = XLALCreateREAL8Vector(nBins);
    dA_theta3 = XLALCreateREAL8Vector(nBins);
    dA_t0 = XLALCreateREAL8Vector(nBins);
    dA_phi = XLALCreateREAL8Vector(nBins);
    dPhase_theta0 = XLALCreateREAL8Vector(nBins);
    dPhase_theta3 = XLALCreateREAL8Vector(nBins);
    dPhase_t0 = XLALCreateREAL8Vector(nBins);
    dPhase_phi = XLALCreateREAL8Vector(nBins);

	printf("About to create amplitude, phase and it's derivatives vectors");

    /* create a frequency vector from fLow to fCut with frequency resolution df */
    if (flow<=fCut){
    for (;k--;) {
        const REAL8 f = flow + k * df;
	
        dPhase_theta0->data[k] = XLALSimIMRPhenomBPhase_Der_theta0(f,theta0,theta3,flow);
        dPhase_theta3->data[k] = XLALSimIMRPhenomBPhase_Der_theta3(f,theta0,theta3,flow);
        dPhase_t0->data[k] = LAL_TWOPI * f;
        dPhase_phi->data[k] = 1.;

        dA_t0->data[k] = 0.;
        dA_phi->data[k] = 0.;

	if (f <= fMerg){
        /* compute the inspiral amplitude of the frequency-domain waveform */
        Amp->data[k] = XLALSimIMRPhenomBAmplitude_Inspiral(f,theta0,theta3,flow);

        /* compute the inspiral waveform deratives with respect to the parameters */
        dA_theta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Inspiral(f,theta0,theta3,flow);
        dA_theta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Inspiral(f,theta0,theta3,flow);
	}

	else if ((fMerg<f)&&(f<=fRing)){
        /* compute the merger amplitude of the frequency-domain waveform */
        Amp->data[k] = XLALSimIMRPhenomBAmplitude_Merger(f,theta0,theta3,flow);

        /* compute the merger waveform deratives with respect to the parameters */
        dA_theta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Merger(f,theta0,theta3,flow);
        dA_theta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Merger(f,theta0,theta3,flow);
	}

	else{
        /* compute the ringdown amplitude of the frequency-domain waveform */
        Amp->data[k] = XLALSimIMRPhenomBAmplitude_Ringdown(f,theta0,theta3,flow);

        /* compute the ringdown waveform deratives with respect to the parameters */
        dA_theta0->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta0_Ringdown(f,theta0,theta3,flow);
        dA_theta3->data[k] = XLALSimIMRPhenomBAmplitude_Der_theta3_Ringdown(f,theta0,theta3,flow);
	}
	/* compute the square of the template norm */
        hSqr += Amp->data[k] * Amp->data[k] / Shdata.data[k];
    }}
    hSqr *= 4 * df;
	printf("Finished Creating amplitude, phase and it's derivatives vectors. Norm is %f", hSqr);

    /* allocate memory, and initialize the Fisher matrix */
    gsl_matrix * g = gsl_matrix_calloc (4, 4);

    /* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
    gsl_matrix_set (g, 0,0, MetricCoeffs(Amp, dPhase_theta0, dPhase_theta0, dA_theta0, dA_theta0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,1, MetricCoeffs(Amp, dPhase_theta0, dPhase_theta3, dA_theta0, dA_theta3, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,2, MetricCoeffs(Amp, dPhase_theta0, dPhase_t0, dA_theta0, dA_t0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,3, MetricCoeffs(Amp, dPhase_theta0, dPhase_phi, dA_theta0, dA_phi, &Shdata, hSqr, df));


    gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
    gsl_matrix_set (g, 1,1, MetricCoeffs(Amp, dPhase_theta3, dPhase_theta3, dA_theta3, dA_theta3, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,2, MetricCoeffs(Amp, dPhase_theta3, dPhase_t0, dA_theta3, dA_t0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,3, MetricCoeffs(Amp, dPhase_theta3, dPhase_phi, dA_theta3, dA_phi, &Shdata, hSqr, df));

    gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
    gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
    gsl_matrix_set (g, 2,2, MetricCoeffs(Amp, dPhase_t0, dPhase_t0, dA_t0, dA_t0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 2,3, MetricCoeffs(Amp, dPhase_t0, dPhase_phi, dA_t0, dA_phi, &Shdata, hSqr, df));

    gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
    gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
    gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
    gsl_matrix_set (g, 3,3, MetricCoeffs(Amp, dPhase_phi, dPhase_phi, dA_phi,dA_phi, &Shdata, hSqr, df));
	printf("Finished Computing Fisher Matrix");
    /* free the memory */
    XLALDestroyREAL8Vector(Amp);
    XLALDestroyREAL8Vector(dA_theta0);
    XLALDestroyREAL8Vector(dA_theta3);
    XLALDestroyREAL8Vector(dA_t0);
    XLALDestroyREAL8Vector(dA_phi);
    XLALDestroyREAL8Vector(dPhase_theta0);
    XLALDestroyREAL8Vector(dPhase_theta3);
    XLALDestroyREAL8Vector(dPhase_t0);
    XLALDestroyREAL8Vector(dPhase_phi);

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

