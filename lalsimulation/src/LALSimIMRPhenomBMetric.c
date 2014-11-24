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



static REAL8 ChirpTime_theta0(const REAL8 mass,	/**< Total Mass of system */
							  const REAL8 eta,	/**< Symmetric mass ratio of system */
							  const REAL8 flow	/**< Lower Frequency Cut-off */

){
	REAL8 pi;
	pi = PI;
	return 5.0/(128.0*eta*pow(pi*mass*LAL_MTSUN_SI*flow,5.0/3.0));
}

static REAL8 ChirpTime_theta3(const REAL8 mass,	/**< Total Mass of system */
							  const REAL8 eta,	/**< Symmetric mass ratio of system */
							  const REAL8 flow	/**< Lower Frequency Cut-off */

){
	REAL8 pi;
	pi = PI;
	return pow(pi,1.0/3.0)/(4.0*eta*pow(mass*LAL_MTSUN_SI*flow,2.0/3.0));
}


static REAL8 ChirpTime_theta3S(const REAL8 mass,	/**< Total Mass of system */
							   const REAL8 eta,	/**< Symmetric mass ratio of system */
							   const REAL8 chi, /**< Reduced spin parameter value	*/
							   const REAL8 flow	/**< Lower Frequency Cut-off */
){
	REAL8 pi;
	pi = PI;
	return (pow(pi,1.0/3.0)/(4.0*eta*pow(mass*LAL_MTSUN_SI*flow,2.0/3.0)))*(17022.*eta - 9565.9*chi*eta) ;
}



static REAL8 TransitionFrequencies_fmerg(
										 const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										 const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
										 const REAL8 theta3S,
										 const REAL8 flow	/**< Lower Frequency Cut-off */
)
{	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_3 = theta0*theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);
	/* Some expressions which occur multiple times */
	
	
	REAL8	expr1 = -0.7794457395540408 + (0.00001800104500292942*(theta3_pow_third*theta3_pow_third)*theta3S)/(theta0_pow_third*theta0_pow_third);
	
	
	return  (-27925.658030729744*flow*theta0_pow_3)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) - (4787.010536911997*flow*pow(theta0,2.3333333333333335))/(theta3_pow_third*theta3_pow_2*theta3_pow_2) + (146.91743091386425*flow*(theta0_pow_third*theta0_pow_third*theta0))/pow(theta3,2.6666666666666665) + (20.106192982974676*flow*theta0)/theta3 + (0.04802651168014136*flow*(theta0_pow_third*theta0_pow_third*theta0)*theta3S)/pow(theta3,3.6666666666666665) + (0.0002862373302760312*flow*theta0*theta3S)/theta3_pow_2 - (1.023947441860004e-8*flow*theta0_pow_third*theta3S_pow_2)/(theta3_pow_third*theta3) - (89.57308973915218*flow*theta0*pow(expr1,0.217))/theta3 + (70.79390549305383*flow*theta0*pow(expr1,0.26))/theta3 ;
	
	
}

static REAL8 TransitionFrequencies_fring(
										 const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										 const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
										 const REAL8 theta3S,
										 const REAL8 flow	/**< Lower Frequency Cut-off */
)
{	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_3 = theta0*theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);
	/* Some expressions which occur multiple times */
	
	/* Some combination of above expressions. For optimization purposes, they are calculated here only once and used in the end expression. */
	
	REAL8	t0_term2	=	theta0_pow_third*theta0_pow_third;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0;
	REAL8	t3_term1	=	theta3_pow_third*theta3_pow_third;
	REAL8	t3_term2	=	theta3_pow_third*theta3_pow_2*theta3_pow_2;
	
	
	
	
	
	REAL8	expr1 = -0.7794457395540408 + (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2);
	
	return (9156.289138283717*flow*theta0_pow_3)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) + (188.3965655180371*flow*pow(theta0,2.3333333333333335))/(t3_term2) - (18.010617308324708*flow*(t0_term6))/pow(theta3,2.6666666666666665) + (10.053096491487338*flow*theta0)/theta3 - (0.002076640683206031*flow*(t0_term6)*theta3S)/pow(theta3,3.6666666666666665) + (0.000453297950563177*flow*theta0*theta3S)/theta3_pow_2 - (9.871711453116567e-10*flow*theta0_pow_third*theta3S_pow_2)/(theta3_pow_third*theta3) - (6.333450789637023*flow*theta0*pow(expr1,0.3))/theta3;
	
}

static REAL8 TransitionFrequencies_fcut(
										const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
										const REAL8 theta3, /**< Theta3 component of Chirp-Time Co-ordinate system*/
										const REAL8 theta3S,
										const REAL8 flow	/**< Lower Frequency Cut-off */
)
{	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_3 = theta0*theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = pow(theta0,1.0/3.0);
	REAL8	theta3_pow_third = pow(theta3,1.0/3.0);
	
	REAL8	t3_term2	=	theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0 ;
	
	
	return (19382.098373745248*flow*theta0_pow_3)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) - (29.730187648067748*flow*pow(theta0,2.3333333333333335))/(t3_term2) + (21.13013741081529*flow*(t0_term6))/pow(theta3,2.6666666666666665) + (9.114273716814719*flow*theta0)/theta3 - (0.0015609287637006245*flow*(t0_term6)*theta3S)/pow(theta3,3.6666666666666665) - (0.0009137763972796158*flow*theta0*theta3S)/theta3_pow_2 - (0.0000350525661224564*flow*theta0_pow_third*theta3S)/theta3_pow_third + (5.490716208946671e-9*flow*theta0_pow_third*theta3S_pow_2)/(theta3_pow_third*theta3) + (8.771363873581556e-11*flow*theta3_pow_third*theta3S_pow_2)/theta0_pow_third;
}

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/





/*****************************************************************************************************************************
	Define the normalization constants and their derivatives.
 
 ******************************************************************************************************************************/




static REAL8 * XLALSimIMRPhenomBNormalization(const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
											  const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
											  const REAL8 theta3S,	/**< Theta3S component of new co-ordinate system */
											  const REAL8 flow		/**< Lower Frequency Cut-off */
){
	static REAL8	normalization_list[8];
	
	REAL8	norm_merg, norm_ring, norm_merg_t0, norm_merg_t3, norm_merg_t3S, norm_ring_t0, norm_ring_t3, norm_ring_t3S;
	REAL8	norm_ring_t0_pre, norm_ring_t3_pre, norm_ring_t3S_pre ;
	/* Normalization constants and their derivatives w.r.t the co-ordinates. The normalization constants are to ensure smooth transition between the amplitudes.*/
	
	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_2 = theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);
	REAL8	pi = PI;
	
	/* Some combination of above expressions. For optimization purposes, they are calculated here only once and used in the end expression. */
	REAL8	t0_term2	=	theta0_pow_third*theta0_pow_third;
	REAL8	t0_term3	=	theta0_pow_third*theta0;
	REAL8	t0_term4	=	theta0_pow_third*theta0_pow_2;
	REAL8	t0_term5	=	theta0_pow_third*theta0_pow_third*theta0_pow_2;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0;
	
	REAL8	t3_term1	=	theta3_pow_third*theta3_pow_third;
	REAL8	t3_term2	=	theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term3	=	theta3_pow_third*theta3_pow_third*theta3;
	REAL8	t3_term4	=	theta3_pow_third*theta3_pow_2;
	REAL8	t3_term5	=	theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term6	=	theta3_pow_2*theta3_pow_2;
	REAL8	t3_term7	=	theta3_pow_2*theta3_pow_2*theta3;
	REAL8	t3_term8	=	theta3_pow_third*theta3_pow_2*theta3;
	REAL8	t3_term9	=	theta3_pow_third*theta3_pow_third*theta3_pow_2;
	REAL8	t3_term10	=	theta3_pow_third*theta3_pow_third*t3_term9 ;
	REAL8	t3_term11	=	theta3_pow_third*theta3;
	REAL8	t3_term12	=	pow(theta3,5.666666666666667);
	REAL8	t3_term13	=	pow(theta3,2.6666666666666665);
	REAL8	t3_term14	=	pow(theta3,3.6666666666666665);
	REAL8	t3_term15	=	pow(theta3,5.333333333333333);

	
	
	
	/* Some expressions which occur multiple times. For optimization purposes, they are calculated here only once and used in the end expression. */
	
	REAL8	expr1 = -0.7794457395540408 + (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2);
	REAL8	expr1_pow_217 = pow(expr1,0.217);
	REAL8	expr1_pow_26  = pow(expr1,0.26);
	REAL8	expr1_pow_3   = pow(expr1,0.3);
	REAL8	expr1_pow_7   = pow(expr1,0.7);
	REAL8	expr1_pow_152 = pow(expr1,1.5230000000000001);
	REAL8	expr1_pow_74  = pow(expr1,0.7400000000000001);
	REAL8	expr1_pow_78  = pow(expr1,0.7830000000000001);
	REAL8	expr1_pow_174 = pow(expr1,1.74);
	REAL8	expr1_pow_178 = pow(expr1,1.783);
	REAL8	expr1_pow_45  = pow(expr1,0.45);
	REAL8	expr1_pow_75  = pow(expr1,0.75);
	REAL8	expr1_pow_216 = pow(expr1,0.21699999999999986);




	
	REAL8	expr2	=	(1.7794457395540408 - (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2));
	REAL8	expr2_term1	=	expr2*expr2;
	
	REAL8	expr3	=	 pow(theta3/(flow*theta0),0.3333333333333333);
	REAL8	expr3_pow_2 = expr3*expr3;
	
	REAL8	expr4	=	pow((-27925.658030729744*(t0_term5) - 4787.010536911997*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(146.91743091386425*theta3 + 0.04802651168014136*theta3S) - 1.0239474418600038e-8*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(0.00028623733027603087*theta3S + theta3*(20.106192982974676 - 89.57308973915218*expr1_pow_217 + 70.79390549305383*expr1_pow_26)))/((t0_term2)*(t3_term7)),0.3333333333333333);
	
	REAL8	expr5	=	(-27925.658030729744*(t0_term5) - 4787.010536911997*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(146.91743091386425*theta3 + 0.04802651168014137*theta3S) - 1.023947441860004e-8*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(0.0002862373302760312*theta3S + theta3*(20.106192982974676 - 89.57308973915218*expr1_pow_217 + 70.79390549305383*expr1_pow_26)))/((t0_term2)*(t3_term7));
	REAL8	expr5_term1	=	cbrt(expr5);
	
	REAL8	expr6	=	(20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*expr2)/(t3_term8) - (0.7131980994477979*(t0_term2)*expr2)/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3;
	REAL8	expr6_term1	=	cbrt(expr6);
	
	REAL8	expr7	=	cbrt((20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*expr2)/(t3_term8) + (4.802428957897166*(t0_term2)*expr2)/(t3_term3) - (1.5716375022681992*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3);
	
	REAL8	expr8	=	cbrt((20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*(expr2))/(t3_term8) - (0.7131980994477979*(t0_term2)*(expr2))/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3);
	
	REAL8	expr9	=	1 + (0.25700804969411556*(1.*(t0_term2) - 0.000037469780438674206*(t3_term1)*theta3S)*(expr4))/(t0_term2) - (0.21294320317109386*(1.*(t0_term2) - 0.000020753711643020906*(t3_term1)*theta3S)*(expr4*expr4))/(t0_term2);
	
	REAL8	expr10	=	cbrt((flow*(1.*(theta0_pow_2*theta0) + 0.020575646167651576*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0019670214686668005*theta3 - 2.2679937820260693e-7*theta3S) - 1.0781345263379211e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(0.0010979444117217687*theta3 + 4.950673179027029e-8*theta3S - 0.0006917049793847144*theta3*expr1_pow_3)))/(theta3_pow_2*theta3_pow_2*theta3_pow_2));
	
	REAL8	expr11	=	pow((20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*(expr2))/(t3_term8) + (4.802428957897166*(t0_term2)*(expr2))/(t3_term3) - (1.5716375022681992*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3,0.3333333333333333);
	
	REAL8	expr12 = ((20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*expr2)/(t3_term8) - (0.20456897851272804*(t0_term2)*expr2)/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3);
	
	REAL8	expr13 = pow((-27925.658030729744*(t0_term5) - 4787.010536911997*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(146.91743091386425*theta3 + 0.04802651168014136*theta3S) - 1.0239474418600038e-8*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(0.00028623733027603087*theta3S + theta3*(20.106192982974676 - 89.57308973915218*expr1_pow_217 + 70.79390549305383*expr1_pow_26)))/((t0_term2)*(t3_term7)),1.3333333333333333);
	
	/* Normalization constant at inspiral-merger transition. norm_merg_t0, norm_merg_t3, norm_merg_t3S denotes derivative of the constant w.r.t theta0, theta3 and theta3S co-ordinate.*/
	
	
	norm_merg = ((1 + (26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/((t0_term3)*(t3_term10*t3_term10)) - (0.19501496506452742*(-10.811580890210543*(t0_term2) + 1.*(t3_term3))*(expr4*expr4))/(t3_term3))/(expr9));

	
	norm_merg_t0 = ((-(((t3_term3)*(1 + (26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/((t0_term3)*(t3_term10*t3_term10)) - (0.19501496506452742*(-10.811580890210543*(t0_term2) + 1.*(t3_term3))*(expr4*expr4))/(t3_term3))*(-4784.745937935112*pow(theta0,3.3333333333333335) + theta0_pow_2*(t3_term4)*(8.390880530053826*theta3 + 0.023231413328905076*theta3S) + (t0_term5)*(t3_term1)*(-546.8001075364147*theta3 + 0.17928337974926692*theta3S) + t3_term15*theta3S_pow_2*(-2.191255533534397e-14*theta3S - (7.48773329722063e-10*theta3)/expr1_pow_78 + (7.090588759353495e-10*theta3)/expr1_pow_74) + (t0_term3)*(theta3_pow_2*theta3)*(-0.00031440445114826307*theta3*theta3S - 1.027771105949522e-7*theta3S_pow_2 - 1.2174346849468309e-15*theta3_pow_2*expr1_pow_216) + (t0_term2)*(t3_term5)*theta3S*(5.848060778260404e-10*theta3S + theta3*(0.000019983392508732752/expr1_pow_78 - 0.000018923486277050577/expr1_pow_74 + 4.5617010343384256e-20*expr1_pow_216)) - 0.17928337974926692*(t3_term1)*theta3S*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) + (7928.77209563223*(1.*(t0_term2) - 0.000020753711643020906*(t3_term1)*theta3S)*(1.*(t0_term5)*expr1_pow_152 + 0.11427986242721798*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0017536731602670974*(t0_term3)*t3_term10*expr1_pow_152 - 5.732662488775536e-7*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 2.9118441129277097e-24*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 - 1.2222301568605703e-13*(t3_term5)*theta3S_pow_2*expr1_pow_152 + 2.544408210464405e-19*(t0_term2)*(t3_term7)*expr1_pow_174 + (theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2*theta3)*theta3S*(-4.176479329925867e-9*expr1_pow_74 + 3.954961563793527e-9*expr1_pow_78))*(expr4))/expr1_pow_152 - 2.9462412233027313e-6*(t0_term2)*(theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2*theta3)*theta3S*expr13))/(expr4*expr4)) + (expr9)*((78940.06029149327*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5)*expr1_pow_152 + 0.13332650616508768*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0029227886004451627*(t0_term3)*t3_term10*expr1_pow_152 - 9.554437481292567e-7*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 3.4166587344280525e-9*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 + 4.074100522868567e-14*(t3_term5)*theta3S_pow_2*expr1_pow_152 + (theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2*theta3)*theta3S*(-2.7843195532839116e-9*expr1_pow_74 + 2.636641042529018e-9*expr1_pow_78) + (t0_term2)*(t3_term7)*(-0.0002399966482538456*expr1_pow_152 + 0.0010691850679708822*expr1_pow_174 - 0.0008450281985017903*expr1_pow_178)))/expr1_pow_152 + 0.17745895018564842*(t3_term1)*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*theta3S*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) + 17542.235620331838*(t0_term2)*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) - 26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) - (78505.3571315806*(t0_term2)*(1.*(t0_term2) - 0.09249341147745195*(t3_term3))*(1.*(t0_term5)*expr1_pow_152 + 0.11427986242721798*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0017536731602670974*(t0_term3)*t3_term10*expr1_pow_152 - 5.732662488775536e-7*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 2.9118441129277097e-24*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 - 1.2222301568605703e-13*(t3_term5)*theta3S_pow_2*expr1_pow_152 + 2.544408210464405e-19*(t0_term2)*(t3_term7)*expr1_pow_174 + (theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2*theta3)*theta3S*(-4.176479329925867e-9*expr1_pow_74 + 3.954961563793527e-9*expr1_pow_78)))/(expr1_pow_152*expr5_term1) + 1.4056133797311474*theta0_pow_2*(t3_term7)*(expr5_term1*expr5_term1)))/((t0_term4)*(t3_term10*t3_term10)*(expr9*expr9))) ;
	
	norm_merg_t3 = ((-(((1 + (26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/((t0_term3)*(t3_term10*t3_term10)) - (0.19501496506452742*(-10.811580890210543*(t0_term2) + 1.*(t3_term3))*(expr4*expr4))/(t3_term3))*(11961.86484483778*theta0_pow_2 + (t0_term3)*(t3_term1)*(1367.000268841037*theta3 - 0.44820844937316734*theta3S) + (t0_term2)*(t3_term4)*(-20.977201325134562*theta3 - 0.06219293335577098*theta3S) + (t3_term15*theta3S_pow_2*(-1.0956277667671992e-14*theta3S + (7.487733297220631e-10*theta3)/expr1_pow_78 - (7.090588759353495e-10*theta3)/expr1_pow_74))/(t0_term3) + (theta3_pow_2*theta3)*(0.0007614893618693598*theta3*theta3S + 4.111084423798093e-7*theta3S_pow_2 + 1.217434684946831e-15*theta3_pow_2*expr1_pow_216) + ((t3_term5)*theta3S*(1.2112282269501943e-9*theta3S + theta3*(-0.000019983392508732752/expr1_pow_78 + 0.000018923486277050577/expr1_pow_74 - 4.561701034338426e-20*expr1_pow_216)))/(t0_term2) + (0.17928337974926692*(t3_term1)*theta3S*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/(t0_term3) - (19821.930239080575*(1.*(t0_term2) - 0.000020753711643020906*(t3_term1)*theta3S)*(1.*(t0_term5)*expr1_pow_152 + 0.114279862427218*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0017536731602670972*(t0_term3)*t3_term10*expr1_pow_152 - 9.172259982040867e-7*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 2.0499952406568365e-9*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 + 2.444460313721142e-14*(t3_term5)*theta3S_pow_2*expr1_pow_152 + 1.017763284185762e-19*(t0_term2)*(t3_term7)*expr1_pow_174 + t3_term12*theta3S*(-1.6705917319703469e-9*expr1_pow_74 + 1.5819846255174108e-9*expr1_pow_78))*(expr4))/((t0_term3)*expr1_pow_152) + (2.9462412233027313e-6*t3_term12*theta3S*expr13)/(t0_term2)))/((theta3_pow_2*theta3_pow_2*theta3_pow_2)*(expr4*expr4))) + ((expr9)*((-157880.12058298654*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5)*expr1_pow_152 + 0.12380318429615285*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0023382308803561297*(t0_term3)*t3_term10*expr1_pow_152 - 1.0509881229421826e-6*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 3.416658734428056e-9*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 + 8.148201045737137e-14*(t3_term5)*theta3S_pow_2*expr1_pow_152 + t3_term12*theta3S*(-1.3921597766419558e-9*expr1_pow_74 + 1.318320521264509e-9*expr1_pow_78) + (t0_term2)*(t3_term7)*(-0.0001199983241269228*expr1_pow_152 + 0.0005345925339854411*expr1_pow_174 - 0.00042251409925089513*expr1_pow_178)))/expr1_pow_152 - 0.17745895018564842*(t3_term1)*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*theta3S*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) - 43855.58905082959*(t0_term2)*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) + 26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) + (196263.3928289515*(t0_term2)*(1.*(t0_term2) - 0.09249341147745195*(t3_term3))*(1.*(t0_term5)*expr1_pow_152 + 0.114279862427218*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.0017536731602670972*(t0_term3)*t3_term10*expr1_pow_152 - 9.172259982040867e-7*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 2.0499952406568365e-9*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 + 2.444460313721142e-14*(t3_term5)*theta3S_pow_2*expr1_pow_152 + 1.017763284185762e-19*(t0_term2)*(t3_term7)*expr1_pow_174 + t3_term12*theta3S*(-1.6705917319703469e-9*expr1_pow_74 + 1.5819846255174108e-9*expr1_pow_78)))/(expr1_pow_152*expr5_term1) - 3.514033449327869*theta0_pow_2*(t3_term7)*(expr5_term1*expr5_term1)))/((t0_term3)*pow(theta3,7.666666666666667)))/(expr9*expr9));
	
	norm_merg_t3S = ((-(((t3_term1)*(1 + (26313.353430497755*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(t0_term2) - 0.000010116096604013777*(t3_term1)*theta3S)*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/((t0_term3)*(t3_term10*t3_term10)) - (0.19501496506452742*(-10.811580890210543*(t0_term2) + 1.*(t3_term3))*(expr4*expr4))/(t3_term3))*((0.12258447263984441*(t0_term2)*(0.6988597173292637 - (0.00002618612016576143*(t3_term1)*theta3S)/(t0_term2))*(0.04802651168014136*(t0_term3) + 0.0002862373302760312*(t0_term2)*(t3_term3) + (t3_term4)*(-2.047894883720008e-8*theta3S - (0.0003498928006197632*theta3)/expr1_pow_78 + (0.00033133471246553445*theta3)/expr1_pow_74)))/t3_term10 + (0.2689250696239004*(t0_term5) + 0.04609907993977405*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0014148200301671839*theta3 - 4.624969976772854e-7*theta3S) + 9.860649900904788e-14*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-2.756475564111513e-9*theta3S + theta3*(-0.00019362334602350447 + 0.0008625920065347126*expr1_pow_217 - 0.0006817478013487593*expr1_pow_26)))/(t3_term7) - (0.006817946156202168*(1.*(t0_term2) - 0.000020753711643020906*(t3_term1)*theta3S)*(1.*(t0_term3)*expr1_pow_152 + 0.005959985854945788*(t0_term2)*(t3_term3)*expr1_pow_152 - 4.264092502405914e-7*(t3_term4)*theta3S*expr1_pow_152 + t3_term10*(-0.007285409420323183*expr1_pow_74 + 0.006898996010208653*expr1_pow_78))*(expr4))/(t3_term10*expr1_pow_152) + 4.419361834954097e-6*(t0_term2)*expr13))/((t0_term3)*(expr4*expr4))) + ((expr9)*(0.04973591971621729*theta0*(3.375 - (10.646770205908554*(t0_term2))/(t3_term3))*t3_term10*(expr2)*(0.04802651168014136*(t0_term3) + 0.0002862373302760312*(t0_term2)*(t3_term3) + (t3_term4)*(-2.047894883720008e-8*theta3S - (0.0003498928006197632*theta3)/expr1_pow_78 + (0.00033133471246553445*theta3)/expr1_pow_74)) - 0.26618842527847253*(1.*(t0_term2) - 0.31699754336080277*(t3_term3))*(1.*(theta0_pow_2*theta0) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.6666904705817113e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))) + (0.06750670739942093*theta0*(t3_term3)*(1.*(t0_term2) - 0.09249341147745195*(t3_term3))*(1.*(t0_term3)*expr1_pow_152 + 0.005959985854945788*(t0_term2)*(t3_term3)*expr1_pow_152 - 4.264092502405914e-7*(t3_term4)*theta3S*expr1_pow_152 + t3_term10*(-0.007285409420323182*expr1_pow_74 + 0.006898996010208652*expr1_pow_78)))/(expr1_pow_152*expr5_term1)))/((t0_term6)*(theta3_pow_2*theta3_pow_2*theta3_pow_2)))/(expr9*expr9));
	
	
	/* Normalization constant at merger-ringdown transition. norm_ring_t0, norm_ring_t3, norm_ring_t3S denotes derivative of the constant w.r.t theta0, theta3 and theta3S co-ordinate.*/
	norm_ring = ((-3.150917089945074e6*flow*(1.*(t0_term5) - 0.10756774166639459*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(0.0015833529729131952*theta3 - 2.1784821538952197e-8*theta3S) - 3.3754960608409787e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(6.018341074822566e-8*theta3S - 0.00044479412261054025*theta3*expr1_pow_45 + 0.0002802202972446404*theta3*expr1_pow_75))*pow(-((flow*(1.*(theta0_pow_2*theta0) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.6666904705817113e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/(theta3_pow_2*theta3_pow_2*theta3_pow_2)),0.6666666666666666)*(-0.009083267631449107*(t3_term1)*theta3S*(-0.10415860320817937*expr3*expr10 + 1.*expr3_pow_2*(expr10*expr10)) + (t0_term2)*(-4.696087900943841 - 25.2496934324519*expr3*expr10 + 437.6695497985125*expr3_pow_2*(expr10*expr10)))*(311.7639249918997*(theta0_pow_2*theta0_pow_2) + pow(theta0,3.3333333333333335)*(t3_term1)*(-45.385890644188066*theta3 - 0.0031538339828645627*theta3S) + (t0_term5)*(t3_term4)*(-18.581339730539014*theta3 - 0.00007704315323053715*theta3S) + 3.6658017555104676e-16*pow(theta3,7)*(theta3S*theta3S*theta3S) + (t0_term2)*t3_term15*theta3S*(-4.6484806500890864e-11*theta3*theta3S - 1.1564133010766257e-15*theta3S_pow_2 + theta3_pow_2*(-7.198154370085134e-7 + 3.2067777718729266e-6*expr1_pow_217 - 2.5344701537069753e-6*expr1_pow_26)) + (t0_term3)*(t3_term5)*(-1.5727410053998532e-9*theta3S_pow_2 + theta3*theta3S*(-1.9760268536044293e-6 - 0.000010116096604013777*expr1_pow_217 + 7.995235946741304e-6*expr1_pow_26) + theta3_pow_2*(0.08300357898163133 - 0.31699754336080277*expr1_pow_217 + 0.25053835020726967*expr1_pow_26 - 0.0023105618636280823*(expr5_term1*expr5_term1))) + theta0_pow_2*(theta3_pow_2*theta3)*(0.00035474001051521043*theta3*theta3S + 5.4239597308179735e-9*theta3S_pow_2 + theta3_pow_2*(0.29547123774198014 + 1.*expr1_pow_217 - 0.790347923681257*expr1_pow_26 + 0.024980826490450633*(expr5_term1*expr5_term1)))))/(theta0*pow(theta3,12.666666666666666)*pow((flow*(1.*(theta0_pow_2*theta0) + 0.020575646167651573*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0019670214686668*theta3 - 2.2679937820260699e-7*theta3S) - 1.0781345263379211e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(4.950673179027029e-8*theta3S + theta3*(0.0010979444117217689 - 0.0006917049793847144*expr1_pow_3))))/(theta3_pow_2*theta3_pow_2*theta3_pow_2),0.6666666666666666)*(-0.000020753711643020906*(t3_term1)*theta3S*(-2.1790556086273494*expr5_term1 + 1.*(expr5_term1*expr5_term1)) + (t0_term2)*(-4.696087900943841 - 1.2069323926137094*expr5_term1 + 1.*(expr5_term1*expr5_term1)))));
	
	norm_ring_t0_pre = ((pi*(expr7*expr7)*((20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*expr2)/(t3_term8) - (0.20456897851272804*(t0_term2)*expr2)/(t3_term3) + (0.5854949302689275*(t0_term2)*(expr2*expr2))/(t3_term3)))/theta3 + (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3)*((0.12258447263984441*expr3*(-1.8897 + 1.4547*expr2)*((0.00002280174653761898*flow*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((910.7929229603028*theta0)/(t3_term7) - (1.119672406156883*theta0_pow_third)/(t3_term8) + 0.5687311066356242/(theta0_pow_third*(t3_term3)) + (0.00006885575619957719*theta3S)/(theta0_pow_third*t3_term9) - (8.558874056109027e-6*theta3S)/(theta0*theta3) + (7.650195440139378*theta0_pow_third*expr2)/(t3_term8) - (0.47546539963186524*expr2)/(theta0_pow_third*(t3_term3)) - (3.6366677468926088e-6*theta3S*expr2)/(theta0*theta3) - (0.10101268416085821*expr2_term1)/(theta0_pow_third*(t3_term3))))/theta3 + (20.106192982974676*flow*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*expr2)/(t3_term8) - (0.7131980994477979*(t0_term2)*expr2)/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (10.053096491487338*flow*(1 - 0.63*expr1_pow_3))/theta3))/(expr6_term1*expr6_term1) + (0.09016171759433259*expr3_pow_2*(1.6557 - 1.8153*expr2)*((0.00002280174653761898*flow*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((910.7929229603028*theta0)/(t3_term7) - (1.119672406156883*theta0_pow_third)/(t3_term8) + 0.5687311066356242/(theta0_pow_third*(t3_term3)) + (0.00006885575619957719*theta3S)/(theta0_pow_third*t3_term9) - (8.558874056109027e-6*theta3S)/(theta0*theta3) + (7.650195440139378*theta0_pow_third*expr2)/(t3_term8) - (0.47546539963186524*expr2)/(theta0_pow_third*(t3_term3)) - (3.6366677468926088e-6*theta3S*expr2)/(theta0*theta3) - (0.10101268416085821*expr2_term1)/(theta0_pow_third*(t3_term3))))/theta3 + (20.106192982974676*flow*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*expr2)/(t3_term8) - (0.7131980994477979*(t0_term2)*expr2)/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (10.053096491487338*flow*(1 - 0.63*expr1_pow_3))/theta3))/expr6_term1 + (6.4200234620069195e-6*(t3_term1)*expr3*theta3S*expr6_term1)/(t0_term6) - (0.12258447263984441*theta3*(-1.8897 + 1.4547*expr2)*expr6_term1)/(flow*theta0_pow_2*expr3_pow_2) - (2.9462412233027313e-6*(t3_term1)*expr3_pow_2*theta3S*(expr6_term1*expr6_term1))/(t0_term6) - (0.09016171759433259*theta3*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1))/(flow*theta0_pow_2*expr3)))/(2.*(expr6_term1*expr6_term1)) + (pi *(expr7*expr7)*((0.00001140087326880949*flow*theta3S)/((t0_term2)*theta3_pow_third*pow(expr1,0.24999999999999994)) + (20.106192982974676*flow*theta0*((-1124.1155729879051*theta0)/(t3_term7) + (82.22622336941672*theta0_pow_third)/(t3_term8) - 1.5865623383204819/(theta0_pow_third*(t3_term3)) - (8.162885715566168e-6*theta3S)/(theta0_pow_third*t3_term9) - (2.4549702589406124e-6*theta3S)/(theta0*theta3) - (0.9069346489870753*theta0_pow_third*expr2)/(t3_term8) - (0.13637931900848535*expr2)/(theta0_pow_third*(t3_term3)) + (0.000014052694118343981*theta3S*expr2)/(theta0*theta3) + (0.39032995351261834*expr2_term1)/(theta0_pow_third*(t3_term3))))/theta3 + (20.106192982974676*flow*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*expr2)/(t3_term8) - (0.20456897851272804*(t0_term2)*expr2)/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 - (0.000027144936354308312*flow*theta3S*(1 - 0.63*expr1_pow_3))/((t0_term2)*theta3_pow_third*pow(expr1,0.55)) + (5.026548245743669*flow*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3)*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(2.*(expr6_term1*expr6_term1)) - (pi *(expr7*expr7)*((0.00002280174653761898*flow*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((910.7929229603028*theta0)/(t3_term7) - (1.119672406156883*theta0_pow_third)/(t3_term8) + 0.5687311066356242/(theta0_pow_third*(t3_term3)) + (0.00006885575619957719*theta3S)/(theta0_pow_third*t3_term9) - (8.558874056109027e-6*theta3S)/(theta0*theta3) + (7.650195440139378*theta0_pow_third*expr2)/(t3_term8) - (0.47546539963186524*expr2)/(theta0_pow_third*(t3_term3)) - (3.6366677468926088e-6*theta3S*expr2)/(theta0*theta3) - (0.10101268416085821*expr2_term1)/(theta0_pow_third*(t3_term3))))/theta3 + (20.106192982974676*flow*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*expr2)/(t3_term8) - (0.7131980994477979*(t0_term2)*expr2)/(t3_term3) - (0.15151902624128732*(t0_term2)*(expr2*expr2))/(t3_term3)))/theta3 + (10.053096491487338*flow*(1 - 0.63*expr1_pow_3))/theta3)*expr12*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(3.*pow(expr6,1.6666666666666667)) + (pi *((20.106192982974676*flow*theta0*((-2777.8165716778267*theta0)/(t3_term7) - (2.6178750257687318*theta0_pow_third)/(t3_term8) + 2.491924961102924/(theta0_pow_third*(t3_term3)) - (0.0015924284861156552*theta3S)/(theta0_pow_third*t3_term9) + (0.000057632493196318874*theta3S)/(theta0*theta3) - (176.92622687810731*theta0_pow_third*expr2)/(t3_term8) + (3.201619305264777*expr2)/(theta0_pow_third*(t3_term3)) - (0.00003772148987549525*theta3S*expr2)/(theta0*theta3) - (1.047758334845466*expr2_term1)/(theta0_pow_third*(t3_term3))))/theta3 + (20.106192982974676*flow*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*expr2)/(t3_term8) + (4.802428957897166*(t0_term2)*expr2)/(t3_term3) - (1.5716375022681992*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (20.106192982974676*flow*theta0*((0.000011601493493937981*(t3_term1)*theta3S)/((t0_term6)*expr1_pow_78) - (0.000010986157772254511*(t3_term1)*theta3S)/((t0_term6)*expr1_pow_74)))/theta3 + (20.106192982974676*flow*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3)*expr12*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(3.*expr7*(expr6_term1*expr6_term1)));
	norm_ring_t0 = norm_merg*norm_ring_t0_pre + (norm_ring/norm_merg)*norm_merg_t0;
	
	norm_ring_t3_pre = ((PI*(expr11*expr11)*((20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*(expr2))/(t3_term8) - (0.20456897851272804*(t0_term2)*(expr2))/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3)*((0.12258447263984441*expr3*(-1.8897 + 1.4547*(expr2))*((-0.00002280174653761898*flow*theta0_pow_third*theta3S)/(t3_term11*expr1_pow_7) + (20.106192982974676*flow*theta0*((-2276.9823074007572*theta0_pow_2)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) + (2.799181015392208*(t0_term3))/(t3_term2) - (1.4218277665890606*(t0_term2))/t3_term9 - (0.00006885575619957719*(t0_term2)*theta3S)/t3_term14 + (8.558874056109027e-6*theta3S)/theta3_pow_2 - (19.125488600348447*(t0_term3)*(expr2))/(t3_term2) + (1.1886634990796632*(t0_term2)*(expr2))/t3_term9 + (3.6366677468926088e-6*theta3S*(expr2))/theta3_pow_2 + (0.25253171040214556*(t0_term2)*expr2_term1)/t3_term9))/theta3 - (20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*(expr2))/(t3_term8) - (0.7131980994477979*(t0_term2)*(expr2))/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3_pow_2 - (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3_pow_2))/(expr8*expr8) + (0.09016171759433259*expr3_pow_2*(1.6557 - 1.8153*(expr2))*((-0.00002280174653761898*flow*theta0_pow_third*theta3S)/(t3_term11*expr1_pow_7) + (20.106192982974676*flow*theta0*((-2276.9823074007572*theta0_pow_2)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) + (2.799181015392208*(t0_term3))/(t3_term2) - (1.4218277665890606*(t0_term2))/t3_term9 - (0.00006885575619957719*(t0_term2)*theta3S)/t3_term14 + (8.558874056109027e-6*theta3S)/theta3_pow_2 - (19.125488600348447*(t0_term3)*(expr2))/(t3_term2) + (1.1886634990796632*(t0_term2)*(expr2))/t3_term9 + (3.6366677468926088e-6*theta3S*(expr2))/theta3_pow_2 + (0.25253171040214556*(t0_term2)*expr2_term1)/t3_term9))/theta3 - (20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*(expr2))/(t3_term8) - (0.7131980994477979*(t0_term2)*(expr2))/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3_pow_2 - (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3_pow_2))/expr8 - (6.4200234620069195e-6*expr3*theta3S*expr8)/((t0_term2)*theta3_pow_third) + (0.12258447263984441*(-1.8897 + 1.4547*(expr2))*expr8)/(flow*theta0*expr3_pow_2) + (2.9462412233027313e-6*expr3_pow_2*theta3S*(expr8*expr8))/((t0_term2)*theta3_pow_third) + (0.09016171759433259*(1.6557 - 1.8153*(expr2))*(expr8*expr8))/(flow*theta0*expr3)))/(2.*(expr8*expr8)) + (PI*(expr11*expr11)*((-0.00001140087326880949*flow*theta0_pow_third*theta3S)/(t3_term11*pow(expr1,0.24999999999999994)) + (20.106192982974676*flow*theta0*((2810.2889324697626*theta0_pow_2)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) - (205.56555842354183*(t0_term3))/(t3_term2) + (3.966405845801205*(t0_term2))/t3_term9 + (8.162885715566168e-6*(t0_term2)*theta3S)/t3_term14 + (2.4549702589406124e-6*theta3S)/theta3_pow_2 + (2.2673366224676883*(t0_term3)*(expr2))/(t3_term2) + (0.3409482975212134*(t0_term2)*(expr2))/t3_term9 - (0.000014052694118343981*theta3S*(expr2))/theta3_pow_2 - (0.9758248837815459*(t0_term2)*expr2_term1)/t3_term9))/theta3 - (20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*(expr2))/(t3_term8) - (0.20456897851272804*(t0_term2)*(expr2))/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3_pow_2 + (0.000027144936354308312*flow*theta0_pow_third*theta3S*(1 - 0.63*expr1_pow_3))/(t3_term11*pow(expr1,0.55)) - (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3_pow_2)*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*(expr2))*expr8 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*(expr2))*(expr8*expr8)))/(2.*(expr8*expr8)) - (PI*(expr11*expr11)*((-0.00002280174653761898*flow*theta0_pow_third*theta3S)/(t3_term11*expr1_pow_7) + (20.106192982974676*flow*theta0*((-2276.9823074007572*theta0_pow_2)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) + (2.799181015392208*(t0_term3))/(t3_term2) - (1.4218277665890606*(t0_term2))/t3_term9 - (0.00006885575619957719*(t0_term2)*theta3S)/t3_term14 + (8.558874056109027e-6*theta3S)/theta3_pow_2 - (19.125488600348447*(t0_term3)*(expr2))/(t3_term2) + (1.1886634990796632*(t0_term2)*(expr2))/t3_term9 + (3.6366677468926088e-6*theta3S*(expr2))/theta3_pow_2 + (0.25253171040214556*(t0_term2)*expr2_term1)/t3_term9))/theta3 - (20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*(expr2))/(t3_term8) - (0.7131980994477979*(t0_term2)*(expr2))/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3_pow_2 - (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3_pow_2)*((20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*(expr2))/(t3_term8) - (0.20456897851272804*(t0_term2)*(expr2))/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3)*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*(expr2))*expr8 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*(expr2))*(expr8*expr8)))/(3.*pow((20.106192982974676*flow*theta0*((455.3964614801514*theta0_pow_2)/(t3_term7) - (0.8397543046176623*(t0_term3))/(t3_term8) + (0.8530966599534363*(t0_term2))/(t3_term3) + (5.737646580104534*(t0_term3)*(expr2))/(t3_term8) - (0.7131980994477979*(t0_term2)*(expr2))/(t3_term3) - (0.15151902624128732*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (10.053096491487338*flow*theta0*(1 - 0.63*expr1_pow_3))/theta3,1.6666666666666667)) + (PI*((20.106192982974676*flow*theta0*((6944.541429194567*theta0_pow_2)/(theta3_pow_2*theta3_pow_2*theta3_pow_2) + (6.54468756442183*(t0_term3))/(t3_term2) - (6.22981240275731*(t0_term2))/t3_term9 + (0.0015924284861156552*(t0_term2)*theta3S)/t3_term14 - (0.000057632493196318874*theta3S)/theta3_pow_2 + (442.31556719526833*(t0_term3)*(expr2))/(t3_term2) - (8.004048263161943*(t0_term2)*(expr2))/t3_term9 + (0.00003772148987549525*theta3S*(expr2))/theta3_pow_2 + (2.6193958371136654*(t0_term2)*expr2_term1)/t3_term9))/theta3 - (20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*(expr2))/(t3_term8) + (4.802428957897166*(t0_term2)*(expr2))/(t3_term3) - (1.5716375022681992*(t0_term2)*expr2_term1)/(t3_term3)))/theta3_pow_2 + (20.106192982974676*flow*theta0*((-0.000011601493493937981*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_78) + (0.000010986157772254511*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_74)))/theta3 - (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3_pow_2)*((20.106192982974676*flow*theta0*((-562.0577864939526*theta0_pow_2)/(t3_term7) + (61.66966752706254*(t0_term3))/(t3_term8) - (2.379843507480723*(t0_term2))/(t3_term3) - (0.6802009867403065*(t0_term3)*(expr2))/(t3_term8) - (0.20456897851272804*(t0_term2)*(expr2))/(t3_term3) + (0.5854949302689275*(t0_term2)*expr2_term1)/(t3_term3)))/theta3 + (5.026548245743669*flow*theta0*expr1_pow_45*(1 - 0.63*expr1_pow_3))/theta3)*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*(expr2))*expr8 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*(expr2))*(expr8*expr8)))/(3.*expr11*(expr8*expr8)));
	norm_ring_t3 = norm_merg*norm_ring_t3_pre + (norm_ring/norm_merg)*norm_merg_t3 ;
	
	norm_ring_t3S_pre = ((PI*(expr7*expr7)*expr12*((0.12258447263984441*expr3*(-1.8897 + 1.4547*expr2)*((-0.000034202619806428474*flow*theta0_pow_third)/(theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((-0.00010328363429936579*(t0_term2))/t3_term13 + 0.00001283831108416354/theta3 + (5.455001620338913e-6*expr2)/theta3))/theta3))/(expr6_term1*expr6_term1) + (0.09016171759433259*expr3_pow_2*(1.6557 - 1.8153*expr2)*((-0.000034202619806428474*flow*theta0_pow_third)/(theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((-0.00010328363429936579*(t0_term2))/t3_term13 + 0.00001283831108416354/theta3 + (5.455001620338913e-6*expr2)/theta3))/theta3))/expr6_term1 - (9.63003519301038e-6*(t3_term1)*expr3*expr6_term1)/(t0_term2) + (4.419361834954097e-6*(t3_term1)*expr3_pow_2*(expr6_term1*expr6_term1))/(t0_term2)))/(2.*(expr6_term1*expr6_term1)) + (PI*(expr7*expr7)*((-0.000017101309903214237*flow*theta0_pow_third)/(theta3_pow_third*pow(expr1,0.24999999999999994)) + (20.106192982974676*flow*theta0*((0.000012244328573349254*(t0_term2))/t3_term13 + 3.682455388410919e-6/theta3 - (0.000021079041177515973*expr2)/theta3))/theta3 + (0.00004071740453146247*flow*theta0_pow_third*(1 - 0.63*expr1_pow_3))/(theta3_pow_third*pow(expr1,0.55)))*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(2.*(expr6_term1*expr6_term1)) - (PI*((-0.000034202619806428474*flow*theta0_pow_third)/(theta3_pow_third*expr1_pow_7) + (20.106192982974676*flow*theta0*((-0.00010328363429936579*(t0_term2))/t3_term13 + 0.00001283831108416354/theta3 + (5.455001620338913e-6*expr2)/theta3))/theta3)*(expr7*expr7)*expr12*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(3.*pow(expr6,1.6666666666666667)) + (PI*((20.106192982974676*flow*theta0*((0.0023886427291734827*(t0_term2))/t3_term13 - 0.00008644873979447832/theta3 + (0.00005658223481324288*expr2)/theta3))/theta3 + (20.106192982974676*flow*theta0*((-0.00001740224024090697*(t3_term1))/((t0_term2)*expr1_pow_78) + (0.000016479236658381764*(t3_term1))/((t0_term2)*expr1_pow_74)))/theta3)*expr12*(1 + 0.36775341791953325*expr3*(-1.8897 + 1.4547*expr2)*expr6_term1 + 0.1352425763914989*expr3_pow_2*(1.6557 - 1.8153*expr2)*(expr6_term1*expr6_term1)))/(3.*expr7*(expr6_term1*expr6_term1)));
	norm_ring_t3S = norm_merg*norm_ring_t3S_pre + (norm_ring/norm_merg)*norm_merg_t3S ;
	
	
	normalization_list[0] = norm_merg ;
	normalization_list[1] = norm_ring ;
	normalization_list[2] = norm_merg_t0 ;
	normalization_list[3] = norm_merg_t3 ;
	normalization_list[4] = norm_merg_t3S ;
	normalization_list[5] = norm_ring_t0 ;
	normalization_list[6] = norm_ring_t3 ;
	normalization_list[7] = norm_ring_t3S ;
	
	return normalization_list ;
	
	
}



/*****************************************************************************************************************************

															AMPLITUDES
 ******************************************************************************************************************************/

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/






/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/





/*****************************************************************************************************************************
							Definition and derivatives of inspiral phase of amplitude w.r.t the co-ordinates
 
 
 ******************************************************************************************************************************/


static REAL8 *XLALSimIMRPhenomBAmplitude_Inspiral(
													  const REAL8 f,		/**<Fourier Frequency*/
													  const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
													  const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
													  const REAL8 theta3S,
													  const REAL8 flow	/**< Lower Frequency Cut-off */
){
	static REAL8	amplitude_inspiral_list[4];
	
	REAL8	theta0_pow_2	=	theta0*theta0;
	REAL8	theta3_pow_2	=	theta3*theta3;
	REAL8	theta3_pow_third=	cbrt(theta3);
	REAL8	theta0_sqrt		=	sqrt(theta0);
	
	
	REAL8	flow_term1 = pow(flow,0.8333333333333334);
	REAL8	flow_term2 = pow(flow,1.8333333333333333);
	REAL8	flow_term3 = flow_term1*flow;
	
	
	REAL8	t3_term1	=	cbrt(theta3*theta3);
	REAL8	t3_term3	=	t3_term1*theta3;
	REAL8	t0_term2	=	pow(theta0,2.1666666666666665);
	REAL8	t0_term3	=	pow(theta0,1.1666666666666667);
	
	
	REAL8 insp_coef1, insp_coef2, insp_coef3, insp_amp;
	REAL8 insp_coef1_t0, insp_coef2_t0, insp_coef3_t0, insp_amp_t0;
	REAL8 insp_coef1_t3, insp_coef2_t3, insp_coef3_t3, insp_amp_t3;
	REAL8 insp_coef1_t3S, insp_coef2_t3S, insp_coef3_t3S, insp_amp_t3S;
	
	
	/******************** Inspiral Phase of waveform ******************************/
	
	
	
	insp_coef1	=	(0.016200730203430484/(flow_term1*theta0_sqrt));
	
	insp_coef2	=	(0.03415794470303461/(flow*sqrt(flow)*theta0_sqrt*theta3) - (0.00315938483464183*(t3_term1))/((flow*sqrt(flow))*t0_term3)) ;
	
	insp_coef3	=	(-0.015265371337209142/(flow_term3*pow(theta0,0.8333333333333334)*(t3_term1)) + (0.004839085212385711*theta3)/(flow_term3*theta0*theta0_sqrt) + (1.5442597114335067e-7*theta3S)/(flow_term3*theta0*theta0_sqrt) - (4.8952653483548386e-8*(t3_term3)*theta3S)/(flow_term3*t0_term2));
	
	
	insp_amp		=	insp_coef1*cbrt(sqrt(1./(f*f*f*f*f*f*f))) + insp_coef2*sqrt(1./(f)) + insp_coef3*cbrt(sqrt(1./(f))) ;
	
	
	
	
	/******************** Derivative of inspiral phase w.r.t theta0 co-ordinate ******************************/
	
	insp_coef1_t0	=	-0.008100365101715244/(flow_term1*theta0*theta0_sqrt);
	
	insp_coef2_t0	=	(1/(flow*sqrt(flow)))*(-0.017078972351517306/(theta0*theta0_sqrt*theta3) + (0.003685948973748802*(t3_term1))/(t0_term2)) ;
	
	insp_coef3_t0	=	(1./flow_term2)*(0.01272114278100762/(pow(theta0,1.8333333333333333)*(t3_term1)) - (0.007258627818578568*theta3)/(theta0_pow_2*theta0_sqrt) - (2.3163895671502604e-7*theta3S)/(theta0_pow_2*theta0_sqrt) + (1.0606408254768816e-7*(t3_term3)*theta3S)/(pow(theta0,3.1666666666666665)));
	
	
	insp_amp_t0		=	insp_coef1_t0*cbrt(sqrt(1./(f*f*f*f*f*f*f))) + insp_coef2_t0*sqrt(1./(f)) + insp_coef3_t0*cbrt(sqrt(1./(f))) ;
	
	
	
	
	/******************** Derivative of inspiral phase w.r.t theta3 co-ordinate ******************************/
	
	
	insp_coef1_t3	=	0.;
	
	insp_coef2_t3	=	-0.03415794470303461/(pow(flow,1.5)*sqrt(theta0)*theta3_pow_2) - 0.0021062565564278868/(pow(flow,1.5)*t0_term3*theta3_pow_third);
	
	insp_coef3_t3	=	0.004839085212385711/(flow_term2*(theta0*sqrt(theta0))) + 0.010176914224806093/(flow_term2*pow(theta0,0.8333333333333334)*(t3_term3)) + (2.6469779601696886e-23*theta3S)/(flow_term2*(theta0*sqrt(theta0))*theta3) - (8.158775580591398e-8*(t3_term1)*theta3S)/(flow_term2*t0_term2);
	
	
	insp_amp_t3		=	 insp_coef2_t3*sqrt(1./(f)) + insp_coef3_t3*cbrt(sqrt(1./(f))) ;
	
	
	
	/******************** Derivative of inspiral phase w.r.t theta3S co-ordinate ******************************/
	
	insp_coef1_t3S	=	0.;
	
	insp_coef2_t3S	=	0.;
	
	insp_coef3_t3S	=	1.5442597114335065e-7/(flow_term2*theta0*theta0_sqrt) - (4.895265348354837e-8*(t3_term3))/(flow_term2*t0_term2);
	
	
	insp_amp_t3S		=	 insp_coef2_t3S*sqrt(1./(f)) + insp_coef3_t3S*cbrt(sqrt(1./(f))) ;
	
	amplitude_inspiral_list[0] = insp_amp;
	amplitude_inspiral_list[1] = insp_amp_t0;
	amplitude_inspiral_list[2] = insp_amp_t3;
	amplitude_inspiral_list[3] = insp_amp_t3S;
	
	return amplitude_inspiral_list;
	
}

/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/






/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/





/*****************************************************************************************************************************
						Definition and derivatives of Merger phase of amplitude w.r.t the co-ordinates
 
 
 ******************************************************************************************************************************/


static REAL8 *XLALSimIMRPhenomBAmplitude_Merger(
														  const REAL8 f,		/**<Fourier Frequency*/
														  const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
														  const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
														  const REAL8 theta3S,
														  const REAL8 norm_merg,
														  const REAL8 norm_merg_t0,
														  const REAL8 norm_merg_t3,
													      const REAL8 norm_merg_t3S,
														  const REAL8 flow	/**< Lower Frequency Cut-off */
){
	static REAL8	amplitude_merger_list[4] ;
	
	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_2 = theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);
	REAL8	theta3_pow_6	=	theta3_pow_2*theta3_pow_2*theta3_pow_2;
	
	/* Some combination of above expressions. For optimization purposes, they are calculated here only once and used in the end expression. */
	REAL8	t0_term1	=	theta0_pow_2*theta0 ;
	REAL8	t0_term2	=	theta0_pow_third*theta0_pow_third;
	REAL8	t0_term3	=	theta0_pow_third*theta0;
	REAL8	t0_term4	=	theta0_pow_third*theta0_pow_2;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0;
	REAL8	t0_term5	=	theta0_pow_third*theta0_pow_third*theta0_pow_2;
	REAL8	t0_term7	=	theta0_pow_third*t0_term1;
	
	
	REAL8	t3_term1	=	cbrt(theta3*theta3);
	REAL8	t3_term2	=	t3_term1*t3_term1;
	REAL8	t3_term3	=	t3_term1*theta3;
	REAL8	t3_term4	=	theta3_pow_third*theta3_pow_2;
	REAL8	t3_term5	=	theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term6	=	theta3_pow_2*theta3_pow_2;
	REAL8	t3_term7	=	theta3_pow_2*theta3_pow_2*theta3;
	REAL8	t3_term8	=	theta3_pow_third*theta3_pow_2*theta3;
	REAL8	t3_term9	=	pow(theta3,5.666666666666667);
	REAL8	t3_term10	=	pow(theta3,5.333333333333333);
	
	
	REAL8	expr1		=	-0.7794457395540408 + (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2);
	REAL8	expr1_pow_217 = pow(expr1,0.217);
	REAL8	expr1_pow_26  = pow(expr1,0.26);
	REAL8	expr1_pow_152 = pow(expr1,1.5230000000000001);
	REAL8	expr1_pow_174 = pow(expr1,1.74);
	REAL8	expr1_pow_74 = pow(expr1,0.74);
	REAL8	expr1_pow_78 = pow(expr1,0.783);
	REAL8	expr1_pow_178 = pow(expr1,1.783);
	

	REAL8	flow_term1	=	cbrt(flow);
	
	REAL8 merg_coef1, merg_coef2, merg_coef3, amp_merg;
	REAL8 merg_coef1_t0, merg_coef2_t0, merg_coef3_t0, amp_merg_t0;	/*  */
	REAL8 merg_coef1_t3, merg_coef2_t3, merg_coef3_t3, amp_merg_t3;	/*  */
	REAL8 merg_coef1_t3S, merg_coef2_t3S, merg_coef3_t3S, amp_merg_t3S;	/*  */


	
	/******************** Merger phase of waveform	******************************/
	
	merg_coef1	=	(-0.00002064413376646512*(t0_term2)*(t3_term1) + 4.2844239930916823e-10*t3_term2*theta3S)/(sqrt(flow*flow*flow)*pow(theta0,1.8333333333333333)*sqrt(-((flow*(1.*(t0_term1) + 0.17141979364082702*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615367 + 0.0032075552039126467*pow(expr1,0.217) - 0.002535084595505371*pow(expr1,0.26)))))/(theta3_pow_2*theta3_pow_2*theta3_pow_2))));
	
	merg_coef2	=	merg_coef1*((-4.696087900943841*flow_term1*flow_term1*(t0_term3))/(1.*(t0_term2)*(t3_term1) - 0.00002075371164302091*t3_term2*theta3S));
	
	merg_coef3	=	merg_coef2*((0.25700804969411556*(t0_term2)*theta3_pow_third - 9.630035193010379e-6*theta3*theta3S)/(flow_term1*theta0));

	amp_merg		=	norm_merg*(merg_coef1 + merg_coef2*cbrt(1./(f*f)) + merg_coef3*cbrt(1./(f))) ;
	
	
	
	
	
	/******************** Derivative of merger phase w.r.t theta0 co-ordinate	******************************/
	
				/* Expression which occurs in all three coefficients */
	REAL8	expr_t0		=	(311.7639249918997*(t0_term5) + 53.44250768676573*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-1.6401960827934574*theta3 - 0.0005361712074463486*theta3S) + 1.1431418128389502e-10*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-3.1955728122094408e-6*theta3S + theta3*(-0.22446689113355778 + 1.*expr1_pow_217 - 0.790347923681257*expr1_pow_26)))*sqrt(-((flow*(1.*(t0_term1) + 0.17141979364082702*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615367 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/theta3_pow_6));

	merg_coef1_t0	=	(pow(theta3,6.333333333333333)*theta3S*(2.5454512977543498e-11*(t0_term2) - 5.282756223494714e-16*(t3_term1)*theta3S)*expr1_pow_78 + (t0_term2)*t3_term9*(0.000034406889610775205*(t0_term2) - 9.996989317213926e-10*(t3_term1)*theta3S)*expr1_pow_174 + (t0_term2)*t3_term9*(-0.00002719341376420639*(t0_term2) + 7.901099749923734e-10*(t3_term1)*theta3S)*expr1_pow_178 + expr1_pow_74*(-2.6880222623976484e-11*(t0_term2)*pow(theta3,6.333333333333333)*theta3S + 5.578643892382148e-16*(theta3_pow_2*theta3_pow_2*theta3_pow_2*theta3)*theta3S_pow_2) + (t3_term1)*expr1_pow_152*(0.017162923122909273*t0_term7 + (t0_term5)*(t3_term1)*(0.0025743066475021753*theta3 - 4.4524294680524365e-7*theta3S) + theta0_pow_2*(t3_term4)*(-0.00006772085467284047*theta3 - 9.082868891174881e-8*theta3S) + (t0_term2)*(t3_term5)*(2.2439931127303988e-10*theta3 + 6.3411670597086005e-15*theta3S)*theta3S - 9.79540842086704e-20*t3_term10*(theta3S*theta3S*theta3S) + (t0_term3)*(theta3_pow_2*theta3)*(-7.723207544506218e-6*theta3_pow_2 + 1.7639957324725234e-9*theta3*theta3S + 6.12582609490153e-13*theta3S_pow_2)))/((flow*sqrt(flow))*pow(theta0,2.8333333333333335)*expr1_pow_152*expr_t0);
	
	merg_coef2_t0	=	(-0.06044894671674265*(t0_term5)*expr1_pow_152 - 0.008635121643324476*theta0_pow_2*(t3_term3)*expr1_pow_152 + 0.0002120153908471349*(t0_term3)*(t3_term8)*expr1_pow_152 + 6.930668186581239e-8*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 + 3.098001326801022e-10*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 - 7.388252562766063e-15*(t3_term5)*theta3S_pow_2*expr1_pow_152 + t3_term9*theta3S*(1.262318882391329e-10*expr1_pow_74 - 1.1953663041826004e-10*expr1_pow_78) + (t0_term2)*(t3_term7)*(0.000021761316903740307*expr1_pow_152 - 0.00009694666680616308*expr1_pow_174 + 0.00007662159681806963*expr1_pow_178))/(pow(flow,0.8333333333333334)*(theta0*sqrt(theta0))*expr1_pow_152*expr_t0);
	
	merg_coef3_t0	=	(-0.00003322143168026296*(t0_term2)*expr1_pow_174*(1.*(t0_term2)*t3_term10 - 0.000056204670658011316*theta3_pow_6*theta3S) + 0.000026256489550214566*(t0_term2)*expr1_pow_178*(1.*(t0_term2)*t3_term10 - 0.00005620467065801129*theta3_pow_6*theta3S) + theta3_pow_6*theta3S*(1.*(t0_term2) - 0.000037469780438674206*(t3_term1)*theta3S)*(3.2442611405545104e-11*expr1_pow_74 - 3.0721876250803296e-11*expr1_pow_78) - 0.018125176885355792*expr1_pow_152*(1.*t0_term7*theta3_pow_third + (t0_term5)*theta3*(0.14693125169213742*theta3 - 0.00004817543199258112*theta3S) + theta0_pow_2*pow(theta3,2.6666666666666665)*(-0.0037578710577152088*theta3 - 8.569069996811477e-6*theta3S) + (t0_term2)*(t3_term7)*(2.31238844128011e-8*theta3 + 4.863418983152387e-13*theta3S)*theta3S - 9.813577633519693e-18*t3_term9*(theta3S*theta3S*theta3S) + (t0_term3)*(t3_term8)*(-0.00041142282557802103*theta3_pow_2 + 1.9127211557019266e-7*theta3*theta3S + 6.444048143503295e-11*theta3S_pow_2)))/(pow(flow,1.1666666666666667)*(theta0_pow_2*sqrt(theta0))*expr1_pow_152*expr_t0);
	
	amp_merg_t0 = ( norm_merg*(merg_coef1_t0 + merg_coef2_t0*cbrt(1./(f*f)) + merg_coef3_t0*cbrt(1./(f))) + norm_merg_t0*amp_merg/norm_merg );
	
	
	
	
	
	
/******************** Derivative of merger phase w.r.t theta3 co-ordinate	******************************/
	
	REAL8	expr_t3_1		=	pow(-((flow*(1.*(t0_term1) + 0.17141979364082702*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615367 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/(theta3_pow_2*theta3*theta3_pow_2*theta3)),0.16666666666666666);
	
	REAL8	expr_t3_2	=	pow(-((theta3_pow_2*theta3*theta3_pow_2*theta3)/(flow*(1.*(t0_term1) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))),0.3333333333333333);
	
	
	merg_coef1_t3	=	(0.005984092419836408*sqrt((t0_term2)/(t3_term3))*sqrt(theta3/(flow*theta0))*expr_t3_2*(-0.012649396405173638*pow(theta0,3.3333333333333335)*expr1_pow_152 + theta0_pow_2*(t3_term4)*(0.00003629930231345255*theta3 + 5.778849005804924e-8*theta3S)*expr1_pow_152 + (t0_term5)*(t3_term1)*(-0.001675548530215911*theta3 + 3.102531846241915e-7*theta3S)*expr1_pow_152 + t3_term10*theta3S_pow_2*(-2.9902292599075614e-16*theta3*expr1_pow_74 + 2.831629431307418e-16*theta3*expr1_pow_78 + 5.2504725947499286e-20*theta3S*expr1_pow_152) + (t0_term2)*(t3_term5)*theta3S*(-3.3989520839993703e-15*theta3S*expr1_pow_152 + theta3*(1.4408166169703531e-11*expr1_pow_74 - 1.3643966342086296e-11*expr1_pow_78 - 9.45065733021649e-11*expr1_pow_152 + 4.2102678406114467e-10*expr1_pow_174 - 3.3275764459692264e-10*expr1_pow_178)) + (t0_term3)*(theta3_pow_2*theta3)*(-9.4552578643607e-10*theta3*theta3S*expr1_pow_152 - 3.8991875857900215e-13*theta3S_pow_2*expr1_pow_152 + theta3_pow_2*(2.89782125150883e-6*expr1_pow_152 - 0.000012909793675471837*expr1_pow_174 + 0.000010203228626562591*expr1_pow_178))))/(flow*(t0_term6)*expr1_pow_152*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284163e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26)))*expr_t3_1);
	
	merg_coef2_t3	=	(0.005984092419836408*sqrt((t0_term2)/(t3_term3))*pow(theta3/(flow*theta0),0.8333333333333334)*(15.15240970786892*(t0_term5)*expr1_pow_152 + 1.8759165715941113*theta0_pow_2*(t3_term3)*expr1_pow_152 - 0.03542983229074711*(t0_term3)*(t3_term8)*expr1_pow_152 - 0.000015925002636924063*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 - 5.177061297602283e-8*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 + 1.2346488062709506e-12*(t3_term5)*theta3S_pow_2*expr1_pow_152 + t3_term9*theta3S*(-2.1094575314494196e-8*expr1_pow_74 + 1.9975732664491157e-8*expr1_pow_78) + (t0_term2)*(t3_term7)*(-0.001818263771428786*expr1_pow_152 + 0.008100365101715242*expr1_pow_174 - 0.006402106739200755*expr1_pow_178))*expr_t3_2)/(theta3*expr1_pow_152*(311.7639249918997*(t0_term5) + 53.44250768676573*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-1.6401960827934574*theta3 - 0.0005361712074463486*theta3S) + 1.1431418128389502e-10*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-3.1955728122094408e-6*theta3S + theta3*(-0.22446689113355778 + 1.*expr1_pow_217 - 0.790347923681257*expr1_pow_26)))*expr_t3_1);
	
	
	merg_coef3_t3	=	(0.005984092419836408*sqrt((t0_term2)/(t3_term3))*pow(theta3/(flow*theta0),0.16666666666666666)*expr_t3_2*(0.01387906024401407*pow(theta0,3.3333333333333335)*expr1_pow_152 + (t0_term5)*(t3_term1)*(0.001784359232218124*theta3 - 6.240544080460071e-7*theta3S)*expr1_pow_152 + theta0_pow_2*(t3_term4)*(-0.00003650900315948638*theta3 - 1.0020373539744504e-7*theta3S)*expr1_pow_152 + t3_term10*theta3S_pow_2*(6.51587583988318e-16*theta3*expr1_pow_74 - 6.170277993844701e-16*theta3*expr1_pow_78 - 9.534226462945023e-20*theta3S*expr1_pow_152) + (t0_term3)*(theta3_pow_2*theta3)*(1.858274050499108e-9*theta3*theta3S*expr1_pow_152 + 7.602173253915788e-13*theta3S_pow_2*expr1_pow_152 + theta3_pow_2*(-2.498195954607433e-6*expr1_pow_152 + 0.000011129462977776115*expr1_pow_174 - 8.79614795617277e-6*expr1_pow_178)) + (t0_term2)*(t3_term5)*theta3S*(4.7249779541337574e-15*theta3S*expr1_pow_152 + theta3*(-1.738968246837619e-11*expr1_pow_74 + 1.6467344941995137e-11*expr1_pow_78 + 1.6849233704146434e-10*expr1_pow_152 - 7.506333615197236e-10*expr1_pow_174 + 5.932615187229959e-10*expr1_pow_178))))/(flow*(t0_term6)*expr1_pow_152*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284163e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26)))*expr_t3_1);


	amp_merg_t3 =  ( norm_merg*(merg_coef1_t3 + merg_coef2_t3*cbrt(1./(f*f)) + merg_coef3_t3*cbrt(1./(f))) + norm_merg_t3*amp_merg/norm_merg ) ;
	
	
	/******************** Derivative of merger phase w.r.t theta3S co-ordinate	******************************/

	
	REAL8	expr_t3S_1		=	expr1_pow_152*(1.*(t0_term5) + 0.17141979364082705*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-1.0249976203284163e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26)))*pow(-((flow*(1.*(t0_term1) + 0.17141979364082702*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615367 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/theta3_pow_6),0.16666666666666666);
	
	
	merg_coef1_t3S	=	(0.005984092419836408*pow(theta3/(flow*theta0),1.5)*pow(-(theta3_pow_6/(flow*(1.*(t0_term1) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))),0.3333333333333333)*(7.159688875942881e-8*(t0_term5)*expr1_pow_152 + 9.306612612383264e-9*theta0_pow_2*(t3_term3)*expr1_pow_152 + (t0_term3)*(t3_term4)*(-3.943529918197914e-10*theta3 - 6.156611977563193e-14*theta3S)*expr1_pow_152 + (t3_term5)*theta3S*(-4.485343889861342e-16*theta3*expr1_pow_74 + 4.2474441469611246e-16*theta3*expr1_pow_78 + 3.009265538105056e-36*theta3S*expr1_pow_152) + (t0_term2)*(t3_term6)*(8.980146494695381e-16*theta3S*expr1_pow_152 + theta3*(2.1612249254555295e-11*expr1_pow_74 - 2.0465949513129443e-11*expr1_pow_78 - 5.154903998299904e-11*expr1_pow_152 + 2.2965097312426077e-10*expr1_pow_174 - 1.8150416978013964e-10*expr1_pow_178))))/(sqrt((t0_term2)/(t3_term3))*theta3*expr_t3S_1);
	
	merg_coef2_t3S	=	(0.005984092419836408*sqrt((t0_term2)/(t3_term3))*(t3_term4)*pow(theta3/(flow*theta0),0.8333333333333334)*(1.3930997749196833e-8*(t0_term3)*expr1_pow_152 + 8.302854953049474e-11*(t0_term2)*(t3_term3)*expr1_pow_152 - 5.940306305338388e-15*(t3_term4)*theta3S*expr1_pow_152 + (t3_term8)*(-1.0149302223649968e-10*expr1_pow_74 + 9.610989788993466e-11*expr1_pow_78))*pow(-(theta3_pow_6/(flow*(1.*(t0_term1) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))),0.3333333333333333))/(expr_t3S_1);
	
	
	merg_coef3_t3S	=	(0.005984092419836408*pow(theta3/(flow*theta0),1.1666666666666667)*pow(-(theta3_pow_6/(flow*(1.*(t0_term1) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326626e-6*theta3S) + 3.666690470581712e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284162e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))),0.3333333333333333)*(-1.5601360201150178e-7*(t0_term5)*expr1_pow_152 - 2.3163440900159563e-8*theta0_pow_2*(t3_term3)*expr1_pow_152 + (t0_term3)*(t3_term4)*(8.421296050362544e-10*theta3 + 1.341559985985139e-13*theta3S)*expr1_pow_152 + pow(theta3,5.666666666666667)*theta3S*(9.77381375982477e-16*expr1_pow_74 - 9.25541699076705e-16*expr1_pow_78) + (t0_term2)*(t3_term6)*(-7.271386841174064e-16*theta3S*expr1_pow_152 + theta3*(-2.6084523702564282e-11*expr1_pow_74 + 2.4701017412992696e-11*expr1_pow_78 + 1.1232822469430954e-10*expr1_pow_152 - 5.004222410131491e-10*expr1_pow_174 + 3.955076791486639e-10*expr1_pow_178))))/(sqrt((t0_term2)/(t3_term3))*theta3*expr_t3S_1);
	
	amp_merg_t3S = ( norm_merg*(merg_coef1_t3S + merg_coef2_t3S*cbrt(1./(f*f)) + merg_coef3_t3S*cbrt(1./(f))) + norm_merg_t3S*amp_merg/norm_merg ) ;
	
	amplitude_merger_list[0]	=	amp_merg ;
	amplitude_merger_list[1]	=	amp_merg_t0 ;
	amplitude_merger_list[2]	=	amp_merg_t3 ;
	amplitude_merger_list[3]	=	amp_merg_t3S ;
	
	return amplitude_merger_list ;
	
}


/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/






/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/





/*****************************************************************************************************************************
 Definition and derivatives of Ringdown phase of amplitude w.r.t the co-ordinates
 
 
 ******************************************************************************************************************************/



static REAL8 *XLALSimIMRPhenomBAmplitude_Ringdown(
												const REAL8 f,		/**<Fourier Frequency*/
												const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
												const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
												const REAL8 theta3S,
												const REAL8 norm_ring,
												const REAL8 norm_ring_t0,
												const REAL8 norm_ring_t3,
												const REAL8 norm_ring_t3S,
												const REAL8 flow	/**< Lower Frequency Cut-off */
){
	static REAL8	amp_list_ringdown[4];
	
	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_2 = theta0*theta0;
	REAL8	theta3_pow_6 = theta3_pow_2*theta3_pow_2*theta3_pow_2;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);

	
	/* Some combination of above expressions. For optimization purposes, they are calculated here only once and used in the end expression. */
	REAL8	t0_term1	=	theta0_pow_2*theta0 ;
	REAL8	t0_term2	=	theta0_pow_third*theta0_pow_third;
	REAL8	t0_term3	=	theta0_pow_third*theta0;
	REAL8	t0_term4	=	theta0_pow_third*theta0_pow_2;
	REAL8	t0_term5	=	theta0_pow_third*theta0_pow_third*theta0_pow_2;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0;
	
	REAL8	t3_term1	=	cbrt(theta3*theta3);
	REAL8	t3_term2	=	theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term3	=	t3_term1*theta3;
	REAL8	t3_term4	=	theta3_pow_third*theta3_pow_2;
	REAL8	t3_term5	=	theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term6	=	theta3_pow_2*theta3_pow_2;
	REAL8	t3_term7	=	theta3_pow_2*theta3_pow_2*theta3;
	REAL8	t3_term8	=	theta3_pow_third*theta3_pow_2*theta3;
	REAL8	t3_term11	=	pow(theta3,2.6666666666666665);
	REAL8	t3_term12	=	pow(theta3,3.6666666666666665);
	
	
	
	REAL8	flow_term1	=	pow(flow,0.8333333333333334);
	
	REAL8	expr1		=	-0.7794457395540408 + (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2);
	REAL8	expr1_pow_217 = pow(expr1,0.217);
	REAL8	expr1_pow_26  = pow(expr1,0.26);
	REAL8	expr1_pow_152 = pow(expr1,1.5230000000000001);
	REAL8	expr1_pow_74 = pow(expr1,0.74);
	REAL8	expr1_pow_78 = pow(expr1,0.783);
	REAL8	expr1_pow_178 = pow(expr1,1.783);
	REAL8	expr1_pow_55  = pow(expr1,0.55);
	REAL8	expr1_pow_75  = pow(expr1,0.75);
	REAL8	expr1_pow_7	  = pow(expr1,0.7);
	
	REAL8	expr3	=	1.7794457395540408 - (0.00001800104500292942*(t3_term1)*theta3S)/(t0_term2);
	
	REAL8	expr2	=	pow((20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*(expr3))/(t3_term8) + (4.802428957897166*(t0_term2)*(expr3))/(t3_term3) - (1.5716375022681992*(t0_term2)*(expr3*expr3))/(t3_term3)))/theta3 + (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3,2.1666666666666665);

	REAL8	f_merg, f_merg_t0, f_merg_t3, f_merg_t3S;
	REAL8	sigma, sigma_t0, sigma_t3, sigma_t3S ;
	REAL8	lorentzian ,ring_amp, ring_amp_t0, ring_amp_t3, ring_amp_t3S ;
	REAL8	amp_const, amp_const_t0, amp_const_t3, amp_const_t3S;
	
	
	/*******************************	Define the terms required for the lorentzian and it's derivatives	*************************************/
	f_merg		=	TransitionFrequencies_fring(theta0, theta3, theta3S, flow) ;
	
	f_merg_t0	=	(flow*(27468.867414851156*theta0_pow_2 + 439.5919862087532*(t0_term3)*(t3_term3) + (t0_term2)*(t3_term4)*(-30.017695513874518*theta3 - 0.003461067805343385*theta3S) + ((t3_term5)*theta3S*(-3.2905704843721867e-10*theta3S + (0.00002280174653761898*theta3)/expr1_pow_7))/(t0_term2) + (t3_term6)*(0.0004532979505631769*theta3S + theta3*(10.053096491487338 - 6.333450789637024*pow(expr1,0.30000000000000004)))))/theta3_pow_6;
	
	f_merg_t3	=	(flow*theta0_pow_third*(-54937.73482970231*(t0_term5) - 816.3851172448275*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(48.02831282219921*theta3 + 0.007614349171755447*theta3S) + (t3_term5)*theta3S*(1.3162281937488755e-9*theta3S - (0.00002280174653761898*theta3)/expr1_pow_7) + (t0_term2)*(t3_term6)*(-0.0009065959011263539*theta3S + theta3*(-10.053096491487338 + 6.333450789637023*pow(expr1,0.30000000000000004)))))/(theta3_pow_6*theta3);
	
	f_merg_t3S	=	(flow*theta0_pow_third*(-0.002076640683206031*(t0_term3) + 0.000453297950563177*(t0_term2)*(t3_term3) + (t3_term4)*(-1.9743422906233134e-9*theta3S - (0.000034202619806428474*theta3)/expr1_pow_7)))/t3_term12;
	
	
	sigma	=	(-11300.842322830988*flow*(1.*(t0_term1) - 0.10756774166639459*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(0.0015833529729131948*theta3 - 2.1784821538952197e-8*theta3S) - 3.3754960608409787e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(6.018341074822566e-8*theta3S - 0.00044479412261054025*theta3*pow(expr1,0.45) + 0.0002802202972446404*theta3*expr1_pow_75)))/theta3_pow_6;
	
	sigma_t0	=	(flow*(-33902.52696849296*theta0_pow_2 + 2836.4142043881975*(t0_term3)*(t3_term3) + (t0_term2)*(t3_term4)*(-29.82203714712949*theta3 + 0.00041031138873785176*theta3S) + ((t3_term5)*theta3S*(1.2715316248300336e-9*theta3S - (0.000027144936354308312*theta3)/expr1_pow_55 + (0.000028502183172023725*theta3)/pow(expr1,0.2500000000000001)))/(t0_term2) + (t3_term6)*(-0.00068012323531587*theta3S + 5.026548245743669*theta3*pow(expr1,0.44999999999999996) - 3.1667253948185117*theta3*expr1_pow_75)))/theta3_pow_6;
	
	sigma_t3	=	(flow*theta0_pow_third*(67805.05393698592*(t0_term5) - 5267.626379578083*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(47.715259435407205*theta3 - 0.0009026850552232741*theta3S) + (t3_term5)*theta3S*(-5.086126499320134e-9*theta3S + (0.000027144936354308312*theta3)/expr1_pow_55 - (0.00002850218317202372*theta3)/pow(expr1,0.2500000000000001)) + (t0_term2)*(t3_term6)*(0.00136024647063174*theta3S - 5.026548245743669*theta3*pow(expr1,0.44999999999999996) + 3.1667253948185117*theta3*expr1_pow_75)))/(theta3_pow_6*theta3);
	
	sigma_t3S	=	(flow*theta0_pow_third*(0.00024618683324271107*(t0_term3) - 0.00068012323531587*(t0_term2)*(t3_term3) + (t3_term4)*(7.629189748980202e-9*theta3S + (0.00004071740453146247*theta3)/expr1_pow_55 - (0.000042753274758035594*theta3)/pow(expr1,0.25))))/t3_term12;
	
	
	amp_const	=	0.016200730203430484/(flow_term1*sqrt(theta0)*pow((20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*(expr3))/(t3_term8) + (4.802428957897166*(t0_term2)*(expr3))/(t3_term3) - (1.5716375022681992*(t0_term2)*(expr3*expr3))/(t3_term3)))/theta3 + (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3,1.1666666666666667));
	
	amp_const_t0	=	(-9.976298781101935e-11*theta3_pow_6*(-4.1046191752998394e8*(t0_term5)*expr1_pow_152 - 5.667993385588423e7*theta0_pow_2*(t3_term3)*expr1_pow_152 + 1.3196627548208495e6*(t0_term3)*(t3_term8)*expr1_pow_152 + 431.3906002440869*(t0_term3)*(t3_term4)*theta3S*expr1_pow_152 + 1.753010369598633*(t0_term2)*(t3_term6)*theta3S*expr1_pow_152 - 0.00003344526225653086*(t3_term5)*theta3S_pow_2*expr1_pow_152 + pow(theta3,5.666666666666667)*theta3S*(1.*expr1_pow_74 - 0.9469606458853771*expr1_pow_78) + (t0_term2)*(t3_term7)*(123136.85555380318*expr1_pow_152 - 548574.6914921931*pow(expr1,1.74) + 433564.86840494093*expr1_pow_178)))/(pow(flow,1.8333333333333333)*pow(theta0,1.8333333333333333)*expr1_pow_152*pow(311.7639249918997*(t0_term5) + 53.44250768676574*theta0_pow_2*(t3_term3) + (t0_term3)*(t3_term4)*(-1.6401960827934576*theta3 - 0.0005361712074463486*theta3S) + 1.1431418128389501e-10*(t3_term5)*theta3S_pow_2 + (t0_term2)*(t3_term6)*(-3.195572812209438e-6*theta3S + theta3*(-0.2244668911335578 + 1.*expr1_pow_217 - 0.790347923681257*expr1_pow_26)),2)*pow(-((flow*(1.*(t0_term1) + 0.17141979364082705*(t0_term4)*(t3_term3) + (t0_term6)*(t3_term4)*(-0.0052610194808012925*theta3 - 1.7197987466326628e-6*theta3S) + 3.6666904705817113e-13*theta0_pow_third*(t3_term5)*theta3S_pow_2 + theta0*(t3_term6)*(-1.0249976203284152e-8*theta3S + theta3*(-0.0007199899447615368 + 0.0032075552039126467*expr1_pow_217 - 0.002535084595505371*expr1_pow_26))))/theta3_pow_6),0.16666666666666666));

	amp_const_t3	=	0. - (0.018900851904002234*((20.106192982974676*flow*theta0*((6944.541429194567*theta0_pow_2)/theta3_pow_6 + (6.54468756442183*(t0_term3))/(t3_term2) - (6.22981240275731*(t0_term2))/t3_term11 + (0.0015924284861156552*(t0_term2)*theta3S)/t3_term12 - (0.000057632493196318874*theta3S)/theta3_pow_2 + (442.31556719526833*(t0_term3)*(expr3))/(t3_term2) - (8.004048263161943*(t0_term2)*(expr3))/t3_term11 + (0.00003772148987549525*theta3S*(expr3))/theta3_pow_2 + (2.6193958371136654*(t0_term2)*(expr3*expr3))/t3_term11))/theta3 - (20.106192982974676*flow*theta0*((-1388.9082858389133*theta0_pow_2)/(t3_term7) - (1.9634062693265488*(t0_term3))/(t3_term8) + (3.737887441654386*(t0_term2))/(t3_term3) - (132.6946701585805*(t0_term3)*(expr3))/(t3_term8) + (4.802428957897166*(t0_term2)*(expr3))/(t3_term3) - (1.5716375022681992*(t0_term2)*(expr3*expr3))/(t3_term3)))/theta3_pow_2 + (20.106192982974676*flow*theta0*((-0.000011601493493937981*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_78) + (0.000010986157772254511*theta3S)/((t0_term2)*theta3_pow_third*expr1_pow_74)))/theta3 - (20.106192982974676*flow*theta0*(1 - 4.455*expr1_pow_217 + 3.521*expr1_pow_26))/theta3_pow_2))/(flow_term1*sqrt(theta0)*(expr2));
	
	amp_const_t3S	=	(-0.018900851904002234*((20.106192982974676*flow*theta0*((0.0023886427291734827*(t0_term2))/t3_term11 - 0.00008644873979447832/theta3 + (0.00005658223481324288*(expr3))/theta3))/theta3 + (20.106192982974676*flow*theta0*((-0.00001740224024090697*(t3_term1))/((t0_term2)*expr1_pow_78) + (0.000016479236658381764*(t3_term1))/((t0_term2)*expr1_pow_74)))/theta3))/(flow_term1*sqrt(theta0)*(expr2));


	
	
	/*******************************	Define the ringdown phase of amplitude and it's derivatives		*************************************/


	lorentzian	=	1./((f - f_merg)*(f - f_merg) + sigma*sigma*0.25);
	
	ring_amp	=	norm_ring*amp_const*(1./(2*PI))*(sigma*lorentzian);
	
	ring_amp_t0	=	(1./(2*PI))*(norm_ring*amp_const*( sigma_t0*lorentzian + (2*sigma*(f_merg_t0*(f - f_merg) - sigma*sigma_t0*0.25))*(lorentzian*lorentzian) ) + norm_ring_t0*amp_const*lorentzian*sigma + norm_ring*lorentzian*amp_const_t0*sigma) ;
	
	ring_amp_t3	=	(1./(2*PI))*(norm_ring*amp_const*( sigma_t3*lorentzian + (2*sigma*(f_merg_t3*(f - f_merg) - sigma*sigma_t3*0.25))*(lorentzian*lorentzian) ) + norm_ring_t3*amp_const*sigma*lorentzian + norm_ring*sigma*lorentzian*amp_const_t3) ;
	
	ring_amp_t3S	=	(1./(2*PI))*(norm_ring*amp_const*( sigma_t3S*lorentzian + (2*sigma*(f_merg_t3S*(f - f_merg) - sigma*sigma_t3S*0.25))*(lorentzian*lorentzian) ) + norm_ring_t3S*amp_const*sigma*lorentzian + norm_ring*sigma*lorentzian*amp_const_t3S) ;
	

	amp_list_ringdown[0]	=	ring_amp;
	amp_list_ringdown[1]	=	ring_amp_t0;
	amp_list_ringdown[2]	=	ring_amp_t3;
	amp_list_ringdown[3]	=	ring_amp_t3S;
	
	return amp_list_ringdown;
	
	
	
}


/*****************************************************************************************************************************
 
 ******************************************************************************************************************************/





/*****************************************************************************************************************************
	Define the Phase and it's derivatives.
 
 ******************************************************************************************************************************/



static REAL8 *XLALSimIMRPhenomBPhase(
													  const REAL8 f,		/**<Fourier Frequency*/
													  const REAL8 theta0,	/**< Theta0 component of Chirp-Time Co-ordinate system*/
													  const REAL8 theta3,	/**< Theta3 component of Chirp-Time Co-ordinate system*/
													  const REAL8 theta3S,
													  const REAL8 flow	/**< Lower Frequency Cut-off */
){
	REAL8	theta3_pow_2 = theta3*theta3;
	REAL8	theta0_pow_2 = theta0*theta0;
	REAL8	theta3S_pow_2 = theta3S*theta3S;
	REAL8	theta0_pow_third = cbrt(theta0);
	REAL8	theta3_pow_third = cbrt(theta3);
	
	
	/* Some combination of above expressions. For optimization purposes, they are calculated here only once and used in the end expression. */
	REAL8	t0_term1	=	theta0_pow_2*theta0 ;
	REAL8	t0_term2	=	theta0_pow_third*theta0_pow_third;
	REAL8	t0_term3	=	theta0_pow_third*theta0;
	REAL8	t0_term4	=	theta0_pow_third*theta0_pow_2;
	REAL8	t0_term5	=	theta0_pow_third*theta0_pow_third*theta0_pow_2;
	REAL8	t0_term6	=	theta0_pow_third*theta0_pow_third*theta0;
	
	REAL8	t3_term1	=	cbrt(theta3*theta3);
	REAL8	t3_term2	=	theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term3	=	t3_term1*theta3;
	REAL8	t3_term4	=	theta3_pow_third*theta3_pow_2;
	REAL8	t3_term5	=	theta3_pow_third*theta3_pow_third*theta3_pow_2*theta3_pow_2;
	REAL8	t3_term6	=	theta3_pow_2*theta3_pow_2;
	REAL8	t3_term7	=	theta3_pow_2*theta3_pow_2*theta3;
	REAL8	t3_term8	=	theta3_pow_third*theta3_pow_2*theta3;
	REAL8	t3_term9	=	pow(theta3,1.3333333333333333);
	REAL8	t3_term10	=	pow(theta3,2.6666666666666665);
	REAL8	t3_term11	=	pow(theta3,3.6666666666666665);
	
	REAL8	flow_term1	=	pow(flow,0.3333333333333333);
	REAL8	flow_term2	=	flow_term1*flow_term1;


	static REAL8	Phase_list[3];
	
	REAL8	ph_coef1_t0, ph_coef2_t0, ph_coef3_t0, ph_coef4_t0, ph_coef5_t0, ph_coef6_t0 ;
	REAL8	ph_coef1_t3, ph_coef2_t3, ph_coef3_t3, ph_coef4_t3, ph_coef5_t3, ph_coef6_t3 ;
	REAL8	ph_coef1_t3S, ph_coef2_t3S, ph_coef3_t3S, ph_coef4_t3S, ph_coef5_t3S, ph_coef6_t3S ;
	REAL8	Phase_theta0, Phase_theta3, Phase_theta3S;


	REAL8	freq_coef_1	=	cbrt(1./(f*f*f*f*f));
	REAL8	freq_coef_2	=	1./f;
	REAL8	freq_coef_3	=	cbrt(1./(f*f));
	REAL8	freq_coef_4	=	cbrt(1./(f));
	REAL8	freq_coef_5	=	cbrt(f);
	REAL8	freq_coef_6	=	cbrt(f*f);
	



	/******************************		Derivative of Phase w.r.t theta0 co-ordinate		*************************************/



	ph_coef1_t0		=	0.6000000000000002*pow(flow,1.6666666666666667);
	
	ph_coef2_t0		=	(-496796.31530040567*flow*(t0_term3))/(t3_term2) + (22200.676447943464*flow*(t0_term2))/t3_term10 + (180.1882246452362*flow)/theta3 + (0.13291697653291495*flow*(t3_term1))/(t0_term2) + (0.051892979771241965*flow*theta3S)/theta3_pow_2 - (0.0027503676945816874*flow*theta3S)/((t0_term2)*theta3_pow_third) - (6.872987910440406e-9*flow*theta3_pow_third*theta3S_pow_2)/(t0_term3);
	
	ph_coef3_t0		=	(2.788900203576317e6*flow_term2*theta0)/(t3_term6) - (113312.12594094372*flow_term2*theta0_pow_third)/(t3_term4) - (798.2717109354624*flow_term2)/(theta0_pow_third*(t3_term1)) - (5.551115123125783e-17*flow_term2*theta3)/theta0 - (0.25063441977539375*flow_term2*theta3S)/(theta0_pow_third*(t3_term3)) + (0.000013489164504710593*flow_term2*(t3_term3)*theta3S)/(t0_term6) + (8.169146323050617e-8*flow_term2*(t3_term1)*theta3S_pow_2)/(t0_term6);
	
	ph_coef4_t0		=	(-6.066903031167111e6*flow_term1*(t0_term2))/t3_term11 + (214157.75437429963*flow_term1)/theta3_pow_2 + (1073.0058727005821*flow_term1)/((t0_term2)*theta3_pow_third) + (0.4762758725280963*flow_term1*t3_term9)/(t0_term3) + (0.3680515121677991*flow_term1*theta3S)/((t0_term2)*t3_term9) + (0.0469184914582182*flow_term1*theta3_pow_third*theta3S)/(t0_term3) - (0.000035592264540930435*flow_term1*theta3_pow_2*theta3S)/theta0_pow_2 - (2.763573554553112e-7*flow_term1*theta3*theta3S_pow_2)/theta0_pow_2 + (3.000456553763886e-10*flow_term1*t3_term10*theta3S_pow_2)/(t0_term5);
	
	ph_coef5_t0		=	-3.2789137589417053e6/(flow_term1*(theta3_pow_2*theta3)) + 57783.760772178706/(flow_term1*(t0_term2)*t3_term9) - (1130.7202644890722*theta3_pow_third)/(flow_term1*(t0_term3)) - (0.424928600282445*theta3S)/(flow_term1*(t0_term3)*(t3_term1)) + (0.12588815603440986*theta3*theta3S)/(flow_term1*theta0_pow_2) - (2.3589270756953156e-7*(t3_term3)*theta3S_pow_2)/(flow_term1*(t0_term5));
	
	ph_coef6_t0		=	776165.1005429978/(flow_term2*theta0_pow_third*t3_term10) - 1.4551915228366852e-11/(flow_term2*theta0*theta3) + (886.7368935774405*(t3_term1))/(flow_term2*(t0_term6)) + (0.3345622241910803*theta3S)/(flow_term2*(t0_term6)*theta3_pow_third) - (0.05919190315883502*t3_term9*theta3S)/(flow_term2*(t0_term4)) + (6.179550113157154e-8*theta3_pow_2*theta3S_pow_2)/(flow_term2*(t0_term1));

	Phase_theta0 = ph_coef1_t0*freq_coef_1 + ph_coef2_t0*freq_coef_2 + ph_coef3_t0*freq_coef_3 + ph_coef4_t0*freq_coef_4 + ph_coef5_t0*freq_coef_5 + ph_coef6_t0*freq_coef_6 ;



	/******************************		Derivative of Phase w.r.t theta3 co-ordinate		*************************************/
	
	
	
	ph_coef1_t3		=	0. ;
	
	ph_coef2_t3		=	(922621.728415039*flow*(t0_term4))/pow(theta3,5.333333333333333) - (35521.08231670954*flow*(t0_term6))/t3_term11 - (180.18822464523618*flow*theta0)/theta3_pow_2 + (0.2658339530658298*flow*theta0_pow_third)/theta3_pow_third - (0.10378595954248393*flow*theta0*theta3S)/(theta3_pow_2*theta3) + (0.0027503676945816874*flow*theta0_pow_third*theta3S)/t3_term9 + (6.872987910440398e-9*flow*theta3S_pow_2)/(theta0_pow_third*(t3_term1));

	ph_coef3_t3		=	0.5001535719852379*flow_term2 - (5.577800407152632e6*flow_term2*theta0_pow_2)/(t3_term7) + (198296.2203966515*flow_term2*(t0_term3))/(t3_term8) + (798.2717109354612*flow_term2*(t0_term2))/(t3_term3) + (0.6265860494384846*flow_term2*(t0_term2)*theta3S)/t3_term10 + (1.3877787807814457e-17*flow_term2*theta3S)/theta3 - (0.00003372291126177649*flow_term2*(t3_term1)*theta3S)/(t0_term2) - (8.169146323050617e-8*flow_term2*theta3S_pow_2)/((t0_term2)*theta3_pow_third);


	ph_coef4_t3		=	(1.3347186668567643e7*flow_term1*(t0_term6))/(t3_term5) - (428315.5087485993*flow_term1*theta0)/(theta3_pow_2*theta3) - (1073.0058727005826*flow_term1*theta0_pow_third)/t3_term9 - (1.9051034901123853*flow_term1*theta3_pow_third)/theta0_pow_third - (1.4722060486711963*flow_term1*theta0_pow_third*theta3S)/(t3_term4) - (0.04691849145821825*flow_term1*theta3S)/(theta0_pow_third*(t3_term1)) + (0.00007118452908186087*flow_term1*theta3*theta3S)/theta0 + (2.7635735545531124e-7*flow_term1*theta3S_pow_2)/theta0 - (4.800730486022218e-10*flow_term1*(t3_term3)*theta3S_pow_2)/(t0_term6);

	
	ph_coef5_t3		=	(9.836741276825115e6*theta0)/(flow_term1*(t3_term6)) - (231135.04308871482*theta0_pow_third)/(flow_term1*(t3_term4)) + 1130.7202644890722/(flow_term1*theta0_pow_third*(t3_term1)) - (0.12588815603440984*theta3S)/(flow_term1*theta0) - (0.8498572005648914*theta3S)/(flow_term1*theta0_pow_third*(t3_term3)) + (2.3589270756953154e-7*(t3_term1)*theta3S_pow_2)/(flow_term1*(t0_term6));

	
	ph_coef6_t3		=	(-3.1046604021719913e6*(t0_term2))/(flow_term2*t3_term11) + 57854.32031258002/(flow_term2*theta3_pow_2) - 886.7368935774405/(flow_term2*(t0_term2)*theta3_pow_third) + (0.16728111209553997*theta3S)/(flow_term2*(t0_term2)*t3_term9) + (0.05919190315883502*theta3_pow_third*theta3S)/(flow_term2*(t0_term3)) - (6.179550113157154e-8*theta3*theta3S_pow_2)/(flow_term2*theta0_pow_2);

	
	Phase_theta3 = ph_coef1_t3*freq_coef_1 + ph_coef2_t3*freq_coef_2 + ph_coef3_t3*freq_coef_3 + ph_coef4_t3*freq_coef_4 + ph_coef5_t3*freq_coef_5 + ph_coef6_t3*freq_coef_6 ;

	
	
	/******************************		Derivative of Phase w.r.t theta3S co-ordinate		*************************************/
	
	
	
	ph_coef1_t3S		=	0. ;
	
	ph_coef2_t3S		=	(0.05189297977124195*flow*theta0)/theta3_pow_2 - (0.008251103083745059*flow*theta0_pow_third)/theta3_pow_third + (4.123792746264244e-8*flow*theta3_pow_third*theta3S)/theta0_pow_third;
	
	ph_coef3_t3S		=	0.05406773306195388*flow_term2 - (0.37595162966309076*flow_term2*(t0_term2))/(t3_term3) - (0.00002023374675706589*flow_term2*(t3_term3))/(t0_term2) - (2.450743896915185e-7*flow_term2*(t3_term1)*theta3S)/(t0_term2);
	
	
	ph_coef4_t3S		=	(1.1041545365033965*flow_term1*theta0_pow_third)/t3_term9 - (0.14075547437465458*flow_term1*theta3_pow_third)/theta0_pow_third + (0.00003559226454093043*flow_term1*theta3_pow_2)/theta0 + (5.527147109106224e-7*flow_term1*theta3*theta3S)/theta0 - (3.6005478645166633e-10*flow_term1*t3_term10*theta3S)/(t0_term6);
	
	
	ph_coef5_t3S		=	1.2747858008473352/(flow_term1*theta0_pow_third*(t3_term1)) - (0.12588815603440986*theta3)/(flow_term1*theta0) + (2.8307124908343786e-7*(t3_term3)*theta3S)/(flow_term1*(t0_term6));
	
	
	ph_coef6_t3S		=	-0.5018433362866205/(flow_term2*(t0_term2)*theta3_pow_third) + (0.04439392736912627*t3_term9)/(flow_term2*(t0_term3)) - (6.179550113157156e-8*theta3_pow_2*theta3S)/(flow_term2*theta0_pow_2);
	
	
	Phase_theta3S = ph_coef1_t3S*freq_coef_1 + ph_coef2_t3S*freq_coef_2 + ph_coef3_t3S*freq_coef_3 + ph_coef4_t3S*freq_coef_4 + ph_coef5_t3S*freq_coef_5 + ph_coef6_t3S*freq_coef_6 ;
	

	Phase_list[0]	=	Phase_theta0;
	Phase_list[1]	=	Phase_theta3;
	Phase_list[2]	=	Phase_theta3S;

	return Phase_list ;

}



/**
 *******************************************************************************************************************************
 */

/**
 * Function to compute the metric elements using waveform derivatives
 */
static REAL8 MetricCoeffs(REAL8Vector *Amp, REAL8Vector *dPsii, REAL8Vector *dPsij,
        REAL8Vector *dAi, REAL8Vector*dAj, REAL8Vector *Sh, REAL8 hSqr, REAL8 df) {
	size_t k = Amp->length;
	REAL8 gij   = 0.;
	REAL8 normalise = df*(1./(hSqr));
	for (;k--;) {
				gij += (Amp->data[k]*Amp->data[k]*dPsii->data[k]*dPsij->data[k]
						+ dAi->data[k]*dAj->data[k])/(2.0*Sh->data[k]);
	}
	return gij*normalise;
	
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

int XLALSimIMRPhenomBMetricTheta0Theta3Theta3S(
	REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
	REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
	REAL8 *gamma02,  /**< template metric coeff. 01/10 PN Chirp Time */
	REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
	REAL8 *gamma12,  /**< template metric coeff. 01/10 PN Chirp Time */
	REAL8 *gamma22,  /**< template metric coeff. 01/10 PN Chirp Time */
	const REAL8 Mass,     /**< Total Mass of the system */
	const REAL8 eta,    /**< Symmetric mass ratio */
	const REAL8	chi,	/** Reduced Spin Parameter of the system **/
	const REAL8 flow,   /**< low-frequency cutoff (Hz) */
	const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
) {
	REAL8Vector *Amp=NULL, *dATheta0=NULL, *dATheta3=NULL, *dATheta3S=NULL;
	REAL8Vector *dAT0=NULL, *dAPhi=NULL, *dPhaseTheta0=NULL;
	REAL8Vector *dPhaseTheta3=NULL, *dPhaseTheta3S=NULL, *dPhaseT0=NULL, *dPhasePhi=NULL;
	REAL8 *normalization;
	
	
	/* compute the chirp-time co-ordinates */
	const REAL8 theta0	=	ChirpTime_theta0(Mass,eta, flow);
	const REAL8 theta3	=	ChirpTime_theta3(Mass,eta, flow);
	const REAL8	theta3S	=	ChirpTime_theta3S(Mass,eta, chi, flow);
	
	/* Compute the transition frequencies */
	
	const REAL8 fMerg  = TransitionFrequencies_fmerg(theta0,theta3, theta3S, flow);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
	
	const REAL8 fRing  = TransitionFrequencies_fring(theta0,theta3,theta3S, flow);	/**Frequency at which merger part transitions to ringdown part of the waveform*/
	
	const REAL8 fCut   = TransitionFrequencies_fcut(theta0,theta3, theta3S, flow);	/**Frequency at which ringdown part of the waveform is terminated*/
	
	/*Compute the normalizations and their derivatives*/
	normalization						=		XLALSimIMRPhenomBNormalization(theta0,theta3, theta3S, flow);
	const REAL8	norm_merg				=		normalization[0];
	const REAL8	norm_ringdown			=		normalization[1];
	const REAL8	norm_merg_theta0		=		normalization[2];
	const REAL8	norm_merg_theta3		=		normalization[3];
	const REAL8 norm_merg_theta3S		=		normalization[4];
	const REAL8	norm_ringdown_theta0	=		normalization[5];
	const REAL8	norm_ringdown_theta3	=		normalization[6];
	const REAL8	norm_ringdown_theta3S	=		normalization[7];

	
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
		REAL8 *amplitude_inspiral, *amplitude_merger, *amplitude_ringdown, *phase;

		const REAL8 f = flow + k * df;
		
		phase	=	XLALSimIMRPhenomBPhase(f, theta0,theta3, theta3S, flow);
		
		
		dPhaseTheta0->data[k]	=	phase[0];
		dPhaseTheta3->data[k]	=	phase[1];
		dPhaseTheta3S->data[k]	=	phase[2];
		dPhaseT0->data[k]		=	LAL_TWOPI * f;
		dPhasePhi->data[k]		=	1.;
		

		
		if (f <= fMerg){
			amplitude_inspiral	=	XLALSimIMRPhenomBAmplitude_Inspiral(f, theta0,theta3, theta3S, flow);
			
			/* inspiral amplitude of the waveform */
			Amp->data[k] = amplitude_inspiral[0];
			
			/* inspiral waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_inspiral[1];
			dATheta3->data[k]	=	amplitude_inspiral[2];
			dATheta3S->data[k]	=	amplitude_inspiral[3];
			

		}
		else if ((fMerg<f) && (f<=fRing)){
			amplitude_merger	=	XLALSimIMRPhenomBAmplitude_Merger(f, theta0, theta3, theta3S, norm_merg, norm_merg_theta0, norm_merg_theta3, norm_merg_theta3S, flow);
			
			/* merger amplitude of the frequency-domain waveform */
			Amp->data[k]	=	amplitude_merger[0];
			
			/* merger waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_merger[1];
			dATheta3->data[k]	=	amplitude_merger[2];
			dATheta3S->data[k]	=	amplitude_merger[3];
			


		}
		
		else{
			
			amplitude_ringdown	=	XLALSimIMRPhenomBAmplitude_Ringdown(f, theta0, theta3, theta3S, norm_ringdown, norm_ringdown_theta0, norm_ringdown_theta3, norm_ringdown_theta3S, flow);
			
			/* ringdown amplitude of the frequency-domain waveform */
			Amp->data[k]	=	amplitude_ringdown[0];

			/* ringdown waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_ringdown[1];
			dATheta3->data[k]	=	amplitude_ringdown[2];
			dATheta3S->data[k]	=	amplitude_ringdown[3];
			
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
 *******************************************************************************************************************************
 */







/**
 *******************************************************************************************************************************
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
	REAL8 *normalization;

	
	/* compute the chirp-time co-ordinates */
	const REAL8 theta0	=	ChirpTime_theta0(Mass,eta, flow);
	const REAL8 theta3	=	ChirpTime_theta3(Mass,eta, flow);
	const REAL8	theta3S	=	ChirpTime_theta3S(Mass,eta, 0., flow);
	
	/* Compute the transition frequencies */
	
	const REAL8 fMerg  = TransitionFrequencies_fmerg(theta0,theta3, theta3S, flow);	/**Frequency at which inspiral part transitions to merger part of the waveform*/
	
	const REAL8 fRing  = TransitionFrequencies_fring(theta0,theta3,theta3S, flow);	/**Frequency at which merger part transitions to ringdown part of the waveform*/
	
	const REAL8 fCut   = TransitionFrequencies_fcut(theta0,theta3, theta3S, flow);	/**Frequency at which ringdown part of the waveform is terminated*/
	
	/*Compute the normalizations and their derivatives*/
	normalization						=		XLALSimIMRPhenomBNormalization(theta0,theta3, theta3S, flow);
	const REAL8	norm_merg				=		normalization[0];
	const REAL8	norm_ringdown			=		normalization[1];
	const REAL8	norm_merg_theta0		=		normalization[2];
	const REAL8	norm_merg_theta3		=		normalization[3];
	const REAL8 norm_merg_theta3S		=		normalization[4];
	const REAL8	norm_ringdown_theta0	=		normalization[5];
	const REAL8	norm_ringdown_theta3	=		normalization[6];
	const REAL8	norm_ringdown_theta3S	=		normalization[7];
	
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
		
		REAL8 *amplitude_inspiral, *amplitude_merger, *amplitude_ringdown, *phase;
		
		const REAL8 f = flow + k * df;
		
		phase	=	XLALSimIMRPhenomBPhase(f, theta0,theta3, theta3S, flow);
		
		
		dPhaseTheta0->data[k]	=	phase[0];
		dPhaseTheta3->data[k]	=	phase[1];
		dPhaseT0->data[k]		=	LAL_TWOPI * f;
		dPhasePhi->data[k]		=	1.;
		
		
		if (f <= fMerg){
			amplitude_inspiral	=	XLALSimIMRPhenomBAmplitude_Inspiral(f, theta0,theta3, theta3S, flow);
			
			/* inspiral amplitude of the waveform */
			Amp->data[k] = amplitude_inspiral[0];
			
			/* inspiral waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_inspiral[1];
			dATheta3->data[k]	=	amplitude_inspiral[2];
			
			
		}
		else if ((fMerg<f) && (f<=fRing)){
			amplitude_merger	=	XLALSimIMRPhenomBAmplitude_Merger(f, theta0, theta3, theta3S, norm_merg, norm_merg_theta0, norm_merg_theta3, norm_merg_theta3S, flow);
			
			/* merger amplitude of the frequency-domain waveform */
			Amp->data[k]	=	amplitude_merger[0];
			
			/* merger waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_merger[1];
			dATheta3->data[k]	=	amplitude_merger[2];
			
			
			
		}
		
		else{
			
			amplitude_ringdown	=	XLALSimIMRPhenomBAmplitude_Ringdown(f, theta0, theta3, theta3S, norm_ringdown, norm_ringdown_theta0, norm_ringdown_theta3, norm_ringdown_theta3S, flow);
			
			/* ringdown amplitude of the frequency-domain waveform */
			Amp->data[k]	=	amplitude_ringdown[0];
			
			/* ringdown waveform deratives with respect to the parameters */
			dATheta0->data[k]	=	amplitude_ringdown[1];
			dATheta3->data[k]	=	amplitude_ringdown[2];
			
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





