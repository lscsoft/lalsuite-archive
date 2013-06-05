/*
 *  Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li, Michalis Agathos
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
 #include <stdlib.h> 
 #include <stdio.h>
 #include <string.h>
 #include <lal/LALSimInspiralEOS.h>

 LALEquationOfState XLALSimEOSfromString(char eos_name[])
 {
   LALEquationOfState eos;
   if (!strcmp("MS1",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_MS1;
   else if (!strcmp("H4",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_H4;
   else if (!strcmp("SQM3",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_SQM3;
   else if (!strcmp("MPA1",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_MPA1;
   else if (!strcmp("GNH3",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_GNH3;
   else if (!strcmp("A",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_A;
   else if (!strcmp("AU",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_AU;
   else if (!strcmp("FPS",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_FPS;
   else if (!strcmp("APR",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_APR;
   else if (!strcmp("UU",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_UU;
   else if (!strcmp("L",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_L;
   else
   {
     fprintf(stderr,"Warning! %s not supported. Equation of state set to NONE\n",eos_name);
     eos = LAL_SIM_INSPIRAL_EOS_NONE;
   }
   return eos;
 }

REAL8 XLALSimInspiralEOSLambda(LALEquationOfState eos_type, REAL8 m_intr_msun){/** this must be fed the INTRINSIC mass */

    /* this is fed the intrinsic masses and then computes the value of \Lambda(m) See Hinderer et al ( http://arxiv.org/abs/0911.3535 ) for details of the EOSes*/
    /* \Lambda(m) is in units of s^-5 */
    REAL8 lambda=0.;
//  printf("EOS number: %d\n", eos_type);
//  printf("mass: %e\n", m_intr_msun);
    switch (eos_type)
    {
        case LAL_SIM_INSPIRAL_EOS_NONE:
            lambda = 0.0;
        break;
    // MS1
        case LAL_SIM_INSPIRAL_EOS_MS1:
           // printf("Using EOS MS1\n");
            lambda = 2.755956E-24*(2.19296 + 20.0273*m_intr_msun - 17.9443*m_intr_msun*m_intr_msun 
            + 5.75129*m_intr_msun*m_intr_msun*m_intr_msun - 0.699095*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // H4
        case LAL_SIM_INSPIRAL_EOS_H4:
            lambda = 2.755956E-24*(0.743351 + 15.8917*m_intr_msun - 14.7348*m_intr_msun*m_intr_msun 
            + 5.32863*m_intr_msun*m_intr_msun*m_intr_msun - 0.942625*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break; 
    // SQM3
        case LAL_SIM_INSPIRAL_EOS_SQM3:
            lambda = 2.755956E-24*(-5.55858 + 20.8977*m_intr_msun - 20.5583*m_intr_msun*m_intr_msun 
            + 9.55465*m_intr_msun*m_intr_msun*m_intr_msun - 1.84933*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // MPA1
    case LAL_SIM_INSPIRAL_EOS_MPA1:
        lambda = 2.755956E-24*(0.276761 + 7.26925*m_intr_msun - 5.72102*m_intr_msun*m_intr_msun
        + 1.51347*m_intr_msun*m_intr_msun*m_intr_msun - 0.152181*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // GNH3
    case LAL_SIM_INSPIRAL_EOS_GNH3:
        lambda = 2.755956E-24*(7.80715 + 0.683549*m_intr_msun + 1.21351*m_intr_msun*m_intr_msun
        - 3.50234*m_intr_msun*m_intr_msun*m_intr_msun + 0.894662*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    default:
        lambda = 0.0;
        break;
    }
//  printf("calculated love number: %e\n", lambda);
    if (lambda<0.0) return 0.0;
    else return lambda;
}

REAL8 XLALSimInspiralEOSQfromLambda(REAL8 lambda) {
    /* Quadrupole-monopole parameter calculated from love number;
       see http://arxiv.org/abs/1303.1528 */
    REAL8 q, loglam;
    REAL8 tolerance = 1E-15;
    if(lambda<tolerance) { //printf("Love number is (nearly) zero; cannot compute QM parameter. Setting to 1.0 (BH value).\n");
                      q = 1.0; } 
    else {
    loglam = log(lambda);
    q =  0.194 + 0.0936*loglam + 0.0474*loglam*loglam;
    q -= 0.00421*loglam*loglam*loglam;
    q += 0.000123*loglam*loglam*loglam*loglam;
    q = exp(q);
    }

//  printf("%e %e\n", l, q); // Testing numerical results from these functions.

    return q;

}

REAL8 XLALSimInspiralEOSqmparameter(LALEquationOfState eos_type, REAL8 m_intr_msun){
  
  REAL8 q = 0.0 ;
  REAL8 m = m_intr_msun ;
  REAL8 m2 = m*m ;
  REAL8 m3 = m2*m ;
  
  switch (eos_type) {
  /*  */
  case LAL_SIM_INSPIRAL_EOS_A:
    q = -6.41414141*m3 + 30.70779221*m2 - 53.37417027*m + 35.62253247 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_AU:
    q = -6.18686869*m3 + 30.15909091*m2 - 52.87806638*m + 35.86616883 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_FPS:
    q = -3.86363636*m3 + 21.03030303*m2 - 42.19448052*m + 32.83722944 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_APR:
    q = -10.55555556*m3 + 49.52380952*m2 - 82.77063492*m + 53.02428571 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_UU:
    q = -8.03030303*m3 + 37.61363636*m2 - 63.48733766*m + 41.75080087 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_L:
    q = -6.59090909*m3 + 33.67424242*m2 - 63.77034632*m + 48.98073593 ;
    break;
  case LAL_SIM_INSPIRAL_EOS_NONE:
    q = 1.0 ;
    break;

  default:
    q = 0.0 ;
    break ;
  }
  
  return q ;
}
