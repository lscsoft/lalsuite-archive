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
 
 #include <lal/LALSimInspiralEOS.h>

LALEquationOfState XLALSimEOSfromString(char eos_name[])
{
    LALEquationOfState eos;
    if (!strcmp("MS1",eos_name)) eos = EOS_MS1;
    else if (!strcmp("H4",eos_name)) eos = EOS_H4;
    else if (!strcmp("SQM3",eos_name)) eos = EOS_SQM3;
    else if (!strcmp("MPA1",eos_name)) eos = EOS_MPA1;
    else if (!strcmp("GNH3",eos_name)) eos = EOS_GNH3;
    else 
    {
        fprintf(stderr,"Warning! %s not supported. Equation of state set to NONE\n",eos_name);
        eos = EOS_NONE;
    }
    return eos;
}

REAL8 LambdaOfM_EOS(LALEquationOfState eos_type, REAL8 m_intr_msun){/** this must be fed the INTRINSIC mass */

    /* this is fed the intrinsic masses and then computes the value of \Lambda(m) See Hinderer et al ( http://arxiv.org/abs/0911.3535 ) for details of the EOSes*/
    
    REAL8 lambda=0.;
    switch (eos_type)
    {
        case EOS_NONE:
            lambda = 0.0;
        break;
    // MS1
        case EOS_MS1:
            lambda = 2.755956E-24*(2.19296 + 20.0273*m_intr_msun - 17.9443*m_intr_msun*m_intr_msun 
            + 5.75129*m_intr_msun*m_intr_msun*m_intr_msun - 0.699095*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // H4
        case EOS_H4:
            lambda = 2.755956E-24*(0.743351 + 15.8917*m_intr_msun - 14.7348*m_intr_msun*m_intr_msun 
            + 5.32863*m_intr_msun*m_intr_msun*m_intr_msun - 0.942625*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break; 
    // SQM3
        case EOS_SQM3:
            lambda = 2.755956E-24*(-5.55858 + 20.8977*m_intr_msun - 20.5583*m_intr_msun*m_intr_msun 
            + 9.55465*m_intr_msun*m_intr_msun*m_intr_msun - 1.84933*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // MPA1
    case EOS_MPA1:
        lambda = 2.755956E-24*(0.276761 + 7.26925*m_intr_msun - 5.72102*m_intr_msun*m_intr_msun
        + 1.51347*m_intr_msun*m_intr_msun*m_intr_msun - 0.152181*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // GNH3
    case EOS_GNH3:
        lambda = 2.755956E-24*(7.80715 + 0.683549*m_intr_msun + 1.21351*m_intr_msun*m_intr_msun
        - 3.50234*m_intr_msun*m_intr_msun*m_intr_msun + 0.894662*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    default:
        lambda = 0.0;
        break;
    }
    if (lambda<0.0) return 0.0;
    else return lambda;
}