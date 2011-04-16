/*
*  Copyright (C) 2007 Bernd Machenschalk, B.S. Sathyaprakash, Thomas Cokelaer, Walter Del Pozzo
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

/*  <lalVerbatim file="LALAdvLIGOPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALAdvLIGOPsd.c}}

Module to calculate the noise power spectral density for the initial LIGO detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALAdvLIGOPsdCP}
\idx{LALAdvLIGOPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it
calculates the noise spectral density (per Hz) $S_{h}(f)$
for that frequency. The noise PSD is based on data provided by
Kip Thorne, and the fit by B.S.Sathyaprakash


\begin{equation}
   S_h(f) = S_0\left\{  \left(\frac{f}{f_0}\right)^{-4.14} - 5\left(\frac{f_0}{f}\right)^2 + 111  \left(\frac{1. -
   \frac{f}{f_0}^2 + 0.5  \frac{f}{f_0}^4}{1. + 0.5\frac{f}{f_0}^2} \right)\right\};
\end{equation}
where, $f_0=215$Hz
The returned value is scaled up by $S_0 = 10^{49}.$

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALAdvLIGOPsdCV}}

</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

NRCSID (LALADVLIGOPSDC,"$Id$");

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*  <lalVerbatim file="LALAdvLIGOPsdCP"> */
void
LALAdvLIGOPsd (LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{ /* </lalVerbatim> */

        REAL8 x2,x;
        x = f/215.;
        x2 = x*x;
        // this is the new fit of the latest curve in the ligo DCC 

//    *psd = 0.5*((60000.0*pow(f/10.0,-9.0))+5.0*pow(f/50.0,-5.25)+1.6*pow(f/100.0,-3.25)+ 
                3.2*pow(f/200.0,-1.25)+0.9*pow(f/300.0,-0.08)+0.85*pow(f/1000.0,0.8)+0.35*pow(f/2000.0,0.85));
    
        *psd = pow(x,-4.14) - 5./x2 + 111. * (1. - x2 + 0.5 * x2*x2)/(1. + 0.5*x2);
}


