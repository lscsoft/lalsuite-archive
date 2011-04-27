/*
*  Copyright (C) 2007 Bernd Machenschalk, B.S. Sathyaprakash, Thomas Cokelaer, Salvatore Vitale
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

/*  <lalVerbatim file="LALAdvVIRGOPsdCV">
Author: Vitale, S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALAdvLIGOPsd.c}}

Module to calculate the noise power spectral density for the Advanced Virgo detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALAdvVIRGOPsdCP}
\idx{LALAdvVIRGOPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it
calculates the noise spectral density (per Hz) $S_{h}(f)$
for that frequency. The noise PSD is based on PhysRev D 82 124065 (2010) (NOTE: REPLACE WITH THE ORIGINAL SOURCE!)


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

\vfill{\footnotesize\input{LALAdvVIRGOPsdCV}}

</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

NRCSID (LALADVVIRGOPSDC,"$Id$");

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*  <lalVerbatim file="LALAdvLIGOPsdCP"> */
void
LALAdvVIRGOPsd (LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{ /* </lalVerbatim> */
        /* The prefactor is 10^-47 */
        
        REAL8 x;
        x = f/720.;
        
        *psd = 2.67/1e+07/pow(x,5.6) + 0.68*pow(LAL_E,-0.73*log(x)*log(x))*pow(x,5.34) + 0.59*pow(LAL_E,log(x)*log(x)*(-3.2 -1.08*log(x) -0.13*log(x)*log(x)))/pow(x,4.1)  + 0.68*pow(LAL_E,-0.73*log(x)*log(x))*pow(x,5.34);
        
}


