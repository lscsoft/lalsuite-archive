/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Steven Caudill, Anand Sengupta, Craig Robinson , Thomas Cokelaer, Chris Van Den Broeck
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

/*  <lalVerbatim file="LALInspiralCreateCoarseBankCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralCreateCoarseBank.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateCoarseBankCP}
\idx{LALInspiralCreateCoarseBank()}

\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters
   \item \texttt{nlist,} Output, the number of templates found by the function
   \item \texttt{coarseIn,} Input, specifies the search space, range of masses, etc.
\end{itemize}

The coarse grid algorithm works in two stages:
After computing the minimum and maximum chirp-times corresponding to the
search space: $(\tau_0^{\rm min}, \tau_0^{\rm max}),$
$(\tau_2^{\rm min}, \tau_2^{\rm max})$ (or\footnote{In what follows
we will only mention $\tau_3$; however, the algorithm is itself valid,
and has been implemented, in the case of $(\tau_0,\tau_2)$ too. However,
we recommend that the space $\tau_0$-$\tau_3$ be used.}
$(\tau_3^{\rm min}, \tau_3^{\rm max})$) the algorithm

\begin{enumerate}
\item chooses a lattice of templates along the equal mass curve and then

\item lays a rectangular grid in the rectangular region defined by
the minimum and maximum values of the chirp-times and retain only
if (a) the point lies in the parameter space, OR (b) one of the
vertices defined by the rectangle lies in the parameter space.
\end{enumerate}

\subsubsection*{Description}
\paragraph*{Templates along the equal mass curve}
The algorithm works in two
stages: In the first stage, templates are built along the equal
mass (that is, $\eta=1/4$) curve starting from the minimum value
of the Newtonian chirp-time and stopping at its maximum value.
Given the $n$~th template at $O$ with parameters $(\tau_0^{(n)},\tau_3^{(n)}),$
and given also the distance between templates in our preferred coordinates
$(D\tau_0^{(n)},D\tau_3^{(n)}),$
consider lines $\tau_0 = \tau_0^{(n)} + D\tau_0^{(n)}$
($QA$ of Fig.~\ref{fig:equal mass algo}) and
$\tau_3 = \tau_3^{(n)} + D\tau_3^{(n)}$
($PB$ of Fig.~\ref{fig:equal mass algo}).
\begin{figure}[h]
\begin{center}
\includegraphics[angle=-90,width=4.0 true in]{LALInspiralBankHequalmass}
\end{center}
\caption{Algorithm sketching the placement of templates along $\eta=1/4$ curve.}
\label{fig:equal mass algo}
\end{figure}
The template next to
$(\tau_0^{(n)},\tau_3^{(n)}),$ on the equal mass curve, must lie
either along $PB$ or along $QA$ (cf. Fig.~\ref{fig:equal mass algo}) in order
that all the signals that may lie on $OAB$
are spanned by at least one of the two templates.  Clearly, if we were
to place the $(n+1)$~th template at $B,$ some of the signals won't
have the required minimal match. However, placing the $(n+1)$~th template at
$A$ suffices our requirement.
(Note, however, that there is
no guarantee that this will always work; it works only if the curve
along which templates are being laid is a slowly varying function.)
To locate the $(n+1)$~th template we
compute the following pairs of coordinates:
\begin{eqnarray}
\tau_0^{(n+1)} = \tau_0^{(n)} + D\tau_0^{(n)}, \ \
\tau_3^{(n+1)} =  4A_3 \left ( \frac{\tau_0^{(n+1)}}{4A_0} \right )^{2/5}
\nonumber \\
\tau_3^{(n+1)} = \tau_3^{(n)} + D\tau_3^{(n)}, \ \
\tau_0^{(n+1)} =  4A_0 \left ( \frac{\tau_3^{(n+1)}}{4A_3} \right )^{5/2},
\end{eqnarray}
where
\begin{equation}
A_0=\frac{5}{256 (\pi f_0)^{8/3}}, \ \ A_3=\frac{\pi}{8 (\pi f_0)^{5/3}}.
\end{equation}
Of the two pairs, the required pair is the one that is closer to the
starting point $(\tau_0^{(n)},\tau_3^{(n)}).$

\paragraph*{Templates in the rest of the parameter space}
In the second stage, the algorithm begins again at the point
$(\tau_0^{\rm min}, \tau_3^{\rm min}),$
corresponding distance between templates
$(D\tau_0^{\rm min}, D\tau_3^{\rm min}),$ and chooses a rectangular lattice
of templates in the rectangular region defined by
$(\tau_0^{\rm min}, \tau_3^{\rm min})$
$(\tau_0^{\rm max}, \tau_3^{\rm min})$
$(\tau_0^{\rm max}, \tau_3^{\rm max})$  and
$(\tau_0^{\rm min}, \tau_3^{\rm max})$.
The implementation of the algorithm along the equal mass curve and
in a rectangular lattice in the rest of the parameter space is shown
plotted in Fig.~\ref{fig:coarse}, where the templates
chosen are represented as points.
\begin{figure}[h]
\begin{center}
\includegraphics[angle=-90,width=4.5 true in]{LALInspiralBankHCoarse2}
\caption{Algorithm sketching the construction of a rectangular lattice of templates.}
\label{fig:coarse}
\end{center}
\end{figure}


\subsubsection*{Algorithm}

The algorithm to lay templates along the equal-mass curve is as follows:
\begin{obeylines}
\texttt{
\hskip 1 true cm Begin at $\tau_0 = \tau_0^{\rm min}$
\hskip 1 true cm do while $(\tau_0 < \tau_0^{\rm max})$
\hskip 1 true cm \{
\hskip 2 true cm $\tau_0^A = \tau_0 + D\tau_0, \ \ \tau_3^A = 4A_3 \left ( {\tau_0^A}/{4A_0} \right )^{2/5}$
\hskip 2 true cm $\tau_3^B = \tau_3 + D\tau_3, \ \ \tau_0^B = 4A_0 \left ( {\tau_3^B}/{4A_3} \right )^{5/2}$
\hskip 2 true cm if ($(\tau_0^A,\tau_3^A)$ is closer $(\tau_0,\tau_3)$ than $(\tau_0^B,\tau_3^B)$)
\hskip 2 true cm \{
\hskip 3 true cm $\tau_0 = \tau_0^A, \tau_3 = \tau_3^A$
\hskip 2 true cm \}
\hskip 2 true cm else
\hskip 2 true cm \{
\hskip 3 true cm $\tau_0 = \tau_0^B, \tau_3 = \tau_3^B$
\hskip 2 true cm \}
\hskip 2 true cm Add $(\tau_0, \tau_3)$ to InspiralTemplateList
\hskip 2 true cm numTemplates++
\hskip 2 true cm Compute metric at $(\tau_0, \tau_3)$
\hskip 2 true cm Compute distance between templates at  new point: $(D\tau_0, D\tau_3)$
\hskip 1 true cm \}
}
\end{obeylines}

The algorithm to lay templates in the rest of the parameter space
is as follows:
\begin{obeylines}
\texttt{
\hskip 1 true cm Begin at $\tau_0 = \tau_0^{\rm min}, \tau_3 = \tau_3^{\rm min}$
\hskip 1 true cm Compute metric at $(\tau_0, \tau_3)$
\hskip 1 true cm Compute distance between templates at  new point: $(D\tau_0, D\tau_3)$
\hskip 1 true cm Add $(\tau_0, \tau_3)$ to InspiralTemplateList
\hskip 1 true cm numTemplates++
\hskip 1 true cm do while ($\tau_3 <= \tau_3^{\rm max}$)
\hskip 1 true cm \{
\hskip 2 true cm do while ($\tau_0 <= \tau_0^{\rm max}$)
\hskip 2 true cm \{
\hskip 3 true cm if ($(\tau_0, \tau_3)$ is inside the parameter space)
\hskip 3 true cm \{
\hskip 4 true cm Compute metric at ($\tau_0, \tau_3$)
\hskip 4 true cm Compute distance between templates at  new point: ($D\tau_0, D\tau_3$)
\hskip 4 true cm Add ($\tau_0, \tau_3$) to InspiralTemplateList
\hskip 4 true cm numTemplates++
\hskip 3 true cm \}
\hskip 3 true cm Increment $\tau_0:$ $\tau_0 = \tau_0 + D\tau_0$
\hskip 2 true cm \}
\hskip 2 true cm Increment $\tau_3:$ $\tau_3 = \tau_3 + D\tau_3$
\hskip 2 true cm Get next template along $\tau_3={\rm const.}$: $(\tau_0, \tau_3)$
\hskip 1 true cm \}
}
\end{obeylines}


\subsubsection*{Uses}
\begin{verbatim}
LALInspiralNextTemplate()
LALInspiralParameterCalc()
LALInspiralSetParams()
LALInspiralSetSearchLimits()
LALInspiralComputeParams()
LALInspiralComputeMetric()
LALInspiralUpdateParams()
LALInspiralValidTemplate()
\end{verbatim}

\subsubsection*{Notes}
\clearpage



\input{LALInspiralCreateFlatBankCP}
\idx{LALInspiralCreateFlatBank()}
\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters
   \item \texttt{bankParams,} Input. It is necessary and sufficient to input
   the eigenvalues of the metric and the angle between the $x_0$ axis and the
   semi-major axis of the ambiguity ellipse, that is,
   \texttt{bankParams.metric.g00, bankParams.metric.g11, bankParams.metric.theta,}
   the minimal match, \texttt{bankParams.minimalMatch} and the range of the two
   coordinates over which templates must be chosen:
({\tt bankParams->x0Min}, {\tt bankParams->x0Max}) and
({\tt bankParams->x1Min}, {\tt bankParams->x1Max}).
\end{itemize}
The code expects {\tt list->vectorLength=2} and allocates just the
requisite amount of memory to {\tt list} and returns the number
of grid points in {\tt list->length.} The data points {\tt list->data[2j],}
{\tt j=1,2,\ldots, list->length,} contain the $x_0$-coordinates of the grid
and data points {\tt list->data[2j+1],} contain the $x_1$-coordinates
of the grid.

\subsubsection*{Description}
Given the {\tt metric} and the {\tt minimalMatch} this routine calls
{\tt bank/LALInspiralUpdateParams} to get the spacings in user coordinates (which are
not necessarily the eigen-directions) and lays a uniform grid of templates in
the range specified in ({\tt bankParams->x0Min}, {\tt bankParams->x0Max}) and
({\tt bankParams->x1Min}, {\tt bankParams->x1Max}).
\subsubsection*{Algorithm}
The algorithm to lay templates is as follows: Given the increments $Dx_0$ and
$Dx_1$ found from calling {\tt bank/LALInspiralUpdateParams} lay a rectangular
grid in the space of $(x_0, x_1).$
\begin{obeylines}
\texttt{
\hskip 1 true cm $x_1 = x_1^{\min}$
\hskip 1 true cm do while ($x_1 <= x_1^{\rm max}$)
\hskip 1 true cm \{
\hskip 2 true cm $x_0 = x_0^{\min}$
\hskip 2 true cm do while ($x_0 <= x_0^{\rm max}$)
\hskip 2 true cm \{
\hskip 3 true cm Add ($x_0, x_1$) to list
\hskip 3 true cm numTemplates++
\hskip 3 true cm Increment $x_0:$ $x_0 = x_0 + Dx_0$
\hskip 2 true cm \}
\hskip 2 true cm Increment $x_1:$ $x_1 = x_1 + Dx_1$
\hskip 1 true cm \}
}
\end{obeylines}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralUpdateParams()
LALRalloc()
\end{verbatim}

\subsubsection*{Notes}



\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralCreateCoarseBankCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>

// Prototypes that don't really belong here...
void 
LALInspiralCreatePNCoarseBankUnphysicalEta(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn); 

void
LALInspiralSetSearchLimitsUnphysicalEta (
    LALStatus            *status,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn  coarseIn);

void 
LALInspiralValidTemplateUnphysicalEta(
  LALStatus            *status,
  INT4                 *valid,
  InspiralBankParams   bankParams, 
  InspiralCoarseBankIn coarseIn);

void LALInspiralValidParamsUnphysicalEta(
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   bankParams, 
    InspiralCoarseBankIn coarseIn);

void LALInspiralParameterCalcUnphysicalEta (
   LALStatus        *status, 
   InspiralTemplate *params);

void 
LALInspiralCreatePNCoarseBankHexaUnphysicalEta(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn);

void 
LALInitHexagonalBankUnphysicalEta(
	LALStatus               *status,
	InspiralCell            **cell, 
	INT4                    id,
	InspiralMomentsEtc      *moments, 
	InspiralTemplate        *paramsIn, 
	HexaGridParam           *gridParam, 
	CellEvolution           *cellEvolution,
	CellList **cellList);

void
LALFindPositionUnphysicalEta(LALStatus       *status, 
		REAL4                   dx0, 
		REAL4                   dx1,
		Position                *position, 
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam);

void
LALPopulateCellUnphysicalEta(
		LALStatus               *status,
		InspiralMomentsEtc      *moments,
		InspiralCell            **cell, 
		INT4                     headId,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam,
		CellEvolution           *cellEvolution, 
		CellList		**cellList
			     );

// End prototypes


NRCSID(LALINSPIRALCREATECOARSEBANKC, "$Id$");


/*  <lalVerbatim file="LALInspiralCreateCoarseBankCP"> */
void
LALInspiralCreateCoarseBank(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    )
{  /*  </lalVerbatim>  */
  INT4 i;

  INITSTATUS( status,
      "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.shf.data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( coarseIn.shf.data->data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( coarseIn.mmCoarse > 0.L, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mmCoarse < 1.L, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.fLower > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.tSampling > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.tSampling >= 2.*coarseIn.fUpper, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  switch( coarseIn.approximant )
  {
    case BCV:
      ASSERT( coarseIn.space == Psi0Psi3, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      LALInspiralCreateBCVBank( status->statusPtr, list, nlist, coarseIn );
      CHECKSTATUSPTR( status );
      break;

    case AmpCorPPN:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case TaylorF1:
    case TaylorF2:
    case Eccentricity:
    case PadeT1:
    case PadeF1:
    case EOB:
    case EOBNR:
    case TaylorEt:
    case TaylorN:
    case FindChirpPTF:
      ASSERT( coarseIn.space == Tau0Tau2 || coarseIn.space == Tau0Tau3, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );

      /* Thomas:: we can use either a square placement or an hexagonal
       * placement. The hexagonal placement is along the eigenvalues but
       * not the square one.*/
      if (coarseIn.gridSpacing == Hexagonal){
	LALInspiralCreatePNCoarseBankHexa( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == HexagonalUnphysical){
	LALInspiralCreatePNCoarseBankHexaUnphysicalEta( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == HybridHexagonal){
	LALInspiralCreatePNCoarseBankHybridHexa( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == SquareNotOriented){
	LALInspiralCreatePNCoarseBank( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == SquareNotOrientedUnphysical){
	LALInspiralCreatePNCoarseBankUnphysicalEta( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else {
        ABORT( status, LALINSPIRALBANKH_EGRIDSPACING, LALINSPIRALBANKH_MSGEGRIDSPACING );
      }

      /* Anand:: Nudge the templates only if using max-total-mass cut */
      if ( coarseIn.massRange == MinComponentMassMaxTotalMass  ||
             coarseIn.massRange == MinMaxComponentTotalMass )
      {
          LALNudgeTemplatesToConstantTotalMassLine( status->statusPtr, list, (*nlist), coarseIn);
          CHECKSTATUSPTR( status );
      }

      break;

    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      break;
  }

  /* record the minimal match of the bank in the template and   */
  /* set up the tmplts as a linked list so that it can be       */
  /* manipulated easily by findchirp                            */
  for ( i = 0; i < *nlist - 1 ; ++i )
  {
    (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
    (*list)[i].params.next     = &((*list)[i+1].params);
    (*list)[i].params.fine     = NULL;
  }
  (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
  (*list)[i].params.next = NULL;
  (*list)[i].params.fine = NULL;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/* Anand:: 26 October 2006
 * This function nudges the templates in the list to
 * the (max-total-mass = constant) line.
 * This is done only for those templates whose total
 * mass exceeds the desired max-total-mass value. The
 * templates are nudged along the metric eigen direction
 * until they lie on the said line.
 */
void
LALNudgeTemplatesToConstantTotalMassLine(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 nlist,
    InspiralCoarseBankIn coarseIn
    )
{
  InspiralTemplate      *tempPars=NULL;
  InspiralMetric        *metric=NULL;
  InspiralMomentsEtc    moments;

  INITSTATUS( status, "LALNudgeTemplatesToConstantTotalMassLine",
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* If there are no templates, return now */
  if ( nlist <= 0 )
  {
      LALWarning( status, "number of templates is <= 0 ! " );

      DETATCHSTATUSPTR(status);
      RETURN (status);
  }

  /* Allocate memory (only required to calculate noise moments) */
  tempPars = (InspiralTemplate *)
          LALCalloc( 1, sizeof(InspiralTemplate) );
  metric   = (InspiralMetric *)
          LALCalloc( 1, sizeof(InspiralMetric) );

  /* Init the tempPars */
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  tempPars->totalMass  = coarseIn.MMax;
  tempPars->eta        = 0.25;
  tempPars->ieta       = 1.L;
  tempPars->fLower     = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta;
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD required in the computation of the metric */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );

  /* Loop over template list and nudge the templates if required */
  {
      INT4   i;
      REAL4  P, Q, M, C, ms; /* t0, t3; */

      M = coarseIn.MMax*LAL_MTSUN_SI;
      P = (5./256.)*pow( (LAL_PI*coarseIn.fLower), -8./3. ) ;
      Q = (LAL_PI/8.)*pow( (LAL_PI*coarseIn.fLower), -5./3. ) ;

      for (i=0; i < nlist; i++)
      {
          /* If the totalMass of this template exceeds max-total-mass
           * then nudge along the metric eigen-direction.
           */
          if ( (*list)[i].params.totalMass > coarseIn.MMax )
          {
             ms = tan( LAL_PI/2. + (*list)[i].metric.theta );
             C  = (*list)[i].params.t3 - ms*((*list)[i].params.t0);

             /* Calculate the new co-ordinates in tau0-tau3 space */
             (*list)[i].params.t3  = C / ( 1. -  (ms*P/(M*Q)) );
             (*list)[i].params.t0  = P*(*list)[i].params.t3/(M*Q);

             /* Calculate the other parameters */
             LALInspiralParameterCalc( status->statusPtr, &(*list)[i].params );
             CHECKSTATUSPTR( status );

             /* Check that the new point has not gone down below the
              * equal mass line. If it has, set it to m1=m2=coarseIn.MMax/2.0
              */
             if ( (*list)[i].params.eta > 0.25L )
             {
                 InputMasses originalMassChoice = (*list)[i].params.massChoice;

                 (*list)[i].params.totalMass = coarseIn.MMax ;
                 (*list)[i].params.eta       = 0.25L;
                 (*list)[i].params.massChoice = totalMassAndEta;

                 LALInspiralParameterCalc( status->statusPtr, &(*list)[i].params );
                 CHECKSTATUSPTR( status );

                 /* Reset the massChoice to whatever it was */
                 (*list)[i].params.massChoice = originalMassChoice;
             }

             /* Recalculate the metric at this new point */
             LALInspiralComputeMetric( status->statusPtr, &((*list)[i].metric),
                     &((*list)[i].params), &moments );
             CHECKSTATUSPTR( status );
          }
      }/* Loop over templates */
  }

  /* Clean up */
  LALFree( tempPars );
  LALFree( metric );

  /* Normal exit */
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

void
LALInspiralCreatePNCoarseBank(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    )
{
  InspiralBankParams bankPars, bankParsOld;
  InspiralTemplate *tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;
  INT4 validPars;
  REAL8 x01, x02, x11, x12, dist1, dist2, ndx1, ndx2, a25;

  INITSTATUS( status, "LALInspiralCreateCoarseBank",
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.mMin > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  ndx1 = 0.0;
  ndx2 = 0.0;
  a25 = 0.0;

  /* Number of templates is nlist */
  *nlist = 0;

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */
  if ( !
      (tempPars = (InspiralTemplate *)LALCalloc( 1, sizeof(InspiralTemplate) ))
      )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  metric.space = coarseIn.space;
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimits( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );
  tempPars->totalMass = coarseIn.MMax;
  tempPars->eta = 0.25;
  tempPars->ieta = 1.L;
  tempPars->fLower = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta;
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric                 */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );

  /* compute the metric at this point, update bankPars and add */
  /* the params to the list                                    */
  LALInspiralComputeMetric( status->statusPtr, &metric, tempPars, &moments );
  CHECKSTATUSPTR( status );
  LALInspiralUpdateParams( status->statusPtr, &bankPars, metric,
      coarseIn.mmCoarse );
  CHECKSTATUSPTR( status );

  /* add the first template to the template list */
  *list = (InspiralTemplateList*)
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
  if ( ! *list )
  {
    LALFree( tempPars );
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

  (*list)[*nlist].ID = *nlist;
  (*list)[*nlist].params = *tempPars;
  (*list)[*nlist].metric = metric;
  ++(*nlist);

  /* First lay templates along the equal mass curve; i.e. eta=1/4.      */
  /* Choose the constant and the index converting the chirp times to    */
  /* one another along the curve depending on whether the templates     */
  /* are laid along the tau0-tau2 or tau0-tau3 space                    */
  switch ( coarseIn.space )
  {
    case Tau0Tau2:
      ndx1 = 0.6L;
      ndx2 = 1.L/ndx1;
      a25 = pow(64.L/5.L, ndx1)*(2435.L/8064.L)/pow(LAL_PI*coarseIn.fLower,.4L);
      break;

    case Tau0Tau3:
      a25 = LAL_PI_2 * pow(64.L/5.L, .4L)/pow(LAL_PI * coarseIn.fLower, .6L);
      ndx1 = 0.4L;
      ndx2 = 2.5L;
      break;

    case Psi0Psi3:
    case PTFIntrinsic:
    case PTFFull:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      break;
  }

  bankParsOld = bankPars;

  while ( bankPars.x0 < bankPars.x0Max )
  {
    x01 = bankPars.x0 + bankPars.dx0;
    x11 = a25 * pow(x01,ndx1);
    x12 = bankPars.x1 + bankPars.dx1;
    x02 = pow(x12/a25,ndx2);
    dist1 = pow(bankPars.x0 - x01,2.L) + pow(bankPars.x1 - x11, 2.L);
    dist2 = pow(bankPars.x0 - x02,2.L) + pow(bankPars.x1 - x12, 2.L);
    if ( dist1 < dist2 )
    {
      bankPars.x0 = x01;
      bankPars.x1 = x11;
    }
    else
    {
      bankPars.x0 = x02;
      bankPars.x1 = x12;
    }

    /* If this is a valid point add it to our list */
    LALInspiralValidTemplate( status->statusPtr,
        &validPars, bankPars, coarseIn );
    CHECKSTATUSPTR( status );

    if ( validPars )
    {
      LALInspiralComputeParams( status->statusPtr,
          tempPars, bankPars, coarseIn);
      CHECKSTATUSPTR( status );
      LALInspiralComputeMetric( status->statusPtr,
          &metric, tempPars, &moments );
      CHECKSTATUSPTR( status );
      LALInspiralUpdateParams( status->statusPtr,
          &bankPars, metric, coarseIn.mmCoarse );
      CHECKSTATUSPTR( status );

      *list = (InspiralTemplateList *)
        LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
      if ( ! *list )
      {
        LALFree( tempPars );
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

      (*list)[*nlist].ID = *nlist;
      (*list)[*nlist].params = *tempPars;
      (*list)[*nlist].metric = metric;
      ++(*nlist);
    }
  }

  /* Begin with the parameters found at the first lattice point */
  bankPars = bankParsOld;

  /* Loop along x1 and x0 coordinates until maximum values are reached */
  while ( bankPars.x1 <= bankPars.x1Max )
  {
    /* step along the tau0 axis until the boundary is reached */
    while ( bankPars.x0 <= bankPars.x0Max )
    {
      /* If this is a valid point add it to our list */
      LALInspiralValidTemplate( status->statusPtr,
          &validPars, bankPars, coarseIn );
      CHECKSTATUSPTR( status );

      if ( validPars )
      {
        LALInspiralComputeParams( status->statusPtr,
            tempPars, bankPars, coarseIn );
        CHECKSTATUSPTR( status );
        LALInspiralComputeMetric( status->statusPtr,
            &metric, tempPars, &moments );
        CHECKSTATUSPTR( status );
        LALInspiralUpdateParams( status->statusPtr,
            &bankPars, metric, coarseIn.mmCoarse );
        CHECKSTATUSPTR( status );

        *list = (InspiralTemplateList *)
          LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
        if ( ! *list )
        {
          LALFree( tempPars );
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

        (*list)[*nlist].ID = *nlist;
        (*list)[*nlist].params = *tempPars;
        (*list)[*nlist].metric = metric;
        ++(*nlist);
      }

      bankPars.x0 += bankPars.dx0;
    }
    bankPars = bankParsOld;
    bankPars.x1 += bankPars.dx1;

    /* Find the t0 coordinate of the next template close to the t2/t3 axis */
    LALInspiralNextTemplate( status->statusPtr, &bankPars, metric );
    CHECKSTATUSPTR( status );

    /* Hop along t0-axis until t0 is inside the region of interest or quit */
    LALInspiralValidTemplate( status->statusPtr,
        &validPars, bankPars, coarseIn );
    CHECKSTATUSPTR( status );
    while ( validPars == 0 && bankPars.x0 < bankPars.x0Max )
    {
      bankPars.x0 += bankPars.dx0;
      LALInspiralValidTemplate( status->statusPtr,
          &validPars, bankPars, coarseIn );
      CHECKSTATUSPTR( status );
    }
    bankParsOld = bankPars;
  }
  LALFree( tempPars );

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}


// UnphysicalEta stuff starts here!
void 
LALInspiralCreatePNCoarseBankUnphysicalEta(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
{  
  InspiralBankParams bankPars, bankParsOld;
  InspiralTemplate *tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;
  INT4 validPars;
  REAL8 x01, x02, x11, x12, dist1, dist2, ndx1, ndx2, a25;

  INITSTATUS( status, "LALInspiralCreateCoarseBank", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.mMin > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  ndx1 = 0.0;
  ndx2 = 0.0;
  a25 = 0.0;

  /* Number of templates is nlist */
  *nlist = 0;

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */ 
  if ( !
      (tempPars = (InspiralTemplate *)LALCalloc( 1, sizeof(InspiralTemplate) ))
      )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  metric.space = coarseIn.space;
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimitsUnphysicalEta( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );
  tempPars->totalMass = coarseIn.MMax;
  tempPars->eta = 0.25;
  tempPars->ieta = 1.L;
  tempPars->fLower = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta; 
  LALInspiralParameterCalcUnphysicalEta( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric                 */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );

  /* compute the metric at this point, update bankPars and add */
  /* the params to the list                                    */
  LALInspiralComputeMetric( status->statusPtr, &metric, tempPars, &moments );
  CHECKSTATUSPTR( status );
  LALInspiralUpdateParams( status->statusPtr, &bankPars, metric, 
      coarseIn.mmCoarse );
  CHECKSTATUSPTR( status );

  /* add the first template to the template list */
  *list = (InspiralTemplateList*) 
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
  if ( ! *list )
  {
    LALFree( tempPars );
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

  (*list)[*nlist].ID = *nlist; 
  (*list)[*nlist].params = *tempPars; 
  (*list)[*nlist].metric = metric; 
  ++(*nlist); 

  /* First lay templates along the equal mass curve; i.e. eta=1/4.      */
  /* Choose the constant and the index converting the chirp times to    */
  /* one another along the curve depending on whether the templates     */
  /* are laid along the tau0-tau2 or tau0-tau3 space                    */
  switch ( coarseIn.space ) 
  {
    case Tau0Tau2:
      ndx1 = 0.6L;
      ndx2 = 1.L/ndx1;
      a25 = pow(64.L/5.L, ndx1)*(2435.L/8064.L)/pow(LAL_PI*coarseIn.fLower,.4L);
      break;
      
    case Tau0Tau3:
      a25 = LAL_PI_2 * pow(64.L/5.L, .4L)/pow(LAL_PI * coarseIn.fLower, .6L)
	* pow(4.0, -0.6);
      ndx1 = 0.4L;
      ndx2 = 2.5L;
      break;

    case Psi0Psi3:
    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      break;
  }

  bankParsOld = bankPars;

  while ( bankPars.x0 < bankPars.x0Max ) 
  {
    x01 = bankPars.x0 + bankPars.dx0;
    x11 = a25 * pow(x01,ndx1);
    x12 = bankPars.x1 + bankPars.dx1;
    x02 = pow(x12/a25,ndx2);
    dist1 = pow(bankPars.x0 - x01,2.L) + pow(bankPars.x1 - x11, 2.L);
    dist2 = pow(bankPars.x0 - x02,2.L) + pow(bankPars.x1 - x12, 2.L);
    if ( dist1 < dist2 )
    {
      bankPars.x0 = x01;
      bankPars.x1 = x11;
    } 
    else 
    {
      bankPars.x0 = x02;
      bankPars.x1 = x12;
    }

    /* If this is a valid point add it to our list */
    LALInspiralValidTemplateUnphysicalEta( status->statusPtr, 
        &validPars, bankPars, coarseIn ); 
    CHECKSTATUSPTR( status );

    if ( validPars )
    {
      LALInspiralComputeParams( status->statusPtr, 
          tempPars, bankPars, coarseIn);
      CHECKSTATUSPTR( status );
      LALInspiralComputeMetric( status->statusPtr, 
          &metric, tempPars, &moments );
      CHECKSTATUSPTR( status );
      LALInspiralUpdateParams( status->statusPtr, 
          &bankPars, metric, coarseIn.mmCoarse );
      CHECKSTATUSPTR( status );

      *list = (InspiralTemplateList *) 
        LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
      if ( ! *list )
      {
        LALFree( tempPars );
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

      (*list)[*nlist].ID = *nlist; 
      (*list)[*nlist].params = *tempPars; 
      (*list)[*nlist].metric = metric; 
      ++(*nlist); 
    }
  }

  /* Begin with the parameters found at the first lattice point */
  bankPars = bankParsOld;

  /* Loop along x1 and x0 coordinates until maximum values are reached */
  while ( bankPars.x1 <= bankPars.x1Max ) 
  {
    /* step along the tau0 axis until the boundary is reached */
    while ( bankPars.x0 <= bankPars.x0Max ) 
    {
      /* If this is a valid point add it to our list */
      LALInspiralValidTemplateUnphysicalEta( status->statusPtr, 
          &validPars, bankPars, coarseIn ); 
      CHECKSTATUSPTR( status );
      
      if ( validPars )
      {
        LALInspiralComputeParams( status->statusPtr, 
            tempPars, bankPars, coarseIn );
        CHECKSTATUSPTR( status );
        LALInspiralComputeMetric( status->statusPtr,
            &metric, tempPars, &moments );
        CHECKSTATUSPTR( status );
        LALInspiralUpdateParams( status->statusPtr, 
            &bankPars, metric, coarseIn.mmCoarse );
        CHECKSTATUSPTR( status );

        *list = (InspiralTemplateList *) 
          LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
        if ( ! *list )
        {
          LALFree( tempPars );
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

        (*list)[*nlist].ID = *nlist; 
        (*list)[*nlist].params = *tempPars; 
        (*list)[*nlist].metric = metric; 
        ++(*nlist); 
      }
      
      bankPars.x0 += bankPars.dx0;
    }
    bankPars = bankParsOld;
    bankPars.x1 += bankPars.dx1;
    
    /* Find the t0 coordinate of the next template close to the t2/t3 axis */
    LALInspiralNextTemplate( status->statusPtr, &bankPars, metric );
    CHECKSTATUSPTR( status );

    /* Hop along t0-axis until t0 is inside the region of interest or quit */
    LALInspiralValidTemplateUnphysicalEta( status->statusPtr, 
        &validPars, bankPars, coarseIn );
    CHECKSTATUSPTR( status );
    while ( validPars == 0 && bankPars.x0 < bankPars.x0Max )
    {
      bankPars.x0 += bankPars.dx0;
      LALInspiralValidTemplateUnphysicalEta( status->statusPtr, 
          &validPars, bankPars, coarseIn );
      CHECKSTATUSPTR( status );
    }
    bankParsOld = bankPars;
  } 
  LALFree( tempPars );

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}





void
LALInspiralSetSearchLimitsUnphysicalEta (
    LALStatus            *status,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn  coarseIn
    )
/* </lalVerbatim> */
{  
   InspiralTemplate *Pars1=NULL, *Pars2=NULL, *Pars3=NULL, *Pars4=NULL;

   INITSTATUS( status, "LALInspiralSetSearchLimits", 
	       LALINSPIRALCREATECOARSEBANKC );

   ATTATCHSTATUSPTR( status );

   ASSERT( bankParams, status, 
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
   ASSERT( coarseIn.space == Tau0Tau2 || coarseIn.space == Tau0Tau3, status,
       LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
   ASSERT( coarseIn.mMin > 0, status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.MMax >= 2. * coarseIn.mMin, status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.mmCoarse > 0., status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.mmCoarse < 1., status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.fLower > 0., status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.tSampling > 0., status, 
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

   Pars1 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars2 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars3 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars4 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );

   if ( ! Pars1 || ! Pars2 || ! Pars3 || !Pars4 )
   {
     ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
   }

   /* Initiate three parameter vectors consistent with the coarseIn structure */
   LALInspiralSetParams(status->statusPtr, Pars1, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars2, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars3, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars4, coarseIn);
   CHECKSTATUSPTR(status);

   Pars1->massChoice = Pars2->massChoice = Pars3->massChoice = totalMassAndEta;
   Pars4->massChoice = totalMassAndEta;

   /* Calculate the value of the parameters at the three corners */
   /* of the search space                                        */
   // NO. Pars1->mass1 = Pars1->mass2 = coarseIn.MMax/2.;
   Pars1->totalMass = coarseIn.MMax;
   Pars1->eta       = 1.00;

   LALInspiralParameterCalcUnphysicalEta( status->statusPtr, Pars1 );
   CHECKSTATUSPTR( status );

   if ( coarseIn.massRange == MinMaxComponentTotalMass )
   {
     // NO. Pars2->mass1 = Pars2->mass2 = coarseIn.MMin/2.;
     Pars2->totalMass = coarseIn.MMin;
     Pars2->eta       = 0.25;
   }
   else
   {   
     // NO. Pars2->mass1 = Pars2->mass2 = coarseIn.mMin;
     Pars2->totalMass = 2 * coarseIn.mMin;
     Pars2->eta       = 0.25;
   }
   LALInspiralParameterCalcUnphysicalEta( status->statusPtr, Pars2 );
   CHECKSTATUSPTR( status );
   
   // NO. Pars3->mass1 = coarseIn.mMin;
   // NO. Pars3->mass2 = coarseIn.MMax - coarseIn.mMin;
   Pars3->totalMass = coarseIn.MMax;
   Pars3->eta       = ((coarseIn.MMax - coarseIn.mMin) * coarseIn.mMin) / pow( coarseIn.MMax, 2.0);

   LALInspiralParameterCalcUnphysicalEta( status->statusPtr, Pars3 ); 
   CHECKSTATUSPTR( status );

   if ( coarseIn.massRange == MinMaxComponentTotalMass )
   {
     // NO. Pars4->mass1 = coarseIn.mMin;
     // NO. Pars4->mass2 = coarseIn.MMin - coarseIn.mMin;

     Pars4->totalMass = coarseIn.MMin;
     Pars4->eta       = (coarseIn.MMin * coarseIn.mMin) / pow( coarseIn.MMin + coarseIn.mMin, 2.0);

     LALInspiralParameterCalcUnphysicalEta( status->statusPtr, Pars4 );
     CHECKSTATUSPTR( status );
   }
   else
   {
     Pars4->t0 = 0.0;
   }

   /* Find the minimum and maximum values of the parameters and set     */
   /* the search space.  (The minimum values of chirp times are those   */
   /* corresponding to m1 = m2 = MMax/2, i.e., Pars1 structure.         */
   bankParams->x0 = bankParams->x0Min = Pars1->t0;
   bankParams->x0Max = (Pars2->t0 > Pars4->t0) ? Pars2->t0 : Pars4->t0;

   switch ( coarseIn.space ) 
   {
     case Tau0Tau2:
       bankParams->x1 = bankParams->x1Min = Pars1->t2;
       bankParams->x1Max = (Pars2->t2 > Pars3->t2) ? Pars2->t2 : Pars3->t2;
       break;
       
     case Tau0Tau3:
       bankParams->x1 = bankParams->x1Min = Pars1->t3;
       bankParams->x1Max = (Pars2->t3 > Pars3->t3) ? Pars2->t3 : Pars3->t3;
       break;
   
     default:
       ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
   }
   
   LALFree( Pars1 );
   LALFree( Pars2 );
   LALFree( Pars3 );
   LALFree( Pars4 );

   DETATCHSTATUSPTR( status );
   RETURN( status );
}



void 
LALInspiralValidTemplateUnphysicalEta(
  LALStatus            *status,
  INT4                 *valid,
  InspiralBankParams   bankParams, 
  InspiralCoarseBankIn coarseIn)
{  /*  </lalVerbatim>  */


  INITSTATUS( status, "LALInspiralValidTemplate", LALINSPIRALCREATECOARSEBANKC );

  ATTATCHSTATUSPTR( status );
  
  ASSERT( coarseIn.fLower > 0, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  *valid = 0;
  if ( bankParams.x0 <=0 || bankParams.x1 <=0 )
  {
    LALInfo( status, "x0 or x1 is less than or equal to zero" );
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  
  /* We have a valid template either if the template itself, or one     */
  /* of the vertices of the 'ambiguity rectangle', is in the region of  */
  /* interest                                                           */

  LALInspiralValidParamsUnphysicalEta( status->statusPtr, valid, bankParams, coarseIn ); 
  CHECKSTATUSPTR( status );

  if ( *valid == 1 ) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  bankParams.x1 = bankParams.x1 - bankParams.dx1/2.;
  
  LALInspiralValidParamsUnphysicalEta( status->statusPtr, valid, bankParams, coarseIn ); 
  CHECKSTATUSPTR( status );
  
  if ( *valid == 1 ) 
  {
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

#if 0
  bankParams.x0 = bankParams.x0 - 2.*bankParams.dx0;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  bankParams.x1 = bankParams.x1 + bankParams.dx1;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  bankParams.x1 = bankParams.x1 - 2.*bankParams.dx1;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
#endif

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void LALInspiralValidParamsUnphysicalEta(
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   bankParams, 
    InspiralCoarseBankIn coarseIn
    )
/* </lalVerbatim> */
{

  InspiralTemplate *Pars=NULL;

  INITSTATUS( status, "LALInspiralValidParams", LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( coarseIn.fLower > 0.L, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  *valid = 0;

  if ( bankParams.x0 <=0 || bankParams.x1 <=0 )
  {
    LALInfo( status, "x0 or x1 are less than or equal to zero" );
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  Pars = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
  if ( ! Pars )
  {
    ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }

  /* First set the chirp times of Pars to be as in bankParams */
  Pars->t0 = bankParams.x0;
  Pars->fLower = coarseIn.fLower;
  switch ( coarseIn.space ) 
  {
    case Tau0Tau2:
      Pars->t2 = bankParams.x1;
      Pars->massChoice = t02;
      break;
    case Tau0Tau3:
      Pars->t3 = bankParams.x1;
      Pars->massChoice = t03;
      break;
    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
  }

  /* Compute all the parameters, including masses, */
  /* corresponding to (t0,t2/t3)                   */
  LALInspiralParameterCalcUnphysicalEta( status->statusPtr, Pars );
  CHECKSTATUSPTR( status );

  /* If the masses are in the correct range accept as valid parameters */
  switch (coarseIn.massRange) 
  {
    case MinComponentMassMaxTotalMass:
      if (
          // Pars->mass1 >= coarseIn.mMin &&
          // Pars->mass2 >= coarseIn.mMin &&
          Pars->totalMass <= coarseIn.MMax &&
          Pars->totalMass >= 2*coarseIn.mMin &&
          Pars->eta <= 1.00 && 
          Pars->eta >= coarseIn.etamin
         ) 
      {
        *valid = 1;
      }
      break;
    
    case MinMaxComponentMass:
      if (
          // Pars->mass1 >= coarseIn.mMin &&
          // Pars->mass2 >= coarseIn.mMin &&
          // Pars->mass1 <= coarseIn.mMax &&
          // Pars->mass2 <= coarseIn.mMax &&
          Pars->totalMass >= 2*coarseIn.mMin &&
          Pars->totalMass <= 2*coarseIn.mMax &&
          Pars->eta <= 1.00 && 
          Pars->eta >= coarseIn.etamin
         ) 
      {
        *valid = 1;
      }
      break;

    case MinMaxComponentTotalMass:
      if (
          // Pars->mass1 >= coarseIn.mMin &&
          // Pars->mass2 >= coarseIn.mMin &&
          Pars->totalMass <= coarseIn.MMax &&
          Pars->totalMass >= coarseIn.MMin &&
          Pars->eta <= 1.00 &&
          Pars->eta >= coarseIn.etamin
         )
      {
        *valid = 1;
      }
      break;
      
    default:
      ABORT(status, 999, "Invalid choice for enum InspiralBankMassRange");
  }

  LALFree( Pars );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void 
LALInspiralParameterCalcUnphysicalEta (
   LALStatus        *status, 
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   REAL8 m1, m2, totalMass, eta, mu, piFl, etamin, tiny, ieta;
   REAL8 x1, x2, A0, A2, A3, A4, B2, B4, C4,v,tN;
   REAL8 theta = -11831.L/9240.L;
   REAL8 lambda = -1987.L/3080.L;
   static REAL8 oneby4;
   void *pars;
   DFindRootIn rootIn;
   EtaTau02In Tau2In;
   EtaTau04In Tau4In;

   INITSTATUS (status, "LALInspiralParameterCalc", LALINSPIRALCREATECOARSEBANKC);
   ATTATCHSTATUSPTR(status);
 
   ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT((INT4)params->massChoice >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->massChoice <= 15, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   totalMass 	= 0.0;
   ieta 	= params->ieta;
   ieta 	= 1.;
   oneby4 	= 1./4.;
   etamin 	= 1.e-10;
   tiny 	= 1.e-10;
   piFl 	= LAL_PI * params->fLower;

   switch(params->massChoice) 
   {
      case massesAndSpin:
      /*case spinOnly:*/
      case minmaxTotalMass:
      case m1Andm2:
      case fixedMasses:

         ASSERT(params->mass1 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mass2 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

         m1 = params->mass1;
         m2 = params->mass2;
         params->totalMass = totalMass = m1+m2;
         params->eta = eta = m1*m2/pow(totalMass,2);
         if (params->eta > oneby4) {
      		 params->eta -= tiny;
         }
         params->mu = mu = m1*m2/totalMass;
         params->chirpMass = pow(mu,0.6)*pow(totalMass,0.4);
         params->psi0 = 3./128./params->eta
	                * 1. * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-5./3.) ;
         params->psi3 = -3./128./params->eta
	                * (16 * LAL_PI) * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-2./3.);

      break;

      case totalMassAndEta:
      case totalMassUAndEta:

         ASSERT(params->totalMass > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->eta > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	 
         // NO!
	 // if (params->eta > oneby4) {
	 // params->eta -= tiny;
	 // }
         // ASSERT(params->eta <= oneby4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	 
         totalMass = params->totalMass;
         eta = params->eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->mu = eta*totalMass;
         params->chirpMass = pow(eta,0.6)*totalMass;

      break;

      case totalMassAndMu:

         ASSERT(params->totalMass > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mu > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mu < params->totalMass, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	 
         totalMass = params->totalMass;
         mu = params->mu;
         eta =  (params->mu)/totalMass;
         if (eta > oneby4) {
		 eta -= tiny;
	 }
            params->eta = eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t02:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t2 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	 
         A0 = 5./ pow(piFl, eightby3)/256.;
         A2 = 3715.0/(64512.0*pow(piFl,2.0));
         B2 = 4620.0/3715 * ieta;
         Tau2In.t2 = params->t2;
         Tau2In.A2 = A2 * pow(params->t0/A0, 0.6);
         Tau2In.B2 = B2;
         
	 pars = (void *) &Tau2In;
         rootIn.function = &LALEtaTau02;
         rootIn.xmax = oneby4+tiny;
         rootIn.xmin = etamin;
         rootIn.xacc = 1.e-8;
         LALEtaTau02(status->statusPtr, &x1, rootIn.xmax, pars);
         CHECKSTATUSPTR(status);
         LALEtaTau02(status->statusPtr, &x2, rootIn.xmin, pars);
         CHECKSTATUSPTR(status);
	 
         if (x1*x2 > 0) {
            params->eta = 0.;
            DETATCHSTATUSPTR(status);
            RETURN(status);
         } else {
            LALDBisectionFindRoot(status->statusPtr, &eta, &rootIn, pars);
            CHECKSTATUSPTR(status);
         }
         if (eta > oneby4) {
		 eta-=tiny;
   	}
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t03:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	 
         A0 = 5./ pow(piFl, eightby3)/256.;
         A3 = LAL_PI / pow(piFl, fiveby3)/8.;
         totalMass = A0 * params->t3/(A3 * params->t0);
         eta = A0/(params->t0 * pow(totalMass, fiveby3));
         
	 if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;
 
      case t04:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t4 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

	 A0 = 5./(256. * pow(piFl, eightby3));
         A4 = 5./(128.0 * pow(piFl,fourby3)) * 3058673./1016064.;
         B4 = 5429./1008 * 1016064./3058673. * ieta;
         C4 = 617./144. * 1016064./3058673. * ieta;
         Tau4In.t4 = params->t4;
         Tau4In.A4 = A4 * pow(params->t0/A0, 0.2);
         Tau4In.B4 = B4;
         Tau4In.C4 = C4;

	 pars = (void *) &Tau4In;
         rootIn.function = &LALEtaTau04;
         rootIn.xmax = oneby4+tiny;
         rootIn.xmin = etamin;
         rootIn.xacc = 1.e-8;
         LALEtaTau04(status->statusPtr, &x1, rootIn.xmax, pars);
         CHECKSTATUSPTR(status);
         LALEtaTau04(status->statusPtr, &x2, rootIn.xmin, pars);
         CHECKSTATUSPTR(status);

	 if (x1*x2 > 0) {
            params->eta = 0.;
            DETATCHSTATUSPTR(status);
            RETURN(status);
         } else {
            LALDBisectionFindRoot(status->statusPtr, &eta, &rootIn, pars);
            CHECKSTATUSPTR(status);
         }
         if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

     
      case psi0Andpsi3:
      if (params->psi0 > 0 && params->psi3 < 0)
      {
	      params->totalMass = totalMass = -params->psi3/(16.L * LAL_PI * LAL_PI * params->psi0)/LAL_MTSUN_SI;
	      params->eta = eta = 3.L/(128.L * params->psi0 * pow (LAL_PI * totalMass*LAL_MTSUN_SI, fiveby3));
		      
	      /* if eta < 1/4 amd M > 0 then physical values*/
	      if (eta <= oneby4) 
	      {
		      params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
		      params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
		      params->mu = eta*totalMass;
		      params->chirpMass = pow(eta,0.6)*totalMass;
	      }
      }
      else 
      {
	      params->eta = 0.;
	      DETATCHSTATUSPTR(status);
	      RETURN(status);
      }
      break;

     default:
      ABORT (status, 999, "Improper choice for massChoice in LALInspiralParameterCalc\n");
      break;
   }
   
   if (params->eta > oneby4) {
	   params->eta-=tiny;
	}
   totalMass 	= totalMass*LAL_MTSUN_SI;

   /* Should use the coefficients from LALInspiraSetup.c to avoid errors. 
    * */
   v = pow(piFl * totalMass, 1.L/3.L);
   tN = 5.L/256.L / eta * totalMass / pow(v,8.L);
   
   params->t0 	= 5.0L/(256.0L*eta*pow(totalMass,fiveby3)*pow(piFl,eightby3));
   params->t2 	= (3715.0L + (4620.0L*ieta*eta))/(64512.0*eta*totalMass*pow(piFl,2.0));
   params->t3 	= LAL_PI/(8.0*eta*pow(totalMass,twoby3)*pow(piFl,fiveby3));
   params->t4 	= (5.0/(128.0*eta*pow(totalMass,oneby3)*pow(piFl,fourby3)))
              	* (3058673./1016064. + 5429.*ieta*eta/1008. +617.*ieta*eta*eta/144.);
   params->t5 	= -5.*(7729./252. - 13./3.*ieta*eta)/(256.*eta*params->fLower);
   /* This is a ddraft. t6 and t7 need to be checked propely*/
   params->t6 =  -10052469856691./23471078400. + 128./3.*LAL_PI*LAL_PI
     +(15335597827.L/15240960.L-451.L/12.L*LAL_PI*LAL_PI+352./3.*theta-2464.L/9.L*lambda)*ieta*eta
     +6848.L/105.L* LAL_GAMMA
     -15211.L/1728.L*ieta*eta*eta+25565.L/1296.L*eta*eta*eta*ieta;
   params->t6 = tN * (params->t6  + 6848.L/105.L*log(4.*v)) * pow(v,6);    
   params->t7 = (-15419335.L/127008.L-75703.L/756.L*ieta*eta+14809.L/378.L*ieta*eta*eta) * LAL_PI * tN * pow(v,7);
     
   params->psi0 = 3.L/(128.L * eta * pow(LAL_PI * totalMass, fiveby3));
   params->psi3 = -3.L * LAL_PI/(8.L * eta * pow(LAL_PI * totalMass, twoby3));
   
   switch (params->order) {
                        
      case LAL_PNORDER_NEWTONIAN :
      case LAL_PNORDER_HALF :
         params->t2=0.0;
/*       params->t3=0.0;*/        
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0;
      break;

      case LAL_PNORDER_ONE:
         params->t3=0.0;
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2;
      break;

      case LAL_PNORDER_ONE_POINT_FIVE:
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3;
      break;

      case LAL_PNORDER_TWO:
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4;
      break;

      case LAL_PNORDER_TWO_POINT_FIVE:
         params->t6 = 0.0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5;
      
      case LAL_PNORDER_THREE:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6;
      
      case LAL_PNORDER_THREE_POINT_FIVE:
      default:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0 and t7=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6 - params->t7;
      break;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);
}




void 
LALInspiralCreatePNCoarseBankHexaUnphysicalEta(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
{  
  INT4                  i; 
  INT4 			firstId = 0;
  REAL4                 piFl;
  REAL4 		A0, A3; 
  InspiralBankParams    bankPars;
  InspiralTemplate      *tempPars;
  InspiralMomentsEtc    moments;
  InspiralCell          *cells;  
  HexaGridParam         gridParam;
  CellEvolution         cellEvolution;
  CellList 		*cellList = NULL;

  // INITSTATUS( status, "LALInspiralHexagonalBank", 
  // LALINSPIRALHEXAGONALBANKC );

  INITSTATUS( status, "LALInspiralHexagonalBank", 
	      LALINSPIRALCREATECOARSEBANKC );

  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.mMin > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */ 
  if ( !(tempPars = (InspiralTemplate *) 
  			LALCalloc( 1, sizeof(InspiralTemplate)))) 
  {
    LALFree(tempPars);
    LALFree(cells);
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );
  
  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimits( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );
  
  tempPars->totalMass   = coarseIn.MMax;
  tempPars->eta         = 0.25;
  tempPars->ieta        = 1.L;
  tempPars->fLower      = coarseIn.fLower;
  tempPars->massChoice  = m1Andm2;
  tempPars->mass1       = coarseIn.mMin;
  tempPars->mass2       = coarseIn.mMax;
  
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );
  
  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric  once for all.   */
  LALGetInspiralMoments( 
  		status->statusPtr, 
  		&moments,
   		&coarseIn.shf,
   	 	tempPars );
  CHECKSTATUSPTR( status );
  
  /* Allocate memory for one cell */
  cells = (InspiralCell*)
      LALCalloc(1,   sizeof(InspiralCell) );

  /*define gridParam*/
  gridParam.mm 			= coarseIn.mmCoarse;
  gridParam.x0Min     	= bankPars.x0Min;
  gridParam.x0Max     	= bankPars.x0Max;
  gridParam.x1Min     	= bankPars.x1Min;
  gridParam.x1Max     	= bankPars.x1Max;
  gridParam.mMin      	= coarseIn.mMin;
  gridParam.mMax      	= coarseIn.mMax;
  gridParam.MMin      	= coarseIn.MMin;
  gridParam.MMax      	= coarseIn.MMax;
  gridParam.etaMin    	= coarseIn.etamin;
  gridParam.space     	= coarseIn.space;
  gridParam.massRange 	= coarseIn.massRange;
  gridParam.gridSpacing = coarseIn.gridSpacing;
  

  cellEvolution.nTemplate 		= 1;
  cellEvolution.nTemplateMax 	= 1;
  cellEvolution.fertile 		= 0;

  /* initialise that first cell */
  tempPars->massChoice  = t03;
  cells[0].t0           = tempPars->t0;
  cells[0].t3           = tempPars->t3;

  /* some aliases */
  piFl  = LAL_PI * tempPars->fLower;
  A0    = 5. / pow(piFl, 8./3.) / 256.;
  A3    = LAL_PI / pow(piFl, 5./3.)/8.;


  /* Initialise the first template */
  LALInitHexagonalBankUnphysicalEta(
  			status->statusPtr, 
		       	&cells, firstId, 
		       	&moments, tempPars,
		       	&gridParam, &cellEvolution, 
		       	&cellList);
  CHECKSTATUSPTR( status );

  {
    INT4 k, kk; /*some indexes*/
    INT4 *mylist 		= NULL;
    CellList *ptr 	= NULL;
    INT4 length 	= 1; /* default size of the bank when we 
    						start the bank generation. */

    /* we re-allocate an array which size equals the 
     * template bank size. */
    if (! (mylist =  LALMalloc(length*sizeof(INT4))))
    {
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }

    /* while there are cells/template which can propagate, we carry on the loop.*/
    while (cellEvolution.fertile) 
    {
      length = LALListLength(cellList);
      /*realloc some memory for the next template*/
      if (! (mylist =  LALRealloc(mylist, length*sizeof(INT4))))
      {
		ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
		/* freeing memory here ? */
      }
      ptr = cellList;
      /* we extract the ids which might change within the LALPopulateCell 
       * function. Indeed the bank might grow and then we will lost track 
       * of ids/bank size and so on. */
      for ( k = 0; k < length; k++)
      {
		mylist[k] = ptr->id;	
		ptr = ptr->next;
      }
      /* look at all the template/ids in the current bank to search for fertile cells */
      for (kk = 0; kk < length; kk++)
	  {	
		k = mylist[kk];
		if ( cells[k].status == Fertile) 
		{
		  LALPopulateCellUnphysicalEta(status->statusPtr, 
					       &moments, &cells,
					       k,  tempPars, &gridParam, 
					       &cellEvolution, &cellList);
		  CHECKSTATUSPTR( status );         	 	  
	  	  /* now the bank might have grown, but we only look at the 
	  	   * template created before this for loop, when we entered 
	  	   * in the while loop
	  	   * */
        }
      }
    }
    LALFree(mylist);
  }
  
  if (cellList != NULL)
    ABORT(status, LALINSPIRALBANKH_EHEXAINIT,LALINSPIRALBANKH_MSGEHEXAINIT);
	/* Here is the current number of template generated. Now, we need 
	 * to clean some of them which might be redundant.
	 * */
  *nlist = cellEvolution.nTemplate;

  {
    INT4 k ;
    INT4 length;
    length = cellEvolution.nTemplate;

    for ( k = 0; k < length; k++)
    {  
      REAL4 a;
      REAL4 b;
      REAL4 x0;
      REAL4 tempA3;
      SFindRootIn input;
      INT4 valid;

      PRIN  prin;
	
      tempA3              = pow(A3, -5./2.)/pow(0.25,-1.5);
      tempPars->t0        = cells[k].t0;
      tempPars->t3        = cells[k].t3;

      /* if non physical parameter i.e below eta=0.25*/
      if(cells[k].RectPosition[0] == Below ) 
      {	  
        INT4 above=0, below=0, in=0, out=0;
      		  
	/*first, we define the line which is along the long semi-axis of the 
	 * ambiguity function, defined by the angle theta and the position of 
	 * the template.
	 * */      		  
	a = tan(cells[k].metric.theta);
	b = cells[k].t3 - a * cells[k].t0;
	/* and call a function to search for a solution along eta=1/4 */
	input.function 	= LALSPAF;
	input.xmin 		= cells[k].t3-1e-3;
	input.xmax 		= 1000;
	input.xacc 		= 1e-6;
	
	prin.ct = a * A0 * tempA3;
	prin.b = b;
	
	LALSBisectionFindRoot(status->statusPtr,
			      &x0, &input, (void *)&prin);
	CHECKSTATUSPTR( status );         
	
	tempPars->t3 = x0 + 1e-3; /* to be sure it is physical */
	tempPars->t0 = (tempPars->t3 - b)/a;
	if (tempPars->t0 > 0) 
	  {
	    LALInspiralParameterCalc(status->statusPtr, tempPars);
	    CHECKSTATUSPTR( status );         		  
	  }
	cells[k].t0  = tempPars->t0;
	cells[k].t3  = tempPars->t3;    
	
	/* update its position values */
	valid = 1;
	GetPositionRectangle(status->statusPtr, &cells, k,  tempPars , 
			     &gridParam, 
			     &cellEvolution, 
			     &cellList, 
			     &valid);
	
	{
	  switch (cells[k].RectPosition[1]){
	      case In:    in    +=1; break;
	      case Below: below +=1; break;
	      case Above: above +=1; break;
	      case Out:   out   +=1; break;
              case Edge:             break;
	  }

	  switch (cells[k].RectPosition[2]){
	      case In:    in    +=1; break;
	      case Below: below +=1; break;
	      case Above: above +=1; break;
	      case Out:   out   +=1; break;
              case Edge:             break;
	  }

	  switch (cells[k].RectPosition[3]){
	      case In:    in    +=1; break;
	      case Below: below +=1; break;
	      case Above: above +=1; break;
	      case Out:   out   +=1; break;
              case Edge:             break;
	  }

	  switch (cells[k].RectPosition[4]){
	      case In:    in    +=1; break;
	      case Below: below +=1; break;
	      case Above: above +=1; break;
	      case Out:   out   +=1; break;
              case Edge:             break;
	  }
	}

	if (above == 2 && cells[k].position == In)
	  {
	    cells[cells[k].child[0]].position = Out;	    
	  }
      }  
    }
  }

  for (i=0; i<cellEvolution.nTemplate; i++) {
    if (cells[i].position == In ) {
      *nlist = *nlist +1; 
    }
  }


  /* allocate appropriate memory and fill the output bank */
  *list = (InspiralTemplateList*) 
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist+1) );
  if ( ! *list )
  {
    LALFree( tempPars );
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );
  {
    *nlist = 0 ;
    for (i=0; i<cellEvolution.nTemplate; i++) 
    {
      if (cells[i].position == In) 
      {
        tempPars->t0  = cells[i].t0;
        tempPars->t3  = cells[i].t3;
        tempPars->massChoice = t03;
        tempPars->fLower = coarseIn.fLower;

        LALInspiralParameterCalc( status->statusPtr, tempPars );
        CHECKSTATUSPTR( status );
	    
        (*list)[*nlist].ID            = *nlist; 
        (*list)[*nlist].params        = *tempPars; 
        (*list)[*nlist].metric        = cells[i].metric; 
        ++(*nlist); 
      }
    }
  }
   
  LALFree( cells );
  LALFree( tempPars );
    
  DETATCHSTATUSPTR( status );
  RETURN ( status );
}



void 
LALInitHexagonalBankUnphysicalEta(
	LALStatus               *status,
	InspiralCell            **cell, 
	INT4                    id,
	InspiralMomentsEtc      *moments, 
	InspiralTemplate        *paramsIn, 
	HexaGridParam           *gridParam, 
	CellEvolution           *cellEvolution,
	CellList **cellList)
{
  INT4          i;
  INT4 		valid;   
  
  //   INITSTATUS( status, "LALInitHexagonalBank", 
  // LALINSPIRALHEXAGONALBANKC );

  INITSTATUS( status, "LALInitHexagonalBank", 
	      LALINSPIRALCREATECOARSEBANKC );

  ATTATCHSTATUSPTR( status );
  
  /* a new cell is created; by default it can create new children, 
     therefore it is fertile */
  cellEvolution->fertile = cellEvolution->fertile + 1;;
  (*cell)[id].status = Fertile;  
  LALListAppend(cellList, id);


  /* all of whom are unset and do not have any id set yet*/
  for (i = 0; i < 6; i++) 
  {
    (*cell)[id].child[i] = -1;
  } 
 
  /* filled some values related to the space */
  (*cell)[id].ID        = id;  
  (*cell)[id].position  = In;
  (*cell)[id].metric.space = gridParam->space;


  /* before any further computation, check that t0, t3 are positive.*/
  if ((*cell)[id].t0 > 0 && (*cell)[id].t3 > 0)
  {
    /* Get the metric at the position of the cell */ 
    paramsIn->t0 = (*cell)[id].t0;
    paramsIn->t3 = (*cell)[id].t3;

    LALInspiralComputeMetric( status->statusPtr, 
			      &((*cell)[id].metric),
			      paramsIn,
			      moments);
    CHECKSTATUSPTR( status );
  
    /* let us store the dx0 and dx3 at that point. */
    (*cell)[id].dx0 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g00 );
    (*cell)[id].dx1 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g11 );

    LALFindPositionUnphysicalEta(status->statusPtr, 
				 (*cell)[id].dx0, (*cell)[id].dx1,
				 &((*cell)[id].RectPosition[0]), 
				 paramsIn, gridParam);
    CHECKSTATUSPTR( status );

    /* if outside, this is a sterile cell which can not propagate */  
    if ((*cell)[id].RectPosition[0] == Out) 
    {
      (*cell)[id].position      = Out;
      for (i = 0; i < 5; i++)
      {
        (*cell)[id].RectPosition[i] = Out;
      }
      (*cell)[id].status = Sterile;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      LALListDelete(cellList, id);      
      
      DETATCHSTATUSPTR(status);
      RETURN(status);
    }
    else
    {
      valid = 1;
      GetPositionRectangle(status->statusPtr, &(*cell), id,  paramsIn , 
			   gridParam, cellEvolution, &(*cellList), &valid);
    }
  }
  else
  {/* if t0 or t3 < 0 , this is not a valid cell*/
    valid = 0;   
  }

  /* If this is not a valid template, we remove it from the bank*/
  if (valid == 0)
  {
    for (i=0; i<5; i++)
    {
      (*cell)[id].RectPosition[i] = Out;
    }
    (*cell)[id].position 		= Out;
    (*cell)[id].status 			= Sterile;
    (cellEvolution->fertile)	=cellEvolution->fertile-1;
    LALListDelete(cellList, id);
  }




#if 1
  if (gridParam->gridSpacing == HybridHexagonal)
  {
    INT4 below=0, above=0;    
    for (i=1; i<=4; i++){
      if ( (*cell)[id].RectPosition[i] == Below) below++;
      if ( (*cell)[id].RectPosition[i] == Above) above++;
    }    
    if (below==2 && above == 2){
      (*cell)[id].status = Edge;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      LALListDelete(cellList, id);
      
    } 
  

  }
#endif  



  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALFindPositionUnphysicalEta(LALStatus       *status, 
		REAL4                   dx0, 
		REAL4                   dx1,
		Position                *position, 
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam
)
{
  REAL8 	mint3;  
  REAL4   	eta;
  REAL4 	totalMass,ieta, oneby4, tiny, piFl, A0, A3;

  INITSTATUS( status, "LALFindPosition", 
	      // LALINSPIRALHEXAGONALBANKC );
	      LALINSPIRALCREATECOARSEBANKC );

  ATTATCHSTATUSPTR( status );

  ieta 	  = 1.;
  oneby4  = 1./4.;
  tiny 	  = 1.e-10;
  piFl 	  = LAL_PI * paramsIn->fLower;
  A0      = 5. / pow(piFl, 8./3.) / 256.;
  A3      = LAL_PI / pow(piFl, 5./3.)/8.;
  
  ASSERT(paramsIn->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(paramsIn->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  /* given t0, t3 we get the totalMass and eta. 
     We do not need to call ParameterCalc again and again here. */
  totalMass     = A0 * paramsIn->t3/(A3 * paramsIn->t0);
  eta           = A0/(paramsIn->t0 * pow(totalMass, fiveby3));
  
  /* be sure eta is inside the space if it is suppose to be */
  // NO
  // if (eta > oneby4) {
  // eta-=tiny;
  // }
  
  /* let us fill the param strucutre now : eta, Mass, mass1, mass2*/
  paramsIn->eta = eta;
  totalMass     = paramsIn->totalMass = totalMass/LAL_MTSUN_SI;
  if (eta <= oneby4) {
    paramsIn->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
    paramsIn->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
  }
  
  /* does t3 positive*/  
  if ((paramsIn->t3-dx1)<0)
  { 
    mint3 = 0;
  }
  else
  {
    mint3 = paramsIn->t3-dx1;
  }
  
  if (    (paramsIn->t0 < gridParam->x0Min - dx0)
       || (paramsIn->t0 > gridParam->x0Max + dx0) 
       || (paramsIn->t3 <= mint3))
  {
    *position = Out;
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }   

  float maxEta = 1.00; // 0.25;

  switch ( gridParam->massRange )
  {
    case MinMaxComponentMass:
      if (
	  // paramsIn->mass1 >= gridParam->mMin &&
	  // paramsIn->mass2 >= gridParam->mMin &&
	  // paramsIn->mass1 <= gridParam->mMax &&
	  // paramsIn->mass2 <= gridParam->mMax &&
	  paramsIn->totalMass <= (gridParam->mMax * 2) &&
	  paramsIn->totalMass >= (gridParam->mMax * 2) &&
          paramsIn->eta <= maxEta && 
          paramsIn->eta >= gridParam->etaMin
          ) 
        {
          *position = In;
        }
      else
        if (paramsIn->eta > maxEta){
          *position = Below; 
        }
        else{
          *position = Above;
        }
      break;

    case MinComponentMassMaxTotalMass:
      if (
	  // paramsIn->mass1 >= gridParam->mMin &&
	  // paramsIn->mass2 >= gridParam->mMin &&
	  paramsIn->totalMass <= (gridParam->mMax * 2) &&
	  paramsIn->totalMass >= (gridParam->mMax * 2) &&

          paramsIn->totalMass <= gridParam->MMax &&
          paramsIn->eta <= maxEta &&
          paramsIn->eta >= gridParam->etaMin
          )
        {
          *position = In;
        }
      else
        if (paramsIn->eta > maxEta){
          *position = Below;
        }
        else{
          *position = Above;
        }
      break;

    case MinMaxComponentTotalMass:
      if (
          // paramsIn->mass1 >= gridParam->mMin &&
          // paramsIn->mass2 >= gridParam->mMin &&
          paramsIn->totalMass <= gridParam->MMax &&
          paramsIn->totalMass >= gridParam->MMin &&
          paramsIn->eta <= maxEta &&
          paramsIn->eta >= gridParam->etaMin
          )
        {
          *position = In;
        }
      else if (paramsIn->eta > maxEta ){
          *position = Below;
        }
      else{
        *position = Above;
        }

      /* Now cut out unnecessary templates */
      if ( paramsIn->totalMass < gridParam->MMin )
      {
        REAL4 totalMass2 = A0 * (paramsIn->t3 - dx1)/(A3 * paramsIn->t0);
        totalMass2 = totalMass2 / LAL_MTSUN_SI;
        totalMass     = A0 * paramsIn->t3/(A3 * (paramsIn->t0 - dx0));
        totalMass = totalMass / LAL_MTSUN_SI;

        if ( totalMass < gridParam->MMin && totalMass2 < gridParam->MMin )
        {
          *position = Out;
        }
      }
      break;

    default:
      ABORT(status, 999, "Invalid choice for enum InspiralBankMassRange"); 
      break;
  }
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}



void
LALPopulateCellUnphysicalEta(
		LALStatus               *status,
		InspiralMomentsEtc      *moments,
		InspiralCell            **cell, 
		INT4                     headId,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam,
		CellEvolution           *cellEvolution, 
		CellList		**cellList
		)
{
  REAL4 dx0, dx1;
  REAL4 newt0, newt3;  
  INT4 i, id1, id2;
  REAL4 theta, ctheta, stheta;
  INT4 offSpring;
  INT4 it;
  INT4 add = 0;

  INITSTATUS( status, "LALPopulateCell", 
	      // LALINSPIRALHEXAGONALBANKC );
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* aliases to get the characteristics of the parent template, 
   * that we refer to its ID (headId) */  
  dx0           = (*cell)[headId].dx0;
  dx1           = (*cell)[headId].dx1;
  theta         = (*cell)[headId].metric.theta;
  ctheta        = cos(theta);
  stheta        = sin(theta);
  offSpring     = cellEvolution->nTemplate;

   /* Around the parent, the offspring can be at most 6 (hexagonal grid). 
   * By default the child are unset. If so it is created and have the 
   * properties of its parents. However, a child migh have been created 
   * earlier. In that case, we do not do anything.  */  
  it = 0 ; 

  for (i = 0; i < 6; i++) 
  {
    if ((*cell)[headId].child[i] == -1) 
    {
      add++;
      /* reallocate memory by set of 1000 cells if needed*/
      if ( (offSpring+add)>cellEvolution->nTemplateMax)
      {
        *cell = (InspiralCell*) 
          LALRealloc( *cell,
           sizeof(InspiralCell) * (cellEvolution->nTemplateMax + 1000) );
        if ( !cell ) {
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        cellEvolution->nTemplateMax +=  1000;
      }
      
      /* creates the child connection if needed. A child heritates the
       * properties of its parent */
      switch ( i ){
      case 0:
	newt0   = dx0 ;
	newt3   = 0 ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 *ctheta + stheta* newt3;
	(*cell)[offSpring + it].t3   += newt0 *stheta - ctheta* newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 1:
	newt0   =   dx0/2. ;
	newt3   =   -dx1 *sqrt(3./2) ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 2:
	newt0   =  -dx0/2 ;
	newt3   =  -dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 3:
	newt0   = -dx0 ;
	newt3   = 0;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 4:
	newt0   =  -dx0/2. ;
	newt3   =  dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 5:
	newt0   = dx0/2. ;
	newt3   = dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALInitHexagonalBankUnphysicalEta(status->statusPtr,  cell,  offSpring+it, 
			     moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      }      
      
      /* Now, tricky part, if a child has been creating, he must have a
       * connection with its parents and vice-versa.  */
      if ((*cell)[offSpring + it].child[(i+3)%6] == -1)
	{
	  (*cell)[offSpring + it].child[(i+3)%6] = (*cell)[headId].ID;
	  (*cell)[headId].child[i] = offSpring+it;
	}
      /* a new cell index */
      it += 1;
    }
  }
  
  cellEvolution->nTemplate +=it;
  
  
  /* Here, the parent has its 6 children set; he become sterile. */
  (*cell)[headId].status 	= Sterile;
  (cellEvolution->fertile) 	= cellEvolution->fertile-1;
  LALListDelete(cellList, headId);
  
  /* what shall we do with that parent. Is he valid ? inside the space,
   * outside since eta > 0.25 but close to the boundary .... */     	
  {
    if ((*cell)[headId].RectPosition[0] == Above && (*cell)[headId].in == 1)
      {
	(*cell)[headId].RectPosition[0]=Out;
      }
  }
  
  /* propagate  connections to the brothers to avoid redundancies */  
  for (i=0; i<6; i++)
    {
      /* for each child*/
      id1 = (*cell)[headId].child[i%6];
      id2 = (*cell)[headId].child[(i+1)%6];
      (*cell)[id1].child[(i+2)%6] = (*cell)[id2].ID;
      (*cell)[id2].child[(i+4+1)%6] = (*cell)[id1].ID;   
    }
  
  /* enfin trouver position[0] (In/out)? of the children. */
  for (i=0; i<6; i++)
    {/* for each child find position[0]*/
      id1 = (*cell)[headId].child[i%6];
      
      if ((*cell)[id1].status == Fertile) 
	{
	  LALSPAValidPosition(status->statusPtr, cell, id1, 
			      moments, cellEvolution, cellList);
	  CHECKSTATUSPTR( status );
	  
	  if ((*cell)[id1].position != In ) 
	    {
	      if ((*cell)[id1].status == Fertile) 
		{
		  (*cell)[id1].status= Sterile;
		  cellEvolution->fertile=cellEvolution->fertile-1;
	  	  LALListDelete(cellList, id1);
		}
	    }
	}
    }
  
  DETATCHSTATUSPTR( status );
  RETURN ( status );
}


