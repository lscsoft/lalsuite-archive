#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Sequence.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALMalloc.h>
#include <lal/LALNoiseModels.h>
#include <lal/MatrixUtils.h>
#include <lal/FindChirpPTF.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALInspiral.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/LALDetectors.h>
#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALError.h>
	
/* PARAMETER STRUCTURE FOR GENERATEPPNAMPCORCONSISTENCY */
typedef struct tagPPNConsistencyParamStruc {
	/* Passed parameters. */
	SkyPosition position; /* location of source on sky */
	REAL4 psi;            /* polarization angle (radians) */
	LIGOTimeGPS epoch;    /* start time of output time series */
	
	/* Input parameters. */
	REAL8 mTot_real8; /* total system mass (Msun) */
	REAL8 eta_real8;  /* mass ratio */
	REAL8 delta;      /* sqrt(1-4eta) */
	REAL4 mTot;       /* total system mass (Msun) */
	REAL4 eta;        /* mass ratio */
	REAL4 d;          /* distance (metres) */
	REAL4 inc;        /* inclination angle (radians) */
	REAL4 cosI;				/* cosine of inclination angle */
	REAL4 sinI;				/* sine of inclination angle */
	REAL4 phi;        /* coalescence phase (radians) */
	REAL8 deltaT;     /* requested sampling interval (s) */
	REAL4 fStartIn;   /* requested start frequency (Hz) */
	REAL4 fStopIn;    /* requested stop frequency (Hz) */
	UINT4 lengthIn;   /* maximum length of waveform */
	REAL4Vector *ppn; /* post-Newtonian selection parameters */
	INT4 ampOrder;    /* PN amplitude selection 0-5 */
	
	/* PN phasing coefficients for use in AmpCorConsistency */
	REAL4 phi0, phi2, phi3, phi4, phi5, phi5l, phi6, phi6l, phi7;
	REAL4 phasePNparams[9];
	
	/* Output parameters. */
	REAL8 tc;         /* time to coalescence from start of waveform */
	REAL4 dfdt;       /* maximum value of df*dt over any timestep */
	REAL4 fStart;     /* actual start frequency (Hz) */
	REAL4 fStop;      /* actual stop frequency (Hz) */
	UINT4 length;     /* length of signal generated */
	INT4 termCode;    /* termination code */
	const CHAR *termDescription; /* description of termination code */
} PPNConsistencyParamStruc;
	
	
	
/* FUNCTION PROTOTYPES */	
void LALGeneratePPNAmpCorConsistency(
                              LALStatus     *stat,
                              CoherentGW    *output,
                              PPNConsistencyParamStruc *params
                            );

/* function to populate phaseParams */                            
void LALPopulatePhasePNparams(PPNConsistencyParamStruc params, INT4 TestParam);
