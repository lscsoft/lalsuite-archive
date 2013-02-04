#include <stdio.h>
#include <stdlib.h>
#include "lal_port.h"

/*
 * ============================================================================
 *
 *                                LEGACY CODE
 *
 * ============================================================================
 *
 * DO NOT USE!!!!!
 *
 * FIXME:  REMOVE AS SOON AS POSSIBLE.
 * 
 * FIXED AGAIN: now in lal_port.c. This is a perfectly fine function, with good use case - create window by numeric identifier.
 * 
 */


#define WINDOWH_EEALLOCATE 4
#define WINDOWH_MSGEEALLOCATE "unable to allocate vector to store window"


static REAL4Window *XLALCreateREAL4Window(UINT4 length, WindowType type, REAL4 beta)
{
       static const char func[] = "XLALCreateREAL4Window";
       REAL4Window *window;

       //XLALPrintDeprecationWarning(func, "XLALCreateRectangularREAL4Window(), XLALCreateHannREAL4Window(), XLALCreateWelchREAL4Window(), XLALCreateBartlettREAL4Window(), XLALCreateParzenREAL4Window(), XLALCreatePapoulisREAL4Window(), XLALCreateHammingREAL4Window(), XLALCreateKaiserREAL4Window(), XLALCreateCreightonREAL4Window(), or XLALCreateTukeyREAL4Window");

       switch (type) {
       case Rectangular:
               window = XLALCreateRectangularREAL4Window(length);
               break;

       case Hann:
               window = XLALCreateHannREAL4Window(length);
               break;

       case Welch:
               window = XLALCreateWelchREAL4Window(length);
               break;

       case Bartlett:
               window = XLALCreateBartlettREAL4Window(length);
               break;

       case Parzen:
               window = XLALCreateParzenREAL4Window(length);
               break;

       case Papoulis:
               window = XLALCreatePapoulisREAL4Window(length);
               break;

       case Hamming:
               window = XLALCreateHammingREAL4Window(length);
               break;

       case Kaiser:
               window = XLALCreateKaiserREAL4Window(length, beta);
               break;

       case Creighton:
               window = XLALCreateCreightonREAL4Window(length, beta);
               break;

       case Tukey:
               window = XLALCreateTukeyREAL4Window(length, beta);
               break;

       default:
               return NULL;
       }

       return window;
}

void LALWindow(LALStatus *status, REAL4Vector *vector, LALWindowParams *parameters)
{
/* emulate legacy code */
REAL4Window *window;

// INITSTATUS(status, "LALWindow", WINDOWC);
// ATTATCHSTATUSPTR(status);

window = XLALCreateREAL4Window(parameters->length, parameters->type, parameters->beta);
if(!window) {
	// ABORT(status, WINDOWH_EEALLOCATE, WINDOWH_MSGEEALLOCATE);
	fprintf(stderr, "*** ERROR: Could not create requested window (length=%d type=%d parameter=%g)\n", parameters->length, parameters->type, parameters->beta);
	exit(-1);
	}
memcpy(vector->data, window->data->data, vector->length * sizeof(*vector->data));
XLALDestroyREAL4Window(window);

// DETATCHSTATUSPTR(status);
// RETURN(status);
}
