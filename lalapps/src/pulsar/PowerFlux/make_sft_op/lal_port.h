#ifndef __LAL_PORT_H__
#define __LAL_PORT_H__
#include <lal/Window.h>


/* 
 * Forward port LAL* functions that were replaced by transition to XLAL
 */


typedef enum {
       Rectangular,
       Hann,
       Welch,
       Bartlett,
       Parzen,
       Papoulis,
       Hamming,
       Kaiser,
       Creighton,
       Tukey
	} WindowType;
	
typedef struct tagLALWindowParams {
       INT4        length;
       WindowType  type;
       REAL4       beta;
	} LALWindowParams;

void LALWindow(LALStatus *, REAL4Vector *, LALWindowParams *);
// Don't need it
// void LALCreateREAL4Window(LALStatus *, REAL4Window **, LALWindowParams *);

#endif