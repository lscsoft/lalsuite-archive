#include <stdio.h>
#include <stderr.h>
#include "lal_port.h"

void LALWindow(LALStatus *, REAL4Vector *, LALWindowParams *)
{
/* emulate legacy code */
REAL4Window *window;

INITSTATUS(status, "LALWindow", WINDOWC);
ATTATCHSTATUSPTR(status);

window = XLALCreateREAL4Window(parameters->length, parameters->type, parameters->beta);
if(!window) {
	ABORT(status, WINDOWH_EEALLOCATE, WINDOWH_MSGEEALLOCATE);
}
memcpy(vector->data, window->data->data, vector->length * sizeof(*vector->data));
XLALDestroyREAL4Window(window);

DETATCHSTATUSPTR(status);
RETURN(status);
}
