#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStdlib.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALInspiralTaylorF2Consistency.h> 

NRCSID (LALINSPIRALWAVEC, "$Id$");

void
LALInspiralTaylorF2Consistency(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{
printf("Hi Dude!\n");


     
    phase = phi0+phi1*(phi2+phi3+phi4+phi5+phi6+phi7)
return;
}

def TaylorF2threePointFivePhase(m1,m2,f):
    gamma = 0.577216 # Euler constant
    mtot = m1+m2
    eta = m1*m2/((m1 + m2)**2)
    Mchirp = mtot*eta**(5./3.)
    v = (numpy.pi*mtot*f)**(1.0/3.0)  
    vlso = (6.**(3./2.))**(-1.)#(6.**(3./2.)*numpy.pi*mtot)**(-1.)

    # calculate the 3.5PN phase. from Buonanno et al. 2009
    
    phi0 = 2.*numpy .pi*tc*f-numpy.pi/4.
    
    phi1 = 3./(128.*eta*v**5) 
    
    phi2 = 1.+(20./9.)*((743./336.)+(11./4.)*v**2.)
    
    phi3 = -16.*numpy.pi*v**3
    
    phi4 = 10.*((3058673.)/(1016064.)+((5429.)/(1008.))*eta + \
            ((617.)/(144.))*eta**2.)*v**4.
    
    phi5 = numpy.pi*((38645.)/(756.)-((65.)/(9.))*eta)*(1. + \
            3.*numpy.log(v/vlso))*v**5.
            
    phi6 = ((11583231236531.)/(4694215680.)-((640.)/(3.))*numpy.pi**2. \
            -((6848.)/(21.))*gamma-((6848.)/(21.))*numpy.log(4.*v)   \
            +(-(15737765635.)/(3048192.)+((2255.)/(12.))*numpy.pi**2.)*eta \
            +((76055.)/(1728.))*eta**2.-((127825.)/(1296.))*eta**3.)*v**6. 
        
    phi7 = numpy.pi*((77096675.)/(254016.)+((378515.)/(1512.))*eta \
            -((74045.)/(756.))*eta**2.)*v**7.