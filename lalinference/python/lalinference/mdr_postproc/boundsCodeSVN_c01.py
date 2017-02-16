from scipy import integrate
from scipy.optimize import newton
import os, sys
import numpy as np
import lal
import pickle as pkl

path0 = os.getcwd()
#LV_resPath='/home/spxas3/public_html/LVC/O2_runs/G268556/LIV/veryPrelim'
#recDirs = ["alpha_0","alpha_1by3","alpha_2by3","alpha_1","alpha_2p5","alpha_3","alpha_4"]

def calculate_redshift(distance,h=0.6790,om=0.3065,ol=0.6935,w0=-1.0):
    """
    Calculate the redshift from the luminosity distance measurement using the
    Cosmology Calculator provided in LAL.
    By default assuming cosmological parameters from arXiv:1502.01589 - 'Planck 2015 results. XIII. Cosmological parameters'
    Using parameters from table 4, column 'TT+lowP+lensing+ext'
    This corresponds to Omega_M = 0.3065, Omega_Lambda = 0.6935, H_0 = 67.90 km s^-1 Mpc^-1
    Returns an array of redshifts
    """
    def find_z_root(z,dl,omega):
        return dl - lal.LuminosityDistance(omega,z)

    omega = lal.CreateCosmologicalParameters(h,om,ol,w0,0.0,0.0) ## initialising with dummy parameters
    lal.SetCosmologicalParametersDefaultValue(omega) ## setting to default Planck 2015 values
    #if isinstance(distance,float):
     #   z = np.array([newton(find_z_root,np.random.uniform(0.0,2.0),args = (distance,omega))])
    #else:
    z = np.array([newton(find_z_root,np.random.uniform(0.0,2.0),args = (d,omega)) for d in distance])
    #return z.reshape(z.shape[0],1)
    return z

def integrand_distance(redshift,nonGR_alpha):
    """
    Calculate D_alpha integral; multiplicative factor put later
    D_alpha = \int ((1+z')^(alpha-2))/sqrt(Omega_m*(1+z')^3 +Omega_lambda) dz' # eq.15 of arxiv 1110.2720
    """
    omega = lal.CreateCosmologicalParameters(0.7,0.3,0.7,-1.0,0.0,0.0) ## these are dummy parameters that are being used for the initialization. they are going to be set to their defaults Planck 2015 values in the next line
    lal.SetCosmologicalParametersDefaultValue(omega) ## setting cosmological parameters to their default Planck 2015 values available on LAL
    omega_m = omega.om # matter density
    omega_l = omega.ol # dark energy density
    #lal.DestroyCosmologicalParameters(omega)
    return (1.0+redshift)**(nonGR_alpha-2.0)/(np.sqrt(omega_m*(1.0+redshift)**3.0 + omega_l))

def DistanceMeasure(redshift,nonGR_alpha):
    """
    D_alpha = ((1+z)^(1-alpha))/H_0 * D_alpha # from eq.15 of arxiv 1110.2720
    D_alpha calculated from integrand in above function
    """
    omega = lal.CreateCosmologicalParameters(0.7,0.3,0.7,-1.0,0.0,0.0) ## these are dummy parameters that are being used for the initialization. they are going to be set to their defaults Planck 2015 values in the next line
    lal.SetCosmologicalParametersDefaultValue(omega) ## setting h, omega_matter and omega_lambda to their default Planck 2015 values available on LAL
    H0 = omega.h*lal.H0FAC_SI ## Hubble constant in SI units
    dist = integrate.quad(integrand_distance, 0, redshift ,args=(nonGR_alpha))[0]
    dist *= (1.0 + redshift)**(1.0 - nonGR_alpha)
    dist /= H0
    #lal.DestroyCosmologicalParameters(omega)
    return dist*lal.C_SI ## returns D_alpha in metres

def lambda_a(redshift, nonGR_alpha, lambda_eff, distance):
    """
    Converting from the effective wavelength-like parameter to \lambda_A:
    \lambda_A = \lambda_{eff}*(D_alpha/D_L)^(1/(2-alpha))*(1/(1+z)^((1-alpha)/(2-alpha)))
    """
    Dfunc = np.vectorize(DistanceMeasure)
    D_alpha = Dfunc(redshift, nonGR_alpha)
    #dl = distance*lal.PC_SI*1e6  ## luminosity distane in metres
    return lambda_eff*(D_alpha/(distance*(1.0+redshift)**(1.0-nonGR_alpha)))**(1./(2.0-nonGR_alpha))

def calc_bounds(sampleFile,distanceFile):
  mpc = 1e6*lal.PC_SI
  bounds={}
  vals=[0,0.5,1,1.5,2.5,3,3.5,4]
  labs=['0','0p5','1','1p5','2p5','3','3p5','4']
  RESCALE=1
  data=pkl.load(open(sampleFile,'r'))
  dist=pkl.load(open(distanceFile,'r')) ## distances are recorded in metres
  for run in ['aplus','aminus']:
    bounds[run]={}
    for (v,l) in zip(vals,labs):
      samps=data[run][v] ## samples for log lambda_eff read off from the file 'samples.pkl'

      ## log lambda_eff converted to log lambda_A as follows:
      distance=dist[run][v] ## luminosity distance in metres read off from file 'distance-samples.pkl'
      lambda_eff = np.power(10,samps) ## lambda_eff converted from it's log
      z = calculate_redshift(distance/mpc) ## converting distance from metres to Mpc
      loglambdaA = np.log10(lambda_a(z, v, lambda_eff, distance))
      ######################

      minLim = loglambdaA.min()
      maxLim = loglambdaA.max()
      #print minLim,maxLim
      #sys.exit(-1)
      pdf,dx = np.histogram(loglambdaA,bins=np.linspace(minLim,maxLim,128),normed=True)
      x = 0.5*(dx[1:]+dx[:-1])

      if RESCALE:
        pdf/=(10**x) 
        #pdf/=(pdf*np.diff(x)[0]).sum()
        cdf = (pdf*np.diff(x)[0]).cumsum()

      if v < 2.0: ## checking values of alpha
        lb = x[np.abs(cdf-0.1).argmin()]
        bounds[run][l]=pow(10,lb)
      else:
        ub = x[np.abs(cdf-0.1).argmax()]
        bounds[run][l] =pow(10,ub)

  return bounds
sampleFile='samples.pkl'
distanceFile='distance-samples.pkl'
LV_bounds = calc_bounds(sampleFile,distanceFile)
print LV_bounds
