from scipy import integrate
from scipy.optimize import newton
import os, sys, optparse, h5py
import numpy as np
import lal
import pickle as pkl

## example command line: 

### python bounds_all_alphas.py -p /home/anuradha/public_html/O2_runs/G268556/LIV/C01/imrphenompv2/plusA/alpha_0p5/lalinferencenest/IMRPhenomPv2threePointFivePN/1167559936.6-0/H1L1/posterior_samples.dat -l plusA_imr_0p5

path0 = os.getcwd()
parser = optparse.OptionParser()

parser.add_option('-p', help='path to posterior file(s) in hdf5 or dat format', dest='pos', action='append',default=None)
parser.add_option('-a', help='mandatory to be passed with hdf5 file as no other way to extract it!', dest='nongralpha', action='store',default=None)
parser.add_option('-l', help='label to check which bound is extracted', dest='lab', action='store',default="alpha_pos")

(opts, args) = parser.parse_args()

if opts.pos is None:
  print "Please input a posterior file in hdf5 or dat format. Exiting..."


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
    z = np.array([newton(find_z_root,0.18,args = (d,omega)) for d in distance])
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
    dl = distance*lal.PC_SI*1e6  ## luminosity distane in metres
    return lambda_eff*(D_alpha/(dl*(1.0+redshift)**(1.0-nonGR_alpha)))**(1./(2.0-nonGR_alpha))

def GetWeightedQuant(dat, qnt, alpha, sample_weight=None,  values_sorted=False, old_style=False):
    #{{{
    # values, quantiles, sample_weight=None, values_sorted=False, old_style=False
    ### taken from here: http://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
   """ Very close to numpy.percentile, but supports weights.
   NOTE: quantiles should be in [0, 1]!
   :param values: numpy.array with data
   :param quantiles: array-like with many quantiles needed
   :param sample_weight: array-like of the same length as `array`
   :param values_sorted: bool, if True, then will avoid sorting of initial array
   :param old_style: if True, will correct output to be consistent with numpy.percentile.
   :return: numpy.array with computed quantiles.
   """
   values = np.copy(dat)
   quantiles = np.array(qnt)
   if sample_weight is None:
       sample_weight = np.ones(len(values))
   #sample_weight = numpy.array(sample_weight)
   assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

   if not values_sorted and alpha < 2.0: ## extract lower bounds for lambda_a when alpha <2
       sorter = np.argsort(values)
       values = values[sorter]
       sample_weight = sample_weight[sorter]
   elif not values_sorted and alpha > 2.0: ## extract upper bounds for lambda_a when alpha>2
       sorter = np.argsort(values)[::-1]
       values = values[sorter]
       sample_weight = sample_weight[sorter]

   weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
   if old_style:
       # To be convenient with numpy.percentile
       weighted_quantiles -= weighted_quantiles[0]
       weighted_quantiles /= weighted_quantiles[-1]
   else:
       weighted_quantiles /= np.sum(sample_weight)
   if alpha <2.0:
     return np.interp(quantiles, weighted_quantiles, values)
   else:
     return np.interp(quantiles, weighted_quantiles[::-1], values[::-1])

def calc_bounds(sampleFiles,lab):
  mpc=1e6*lal.PC_SI
  bounds={}
  bounds[lab]={}
  cdf=0.0
  cdftot=0.0
  RESCALE=1
  for posfile in sampleFiles:
    filename, filext = os.path.splitext(posfile)
    if 'hdf5' in filext:
      ## have to convert from lambdaeff to lambdaA
      if opts.nongralpha is None:
        print "alpha is mandatory to be passed with hdf5 file! Exiting..."
        sys.exit(-1)
      else:
        alpha=float(opts.nongralpha)
        sampObj = h5py.File(posfile)
        loglambda_eff = sampObj['lalinference']['lalinference_nest']['posterior_samples']['log10lambda_eff']
        lambda_eff = np.power(10,loglambda_eff) ## in metres
        logdistance = sampObj['lalinference']['lalinference_nest']['posterior_samples']['logdistance']
        distance = np.exp(logdistance)  ## in Mpc
        z = calculate_redshift(distance)
        loglambdaA = np.log10(lambda_a(z, alpha, lambda_eff, distance))
        minLim = loglambdaA.min()
        maxLim = loglambdaA.max()
    elif 'dat' in filext:
      ## extract lambdaA directly
      with open(posfile) as f:
        params=f.readline().split()
      d = np.loadtxt(posfile,skiprows=1)
      alpha = d[:,params.index("nongr_alpha")][0]
      parid = params.index("log10lambda_a")
      lambda_eff = np.power(10,d[:,params.index("log10lambda_eff")])
      loglambdaA = d[:,parid]
      minLim = loglambdaA.min()
      maxLim = loglambdaA.max()
    else:
      print "file type not understood"
      sys.exit(-1)
      ######################

    pdf,dx = np.histogram(loglambdaA,bins=np.linspace(minLim,maxLim,512),normed=True)
    x = 0.5*(dx[1:]+dx[:-1])
    ##### using method 2 of reweighting samples
    weights = 1./lambda_eff
    CI01 =  GetWeightedQuant(loglambdaA, 0.1, alpha, sample_weight=weights)
    bounds[lab]["nobin"]=pow(10,CI01)
    ###########################################
    if RESCALE:
      pdf/=(10**x) 
      pdf/=(pdf*np.diff(x)[0]).sum()
    #print (pdf*np.diff(x)[0]).sum()
    cdf = (pdf*np.diff(x)[0]).cumsum()
    cdftot += cdf
  if len(sampleFiles)==2:
    cdftot/=2.
  if alpha < 2.0: ## checking values of alpha
    lb = x[np.abs(cdftot-0.1).argmin()]
    bounds[lab]["binning"]=pow(10,lb)
  else:
    ub = x[np.abs(cdftot-0.1).argmax()]
    bounds[lab]["binning"] =pow(10,ub)

  return bounds
LV_bounds = calc_bounds(opts.pos,opts.lab)  ## from uniform prior on mg method

print LV_bounds
