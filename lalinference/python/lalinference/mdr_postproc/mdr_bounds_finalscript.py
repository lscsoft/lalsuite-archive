#/bin/python

'''

Modified Dispersion Relation post-processing module for PE 

'''

import sys
import os
import argparse

import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy import stats
from scipy import integrate
from scipy.optimize import newton 
import random 
import h5py
import lal
import pickle

def MassScale(lambda_A):
    """
    Calculate mass/energy scale in eV as:    
    m c^2 = h c / \lambda_A
    Valid for alpha != 2
    """
    return lal.c_SI*lal.h_SI/lambda_A


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

    z = np.array([newton(find_z_root,np.random.uniform(0.0,2.0),args = (d,omega)) for d in distance])
    return z


def weighted_1dQuantile(quantile, data, weights=None):
  '''Calculate quantiles for confidence intervals on sorted weighted posteriors'''

  if weights is None:
    norm = len(data)
    weights = ones(len(data))
  else:
    norm = sum(weights)

  ### Sort data and weights array
  wdata = array(zip(data, weights) ,dtype=[('data','<f8'),('weights','<f8')])

  swdata = sort(wdata, order='data')
  sdata = swdata['data']
  sweights = swdata['weights']

  cumw = cumsum(sweights)/norm
  iplus = where(cumw > quantile)[0][0]
  x = sdata[iplus]
  # Interpolate if possible
  if iplus > 0:
    x -= (cumw[iplus] - quantile) * (sdata[iplus] - sdata[iplus-1])/(cumw[iplus] - cumw[iplus-1])
  return x


if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################

  parser = argparse.ArgumentParser(description="This is a post-processing script for calculating bounds on parameterized modified dispersion relations.")
  parser.add_argument("-i", "--input", type=str, dest="datafiles", nargs='+', help="list of .hdf5 or .dat file(s) containing data points (one file per chain)", metavar="POSTERIOR FILES")
  parser.add_argument("-l", "--label", type=str, dest="labels", nargs='+', help="source-identifying string(s)", default=None)
  parser.add_argument("-o", "--output", type=str, dest="outputfolder", help="outputfolder", metavar="OUTPUT FOLDER",default=".")
  parser.add_argument("--mprior", action="store_true", dest="mprior", help="use uniform mass prior")
  parser.add_argument("-a", "--alpha", type=float, dest="alphaLIV", help="Exponent of Lorentz invariance violating term", default=0.0)


  args = parser.parse_args()
  
  datafiles = args.datafiles
  labels = args.labels
  outfolder = args.outputfolder
  MASSPRIOR = args.mprior
  alphaLIV = args.alphaLIV

    cosmology = lal.CreateCosmologicalParameters(0.7,0.3,0.7,-1.0,0.0,0.0) ## these are dummy parameters that are being used for the initialization. they are going to be set to their defaults Planck 2015 values in the next line
    lal.SetCosmologicalParametersDefaultValue(cosmology) ## setting h, omega_matter and omega_lambda to their default Pla
    #this_cosmology = LCDM_Cosmology(67.90, 0.3065, 0.6935)

  if not os.path.exists(outfolder):
    os.makedirs(outfolder)
  
  if labels:
    if len(labels)!=len(datafiles):
      print "ERROR: need to give same number of datafiles and labels"
      sys.exit(-1)

  print datafiles
  for (dfile, lab) in zip(datafiles, labels):
    if os.path.splitext(dfile)[1] is '.hdf5':
      #  with h5py.File(dfile, 'r') as h5file:
      h5file = h5py.File(dfile, 'r') 
      ns = h5file['lalinference']['lalinference_nest']['posterior_samples']
      data = array(ns[()])
    else:
      if os.path.splitext(dfile)[1] is not '.dat':
        print 'WARNING: data format seems to be incompatible...'
      data = genfromtxt(dfile, names=True)


    """Converting (log)distance posterior to meters"""
    if "logdistance" in data.dtype.names:
      distdata = exp(data["logdistance"]) * 1e6 * lal.PC_SI
      print "Logarithmic distance parameter detected."
    elif "distance" in data.dtype.names:
      distdata = data["distance"] * 1e6 * lal.PC_SI
      print "Linear distance parameter detected."
    else:
      print "ERROR: No distance posterior! Exiting..."
      sys.exit(-1)

    """Calculating redshifts"""
    zdata = calculate_redshift(distdata)

    logldata = data["logL"]

    """Calculating posteriors for lambda_{eff} parameters"""
    if "log10lambda_a" in data.dtype.names:
      loglamAdata = data["log10lambda_a"]
      lamAdata = pow(10, loglamAdata)
    elif "lambda_a" in data.dtype.names:
      lamAdata = data["lambda_a"]
      loglamAdata = log10(lamAdata)
    elif "lambda_eff" in data.dtype.names:
      leffdata = data["lambda_eff"]
      logleffdata = log10(leffdata)
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, this_cosmology)
      lameff = True
    elif "log10lambda_eff" in data.dtype.names:
      logleffdata = data["log10lambda_eff"]
      leffdata = pow(10, logleffdata)
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, this_cosmology)
      loglamAdata = log10(lamAdata)
      lameff = True
    if alphaLIV == 0.0:
        mgdata = GravitonMass(lamAdata)
    if MASSPRIOR:
        # apply uniform mass prior
        print "Weighing posterior points by 1/\lambda_\mathbb{A}^2"
        print "TODO: Update for arbitrary values of alpha!"
        weights = 1.0/lamAdata**2
        logweights = 1.0/lamAdata
    else:
        weights = None
        logweights = None

    CI_68 = weighted_1dQuantile(0.32, loglamAdata,logweights)
    CI_90 = weighted_1dQuantile(0.1, loglamAdata, logweights)
    CI_95 = weighted_1dQuantile(0.05, loglamAdata, logweights)
    CI_99 = weighted_1dQuantile(0.01, loglamAdata, logweights)

    # print min(loglamAdata)
    print " Summary"
    print "=-=-=-=-="
    print "shape:", shape(loglamAdata), " min:", min(loglamAdata), " max:", max(loglamAdata)
    print "\t68% CI: ", CI_68, "\t90% CI: ", CI_90, "\t95% CI: ", CI_95, "\t99% CI: ", CI_99

    
