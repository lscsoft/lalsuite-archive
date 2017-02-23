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
from pylal.bayespputils import calculate_redshift, DistanceMeasure, lambda_a, amplitudeMeasure

def EnergyScale(lambda_A):
    """
    Calculate mass/energy scale in eV as:    
    m c^2 = h c / \lambda_A
    Valid for alpha != 2
    """
    return 4.135667662E-15*lal.C_SI/lambda_A

def A_LIV(lambda_A, alpha):
    """
    Calculate \mathbb{A}*c^{2-alpha} in (eV)^{2-a} as: 
    A = (m_A*c)^{2-alpha} = (h/lambda_A)^{2-alpha}
    Valid for alpha != 2
    """
    return (4.135667662E-15*lal.C_SI/lambda_A)**(2-alpha)

def lambda_A_of_eff(leffedata, zdata, alphaLIV, cosmology):
    """
    Interface function for lambda_a
    """
    dist = vectorize(lal.LuminosityDistance, excluded=[0])(cosmology, zdata)
    return lambda_a(zdata, alphaLIV, leffedata, dist)

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
  else:
    print "WARNING: Quantile interpolation not possible, edge value used."
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
  lal.SetCosmologicalParametersDefaultValue(cosmology) ## setting h, omega_matter and omega_lambda to their default Planck 2015 values available in LAL

  if not os.path.exists(outfolder):
    os.makedirs(outfolder)
  
  if labels:
    if len(labels)!=len(datafiles):
      print "ERROR: need to give same number of datafiles and labels"
      sys.exit(-1)

  print datafiles
  for (dfile, lab) in zip(datafiles, labels):
    if os.path.splitext(dfile)[1] == '.hdf5':
      #  with h5py.File(dfile, 'r') as h5file:
      h5file = h5py.File(dfile, 'r') 
      ns = h5file['lalinference']['lalinference_nest']['posterior_samples']
      data = array(ns[()])
    else:
      if os.path.splitext(dfile)[1] !=  '.dat':
        print 'WARNING: data format seems to be incompatible...'
      data = genfromtxt(dfile, names=True)


    """Converting (log)distance posterior to meters"""
    if "logdistance" in data.dtype.names:
      distdata = exp(data["logdistance"]) # calculate_redshift needs distances in Mpc. Use * 1e6 * lal.PC_SI to convert to meters
      print "Logarithmic distance parameter detected."
    elif "distance" in data.dtype.names:
      distdata = data["distance"] # calculate_redshift needs distances in Mpc. * 1e6 * lal.PC_SI to convert to meters
      print "Linear distance parameter detected."
    else:
      print "ERROR: No distance posterior! Exiting..."
      sys.exit(-1)

    """Calculating redshifts"""
    zdata = vectorize(calculate_redshift)(distdata, cosmology.h, cosmology.om, cosmology.ol, cosmology.w0)

    #logldata = data["logl"]

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
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, cosmology)
      lameff = True
    elif "log10lambda_eff" in data.dtype.names:
      logleffdata = data["log10lambda_eff"]
      leffdata = pow(10, logleffdata)
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, cosmology)
      loglamAdata = log10(lamAdata)
      lameff = True
    if alphaLIV == 0.0:
        mgdata = EnergyScale(lamAdata)
    if MASSPRIOR:
        # apply uniform mass prior
        print "Weighing lambda_A posterior points by 1/\lambda_\mathbb{A}^2"
        weights = 1.0/lamAdata**2
        print "Weighing loglambda_A posterior points by 1/\lambda_\mathbb{A}"
        logweights = 1.0/lamAdata
    else:
        weights = None
        logweights = None
        

    if alphaLIV < 2.0:
        """Calculating Posterior Quantiles (lower)"""
        PQ_68 = weighted_1dQuantile(0.32, loglamAdata,logweights)
        PQ_90 = weighted_1dQuantile(0.1, loglamAdata, logweights)
        PQ_95 = weighted_1dQuantile(0.05, loglamAdata, logweights)
        PQ_99 = weighted_1dQuantile(0.01, loglamAdata, logweights)
    elif alphaLIV > 2.0:
        """Calculating Posterior Quantiles (upper)"""
        PQ_68 = -weighted_1dQuantile(0.32, -loglamAdata,logweights)
        PQ_90 = -weighted_1dQuantile(0.1, -loglamAdata, logweights)
        PQ_95 = -weighted_1dQuantile(0.05, -loglamAdata, logweights)
        PQ_99 = -weighted_1dQuantile(0.01, -loglamAdata, logweights)
    else:
        print "Cannot handle alpha=2 yet. Exiting..."
        sys.exit(-1)

    # print min(loglamAdata)
    print " Summary"
    print "=-=-=-=-="
    print "shape:", shape(loglamAdata), " min:", min(loglamAdata), " max:", max(loglamAdata)
    print "log(lambda_A)\t68% PQ: ", PQ_68, "\t90% PQ: ", PQ_90, "\t95% PQ: ", PQ_95, "\t99% PQ: ", PQ_99
    print "lambda_A [m]\t68% PQ: ", 10**PQ_68, "\t90% PQ: ", 10**PQ_90, "\t95% PQ: ", 10**PQ_95, "\t99% PQ: ", 10**PQ_99
    print "E_A [eV]\t68% PQ: ", EnergyScale(10**PQ_68), "\t90% PQ: ", EnergyScale(10**PQ_90), "\t95% PQ: ", EnergyScale(10**PQ_95), "\t99% PQ: ", EnergyScale(10**PQ_99)
    print "A [(eV/c)^", str(2-alphaLIV), "]\t68% PQ: ", A_LIV(10**PQ_68, alphaLIV), "\t90% PQ: ", A_LIV(10**PQ_90, alphaLIV), "\t95% PQ: ", A_LIV(10**PQ_95, alphaLIV), "\t99% PQ: ", A_LIV(10**PQ_99, alphaLIV)

    
