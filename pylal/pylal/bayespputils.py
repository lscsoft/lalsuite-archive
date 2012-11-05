# -*- coding: utf-8 -*-
#
#       bayespputils.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#       Salvatore Vitale <salvatore.vitale@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#===============================================================================
# Preamble
#===============================================================================

"""
This module contains classes and functions for post-processing the output
of the Bayesian parameter estimation codes.
"""

#standard library imports
import os
import sys
from math import ceil,floor,sqrt,pi as pi_constant
import xml
from xml.dom import minidom
from operator import itemgetter

#related third party imports
import numpy as np
from numpy import fmod
import matplotlib
from matplotlib import pyplot as plt,cm as mpl_cm,lines as mpl_lines
from scipy import stats
from scipy import special
from scipy import signal
from numpy import linspace
import random

from matplotlib.ticker import FormatStrFormatter,ScalarFormatter,AutoMinorLocator

# Default font properties
fig_width_pt = 246  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 11,
        'text.fontsize':   11,
        'legend.fontsize': 11,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'text.usetex': False,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ['Times','Palatino','New Century Schoolbook','Bookman','Computer Modern Roman','Times New Roman','Liberation Serif'],
        'font.weight':'normal',
        'font.size':11,
        'savefig.dpi': 120
        })


try:
    from xml.etree.cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser
except ImportError:
    #Python < 2.5
    from cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser

#local application/library specific imports
import pylal
from pylal import lalconstants
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import git_version
#C extensions
from _bayespputils import _skyhist_cart,_burnin

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

#===============================================================================
# Constants
#===============================================================================
#Parameters which are not to be exponentiated when found
logParams=['logl','loglh1','loglh2','logll1','loglv1','deltalogl','deltaloglh1','deltalogll1','deltaloglv1','logw']
#Pre-defined ordered list of line styles for use in matplotlib contour plots.
__default_line_styles=['solid', 'dashed', 'dashdot', 'dotted']
#Pre-defined ordered list of matplotlib colours for use in plots.
__default_color_lst=['r','b','y','g','c','m']
#A default css string for use in html results pages.
__default_css_string="""
p,h1,h2,h3,h4,h5
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
}

p
{
font-size:14px;
}

h1
{
font-size:20px;
}

h2
{
font-size:18px;
}

h3
{
font-size:16px;
}



table
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
width:100%;
border-collapse:collapse;
}
td,th
{
font-size:12px;
border:1px solid #B5C1CF;
padding:3px 7px 2px 7px;
}
th
{
font-size:14px;
text-align:left;
padding-top:5px;
padding-bottom:4px;
background-color:#B3CEEF;
color:#ffffff;
}
#postable tr:hover
{
background: #DFF4FF;
}
#covtable tr:hover
{
background: #DFF4FF;
}
#statstable tr:hover
{
background: #DFF4FF;
}

.ppsection
{
border-bottom-style:double;
}

"""

#===============================================================================
# Functions used to parse injection structure.
#===============================================================================
def _inj_m1(inj):
    """
    Return the mapping of (mchirp,eta)->m1; m1>m2 i.e. return the greater of the mass 
    components (m1) calculated from the chirp mass and the symmetric mass ratio.
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attributes 'mchirp' and 'eta'.
    @rtype: number
    """
    (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
    return mass1
def _inj_m2(inj):
    """
    Return the mapping of (mchirp,eta)->m2; m1>m2 i.e. return the lesser of the mass 
    components (m2) calculated from the chirp mass and the symmetric mass ratio.
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attributes 'mchirp' and 'eta'.
    @rtype: number
    """
    (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
    return mass2

def _inj_q(inj):
    """
    Return the mapping of (mchirp,eta)->q; m1>m2 i.e. return the mass ratio q=m2/m1.
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attributes 'mchirp' and 'eta'.
    @rtype: number 
    """
    (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
    return mass2/mass1

def _inj_longitude(inj):
    """
    Return the mapping of longitude found in inj to the interval [0,2*pi).
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attribute 'longitude'.
    @rtype: number
    """
    if inj.longitude>2*pi_constant or inj.longitude<0.0:
        maplong=2*pi_constant*(((float(inj.longitude))/(2*pi_constant)) - floor(((float(inj.longitude))/(2*pi_constant))))
        print "Warning: Injected longitude/ra (%s) is not within [0,2\pi)! Angles are assumed to be in radians so this will be mapped to [0,2\pi). Mapped value is: %s."%(str(inj.longitude),str(maplong))
        return maplong
    else:
        return inj.longitude

def _inj_a1(inj):
    """
    Return the magnitude of the spin 1 vector. Calculates the spin magnitude
    from it's components.
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attribute 'spin1x','spin1y', and 'spin1z' (the spin components).
    @rtype: number
    """
    x = inj.spin1x
    y = inj.spin1y
    z = inj.spin1z
    return sqrt(x*x + y*y + z*z)

def _inj_a2(inj):
    """
    Return the magnitude of the spin 2 vector. Calculates the spin magnitude
    from it's components.
    
    @type inj: glue.ligolw.lsctables.SimInspiral
    @param inj: a custom type with the attribute 'spin2x','spin2y', and 'spin2z' (the spin components).
    @rtype: number
    """
    x = inj.spin2x
    y = inj.spin2y
    z = inj.spin2z
    return sqrt(x*x + y*y + z*z)

def _inj_theta1(inj):
    x = inj.spin1x
    y = inj.spin1y
    z = inj.spin1z
    if x == 0.0 and y == 0.0 and z == 0.0:
        return None
    else:
        return np.arccos( z / sqrt(x*x+y*y+z*z) )

def _inj_theta2(inj):
    x = inj.spin2x
    y = inj.spin2y
    z = inj.spin2z
    if x == 0.0 and y == 0.0 and z == 0.0:
        return None
    else:
        return np.arccos( z / sqrt(x*x+y*y+z*z) )

def _inj_phi1(inj):
    x = inj.spin1x
    y = inj.spin1y
    z = inj.spin1z
    if x == 0.0 and y == 0.0 and z == 0.0:
        return None
    else:
        phi_mpi_to_pi = np.arctan2(y, x)
        if phi_mpi_to_pi < 0.0:
            return phi_mpi_to_pi + 2*pi_constant
        else:
            return phi_mpi_to_pi

def _inj_phi2(inj):
    x = inj.spin2x
    y = inj.spin2y
    z = inj.spin2z
    if x == 0.0 and y == 0.0 and z == 0.0:
        return None
    else:
        phi_mpi_to_pi = np.arctan2(y, x)
        if phi_mpi_to_pi < 0.0:
            return phi_mpi_to_pi + 2*pi_constant
        else:
            return phi_mpi_to_pi

def _inj_tilt1(inj):
    S1  = np.hstack((inj.spin1x,inj.spin1y,inj.spin1z))
    L  = orbital_momentum(inj.f_lower, inj.mchirp, inj.inclination)
    tilt1 = array_ang_sep(L,S1)
    if np.max(S1) == 0.0:
        return None
    else:
        return tilt1

def _inj_tilt2(inj):
    S2  = np.hstack((inj.spin2x,inj.spin2y,inj.spin2z))
    L  = orbital_momentum(inj.f_lower, inj.mchirp, inj.inclination)
    tilt2 = array_ang_sep(L,S2)
    if np.max(S2) == 0.0:
        return None
    else:
        return tilt2

def _inj_thetas(inj):
    L  = orbital_momentum(inj.f_lower, inj.mchirp, inj.inclination)
    S1  = inj.mass1*inj.mass1*np.hstack((inj.spin1x,inj.spin1y,inj.spin1z))
    S2  = inj.mass2*inj.mass2*np.hstack((inj.spin2x,inj.spin2y,inj.spin2z))
    J = L + S1 + S2

    thetas = array_polar_ang(J)
    return thetas
    
def _inj_beta(inj):
    L  = orbital_momentum(inj.f_lower, inj.mchirp, inj.inclination)
    S1  = inj.mass1*inj.mass1*np.hstack((inj.spin1x,inj.spin1y,inj.spin1z))
    S2  = inj.mass2*inj.mass2*np.hstack((inj.spin2x,inj.spin2y,inj.spin2z))
    J = L + S1 + S2
    
    beta  = array_ang_sep(J,L)
    return beta


#===============================================================================
# Class definitions
#===============================================================================

class PosteriorOneDPDF(object):
    """
    A data structure representing one parameter in a chain of posterior samples. 
    The Posterior class generates instances of this class for pivoting onto a given
    parameter (the Posterior class is per-Sampler oriented whereas this class represents
    the same one parameter in successive samples in the chain).
    """
    def __init__(self,name,posterior_samples,injected_value=None,trigger_values=None,prior=None):
        """
        Create an instance of PosteriorOneDPDF based on a table of posterior_samples.

        @param name: A literal string name for the parameter.
        @param posterior_samples: A 1D array of the samples.
        @keyword injected_value: The injected or real value of the parameter.
        @type injected_value: glue.ligolw.lsctables.SimInspiral
        @keyword trigger_values: The trigger values of the parameter (dictionary w/ IFOs as keys).
        @keyword prior: The prior value corresponding to each sample.
        """
        self.__name=name
        self.__posterior_samples=np.array(posterior_samples)

        self.__injval=injected_value
        self.__trigvals=trigger_values
        self.__prior=prior

        return

    def __len__(self):
        """
        Container method. Defined as number of samples.
        """
        return len(self.__posterior_samples)

    def __getitem__(self,idx):
        """
        Container method . Returns posterior containing sample idx (allows slicing).
        """
        return PosteriorOneDPDF(self.__name, self.__posterior_samples[idx], injected_value=self.__injval, trigger_values=self.__trigvals)
        
    @property
    def name(self):
        """
        Return the string literal name of the parameter.
        
        @rtype: string
        """
        return self.__name

    @property
    def mean(self):
        """
        Return the arithmetic mean for the marginal PDF on the parameter.
        
        @rtype: number
        """
        return np.mean(self.__posterior_samples)

    @property
    def median(self):
        """
        Return the median value for the marginal PDF on the parameter.
        
        @rtype: number
        """
        return np.median(self.__posterior_samples)

    @property
    def stdev(self):
        """
        Return the standard deviation of the marginal PDF on the parameter.
        
        @rtype: number
        """
        try:
            stdev = sqrt(np.var(self.__posterior_samples))
            if not np.isfinite(stdev): 
                raise OverflowError
        except OverflowError:
            mean = np.mean(self.__posterior_samples)
            stdev = mean * sqrt(np.var(self.__posterior_samples/mean))
        return stdev

    @property
    def stacc(self):
        """
        Return the 'standard accuracy statistic' (stacc) of the marginal 
        posterior of the parameter. 
        
        stacc is a standard deviant incorporating information about the 
        accuracy of the waveform recovery. Defined as the mean of the sum 
        of the squared differences between the points in the PDF 
        (x_i - sampled according to the posterior) and the true value 
        (x_{\rm true}).  So for a marginalized one-dimensional PDF:
        stacc = \sqrt{\frac{1}{N}\sum_{i=1}^N (x_i-x_{\rm true})2}
        
        @rtype: number
        """
        if self.__injval is None:
            return None
        else:
            return np.sqrt(np.mean((self.__posterior_samples - self.__injval)**2.0))

    @property
    def injval(self):
        """
        Return the injected value set at construction . If no value was set
        will return None .
        
        @rtype: undefined
        """
        return self.__injval

    @property
    def trigvals(self):
        """
        Return the trigger values set at construction. If no value was set
        will return None .
        
        @rtype: undefined
        """
        return self.__trigvals

    #@injval.setter #Python 2.6+
    def set_injval(self,new_injval):
        """
        Set the injected/real value of the parameter.

        @param new_injval: The injected/real value to set.
        @type new_injval: glue.ligolw.lsctables.SimInspiral
        """

        self.__injval=new_injval

    def set_trigvals(self,new_trigvals):
        """
        Set the trigger values of the parameter.

        @param new_trigvals: Dictionary containing trigger values with IFO keys.
        @type new_trigvals: glue.ligolw.lsctables.SnglInspiral
        """

        self.__trigvals=new_trigvals

    @property
    def samples(self):
        """
        Return a 1D numpy.array of the samples.
        
        @rtype: numpy.array
        """
        return self.__posterior_samples

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from posterior, analagous to numpy.delete but opperates in place.

        @param samples: A list containing the indexes of the samples to be removed.
        @type samples: list
        """
        self.__posterior_samples=np.delete(self.__posterior_samples,samples).reshape(-1,1)

    @property
    def gaussian_kde(self):
        """
        Return a SciPy gaussian_kde (representing a Gaussian KDE) of the samples.
        
        @rtype: scipy.stats.kde.gaussian_kde
        """
        from numpy import seterr as np_seterr
        from scipy import seterr as sp_seterr

        np_seterr(under='ignore')
        sp_seterr(under='ignore')
        try:
            return_value=stats.kde.gaussian_kde(np.transpose(self.__posterior_samples))
        except:
            exfile=open('exception.out','w')
            np.savetxt(exfile,self.__posterior_samples)
            exfile.close()
            raise

        return return_value

    def prob_interval(self,intervals):
        """
        Evaluate probability intervals.

        @param intervals: A list of the probability intervals [0-1]
        @type intervals: list
        """
        list_of_ci=[]
        samples_temp=np.sort(np.squeeze(self.samples))

        for interval in intervals:
            if interval<1.0:
                samples_temp
                N=np.size(np.squeeze(self.samples))
                #Find index of lower bound
                lower_idx=int(floor((N/2.0)*(1-interval)))
                if lower_idx<0:
                    lower_idx=0
                #Find index of upper bound
                upper_idx=N-int(floor((N/2.0)*(1-interval)))
                if upper_idx>N:
                    upper_idx=N-1

                list_of_ci.append((float(samples_temp[lower_idx]),float(samples_temp[upper_idx])))
            else:
                list_of_ci.append((None,None))

        return list_of_ci


class Posterior(object):
    """
    Data structure for a table of posterior samples .
    """
    def __init__(self,commonResultsFormatData,SimInspiralTableEntry=None,SnglInpiralList=None,name=None,description=None,votfile=None):
        """
        Constructor.

        @param commonResultsFormatData: A 2D array containing the posterior
            samples and related data. The samples chains form the columns.
        @type commonResultsFormatData: custom type
        @param SimInspiralTableEntry: A SimInspiralTable row containing the injected values.
        @type SimInspiralTableEntry: glue.ligolw.lsctables.SimInspiral
        @param SnglInspiralList: A list of SnglInspiral objects containing the triggers.
        @type SnglInspiralList: list
        
        """
        common_output_table_header,common_output_table_raw =commonResultsFormatData
        self._posterior={}
        self._injection=SimInspiralTableEntry
        self._triggers=SnglInpiralList
        self._loglaliases=['posterior', 'logl','logL','likelihood']
        self._votfile=votfile
        
        common_output_table_header=[i.lower() for i in common_output_table_header]
        
        for one_d_posterior_samples,param_name in zip(np.hsplit(common_output_table_raw,common_output_table_raw.shape[1]),common_output_table_header):
            
            self._posterior[param_name]=PosteriorOneDPDF(param_name.lower(),one_d_posterior_samples,injected_value=self._getinjpar(param_name),trigger_values=self._gettrigpar(param_name))

        if 'mchirp' in common_output_table_header and 'eta' in common_output_table_header \
        and (not 'm1' in common_output_table_header) and (not 'm2' in common_output_table_header):
            try:
                print 'Inferring m1 and m2 from mchirp and eta'
                (m1,m2)=mc2ms(self._posterior['mchirp'].samples, self._posterior['eta'].samples)
                self._posterior['m1']=PosteriorOneDPDF('m1',m1,injected_value=self._getinjpar('m1'),trigger_values=self._gettrigpar('m1'))
                self._posterior['m2']=PosteriorOneDPDF('m2',m2,injected_value=self._getinjpar('m2'),trigger_values=self._gettrigpar('m2'))
            except KeyError:
                print 'Unable to deduce m1 and m2 from input columns'

        
        logLFound=False
        
        for loglalias in self._loglaliases:
        
            if loglalias in common_output_table_header:
                try:
                    self._logL=self._posterior[loglalias].samples
                except KeyError:
                    print "No '%s' column in input table!"%loglalias
                    continue
                logLFound=True
                
        if not logLFound:
            raise RuntimeError("No likelihood/posterior values found!")

        if name is not None:
            self.__name=name
            
        if description is not None:
            self.__description=description

        return

    def bootstrap(self):
        """
        Returns a new Posterior object that contains a bootstrap
        sample of self.
        
        @rtype: Posterior
        """
        names=[]
        samples=[]
        for name,oneDpos in self._posterior.items():
            names.append(name)
            samples.append(oneDpos.samples)

        samplesBlock=np.hstack(samples)

        bootstrapSamples=samplesBlock[:,:]
        Nsamp=bootstrapSamples.shape[0]

        rows=np.vsplit(samplesBlock,Nsamp)

        for i in range(Nsamp):
            bootstrapSamples[i,:]=random.choice(rows)

        return Posterior((names,bootstrapSamples),self._injection,self._triggers)

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from all OneDPosteriors.

        @param samples: The indexes of the samples to be removed.
        @type samples: list
        """
        for name,pos in self:
            pos.delete_samples_by_idx(samples)
        return

    def delete_NaN_entries(self,param_list):
        """
        Remove samples containing NaN in request params.

        @param param_list: The parameters to be checked for NaNs.
        @type param_list: list
        """
        nan_idxs = np.array(())
        nan_dict = {}
        for param in param_list:
            nan_bool_array = np.isnan(self[param].samples).any(1)
            idxs = np.where(nan_bool_array == True)[0]
            if len(idxs) > 0:
                nan_dict[param]=len(idxs)
                nan_idxs = np.append(nan_idxs, idxs)
        total_samps = len(self)
        nan_samps   = len(nan_idxs)
        if nan_samps is not 0:
            print "WARNING: removing %i of %i total samples due to NaNs:"% (nan_samps,total_samps)
            for param in nan_dict.keys():
                print "\t%i NaNs in %s."%(nan_dict[param],param)
            self.delete_samples_by_idx(nan_idxs)
        return

    @property
    def injection(self):
        """
        Return the injected values.
        
        @rtype: glue.ligolw.lsctables.SimInspiral
        """

        return self._injection

    @property
    def triggers(self):
        """
        Return the trigger values .
        
        @rtype: list
        """

        return self._triggers

    def _total_incl_restarts(self, samples):
        total=0
        last=samples[0]
        for x in samples[1:]:
            if x < last:
                total += last
            last = x
        total += samples[-1]
        return total

    def longest_chain_cycles(self):
        """
        Returns the number of cycles in the longest chain
        
        @rtype: number
        """
        samps,header=self.samples()
        header=header.split()
        if not ('cycle' in header):
            raise RuntimeError("Cannot compute number of cycles in longest chain")
        if 'chain' in header:
            chain_col=header.index('chain')
            cycle_col=header.index('cycle')
            chain_indexes=np.unique(samps[:,chain_col])
            max_cycle=0
            for ind in chain_indexes:
                chain_cycle_samps=samps[ samps[:,chain_col] == ind, cycle_col ]
                max_cycle=max(max_cycle, self._total_incl_restarts(chain_cycle_samps))
            return int(max_cycle)
        else:
            return int(self._total_incl_restarts(samps[:,cycle_col]))

    #@injection.setter #Python 2.6+
    def set_injection(self,injection):
        """
        Set the injected values of the parameters.

        @param injection: A SimInspiralTable row object containing the injected parameters.
        @type injection: glue.ligolw.lsctables.SimInspiral
        """
        if injection is not None:
            self._injection=injection
            for name,onepos in self:
                new_injval=self._getinjpar(name)
                if new_injval is not None:
                    self[name].set_injval(new_injval)

    def set_triggers(self,triggers):
        """
        Set the trigger values of the parameters.

        @param triggers: A list of SnglInspiral objects.
        @type triggers: list
        """
        if triggers is not None:
            self._triggers=triggers
            for name,onepos in self:
                new_trigvals=self._gettrigpar(name)
                if new_trigvals is not None:
                    self[name].set_trigvals(new_trigvals)


    _injXMLFuncMap={
                        'mchirp':lambda inj:inj.mchirp,
                        'chirpmass':lambda inj:inj.mchirp,
                        'mc':lambda inj:inj.mchirp,
                        'mass1':_inj_m1,
                        'm1':_inj_m1,
                        'mass2':_inj_m2,
                        'm2':_inj_m2,
                        'eta':lambda inj:inj.eta,
                        'q':_inj_q,
                        'asym_massratio':_inj_q,
                        'time': lambda inj:float(inj.get_end()),
                        'end_time': lambda inj:float(inj.get_end()),
                        'phi0':lambda inj:inj.phi0,
                        'phi_orb': lambda inj: inj.coa_phase,
                        'dist':lambda inj:inj.distance,
                        'distance':lambda inj:inj.distance,
                        'ra':_inj_longitude,
                        'rightascension':_inj_longitude,
                        'long':_inj_longitude,
                        'longitude':_inj_longitude,
                        'dec':lambda inj:inj.latitude,
                        'declination':lambda inj:inj.latitude,
                        'lat':lambda inj:inj.latitude,
                        'latitude':lambda inj:inj.latitude,
                        'psi': lambda inj: np.mod(inj.polarization, np.pi),
                        'iota':lambda inj: inj.inclination,
                        'inclination': lambda inj: inj.inclination,
                        'spinchi': lambda inj: (inj.spin1z + inj.spin2z) + sqrt(1-4*inj.eta)*(inj.spin1z - spin2z),
                        'f_lower': lambda inj: inj.f_lower,
                        'a1': lambda inj: np.sqrt(inj.spin1x**2+inj.spin1y**2+inj.spin1z**2) ,
                        'a2': lambda inj: np.sqrt(inj.spin2x**2+inj.spin2y**2+inj.spin2z**2) ,
                        'theta1':_inj_theta1,
                        'theta2':_inj_theta2,
                        'phi1':_inj_phi1,
                        'phi2':_inj_phi2,
                        'tilt1':_inj_tilt1,
                        'tilt2':_inj_tilt2,
                        'costilt1': lambda inj: np.cos(_inj_tilt1),
                        'costilt2': lambda inj: np.cos(_inj_tilt2),
                        'cos(iota)': lambda inj: np.cos(inj.inclination),
                        'thetas':_inj_thetas,
                        'beta':_inj_beta
                       }

    def _getinjpar(self,paramname):
        """
        Map parameter names to parameters in a SimInspiralTable .
        """
        if self._injection is not None:
            for key,value in self._injXMLFuncMap.items():
                if paramname.lower().strip() == key.lower().strip():
                    return self._injXMLFuncMap[key](self._injection)
        return None

    def _gettrigpar(self,paramname):
        """
        Map parameter names to parameters in a SnglInspiral.
        """
        vals = None
        if self._triggers is not None:
            for key,value in self._injXMLFuncMap.items():
                if paramname.lower().strip() == key.lower().strip():
                    try:
                        vals = dict([(trig.ifo,self._injXMLFuncMap[key](trig)) for trig in self._triggers])
                    except AttributeError:
                        break
        return vals

    def __getitem__(self,key):
        """
        Container method . Returns posterior chain,one_d_pos, with name one_d_pos.name.
        """
        return self._posterior[key.lower()]

    def __len__(self):
        """
        Container method. Defined as number of samples.
        """
        return len(self._logL)

    def __iter__(self):
        """
        Container method. Returns iterator from self.forward for use in
        for (...) in (...) etc.
        """
        return self.forward()

    def forward(self):
        """
        Generate a forward iterator (in sense of list of names) over Posterior
        with name,one_d_pos.
        """
        current_item = 0
        while current_item < self.dim:
            name=self._posterior.keys()[current_item]
            pos=self._posterior[name]
            current_item += 1
            yield name,pos

    def bySample(self):
        """
        Generate a forward iterator over the list of samples corresponding to
        the data stored within the Posterior instance. These are returned as
        ParameterSamples instances.
        """
        current_item=0
        pos_array,header=self.samples
        while current_item < len(self):
            sample_array=(np.squeeze(pos_array[current_item,:]))
            yield PosteriorSample(sample_array, header, header)
            current_item += 1


    @property
    def dim(self):
        """
        Return number of parameters.
        """
        return len(self._posterior.keys())

    @property
    def names(self):
        """
        Return list of parameter names.
        """
        nameslist=[]
        for key,value in self:
            nameslist.append(key)
        return nameslist

    @property
    def means(self):
        """
        Return dict {paramName:paramMean} .
        """
        meansdict={}
        for name,pos in self:
            meansdict[name]=pos.mean
        return meansdict

    @property
    def medians(self):
        """
        Return dict {paramName:paramMedian} .
        """
        mediansdict={}
        for name,pos in self:
            mediansdict[name]=pos.median
        return mediansdict

    @property
    def stdevs(self):
        """
        Return dict {paramName:paramStandardDeviation} .
        """
        stdsdict={}
        for name,pos in self:
            stdsdict[name]=pos.stdev
        return stdsdict

    @property
    def name(self):
        """
        Return qualified string containing the 'name' of the Posterior instance.
        """
        return self.__name

    @property
    def description(self):
        """
        Return qualified string containing a 'description' of the Posterior instance.
        """
        return self.__description

    def append(self,one_d_posterior):
        """
        Container method. Add a new OneDParameter to the Posterior instance.
        """
        self._posterior[one_d_posterior.name]=one_d_posterior
        return

    def pop(self,param_name):
        """
        Container method.  Remove PosteriorOneDPDF from the Posterior instance.
        """
        return self._posterior.pop(param_name)

    def append_mapping(self, new_param_names, func, post_names):
        """
        Append posteriors pos1,pos2,...=func(post_names)
        """
        #1D input
        if isinstance(post_names, str):
            old_post = self[post_names]
            old_inj  = old_post.injval
            old_trigs  = old_post.trigvals
            if old_inj:
                new_inj = func(old_inj)
            else:
                new_inj = None
            if old_trigs:
                new_trigs = {}
                for IFO in old_trigs.keys():
                    new_trigs[IFO] = func(old_trigs[IFO])
            else:
                new_trigs = None

            samps = func(old_post.samples)
            new_post = PosteriorOneDPDF(new_param_names, samps, injected_value=new_inj, trigger_values=new_trigs)
            if new_post.samples.ndim is 0:
                print "WARNING: No posterior calculated for %s ..." % post.name
            else:
                self.append(new_post)
        #MultiD input
        else:
            old_posts = [self[post_name] for post_name in post_names]
            old_injs = [post.injval for post in old_posts]
            old_trigs = [post.trigvals for post in old_posts]
            samps = func(*[post.samples for post in old_posts])
            #1D output
            if isinstance(new_param_names, str):
                if None not in old_injs:
                    inj = func(*old_injs)
                else:
                    inj = None
                if None not in old_trigs:
                    new_trigs = {}
                    for IFO in old_trigs[0].keys():
                        oldvals = [param[IFO] for param in old_trigs]
                        new_trigs[IFO] = func(*oldvals)
                else:
                    new_trigs = None
                new_post = PosteriorOneDPDF(new_param_names, samps, injected_value=inj, trigger_values=new_trigs)
                self.append(new_post)
            #MultiD output
            else:
                if None not in old_injs:
                    injs = func(*old_injs)
                else:
                    injs = [None for name in new_param_names]
                if None not in old_trigs:
                    new_trigs = [{} for param in range(len(new_param_names))]
                    for IFO in old_trigs[0].keys():
                        oldvals = [param[IFO] for param in old_trigs]
                        newvals = func(*oldvals)
                        for param,newval in enumerate(newvals):
                            new_trigs[param][IFO] = newval
                else:
                    new_trigs = [None for param in range(len(new_param_names))]
                new_posts = [PosteriorOneDPDF(new_param_name,samp,injected_value=inj,trigger_values=new_trigs) for (new_param_name,samp,inj,new_trigs) in zip(new_param_names,samps,injs,new_trigs)]
                for post in new_posts: 
                    if post.samples.ndim is 0: 
                        print "WARNING: No posterior calculated for %s ..." % post.name
                    else:
                        self.append(post)
        return

    def _average_posterior(self, samples, post_name):
        """
        Returns the average value of the 'post_name' column of the
        given samples.
        """
        ap = 0.0
        for samp in samples:
            ap = ap + samp[post_name]
        return ap / len(samples)

    def _average_posterior_like_prior(self, samples, logl_name, prior_name, log_bias = 0):
        """
        Returns the average value of the posterior assuming that the
        'logl_name' column contains log(L) and the 'prior_name' column
        contains the prior (un-logged).
        """
        ap = 0.0
        for samp in samples:
            ap += np.exp(samp[logl_name]-log_bias)*samp[prior_name]
        return ap / len(samples)

    def _bias_factor(self):
        """
        Returns a sensible bias factor for the evidence so that
        integrals are representable as doubles.
        """
        return np.mean(self._logL)

    def di_evidence(self, boxing=64):
        """
        Returns the log of the direct-integration evidence for the
        posterior samples.
        """
        allowed_coord_names=["spin1", "spin2", "a1", "phi1", "theta1", "a2", "phi2", "theta2",
                             "iota", "psi", "ra", "dec",
                             "phi_orb", "phi0", "dist", "time", "mc", "mchirp", "chirpmass", "q"]
        samples,header=self.samples()
        header=header.split()
        coord_names=[name for name in allowed_coord_names if name in header]
        coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]
        tree=KDTree(coordinatized_samples)

        if "prior" in header and "logl" in header:
            bf = self._bias_factor()
            return bf + np.log(tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "logl", "prior", bf), boxing))
        elif "prior" in header and "likelihood" in header:
            bf = self._bias_factor()
            return bf + np.log(tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "likelihood", "prior", bf), boxing))
        elif "post" in header:
            return np.log(tree.integrate(lambda samps: self._average_posterior(samps, "post"), boxing))
        elif "posterior" in header:
            return np.log(tree.integrate(lambda samps: self._average_posterior(samps, "posterior"), boxing))
        else:
            raise RuntimeError("could not find 'post', 'posterior', 'logl' and 'prior', or 'likelihood' and 'prior' columns in output to compute direct integration evidence")

    def elliptical_subregion_evidence(self):
        """Returns an approximation to the log(evidence) obtained by
        fitting an ellipse around the highest-posterior samples and
        performing the harmonic mean approximation within the ellipse.
        Because the ellipse should be well-sampled, this provides a
        better approximation to the evidence than the full-domain HM."""
        allowed_coord_names=["spin1", "spin2", "a1", "phi1", "theta1", "a2", "phi2", "theta2",
                             "iota", "psi", "ra", "dec",
                             "phi_orb", "phi0", "dist", "time", "mc", "mchirp", "chirpmass", "q"]
        samples,header=self.samples()
        header=header.split()

        n=int(0.05*samples.shape[0])
        if not n > 1:
            raise IndexError

        coord_names=[name for name in allowed_coord_names if name in header]
        indexes=np.argsort(self._logL[:,0])

        my_samples=samples[indexes[-n:], :] # The highest posterior samples.
        my_samples=np.array([PosteriorSample(sample,header,coord_names).coord() for sample in my_samples])

        mu=np.mean(my_samples, axis=0)
        cov=np.cov(my_samples, rowvar=0)

        d0=None
        for mysample in my_samples:
            d=np.dot(mysample-mu, np.linalg.solve(cov, mysample-mu))
            if d0 is None:
                d0 = d
            else:
                d0=max(d0,d)

        ellipse_logl=[]
        ellipse_samples=[]
        for sample,logl in zip(samples, self._logL):
            coord=PosteriorSample(sample, header, coord_names).coord()
            d=np.dot(coord-mu, np.linalg.solve(cov, coord-mu))

            if d <= d0:
                ellipse_logl.append(logl)
                ellipse_samples.append(sample)
        
        if len(ellipse_samples) > 5*n:
            print 'WARNING: ellpise evidence region encloses significantly more samples than %d'%n

        ellipse_samples=np.array(ellipse_samples)
        ellipse_logl=np.array(ellipse_logl)

        ndim = len(coord_names)
        ellipse_volume=np.pi**(ndim/2.0)*d0**(ndim/2.0)/special.gamma(ndim/2.0+1)*np.sqrt(np.linalg.det(cov))

        try:
            prior_index=header.index('prior')
            pmu=np.mean(ellipse_samples[:,prior_index])
            pstd=np.std(ellipse_samples[:,prior_index])
            if pstd/pmu > 1.0:
                print 'WARNING: prior variation greater than 100\% over elliptical volume.'
            approx_prior_integral=ellipse_volume*pmu
        except KeyError:
            # Maybe prior = 1?
            approx_prior_integral=ellipse_volume

        ll_bias=np.mean(ellipse_logl)
        ellipse_logl = ellipse_logl - ll_bias

        return np.log(approx_prior_integral) - np.log(np.mean(1.0/np.exp(ellipse_logl))) + ll_bias

    def harmonic_mean_evidence(self):
        """
        Returns the log of the harmonic mean evidence for the set of
        posterior samples.
        """
        bf = self._bias_factor()
        return bf + np.log(1/np.mean(1/np.exp(self._logL-bf)))

    def _posMode(self):
        """
        Find the sample with maximum posterior probability. Returns value
        of posterior and index of sample .
        """
        pos_vals=self._logL
        max_i=0
        max_pos=pos_vals[0]
        for i in range(len(pos_vals)):
            if pos_vals[i] > max_pos:
                max_pos=pos_vals[i]
                max_i=i
        return max_pos,max_i

    def _print_table_row(self,name,entries):
        """
        Print a html table row representation of

        name:item1,item2,item3,...
        """

        row_str='<tr><td>%s</td>'%name
        for entry in entries:
            row_str+='<td>%s</td>'%entry
        row_str+='</tr>'
        return row_str

    @property
    def maxL(self):
        """
        Return the maximum posterior probability and the corresponding
        set of parameters.
        """
        maxLvals={}
        max_pos,max_i=self._posMode()
        for param_name in self.names:
            maxLvals[param_name]=self._posterior[param_name].samples[max_i][0]

        return (max_pos,maxLvals)

    def samples(self):
        """
        Return an (M,N) numpy.array of posterior samples; M = len(self);
        N = dim(self) .
        """
        header_string=''
        posterior_table=[]
        for param_name,one_pos in self:
            column=np.array(one_pos.samples)
            header_string+=param_name+' '
            posterior_table.append(column)
        posterior_table=tuple(posterior_table)
        return np.column_stack(posterior_table),header_string

    def write_to_file(self,fname):
        """
        Dump the posterior table to a file in the 'common format'.
        """
        column_list=()

        posterior_table,header_string=self.samples()

        fobj=open(fname,'w')

        fobj.write(header_string+'\n')
        np.savetxt(fobj,posterior_table)
        fobj.close()

        return

    def gelman_rubin(self, pname):
        """
        Returns an approximation to the Gelman-Rubin statistic (see
        Gelman, A. and Rubin, D. B., Statistical Science, Vol 7,
        No. 4, pp. 457--511 (1992)) for the parameter given, accurate
        as the number of samples in each chain goes to infinity.  The
        posterior samples must have a column named 'chain' so that the
        different chains can be separated.
        """
        if "chain" in self.names:
            chains=np.unique(self["chain"].samples)
            chain_index=self.names.index("chain")
            param_index=self.names.index(pname)
            data,header=self.samples()
            chainData=[data[ data[:,chain_index] == chain, param_index] for chain in chains]
            allData=np.concatenate(chainData)
            chainMeans=[np.mean(data) for data in chainData]
            chainVars=[np.var(data) for data in chainData]
            BoverN=np.var(chainMeans)
            W=np.mean(chainVars)
            sigmaHat2=W + BoverN
            m=len(chainData)
            VHat=sigmaHat2 + BoverN/m
            R = VHat/W
            return R
        else:
            raise RuntimeError('could not find necessary column header "chain" in posterior samples')

    def __str__(self):
        """
        Define a string representation of the Posterior class ; returns
        a html formatted table of various properties of posteriors.
        """
        return_val='<table border="1" id="statstable"><tr><th/>'

        column_names=['maxL','stdev','mean','median','stacc','injection value']
        IFOs = []
        if self._triggers is not None:
            IFOs = [trig.ifo for trig in self._triggers]
            for IFO in IFOs:
                column_names.append(IFO+' trigger values')

        for column_name in column_names:
            return_val+='<th>%s</th>'%column_name

        return_val+='</tr>'

        for name,oned_pos in self:

            max_pos,max_i=self._posMode()
            maxL=oned_pos.samples[max_i][0]
            mean=str(oned_pos.mean)
            stdev=str(oned_pos.stdev)
            median=str(np.squeeze(oned_pos.median))
            stacc=str(oned_pos.stacc)
            injval=str(oned_pos.injval)
            trigvals=oned_pos.trigvals

            row = [maxL,stdev,mean,median,stacc,injval]
            if self._triggers is not None:
                for IFO in IFOs:
                    try:
                        row.append(str(trigvals[IFO]))
                    except TypeError:
                        row.append(None)
            return_val+=self._print_table_row(name,row)

        return_val+='</table>'

        parser=XMLParser()
        parser.feed(return_val)
        Estr=parser.close()

        elem=Estr
        rough_string = tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return_val=reparsed.toprettyxml(indent="  ")

        return return_val

    def write_vot_info(self):
      """
      Writes the information stored in the VOTtree if there is one
      """
      target=VOT2HTML()
      parser=XMLParser(target=target)
      parser.feed(self._votfile)
      return parser.close()


class KDTree(object):
    """
    A kD-tree.
    """
    def __init__(self, objects):
        """
        Construct a kD-tree from a sequence of objects.  Each object
        should return its coordinates using obj.coord().
        """
        if len(objects) == 0:
            raise RuntimeError("cannot construct kD-tree out of zero objects---you may have a repeated sample in your list")
        elif len(objects) == 1:
            self._objects = objects[:]
            coord=self._objects[0].coord()
            self._bounds = coord,coord
        elif self._same_coords(objects):
            # All the same coordinates
            self._objects = [ objects[0] ]
            coord=self._objects[0].coord()
            self._bounds = coord,coord
        else:
            self._objects = objects[:]
            self._bounds = self._bounds_of_objects()
            low,high=self._bounds
            self._split_dim=self._longest_dimension()
            longest_dim = self._split_dim
            sorted_objects=sorted(self._objects, key=lambda obj: (obj.coord())[longest_dim])
            N = len(sorted_objects)
            bound=0.5*(sorted_objects[N/2].coord()[longest_dim] + sorted_objects[N/2-1].coord()[longest_dim])
            low = [obj for obj in self._objects if obj.coord()[longest_dim] < bound]
            high = [obj for obj in self._objects if obj.coord()[longest_dim] >= bound]
            if len(low)==0:
                # Then there must be multiple values with the same
                # coordinate as the minimum element of high
                low = [obj for obj in self._objects if obj.coord()[longest_dim]==bound]
                high = [obj for obj in self._objects if obj.coord()[longest_dim] > bound]
            self._left = KDTree(low)
            self._right = KDTree(high)

    def _same_coords(self, objects):
        """
        True if and only if all the given objects have the same
        coordinates.
        """
        if len(objects) <= 1:
            return True
        coords = [obj.coord() for obj in objects]
        c0 = coords[0]
        for ci in coords[1:]:
            if not np.all(ci == c0):
                return False
        return True

    def _bounds_of_objects(self):
        """
        Bounds of the objects contained in the tree.
        """
        low=self._objects[0].coord()
        high=self._objects[0].coord()
        for obj in self._objects[1:]:
            low=np.minimum(low,obj.coord())
            high=np.maximum(high,obj.coord())
        return low,high

    def _longest_dimension(self):
        """
        Longest dimension of the tree bounds.
        """
        low,high = self._bounds
        widths = high-low
        return np.argmax(widths)

    def objects(self):
        """
        Returns the objects in the tree.
        """
        return self._objects[:]

    def __iter__(self):
        """
        Iterator over all the objects contained in the tree.
        """
        return self._objects.__iter__()

    def left(self):
        """
        Returns the left tree.
        """
        return self._left

    def right(self):
        """
        Returns the right tree.
        """
        return self._right

    def split_dim(self):
        """
        Returns the dimension along which this level of the kD-tree
        splits.
        """
        return self._split_dim

    def bounds(self):
        """
        Returns the coordinates of the lower-left and upper-right
        corners of the bounding box for this tree: low_left, up_right
        """
        return self._bounds

    def volume(self):
        """
        Returns the volume of the bounding box of the tree.
        """
        v = 1.0
        low,high=self._bounds
        for l,h in zip(low,high):
            v = v*(h - l)
        return v

    def integrate(self,f,boxing=64):
        """
        Returns the integral of f(objects) over the tree.  The
        optional boxing parameter determines how deep to descend into
        the tree before computing f.
        """
        # if len(self._objects) <= boxing:
        #     return self.volume()*f(self._objects)
        # else:
        #     return self._left.integrate(f, boxing) + self._right.integrate(f, boxing)
        
        def x(tree):
            return tree.volume()*f(tree._objects)
            
        def y(a,b):
            return a+b
        
        return self.operate(x,y,boxing=boxing)
        
    def operate(self,f,g,boxing=64):
        """
        Operates on tree nodes exceeding boxing parameter depth.
        """
        if len(self._objects) <= boxing:
            return f(self)
        else:
            
            return g(self._left.operate(f,g,boxing),self._right.operate(f,g,boxing))


class KDTreeVolume(object):
    """
    A kD-tree suitable for splitting parameter spaces and counting hypervolumes.
    Is modified from the KDTree class so that bounding boxes are stored. This means that
    there are no longer gaps in the hypervolume once the samples have been split into groups.
    """
    def __init__(self, objects,boundingbox,dims=0):
        """
        Construct a kD-tree from a sequence of objects.  Each object
        should return its coordinates using obj.coord().
        the obj should also store the bounds of the hypervolume its found in.
        for non-leaf objects we need the name of the dimension split and value at split.
        """
        self._dimension = dims
        self._bounds = boundingbox
        self._weight = 1
        if len(objects) == 0: #for no objects - something is wrong, i think it can only happen in first call
            raise RuntimeError("cannot construct kD-tree out of zero objects---you may have a repeated sample in your list")
        elif len(objects) == 1: #1 object, have reached leaf of tree
            self._objects = objects[:]
        elif self._same_coords(objects): # When ALL samples have the same coordinates in all dimensions
            self._weight = len(objects)
            self._objects = [ objects[0] ] #need to modify kdtree_bin functions to use _weight to get correct number of samples
            coord=self._objects[0].coord()
        else: #construct next level of tree with multiple samples
            self._objects = objects[:]
            split_dim = self._dimension
            sorted_objects=sorted(self._objects, key=lambda obj: (obj.coord())[split_dim])
            N = len(sorted_objects)
            self._split_value = 0.5*(sorted_objects[N/2].coord()[split_dim] + sorted_objects[N/2-1].coord()[split_dim])
            bound = self._split_value
            low = [obj for obj in self._objects if obj.coord()[split_dim] < bound]
            high = [obj for obj in self._objects if obj.coord()[split_dim] >= bound]
            if len(low)==0:
                # Then there must be multiple values with the same
                # coordinate as the minimum element of 'high'
                low = [obj for obj in self._objects if obj.coord()[split_dim] == bound]
                high = [obj for obj in self._objects if obj.coord()[split_dim] > bound]
            leftBoundingbox = []
            rightBoundingbox = []
            for i in self._bounds:
                leftBoundingbox.append(list(i))
                rightBoundingbox.append(list(i))
            leftBoundingbox[1][split_dim] = bound
            rightBoundingbox[0][split_dim] = bound
            # designate the next dimension to use for split for sub-trees
            # if has got to the end of the list of dimensions then starts
            # again at dimension = 0
            if (split_dim < (len(self._objects[0].coord()) - 1)):
                child_dim = split_dim + 1
            else:
                child_dim = 0
            self._left = KDTreeVolume(low,leftBoundingbox,dims = child_dim)
            # added in a load of messing about incase there are repeated values in the currently checked dimension
            if (len(high) != 0):
                self._right = KDTreeVolume(high,rightBoundingbox,dims = child_dim)
            else:
                self._right = None

    def _same_coords(self, objects):
        """
        True if and only if all the given objects have the same
        coordinates.
        """
        if len(objects) <= 1:
            return True
        coords = [obj.coord() for obj in objects]
        c0 = coords[0]
        for ci in coords[1:]:
            if not np.all(ci == c0):
                return False
        return True

    def objects(self):
        """
        Returns the objects in the tree.
        """
        return self._objects[:]

    def __iter__(self):
        """
        Iterator over all the objects contained in the tree.
        """
        return self._objects.__iter__()

    def left(self):
        """
        Returns the left tree.
        """
        return self._left

    def right(self):
        """
        Returns the right tree.
        """
        return self._right

    def split_dim(self):
        """
        Returns the dimension along which this level of the kD-tree
        splits.
        """
        return self._split_dim

    def bounds(self):
        """
        Returns the coordinates of the lower-left and upper-right
        corners of the bounding box for this tree: low_left, up_right
        """
        return self._bounds

    def volume(self):
        """
        Returns the volume of the bounding box of the tree.
        """
        v = 1.0
        low,high=self._bounds
        for l,h in zip(low,high):
            v = v*(h - l)
        return v

    def integrate(self,f,boxing=64):
        """
        Returns the integral of f(objects) over the tree.  The
        optional boxing parameter determines how deep to descend into
        the tree before computing f.
        """
        def x(tree):
            return tree.volume()*f(tree._objects)

        def y(a,b):
            return a+b

        return self.operate(x,y,boxing=boxing)

    def operate(self,f,g,boxing=64):
        """
        Operates on tree nodes exceeding boxing parameter depth.
        """
        if len(self._objects) <= boxing:
            return f(self)
        else:
            return g(self._left.operate(f,g,boxing),self._right.operate(f,g,boxing))

    def search(self,coordinates,boxing = 64):
        """
        takes a set of coordinates and searches down through the tree untill it gets
        to a box with less than 'boxing' objects in it and returns the box bounds,
        number of objects in the box, and the weighting.
        """
        if len(self._objects) <= boxing:
            return self._bounds,len(self._objects),self._weight
        elif coordinates[self._dimension] < self._split_value:
            return self._left.search(coordinates,boxing)
        else:
            return self._right.search(coordinates,boxing)



class PosteriorSample(object):
    """
    A single parameter sample object, suitable for inclusion in a
    kD-tree.
    """

    def __init__(self, sample_array, headers, coord_names):
        """
        Given the sample array, headers for the values, and the names
        of the desired coordinates, construct a parameter sample
        object.
        """
        self._samples=sample_array[:]
        self._headers=headers
        if not (len(sample_array) == len(self._headers)):
            print "Header length = ", len(self._headers)
            print "Sample length = ", len(sample_array)
            raise RuntimeError("parameter and sample lengths do not agree")
        self._coord_names=coord_names
        self._coord_indexes=[self._headers.index(name) for name in coord_names]

    def __getitem__(self, key):
        """
        Return the element with the corresponding name.
        """
        key=key.lower()
        if key in self._headers:
            idx=self._headers.index(key)
            return self._samples[idx]
        else:
            raise KeyError("key not found in posterior sample: %s"%key)

    def coord(self):
        """
        Return the coordinates for the parameter sample.
        """
        return self._samples[self._coord_indexes]





class AnalyticLikelihood(object):
    """
    Return analytic likelihood values.
    """

    def __init__(self, covariance_matrix_files, mean_vector_files):
        """
        Prepare analytic likelihood for the given parameters.
        """
        # Make sure files names are in a list
        if isinstance(covariance_matrix_files, str):
            covariance_matrix_files = [covariance_matrix_files]
        if isinstance(mean_vector_files, str):
            mean_vector_files = [mean_vector_files]

        covarianceMatrices = [np.loadtxt(csvFile, delimiter=',') for csvFile in covariance_matrix_files]
        num_matrices = len(covarianceMatrices)

        if num_matrices != len(mean_vector_files):
            raise RuntimeError('Must give a mean vector list for every covariance matrix')

        param_line = open(mean_vector_files[0]).readline()
        self._params = [param.strip() for param in param_line.split(',')]

        converter=lambda x: eval(x.replace('pi','%.32f'%pi_constant))  # converts fractions w/ pi (e.g. 3.0*pi/2.0)
        self._modes = []
        for i in range(num_matrices):
            CM = covarianceMatrices[i]
            vecFile = mean_vector_files[i]

            param_line = open(vecFile).readline()
            params = [param.strip() for param in param_line.split(',')]
            if set(params)!=set(self._params):
                raise RuntimeError('Parameters do not agree between mean vector files.')

            sigmas = dict(zip(params,np.sqrt(CM.diagonal())))
            colNums = range(len(params))
            converters = dict(zip(colNums,[converter for i in colNums]))
            meanVectors = np.loadtxt(vecFile, delimiter=',', skiprows=1, converters=converters)
            try:
                for vec in meanVectors:
                    means = dict(zip(params,vec))
                    mode = [(param, stats.norm(loc=means[param],scale=sigmas[param])) for param in params]
                    self._modes.append(dict(mode))
            except TypeError:
                means = dict(zip(params,meanVectors))
                mode = [(param, stats.norm(loc=means[param],scale=sigmas[param])) for param in params]
                self._modes.append(dict(mode))
        self._num_modes = len(self._modes)

    def pdf(self, param):
        """
        Return PDF function for parameter.
        """
        pdf = None
        if param in self._params:
            pdf = lambda x: (1.0/self._num_modes) * sum([mode[param].pdf(x) for mode in self._modes])
        return pdf

    def cdf(self, param):
        """
        Return PDF function for parameter.
        """
        cdf = None
        if param in self._params:
            cdf = lambda x: (1.0/self._num_modes) * sum([mode[param].cdf(x) for mode in self._modes])
        return cdf

    @property
    def names(self):
        """
        Return list of parameter names described by analytic likelihood function.
        """
        return self._params



#===============================================================================
# Web page creation classes (wrap ElementTrees)
#===============================================================================

class htmlChunk(object):
    """
    A base class for representing web content using ElementTree .
    """
    def __init__(self,tag,attrib=None,parent=None):

        self._html=Element(tag)#attrib={'xmlns':"http://www.w3.org/1999/xhtml"})
        if attrib:
            for attribname,attribvalue in attrib.items():
                self._html.attrib[attribname]=attribvalue
        if parent:
            parent.append(self._html)


    def toprettyxml(self):
        """
        Return a pretty-printed XML string of the htmlPage.
        """
        elem=self._html
        rough_string = tostring(elem)
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    def __str__(self):
        return self.toprettyxml()

    def write(self,string):
        parser=XMLParser()
        parser.feed(string)
        Estr=parser.close()
        self._html.append(Estr)

    def p(self,pstring):
        Ep=Element('p')
        Ep.text=pstring
        self._html.append(Ep)
        return Ep

    def h1(self,h1string):
        Ep=Element('h1')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep
#
    def h5(self,h1string):
        Ep=Element('h5')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep

    def h2(self,h2string):
        Ep=Element('h2')
        Ep.text=h2string
        self._html.append(Ep)
        return Ep

    def h3(self,h1string):
        Ep=Element('h3')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep

    def br(self):
        Ebr=Element('br')
        self._html.append(Ebr)
        return Ebr

    def hr(self):
        Ehr=Element('hr')
        self._html.append(Ehr)
        return Ehr

    def a(self,url,linktext):
        Ea=Element('a')
        Ea.attrib['href']=url
        Ea.text=linktext
        self._html.append(Ea)
        return Ea

    def append(self,element):
        self._html.append(element)


#
class htmlPage(htmlChunk):
    """
    A concrete class for generating an XHTML(1) document. Inherits from htmlChunk.
    """
    def __init__(self,title=None,css=None,toc=False):
        htmlChunk.__init__(self,'html',attrib={'xmlns':"http://www.w3.org/1999/xhtml"})
        self.doctype_str='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'

        self._head=SubElement(self._html,'head')
        Etitle=SubElement(self._head,'title')
        self._body=SubElement(self._html,'body')
        self._css=None
        if title is not None:
            Etitle.text=str(title)
            self._title=SubElement(self._body,'h1')
            self._title.text=title

        if css is not None:
            self._css=SubElement(self._head,'style')
            self._css.attrib['type']="text/css"
            self._css.text=str(css)

    def __str__(self):
        return self.doctype_str+'\n'+self.toprettyxml()

    def add_section(self,section_name):
        newSection=htmlSection(section_name)
        self._body.append(newSection._html)
        return newSection

    @property
    def body():
        return self._body

    @property
    def head():
        return self._head


class htmlSection(htmlChunk):
    """
    Represents a block of html fitting within a htmlPage. Inherits from htmlChunk.
    """
    def __init__(self,section_name,htmlElement=None):
        htmlChunk.__init__(self,'div',attrib={'class':'ppsection'},parent=htmlElement)

        self.h2(section_name)

#===============================================================================
# Internal module functions
#===============================================================================

def _skyhist_cart_slow(skycarts,sky_samples):
    """
    @deprecated: This is a pure python version of the C extension function
        pylal._bayespputils._skyhist_cart .
    """

    N=len(skycarts)
    print 'operating on %d sky points'%(N)
    bins=np.zeros(N)
    for RAsample,decsample in sky_samples:
        sampcart=pol2cart(RAsample,decsample)
        maxdx=-1
        maxvalue=-1
        for i in xrange(0,N):
            dx=np.dot(sampcart,skycarts[i])
            if dx>maxvalue:
                    maxdx=i
                    maxvalue=dx

        bins[maxdx]+=1
    return bins
#
def _sky_hist(skypoints,samples):
    """
    @deprecated: This is an old pure python version of the C extension function
        pylal._bayespputils._skyhist_cart .
    """
    N=len(skypoints)
    print 'operating on %d sky points' % (N)
    bins=zeros(N)
    j=0
    for sample in samples:
        seps=map(lambda s: ang_dist(sample[RAdim],sample[decdim],s[1],s[0]),skypoints)
        minsep=math.pi
        for i in range(0,N):
            if seps[i]<minsep:
                minsep=seps[i]
                mindx=i
        bins[mindx]=bins[mindx]+1
        j=j+1
        print 'Done %d/%d iterations, minsep=%f degrees'\
            %(j,len(samples),minsep*(180.0/3.1415926))
    return (skypoints,bins)
#

def _calculate_confidence_levels(hist, points, injBin, NSamples):
    """
    Returns (injectionconf, toppoints), where injectionconf is the
    confidence level of the injection, contained in the injBin and
    toppoints is a list of (pointx, pointy, ptindex, frac), with
    pointx and pointy the (x,y) coordinates of the corresponding
    element of the points array, ptindex the index of the point in the
    array, and frac the cumulative fraction of points with larger
    posterior probability.

    The hist argument should be a one-dimensional array that contains
    counts of sample points in each bin.

    The points argument should be a 2-D array storing the sky location
    associated with each bin; the first index runs from 0 to NBins -
    1, while the second index runs from 0 to 1.

    The injBin argument gives the bin index in which the injection is
    found.

    The NSamples argument is used to normalize the histogram counts
    into fractional probability.
    """

    histIndices=np.argsort(hist)[::-1]  # In decreasing order

    toppoints=[]
    frac=0.0
    injConf=None
    for i in histIndices:
        frac+=float(hist[i])/float(NSamples)
        toppoints.append((points[i,0], points[i,1], i, frac))
        if i == injBin:
            injConf=frac
            print 'Injection found at confidence level %g'%injConf

    return (injConf, toppoints)

def _greedy_bin(greedyHist,greedyPoints,injection_bin_index,bin_size,Nsamples,confidence_levels):
    """
    An interal function representing the common, dimensionally-independent part of the
    greedy binning algorithms.
    """

    #Now call confidence level C extension function to determine top-ranked pixels
    (injectionconfidence,toppoints)=_calculate_confidence_levels(greedyHist, greedyPoints, injection_bin_index, Nsamples)

    #Determine interval/area contained within given confidence intervals
    nBins=0
    confidence_levels.sort()
    reses={}
    toppoints=np.array(toppoints)
    for printcl in confidence_levels:
        nBins=np.searchsorted(toppoints[:,3], printcl) + 1

        if nBins >= len(toppoints):
            nBins=len(toppoints)-1

        accl=toppoints[nBins-1,3]

        reses[printcl]=nBins*bin_size

    #Find area
    injection_area=None
    if injection_bin_index and injectionconfidence:
        i=list(np.nonzero(np.asarray(toppoints)[:,2]==injection_bin_index))[0]
        injection_area=bin_size*(i+1)

    return toppoints,injectionconfidence,reses,injection_area
#


#
#===============================================================================
# Public module functions
#===============================================================================

def kdtree_bin_sky_volume(posterior,confidence_levels):
    
    confidence_levels.sort()
    
    class Harvester(list):
        
        def __init__(self):
            list.__init__(self)
            self.unrho=0.
            
        def __call__(self,tree):
            number_density=float(len(tree.objects()))/float(tree.volume())
            self.append([number_density,tree.volume(),tree.bounds()])
            self.unrho+=number_density
            
        def close_ranks(self):
            
            for i in range(len(self)):
                self[i][0]/=self.unrho
            
            return sorted(self,key=itemgetter(0))
    
    def h(a,b):
        pass
    
    samples,header=posterior.samples()
    header=header.split()
    coord_names=["ra","dec","dist"]
    coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]
    tree=KDTree(coordinatized_samples)
    
    a=Harvester()
    samples_per_bin=10
    tree.operate(a,h,boxing=samples_per_bin)
    
    b=a.close_ranks()
    b.reverse()
    
    acc_rho=0.
    acc_vol=0.
    cl_idx=0
    confidence_intervals={}
    for rho,vol,bounds in b:
        acc_rho+=rho
        acc_vol+=vol
    
        if acc_rho>confidence_levels[cl_idx]:
            confidence_intervals[acc_rho]=acc_vol
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break
    
    return confidence_intervals


def kdtree_bin_sky_area(posterior,confidence_levels,samples_per_bin=10):
    """
    takes samples and applies a KDTree to them to return confidence levels
    returns confidence_intervals - dictionary of user_provided_CL:calculated_area
            b - ordered list of KD leaves
            injInfo - if injection values provided then returns
                      [Bounds_of_inj_kd_leaf ,number_samples_in_box, weight_of_box,injection_CL ,injection_CL_area]
    Not quite sure that the repeated samples case is fixed, posibility of infinite loop.
    """
    confidence_levels.sort()
    from math import cos, pi
    class Harvester(list):
        """
        when called by kdtree.operate will be used to calculate the density of each bin (sky area)
        """
        def __init__(self):
            list.__init__(self)
            self.unrho=0.

        def __call__(self,tree):
            #area = (cos(tree.bounds()[0][1])-cos(tree.bounds()[1][1]))*(tree.bounds()[1][0] - tree.bounds()[0][0])
            area = - (cos(pi/2. - tree.bounds()[0][1])-cos(pi/2. - tree.bounds()[1][1]))*(tree.bounds()[1][0] - tree.bounds()[0][0])
            number_density=float(len(tree.objects()))/float(area)
            self.append([number_density,len(tree.objects()),area,tree.bounds()])
            self.unrho+=number_density

        def close_ranks(self):

            for i in range(len(self)):
                self[i][0]/=self.unrho

            return sorted(self,key=itemgetter(0))

    def h(a,b):
        pass

    peparser=PEOutputParser('common')

    samples,header=posterior.samples()
    header=header.split()
    coord_names=["ra","dec"]
    initial_dimensions = [[0.,-pi/2.],[2.*pi, pi/2.]]
    coordinatized_samples=[ParameterSample(row, header, coord_names) for row in samples]
    tree=KDTreeVolume(coordinatized_samples,initial_dimensions)

    a=Harvester()
    tree.operate(a,h,boxing=samples_per_bin)
    totalSamples = len(tree.objects())
    b=a.close_ranks()
    b.reverse()
    samplecounter=0.0
    for entry in b:
        samplecounter += entry[1]
        entry[1] = float(samplecounter)/float(totalSamples)

    acc_rho=0.
    acc_vol=0.
    cl_idx=0

    #checks for injection and extract details of the node in the tree that the injection is found
    if posterior['ra'].injval is not None and posterior['dec'].injval is not None:
        injBound,injNum,injWeight = tree.search([posterior['ra'].injval,posterior['dec'].injval],boxing = samples_per_bin)
        injInfo = [injBound,injNum,injWeight]
        inj_area = - (cos(pi/2. - injBound[0][1])-cos(pi/2. - injBound[1][1]))*(injBound[1][0] - injBound[0][0])
        inj_number_density=float(injNum)/float(inj_area)
        inj_rho = inj_number_density / a.unrho
    else:
        injInfo = None
        inj_area = None
        inj_number_density=None
        inj_rho = None

    #finds the volume contained within the confidence levels requested by user
    confidence_intervals={}
    for rho,confidence_level,vol,bounds in b:
        acc_vol+=vol

        if confidence_level>confidence_levels[cl_idx]:
            print str(confidence_level)
            print acc_vol
            confidence_intervals[confidence_levels[cl_idx]]=acc_vol
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break

    acc_vol = 0.
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol
    print 'total area: ' + str(acc_vol)

    #finds the confidence level of the injection and the volume of the associated contained region
    inj_confidence = None
    inj_confidence_area = None
    if inj_rho is not None:
        acc_vol=0.
        for rho,confidence_level,vol,bounds in b:
            acc_vol+=vol
            if rho <= inj_rho:
                inj_confidence = confidence_level
                inj_confidence_area = acc_vol
                injInfo.append(inj_confidence)
                injInfo.append(inj_confidence_area)
                print 'inj ' +str(vol)
                break
    return confidence_intervals, b, injInfo

def kdtree_bin(posterior,coord_names,confidence_levels,initial_boundingbox = None,samples_per_bin = 10):
    """
    takes samples and applies a KDTree to them to return confidence levels
    returns confidence_intervals - dictionary of user_provided_CL:calculated_volume
            b - ordered list of KD leaves
            initial_boundingbox - list of lists [upperleft_coords,lowerright_coords]
            injInfo - if injection values provided then returns
                      [Bounds_of_inj_kd_leaf ,number_samples_in_box, weight_of_box,injection_CL ,injection_CL_volume]
    Not quite sure that the repeated samples case is fixed, posibility of infinite loop.
    """
    confidence_levels.sort()
    print confidence_levels
    class Harvester(list):
        """
        when called by kdtree.operate will be used to calculate the density of each bin
        """
        def __init__(self):
            list.__init__(self)
            self.unrho=0.

        def __call__(self,tree):
            number_density=float(len(tree.objects()))/float(tree.volume())
            self.append([number_density,len(tree.objects()),tree.volume(),tree.bounds()])
            self.unrho+=number_density

        def close_ranks(self):

            for i in range(len(self)):
                self[i][0]/=self.unrho

            return sorted(self,key=itemgetter(0))

    def h(a,b):
        pass

    peparser=PEOutputParser('common')

    samples,header=posterior.samples()
    header=header.split()
    coordinatized_samples=[ParameterSample(row, header, coord_names) for row in samples]

    #if initial bounding box is not provided, create it using max/min of sample coords.
    if initial_boundingbox is None:
        low=coordinatized_samples[0].coord()
        high=coordinatized_samples[0].coord()
        for obj in coordinatized_samples[1:]:
            low=np.minimum(low,obj.coord())
            high=np.maximum(high,obj.coord())
        initial_boundingbox = [low,high]

    tree=KDTreeVolume(coordinatized_samples,initial_boundingbox)

    a=Harvester()
    tree.operate(a,h,boxing=samples_per_bin)

    b=a.close_ranks()
    b.reverse()
    totalSamples = len(tree.objects())
    samplecounter=0.0
    for entry in b:
        samplecounter += entry[1]
        entry[1] = float(samplecounter)/float(totalSamples)

    acc_rho=0.
    acc_vol=0.
    cl_idx=0

    #check that there is an injection value for all dimension names
    def checkNone(listoParams):
        for param in listoParams:
            if posterior[param].injval is None:
                return False
        return True

    #checks for injection and extract details of the lnode in the tree that the injection is found
    if checkNone(coord_names):
        injBound,injNum,injWeight = tree.search([posterior[x].injval for x in coord_names],boxing = samples_per_bin)
        injInfo = [injBound,injNum,injWeight]
        #calculate volume of injections bin
        inj_volume = 1.
        low = injBound[1]
        high = injBound[0]
        for aCoord,bCoord in zip(low,high):
            inj_volume = inj_volume*(bCoord - aCoord)
        inj_number_density=float(injNum)/float(inj_volume)
        inj_rho = inj_number_density / a.unrho
    else:
        injInfo = None
        inj_area = None
        inj_number_density=None
        inj_rho = None

    #finds the volume contained within the confidence levels requested by user
    confidence_intervals={}
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol

        if sample_number>confidence_levels[cl_idx]:
            confidence_intervals[confidence_levels[cl_idx]]=(acc_vol,sample_number)
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break

    acc_vol = 0.
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol

    #finds the confidence level of the injection and the volume of the associated contained region
    inj_confidence = None
    inj_confidence_area = None
    if inj_rho is not None:
        print 'calculating cl'
        acc_vol=0.
        for rho,confidence_level,vol,bounds in b:
            acc_vol+=vol
            if rho <= inj_rho:
                inj_confidence = confidence_level
                inj_confidence_area = acc_vol
                injInfo.append(inj_confidence)
                injInfo.append(inj_confidence_area)
                break

    return confidence_intervals, b, initial_boundingbox,injInfo

def greedy_bin_two_param(posterior,greedy2Params,confidence_levels):
    """
    Determine the 2-parameter Bayesian Confidence Intervals using a greedy
    binning algorithm.

    @param posterior: an instance of the Posterior class.

    @param greedy2Params: a dict - {param1Name:param1binSize,param2Name:param2binSize} .

    @param confidence_levels: A list of floats of the required confidence intervals [(0-1)].
    """

    #Extract parameter names
    par1_name,par2_name=greedy2Params.keys()

    #Set posterior array columns
    par1pos=posterior[par1_name.lower()].samples
    par2pos=posterior[par2_name.lower()].samples

    #Extract bin sizes
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]

    #Extract injection information
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval

    #Create 2D bin array
    par1pos_min=min(par1pos)[0]
    par2pos_min=min(par2pos)[0]

    par1pos_max=max(par1pos)[0]
    par2pos_max=max(par2pos)[0]

    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1

    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

    greedyHist = np.zeros(par1pos_Nbins*par2pos_Nbins,dtype='i8')
    greedyPoints = np.zeros((par1pos_Nbins*par2pos_Nbins,2))

    #Fill bin values
    par1_point=par1pos_min
    par2_point=par2pos_min
    for i in range(par2pos_Nbins):

        par1_point=par1pos_min
        for j in range(par1pos_Nbins):

            greedyPoints[j+par1pos_Nbins*i,0]=par1_point
            greedyPoints[j+par1pos_Nbins*i,1]=par2_point
            par1_point+=par1_bin
        par2_point+=par2_bin


    #If injection point given find which bin its in...
    injbin=None
    if par1_injvalue is not None and par2_injvalue is not None:

        par1_binNumber=int(floor((par1_injvalue-par1pos_min)/par1_bin))
        par2_binNumber=int(floor((par2_injvalue-par2pos_min)/par2_bin))

        injbin=int(par1_binNumber+par2_binNumber*par1pos_Nbins)
    elif par1_injvalue is None and par2_injvalue is not None:
        print "Injection value not found for %s!"%par1_name

    elif par1_injvalue is not None and par2_injvalue is None:
        print "Injection value not found for %s!"%par2_name

    #Bin posterior samples
    for par1_samp,par2_samp in zip(par1pos,par2pos):
        par1_samp=par1_samp[0]
        par2_samp=par2_samp[0]
        par1_binNumber=int(floor((par1_samp-par1pos_min)/par1_bin))
        par2_binNumber=int(floor((par2_samp-par2pos_min)/par2_bin))
        try:
            greedyHist[par1_binNumber+par2_binNumber*par1pos_Nbins]+=1
        except:
            raise RuntimeError("Problem binning samples: %i,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f .")%(par1_binNumber,par2_binNumber,par1pos_Nbins,par2pos_Nbins,par1_binNumber+par2_binNumber*par1pos_Nbins,par1_samp,par1pos_min,par1_bin,par1_samp,par2pos_min,par2_bin)
    #Call greedy bins routine
    toppoints,injection_cl,reses,injection_area=\
                                _greedy_bin(
                                                greedyHist,
                                                greedyPoints,
                                                injbin,
                                                float(par1_bin*par2_bin),
                                                int(len(par1pos)),
                                                confidence_levels
                                            )

    return toppoints,injection_cl,reses,injection_area

def pol2cart(long,lat):
    """
    Utility function to convert longitude,latitude on a unit sphere to
    cartesian co-ordinates.
    """

    x=np.cos(lat)*np.cos(long)
    y=np.cos(lat)*np.sin(long)
    z=np.sin(lat)
    return np.array([x,y,z])
#

def sph2cart(r,theta,phi):
    """
    Utiltiy function to convert r,theta,phi to cartesian co-ordinates.
    """
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z


def cart2sph(x,y,z):
    """
    Utility function to convert cartesian coords to r,theta,phi.
    """
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    phi = np.fmod(2*pi_constant + np.arctan2(y,x), 2*pi_constant)
    
    return r,theta,phi


def greedy_bin_sky(posterior,skyres,confidence_levels):
    """
    Greedy bins the sky posterior samples into a grid on the sky constructed so that
    sky boxes have roughly equal size (determined by skyres).

    @param posterior: Posterior class instance containing ra and dec samples.

    @param skyres: Desired approximate size of sky pixel on one side.

    @param confidence_levels: List of desired confidence levels [(0-1)].
    """

    from pylal import skylocutils

    np.seterr(under='ignore')

    if 'ra' in posterior.names and 'dec' in posterior.names:
        skypos=np.column_stack([posterior['ra'].samples,posterior['dec'].samples])
        raname='ra'
        decname='dec'
    elif 'rightascension' in posterior.names and 'declination' in posterior.names:
        skypos=np.column_stack([posterior['rightascension'].samples,posterior['declination'].samples])
        raname='rightascension'
        decname='declination'
    else:
        raise RuntimeError('could not find ra and dec or rightascention and declination in column names for sky position')

    injvalues=None

    sky_injpoint=(posterior[raname].injval,posterior[decname].injval)

    skypoints=np.array(skylocutils.gridsky(float(skyres)))
    skycarts=map(lambda s: pol2cart(s[1],s[0]),skypoints)
    skyinjectionconfidence=None

    shist=_skyhist_cart(np.array(skycarts),skypos)

    #shist=skyhist_cart(skycarts,list(pos))
    bins=skycarts

    # Find the bin of the injection if available
    injbin=None
    if None not in sky_injpoint:
        injhist=_skyhist_cart_slow(skycarts,np.array([sky_injpoint]))
        injbin=injhist.tolist().index(1)
        print 'Found injection in bin %d with co-ordinates %f,%f .'%(
                                                                     injbin,
                                                                     skypoints[injbin,0],
                                                                     skypoints[injbin,1]
                                                                     )

    return _greedy_bin(shist,skypoints,injbin,float(skyres)*float(skyres),len(skypos),confidence_levels)


def plot_sky_map(inj_pos,top_ranked_pixels,outdir):
    """
    Plots a sky map using the Mollweide projection in the Basemap package.

    @inj_pos: injected position in the sky in the form [dec,ra]

    @param top_ranked_pixels: the top-ranked sky pixels as determined by greedy_bin_sky.

    @param outdir: Output directory in which to save skymap.png image.
    """
    from mpl_toolkits.basemap import Basemap

    np.seterr(under='ignore')

    myfig=plt.figure()
    plt.clf()
    m=Basemap(projection='moll',lon_0=180.0,lat_0=0.0)
    
    # Plot an X on the injected position
    if (inj_pos is not None and inj_pos[1] is not None and inj_pos[0] is not None):
        ra_inj_rev=2*pi_constant - inj_pos[1]*57.296
        inj_plx,inj_ply=m(ra_inj_rev, inj_pos[0]*57.296)
        plt.plot(inj_plx,inj_ply,'wx',linewidth=12, markersize=22,mew=2,alpha=0.6)
    
    ra_reverse = 2*pi_constant - np.asarray(top_ranked_pixels)[::-1,1]*57.296

    plx,ply=m(
              ra_reverse,
              np.asarray(top_ranked_pixels)[::-1,0]*57.296
              )

    cnlevel=[1-tp for tp in np.asarray(top_ranked_pixels)[::-1,3]]
    plt.scatter(plx,ply,s=5,c=cnlevel,faceted=False,cmap=mpl_cm.jet)
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,45.),labels=[1,0,0,0],labelstyle='+/-')
    # draw parallels
    m.drawmeridians(np.arange(0.,360.,90.),labels=[0,0,0,1],labelstyle='+/-')
    # draw meridians
    plt.title("Skymap") # add a title
    plt.colorbar()
    myfig.savefig(os.path.join(outdir,'skymap.png'))
    plt.clf()

    #Save skypoints
    
    fid = open( os.path.join(outdir,'ranked_sky_pixels.dat'), 'w' ) 
    fid.write( 'dec(deg.)\tra(h.)\tprob.\tcumul.\n' ) 
    np.savetxt(
               fid,
               #os.path.join(outdir,'ranked_sky_pixels.dat'),
               np.column_stack(
                               [
                                np.asarray(top_ranked_pixels)[:,0]*57.296,
                                np.asarray(top_ranked_pixels)[:,1]*3.820,
                                np.append(np.asarray(top_ranked_pixels)[0,3],np.asarray(top_ranked_pixels)[1:,3]-np.asarray(top_ranked_pixels)[:-1,3]),
                                np.asarray(top_ranked_pixels)[:,3]
                                ]
                               ),
               fmt='%.4f',
               delimiter='\t'
               )
    fid.close() 

    return myfig
#

def mc2ms(mc,eta):
    """
    Utility function for converting mchirp,eta to component masses. The
    masses are defined so that m1>m2. The rvalue is a tuple (m1,m2).
    """
    root = np.sqrt(0.25-eta)
    fraction = (0.5+root) / (0.5-root)
    invfraction = 1/fraction

    m2= mc * np.power((1+fraction),0.2) / np.power(fraction,0.6)

    m1= mc* np.power(1+invfraction,0.2) / np.power(invfraction,0.6)
    return (m1,m2)
#
#

def q2ms(mc,q):
    """
    Utility function for converting mchirp,q to component masses. The
    masses are defined so that m1>m2. The rvalue is a tuple (m1,m2).
    """
    factor = mc * np.power(1+q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, 2.0/5.0);
    return (m1,m2)
#
#

def q2eta(mc,q):
    """
    Utility function for converting mchirp,q to eta. The
    rvalue is eta.
    """
    m1,m2 = q2ms(mc,q)
    eta = m1*m2/( (m1+m2)*(m1+m2) )
    return eta
#
#

def mc2q(mc,eta):
    """
    Utility function for converting mchirp,eta to new mass ratio q (m2/m1).
    """
    m1,m2 = mc2ms(mc,eta)
    q = m2/m1
    return q
#
#

def ang_dist(long1,lat1,long2,lat2):
    """
    Find the angular separation of (long1,lat1) and (long2,lat2), which are
        specified in radians.
    """

    x1=np.cos(lat1)*np.cos(long1)
    y1=np.cos(lat1)*np.sin(long1)
    z1=np.sin(lat1)
    x2=np.cos(lat2)*np.cos(long2)
    y2=np.cos(lat2)*np.sin(long2)
    z2=np.sin(lat2)
    sep=math.acos(x1*x2+y1*y2+z1*z2)
    return(sep)
#
#

def array_dot(vec1, vec2):
    """
    Calculate dot products between vectors in rows of numpy arrays.
    """
    if vec1.ndim==1:
        product = (vec1*vec2).sum()
    else:
        product = (vec1*vec2).sum(axis=1).reshape(-1,1)
    return product
#
#

def array_ang_sep(vec1, vec2):
    """
    Find angles between vectors in rows of numpy arrays.
    """
    vec1_mag = np.sqrt(array_dot(vec1, vec1))
    vec2_mag = np.sqrt(array_dot(vec2, vec2))
    return np.arccos(array_dot(vec1, vec2)/(vec1_mag*vec2_mag))
#
#

def array_polar_ang(vec):
    """
    Find polar angles of vectors in rows of a numpy array.
    """
    if vec.ndim==1:
        z = vec[2]
    else:
        z = vec[:,2].reshape(-1,1)
    norm = np.sqrt(array_dot(vec,vec))
    return np.arccos(z/norm)
#
#

def rotation_matrix(angle, direction):
    """
    Compute general rotation matrices for a given angles and direction vectors.
    """
    cosa = np.cos(angle)
    sina = np.sin(angle)
    direction /= np.sqrt(array_dot(direction,direction))
    #Assume calculating array of rotation matrices.
    try:
        nSamps = len(angle)
        R = np.array( [np.diag([i,i,i]) for i in cosa.flat] )
        R += np.array( [np.outer(direction[i],direction[i])*(1.0-cosa[i]) for i in range(nSamps)] )
        R += np.array( [np.array(   [[ 0.0,            -direction[i,2],    direction[i,1]],
                                     [ direction[i,2],  0.0,              -direction[i,0]],
                                     [-direction[i,1],  direction[i,0],    0.0          ]] ) * sina[i] for i in range(nSamps)] )
    #Only computing one rotation matrix.
    except TypeError:
        R = np.diag([cosa,cosa,cosa])
        R += np.outer(direction,direction) * (1.0 - cosa)
        R += np.array(   [[ 0.0,            -direction[2],    direction[1]],
                          [ direction[2],  0.0,              -direction[0]],
                          [-direction[1],  direction[0],    0.0          ]] ) * sina
    return R
#
#

def rotate_vector(R, vec):
    """
    Rotate vectors using the given rotation matrices.
    """
    if vec.ndim == 1:
        newVec = np.dot(R[i],vec[i])
    else:
        newVec = np.array( [np.dot(R[i],vec[i]) for i in range(len(vec))] )
    return newVec
#
#

def orbital_momentum(f_lower, mc, inclination):
    """
    Calculate orbital angular momentum vector.
    """
    Lmag = np.power(mc, 5.0/3.0) / np.power(pi_constant * lalconstants.LAL_MTSUN_SI * f_lower, 1.0/3.0)
    Lx, Ly, Lz = sph2cart(Lmag, inclination, 0.0)
    return np.hstack((Lx,Ly,Lz))
#
#

def component_momentum(m, a, theta, phi):
    """
    Calculate BH angular momentum vector.
    """
    Sx, Sy, Sz = sph2cart(m**2 * a, theta, phi)
    return np.hstack((Sx,Sy,Sz))
#
#

def spin_angles(f_lower,mc,eta,incl,a1,theta1,phi1,a2=None,theta2=None,phi2=None):
    """
    Calculate physical spin angles.
    """
    singleSpin = None in (a2,theta2,phi2)
    m1, m2 = mc2ms(mc,eta)
    L  = orbital_momentum(f_lower, mc, incl)
    S1 = component_momentum(m1, a1, theta1, phi1)
    if not singleSpin:
        S2 = component_momentum(m2, a2, theta2, phi2)
    else:
        S2 = 0.0
    J = L + S1 + S2
    tilt1 = array_ang_sep(L,S1)
    if not singleSpin:
        tilt2 = array_ang_sep(L,S2)
    else:
        tilt2 = None
    thetas = array_polar_ang(J)
    beta  = array_ang_sep(J,L)
    return tilt1, tilt2, thetas, beta
#
#

def plot_one_param_pdf_kde(fig,onedpos):

    from scipy import seterr as sp_seterr

    np.seterr(under='ignore')
    sp_seterr(under='ignore')
    pos_samps=onedpos.samples
    try:
        gkde=onedpos.gaussian_kde
    except np.linalg.linalg.LinAlgError:
        print 'Singular matrix in KDE. Skipping'
    else:
        ind=np.linspace(np.min(pos_samps),np.max(pos_samps),101)
        kdepdf=gkde.evaluate(ind)
        plt.plot(ind,kdepdf,color='green')
    return

def plot_one_param_pdf_line_hist(fig,pos_samps):
    plt.hist(pos_samps,kdepdf)


def plot_one_param_pdf(posterior,plot1DParams,analyticPDF=None,analyticCDF=None,plotkde=False):
    """
    Plots a 1D histogram and (gaussian) kernel density estimate of the
    distribution of posterior samples for a given parameter.

    @param posterior: an instance of the Posterior class.

    @param plot1DParams: a dict; {paramName:Nbins}

    @param analyticPDF: an analytic probability distribution function describing the distribution.

    @param analyticCDF: an analytic cumulative distribution function describing the distribution.

    """
    
    # matplotlib.rcParams['text.usetex']=True
    
    param=plot1DParams.keys()[0].lower()
    histbins=plot1DParams.values()[0]

    pos_samps=posterior[param].samples
    injpar=posterior[param].injval
    trigvals=posterior[param].trigvals

    myfig=plt.figure(figsize=(4,3.5),dpi=200)
    axes=plt.Axes(myfig,[0.2, 0.2, 0.7,0.7])
    myfig.add_axes(axes)
    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.6g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.6g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)
    
    if param.find('time')!=-1:
      offset=floor(min(pos_samps))
      pos_samps=pos_samps-offset
      ax1_name=param+' + %i'%(int(offset))
    else: ax1_name=param

    (n, bins, patches)=plt.hist(pos_samps,histbins,normed='true',facecolor='grey')
    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
        Nticks=3
    elif Nchars>5:
        Nticks=4
    elif Nchars>4:
        Nticks=6
    else:
        Nticks=6
    locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks)
    xmin,xmax=plt.xlim()
    if param=='rightascension' or param=='ra':
        locatorX=RALocator(min=xmin,max=xmax)
        majorFormatterX=RAFormatter()
    if param=='declination' or param=='dec':
        locatorX=DecLocator(min=xmin,max=xmax)
        majorFormatterX=DecFormatter()
    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)

    locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)
    if plotkde:  plot_one_param_pdf_kde(myfig,posterior[param])
    histbinSize=bins[1]-bins[0]
    if analyticPDF:
        (xmin,xmax)=plt.xlim()
        x = np.linspace(xmin,xmax,2*len(bins))
        plt.plot(x, analyticPDF(x), color='r', linewidth=2, linestyle='dashed')
        if analyticCDF:
            D,p = stats.kstest(pos_samps.flatten(), analyticCDF)
            plt.title("%s: ks p-value %.3f"%(param,p))

    rbins=None

    if injpar is not None:
        if min(pos_samps)<injpar and max(pos_samps)>injpar:
            plt.axvline(injpar, color='b', linestyle='-.')

            #rkde=gkde.integrate_box_1d(min(pos[:,i]),getinjpar(injection,i))
            #print "r of injected value of %s (kde) = %f"%(param,rkde)

            #Find which bin the true value is in
            bins_to_inj=(injpar-bins[0])/histbinSize
            injbinh=int(floor(bins_to_inj))
            injbin_frac=bins_to_inj-float(injbinh)

            #Integrate over the bins
            rbins=(sum(n[0:injbinh-1])+injbin_frac*n[injbinh])*histbinSize

    if trigvals is not None:
        for IFO in [IFO for IFO in trigvals.keys()]:
            trigval = trigvals[IFO]
            if min(pos_samps)<trigval and max(pos_samps)>trigval:
                if IFO=='H1': color = 'r'
                elif IFO=='L1': color = 'g'
                elif IFO=='V1': color = 'm'
                else: color = 'c'
                plt.axvline(trigval, color=color, linestyle='-.')
    #
    plt.grid()
    plt.xlabel(ax1_name)
    plt.ylabel('Probability Density')

    # For RA and dec set custom labels and for RA reverse
    if(param.lower()=='ra' or param.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        plt.xlim(xmax,xmin)
    #if(param.lower()=='ra' or param.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(param.lower()=='dec' or param.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return rbins,myfig#,rkde
#

class RALocator(matplotlib.ticker.MultipleLocator):
    """
    RA tick locations with some intelligence
    """
    def __init__(self,min=0.0,max=2.0*pi_constant):
      hour=pi_constant/12.0
      if(max-min)>12.0*hour:
        base=3.0*hour
      elif(max-min)>6.0*hour:
        base=2.0*hour
      # Put hour ticks if there are more than 3 hours displayed
      elif (max-min)>3.0*pi_constant/12.0:
        base=hour
      elif (max-min)>hour:
        base=hour/2.0
      else:
        base=hour/4.0
         
      matplotlib.ticker.MultipleLocator.__init__(self,base=base)

class DecLocator(matplotlib.ticker.MultipleLocator):
    """
    Dec tick locations with some intelligence
    """
    def __init__(self, min=-pi_constant/2.0,max=pi_constant/2.0):
      deg=pi_constant/180.0
      if (max-min)>60*deg:
        base=30.0*deg
      elif (max-min)>20*deg:
        base=10*deg
      elif (max-min)>10*deg:
        base=5*deg
      else:
        base=deg
      matplotlib.ticker.MultipleLocator.__init__(self,base=base)

class RAFormatter(matplotlib.ticker.FuncFormatter):
    def __init__(self,accuracy='auto'):
      matplotlib.ticker.FuncFormatter.__init__(self,getRAString)

class DecFormatter(matplotlib.ticker.FuncFormatter):
    def __init__(self,accuracy='auto'):
      matplotlib.ticker.FuncFormatter.__init__(self,getDecString)

def formatRATicks(locs, accuracy='auto'):
    """
    Format locs, ticks to RA angle with given accuracy
    accuracy can be 'hour', 'min', 'sec', 'all'
    returns (locs, ticks)
    'all' does no rounding, just formats the tick strings
    'auto' will use smallest appropriate units
    """
    newmax=max(locs)
    newmin=min(locs)
    if(accuracy=='auto'):
        acc='hour'
        if abs(newmax-newmin)<pi_constant/12.:
            acc='min'
        if abs(newmax-newmin)<pi_constant/(12.*60.):
            acc='sec'
    else:
        acc=accuracy
    
    if max(locs)>2*pi_constant: newmax=2.0*pi_constant
    if min(locs)<0.0: newmin=0.0
    locs=linspace(newmin,newmax,len(locs))
    
    roundlocs=map(lambda a: roundRadAngle(a, accuracy=acc), locs)
    
    newlocs=filter(lambda a:a>=0 and a<=2.0*pi_constant, roundlocs)
    return (list(newlocs), map(getRAString, list(newlocs) ) )

def formatDecTicks(locs, accuracy='auto'):
    """
    Format locs to Dec angle with given accuracy
    accuracy can be 'deg', 'arcmin', 'arcsec', 'all'
    'all' does no rounding, just formats the tick strings
    """
    newmax=max(locs)
    newmin=min(locs)
    if(accuracy=='auto'):
        acc='deg'
        if abs(newmax-newmin)<pi_constant/360.:
            acc='arcmin'
        if abs(newmax-newmin)<pi_constant/(360.*60.):
            acc='arcsec'
    else:
        acc=accuracy
    if newmax>0.5*pi_constant: newmax=0.5*pi_constant
    if newmin<-0.5*pi_constant: newmin=-0.5*pi_constant
    locs=linspace(newmin,newmax,len(locs))
    
    roundlocs=map(lambda a: roundRadAngle(a, accuracy=acc), locs)
    newlocs=filter(lambda a:a>=-pi_constant/2.0 and a<=pi_constant/2.0, roundlocs)
    return (list(newlocs), map(getDecString, list(newlocs) ) )

def roundRadAngle(rads,accuracy='all'):
    """
    round given angle in radians to integer hours, degrees, mins or secs
    accuracy can be 'hour'. 'deg', 'min', 'sec', 'all', all does nothing
    'arcmin', 'arcsec'
    """
    if accuracy=='all': return locs
    if accuracy=='hour': mult=24
    if accuracy=='deg': mult=360
    if accuracy=='min': mult=24*60
    if accuracy=='sec': mult=24*60*60
    if accuracy=='arcmin': mult=360*60
    if accuracy=='arcsec': mult=360*60*60
    mult=mult/(2.0*pi_constant)
    return round(rads*mult)/mult

def getRAString(radians,accuracy='auto'):
    secs=radians*12.0*3600/pi_constant
    hours, rem = divmod(secs, 3600 )
    mins,rem = divmod(rem, 60 )
    secs = rem
    if secs>=59.5:
        secs=secs-60
        mins=mins+1
    if mins>=59.5:
        mins=mins-60
        hours=hours+1
    if accuracy=='hour': return ur'%ih'%(hours)
    if accuracy=='min': return ur'%ih%im'%(hours,mins)
    if accuracy=='sec': return ur'%ih%im%2.0fs'%(hours,mins,secs)
    else:
        if abs(fmod(secs,60.0))>=0.5: return(getRAString(radians,accuracy='sec'))
        if abs(fmod(mins,60.0))>=0.5: return(getRAString(radians,accuracy='min'))
        else: return(getRAString(radians,accuracy='hour'))
        
def getDecString(radians,accuracy='auto'):
    # LaTeX doesn't like unicode degree symbols etc
    if matplotlib.rcParams['text.usetex']:
        degsymb='^\circ'
        minsymb="'"
        secsymb="''"
    else:
        degsymb=u'\u00B0'
        minsymb=u'\u0027'
        secsymb=u'\u2033'
    if(radians<0):
        radians=-radians
        sign=-1
    else: sign=+1
    deg,rem=divmod(radians,pi_constant/180.0)
    mins, rem = divmod(rem, pi_constant/(180.0*60.0))
    secs = rem * (180.0*3600.0)/pi_constant
    #if secs>=59.5:
    #    secs=secs-60.0
    #    mins=mins+1
    #if mins>=59.5:
    #    mins=mins-60.0
    #    deg=deg+1
    if (accuracy=='arcmin' or accuracy=='deg') and secs>30: mins=mins+1
    if accuracy=='deg' and mins>30: deg=deg+1
    if accuracy=='deg': return ur'%i'%(sign*deg)+degsymb
    if accuracy=='arcmin': return ur'%i%s%i%s'%(sign*deg,degsymb,mins,minsymb)
    if accuracy=='arcsec': return ur'%i%s%i%s%2.0f%s'%(sign*deg,degsymb,mins,minsymb,secs,secsymb)
    else:
    #    if abs(fmod(secs,60.0))>=0.5 and abs(fmod(secs,60)-60)>=0.5 : return(getDecString(sign*radians,accuracy='arcsec'))
    #    if abs(fmod(mins,60.0))>=0.5 and abs(fmod(mins,60)-60)>=0.5: return(getDecString(sign*radians,accuracy='arcmin'))
    #    else: return(getDecString(sign*radians,accuracy='deg'))
      return(getDecString(sign*radians,accuracy='deg'))

def plot_two_param_kde(posterior,plot2DkdeParams):
    """
    Plots a 2D kernel density estimate of the 2-parameter marginal posterior.

    @param posterior: an instance of the Posterior class.

    @param plot2DkdeParams: a dict {param1Name:Nparam1Bins,param2Name:Nparam2Bins}
    """

    from scipy import seterr as sp_seterr

    par1_name,par2_name=plot2DkdeParams.keys()
    Nx=plot2DkdeParams[par1_name]
    Ny=plot2DkdeParams[par2_name]

    xdat=posterior[par1_name].samples
    ydat=posterior[par2_name].samples

    par_injvalue1=posterior[par1_name].injval
    par_injvalue2=posterior[par2_name].injval

    par_trigvalues1=posterior[par1_name].trigvals
    par_trigvalues2=posterior[par2_name].trigvals

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    myfig=plt.figure(1,figsize=(6,4),dpi=200)
    myfig.add_axes(plt.Axes(myfig,[0.20,0.20,0.75,0.7]))
    plt.clf()

    xax=np.linspace(min(xdat),max(xdat),Nx)
    yax=np.linspace(min(ydat),max(ydat),Ny)
    x,y=np.meshgrid(xax,yax)

    samp=np.transpose(np.column_stack((xdat,ydat)))

    kde=stats.kde.gaussian_kde(samp)
    grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)

    z = kde(grid_coords.T)
    z = z.reshape(Nx,Ny)
    asp=xax.ptp()/yax.ptp()
#    if(asp<0.8 or asp > 1.6): asp=1.4
    plt.imshow(z,extent=(xax[0],xax[-1],yax[0],yax[-1]),aspect=asp,origin='lower')
    plt.colorbar()

    if par_injvalue1 is not None and par_injvalue2 is not None:
        plt.plot([par_injvalue1],[par_injvalue2],'bo',scalex=False,scaley=False)

    if par_trigvalues1 is not None and par_trigvalues2 is not None:
        par1IFOs = set([IFO for IFO in par_trigvalues1.keys()])
        par2IFOs = set([IFO for IFO in par_trigvalues2.keys()])
        IFOs = par1IFOs.intersection(par2IFOs)
        for IFO in IFOs:
            if IFO=='H1': color = 'r'
            elif IFO=='L1': color = 'g'
            elif IFO=='V1': color = 'm'
            else: color = 'c'
            plt.plot([par_trigvalues1[IFO]],[par_trigvalues2[IFO]],color=color,marker='o',scalex=False,scaley=False)

    plt.xlabel(par1_name)
    plt.ylabel(par2_name)
    plt.grid()

    # For RA and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            xmin,xmax=plt.xlim()
            plt.xlim(xmax,xmin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.xticks()
    #        (newlocs, newticks)=formatRATicks(locs, ticks)
    #        plt.xticks(newlocs,newticks,rotation=45)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.xticks()
    #        newlocs, newticks = formatDecTicks(locs)
    #        plt.xticks(newlocs,newticks,rotation=45)

    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    ymin,ymax=plt.ylim()
    #    plt.ylim(ymax,ymin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.yticks()
    #    newlocs,newticks=formatRATicks(locs)
    #    plt.yticks(newlocs,newticks)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.yticks()
    #    newlocs,newticks=formatDecTicks(locs)
    #    plt.yticks(newlocs,newticks)

    return myfig
#

def get_inj_by_time(injections,time):
    """
    Filter injections to find the injection with end time given by time +/- 0.1s
    """
    import itertools
    injection = itertools.ifilter(lambda a: abs(float(a.get_end()) - time) < 0.1, injections).next()
    return injection

def histogram2D(posterior,greedy2Params,confidence_levels):
    """
    Returns a 2D histogram and edges for the two parameters passed in greedy2Params, plus the actual discrete confidence levels
    imposed by the finite number of samples.
       H,xedges,yedges,Hlasts = histogram2D(posterior,greedy2Params,confidence_levels)
    @param posterior: Posterior instance
    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}
    @param confidence_levels: a list of the required confidence levels to plot on the contour map.
    """

    par1_name,par2_name=greedy2Params.keys()
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval
    a=np.squeeze(posterior[par1_name].samples)
    b=np.squeeze(posterior[par2_name].samples)
    par1pos_min=a.min()
    par2pos_min=b.min()
    par1pos_max=a.max()
    par2pos_max=b.max()
    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1
    H, xedges, yedges = np.histogram2d(a,b, bins=(par1pos_Nbins, par2pos_Nbins),normed=True)
    temp=np.copy(H)
    temp=temp.ravel()
    confidence_levels.sort()
    Hsum=0
    Hlasts=[]
    idxes=np.argsort(temp)
    j=len(idxes)-1
    for cl in confidence_levels:
        while float(Hsum/np.sum(H))<cl:
            #ind = np.argsort(temp)
            max_i=idxes[j]
            j-=1
            val = temp[max_i]
            Hlast=val
            Hsum+=val
            temp[max_i]=0
        Hlasts.append(Hlast)
    return (H,xedges,yedges,Hlasts)

def plot_two_param_greedy_bins_contourf(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,figsize=(7,6),dpi=120,figposition=[0.3,0.3,0.5,0.5],legend='right',hatches_by_name=None):
    """
    @param posteriors_by_name A dictionary of posterior objects indexed by name
    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}
    @param confidence_levels: a list of the required confidence levels to plot on the contour map. 
    
    """
    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()
    fig.add_axes(figposition)
    CSlst=[]
    name_list=[]
    par1_name,par2_name=greedy2Params.keys()
    for name,posterior in posteriors_by_name.items():
        name_list.append(name)
        H,xedges,yedges,Hlasts=histogram2D(posterior,greedy2Params,confidence_levels+[0.99999999999999])
        extent= [xedges[0], yedges[-1], xedges[-1], xedges[0]]
        CS2=plt.contourf(yedges[:-1],xedges[:-1],H,Hlasts,extend='max',colors=[colors_by_name[name]] ,alpha=0.3 )
        CS=plt.contour(yedges[:-1],xedges[:-1],H,Hlasts,extend='max',colors=[colors_by_name[name]] )
        CSlst.append(CS)
    
    plt.title("%s-%s confidence contours (greedy binning)"%(par1_name,par2_name)) # add a title
    plt.xlabel(par2_name)
    plt.ylabel(par1_name)
    if len(name_list)!=len(CSlst):
        raise RuntimeError("Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name.")
    full_name_list=[]
    dummy_lines=[]
    for plot_name in name_list:
        full_name_list.append(plot_name)
        for cl in confidence_levels+[1]:
            dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),color='k'))
            full_name_list.append('%s%%'%str(int(cl*100)))
        fig_actor_lst = [cs.collections[0] for cs in CSlst]
        fig_actor_lst.extend(dummy_lines)
    if legend is not None: twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')
    for text in twodcontour_legend.get_texts():
        text.set_fontsize('small')
    fix_axis_names(plt,par1_name,par2_name)
    return fig

def fix_axis_names(plt,par1_name,par2_name):
    """
    Fixes names of axes
    """
    return

    # For ra and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            ymin,ymax=plt.ylim()
            if(ymin<0.0): ylim=0.0
            if(ymax>2.0*pi_constant): ymax=2.0*pi_constant
            plt.ylim(ymax,ymin)
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            locs, ticks = plt.yticks()
            newlocs, newticks = formatRATicks(locs)
            plt.yticks(newlocs,newticks)
    if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
            locs, ticks = plt.yticks()
            newlocs,newticks=formatDecTicks(locs)
            plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        locs, ticks = plt.xticks()
        newlocs, newticks = formatRATicks(locs)
        plt.xticks(newlocs,newticks,rotation=45)
    if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
        locs, ticks = plt.xticks()
        newlocs, newticks = formatDecTicks(locs)
        plt.xticks(newlocs,newticks,rotation=45)
    return plt

def plot_two_param_greedy_bins_contour(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,line_styles=__default_line_styles,figsize=(7,6),dpi=250,figposition=[0.2,0.2,0.48,0.75],legend='right'):
    """
    Plots the confidence level contours as determined by the 2-parameter
    greedy binning algorithm.

    @param posteriors_by_name: A dict containing Posterior instances referenced by some id.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}

    @param confidence_levels: a list of the required confidence levels to plot on the contour map.

    @param colors_by_name: A dict of colors cross-referenced to the above Posterior ids.

    @param legend: Argument for legend placement or None for no legend ('right', 'upper left', 'center' etc)

    """

    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()

    axes=fig.add_axes(figposition)

    #This fixes the precedence of line styles in the plot
    if len(line_styles)<len(confidence_levels):
        raise RuntimeError("Error: Need as many or more line styles to choose from as confidence levels to plot!")

    CSlst=[]
    name_list=[]
    for name,posterior in posteriors_by_name.items():

        name_list.append(name)
        #Extract parameter names
        par1_name,par2_name=greedy2Params.keys()
        #Extract bin sizes
        par1_bin=greedy2Params[par1_name]
        par2_bin=greedy2Params[par2_name]

        #Extract injection information
        par1_injvalue=posterior[par1_name.lower()].injval
        par2_injvalue=posterior[par2_name.lower()].injval

        #Extract trigger information
        par1_trigvalues=posterior[par1_name.lower()].trigvals
        par2_trigvalues=posterior[par2_name.lower()].trigvals

        a=np.squeeze(posterior[par1_name].samples)
        b=np.squeeze(posterior[par2_name].samples)

        #Create 2D bin array
        par1pos_min=a.min()
        par2pos_min=b.min()

        par1pos_max=a.max()
        par2pos_max=b.max()

        par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
        par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

        if par1_name.find('time')!=-1:
          offset=floor(min(a))
          a=a-offset
          ax1_name=par1_name+' + %i'%(int(offset))
        else: ax1_name=par1_name

        if par2_name.find('time')!=-1:
          offset=floor(min(b))
          b=b-offset
          ax2_name=par2_name+' + %i'%(int(offset))
        else: ax2_name=par2_name


        majorFormatterX=ScalarFormatter(useMathText=True)
        majorFormatterX.format_data=lambda data:'%.4g'%(data)
        majorFormatterY=ScalarFormatter(useMathText=True)
        majorFormatterY.format_data=lambda data:'%.4g'%(data)
        majorFormatterX.set_scientific(True)
        majorFormatterY.set_scientific(True)
        axes.xaxis.set_major_formatter(majorFormatterX)
        axes.yaxis.set_major_formatter(majorFormatterY)
        
        H, xedges, yedges = np.histogram2d(a,b, bins=(par1pos_Nbins, par2pos_Nbins),normed=True)

        extent = [xedges[0], yedges[-1], xedges[-1], xedges[0]]

        temp=np.copy(H)
        temp=temp.ravel()
        confidence_levels.sort()
        Hsum=0
        Hlasts=[]
        idxes=np.argsort(temp)
        j=len(idxes)-1
        for cl in confidence_levels:
            while float(Hsum/np.sum(H))<cl:
                #ind = np.argsort(temp)
                max_i=idxes[j]
                j-=1
                val = temp[max_i]
                Hlast=val
                Hsum+=val
                temp[max_i]=0
            Hlasts.append(Hlast)

        CS=plt.contour(yedges[:-1],xedges[:-1],H,Hlasts,colors=[colors_by_name[name]],linestyles=line_styles)
        plt.grid()
        if(par1_injvalue is not None and par2_injvalue is not None):
            plt.plot([par2_injvalue],[par1_injvalue],'b*',scalex=False,scaley=False)
        if(par1_trigvalues is not None and par2_trigvalues is not None):
            par1IFOs = set([IFO for IFO in par1_trigvalues.keys()])
            par2IFOs = set([IFO for IFO in par2_trigvalues.keys()])
            IFOs = par1IFOs.intersection(par2IFOs)
            if IFO=='H1': color = 'r'
            elif IFO=='L1': color = 'g'
            elif IFO=='V1': color = 'm'
            else: color = 'c'
            plt.plot([par2_trigvalues[IFO]],[par1_trigvalues[IFO]],color=color,marker='*',scalex=False,scaley=False)
        CSlst.append(CS)

    	Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    	if Nchars>8:
      		Nticks=3
    	elif Nchars>5:
      		Nticks=4
    	elif Nchars>4:
      		Nticks=5
    	else:
      		Nticks=6
    	locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks-1)
        if par2_name=='rightascension' or par2_name=='ra':
            (ramin,ramax)=plt.xlim()
            locatorX=RALocator(min=ramin,max=ramax)
            majorFormatterX=RAFormatter()
        if par2_name=='declination' or par2_name=='dec':
            (decmin,decmax)=plt.xlim()
            locatorX=DecLocator(min=decmin,max=decmax)
            majorFormatterX=DecFormatter()
        axes.xaxis.set_major_formatter(majorFormatterX)
        if par1_name=='rightascension' or par1_name=='ra':
            (ramin,ramax)=plt.ylim()
            locatorY=RALocator(ramin,ramax)
            axes.yaxis.set_major_locator(locatorY)
            majorFormatterY=RAFormatter()
        if par1_name=='declination' or par1_name=='dec':
            (decmin,decmax)=plt.ylim()
            locatorY=DecLocator(min=decmin,max=decmax)
            majorFormatterY=DecFormatter()
            axes.yaxis.set_major_locator(locatorY)

        axes.yaxis.set_major_formatter(majorFormatterY)
    	#locatorX.view_limits(bins[0],bins[-1])
    	axes.xaxis.set_major_locator(locatorX)

    plt.title("%s-%s confidence contours (greedy binning)"%(par1_name,par2_name)) # add a title
    plt.xlabel(ax2_name)
    plt.ylabel(ax1_name)

    if len(name_list)!=len(CSlst):
        raise RuntimeError("Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name.")
    full_name_list=[]
    dummy_lines=[]

    for plot_name in name_list:
        full_name_list.append(plot_name)
    for ls_,cl in zip(line_styles[0:len(confidence_levels)],confidence_levels):
        dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),ls=ls_,color='k'))
        full_name_list.append('%s%%'%str(int(cl*100)))

    fig_actor_lst = [cs.collections[0] for cs in CSlst]

    fig_actor_lst.extend(dummy_lines)

    if legend is not None: twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')

    for text in twodcontour_legend.get_texts():
        text.set_fontsize('small')


    # For ra and dec set custom labels and for RA reverse
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        ymin,ymax=plt.ylim()
    #        if(ymin<0.0): ylim=0.0
    #        if(ymax>2.0*pi_constant): ymax=2.0*pi_constant
    #        plt.ylim(ymax,ymin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatRATicks(locs)
    #        plt.yticks(newlocs,newticks)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.yticks()
    #        newlocs,newticks=formatDecTicks(locs)
    #        plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return fig
#

def plot_two_param_greedy_bins_hist(posterior,greedy2Params,confidence_levels):
    """
    Histograms of the ranked pixels produced by the 2-parameter greedy
    binning algorithm colured by their confidence level.

    @param toppoints: Nx2 array of 2-parameter posterior samples.

    @param posterior: an instance of the Posterior class.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}
    """

    from scipy import seterr as sp_seterr

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    #Extract parameter names
    par1_name,par2_name=greedy2Params.keys()
    #Extract bin sizes
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]

    a=np.squeeze(posterior[par1_name].samples)
    b=np.squeeze(posterior[par2_name].samples)

    #Create 2D bin array
    par1pos_min=a.min()
    par2pos_min=b.min()

    par1pos_max=a.max()
    par2pos_max=b.max()

    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

    # Adjust for time parameter
    if par1_name.find('time')!=-1:
      offset=floor(min(a))
      a=a-offset
      ax1_name=par1_name+' + %i'%(int(offset))
    else: ax1_name=par1_name

    if par2_name.find('time')!=-1:
      offset=floor(min(b))
      b=b-offset
      ax2_name=par2_name+' + %i'%(int(offset))
    else: ax2_name=par2_name

    #Extract injection information
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval

    #Extract trigger information
    par1_trigvalues=posterior[par1_name.lower()].trigvals
    par2_trigvalues=posterior[par2_name.lower()].trigvals

    myfig=plt.figure()
    axes=plt.Axes(myfig,[0.3,0.3,0.95-0.3,0.90-0.3])
    myfig.add_axes(axes)
    
    #plt.clf()
    plt.xlabel(ax2_name)
    plt.ylabel(ax1_name)

    #bins=(par1pos_Nbins,par2pos_Nbins)
    bins=(50,50) # Matches plot_one_param_pdf
    
    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.4g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.4g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)
    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)
    H, xedges, yedges = np.histogram2d(a,b, bins,normed=False)

      
    #Replace H with greedy bin confidence levels at each pixel...
    temp=np.copy(H)
    temp=temp.flatten()

    Hsum=0
    Hsum_actual=np.sum(H)
    
    idxes=np.argsort(temp)
    j=len(idxes)-1
    while Hsum<Hsum_actual:
        #ind = np.argsort(temp)
        max_i=idxes[j]
        j-=1
        val = temp[max_i]
        Hsum+=int(val)
        temp[max_i]=0

        #print Hsum,Hsum_actual
        H.flat[max_i]=1-float(Hsum)/float(Hsum_actual)

    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    plt.imshow(np.flipud(H), axes=axes, aspect='auto', extent=extent, interpolation='nearest',cmap='gray_r')
    plt.gca().autoscale_view()
    plt.colorbar()
    
    #plt.hexbin(a,b,cmap='gray_r',axes=axes )
    
    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
      Nticks=3
    elif Nchars>5:
      Nticks=4
    elif Nchars>4:
      Nticks=5
    else:
      Nticks=6
    locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks-1)
    (xmin,xmax)=plt.xlim()
    (ymin,ymax)=plt.ylim()
    if par2_name=='rightascension' or par2_name=='ra':
        locatorX=RALocator(min=xmin,max=xmax)
        majorFormatterX=RAFormatter()
    if par2_name=='declination' or par2_name=='dec':
        locatorX=DecLocator(min=xmin,max=xmax)
        majorFormatterX=DecFormatter()
    if par1_name=='rightascension' or par1_name=='ra':
        locatorY=RALocator(min=ymin,max=ymax)
        axes.yaxis.set_major_locator(locatorY)
        majorFormatterY=RAFormatter()
    if par1_name=='declination' or par1_name=='dec':
        locatorY=DecLocator(min=ymin,max=ymax)
        axes.yaxis.set_major_locator(locatorY)
        majorFormatterY=DecFormatter()

    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)
    #locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)

    if par1_injvalue is not None and par2_injvalue is not None:
        plt.plot([par1_injvalue],[par2_injvalue],'bo',scalex=False,scaley=False)

    if par1_trigvalues is not None and par2_trigvalues is not None:
        par1IFOs = set([IFO for IFO in par1_trigvalues.keys()])
        par2IFOs = set([IFO for IFO in par2_trigvalues.keys()])
        IFOs = par1IFOs.intersection(par2IFOs)
        if IFO=='H1': color = 'r'
        elif IFO=='L1': color = 'g'
        elif IFO=='V1': color = 'm'
        else: color = 'c'
        plt.plot([par1_trigvalues[IFO]],[par2_trigvalues[IFO]],color=color,marker='o',scalex=False,scaley=False)

    # For RA and dec set custom labels and for RA reverse
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        ymin,ymax=plt.ylim()
    #        plt.ylim(ymax,ymin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatRATicks(locs)
    #        plt.yticks(newlocs,newticks)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatDecTicks(locs)
    #        plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin)<0.0: xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return myfig

def greedy_bin_one_param(posterior,greedy1Param,confidence_levels):
    """
    Determine the 1-parameter Bayesian Confidence Interval using a greedy
    binning algorithm.

    @param posterior: an instance of the posterior class.

    @param greedy1Param: a dict; {paramName:paramBinSize}.

    @param confidence_levels: A list of floats of the required confidence intervals [(0-1)].
    """

    paramName=greedy1Param.keys()[0]
    par_bin=greedy1Param.values()[0]
    par_samps=posterior[paramName.lower()].samples

    parpos_min=min(par_samps)[0]
    parpos_max=max(par_samps)[0]

    par_point=parpos_min

    parpos_Nbins= int(ceil((parpos_max - parpos_min)/par_bin))+1

    greedyPoints=np.zeros((parpos_Nbins,2))
    # ...NB 2D so it can be put through same confidence level function
    greedyHist=np.zeros(parpos_Nbins,dtype='i8')

    #Bin up
    for i in range(parpos_Nbins):
        greedyPoints[i,0]=par_point
        greedyPoints[i,1]=par_point
        par_point+=par_bin

    for par_samp in par_samps:
        par_samp=par_samp[0]
        par_binNumber=int(floor((par_samp-parpos_min)/par_bin))
        try:
            greedyHist[par_binNumber]+=1
        except IndexError:
            print "IndexError: bin number: %i total bins: %i parsamp: %f "\
                %(par_binNumber,parpos_Nbins,par_samp)

    #Find injection bin
    injbin=None
    par_injvalue=posterior[paramName].injval
    if par_injvalue:
        par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
        injbin=par_binNumber

    toppoints,injectionconfidence,reses,injection_area=_greedy_bin(greedyHist,greedyPoints,injbin,float(par_bin),int(len(par_samps)),confidence_levels)
    cl_intervals=[]
    confidence_levels.sort()
    for cl in confidence_levels:
        ind=np.nonzero(toppoints[:,-1]<cl)

        if len(ind[0]) > 1:
            cl_intervals.append((np.min(toppoints[ind,0]),np.max(toppoints[ind,0])))

        else:

            cl_intervals.append((toppoints[ind[0],0],toppoints[ind[0],0]))

    return toppoints,injectionconfidence,reses,injection_area,cl_intervals

#
def contigious_interval_one_param(posterior,contInt1Params,confidence_levels):
    """
    Calculates the smallest contigious 1-parameter confidence interval for a
    set of given confidence levels.

    @param posterior: an instance of the Posterior class.

    @param contInt1Params: a dict {paramName:paramBinSize}.

    @param confidence_levels: Required confidence intervals.

    """
    oneDContCL={}
    oneDContInj={}

    paramName=contInt1Params.keys()[0]
    par_bin=contInt1Params.values()[0]

    par_injvalue=posterior[paramName].injval

    par_samps=posterior[paramName].samples

    parpos_min=min(par_samps)
    parpos_max=max(par_samps)

    par_point=parpos_min
    parpos_Nbins= int(ceil((parpos_max - parpos_min)/par_bin))+1

    greedyHist=np.zeros(parpos_Nbins,dtype='i8')

    for par_samp in par_samps:
        par_binNumber=int(floor((par_samp-parpos_min)/par_bin))
        try:
            greedyHist[par_binNumber]+=1
        except IndexError:
            print "IndexError: bin number: %i total bins: %i parsamp: %f bin: %f - %f"\
                %(
                  par_binNumber,
                  parpos_Nbins,
                  par_samp,
                  greedyPoints[par_binNumber-1,0],
                  greedyPoints[par_binNumber-1,0]+par_bin
                  )

    injbin=None
    #Find injection bin
    if par_injvalue:
        par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
        injbin=par_binNumber

    j=0
    #print "Calculating contigious confidence intervals for %s..."%par_name
    len_par_samps=len(par_samps)

    injinterval=None

    #Determine smallest contigious interval for given confidence levels (brute force)
    while j < len(confidence_levels):
        confidence_level=confidence_levels[j]
        #Loop over size of interval
        max_left=0
        max_right=0

        for i in range(len(greedyHist)):

            max_frac=None
            left=0
            right=i

            #Slide interval
            while right<len(greedyHist):
                Npoints=sum(greedyHist[left:right])
                frac=float(Npoints)/float(len_par_samps)
                #print "left %i , right %i , frac %f"%(left,right,frac)

                if max_frac is None:
                    max_frac=frac
                    max_left=left
                    max_right=right
                else:
                    if frac>max_frac:
                        max_frac=frac
                        max_left=left
                        max_right=right

                left+=1
                right+=1

            if injbin is not None and injinterval is None:
                if injbin in range(max_left,max_right):
                    injinterval=(max_right-max_left)*par_bin
                    oneDContInj['interval']=injinterval
                    oneDContInj['confidence']=1-frac
            if max_frac > confidence_level:
                break

            max_frac=None

        if max_frac is None:
            print "Cant determine intervals at %f confidence!"%confidence_level
        else:

            oneDContCL['left']=max_left*par_bin
            oneDContCL['right']=max_right*par_bin
            oneDContCL['width']=(max_right-max_left)*par_bin
            k=j
            while k+1<len(confidence_levels) :
                if confidence_levels[k+1]<max_frac:
                    j+=1
                k+=1
        j+=1

    return oneDContCL,oneDContInj
#
def burnin(data,spin_flag,deltaLogL,outputfile):

    pos,bayesfactor=_burnin(data,spin_flag,deltaLogL,outputfile)

    return pos,bayesfactor


class ACLError(StandardError):
    def __init__(self, *args):
        super(ACLError, self).__init__(*args)


def autocorrelation(series):
    """Returns an estimate of the autocorrelation function of a given
    series.  Returns only the positive-lag portion of the ACF,
    normalized so that the zero-th element is 1."""
    x=series-np.mean(series) 
    y=np.conj(x[::-1])

    acf=np.fft.ifftshift(signal.fftconvolve(y,x,mode='full'))

    N=series.shape[0]

    acf = acf[0:N]

    return acf/acf[0]


def autocorrelation_length_estimate(series, acf=None, M=5, K=2):
    """Attempts to find a self-consistent estimate of the
    autocorrelation length of a given series.  

    If C(tau) is the autocorrelation function (normalized so C(0) = 1,
    for example from the autocorrelation procedure in this module),
    then the autocorrelation length is the smallest s such that

    1 + 2*C(1) + 2*C(2) + ... + 2*C(M*s) < s

    In words: the autocorrelation length is the shortest length so
    that the sum of the autocorrelation function is smaller than that
    length over a window of M times that length.

    The maximum window length is restricted to be len(series)/K as a
    safety precaution against relying on data near the extreme of the
    lags in the ACF, where there is a lot of noise.  Note that this
    implies that the series must be at least M*K*s samples long in
    order to get a reliable estimate of the ACL.

    If no such s can be found, raises ACLError; in this case it is
    likely that the series is too short relative to its true
    autocorrelation length to obtain a consistent ACL estimate."""

    if acf is None:
      acf=autocorrelation(series)
    acf[1:] *= 2.0

    imax=int(acf.shape[0]/K)
    
    # Cumulative sum and ACL length associated with each window
    cacf=np.cumsum(acf)
    s=np.arange(1, cacf.shape[0]+1)/float(M)

    # Find all places where cumulative sum over window is smaller than
    # associated ACL.
    estimates=np.flatnonzero(cacf[:imax] < s[:imax])

    if estimates.shape[0] > 0:
        # Return the first index where cumulative sum is smaller than
        # ACL associated with that index's window
        return s[estimates[0]]
    else:
        # Cannot find self-consistent ACL estimate.
        raise ACLError('autocorrelation length too short for consistent estimate')


def effectiveSampleSize(samples, Nskip=1):
    """
    Compute the effective sample size, calculating the ACL using only
    the second half of the samples to avoid ACL overestimation due to
    chains equilibrating after adaptation.
    """
    N = len(samples)
    acf = autocorrelation(samples[N/2:])
    try:
      acl = autocorrelation_length_estimate(samples[N/2:], acf=acf)
    except ACLError:
      acl = N
    Neffective = floor(N/acl)
    acl *= Nskip
    return (Neffective, acl, acf)


def readCoincXML(xml_file, trignum):
    triggers=None

    coincXML = utils.load_filename(xml_file)
    coinc = lsctables.getTablesByType(coincXML, lsctables.CoincTable)[0]
    coincMap = lsctables.getTablesByType(coincXML, lsctables.CoincMapTable)[0]
    snglInsps = lsctables.getTablesByType(coincXML, lsctables.SnglInspiralTable)[0]

    if (trignum>len(coinc)):
        raise RuntimeError("Error: You asked for trigger %d, but %s contains only %d triggers" %(trignum,coincfile,len(tiggers)))
    else:
        coincEventID = coinc.getColumnByName('coinc_event_id')[trignum]
        eventIDs = [row.event_id for row in coincMap if row.coinc_event_id == coincEventID]
        triggers = [row for row in snglInsps if row.event_id in eventIDs]
    return triggers

#===============================================================================
# Parameter estimation codes results parser
#===============================================================================

class PEOutputParser(object):
    """
    A parser for the output of Bayesian parameter estimation codes.

    TODO: Will be abstract class when LDG moves over to Python >2.6,
    inherited by each method .
    """
    def __init__(self,inputtype):
        if inputtype is 'mcmc_burnin':
            self._parser=self._mcmc_burnin_to_pos
        elif inputtype is 'ns':
            self._parser=self._ns_to_pos
        elif inputtype is 'common':
            self._parser=self._common_to_pos
        elif inputtype is 'fm':
            self._parser=self._followupmcmc_to_pos
        elif inputtype is "inf_mcmc":
            self._parser=self._infmcmc_to_pos
        elif inputtype is "xml":
            self._parser=self._xml_to_pos

    def parse(self,files,**kwargs):
        """
        Parse files.
        """
        return self._parser(files,**kwargs)

    def _infmcmc_to_pos(self,files,outdir=None,deltaLogL=None,fixedBurnins=None,nDownsample=None,oldMassConvention=False,**kwargs):
        """
        Parser for lalinference_mcmcmpi output.
        """
        if not (fixedBurnins is None):
            if not (deltaLogL is None):
                print "Warning: using deltaLogL criteria in addition to fixed burnin"
            if len(fixedBurnins) == 1 and len(files) > 1:
                print "Only one fixedBurnin criteria given for more than one output.  Applying this to all outputs."
                fixedBurnins = np.ones(len(files),'int')*fixedBurnins[0]
            elif len(fixedBurnins) != len(files):
                raise RuntimeError("Inconsistent number of fixed burnin criteria and output files specified.")
            print "Fixed burning criteria: ",fixedBurnins
        else:
            fixedBurnins = np.zeros(len(files))
        logLThreshold=-1e200 # Really small?
        if not (deltaLogL is None):
            logLThreshold=self._find_max_logL(files) - deltaLogL
            print "Eliminating any samples before log(Post) = ", logLThreshold
        nskips=self._find_ndownsample(files, logLThreshold, fixedBurnins, nDownsample)
        if nDownsample is None:
            print "Downsampling to take only uncorrelated posterior samples from each file."
            if len(nskips) == 1 and np.isnan(nskips[0]):
                print "WARNING: All samples in chain are correlated.  Downsampling to 10000 samples for inspection!!!"
                nskips=self._find_ndownsample(files, logLThreshold, fixedBurnins, 10000)
            else:
                for i in range(len(nskips)):
                    if np.isnan(nskips[i]):
                        print "%s eliminated since all samples are correlated."
                    else:
                        print "Downsampling by a factor of ", nskips[0], " to achieve approximately ", nDownsample, " posterior samples"
        if outdir is None:
            outdir=''
        runfileName=os.path.join(outdir,"lalinfmcmc_headers.dat")
        postName="posterior_samples.dat"
        runfile=open(runfileName, 'w')
        outfile=open(postName, 'w')
        try:
            self._infmcmc_output_posterior_samples(files, runfile, outfile, logLThreshold, fixedBurnins, nskips, oldMassConvention)
        finally:
            runfile.close()
            outfile.close()
        return self._common_to_pos(open(postName,'r'))


    def _infmcmc_output_posterior_samples(self, files, runfile, outfile, logLThreshold, fixedBurnins, nskips=None, oldMassConvention=False):
        """
        Concatenate all the samples from the given files into outfile.
        For each file, only those samples past the point where the
        log(L) > logLThreshold are concatenated after eliminating
        fixedBurnin.
        """
        nRead=0
        outputHeader=False
        acceptedChains=0
        if nskips is None:
            nskips = np.ones(len(files),'int')
        for infilename,i,nskip,fixedBurnin in zip(files,range(1,len(files)+1),nskips,fixedBurnins):
            infile=open(infilename,'r')
            try:
                print "Writing header of %s to %s"%(infilename,runfile.name)
                runInfo,header=self._clear_infmcmc_header(infile)
                runfile.write('Chain '+str(i)+':\n')
                runfile.writelines(runInfo)
                print "Processing file %s to %s"%(infilename,outfile.name)
                f_lower=self._find_infmcmc_f_lower(runInfo)
                if oldMassConvention:
                    # Swap #1 for #2 because our old mass convention
                    # has m2 > m1, while the common convention has m1
                    # > m2
                    header=[self._swaplabel12(label) for label in header]
                if not outputHeader:
                    for label in header:
                        outfile.write(label)
                        outfile.write(" ")
                    outfile.write("f_lower")
                    outfile.write(" ")
                    outfile.write("chain")
                    outfile.write("\n")
                    outputHeader=header
                iterindex=header.index("cycle")
                loglindex=header.index("logpost")
                output=False
                for line in infile:
                    line=line.lstrip()
                    lineParams=line.split()
                    iter=int(lineParams[iterindex])
                    logL=float(lineParams[loglindex])
                    if (iter > fixedBurnin) and (logL >= logLThreshold):
                        output=True
                    if output:
                        if nRead % nskip == 0:
                            for label in outputHeader:
                                # Note that the element "a1" in the
                                # *header* actually already
                                # corresponds to the "a2" *column* of
                                # the input because we switched the
                                # names above
                                outfile.write(lineParams[header.index(label)])
                                outfile.write(" ")
                            outfile.write(f_lower)
                            outfile.write(" ")
                            outfile.write(str(i))
                            outfile.write("\n")
                        nRead=nRead+1
                if output: acceptedChains += 1
            finally:
                infile.close()
        print "%i of %i chains accepted."%(acceptedChains,len(files))

    def _swaplabel12(self, label):
        if label[-1] == '1':
            return label[0:-1] + '2'
        elif label[-1] == '2':
            return label[0:-1] + '1'
        else:
            return label[:]

    def _find_max_logL(self, files):
        """
        Given a list of files, reads them, finding the maximum log(L)
        """
        maxLogL = -1e200  # Really small, I hope!
        for inpname in files:
            infile=open(inpname, 'r')
            try:
                runInfo,header=self._clear_infmcmc_header(infile)
                loglindex=header.index("logpost")
                for line in infile:
                    line=line.lstrip().split()
                    logL=float(line[loglindex])
                    if logL > maxLogL:
                        maxLogL=logL
            finally:
                infile.close()
        print "Found max log(Post) = ", maxLogL
        return maxLogL

    def _find_ndownsample(self, files, logLthreshold, fixedBurnins, nDownsample):
        """
        Given a list of files, threshold value, and a desired
        number of outputs posterior samples, return the skip number to
        achieve the desired number of posterior samples.
        """
        nfiles = len(files)
        ntots=[]
        nEffectives = []
        for inpname,fixedBurnin in zip(files,fixedBurnins):
            infile = open(inpname, 'r')
            try:
                runInfo,header = self._clear_infmcmc_header(infile)
                header = [name.lower() for name in header]
                loglindex = header.index("logpost")
                iterindex = header.index("cycle")
                deltaLburnedIn = False
                fixedBurnedIn  = False
                lines=[]
                ntot=0
                for line in infile:
                    line = line.lstrip().split()
                    iter = int(line[iterindex])
                    logL = float(line[loglindex])
                    if iter > fixedBurnin:
                        fixedBurnedIn = True
                    if logL > logLthreshold:
                        deltaLburnedIn = True
                    if fixedBurnedIn and deltaLburnedIn:
                        ntot += 1
                        lines.append(line)
                ntots.append(ntot)
                if nDownsample is None:
                    try:
                        nonParams = ["logpost", "cycle", "logprior", "logl", "loglh1", "logll1", "loglv1"]
                        nonParamsIdxs = [header.index(name) for name in nonParams if name in header]
                        paramIdxs = [i for i in range(len(header)) if i not in nonParamsIdxs]
                        samps = np.array(lines).astype(float)
                        nEffectives.append(min([effectiveSampleSize(samps[:,i])[0] for i in paramIdxs]))
                    except:
                        nEffectives.append(None)
                        print "Error computing effective sample size of %s!"%inpname

            finally:
                infile.close()
        nskips = np.ones(nfiles)
        ntot = sum(ntots)
        if nDownsample is not None:
            if ntot > nDownsample:
                nskips *= floor(ntot/nDownsample)

        else:
            for i in range(nfiles):
                nEff = nEffectives[i]
                ntot = ntots[i]
                if nEff > 1:
                    if ntot > nEff:
                        nskips[i] = ceil(ntot/nEff)
                else:
                    nskips[i] = None
        return nskips

    def _find_infmcmc_f_lower(self, runInfo):
        """
        Searches through header to determine starting frequency of waveforms.
        Assumes same for all IFOs.
        """
        runInfo = iter(runInfo)
        for line in runInfo:
            headers=line.lstrip().lower().split()
            try:
                flowColNum = headers.index('f_low')
                IFOinfo = runInfo.next().lstrip().lower().split()
                f_lower = IFOinfo[flowColNum]
                break
            except ValueError:
                continue
        return f_lower

    def _clear_infmcmc_header(self, infile):
        """
        Reads lalinference_mcmcmpi file given, returning the run info and 
        common output header information.
        """
        runInfo = []
        for line in infile:
            runInfo.append(line)
            headers=line.lstrip().lower().split()
            try:
                headers.index('cycle')
                break
            except ValueError:
                continue
        else:
            raise RuntimeError("couldn't find line with 'cycle' in LALInferenceMCMC input")
        return runInfo[:-1],headers


    def _mcmc_burnin_to_pos(self,files,spin=False,deltaLogL=None):
        """
        Parser for SPINspiral output .
        """
        raise NotImplementedError
        if deltaLogL is not None:
            pos,bayesfactor=burnin(data,spin,deltaLogL,"posterior_samples.dat")
            return self._common_to_pos(open("posterior_samples.dat",'r'))

    def _ns_to_pos(self,files,Nlive=None,Npost=10000):
        """
        Parser for nested sampling output.
        files : list of input NS files
        Nlive : Number of live points
        Npost : Desired number of posterior samples
        """
        try:
            from lalapps.nest2pos import draw_N_posterior_many,draw_posterior_many
        except ImportError:
            print "Need lalapps.nest2pos to convert nested sampling output!"
            raise

        if Nlive is None:
            raise RuntimeError("Need to specify number of live points in positional arguments of parse!")
                       
        posfilename='posterior_samples.dat'
       
        #posfile.write('mchirp \t eta \t time \t phi0 \t dist \t RA \t dec \t
        #psi \t iota \t likelihood \n')
        # get parameter list
        it = iter(files)
        
        # check if there's a file containing the parameter names
        parsfilename = it.next()+'_params.txt'
        
        if os.path.isfile(parsfilename):
            print 'Looking for '+parsfilename

            if os.access(parsfilename,os.R_OK):

                with open(parsfilename,'r') as parsfile: 
                    outpars=parsfile.readline()+'\n'
            else:
                raise RuntimeError('Cannot open parameters file %s!'%(parsfilename))

        else: # Use hardcoded CBC parameter names
            outpars='mchirp \t eta \t time \t phi0 \t dist \t RA \t \
            dec \t psi \t iota \t logl \n'

        # Find the logL column
        parsvec=outpars.split()
        logLcol=-1
        for i in range(len(parsvec)):
            if parsvec[i].lower()=='logl':
                logLcol=i
        if logLcol==-1:
            print 'Error! Could not find logL column in parameter list: %s'%(outpars)
            raise RuntimeError

        inarrays=map(np.loadtxt,files)
        if Npost is None:
            pos=draw_posterior_many(inarrays,[Nlive for f in files],logLcol=logLcol)
        else:
            pos=draw_N_posterior_many(inarrays,[Nlive for f in files],Npost,logLcol=logLcol)

        with open(posfilename,'w') as posfile:
            
            posfile.write(outpars)
        
            for row in pos:
                for i in row:
                    posfile.write('%.12g\t' %(i))
                posfile.write('\n')
        
        with open(posfilename,'r') as posfile:
            return_val=self._common_to_pos(posfile)
        
        return return_val

    def _followupmcmc_to_pos(self,files):
        """
        Parser for followupMCMC output.
        """
        return self._common_to_pos(open(files[0],'r'),delimiter=',')


    def _multinest_to_pos(self,files):
        """
        Parser for MultiNest output.
        """
        return self._common_to_pos(open(files[0],'r'))

    def _xml_to_pos(self,infile):
        """
        Parser for VOTable XML Using
        """
        from xml.etree import ElementTree as ET
        xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
        try:
                register_namespace=ET.register_namespace
        except AttributeError:
                def register_namespace(prefix,uri):
                    ET._namespace_map[uri]=prefix
        register_namespace('vot',xmlns)
        tree = ET.ElementTree()
        
        tree.parse(infile)
        # Find the posterior table
        tables = tree.findall('.//{%s}TABLE'%(xmlns))
        for table in tables:
            if table.get('utype')=='lalinference:results:posteriorsamples':
                return(self._VOTTABLE2pos(table))
        for table in tables:
          if table.get('utype')=='lalinference:results:nestedsamples':
            nsresource=[node for node in tree.findall('{%s}RESOURCE'%(xmlns)) if node.get('utype')=='lalinference:results'][0]
            return(self._VOTTABLE2pos(vo_nest2pos(nsresource)))
        raise RuntimeError('Cannot find "Posterior Samples" TABLE element in XML input file %s'%(infile))
        
    def _VOTTABLE2pos(self,table):
        """
        Parser for a VOT TABLE element with FIELDs and TABLEDATA elements
        """
        from xml.etree import ElementTree as ET
        xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
        try:
            register_namespace=ET.register_namespace
        except AttributeError:
            def register_namespace(prefix,uri):
                ET._namespace_map[uri]=prefix
                register_namespace('vot',xmlns)

        header=[]
        for field in table.findall('./{%s}FIELD'%(xmlns)):
            header.append(field.attrib['name'])
        if(len(header)==0):
            raise RuntimeError('Unable to find FIELD nodes for table headers in XML table')
        data=table.findall('./{%s}DATA'%(xmlns))
        tabledata=data[0].find('./{%s}TABLEDATA'%(xmlns))
        llines=[]
        for row in tabledata:
            llines.append(np.array(map(lambda a:float(a.text),row)))
        flines=np.array(llines)
        for i in range(0,len(header)):
            if header[i].lower().find('log')!=-1 and header[i].lower() not in logParams:
                print 'exponentiating %s'%(header[i])

                flines[:,i]=np.exp(flines[:,i])

                header[i]=header[i].replace('log','')
            if header[i].lower().find('sin')!=-1:
                print 'asining %s'%(header[i])
                flines[:,i]=np.arcsin(flines[:,i])
                header[i]=header[i].replace('sin','')
            if header[i].lower().find('cos')!=-1:
                print 'acosing %s'%(header[i])
                flines[:,i]=np.arccos(flines[:,i])
                header[i]=header[i].replace('cos','')
            header[i]=header[i].replace('(','')
            header[i]=header[i].replace(')','')
        print 'Read columns %s'%(str(header))
        return header,flines


    def _common_to_pos(self,infile,info=[None,None]):
        """
        Parse a file in the 'common format' and return an array of posterior
        samples and list of parameter names. Will apply inverse functions to
        columns with names containing sin,cos,log.
        """
        
        [headerfile,delimiter]=info

        if headerfile==None:
        	formatstr=infile.readline().lstrip()
        else:
        	hf=open(headerfile,'r')
        	formatstr=hf.readline().lstrip()
        	hf.close()
        
        formatstr=formatstr.replace('#','')
        formatstr=formatstr.replace('"','')

        header=formatstr.split(delimiter)
        header[-1]=header[-1].rstrip('\n')
        nparams=len(header)
        llines=[]
        import re
        dec=re.compile(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$|^inf$')
        
        for line_number,line in enumerate(infile):
            sline=line.split(delimiter)
            if sline[-1] == '\n':
                del(sline[-1])
            proceed=True
            if len(sline)<1:
                print 'Ignoring empty line in input file: %s'%(sline)
                proceed=False
            elif len(sline)!=nparams:
                sys.stderr.write('WARNING: Malformed row %i, read %i elements but there is meant to be %i\n'%(line_number,len(sline),nparams))
                proceed=False

            for elemn,st in enumerate(sline):
                s=st.replace('\n','')
                if dec.search(s) is None:
                    print 'Warning! Ignoring non-numeric data after the header: %s. Row = %i,Element=%i'%(s,line_number,elemn)
                    proceed=False
                elif s is '\n':
                    proceed=False

            if proceed:
                llines.append(map(float,sline))	

        flines=np.array(llines)

        if not flines.any():
            raise RuntimeError("ERROR: no lines read in!")	
           

        for i in range(0,len(header)):
            if header[i].lower().find('log')!=-1 and header[i].lower() not in logParams:
                print 'exponentiating %s'%(header[i])

                flines[:,i]=np.exp(flines[:,i])

                header[i]=header[i].replace('log','')
            if header[i].lower().find('sin')!=-1:
                print 'asining %s'%(header[i])
                flines[:,i]=np.arcsin(flines[:,i])
                header[i]=header[i].replace('sin','')
            if header[i].lower().find('cos')!=-1:
                print 'acosing %s'%(header[i])
                flines[:,i]=np.arccos(flines[:,i])
                header[i]=header[i].replace('cos','')
            header[i]=header[i].replace('(','')
            header[i]=header[i].replace(')','')
        print 'Read columns %s'%(str(header))
        return header,flines
#

def parse_converge_output_section(fo):
        result={}
        lines=fo.split('\n')
        chain_line=False
        for line in lines:
            
            if '[1]' in line:
                key=line.replace('[1]','').strip(' ').strip('"')
                result[key]={}
                out=result[key]
                continue
            if result is not {}:
                if 'chain' in line:
                    chain_line=True
                    continue
                if chain_line:
                    chain_line=False
                    key=line.strip('"').split()[1]
                    out[key]=[]
                    out2=out[key]
                else:
                    try:
                        newline=line.strip('"').split()
                        if newline is not []:
                            out2.append(line.strip('"').split())
                    except:
                        pass

        return result
#

def vo_nest2pos(nsresource,Nlive=None):
    """
    Parse a VO Table RESOURCE containing nested sampling output and
    return a VOTable TABLE element with posterior samples in it.
    This can be added to an existing tree by the user.
    Nlive will be read from the nsresource, unless specified
    """
    from xml.etree import ElementTree as ET
    import copy
    from math import log, exp
    xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
    try:
        register_namespace=ET.register_namespace
    except AttributeError:
        def register_namespace(prefix,uri):
            ET._namespace_map[uri]=prefix
    register_namespace('vot',xmlns)
    
    postable=ET.Element("{%s}TABLE"%(xmlns),attrib={'name':'Posterior Samples','utype':'lalinference:results:posteriorsamples'})
    i=0
    nstables=[resource for resource in nsresource.findall("./{%s}TABLE"%(xmlns)) if resource.get("utype")=="lalinference:results:nestedsamples"]

    nstable=nstables[0]
    if Nlive is None:
        runstateResource = [resource for resource in nsresource.findall("./{%s}RESOURCE"%(xmlns)) if resource.get("utype")=="lalinference:state"][0]
        algTable = [table for table in runstateResource.findall("./{%s}TABLE"%(xmlns)) if table.get("utype")=="lalinference:state:algorithmparams"][0]
        Nlive = int ([param for param in algTable.findall("./{%s}PARAM"%(xmlns)) if param.get("name")=='Nlive'][0].get('value'))
        print 'Found Nlive %i'%(Nlive)
    if Nlive is None:
        raise RuntimeError("Cannot find number of live points in XML table, please specify")
    logLcol = None
    for fieldnode in nstable.findall('./{%s}FIELD'%xmlns):
        if fieldnode.get('name') == 'logL':
            logLcol=i
        i=i+1
        postable.append(copy.deepcopy(fieldnode))
    for paramnode in nstable.findall('./{%s}PARAM'%(xmlns)):
        postable.append(copy.deepcopy(paramnode))
    if logLcol is None:
        RuntimeError("Unable to find logL column")
    posdataNode=ET.Element("{%s}DATA"%(xmlns))
    postabledataNode=ET.Element("{%s}TABLEDATA"%(xmlns))
    postable.append(posdataNode)
    posdataNode.append(postabledataNode)
    nstabledata=nstable.find('./{%s}DATA/{%s}TABLEDATA'%(xmlns,xmlns))
    logw=log(1.0 - exp(-1.0/float(Nlive)))
    weights=[]
    for row in nstabledata:
        logL=float(row[logLcol].text)
        weights.append(logL-logw)
        logw=logw-1.0/float(Nlive)
    mw=max(weights)
    weights = [w - mw for w in weights]
    for (row,weight) in zip(nstabledata,weights):
        if weight > log(random.random()):
            postabledataNode.append(copy.deepcopy(row))
    return postable

xmlns='http://www.ivoa.net/xml/VOTable/v1.1'

class VOT2HTML:
    def __init__(self):
        self.html=htmlSection("VOTable information")
        self.skiptable=0
        
    def start(self,tag,attrib):
        if tag=='{%s}TABLE'%(xmlns):
            if attrib['utype']=='lalinference:results:nestedsamples'\
            or attrib['utype']=='lalinference:results:posteriorsamples':
                self.skiptable=1
            else:
                self.skiptable=0
            self.tableouter=htmlChunk('div')
            self.tableouter.h2(attrib['name'])
            try:
                self.tableouter.p(attrib['utype'])
            except KeyError:
                pass
            self.fixedparams=htmlChunk('table',attrib={'class':'statstable'},parent=self.tableouter)
            self.table=htmlChunk('table',attrib={'class':'statstable'},parent=self.tableouter)
            self.tabheader=htmlChunk('tr',parent=self.table)
        if tag=='{%s}FIELD'%(xmlns):
            self.field=htmlChunk('th',{'name':attrib['name']},parent=self.tabheader)
        if tag=='{%s}TR'%(xmlns):
            self.tabrow=htmlChunk('tr',parent=self.table)
        if tag=='{%s}TD'%(xmlns):
            self.td=htmlChunk('td',parent=self.tabrow)
        if tag=='{%s}PARAM'%(xmlns):
            pnode=htmlChunk('tr',parent=self.fixedparams)
            namenode=htmlChunk('td',parent=pnode)
            namenode.p(attrib['name'])
            valnode=htmlChunk('td',parent=pnode)
            valnode.p(attrib['value'])
        
    def end(self,tag):
        if tag=='{%s}TABLE'%(xmlns):
            if not self.skiptable:
                self.html.append(self.tableouter._html)
        if tag=='{%s}FIELD'%(xmlns):
            self.field.p(self.data)
        if tag=='{%s}TD'%(xmlns):
            self.td.p(self.data)

    def data(self,data):
        self.data=data
  
    def close(self):
        return self.html.toprettyxml()
