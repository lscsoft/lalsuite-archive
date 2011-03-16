# -*- coding: utf-8 -*-
#
#       bayespputils.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#
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
from math import ceil,floor,sqrt,pi as pi_constant
import xml
from xml.dom import minidom

#related third party imports
import numpy as np
from matplotlib import pyplot as plt,cm as mpl_cm,lines as mpl_lines
from scipy import stats

import random

try:
    from xml.etree.cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser
except ImportError:
    #Python < 2.5
    from cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser

#local application/library specific imports
import pylal
from pylal import git_version
#C extensions
from _bayespputils import _skyhist_cart,_calculate_confidence_levels,_burnin

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

#===============================================================================
# Constants
#===============================================================================

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
    Function mapping (mchirp,eta)->m1; m1>m2 .
    """
    (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
    return mass1
def _inj_m2(inj):
    """
    Function mapping (mchirp,eta)->m2; m1>m2 .
    """
    (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
    return mass2

def _inj_longitude(inj):
    """
    Map the value of the longitude found in inj to an interval [0,2*pi).
    """
    if inj.longitude>2*pi_constant or inj.longitude<0.0:
        maplong=2*pi_constant*(((float(inj.longitude))/(2*pi_constant)) - floor(((float(inj.longitude))/(2*pi_constant))))
        print "Warning: Injected longitude/ra (%s) is not within [0,2\pi)! Angles are assumed to be in radians so this will be mapped to [0,2\pi). Mapped value is: %s."%(str(inj.longitude),str(maplong))
        return maplong
    else:
        return inj.longitude

def _inj_a1(inj):
    x = inj.spin1x
    y = inj.spin1y
    z = inj.spin1z
    return sqrt(x*x + y*y + z*z)

def _inj_a2(inj):
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
    Sx  = inj.spin1x
    Sy  = inj.spin1y
    Sz  = inj.spin1z
    Lnx = np.arcsin(inj.inclination)
    Lny = 0.0
    Lnz = np.arccos(inj.inclination)
    if Sx == 0.0 and Sy == 0.0 and Sz == 0.0:
        return None
    else:
        return np.arccos(Sx*Lnx + Sy*Lny + Sz*Lnz)

def _inj_tilt2(inj):
    Sx  = inj.spin2x
    Sy  = inj.spin2y
    Sz  = inj.spin2z
    Lnx = np.arcsin(inj.inclination)
    Lny = 0.0
    Lnz = np.arccos(inj.inclination)
    if Sx == 0.0 and Sy == 0.0 and Sz == 0.0:
        return None
    else:
        return np.arccos(Sx*Lnx + Sy*Lny + Sz*Lnz)

def _inj_thetas(inj):
    mtsun = 4.92549095e-06  #Msol in seconds
    f_inj = 40.0            #Assume starting frequency is 40Hz TODO: not assume

    Lmag = np.power(inj.mchirp,5.0/3.0) / np.power(pi_constant * mtsun * f_inj,1.0/3.0)
    Lx = Lmag * np.arcsin(inj.inclination)
    Ly = 0.0
    Lz = Lmag * np.arccos(inj.inclination)
    
    S1x  = inj.m1*inj.m1*inj.spin1x
    S1y  = inj.m1*inj.m1*inj.spin1y
    S1z  = inj.m1*inj.m1*inj.spin1z
    
    S2x  = inj.m2*inj.m2*inj.spin2x
    S2y  = inj.m2*inj.m2*inj.spin2y
    S2z  = inj.m2*inj.m2*inj.spin2z

    Jx = Lx + S1x + S2x
    Jy = Ly + S1y + S2y
    Jz = Lz + S1z + S2z
    Jmag = np.sqrt(Jx*Jx + Jy*Jy + Jz*Jz)

    return np.arccos(Jz/Jmag)
    
def _inj_beta(inj):
    mtsun = 4.92549095e-06  #Msol in seconds
    f_inj = 40.0            #Assume starting frequency is 40Hz TODO: not assume

    Lmag = np.power(inj.mchirp,5.0/3.0) / np.power(pi_constant * mtsun * f_inj,1.0/3.0)
    Lx = Lmag * np.arcsin(inj.inclination)
    Ly = 0.0
    Lz = Lmag * np.arccos(inj.inclination)
    
    S1x  = inj.m1*inj.m1*inj.spin1x
    S1y  = inj.m1*inj.m1*inj.spin1y
    S1z  = inj.m1*inj.m1*inj.spin1z
    
    S2x  = inj.m2*inj.m2*inj.spin2x
    S2y  = inj.m2*inj.m2*inj.spin2y
    S2z  = inj.m2*inj.m2*inj.spin2z

    Jx = Lx + S1x + S2x
    Jy = Ly + S1y + S2y
    Jz = Lz + S1z + S2z
    Jmag = np.sqrt(Jx*Jx + Jy*Jy + Jz*Jz)

    return np.arccos((Jx*Lx + Jy*Ly + Jz*Lz)/(Jmag*Lmag))


#===============================================================================
# Class definitions
#===============================================================================

class OneDPosterior(object):
    """
    A data structure for a single chain of posterior samples.
    """
    def __init__(self,name,posterior_samples,injected_value=None,prior=None):
        """
        Constructor.

        @param name: A literal string name for the parameter.
        @param posterior_samples: A 1D array of the samples.
        @keyword injected_value: The injected or real value of the parameter.
        @keyword prior: The prior value corresponding to each sample.
        """
        self.__name=name
        self.__posterior_samples=np.array(posterior_samples)

        self.__injval=injected_value
        self.__prior=prior

        return

    @property
    def name(self):
        """
        Return the literal name of the parameter.
        """
        return self.__name

    @property
    def mean(self):
        """
        Calculate the arithmetic mean of the 1D samples.
        """
        return np.mean(self.__posterior_samples)

    @property
    def median(self):
        """
        Find the median value of the 1D samples.
        """
        return np.median(self.__posterior_samples)

    @property
    def stdev(self):
        """
        Return the standard deviation of the 1D samples.
        """
        return sqrt(np.var(self.__posterior_samples))

    @property
    def stacc(self):
        """
         The 'standard accuracy statistic' (stacc) - a standard deviant incorporating
        information about the accuracy of the waveform recovery. Defined as
        the mean of the sum of the squared differences between the points
        in the PDF (x_i - sampled according to the posterior) and the
        true value (x_{\rm true}).
        So for a marginalized one-dimensional PDF:
        stacc = \sqrt{\frac{1}{N}\sum_{i=1}^N (x_i-x_{\rm true})2}
        """
        if self.__injval is None:
            return None
        else:
            return sqrt(np.var(self.__posterior_samples)+pow((np.mean(self.__posterior_samples)-self.__injval),2) )

    @property
    def injval(self):
        """
        Return the injected value set at construction . If no value was set
        will return None .
        """
        return self.__injval

    #@injval.setter #Python 2.6+
    def set_injval(self,new_injval):
        """
        Set the injected/real value of the parameter.

        @param new_injval: The injected/real value to set.
        """

        self.__injval=new_injval

    @property
    def samples(self):
        """
        Return a 1D numpy.array of the samples.
        """
        return self.__posterior_samples

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from posterior, analagous to numpy.delete but opperates in place.

        @param samples: A list of the indexes of the samples to remove.
        """
        self.__posterior_samples=np.delete(self.__posterior_samples,samples).reshape(-1,1)

    @property
    def gaussian_kde(self):
        """
        Return a gaussian kde of the samples.
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
    def __init__(self,commonResultsFormatData,SimInspiralTableEntry=None,name=None,description=None):
        """
        Constructor.

        @param commonResultsFormatData: A 2D array containing the posterior
            samples and related data. The samples chains form the columns.

        """
        common_output_table_header,common_output_table_raw =commonResultsFormatData
        self._posterior={}
        self._injection=SimInspiralTableEntry
        self._loglaliases=['logl','logL','likelihood','posterior']

        common_output_table_header=[i.lower() for i in common_output_table_header]
        
        for one_d_posterior_samples,param_name in zip(np.hsplit(common_output_table_raw,common_output_table_raw.shape[1]),common_output_table_header):
            
            self._posterior[param_name]=OneDPosterior(param_name.lower(),one_d_posterior_samples,injected_value=self._getinjpar(param_name))

        if 'mchirp' in common_output_table_header and 'eta' in common_output_table_header \
        and (not 'm1' in common_output_table_header) and (not 'm2' in common_output_table_header):
            try:
                print 'Inferring m1 and m2 from mchirp and eta'
                (m1,m2)=mc2ms(self._posterior['mchirp'].samples, self._posterior['eta'].samples)
                self._posterior['m1']=OneDPosterior('m1',m1,injected_value=self._getinjpar('m1'))
                self._posterior['m2']=OneDPosterior('m2',m2,injected_value=self._getinjpar('m2'))
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

        return Posterior((names,bootstrapSamples),self._injection)

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from all OneDPosteriors.

        @param samples: The indixes of the samples to remove.
        """
        for name,pos in self:
            pos.delete_samples_by_idx(samples)

    @property
    def injection(self):
        """
        Return the injected values .
        """

        return self._injection

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

        @param injection: A SimInspiralTable row object.
        """
        if injection is not None:
            self._injection=injection
            for name,onepos in self:
                new_injval=self._getinjpar(name)
                if new_injval is not None:
                    self[name].set_injval(new_injval)


    _injXMLFuncMap={
                        'mchirp':lambda inj:inj.mchirp,
                        'mc':lambda inj:inj.mchirp,
                        'mass1':_inj_m1,
                        'm1':_inj_m1,
                        'mass2':_inj_m2,
                        'm2':_inj_m2,
                        'eta':lambda inj:inj.eta,
                        'time': lambda inj:float(inj.get_end()),
                        'end_time': lambda inj:float(inj.get_end()),
                        'phi0':lambda inj:inj.phi0,
                        'phi_orb': lambda inj: inj.phi0,
                        'dist':lambda inj:inj.distance,
                        'distance':lambda inj:inj.distance,
                        'ra':_inj_longitude,
                        'long':_inj_longitude,
                        'longitude':_inj_longitude,
                        'dec':lambda inj:inj.latitude,
                        'lat':lambda inj:inj.latitude,
                        'latitude':lambda inj:inj.latitude,
                        'psi': lambda inj: inj.polarization,
                        'iota':lambda inj: inj.inclination,
                        'inclination': lambda inj: inj.inclination,
                        'spinchi': lambda inj: (inj.spin1z + inj.spin2z) + sqrt(1-4*inj.eta)*(inj.spin1z - spin2z),
                        'a1':_inj_a1,
                        'a2':_inj_a2,
                        'theta1':_inj_theta1,
                        'theta2':_inj_theta2,
                        'phi1':_inj_phi1,
                        'phi2':_inj_phi2,
                        'tilt1':_inj_tilt1,
                        'tilt2':_inj_tilt2,
                        'cos(iota)': lambda inj: np.cos(inj.inclination),
                        'theta_s':_inj_thetas,
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
            yield ParameterSample(sample_array, header, header)
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

    def _average_posterior(self, samples, post_name):
        ap = 0.0
        for samp in samples:
            ap = ap + samp[post_name]
        return ap / len(samples)

    def _average_posterior_like_prior(self, samples, logl_name, prior_name):
        ap = 0.0
        for samp in samples:
            ap += np.exp(samp[logl_name])*samp[prior_name]
        return ap / len(samples)

    def di_evidence(self, boxing=64):
        """
        Returns the direct-integration evidence for the posterior
        samples.
        """
        allowed_coord_names=["a1", "phi1", "theta1", "a2", "phi2", "theta2",
                             "iota", "psi", "ra", "dec",
                             "phi_orb", "phi0", "dist", "time", "mc", "mchirp", "eta"]
        samples,header=self.samples()
        header=header.split()
        coord_names=[name for name in allowed_coord_names if name in header]
        coordinatized_samples=[ParameterSample(row, header, coord_names) for row in samples]
        tree=KDTree(coordinatized_samples)

        if "post" in header:
            return tree.integrate(lambda samps: self._average_posterior(samps, "post"), boxing)
        elif "posterior" in header:
            return tree.integrate(lambda samps: self._average_posterior(samps, "posterior"), boxing)
        elif "prior" in header and "logl" in header:
            return tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "logl", "prior"), boxing)
        elif "prior" in header and "likelihood" in header:
            return tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "likelihood", "prior"), boxing)
        else:
            raise RuntimeError("could not find 'post', 'posterior', 'logl' and 'prior', or 'likelihood' and 'prior' columns in output to compute direct integration evidence")



    def harmonic_mean_evidence(self):
        """
        Returns the harmonic mean evidence for the set of posterior
        samples.
        """
        return 1/np.mean(1/np.exp(self._logL))

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

            return_val+=self._print_table_row(name,[maxL,stdev,mean,median,stacc,injval])

        return_val+='</table>'

        parser=XMLParser()
        parser.feed(return_val)
        Estr=parser.close()

        elem=Estr
        rough_string = tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return_val=reparsed.toprettyxml(indent="  ")

        return return_val

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

    def integrate(self, f, boxing=64):
        """
        Returns the integral of f(objects) over the tree.  The
        optional boxing parameter determines how deep to descend into
        the tree before computing f.
        """
        if len(self._objects) <= boxing:
            return self.volume()*f(self._objects)
        else:
            return self._left.integrate(f, boxing) + self._right.integrate(f, boxing)

class ParameterSample(object):
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

def _calculate_sky_confidence_slow(
                                shist,
                                skypoints,
                                injbin,
                                skyres_,
                                confidence_levels,
                                lenpos):
    """
    @deprecated: This is a pure python version of the C extension function
        pylal._bayespputils._calculate_confidence_levels.
    """
    frac=0
    Nbins=0
    injectionconfidence=None
    #print "lenpos : %i"%lenpos
    #toppoints=[(None,None,None)]*lenpos
    toppoints=[]

    skyreses=[]
    lenbins=len(shist)
    range_lenbins=range(0,lenbins)
    for confidence_level in confidence_levels:
        while(frac<confidence_level):
            maxbin=0
            for i in range_lenbins:
                if shist[i]>maxbin:
                    maxbin=shist[i]
                    maxpos=i

            shist[maxpos]=0
            frac=frac+(float(maxbin)/(lenpos))

            Nbins=Nbins+1
            toppoints.append((skypoints[maxpos,0],skypoints[maxpos,1],maxpos,frac))
            if injbin is not None:
                if (injbin==maxpos):
                    injectionconfidence=frac
                    print 'Injection sky point found at confidence %f'%(frac)

        print '%f confidence region: %f square degrees'%(frac,Nbins*float(skyres_)*float(skyres_))

        skyreses.append((frac,Nbins*float(skyres_)*float(skyres_)))
        toppoints=toppoints[:Nbins]
    return injectionconfidence,toppoints,skyreses

def _greedy_bin(greedyHist,greedyPoints,injection_bin_index,bin_size,Nsamples,confidence_levels):
    """
    An interal function representing the common, dimensionally-independent part of the
    greedy binning algorithms.
    """

    #Now call confidence level C extension function to determine top-ranked pixels
    (injectionconfidence,toppoints)=_calculate_confidence_levels(
                                                                    greedyHist,
                                                                    greedyPoints,
                                                                    injection_bin_index,
                                                                    bin_size,
                                                                    Nsamples
                                                                    )

    #Determine interval/area contained within given confidence intervals
    nBins=0
    confidence_levels.sort()
    reses={}
    toppoints=np.array(toppoints)
    for printcl in confidence_levels:
        nBins=1
        #Start at top of list of ranked pixels...
        accl=toppoints[0,3]

        #Loop over next significant pixels and their confidence levels

        while accl<printcl and nBins<=len(toppoints):
            nBins=nBins+1
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
            print par1_binNumber,par2_binNumber,par1pos_Nbins,par2pos_Nbins,par1_binNumber+par2_binNumber*par1pos_Nbins,par1_samp,par1pos_min,par1_bin,par1_samp,par2pos_min,par2_bin
            exit(1)
    #Call greedy bins routine
    toppoints,injection_cl,reses,injection_area=\
                                _greedy_bin(
                                                greedyHist,
                                                greedyPoints,
                                                injbin,
                                                float(sqrt(par1_bin*par2_bin)),
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

    return _greedy_bin(shist,skypoints,injbin,float(skyres),len(skypos),confidence_levels)


def plot_sky_map(top_ranked_pixels,outdir):
    """
    Plots a sky map using the Mollweide projection in the Basemap package.

    @param top_ranked_pixels: the top-ranked sky pixels as determined by greedy_bin_sky.

    @param outdir: Output directory in which to save skymap.png image.
    """
    from mpl_toolkits.basemap import Basemap

    np.seterr(under='ignore')

    myfig=plt.figure()
    plt.clf()
    m=Basemap(projection='moll',lon_0=180.0,lat_0=0.0)
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
    np.savetxt(
               os.path.join(outdir,'ranked_sky_pixels.dat'),
               np.column_stack(
                               [
                                np.asarray(top_ranked_pixels)[:,0:1],
                                np.asarray(top_ranked_pixels)[:,1],
                                np.asarray(top_ranked_pixels)[:,3]
                                ]
                               )
               )

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

def plot_one_param_pdf_kde(fig,onedpos):

    from scipy import seterr as sp_seterr

    np.seterr(under='ignore')
    sp_seterr(under='ignore')
    pos_samps=onedpos.samples
    gkde=onedpos.gaussian_kde

    ind=np.linspace(np.min(pos_samps),np.max(pos_samps),101)
    kdepdf=gkde.evaluate(ind)
    plt.plot(ind,kdepdf)

    return

def plot_one_param_pdf_line_hist(fig,pos_samps):
    plt.hist(pos_samps,kdepdf)


def plot_one_param_pdf(posterior,plot1DParams):
    """
    Plots a 1D histogram and (gaussian) kernel density estimate of the
    distribution of posterior samples for a given parameter.

    @param posterior: an instance of the Posterior class.

    @param plot1DParams: a dict; {paramName:Nbins}

    """

    param=plot1DParams.keys()[0].lower()
    histbins=plot1DParams.values()[0]

    pos_samps=posterior[param].samples
    injpar=posterior[param].injval

    myfig=plt.figure(figsize=(4,3.5),dpi=200)
    axes=plt.Axes(myfig,[0.2, 0.2, 0.7,0.7])
    myfig.add_axes(axes)

    (n, bins, patches)=plt.hist(pos_samps,histbins,normed='true')
    histbinSize=bins[1]-bins[0]

    plot_one_param_pdf_kde(myfig,posterior[param])

    rbins=None

    if injpar:
        if min(pos_samps)<injpar and max(pos_samps)>injpar:
            plt.axvline(injpar, color='r', linestyle='-.')

            #rkde=gkde.integrate_box_1d(min(pos[:,i]),getinjpar(injection,i))
            #print "r of injected value of %s (kde) = %f"%(param,rkde)

            #Find which bin the true value is in
            bins_to_inj=(injpar-bins[0])/histbinSize
            injbinh=int(floor(bins_to_inj))
            injbin_frac=bins_to_inj-float(injbinh)

            #Integrate over the bins
            rbins=(sum(n[0:injbinh-1])+injbin_frac*n[injbinh])*histbinSize

    #
    plt.grid()
    plt.xlabel(param)
    plt.ylabel('Probability Density')

    # For RA and dec set custom labels and for RA reverse
    if(param.lower()=='ra' or param.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        plt.xlim(xmax,xmin)
    if(param.lower()=='ra' or param.lower()=='rightascension'):
        locs, ticks = plt.xticks()
        strticks=map(getRAString,locs)
        plt.xticks(locs,strticks,rotation=45)
    if(param.lower()=='dec' or param.lower()=='declination'):
        locs, ticks = plt.xticks()
        strticks=map(getDecString,locs)
        plt.xticks(locs,strticks,rotation=45)

    return rbins,myfig#,rkde
#

def getRAString(radians):
    hours = floor(radians*(12.0/pi_constant))
    rem = radians-hours*(pi_constant/12.0)
    mins = floor(rem*((12*60)/pi_constant))
    rem = rem - mins*(pi_constant/(12*60))
    secs = rem*(12*3600/pi_constant)
    return '$%i\mathrm{h}%i^{\'}%2.0f^{\'\'}$'%(hours,mins,secs)

def getDecString(radians):
    if(radians<0):
        round = ceil
        sign=-1
    else:
        round = floor
        sign=+1
    deg = round(radians*(180.0/pi_constant))
    rem = radians - deg*(pi_constant/180.0)
    mins = round(rem*((180.0*60.0)/pi_constant))
    rem = rem - mins*(pi_constant/(180.0*60.0))
    secs = rem * (180.0*60.0*60.0)/pi_constant
    return '$%i^\circ%i^{\'}%2.0f^{\'\'}$'%(deg,sign*mins,sign*secs)

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

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    myfig=plt.figure(1,figsize=(6,4),dpi=200)
    myfig.add_axes(plt.Axes(myfig,[0.2,0.25,0.75,0.7]))
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
        plt.plot([par_injvalue1],[par_injvalue2],'go',scalex=False,scaley=False)

    plt.xlabel(par1_name)
    plt.ylabel(par2_name)
    plt.grid()

    # For RA and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            xmin,xmax=plt.xlim()
            plt.xlim(xmax,xmin)
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            locs, ticks = plt.xticks()
            strticks=map(getRAString,locs)
            plt.xticks(locs,strticks,rotation=45)
    if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
            locs, ticks = plt.xticks()
            strticks=map(getDecString,locs)
            plt.xticks(locs,strticks,rotation=45)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        ymin,ymax=plt.ylim()
        plt.ylim(ymax,ymin)
    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        locs, ticks = plt.yticks()
        strticks=map(getRAString,locs)
        plt.yticks(locs,strticks)
    if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
        locs, ticks = plt.yticks()
        strticks=map(getDecString,locs)
        plt.yticks(locs,strticks)

    return myfig
#

def get_inj_by_time(injections,time):
    """
    Filter injections to find the injection with end time given by time +/- 0.1s
    """
    import itertools
    injection = itertools.ifilter(lambda a: abs(float(a.get_end()) - time) < 0.1, injections).next()
    return injection

def plot_two_param_greedy_bins_contour(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,line_styles=__default_line_styles,figsize=(7,6),dpi=250,figposition=[0.2,0.2,0.48,0.75]):
    """
    Plots the confidence level contours as determined by the 2-parameter
    greedy binning algorithm.

    @param posteriors_by_name: A dict containing Posterior instances referenced by some id.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}

    @param confidence_levels: a list of the required confidence levels to plot on the contour map.

    @param colors_by_name: A dict of colors cross-referenced to the above Posterior ids.

    """

    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()

    fig.add_axes(figposition)

    #This fixes the precedence of line styles in the plot
    if len(line_styles)<len(confidence_levels):
        print "Error: Need as many or more line styles to choose from as confidence levels to plot!"
        exit(0)

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

        a=np.squeeze(posterior[par1_name].samples)
        b=np.squeeze(posterior[par2_name].samples)

        #Create 2D bin array
        par1pos_min=a.min()
        par2pos_min=b.min()

        par1pos_max=a.max()
        par2pos_max=b.max()

        par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
        par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

        H, xedges, yedges = np.histogram2d(a,b, bins=(par1pos_Nbins, par2pos_Nbins),normed=True)

        extent = [xedges[0], yedges[-1], xedges[-1], xedges[0]]

        temp=np.copy(H)
        temp=temp.ravel()
        confidence_levels.sort()
        Hsum=0
        Hlasts=[]
        for cl in confidence_levels:
            while float(Hsum/np.sum(H))<cl:
                ind = np.argsort(temp)
                max_i=ind[-1:]
                val = temp[max_i]
                Hlast=val[0]
                Hsum+=val
                temp[max_i]=0
            Hlasts.append(Hlast)

        CS=plt.contour(yedges[:-1],xedges[:-1],H,Hlasts,colors=[colors_by_name[name]],linestyles=line_styles)
        plt.grid()

        CSlst.append(CS)


    plt.title("%s-%s confidence contours (greedy binning)"%(par1_name,par2_name)) # add a title
    plt.xlabel(par2_name)
    plt.ylabel(par1_name)

    if len(name_list)!=len(CSlst):
        print "Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name."
        exit(0)
    full_name_list=[]
    dummy_lines=[]

    for plot_name in name_list:
        full_name_list.append(plot_name)
    for ls_,cl in zip(line_styles[0:len(confidence_levels)],confidence_levels):
        dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),ls=ls_,color='k'))
        full_name_list.append('%s%%'%str(int(cl*100)))

    fig_actor_lst = [cs.collections[0] for cs in CSlst]

    fig_actor_lst.extend(dummy_lines)

    twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')

    for text in twodcontour_legend.get_texts():
        text.set_fontsize('small')


    # For ra and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            ymin,ymax=plt.ylim()
            if(ymin<0.0): ylim=0.0
            if(ymax>2.0*pi_constant): ymax=2.0*pi_constant
            plt.ylim(ymax,ymin)
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            locs, ticks = plt.yticks()
            strticks=map(getRAString,locs)
            plt.yticks(locs,strticks)
    if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
            locs, ticks = plt.yticks()
            strticks=map(getDecString,locs)
            plt.yticks(locs,strticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        locs, ticks = plt.xticks()
        strticks=map(getRAString,locs)
        plt.xticks(locs,strticks,rotation=45)
    if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
        locs, ticks = plt.xticks()
        strticks=map(getDecString,locs)
        plt.xticks(locs,strticks,rotation=45)

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

    #Extract injection information
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval

    myfig=plt.figure(1,figsize=(10,8),dpi=300)
    myfig.add_axes([0.2,0.2,0.8,0.8])
    plt.clf()
    plt.xlabel(par2_name)
    plt.ylabel(par1_name)

    bins=(100,100)

    H, xedges, yedges = np.histogram2d(a,b, bins,normed=False)

    #Replace H with greedy bin confidence levels at each pixel...
    temp=np.copy(H)
    temp=temp.flatten()

    Hsum=0
    Hsum_actual=np.sum(H)

    while Hsum<Hsum_actual:
        ind = np.argsort(temp)
        max_i=ind[-1:]
        val = temp[max_i]
        Hsum+=int(val)
        temp[max_i]=0

        #print Hsum,Hsum_actual
        H.flat[max_i]=1-float(Hsum)/float(Hsum_actual)

    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    plt.imshow(np.flipud(H), aspect='auto', extent=extent, interpolation='nearest')
    plt.gca().autoscale_view()
    plt.colorbar()

    if par1_injvalue is not None and par2_injvalue is not None:
        plt.plot([par1_injvalue],[par2_injvalue],'go',scalex=False,scaley=False)


    # For RA and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            ymin,ymax=plt.ylim()
            plt.ylim(ymax,ymin)
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
            locs, ticks = plt.yticks()
            strticks=map(getRAString,locs)
            plt.yticks(locs,strticks)
    if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
            locs, ticks = plt.yticks()
            strticks=map(getDecString,locs)
            plt.yticks(locs,strticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin)<0.0: xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        locs, ticks = plt.xticks()
        strticks=map(getRAString,locs)
        plt.xticks(locs,strticks,rotation=45)
    if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
        locs, ticks = plt.xticks()
        strticks=map(getDecString,locs)
        plt.xticks(locs,strticks,rotation=45)

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

    toppoints,injectionconfidence,reses,injection_area=_greedy_bin(greedyHist,greedyPoints,injbin,float(sqrt(par_bin*par_bin)),int(len(par_samps)),confidence_levels)
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

    def parse(self,files,**kwargs):
        """
        Parse files.
        """
        return self._parser(files,**kwargs)

    def _infmcmc_to_pos(self,files,deltaLogL=None,nDownsample=None,oldMassConvention=False,**kwargs):
        """
        Parser for lalinference_mcmcmpi output.
        """
        logLThreshold=-1e200 # Really small?
        if not (deltaLogL is None):
            logLThreshold=self._find_max_logL(files) - deltaLogL
            print "Eliminating any samples before log(Post) = ", logLThreshold
        nskip=1
        if not (nDownsample is None):
            nskip=self._find_ndownsample(files, logLThreshold, nDownsample)
            print "Downsampling by a factor of ", nskip, " to achieve approximately ", nDownsample, " posterior samples"
        postName="posterior_samples.dat"
        outfile=open(postName, 'w')
        try:
            self._infmcmc_output_posterior_samples(files, outfile, logLThreshold, nskip, oldMassConvention)
        finally:
            outfile.close()
        return self._common_to_pos(open(postName,'r'))


    def _infmcmc_output_posterior_samples(self, files, outfile, logLThreshold, nskip=1, oldMassConvention=False):
        """
        Concatenate all the samples from the given files into outfile.
        For each file, only those samples past the point where the
        log(L) > logLThreshold are concatenated.
        """
        print oldMassConvention
        allowedCols=["cycle", "logl", "logpost", "logprior",
                     "a1", "theta1", "phi1",
                     "a2", "theta2", "phi2",
                     "mc", "eta", "time",
                     "phi_orb", "iota", "psi",
                     "ra", "dec",
                     "dist"]
        nRead=0
        outputHeader=False
        for infilename,i in zip(files,range(1,len(files)+1)):
            infile=open(infilename,'r')
            try:
                print "Processing file %s to posterior_samples.dat"%infilename
                header=self._clear_infmcmc_header(infile)
                # Remove unwanted columns, and accound for 1 <--> 2 reversal of convention in lalinference.
                if oldMassConvention:
                    header=[self._swaplabel12(label) for label in header if label in allowedCols]
                else:
                    header=[label for label in header if label in allowedCols]
                if not outputHeader:
                    for label in header:
                        outfile.write(label)
                        outfile.write(" ")
                    outfile.write("chain")
                    outfile.write("\n")
                    outputHeader=header
                loglindex=header.index("logpost")
                output=False
                for line in infile:
                    line=line.lstrip()
                    lineParams=line.split()
                    logL=float(lineParams[loglindex])
                    if logL >= logLThreshold:
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
                            outfile.write(str(i))
                            outfile.write("\n")
                        nRead=nRead+1
            finally:
                infile.close()

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
                header=self._clear_infmcmc_header(infile)
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

    def _find_ndownsample(self, files, logLthreshold, nDownsample):
        """
        Given a list of files, and threshold value, and a desired
        number of outputs posterior samples, return the skip number to
        achieve the desired number of posterior samples.
        """
        ntot=0
        for inpname in files:
            infile=open(inpname, 'r')
            try:
                header=self._clear_infmcmc_header(infile)
                loglindex=header.index("logpost")
                burnedIn=False
                for line in infile:
                    line=line.lstrip().split()
                    logL=float(line[loglindex])
                    if logL > logLthreshold:
                        burnedIn=True
                    if burnedIn:
                        ntot=ntot+1
            finally:
                infile.close()
        if ntot < nDownsample:
            return 1
        else:
            return floor(ntot/nDownsample)

    def _clear_infmcmc_header(self, infile):
        """
        Reads past the header information from the
        lalinference_mcmcmpi file given, returning the common output
        header information.
        """
        for line in infile:
            headers=line.lstrip().lower().split()
            if len(headers) is 0:
                continue
            if "cycle" in headers[0]:
                break
        else:
            raise RuntimeError("couldn't find line beginning with 'cycle' in LALInferenceMCMC input")
        return headers


    def _mcmc_burnin_to_pos(self,files,spin=False,deltaLogL=None):
        """
        Parser for SPINspiral output .
        """
        raise NotImplementedError
        if deltaLogL is not None:
            pos,bayesfactor=burnin(data,spin,deltaLogL,"posterior_samples.dat")
            return self._common_to_pos(open("posterior_samples.dat",'r'))

    def _ns_to_pos(self,files,Nlive=None,xflag=False):
        """
        Parser for nested sampling output.
        """
        try:
            from lalapps.combine_evidence import combine_evidence
        except ImportError:
            print "Need lalapps.combine_evidence to convert nested sampling output!"
            exit(1)

        if Nlive is None:
            print "Need to specify number of live points in positional arguments of parse!"
            exit(1)

        pos,d_all,totalBayes,ZnoiseTotal=combine_evidence(files,False,Nlive)

        posfilename='posterior_samples.dat'
        posfile=open(posfilename,'w')
        posfile.write('mchirp \t eta \t time \t phi0 \t dist \t RA \t dec \t psi \t iota \t likelihood \n')
        for row in pos:
            for i in row:
                posfile.write('%f\t' %(i))
            posfile.write('\n')
        posfile.close()

        posfile=open(posfilename,'r')
        return_val=self._common_to_pos(posfile)
        posfile.close()

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

    def _common_to_pos(self,infile,delimiter=None):
        """
        Parse a file in the 'common format' and return an array of posterior
        samples and list of parameter names. Will apply inverse functions to
        columns with names containing sin,cos,log.
        """

        formatstr=infile.readline().lstrip()
        formatstr=formatstr.replace('#','')
        formatstr=formatstr.replace('"','')

        header=formatstr.split(delimiter)
        header[-1]=header[-1].rstrip('\n')

        llines=[]
        import re
        dec=re.compile(r'[^Ee+\d.-]+')
        line_count=0
        for line in infile:
            sline=line.split(delimiter)
            if sline[-1] == '\n':
                del(sline[-1])
            proceed=True
            if len(sline)<1:
                print 'Ignoring empty line in input file: %s'%(sline)
                proceed=False

            for st in sline:
                s=st.replace('\n','')
                if dec.search(s) is not None:
                    print 'Warning! Ignoring non-numeric data after the header: %s'%s
                    proceed=False
                if s is '\n':
                    proceed=False

            if proceed:
                llines.append(np.array(map(float,sline)))

        flines=np.array(llines)
        for i in range(0,len(header)):
            if header[i].lower().find('log')!=-1 and header[i].lower()!='logl':
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

def convergenceTests(posterior,gelman=True,geweke=True,geweke_frac1=0.1, geweke_frac2=0.5,heidelberger=True,effectiveSize=True,warnings=False):
    """
    This function spawns a process to run an R script to use the coda package
    to calculate the following convergence diagnostics for the chains contained in the input
    Posterior instance:
    
    Gelman-Rubin
    Geweke
    Heidelberger & Welch
    Effective sample sizes

    These can be turned on and off using the flags.
    
    Python scrapes the output produced by the R script and returns a (nested) dictionary:
    {'test':{'chainN':[['param':'test_value'],...],...},...}

    You will need R with the coda and lattice libraries for this to work . 
    """
    from subprocess import Popen,PIPE

    print "Calculating convergence diagnostics using R coda library..."

    posterior.write_to_file('tmp.dat')

    R_process=Popen(["R","-q","--vanilla","--slave"],stdin=PIPE,stdout=PIPE)
    script=convergenceTests_R%(str(gelman).upper(),str(geweke).upper(),geweke_frac1,geweke_frac2,str(heidelberger).upper(),str(effectiveSize).upper(),str(warnings).upper())
    #print script
    rp_stdout,rp_stderr=R_process.communicate(input=script)

    outdict=parse_converge_output_section(rp_stdout)
    
    #Check outdict
    if outdict is {}:
        print "No results found! Check for any errors above. You can turn on the R warning messages with warnings=True ."
        return None
    else:
        test_pass=False
        for test,test_result in outdict.items():
            
            if not test_result:
                print "Convergence test failed : %s ! Check for any errors above. You can turn on the R warning messages with warnings=True ."%test
            else:
                test_pass=True

        if test_pass:
            return outdict

        else:
            return None
    
#

#R convergenceTest script
convergenceTests_R="""
convergenceTests <- function(data,
                             gelman=TRUE,
                             geweke=TRUE, geweke.frac1=0.1, geweke.frac2=0.5,
                             heidelberger=TRUE,
                             effectiveSize=TRUE)
# The argument `data' is a matrix (or data frame)
# of (supposedly post burn-in) posterior samples.
# Rows correspond to samples, columns correspond to variables.
# If `data' contains a variable named "chain",
# this is assumed to indicate chains from independent runs.
# By default, Gelman's, Geweke's and Heidelberger & Welch's convergence
# diagnostics, and effectice sample sizes are computed.
# Geweke's diagnostic is directly converted to p-values,
# and the effective sample size is divided by the number of MCMC samples
# to give /relative/ sample sizes.
# In order to suppress computation of any of the 3 figures, supply
# an additional (e.g.)  `geweke=FALSE'  or  `geweke=F'  argument.
{
  # check whether "coda" library is installed:
  codaInstalled <- require("coda")
  if (codaInstalled) {
    # check whether a "chain" variable indicates presence of multiple chains:
    if (is.element("chain",colnames(data))) {
      chainIndicator <- data[,"chain"]
      # which chains present?:
      chainLabels <- sort(unique(data[,"chain"]))
      # how many samples of each chain?:
      chainFrequencies <- table(data[,"chain"])
      # drop "chain" column from data:
      data <- data[,colnames(data)!="chain"]
    }
    else { # (no "chain" indicator means a single chain)
      chainIndicator <- rep(1, nrow(data))
      chainLabels <- 1
      chainFrequencies <- c("1"=nrow(data))
    }

    
    # Gelman diagnostic - need at least 2 chains with at least 2 samples:
    if (gelman && (length(chainFrequencies>=2) >= 2)) {
      codaList <- mcmc.list()
      for (i in 1:sum(chainFrequencies>=2)){
        # filter out i-th chain:
        finalSamples <- which(chainIndicator==chainLabels[chainFrequencies>=2][i])
        # filter out final samples:
        finalSamples <- finalSamples[length(finalSamples)+((-(min(chainFrequencies[chainFrequencies>=2])-1)):0)]
        # add to list:
        codaList[[i]] <- mcmc(data=data[finalSamples,])
      }
      GD <- try(gelman.diag(codaList), silent=TRUE)
      rm("codaList")
      if (class(GD) != "try-error") RHatP <- c("RHatP"=Re(GD$mpsrf))
      else RHatP <- c("RHatP"=NA)
    }
    else RHatP <- c("RHatP"=NA)

    
    # Geweke diagnostic:
    if (geweke) {
      geweke.p <- geweke.z <- matrix(NA, nrow=ncol(data), ncol=length(chainLabels),
                                     dimnames=list("variable"=colnames(data), "chain"=chainLabels))
      # compute diagnostic for each variable and chain:
      for (i in 1:length(chainLabels)) {
        GD <- try(geweke.diag(mcmc(data[chainIndicator==chainLabels[i],]),
                              frac1=geweke.frac1, frac2=geweke.frac2))
        if (class(GD) != "try-error") zScores <- GD$z
        else zScores <- rep(NA, ncol(data))
        geweke.z[,i] <- zScores
      }
      # compute corresponding matrix of p-values:
      geweke.p <- (1.0 - pnorm(abs(geweke.z))) / 2.0
    }
    else geweke.p <- NA

    
    # Heidelberger & Welch diagnostic:
    if (heidelberger) {
      heidelberger.p  <- matrix(NA, nrow=ncol(data), ncol=length(chainLabels),
                                    dimnames=list("variable"=colnames(data), "chain"=chainLabels))
      # compute diagnostic for each variable and chain:
      for (i in 1:length(chainLabels)) {
        HWD <- try(heidel.diag(mcmc(data[chainIndicator==chainLabels[i],])))
        if (class(HWD) != "try-error") pvalues <- HWD[,"pvalue"]
        else pvalues <- rep(NA, ncol(data))
        heidelberger.p[,i] <- pvalues
      }
    }
    else heidelberger.p <- NA

    
    # effective sample size:
    if (effectiveSize) {
      Neff <- matrix(NA, nrow=ncol(data), ncol=length(chainLabels),
                     dimnames=list("variable"=colnames(data), "chain"=chainLabels))
      # compute N_eff for each variable and chain:
      for (i in 1:length(chainLabels)) {
        ES <- try(effectiveSize(mcmc(data[chainIndicator==chainLabels[i],])))
        if (class(ES) != "try-error") sizes <- ES
        else sizes <- rep(NA, ncol(data))
        Neff[,i] <- sizes
      }
      # normalize (to /relative/ sample sizes):
      Reff <- Neff / matrix(chainFrequencies, nrow=ncol(data), ncol=length(chainLabels), byrow=TRUE)
    }
    else Reff <- NA

    
    # assemble eventual results:
    result <- list("gelman"       = RHatP,          # "multivariate scale reduction factor", $\hat{R}^p$
                   "geweke"       = geweke.p,       # matrix of p-values
                   "heidelberger" = heidelberger.p, # matrix of p-values
                   "effectiveSize"     = Reff)           # matrix of /relative/ effective sample sizes
  }
  else {
    warning("coda library not installed.")
    result <- list("gelman"       = NA,
                   "geweke"       = NA,
                   "heidelberger" = NA,
                   "effectiveSize"     = NA)
  }
  return(result)
} 

A <- read.table("tmp.dat",header=TRUE)
result <- convergenceTests(A,gelman=%s,geweke=%s,geweke.frac1=%f, geweke.frac2=%f,heidelberger=%s,effectiveSize=%s)
for (name in names(result)) {
    print(name)
    print(result[[name]])
    
}
warnings_flag <- %s
if (warnings_flag){
    warnings()
}

"""
#
