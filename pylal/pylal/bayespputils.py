# -*- coding: utf-8 -*-
#
#       bayespputils.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
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

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

#===============================================================================
# Constants
#===============================================================================

#Pre-defined ordered list of line styles for use in matplotlib contour plots.
__default_line_styles=['solid', 'dashed', 'dashdot', 'dotted']
#Pre-defined ordered list of matplotlib colours for use in plots.
__default_color_lst=['r','b','y','g','k']
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
    def __init__(self,commonResultsFormatData,SimInspiralTableEntry=None):
        """
        Constructor.
        
        @param commonResultsFormatData: A 2D array containing the posterior 
            samples and related data. The samples chains form the columns.
                
        """
        common_output_table_header,common_output_table_raw =commonResultsFormatData
        self._posterior={}
        self._injection=SimInspiralTableEntry
        for one_d_posterior_samples,param_name in zip(np.hsplit(common_output_table_raw,common_output_table_raw.shape[1]),common_output_table_header):
            param_name=param_name.lower()
            self._posterior[param_name]=OneDPosterior(param_name.lower(),one_d_posterior_samples,injected_value=self._getinjpar(param_name))

        if 'logl' in common_output_table_header:
            try:
                self._logL=self._posterior['logl'].samples

            except KeyError:
                print "No 'logl' column in input table!"
                raise
        elif 'likelihood' in common_output_table_header:
            try:
                self._logL=self._posterior['likelihood'].samples

            except KeyError:
                print "No 'logl' column in input table!"
                raise

        elif 'post' in common_output_table_header:
            try:
                self._logL=self._posterior['post'].samples

            except KeyError:
                print "No 'post' column in input table!"
                raise

        elif 'posterior' in common_output_table_header:
            try:
                self._logL=self._posterior['posterior'].samples

            except KeyError:
                print "No 'posterior' column in input table!"
                raise

        else:
            print "No likelihood/posterior values found!"
            import sys
            sys.exit(1)

        return

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from all OneDPosteriors.
        
        @param samples: The indixes of the samples to remove.
        """
        for name,pos in self:
            pos.delete(samples)

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

    def _inj_mchirp(inj):

        return inj.mchirp

    def _inj_eta(inj):
        return inj.eta

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
                        'spinchi': lambda inj: (inj.spin1z + inj.spin2z) + sqrt(1-4*inj.eta)*(inj.spin1z - spin2z)
                       }

    def _getinjpar(self,paramname):
        """
        Map parameter names to parameters in a SimInspiralTable .
        """
        if self._injection is not None:
            for key,value in self._injXMLFuncMap.items():
                if paramname in key:
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
            sample_array=(numpy.squeeze(pos_array[current_item,:]))
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
            median=str(oned_pos.median)
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
    def __init__(self,title=None,css=None):
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

def _histN(mat,N):
    """
    @deprecated: UNUSED .
    """
    Nd=size(N)
    histo=zeros(N)
    scale=array(map(lambda a,b:a/b,map(lambda a,b:(1*a)-b,map(max,mat),map(min,mat)),N))
    axes=array(map(lambda a,N:linspace(min(a),max(a),N),mat,N))
    bins=floor(map(lambda a,b:a/b , map(lambda a,b:a-b, mat, map(min,mat) ),scale*1.01))

    hbins=reshape(map(int,bins.flat),bins.shape)
    for co in transpose(hbins):
        t=tuple(co)
        histo[t[::-1]]=histo[t[::-1]]+1
    return (axes,histo)
#

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

    skypos=np.column_stack([posterior['ra'].samples,posterior['dec'].samples])

    injvalues=None

    sky_injpoint=(posterior['ra'].injval,posterior['dec'].injval)

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
    plx,ply=m(
              np.asarray(top_ranked_pixels)[::-1,1]*57.296,
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
    return '%ih%im%2.0fs'%(hours,mins,secs)

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
    return '%ideg%im%2.0fs'%(deg,sign*mins,sign*secs)

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

def plot_two_param_greedy_bins_contour(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,line_styles=__default_line_styles):
    """
    Plots the confidence level contours as determined by the 2-parameter
    greedy binning algorithm.

    @param posteriors_by_name: A dict containing Posterior instances referenced by some id.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}

    @param confidence_levels: a list of the required confidence levels to plot on the contour map.

    @param colors_by_name: A dict of colors cross-referenced to the above Posterior ids.

    """

    fig=plt.figure(1,figsize=(30,20),dpi=150)
    plt.clf()

    fig.add_axes([0.1,0.1,0.58,0.85])

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
    plt.xlabel(par1_name)
    plt.ylabel(par2_name)

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
    fig.savefig('test.png')
    twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')

    for text in twodcontour_legend.get_texts():
        text.set_fontsize('small')


    # For ra and dec set custom labels and for RA reverse
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
    plt.clf()

    #bins=(par1pos_Nbins, par2pos_Nbins)
    bins=(100,100)

    H, xedges, yedges = np.histogram2d(a,b, bins,normed=True)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    plt.imshow(H, aspect='equal', extent=None, interpolation='nearest')
    plt.colorbar()

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

    def _infmcmc_to_pos(self,files,deltaLogL=None,nDownsample=None,**kwargs):
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
            self._infmcmc_output_posterior_samples(files, outfile, logLThreshold, nskip)
        finally:
            outfile.close()
        return self._common_to_pos(open(postName,'r'))


    def _infmcmc_output_posterior_samples(self, files, outfile, logLThreshold, nskip=1):
        """
        Concatenate all the samples from the given files into outfile.
        For each file, only those samples past the point where the
        log(L) > logLThreshold are concatenated.
        """
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
                header=self._clear_infmcmc_header(infile)
                # Remove unwanted columns, and accound for 1 <--> 2 reversal of convention in lalinference.
                header=[self._swaplabel12(label) for label in header if label in allowedCols]
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
