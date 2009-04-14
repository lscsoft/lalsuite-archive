#!/usr/bin/env python
#
# Copyright (C) 2008  Nickolas Fotopoulos and Alexander Dietz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
import os
import sys
import random
itertools = __import__("itertools")

import numpy as np
from scipy import stats

from glue import iterutils
from pylal import plotutils
from pylal import rate

class PopStatement(object):

    def __init__(self, grb_data, name_suffix):
        """
        Initializes the class with the background population
        """

        self.off_lik_by_grb = []
        self.off_ifar_by_grb = []
        self.on_lik_by_grb = []
        self.on_ifar_by_grb = []

        self.off_lik = []
        self.off_ifar = []        
        self.use_lik = []
        self.use_ifar = []

        
        self.list_grbs = []
        self.grb_data = grb_data
        self.name_suffix = name_suffix

        
        self.p_one_sided = None
        self.p_two_sided = None
        self.u = None
        self.z = None
        
        self.colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
 

    def add_background_by_trial(self, grb_name, pop_background_by_trial):
        """
        Adding background values from a different GRB to
        this class.
        """

        self.off_lik_by_grb.append(list(pop_background_by_trial))
        vector_ifar = self.calculate_off_to_ifar(pop_background_by_trial,\
                                                 pop_background_by_trial)
        self.off_ifar_by_grb.append(vector_ifar)

        self.list_grbs.append(grb_name)
        
    def add_foreground(self, grb_name, value_foreground):
        """
        Adding foreground value from a different GRB to
        this class.
        """
        self.on_lik_by_grb.append(value_foreground)


    def finalize(self):
        """
        Finalize the values, recalculates the FAR of the
        foreground values
        """

        #self.on_lik_by_grb = np.asarray(self.on_lik_by_grb)

        # create the IFAR values for the onsource
       ##  for on_lik, off_lik in zip(self.on_lik_by_grb,self.off_lik_by_grb):
##             n_sample = float(len(off_lik))
##             ifar = self.calculate_ifar(on_lik, off_lik, n_sample)
##             self.on_ifar_by_grb.append(ifar)

        self.on_ifar_by_grb = map(self.calculate_ifar, \
                                      self.on_lik_by_grb, self.off_lik_by_grb)

        # create the combined/complete offsource sample
        self.off_lik = list(iterutils.flatten(self.off_lik_by_grb))
        self.off_ifar = list(iterutils.flatten(self.off_ifar_by_grb))
        
    def calculate_off_to_ifar(self, sample, sample_ref):
        """
        Calculates a list of FAR given the items in the sample
        """
        vector_ifar = []
        for item in sample:
            ifar = self.calculate_ifar(item, sample, len(sample))
            vector_ifar.append(ifar)

        return vector_ifar


    def calculate_ifar(self, value, sample):

        n_sample = float(len(sample))
        count_louder = (sample >= value).sum(axis=0)

        count_louder = max(count_louder, 1)
        return n_sample/count_louder


   ##  def remove_infs_from_list(self,list):
##         """
##         Removing infinite values
##         """        
##         return [l for l in list if l != -np.inf and l != +np.inf]

    def get_min_max(self, list_of_lists_unequal, use_infs = False):
        """
        Get the min and max values from a list of lists,
        in which the sub-lists have different dimensions.
        Using Nick algorithm, email 3 Apr 2009.
        """

        val_min = +np.infty
        val_min = -np.infty        
        for ilist in list_of_lists_unequal:                        
            val_min = min(ilist + [val_min])
            val_max = max(ilist + [val_min])            
            
            ## val_max = max(ilist + [val_min])
##             # NB: We want to ignore all +/-inf values, so
##             # set them all to be +inf for min and -inf for max.
##             tmp_arr = np.array(orig_list, dtype=float)
##             tmp_arr = tmp_arr [~np.isinf(tmp_arr)]
##             if len(tmp_arr) == 0:
##                 # NB: This also catches the case where len(orig_list) == 0.
##                 val_min = np.inf
##                 val_max = -np.inf
##                 continue
##             val_min = tmp_arr.min()
##             val_max = tmp_arr.max()

        return val_min, val_max
    
    def check_off_distribution(self, off_by_grb, far = False):

        
        # prepare the plot
        if far:
            tag = 'log(IFAR)'
            sname = 'far'
        else:
            tag = 'Likelihood'
            sname = 'lik'
        plot = plotutils.SimplePlot(tag, r"cumulative sum",\
                                    r"Cumulative distribution offsource")
        
        # create the hist data in a consistent binning
        pop_min, pop_max = self.get_min_max(off_by_grb)
        nbins = 20
        bins = rate.LinearBins(pop_min, pop_max, nbins)
        px = bins.lower()
        
        for grb_name, tmp_pop in zip(self.list_grbs, off_by_grb):

            tmp_arr = np.array(tmp_pop, dtype=float)
            off_pop = tmp_arr[~np.isinf(tmp_arr)]            
            #off_pop = self.remove_infs_from_list(off_pop)
            off_pop.sort()
            py = range(len(off_pop), 0, -1)
            if far:
                off_pop = np.log10(off_pop)


            ifos = self.grb_data[grb_name]['ifos']
            if ifos=='H1L1':
                linestyle = '-'
            elif ifos == 'H1H2':
                linestyle = '-.'
            else:
                linestyle = ':'
            
            
            # add content to the plot
            plot.add_content(off_pop, py, color = self.colors.next(),\
                             linestyle = linestyle, label=grb_name)
        plot.finalize()
        plot.ax.set_yscale("log")
        return plot
        
   

    def check_off_distribution_lik(self):
        return self.check_off_distribution(self.off_lik_by_grb, far = False)

    def check_off_distribution_far(self):        
        return self.check_off_distribution(self.off_ifar_by_grb, far = True)


    def select_onsource(self, type):
        """
        Selecting fake trials from the set of offsource
        trials for each GRB. 
        """

        # get the number of onsource values
        n = len(self.on_lik_by_grb)

        # delete the old items
        self.use_lik = []
        self.use_ifar = []

        if type=='box':
            self.use_lik  = self.on_lik_by_grb
            self.use_ifar = self.on_ifar_by_grb
            return

        
        for counter, (off_lik, off_ifar) in \
                enumerate(zip(self.off_lik_by_grb, self.off_ifar_by_grb)):

            if len(off_lik) != len(off_ifar):
                print "ERROR: Lengths different..."
                
            if type=='random' or (type=='single' and counter>0):
                index = random.randrange(len(off_lik))
                self.use_lik.append(off_lik[index])
                self.use_ifar.append(off_ifar[index])
            elif type=='max' or (type=='single' and counter==0):
                self.use_lik.append(max(off_lik))
                self.use_ifar.append(max(off_ifar))
                
    
    def mannwhitney_u(self, x, y):
        """
        Return the Mann-Whitney U statistic on the provided scores.  Copied from
        scipy.stats.mannwhitneyu except that we only return the U such that
        large U means that population x was systematically larger than population
        y, rather than the smaller U between x and y.  The two possible U values
        one can report are related by U' = n1*n2 - U.
        """
        x = np.asarray(x)
        y = np.asarray(y)
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError, "populations must be rank 1 collections"
        n1 = len(x)
        n2 = len(y)

        ranked = stats.rankdata(np.concatenate((x,y)))
        rankx = ranked[0:n1]  # get the x-ranks
        u1 = n1 * n2 + (n1 * (n1 + 1)) / 2.0 - rankx.sum()  # calc U for x
        self.u =  n1 * n2 - u1  # return U for y
        return self.u

    def mannwhitney_u_zscore(self, pop_test, pop_ref):
        """
        Return the z-score of a given sample.
        Not appropriate for n1 + n2 < ~20.
        """

        n1 = len(pop_test)
        n2 = len(pop_ref)            
        u_value = self.mannwhitney_u(pop_test, pop_ref)
        mean_U = n1 * n2 / 2
        stdev_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        self.z = (u_value - mean_U) / stdev_U
        
        return self.z

    def compute_wmu(self):
        """
        Computes the WMU z-score for the both cases
        using likelihood and the IFAR
        """

        z_ifar = self.mannwhitney_u_zscore(self.use_ifar, self.off_ifar)        
        z_lik = self.mannwhitney_u_zscore(self.use_lik, self.off_lik)

        # sf = 1 - cdf
        self.p_one_sided = stats.distributions.norm.sf(z_lik)
         
        # erfc = 1 - erf
        self.p_two_sided = stats.erfc(abs(z_lik) / np.sqrt(2.))  
         
        return z_lik, z_ifar

    def float_to_latex(self, x, format="%g"):
        """
        Convert a floating point number to a latex representation.  In particular,
        scientific notation is handled gracefully: e -> 10^
        """
        base_str = format % x
        if "e" not in base_str:
            return base_str
        mantissa, exponent = base_str.split("e")
        exponent = str(int(exponent))  # remove leading 0 or +
    
    def create_plot_hist(self):
                
        #
        # Create the histogram comparison
        #
        plot_title = r"$m_2 \in [%s), U=%d, z_U=%s, p_1=%s, p_2=%s$" \
            % (self.name_suffix , int(self.u), self.float_to_latex(self.z, "%5.2g"),
               self.float_to_latex(self.p_one_sided, "%5.2g"),
               self.float_to_latex(self.p_two_sided, "%5.2g"))    
        plot = plotutils.VerticalBarHistogram(r"$IFAR(m_2 \in [%s))$" %\
                                              self.name_suffix, "PDF", plot_title)
        plot.add_content(self.use_lik, color='r', label = r'On source', bottom = 1.0e-4)
        plot.add_content(self.off_lik, color='b', label = r'Off source', bottom = 1.0e-4)
        plot.finalize(normed=True)
        plot.ax.set_yscale('log')
        return plot


    def create_plot_qq(self):
        
        #
        # Create the QQ plot
        #
        plot_title = r"$m_2 \in [%s), U=%d, z_U=%s, p_1=%s, p_2=%s$" \
                     % (self.name_suffix , int(self.u), self.float_to_latex(self.z, "%5.2g"),
                        self.float_to_latex(self.p_one_sided, "%5.2g"),
                        self.float_to_latex(self.p_two_sided, "%5.2g")) 
        plot = plotutils.QQPlot(r"self quantile", "combined quantile", plot_title)
        plot.add_bg(self.off_lik, linewidth = 3, label="\"Off source\"")
        plot.add_fg(self.use_lik, color='r', marker = 'o',label = r'On source',\
                         linestyle='None',markersize=10)    
        plot.finalize()
        return plot

