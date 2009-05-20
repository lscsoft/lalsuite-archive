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
from __future__ import division

import os
import sys
import random
itertools = __import__("itertools")

import numpy as np
from scipy import stats

from glue import iterutils
from pylal import plotutils
from pylal import rate

####################################################
## class GRBdata
####################################################
class GRBdata(object):

    # -------------------------------------------------
    def __init__(self, grb_name, inj_by_trial, off_by_trial, on_source):
        """
        Initializes this object by passing the injection/offsource
        and onource data
        """

        # store the grb name
        self.grb_name = grb_name
        
        # store the likelihood data
        self.inj_by_trial = inj_by_trial
        self.off_by_trial = off_by_trial
        self.onsource = on_source
        
        # sort the given data
        self.inj_by_trial.sort()
        self.off_by_trial.sort()

    # -------------------------------------------------
    def choice(self, type):
      """
      Returns a particluar number from any sample
      """ 

      if type=='box':
        return self.onsource
      elif type=='random':
        index = random.randrange(len(self.off_by_trial))
        return self.off_by_trial[index]
      elif type=='max':
        return max(self.off_by_trial)
      elif type=='injectMax':
        return max(self.inj_by_trial)
      elif type=='injectUpper':
        n = len(self.inj_by_trial)
        index = random.randrange(int(0.7*n), n)
        return self.inj_by_trial[index]
      elif type=='injectMedian':
        n = len(self.inj_by_trial)
        index = random.randrange(int(0.3*n), int(0.5*n))
        return self.inj_by_trial[index]
      elif type=='injectLower':
        n = len(self.inj_by_trial)
        index = random.randrange(0, int(0.4*n))
        return self.inj_by_trial[index]
      elif type=='injectRange':
        pseudo_list = [x for x in self.inj_by_trial if x>3.5 and x<5.5]
        return random.choice(pseudo_list)
        
      elif type=='injectRandom':
        index = random.randrange(len(self.inj_by_trial))
        return self.inj_by_trial[index]
      else:
        raise ValueError, "Unknown type %s for trigger selection. " % type 
       
        
####################################################
#
####################################################
class PopStatement(object):

    # -------------------------------------------------
    def __init__(self, grb_data, name_suffix):
        """
        Initializes the class with the GRB data
        """
               
        # the data 
        self.off_list = []
        self.on_list = []

        # the main lists containing the data objects
        self.lik_by_grb = []
        self.ifar_by_grb = []

        # sets the statistic to be used: lik or ifar
        self.stat = None

        # sets the type of data to be used:
        # box or offsource/injection data
        self.type = None

        # a list of used GRBs (names only)
        self.list_grbs = []

        # a list of all GRBs
        self.grb_data = grb_data

        # two counters
        self.n_grb = 0
        self.n_off = 0

        # a name suffix for plot creation
        self.name_suffix = name_suffix

        # the statistical results
        self.p_one_sided = None
        self.p_two_sided = None
        self.u = None
        self.z = None
        
        self.colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

    # -------------------------------------------------
    def add_data(self,  grb_name, pop_injections_by_trial, \
                 pop_background_by_trial, pop_onsource):

        # append the GRB to the list
        self.list_grbs.append(grb_name)
     
        # create an instance for the likelihood data
        data_lik = GRBdata(grb_name, pop_injections_by_trial, \
                           pop_background_by_trial, pop_onsource)
        
        # calculate the IFAR as well
        off_ifar = self.calculate_sample_to_ifar(pop_background_by_trial,\
                                                 pop_background_by_trial)
        inj_ifar = [0] #FIXME: fix me
        on_ifar = self.calculate_sample_to_ifar([pop_onsource], \
                                                pop_background_by_trial)
        
        # and an instance for the IFAR data 
        data_ifar = GRBdata(grb_name, inj_ifar, off_ifar, on_ifar)
        
        # store the samples
        self.lik_by_grb.append(data_lik)
        self.ifar_by_grb.append(data_ifar)

    # -------------------------------------------------  
    def finalize(self):
        """
        Just sets some numbers
        """

        self.n_grb = len(self.lik_by_grb)

        self.n_off = 0
        for data_list in self.lik_by_grb:
            self.n_off += len(data_list.off_by_trial)
        
    # -------------------------------------------------  
    def calculate_sample_to_ifar(self, sample, sample_ref):
        """
        Calculates a list of FAR given the items in the sample
        """

        vector_ifar = []
        for item in sample:
            ifar = self.calculate_ifar(item, sample)
            vector_ifar.append(ifar)

        return vector_ifar
    
    # -------------------------------------------------=
    def calculate_ifar(self, value, sample):

        n_sample = float(len(sample))
        count_louder = (sample >= value).sum(axis=0)

        count_louder = max(count_louder, 1)
        return n_sample/count_louder  

    # -------------------------------------------------
    def check_off_distribution(self, pop_min, pop_max, far = False):

        data_list_by_grb = self.lik_by_grb
        if far:
           data_list_by_grb = self.ifar_by_grb
           
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
        nbins = 20
        bins = rate.LinearBins(pop_min, pop_max, nbins)
        px = bins.lower()
        
        for data_list in data_list_by_grb:

            grb_name = data_list.grb_name
            
            tmp_pop = data_list.off_by_trial
            tmp_arr = np.array(tmp_pop, dtype=float)
            off_pop = tmp_arr[~np.isinf(tmp_arr)]            
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
           

    # -------------------------------------------------
    def check_off_distribution_lik(self):
        return self.check_off_distribution(0.0, 25.0, far = False)

    # -------------------------------------------------
    def check_off_distribution_far(self):        
        return self.check_off_distribution(0.0, 3.0, far = True)

    # -------------------------------------------------
    def set_stat(self, stat):
        """
        Choose your statistics, either lik or ifar
        """
        self.stat = stat

    # -------------------------------------------------
    def select_onsource(self, type):
        """
        Selecting fake trials from the set of offsource
        trials for each GRB. 
        """

        # select the statistic value to be used 
        if self.stat=='ifar':
            data_list_by_grb = self.ifar_by_grb
        elif self.stat=='lik':
            data_list_by_grb = self.lik_by_grb               
        else:
            raise ValueError, "Unknown statistics %s in select_onsource" %\
                  self.stat

        # delete the old items
        self.on_list = []
        self.off_list = []
        self.type = type

        for counter, data_list in \
                enumerate(data_list_by_grb):
            
            if type=='box':
                value = data_list.choice(type) 

            elif 'inject' in type:
                number_inj = int(type[-1])
                if counter<number_inj:
                    if 'Max' in type:  
                        value = data_list.choice('injectMax')
                        print "injectMax lik: ", value
                    elif 'Random' in type:
                        value = data_list.choice('injectRandom')
                        print "injectRandom lik: ", value
                    elif 'Upper' in type:
                        value = data_list.choice('injectUpper')
                        print "injectUpper lik: ", value
                    elif 'Median' in type:
                        value = data_list.choice('injectMedian')
                        print "injectMedian lik: ", value
                    elif 'Lower' in type:
                        value = data_list.choice('injectLower')
                        print "injectLower lik: ", value
                    elif 'Range' in type:
                        value = data_list.choice('injectRange')
                        print "injectRange lik: ", value
                    else:
                        raise ValueError,"Unknown injection type: ", type
                else:
                    value = data_list.choice('box')
                    print "box lik: ", value
                 
            elif type=='random' or (type=='single' and counter>0):
                value = data_list.choice('random')
                
            elif type=='max' or (type=='single' and counter==0):
                value = data_list.choice('max')
                
            else:
                raise ValueError,"Unknown selection type: ", type

            # fill the data to be used for the test
            self.on_list.append(value)
            self.off_list.extend(data_list.off_by_trial) 

    # -------------------------------------------------
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

    # -------------------------------------------------
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

    # -------------------------------------------------
    def compute_wmu(self):
        """
        Computes the WMU z-score for the both cases
        using likelihood and the IFAR
        """
        # compute the z-value
        z_value = self.mannwhitney_u_zscore(self.on_list, self.off_list)

        # sf = 1 - cdf
        self.p_one_sided = stats.distributions.norm.sf(z_value)
         
        # erfc = 1 - erf
        self.p_two_sided = stats.erfc(abs(z_value) / np.sqrt(2.))  
         
        return z_value

    # -------------------------------------------------
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

    # -------------------------------------------------    
    def create_plot_hist(self):
                
        #
        # Create the histogram comparison
        #
        plot_title = r"$m_2 \in [%s), U=%d, z_U=%s, p_1=%s, p_2=%s$" \
            % (self.name_suffix , int(self.u), self.float_to_latex(self.z,\
                                                                   "%5.2g"),
               self.float_to_latex(self.p_one_sided, "%5.2g"),
               self.float_to_latex(self.p_two_sided, "%5.2g"))
        plot = plotutils.VerticalBarHistogram(r"$IFAR(m_2 \in [%s))$" %\
                                              self.name_suffix, "PDF", \
                                              plot_title)
        plot.add_content(self.on_list, color='r', label = r'On source', \
                         bottom = 1.0e-4)
        plot.add_content(self.off_list, color='b', label = r'Off source', \
                         bottom = 1.0e-4)
        plot.finalize(normed=True)
        plot.ax.set_yscale('log')
        return plot

    # ------------------------------------------------- 
    def create_plot_qq(self):
        
        #
        # Create the QQ plot
        #
        plot_title = r"$m_2 \in [%s), U=%d, z_U=%s, p_1=%s, p_2=%s$" \
                     % (self.name_suffix , int(self.u), \
                        self.float_to_latex(self.z, "%5.2g"),
                        self.float_to_latex(self.p_one_sided, "%5.2g"),
                        self.float_to_latex(self.p_two_sided, "%5.2g"))
        if self.type=="box":
            plot_title = r"$ $"
        plot = plotutils.QQPlot(r"self quantile", "combined quantile", \
                                plot_title)
        plot.add_bg(self.off_list, linewidth = 3, label="\"Off source\"")
        plot.add_fg(self.on_list, color='r', marker = 'o',\
                    label = r'On source',\
                    linestyle='None',markersize=10)    
        plot.finalize()
        return plot
    
    # -------------------------------------------------
    def create_cdf_plot(self):

        # create the cumulative data set
        self.off_list.sort()
        on_data = np.asarray(self.on_list)
        y_on = []
        y_off = []
        x_onoff = []
        counter = 0
        for value in self.off_list[::-1]:
            counter += 1
            x_onoff.append(value)
            y_off.append(counter)
            y_on.append((on_data>value).sum())

        # replace the infinite values by 0.0 for plotting purposes
        x_onoff = np.asarray(x_onoff)
        y_on = np.asarray(y_on)
        y_off = np.asarray(y_off)                    
        inf_ind = np.isinf(x_onoff)
        x_onoff[inf_ind] = 0.0

        # errors on the background counting?
        #y_err = np.asarray([np.sqrt(n) for n in y_off])
        
        # scale the values correspondingly
        scale = y_on.max()/y_off.max()
        y_off_scale = y_off*scale
        #y_err_scale = y_err*scale

        # create the plot
        plot = plotutils.SimplePlot(r"Likelihood", r"Cumulative sum",\
                                    r"Cumulative distribution")
        plot.add_content(x_onoff, y_off_scale, color = 'r', \
                         linewidth = 2, label = 'off-source')
        plot.add_content(x_onoff, y_on, color = 'b',\
                         linewidth = 3, label = 'on-source')
        ## for x, y, err in zip(x_onoff, y_off_scale, y_err_scale):
##             plot.ax.plot([x,x],[y-err, y+err], 'g')
            
        plot.finalize()
        
        # make the plot nice
        plot.ax.set_xscale('log')
        plot.ax.axis([2, 20, 0.0, 22.0])
        plot.ax.set_xticks([2,3,5,10])
        plot.ax.set_xticklabels(['2','3','5','10'])

        return plot

    # -------------------------------------------------
    def create_pdf_plot(self):
        """
        Create a pdf plot of the used data.
        """
    
        # get the data and sort them
        data_off = self.off_list
        data_on = self.on_list
        data_off.sort(reverse=True)
        data_on.sort(reverse=True)

        # transform to array and remove infinities
        data_off = np.asarray(data_off)
        data_on = np.asarray(data_on)                
        inf_ind = np.isinf(data_on)
        data_on[inf_ind] = 0.0
        inf_ind = np.isinf(data_off)
        data_off[inf_ind] = 0.0

        # prepare lists
        pdf_x = []
        pdf_off_y = []
        pdf_on_y = []
        pdf_on_r = []

        # prepare the main loop
        index = 0
        sum_off = 0.0
        sum_tot = 0
        for d in data_off:

            # add one counter to the offsource sum
            sum_off += 1
            sum_tot +=1

            # have we reached an onsource statistics value?
            if d<=data_on[index]:

                # compute the values to be stored
                sum_scale = sum_off/self.n_off
                y = 1.0/self.n_grb
                r = y-sum_scale
                
                # store these values in the lists
                pdf_x.append(d)
                pdf_off_y.append(sum_scale)
                pdf_on_y.append(y)        
                pdf_on_r.append(r)
                
                # check breaking condition
                if index==self.n_grb-1:
                    break

                # reset the sum and increase the index                
                index += 1
                sum_off = 0.0
        
        # double check the data
        remain_data = data_off[data_off<=d]
        sum_off = sum(pdf_off_y)+len(remain_data)/self.n_off
        sum_on = sum(pdf_on_y)
        print "Sum of the onsource: %.1f  offsource: %.1f" % \
              (sum_on, sum_off)
                
        # create the plot
        plot = plotutils.SimplePlot(r"Likelihood", r"pdf",\
                                    r"PDF Distribution")
        plot.add_content(pdf_x, pdf_off_y, color = 'r', \
                         linewidth = 2, label = 'off-source')
        plot.add_content(pdf_x, pdf_on_y, color = 'b', marker = 'o',\
                         markersize = 10.0, label = 'on-source')
        plot.finalize()
        #plot.ax.axis([2, 6.0, 0.0, 0.25])

        return plot
    
    # -------------------------------------------------
    def create_hist_plot(self, n_bin, range = None):

        def create_area(x, dx, y, dy):
            px = [x-dx/2, x+dx/2, x+dx/2, x-dx/2, x-dx/2]
            py = [y-dy/2, y-dy/2, y+dy/2, y+dy/2, y-dy/2]
            return px, py
        
        def draw_error_boxes(plot, x, dx, y, dy, col):
            
            bx, by = create_area(x, dx, y, dy )
            plot.ax.fill(bx, by, ec='w', fc = col, alpha = 0.2)
            
            bx, by = create_area(x, dx, y, 2*dy )
            plot.ax.fill(bx, by, ec='w', fc = col, alpha = 0.2)
            
            bx, by = create_area(x, dx, y, 3*dy )
            plot.ax.fill(bx, by, ec='k', fc = col, alpha = 0.2)

            return plot
        

        # set the surroundings of the parameter space
        if range is None:
            data_on = np.asarray(self.on_list)                
            inf_ind = np.isinf(data_on)
            val_min = data_on[~inf_ind].min()
            val_max = data_on.max()
        else:
            val_min = range[0]
            val_max = range[1]

        # create the hists
        hist_on = np.zeros(n_bin)
        hist_off = np.zeros(n_bin)   
        
        # create the rate bins
        lik_bins = rate.LinearBins(val_min, val_max, n_bin) 

        # and fill the histograms
        for x in self.off_list:
            if x>=val_min and x<=val_max:            
                hist_off[lik_bins[x]] += 1
            
        for x in self.on_list:
            if x>=val_min and x<=val_max:            
                hist_on[lik_bins[x]] += 1

        # get the centres
        px = lik_bins.centres()
        dx = px[1]-px[0]

        # norm the histograms
        norm = self.n_grb/self.n_off
        hist_on_norm = hist_on
        hist_off_norm = norm*hist_off

        # create the plot
        plot = plotutils.SimplePlot(r"Likelihood", r"counts",\
                                    r"Histogramm of on/offsource with %d bins"%\
                                    (n_bin))

        # the norm of the normed histograms: 1.0
        plot.add_content(px, hist_on, color = 'r', marker = 'o',\
                 markersize = 10.0, label = 'on-source')
        plot.add_content(px, hist_off_norm, color = 'b',marker = 's', \
                 markersize = 10, label = 'off-source')

        # add the error shading
        for x,n in zip(px, hist_on):
            
            # calculate the error (just the statistic error)
            dn = np.sqrt(n)
            plot = draw_error_boxes(plot, x, dx, n, dn, 'r')

        plot.finalize()           
        plot.ax.axis([val_min, val_max, 0.0, 20.0])
        #return plot
            
        # insert the lines of the data itself
        for x in self.on_list:
            if x>=val_min and x<=val_max: 
                plot.ax.plot([x,x],[0.0, 1.0],'k')            
        plot.ax.axis([val_min, val_max, 0.0, 20.0])
        
        return plot


