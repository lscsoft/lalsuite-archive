

import os
import sys
import random
itertools = __import__("itertools")

import numpy as np
from scipy import stats

from pylal import plotutils
from pylal import rate

class PopStatement:

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
        
        # some combinations of colors/styles
        linestyles = []
        colors = []
        for linestyle in ['-','-.',':']:
            for color in ['b', 'g', 'r', 'c', 'm', 'y', 'k']:
                linestyles.append(linestyle)
                colors.append(color)
        self.linestyles = itertools.cycle(linestyles)
        self.colors = itertools.cycle(colors)
 

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

        self.on_lik_by_grb = np.asarray(self.on_lik_by_grb)

        # create the IFAR values for the onsource
        for on_lik, off_lik in zip(self.on_lik_by_grb,self.off_lik_by_grb):
            n_sample = float(len(off_lik))
            ifar = self.calculate_ifar(on_lik, off_lik, n_sample)
            self.on_ifar_by_grb.append(ifar)

        # create the combined/complete offsource sample
        self.off_lik = []
        self.off_ifar = []
        for off_lik, off_ifar in zip(self.off_lik_by_grb, self.off_ifar_by_grb):
            self.off_lik.extend(off_lik)
            self.off_ifar.extend(off_ifar)            


    def calculate_off_to_ifar(self, sample, sample_ref):
        """
        Calculates a list of FAR given the items in the sample
        """
        vector_ifar = []
        n_sample = float(len(sample))
        for item in sample:
            ifar = self.calculate_ifar(item, sample, n_sample)
            vector_ifar.append(ifar)
                        
                
        return vector_ifar

    def calculate_ifar(self, value, sample, n_sample):
        count_louder = (sample >= value).sum(axis=0)

        count_louder = max(count_louder, 1)
        return n_sample/count_louder
        #return count_louder/n_sample


    def remove_infs_from_list(self,list):
        """
        Just removes any inf's from a given list
        """        
        return [l for l in list if l != -np.inf and l != +np.inf]

    def get_min_max(self, list_of_lists_unequal, use_infs = False):
        """
        Get the min and max values from a list of lists,
        in which the sub-lists have different dimensions.
        Using Nick algorithm, email 3 Apr 2009.
        """

        for orig_list in list_of_lists_unequal:
            # NB: We want to ignore all +/-inf values, so
            # set them all to be +inf for min and -inf for max.
            tmp_arr = np.array(orig_list, dtype=float)
            tmp_arr = tmp_arr [~np.isinf(tmp_arr)]
            if len(tmp_arr) == 0:
                # NB: This also catches the case where len(orig_list) == 0.
                val_min = np.inf
                val_max = -np.inf
                continue
            val_min = tmp_arr.min()
            val_max = tmp_arr.max()

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
        
        for grb_name, off_pop in zip(self.list_grbs, off_by_grb):

            off_pop = self.remove_infs_from_list(off_pop)
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
        return n1 * n2 - u1  # return U for y


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
        z_val = (u_value - mean_U) / stdev_U
        
        return z_val

    def compute_wmu(self):
        """
        Computes the WMU z-score for the both cases
        using likelihood and the IFAR
        """
        
        z_lik = self.mannwhitney_u_zscore(self.use_lik, self.off_lik)
        z_ifar = self.mannwhitney_u_zscore(self.use_ifar, self.off_ifar)

        return z_lik, z_ifar

        

