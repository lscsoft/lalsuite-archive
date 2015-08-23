#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPostProc.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#       Vivien Raymond <vivien.raymond@ligo.org>
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

#standard library imports
import sys
import os

from math import ceil,floor
import cPickle as pickle

from time import strftime

#related third party imports
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack,cov,unique,hsplit,correlate,log,dot,power,squeeze,sort
from scipy import stats

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

# Default font properties
fig_width_pt = 246  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 16,
        'font.size':       16,
        'legend.fontsize': 16,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'text.usetex': False,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ['Times','Palatino','New Century Schoolbook','Bookman','Computer Modern Roman','Times New Roman','Liberation Serif'],
        'font.weight':'normal',
        'font.size':16,
        'savefig.dpi': 120
        })

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

from glue.ligolw import table
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def email_notify(address,path):
    import smtplib
    import subprocess
    address=address.split(',')
    SERVER="localhost"
    USER=os.environ['USER']
    HOST=os.environ['HOSTNAME']
    FROM=USER+'@'+HOST
    SUBJECT="LALInference result is ready at "+HOST+"!"
    # Guess the web space path for the clusters
    fslocation=os.path.abspath(path)
    webpath='posplots.html'
    if 'public_html' in fslocation:
        k='public_html/'
    elif 'WWW' in fslocation:
        k='WWW/'
    elif 'www_html' in fslocation:
        k='www_html/'
    else:
        k=None
    if k is not None:
        (a,b)=fslocation.split(k)
        webpath=os.path.join('~%s'%(USER),b,webpath)
        onweb=True
    else:
        webpath=os.path.join(fslocation,'posplots.html')
        onweb=False

    if 'atlas.aei.uni-hannover.de' in HOST:
        url="https://atlas1.atlas.aei.uni-hannover.de/"
    elif 'ligo.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo.caltech.edu/"
    elif 'ligo-wa.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo-wa.caltech.edu/"
    elif 'ligo-la.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo-la.caltech.edu/"
    elif 'phys.uwm.edu' in HOST:
        url="https://ldas-jobs.phys.uwm.edu/"
    elif 'phy.syr.edu' in HOST:
        url="https://sugar-jobs.phy.syr.edu/"
    elif 'arcca.cf.ac.uk' in HOST:
        url="https://geo2.arcca.cf.ac.uk/"
    else:
        if onweb:
          url="http://%s/"%(HOST)
        else:
          url=HOST+':'
    url=url+webpath

    TEXT="Hi "+USER+",\nYou have a new parameter estimation result on "+HOST+".\nYou can view the result at "+url+"\n"
    message="From: %s\nTo: %s\nSubject: %s\n\n%s"%(FROM,', '.join(address),SUBJECT,TEXT)
    server=smtplib.SMTP(SERVER)
    server.sendmail(FROM,address,message)
    server.quit()

#Import content handler 
from pylal.SimInspiralUtils import ExtractSimInspiralTableLIGOLWContentHandler
lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)


def pickle_to_file(obj,fname):
    """
    Pickle/serialize 'obj' into 'fname'.
    """
    filed=open(fname,'w')
    pickle.dump(obj,filed)
    filed.close()

def oneD_dict_to_file(dict,fname):
    filed=open(fname,'w')
    for key,value in dict.items():
        filed.write("%s %s\n"%(str(key),str(value)) )

def multipleFileCB(opt, opt_str, value, parser):
    args=[]

    def floatable(str):
      try:
        float(str)
        return True
      except ValueError:
        return False

    for arg in parser.rargs:
      # stop on --foo like options
      if arg[:2] == "--" and len(arg) > 2:
        break
      # stop on -a, but not on -3 or -3.0
      if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
        break
      args.append(arg)

    del parser.rargs[:len(args)]
    #Append new files to list if some already specified
    if getattr(parser.values, opt.dest):
        oldargs = getattr(parser.values, opt.dest)
        oldargs.extend(args)
        args = oldargs
    setattr(parser.values, opt.dest, args)

def cbcBayesPostProc(
                        outdir,data,oneDMenu,twoDGreedyMenu,GreedyRes,
                        confidence_levels,twoDplots,
                        #misc. optional
                        injfile=None,eventnum=None,
                        trigfile=None,trignum=None,
                        skyres=None,
                        #direct integration evidence
                        dievidence=False,boxing=64,difactor=1.0,
                        #elliptical evidence
                        ellevidence=False,
                        #manual input of bayes factors optional.
                        bayesfactornoise=None,bayesfactorcoherent=None,
                        #manual input for SNR in the IFOs, optional.
                        snrfactor=None,
                        #nested sampling options
                        ns_flag=False,ns_Nlive=None,
                        #spinspiral/mcmc options
                        ss_flag=False,ss_spin_flag=False,
                        #lalinferenceMCMC options
                        li_flag=False,deltaLogL=None,fixedBurnins=None,nDownsample=None,oldMassConvention=False,
                        #followupMCMC options
                        fm_flag=False,
                        #spin frame for the injection
                        inj_spin_frame='OrbitalL',
                        # on ACF?
                        noacf=False,
                        #Turn on 2D kdes
                        twodkdeplots=False,
                        #Turn on R convergence tests
                        RconvergenceTests=False,
                        # Save PDF figures?
                        savepdfs=True,
                        #List of covariance matrix csv files used as analytic likelihood
                        covarianceMatrices=None,
                        #List of meanVector csv files used, one csv file for each covariance matrix
                        meanVectors=None,
                        #header file
                        header=None,
                        psd_files=None,
                        greedy=True ## If true will use greedy bin for 1-d credible regions. Otherwise use 2-steps KDE
                    ):
    """
    This is a demonstration script for using the functionality/data structures
    contained in pylal.bayespputils . It will produce a webpage from a file containing
    posterior samples generated by the parameter estimation codes with 1D/2D plots
    and stats from the marginal posteriors for each parameter/set of parameters.
    """
    votfile=None
    if eventnum is not None and injfile is None:
        print "You specified an event number but no injection file. Ignoring!"

    if trignum is not None and trigfile is None:
        print "You specified a trigger number but no trigger file. Ignoring!"

    if trignum is None and trigfile is not None:
        print "You specified a trigger file but no trigger number. Taking first entry (the case for GraceDB events)."
        trignum=0

    if data is None:
        raise RuntimeError('You must specify an input data file')
    #
    if outdir is None:
        raise RuntimeError("You must specify an output directory.")

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #
    if fm_flag:
        peparser=bppu.PEOutputParser('fm')
        commonResultsObj=peparser.parse(data)

    elif ns_flag and not ss_flag:
        peparser=bppu.PEOutputParser('ns')
        commonResultsObj=peparser.parse(data,Nlive=ns_Nlive)

    elif ss_flag and not ns_flag:
        peparser=bppu.PEOutputParser('mcmc_burnin')
        commonResultsObj=peparser.parse(data,spin=ss_spin_flag,deltaLogL=deltaLogL)

    elif li_flag:
        peparser=bppu.PEOutputParser('inf_mcmc')
        commonResultsObj=peparser.parse(data,outdir=outdir,deltaLogL=deltaLogL,fixedBurnins=fixedBurnins,nDownsample=nDownsample,oldMassConvention=oldMassConvention)

    elif ss_flag and ns_flag:
        raise RuntimeError("Undefined input format. Choose only one of:")

    elif '.xml' in data[0]:
        peparser=bppu.PEOutputParser('xml')
        commonResultsObj=peparser.parse(data[0])
        thefile=open(data[0],'r')
        votfile=thefile.read()
    else:
        peparser=bppu.PEOutputParser('common')
        commonResultsObj=peparser.parse(open(data[0],'r'),info=[header,None])
    
    #Extract f_ref from CRO if present.  This is needed to calculate orbital angular momentum
    #  when converting spin parameters.  Ideally this info will be provided in the
    #  SimInspiralTable in the near future.
    ps,samps = commonResultsObj
    try:
        f_refIdx = ps.index('f_ref')
        injFref = unique(samps[:,f_refIdx])
        if len(injFref) > 1:
            print "ERROR: Expected f_ref to be constant for all samples.  Can't tell which value was injected!"
            print injFref
            injFref = None
        else:
            injFref = injFref[0]
    except ValueError:
        injFref = None

    #Select injections using tc +/- 0.1s if it exists or eventnum from the injection file
    injection=None
    if injfile and eventnum is not None:
        print 'Looking for event %i in %s\n'%(eventnum,injfile)
        xmldoc = utils.load_filename(injfile,contenthandler=ExtractSimInspiralTableLIGOLWContentHandler)
        siminspiraltable=lsctables.table.get_table(xmldoc,lsctables.SimInspiralTable.tableName)
        injection=siminspiraltable[eventnum]
	#injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])
	#if(len(injections)!=1): raise RuntimeError('Error: something unexpected happened while loading the injection file!\n')
        #injection=injections[0]

    #Get trigger
    triggers = None
    if trigfile is not None and trignum is not None:
        triggers = bppu.readCoincXML(trigfile, trignum)

    ## Load Bayes factors ##
    # Add Bayes factor information to summary file #
    if bayesfactornoise is not None:
        bfile=open(bayesfactornoise,'r')
        BSN=bfile.read()
        bfile.close()
        if(len(BSN.split())!=1):
          BSN=BSN.split()[0]
        print 'BSN: %s'%BSN
    if bayesfactorcoherent is not None:
        bfile=open(bayesfactorcoherent,'r')
        BCI=bfile.read()
        bfile.close()
        print 'BCI: %s'%BCI

    if snrfactor is not None:
        if not os.path.isfile(snrfactor):
            print "No snr file provided or wrong path to snr file\n"
            snrfactor=None
        else:
            snrstring=""
            snrfile=open(snrfactor,'r')
            snrs=snrfile.readlines()
            snrfile.close()
            for snr in snrs:
                if snr=="\n":
                    continue
                snrstring=snrstring +" "+str(snr[0:-1])+" ,"
            snrstring=snrstring[0:-1]

    #Create an instance of the posterior class using the posterior values loaded
    #from the file and any injection information (if given).
    pos = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injection,inj_spin_frame=inj_spin_frame, injFref=injFref,SnglInpiralList=triggers,votfile=votfile)
  
    #Create analytic likelihood functions if covariance matrices and mean vectors were given
    analyticLikelihood = None
    if covarianceMatrices and meanVectors:
        analyticLikelihood = bppu.AnalyticLikelihood(covarianceMatrices, meanVectors)

        #Plot only analytic parameters
        oneDMenu = analyticLikelihood.names
        if twoDGreedyMenu:
            twoDGreedyMenu = []
            for i in range(len(oneDMenu)):
                for j in range(i+1,len(oneDMenu)):
                    twoDGreedyMenu.append([oneDMenu[i],oneDMenu[j]])
        twoDplots = twoDGreedyMenu

    if eventnum is None and injfile is not None:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])

        if(len(injections)<1):
            try:
                print 'Warning: Cannot find injection with end time %f' %(pos['time'].mean)
            except KeyError:
                print "Warning: No 'time' column!"

        else:
            try:
                injection = itertools.ifilter(lambda a: abs(float(a.get_end()) - pos['time'].mean) < 0.1, injections).next()
                pos.set_injection(injection)
            except KeyError:
                print "Warning: No 'time' column!"

    try:
      pos.extend_posterior()
    except:
      pass

    #Perform necessary mappings
    functions = {'cos':cos,'sin':sin,'exp':exp,'log':log}
    for pos_name in oneDMenu:
        if pos_name not in pos.names:
            for func in functions.keys():
                old_pos_name = pos_name.replace(func,'')
                if pos_name.find(func)==0 and old_pos_name in pos.names:
                    print "Taking %s of %s ..."% (func,old_pos_name)
                    pos.append_mapping(pos_name,functions[func],old_pos_name)

    #Remove samples with NaNs in requested params
    requested_params = set(pos.names).intersection(set(oneDMenu))
    pos.delete_NaN_entries(requested_params)

    #Remove non-analytic parameters if analytic likelihood is given:
    if analyticLikelihood:
        dievidence_names = ['post','posterior','logl','prior','likelihood','cycle','chain']
        [pos.pop(param) for param in pos.names if param not in analyticLikelihood.names and param not in dievidence_names]

    ##Print some summary stats for the user...##
    #Number of samples
    print "Number of posterior samples: %i"%len(pos)
    # Means
    print 'Means:'
    print str(pos.means)
    #Median
    print 'Median:'
    print str(pos.medians)
    #maxL
    print 'maxL:'
    max_pos,max_pos_co=pos.maxL
    print max_pos_co

    #==================================================================#
    #Create web page
    #==================================================================#

    html=bppu.htmlPage('Posterior PDFs',css=bppu.__default_css_string,javascript=bppu.__default_javascript_string)

    #Create a section for meta-data/run information
    html_meta=html.add_section('Summary')
    table=html_meta.tab()
    row=html_meta.insert_row(table,label='thisrow')
    td=html_meta.insert_td(row,'',label='Samples')
    SampsStats=html.add_section_to_element('Samples',td)
    SampsStats.p('Produced from '+str(len(pos))+' posterior samples.')
    if 'chain' in pos.names:
        acceptedChains = unique(pos['chain'].samples)
        acceptedChainText = '%i of %i chains accepted: %i'%(len(acceptedChains),len(data),acceptedChains[0])
        if len(acceptedChains) > 1:
            for chain in acceptedChains[1:]:
                acceptedChainText += ', %i'%(chain)
        SampsStats.p(acceptedChainText)
    if 'cycle' in pos.names:
        SampsStats.p('Longest chain has '+str(pos.longest_chain_cycles())+' cycles.')
    filenames='Samples read from %s'%(data[0])
    if len(data) > 1:
        for fname in data[1:]:
            filenames+=', '+str(fname)
    SampsStats.p(filenames)
    td=html_meta.insert_td(row,'',label='SummaryLinks')
    legend=html.add_section_to_element('Sections',td)
    
    #Create a section for model selection results (if they exist)
    if bayesfactornoise is not None:
        html_model=html.add_section('Model selection',legend=legend)
        html_model.p('log Bayes factor ( coherent vs gaussian noise) = %s, Bayes factor=%f'%(BSN,exp(float(BSN))))
        if bayesfactorcoherent is not None:
            html_model.p('log Bayes factor ( coherent vs incoherent OR noise ) = %s, Bayes factor=%f'%(BCI,exp(float(BCI))))

    if dievidence:
        html_model=html.add_section('Direct Integration Evidence',legend=legend)
        log_ev = log(difactor) + pos.di_evidence(boxing=boxing)
        ev=exp(log_ev)
        evfilename=os.path.join(outdir,"evidence.dat")
        evout=open(evfilename,"w")
        evout.write(str(ev))
        evout.write(" ")
        evout.write(str(log_ev))
        evout.close()
        print "Computing direct integration evidence = %g (log(Evidence) = %g)"%(ev, log_ev)
        html_model.p('Direct integration evidence is %g, or log(Evidence) = %g.  (Boxing parameter = %d.)'%(ev,log_ev,boxing))
        if 'logl' in pos.names:
            log_ev=pos.harmonic_mean_evidence()
            html_model.p('Compare to harmonic mean evidence of %g (log(Evidence) = %g).'%(exp(log_ev),log_ev))

    if ellevidence:
        try:
            html_model=html.add_section('Elliptical Evidence',legend=legend)
            log_ev = pos.elliptical_subregion_evidence()
            ev = exp(log_ev)
            evfilename=os.path.join(outdir, 'ellevidence.dat')
            evout=open(evfilename, 'w')
            evout.write(str(ev) + ' ' + str(log_ev))
            evout.close()
            print 'Computing elliptical region evidence = %g (log(ev) = %g)'%(ev, log_ev)
            html_model.p('Elliptical region evidence is %g, or log(Evidence) = %g.'%(ev, log_ev))

            if 'logl' in pos.names:
                log_ev=pos.harmonic_mean_evidence()
                html_model.p('Compare to harmonic mean evidence of %g (log(Evidence = %g))'%(exp(log_ev), log_ev))
        except IndexError:
            print 'Warning: Sample size too small to compute elliptical evidence!'

    #Create a section for SNR, if a file is provided
    if snrfactor is not None:
        html_snr=html.add_section('Signal to noise ratio(s)',legend=legend)
        html_snr.p('%s'%snrstring)

    #Create a section for summary statistics
    # By default collapse section are collapsed when the page is opened or reloaded. Use start_closed=False option as here below to change this
    tabid='statstable'
    html_stats=html.add_collapse_section('Summary statistics',legend=legend,innertable_id=tabid,start_closed=False)
    html_stats.write(str(pos))
    statfilename=os.path.join(outdir,"summary_statistics.dat")
    statout=open(statfilename,"w")
    statout.write("\tmaP\tmaxL\tstdev\tmean\tmedian\tstacc\tinjection\tvalue\n")
    
    for statname,statoned_pos in pos:

      statmax_pos,max_i=pos._posMaxL()
      statmaxL=statoned_pos.samples[max_i][0]
      statmax_pos,max_j=pos._posMap()
      statmaxP=statoned_pos.samples[max_j][0]
      statmean=str(statoned_pos.mean)
      statstdev=str(statoned_pos.stdev)
      statmedian=str(squeeze(statoned_pos.median))
      statstacc=str(statoned_pos.stacc)
      statinjval=str(statoned_pos.injval)
      
      statarray=[str(i) for i in [statname,statmaxP,statmaxL,statstdev,statmean,statmedian,statstacc,statinjval]]
      statout.write("\t".join(statarray))
      statout.write("\n")
      
    statout.close()

    #==================================================================#
    #Generate sky map, WF, and PSDs
    #==================================================================#
   
    skyreses=None
    sky_injection_cl=None
    inj_position=None
    tabid='skywftable'
    html_wf=html.add_collapse_section('Sky Localization and Waveform',innertable_id=tabid)
    
    table=html_wf.tab(idtable=tabid)
    row=html_wf.insert_row(table,label='SkyandWF')
    skytd=html_wf.insert_td(row,'',label='SkyMap',legend=legend)
    html_sky=html.add_section_to_element('SkyMap',skytd)
    #If sky resolution parameter has been specified try and create sky map...
    if skyres is not None and \
       (('ra' in pos.names and 'dec' in pos.names) or \
        ('rightascension' in pos.names and 'declination' in pos.names)):
        #Greedy bin sky samples (ra,dec) into a grid on the sky which preserves
        #?
        if ('ra' in pos.names and 'dec' in pos.names):
            inj_position=[pos['dec'].injval,pos['ra'].injval]
        elif ('rightascension' in pos.names and 'declination' in pos.names):
            inj_position=[pos['declination'].injval,pos['rightascension'].injval]
        top_ranked_sky_pixels,sky_injection_cl,skyreses,injection_area=bppu.greedy_bin_sky(pos,skyres,confidence_levels)
        print "BCI for sky area:"
        print skyreses
        #Create sky map in outdir
        bppu.plot_sky_map(inj_position,top_ranked_sky_pixels,outdir)
        
        #Create a web page section for sky localization results/plots (if defined)

        #html_sky=html.add_section('Sky Localization',legend=legend)
        if injection:
            if sky_injection_cl:
                html_sky.p('Injection found at confidence interval %f in sky location'%(sky_injection_cl))
            else:
                html_sky.p('Injection not found in posterior bins in sky location!')
        html_sky.write('<a href="skymap.png" target="_blank"><img src="skymap.png"/></a>')

        html_sky_write='<table border="1" id="statstable"><tr><th>Confidence region</th><th>size (sq. deg)</th></tr>'

        fracs=skyreses.keys()
        fracs.sort()

        skysizes=[skyreses[frac] for frac in fracs]
        for frac,skysize in zip(fracs,skysizes):
            html_sky_write+=('<tr><td>%f</td><td>%f</td></tr>'%(frac,skysize))
        html_sky_write+=('</table>')

        html_sky.write(html_sky_write)
    else:
        html_sky.write('<b>No skymap generated!</b>')
    wfdir=os.path.join(outdir,'Waveform')
    if not os.path.isdir(wfdir):
        os.makedirs(wfdir)
    try:
        wfpointer= bppu.plot_waveform(pos=pos,siminspiral=injfile,event=eventnum,path=wfdir)
    except:
        wfpointer = None
    wftd=html_wf.insert_td(row,'',label='Waveform',legend=legend)
    wfsection=html.add_section_to_element('Waveforms',wftd)
    if wfpointer:
      wfsection.write('<a href="Waveform/WF_DetFrame.png" target="_blank"><img src="Waveform/WF_DetFrame.png"/></a>')
    else:
      wfsection.write("<b>No Waveform generated!</b>")
      
    wftd=html_wf.insert_td(row,'',label='PSDs',legend=legend)
    wfsection=html.add_section_to_element('PSDs',wftd)
    psd_pointer=None
    if psd_files is not None:
      psd_files=list(psd_files.split(','))
      psddir=os.path.join(outdir,'PSDs')
      if not os.path.isdir(psddir):
        os.makedirs(psddir)
      try:
        psd_pointer=bppu.plot_psd(psd_files,outpath=psddir)    
      except:
        psd_pointer=None
    if psd_pointer:
      wfsection.write('<a href="PSDs/PSD.png" target="_blank"><img src="PSDs/PSD.png"/></a>')
    else:
      wfsection.write("<b>No PSD file found!</b>")
    #==================================================================#
    #1D posteriors
    #==================================================================#

    #Loop over each parameter and determine the contigious and greedy
    #confidence levels and some statistics.

    #Add section for 1D marginal PDFs and sample plots
    tabid='onedmargtable'
    html_ompdf=html.add_collapse_section('1D marginal posterior PDFs',legend=legend,innertable_id=tabid)
    #Table matter
    if not noacf:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th><th>Autocorrelation</th></tr>'%tabid
    else:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th></tr>'%tabid

#Add section for 1D confidence intervals
    tabid='onedconftable'
    html_ogci=html.add_collapse_section('1D confidence intervals (greedy binning)',legend=legend,innertable_id=tabid)
    html_ogci_write='<table id="%s" border="1"><tr><th/>'%tabid
    clasciiout="#parameter \t"
    confidence_levels.sort()
    for cl in confidence_levels:
        html_ogci_write+='<th>%f</th>'%cl
        clasciiout+="%s\t"%('%.02f'%cl)
    if injection:
        html_ogci_write+='<th>Injection Confidence Level</th>'
        html_ogci_write+='<th>Injection Confidence Interval</th>'
        clasciiout+="Injection_Confidence_Level\t"
        clasciiout+="Injection_Confidence_Interval"
    clasciiout+='\n'
    html_ogci_write+='</tr>'

    onepdfdir=os.path.join(outdir,'1Dpdf')
    if not os.path.isdir(onepdfdir):
        os.makedirs(onepdfdir)

    sampsdir=os.path.join(outdir,'1Dsamps')
    if not os.path.isdir(sampsdir):
        os.makedirs(sampsdir)
    Nskip=0
    if 'chain' in pos.names:
        data,header=pos.samples()
        par_index=pos.names.index('cycle')
        chain_index=pos.names.index("chain")
        chains=unique(pos["chain"].samples)
        chainCycles = [sort(data[ data[:,chain_index] == chain, par_index ]) for chain in chains]
        chainNcycles = [cycles[-1]-cycles[0] for cycles in chainCycles]
        chainNskips = [cycles[1] - cycles[0] for cycles in chainCycles]
    elif 'cycle' in pos.names:
        cycles = sort(pos['cycle'].samples)
        Ncycles = cycles[-1]-cycles[0]
        Nskip = cycles[1]-cycles[0]
    
    printed=0
    for par_name in oneDMenu:
        par_name=par_name.lower()
        try:
            pos[par_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par_name
            continue
        try:
            par_bin=GreedyRes[par_name]
        except KeyError:
            print "Bin size is not set for %s, skipping binning."%par_name
            continue

        #print "Binning %s to determine confidence levels ..."%par_name
        binParams={par_name:par_bin}
        injection_area=None
        injection_area=None
        if greedy:
          if printed==0:
            print "Using greedy 1-d binning credible regions\n"
            printed=1
          toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(pos,binParams,confidence_levels)
        else:
          if printed==0:
            print "Using 2-step KDE 1-d credible regions\n"
            printed=1
          if pos[par_name].injval is None:
            injCoords=None
          else:
            injCoords=[pos[par_name].injval]
          _,reses,injstats=bppu.kdtree_bin2Step(pos,[par_name],confidence_levels,injCoords=injCoords)
          if injstats is not None:
            injectionconfidence=injstats[3]
            injection_area=injstats[4]
             
        #oneDContCL,oneDContInj = bppu.contigious_interval_one_param(pos,binParams,confidence_levels)

        #Generate new BCI html table row
        BCItableline='<tr><td>%s</td>'%(par_name)
        clasciiout+="%s\t"%par_name
        cls=reses.keys()
        cls.sort()

        for cl in cls:
            BCItableline+='<td>%f</td>'%reses[cl]
            clasciiout+="%f\t"%reses[cl]
        if injection is not None:
            if injectionconfidence is not None and injection_area is not None:

                BCItableline+='<td>%f</td>'%injectionconfidence
                BCItableline+='<td>%f</td>'%injection_area
                clasciiout+="%f\t"%injectionconfidence
                clasciiout+="%f"%injection_area

            else:
                BCItableline+='<td/>'
                BCItableline+='<td/>'
                clasciiout+="nan\t"
                clasciiout+="nan"
        BCItableline+='</tr>'
        clasciiout+="\n"
        #Append new table line to section html
        html_ogci_write+=BCItableline

        #Generate 1D histogram/kde plots
        print "Generating 1D plot for %s."%par_name

        #Get analytic description if given
        pdf=cdf=None
        if analyticLikelihood:
            pdf = analyticLikelihood.pdf(par_name)
            cdf = analyticLikelihood.cdf(par_name)

        oneDPDFParams={par_name:50}
        try:
            rbins,plotFig=bppu.plot_one_param_pdf(pos,oneDPDFParams,pdf,cdf,plotkde=False)
        except:
            print "Failed to produce plot for %s."%par_name
            continue

        figname=par_name+'.png'
        oneDplotPath=os.path.join(onepdfdir,figname)
        plotFig.savefig(oneDplotPath)
        if(savepdfs): plotFig.savefig(os.path.join(onepdfdir,par_name+'.pdf'))
        plt.close(plotFig)

        if rbins:
            print "r of injected value of %s (bins) = %f"%(par_name, rbins)

        ##Produce plot of raw samples
        myfig=plt.figure(figsize=(4,3.5),dpi=200)
        pos_samps=pos[par_name].samples
        if not ("chain" in pos.names):
            # If there is not a parameter named "chain" in the
            # posterior, then just produce a plot of the samples.
            plt.plot(pos_samps,'k,',linewidth=0.0, markeredgewidth=0,figure=myfig)
            maxLen=len(pos_samps)
        else:
            # If there is a parameter named "chain", then produce a
            # plot of the various chains in different colors, with
            # smaller dots.
            data,header=pos.samples()
            par_index=pos.names.index(par_name)
            chain_index=pos.names.index("chain")
            chains=unique(pos["chain"].samples)
            chainData=[data[ data[:,chain_index] == chain, par_index ] for chain in chains]
            chainDataRanges=[range(len(cd)) for cd in chainData]
            maxLen=max([len(cd) for cd in chainData])
            for rng, data in zip(chainDataRanges, chainData):
                plt.plot(rng, data, marker=',',linewidth=0.0, markeredgewidth=0,figure=myfig)
            plt.title("Gelman-Rubin R = %g"%(pos.gelman_rubin(par_name)))

            #dataPairs=[ [rng, data] for (rng,data) in zip(chainDataRanges, chainData)]
            #flattenedData=[ item for pair in dataPairs for item in pair ]
            #maxLen=max([len(data) for data in flattenedData])
            #plt.plot(array(flattenedData),marker=',',linewidth=0.0,figure=myfig)


        injpar=pos[par_name].injval

        if injpar is not None:
            # Allow injection to be 5% outside the posterior plot
            minrange=min(pos_samps)-0.05*(max(pos_samps)-min(pos_samps))
            maxrange=max(pos_samps)+0.05*(max(pos_samps)-min(pos_samps))
            if minrange<injpar and maxrange>injpar:
                plt.axhline(injpar, color='r', linestyle='-.',linewidth=4)
        myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.png')))
        if(savepdfs): myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.pdf')))
        plt.close(myfig)
        acfail=0
        if not (noacf):
            acffig=plt.figure(figsize=(4,3.5),dpi=200)
            if not ("chain" in pos.names):
                data=pos_samps[:,0]
                try:
                    (Neff, acl, acf) = bppu.effectiveSampleSize(data, Nskip)
                    lines=plt.plot(acf, 'k,', marker=',',linewidth=0.0, markeredgewidth=0, figure=acffig)
                    # Give ACL info if not already downsampled according to it
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    elif 'cycle' in pos.names:
                        last_color = lines[-1].get_color()
                        plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                        plt.title('ACL = %i   N = %i'%(acl,Neff))
                except FloatingPointError:
                    # Ignore
                    acfail=1
                    pass
            else:
                try:
                    acls = []
                    Nsamps = 0.0;
                    for rng, data, Nskip, Ncycles in zip(chainDataRanges, chainData, chainNskips, chainNcycles):
                        (Neff, acl, acf) = bppu.effectiveSampleSize(data, Nskip)
                        acls.append(acl)
                        Nsamps += Neff
                        lines=plt.plot(acf,'k,', marker=',',linewidth=0.0, markeredgewidth=0, figure=acffig)
                        # Give ACL info if not already downsampled according to it
                        if nDownsample is not None:
                            last_color = lines[-1].get_color()
                            plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    else:
                        plt.title('ACL = %i  N = %i'%(max(acls),Nsamps))
                except FloatingPointError:
                    # Ignore
                    acfail=1
                    pass

            acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.png')))
            if(savepdfs): acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.pdf')))
            plt.close(acffig)

        if not noacf:
	  if not acfail:
	    acfhtml='<td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png', '_acf.png')+'"/></td>'
	  else:
	    acfhtml='<td>ACF generation failed!</td>'
          html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td>'+acfhtml+'</tr>'
        else:
            html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td></tr>'


    html_ompdf_write+='</table>'
    html_ompdf.write(html_ompdf_write)

    html_ogci_write+='</table>'
    html_ogci.write(html_ogci_write)

    cornerdir=os.path.join(outdir,'corner')
    if not os.path.isdir(cornerdir):
        os.makedirs(cornerdir)

    #===============================#
    # Corner plots
    #===============================#

    massParams=['mtotal','m1','m2','mc']
    distParams=['distance','distMPC','dist']
    incParams=['iota','inclination','theta_jn']
    polParams=['psi','polarisation','polarization']
    skyParams=['ra','rightascension','declination','dec']
    timeParams=['time']
    spinParams=['spin1','spin2','a1','a2','phi1','theta1','phi2','theta2','chi','effectivespin','beta','tilt1','tilt2','phi_jl','theta_jn','phi12']
    intrinsicParams=massParams+spinParams
    extrinsicParams=incParams+distParams+polParams+skyParams
    try:
      myfig=bppu.plot_corner(pos,[0.05,0.5,0.95],parnames=intrinsicParams)
    except:
      myfig=None
    tabid='CornerTable'
    html_corner=''
    got_any=0
    if myfig:
      html_corner+='<tr><td width="100%"><a href="corner/intrinsic.png" target="_blank"><img width="70%" src="corner/intrinsic.png"/></a></td></tr>'
      myfig.savefig(os.path.join(cornerdir,'intrinsic.png'))
      myfig.savefig(os.path.join(cornerdir,'intrinsic.pdf'))
      got_any+=1
    try:
      myfig=bppu.plot_corner(pos,[0.05,0.5,0.95],parnames=extrinsicParams)
    except:
      myfig=None

    if myfig:
      myfig.savefig(os.path.join(cornerdir,'extrinsic.png'))
      myfig.savefig(os.path.join(cornerdir,'extrinsic.pdf'))
      html_corner+='<tr><td width="100%"><a href="corner/extrinsic.png" target="_blank"><img width="70%" src="corner/extrinsic.png"/></a></td></tr>'
      got_any+=1
    if got_any>0:
      html_corner='<table id="%s" border="1">'%tabid+html_corner
      html_corner+='</table>'
    if html_corner!='':
      html_co=html.add_collapse_section('Corner plots',legend=legend,innertable_id=tabid)
      html_co.write(html_corner)
    if clasciiout:
      fout=open(os.path.join(outdir,'confidence_levels.txt'),'w')
      fout.write(clasciiout)
      fout.close()
    #==================================================================#
    #2D posteriors
    #==================================================================#

    #Loop over parameter pairs in twoDGreedyMenu and bin the sample pairs
    #using a greedy algorithm . The ranked pixels (toppoints) are used
    #to plot 2D histograms and evaluate Bayesian confidence intervals.

    #Make a folder for the 2D kde plots
    margdir=os.path.join(outdir,'2Dkde')
    if not os.path.isdir(margdir):
        os.makedirs(margdir)

    twobinsdir=os.path.join(outdir,'2Dbins')
    if not os.path.isdir(twobinsdir):
        os.makedirs(twobinsdir)

    greedytwobinsdir=os.path.join(outdir,'greedy2Dbins')
    if not os.path.isdir(greedytwobinsdir):
        os.makedirs(greedytwobinsdir)

    #Add a section to the webpage for a table of the confidence interval
    #results.
    tabid='2dconftable'
    html_tcig=html.add_collapse_section('2D confidence intervals (greedy binning)',legend=legend,innertable_id=tabid)
    #Generate the top part of the table
    html_tcig_write='<table id="%s" border="1"><tr><th/>'%tabid
    confidence_levels.sort()
    for cl in confidence_levels:
        html_tcig_write+='<th>%f</th>'%cl
    if injection:
        html_tcig_write+='<th>Injection Confidence Level</th>'
        html_tcig_write+='<th>Injection Confidence Interval</th>'
    html_tcig_write+='</tr>'


    #=  Add a section for a table of 2D marginal PDFs (kde)
    twodkdeplots_flag=twodkdeplots
    if twodkdeplots_flag:
        tabid='2dmargtable'
        html_tcmp=html.add_collapse_section('2D Marginal PDFs',legend=legend,innertable_id=tabid)
        #Table matter
        html_tcmp_write='<table border="1" id="%s">'%tabid

    tabid='2dgreedytable'
    html_tgbh=html.add_collapse_section('2D Greedy Bin Histograms',legend=legend,innertable_id=tabid)
    html_tgbh_write='<table border="1" id="%s">'%tabid

    row_count=0
    row_count_gb=0

    for par1_name,par2_name in twoDGreedyMenu:
        par1_name=par1_name.lower()
        par2_name=par2_name.lower()
        # Don't plot a parameter against itself!
        if par1_name == par2_name: continue
        try:
            pos[par1_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par1_name
            continue
        try:
            pos[par2_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par2_name
            continue
        #Bin sizes
        try:
            par1_bin=GreedyRes[par1_name]
        except KeyError:
            print "Bin size is not set for %s, skipping %s/%s binning."%(par1_name,par1_name,par2_name)
            continue
        try:
            par2_bin=GreedyRes[par2_name]
        except KeyError:
            print "Bin size is not set for %s, skipping %s/%s binning."%(par2_name,par1_name,par2_name)
            continue

        # Skip any fixed parameters to avoid matrix inversion problems
        par1_pos=pos[par1_name].samples
        par2_pos=pos[par2_name].samples
        if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
            continue

        #print "Binning %s-%s to determine confidence levels ..."%(par1_name,par2_name)
        #Form greedy binning input structure
        greedy2Params={par1_name:par1_bin,par2_name:par2_bin}

        #Greedy bin the posterior samples
        toppoints,injection_cl,reses,injection_area=\
        bppu.greedy_bin_two_param(pos,greedy2Params,confidence_levels)

        print "BCI %s-%s:"%(par1_name,par2_name)
        print reses

        #Generate new BCI html table row
        BCItableline='<tr><td>%s-%s</td>'%(par1_name,par2_name)
        cls=reses.keys()
        cls.sort()

        for cl in cls:
            BCItableline+='<td>%f</td>'%reses[cl]

        if injection is not None:
            if injection_cl is not None:
                BCItableline+='<td>%f</td>'%injection_cl
                BCItableline+='<td>'+str(injection_area)+'</td>'

            else:
                BCItableline+='<td/>'
                BCItableline+='<td/>'

        BCItableline+='</tr>'

        #Append new table line to section html
        html_tcig_write+=BCItableline


        #= Plot 2D histograms of greedily binned points =#

        #greedy2ContourPlot=bppu.plot_two_param_greedy_bins_contour({'Result':pos},greedy2Params,[0.67,0.9,0.95],{'Result':'k'})
        greedy2ContourPlot=bppu.plot_two_param_kde_greedy_levels({'Result':pos},greedy2Params,[0.67,0.9,0.95],{'Result':'k'})
        greedy2contourpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2contour.png'%(par1_name,par2_name))
        if greedy2ContourPlot is not None:
          greedy2ContourPlot.savefig(greedy2contourpath)
          if(savepdfs): greedy2ContourPlot.savefig(greedy2contourpath.replace('.png','.pdf'))
          plt.close(greedy2ContourPlot)

        greedy2HistFig=bppu.plot_two_param_greedy_bins_hist(pos,greedy2Params,confidence_levels)
        greedy2histpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2.png'%(par1_name,par2_name))
        greedy2HistFig.savefig(greedy2histpath)
        if(savepdfs): greedy2HistFig.savefig(greedy2histpath.replace('.png','.pdf'))
        plt.close(greedy2HistFig)

        greedyFile = open(os.path.join(twobinsdir,'%s_%s_greedy_stats.txt'%(par1_name,par2_name)),'w')

        #= Write out statistics for greedy bins
        for cl in cls:
            greedyFile.write("%lf %lf\n"%(cl,reses[cl]))
        greedyFile.close()

        if [par1_name,par2_name] in twoDplots or [par2_name,par1_name] in twoDplots :
            print 'Generating %s-%s greedy hist plot'%(par1_name,par2_name)

            par1_pos=pos[par1_name].samples
            par2_pos=pos[par2_name].samples

            if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
                continue
            head,figname=os.path.split(greedy2histpath)
            head,figname_c=os.path.split(greedy2contourpath)
            if row_count_gb==0:
                html_tgbh_write+='<tr>'
            html_tgbh_write+='<td width="30%"><img width="100%" src="greedy2Dbins/'+figname+'"/>[<a href="greedy2Dbins/'+figname_c+'">contour</a>]</td>'
            row_count_gb+=1
            if row_count_gb==3:
                html_tgbh_write+='</tr>'
                row_count_gb=0

        #= Generate 2D kde plots =#

        if twodkdeplots_flag is True:
            if [par1_name,par2_name] in twoDplots or [par2_name,par1_name] in twoDplots :
                print 'Generating %s-%s plot'%(par1_name,par2_name)

                par1_pos=pos[par1_name].samples
                par2_pos=pos[par2_name].samples

                if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
                    continue

                plot2DkdeParams={par1_name:50,par2_name:50}
                myfig=bppu.plot_two_param_kde(pos,plot2DkdeParams)

                figname=par1_name+'-'+par2_name+'_2Dkernel.png'
                twoDKdePath=os.path.join(margdir,figname)

                if row_count==0:
                    html_tcmp_write+='<tr>'
                html_tcmp_write+='<td width="30%"><img width="100%" src="2Dkde/'+figname+'"/></td>'
                row_count+=1
                if row_count==3:
                    html_tcmp_write+='</tr>'
                    row_count=0

                myfig.savefig(twoDKdePath)
                if(savepdfs): myfig.savefig(twoDKdePath.replace('.png','.pdf'))
                plt.close(myfig)


    #Finish off the BCI table and write it into the etree
    html_tcig_write+='</table>'
    html_tcig.write(html_tcig_write)

    if twodkdeplots_flag is True:
    #Finish off the 2D kde plot table
        while row_count!=0:
            html_tcmp_write+='<td/>'
            row_count+=1
            if row_count==3:
                row_count=0
                html_tcmp_write+='</tr>'
        html_tcmp_write+='</table>'
        html_tcmp.write(html_tcmp_write)
        #Add a link to all plots
        html_tcmp.a("2Dkde/",'All 2D marginal PDFs (kde)')

    #Finish off the 2D greedy histogram plot table
    while row_count_gb!=0:
        html_tgbh_write+='<td/>'
        row_count_gb+=1
        if row_count_gb==3:
            row_count_gb=0
            html_tgbh_write+='</tr>'
    html_tgbh_write+='</table>'
    html_tgbh.write(html_tgbh_write)
    #Add a link to all plots
    html_tgbh.a("greedy2Dbins/",'All 2D Greedy Bin Histograms')

    if RconvergenceTests is True:
        convergenceResults=bppu.convergenceTests(pos,gelman=False)
        
        if convergenceResults is not None:
            tabid='convtable'
            html_conv_test=html.add_collapse_section('Convergence tests',legend=legend,innertable_id=tabid)
            data_found=False
            for test,test_data in convergenceResults.items():
                
                if test_data:
                    data_found=True
                    html_conv_test.h3(test)
                                       
                    html_conv_table_rows={}
                    html_conv_table_header=''
                    for chain,chain_data in test_data.items():
                        html_conv_table_header+='<th>%s</th>'%chain
                        
                        
                        for data in chain_data:
                            if len(data)==2:
                                try:
                                    html_conv_table_rows[data[0]]+='<td>'+data[1]+'</td>'
                                except KeyError:
                                    html_conv_table_rows[data[0]]='<td>'+data[1]+'</td>'
                                
                    html_conv_table='<table id="%s"><tr><th>Chain</th>'%tabid+html_conv_table_header+'</tr>'
                    for row_name,row in html_conv_table_rows.items():
                        html_conv_table+='<tr><td>%s</td>%s</tr>'%(row_name,row)
                    html_conv_table+='</table>'
                    html_conv_test.write(html_conv_table)
            if data_found is False:
                html_conv_test.p('No convergence diagnostics generated!')

    #Create a section for the covariance matrix
    tabid='covtable'
    html_stats_cov=html.add_collapse_section('Covariance matrix',legend=legend,innertable_id=tabid)
    pos_samples,table_header_string=pos.samples()

    #calculate cov matrix
    try:
        cov_matrix=cov(pos_samples,rowvar=0,bias=1)

    #Create html table
        table_header_list=table_header_string.split()
        cov_table_string='<table border="1" id="%s"><tr><th/>'%tabid
        for header in table_header_list:
            cov_table_string+='<th>%s</th>'%header
        cov_table_string+='</tr>'

        cov_column_list=hsplit(cov_matrix,cov_matrix.shape[1])

        for cov_column,cov_column_name in zip(cov_column_list,table_header_list):
            cov_table_string+='<tr><th>%s</th>'%cov_column_name
            for cov_column_element in cov_column:
                cov_table_string+='<td>%.3e</td>'%(cov_column_element[0])
            cov_table_string+='</tr>'
        cov_table_string+='</table>'
        html_stats_cov.write(cov_table_string)

    except:
        print 'Unable to compute the covariance matrix.'

    #Create a section for run configuration information if it exists
    if pos._votfile is not None:
        html_vot=html.add_section('Run information',legend=legend)
        html_vot.write(pos.write_vot_info())
    
    html_footer=html.add_section('')
    html_footer.p('Produced using cbcBayesPostProc.py at '+strftime("%Y-%m-%d %H:%M:%S")+' .')

    cc_args=''
    for arg in sys.argv:
        cc_args+=arg+' '
        
    html_footer.p('Command line: %s'%cc_args)
    html_footer.p(git_version.verbose_msg)

    #Save results page
    resultspage=open(os.path.join(outdir,'posplots.html'),'w')
    resultspage.write(str(html))

    # Save posterior samples too...
    posfilename=os.path.join(outdir,'posterior_samples.dat')
    pos.write_to_file(posfilename)

    #Close files
    resultspage.close()

USAGE='''%prog [options] datafile.dat [datafile2.dat ...]
Generate a web page displaying results of parameter estimation based on the contents
of one or more datafiles containing samples from one of the bayesian algorithms (MCMC, nested sampling).
Options specify which extra statistics to compute and allow specification of additional information.
'''

if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser(USAGE)
    parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
    parser.add_option("-d","--data",dest="data",action="callback",callback=multipleFileCB,help="datafile")
    #Optional (all)
    parser.add_option("-i","--inj",dest="injfile",help="SimInsipral injection file",metavar="INJ.XML",default=None)
    parser.add_option("-t","--trig",dest="trigfile",help="Coinc XML file",metavar="COINC.XML",default=None)
    parser.add_option("--skyres",dest="skyres",help="Sky resolution to use to calculate sky box size",default=None)
    parser.add_option("--eventnum",dest="eventnum",action="store",default=None,help="event number in SimInspiral file of this signal",type="int",metavar="NUM")
    parser.add_option("--trignum",dest="trignum",action="store",default=None,help="trigger number in CoincTable",type="int",metavar="NUM")
    parser.add_option("--bsn",action="store",default=None,help="Optional file containing the bayes factor signal against noise",type="string")
    parser.add_option("--bci",action="store",default=None,help="Optional file containing the bayes factor coherent signal model against incoherent signal model.",type="string")
    parser.add_option("--snr",action="store",default=None,help="Optional file containing the SNRs of the signal in each IFO",type="string")
    parser.add_option("--dievidence",action="store_true",default=False,help="Calculate the direct integration evidence for the posterior samples")
    parser.add_option("--boxing",action="store",default=64,help="Boxing parameter for the direct integration evidence calculation",type="int",dest="boxing")
    parser.add_option("--evidenceFactor",action="store",default=1.0,help="Overall factor (normalization) to apply to evidence",type="float",dest="difactor",metavar="FACTOR")
    
    parser.add_option('--ellipticEvidence', action='store_true', default=False,help='Estimate the evidence by fitting ellipse to highest-posterior points.', dest='ellevidence')

    parser.add_option("--no2D",action="store_true",default=False,help="Skip 2-D plotting.")
    parser.add_option("--header",action="store",default=None,help="Optional file containing the header line for posterior samples",type="string")
    #NS
    parser.add_option("--ns",action="store_true",default=False,help="(inspnest) Parse input as if it was output from parallel nested sampling runs.")
    parser.add_option("--Nlive",action="store",default=None,help="(inspnest) Number of live points used in each parallel nested sampling run.",type="int")
    parser.add_option("--xflag",action="store_true",default=False,help="(inspnest) Convert x to iota.")
    #SS
    parser.add_option("--ss",action="store_true",default=False,help="(SPINspiral) Parse input as if it was output from SPINspiral.")
    parser.add_option("--spin",action="store_true",default=False,help="(SPINspiral) Specify spin run (15 parameters). ")
    #LALInf
    parser.add_option("--lalinfmcmc",action="store_true",default=False,help="(LALInferenceMCMC) Parse input from LALInferenceMCMC.")
    parser.add_option("--inj-spin-frame",default='OrbitalL', help="The reference frame used for the injection (default: OrbitalL)")
    parser.add_option("--downsample",action="store",default=None,help="(LALInferenceMCMC) approximate number of samples to record in the posterior",type="int")
    parser.add_option("--deltaLogL",action="store",default=None,help="(LALInferenceMCMC) Difference in logL to use for convergence test.",type="float")
    parser.add_option("--fixedBurnin",dest="fixedBurnin",action="callback",callback=multipleFileCB,help="(LALInferenceMCMC) Fixed number of iteration for burnin.")
    parser.add_option("--oldMassConvention",action="store_true",default=False,help="(LALInferenceMCMC) if activated, m2 > m1; otherwise m1 > m2 in PTMCMC.output.*.00")
    #FM
    parser.add_option("--fm",action="store_true",default=False,help="(followupMCMC) Parse input as if it was output from followupMCMC.")
    # ACF plots off?
    parser.add_option("--no-acf", action="store_true", default=False, dest="noacf")
    # Turn on 2D kdes
    parser.add_option("--twodkdeplots", action="store_true", default=False, dest="twodkdeplots")
    # Turn on R convergence tests
    parser.add_option("--RconvergenceTests", action="store_true", default=False, dest="RconvergenceTests")
    parser.add_option("--nopdfs",action="store_false",default=True,dest="nopdfs")
    parser.add_option("-c","--covarianceMatrix",dest="covarianceMatrices",action="append",default=None,help="CSV file containing covariance (must give accompanying mean vector CSV. Can add more than one matrix.")
    parser.add_option("-m","--meanVectors",dest="meanVectors",action="append",default=None,help="Comma separated list of locations of the multivariate gaussian described by the correlation matrix.  First line must be list of params in the order used for the covariance matrix.  Provide one list per covariance matrix.")
    parser.add_option("--email",action="store",default=None,type="string",metavar="user@ligo.org",help="Send an e-mail to the given address with a link to the finished page.")
    parser.add_option("--archive",action="store",default=None,type="string",metavar="results.tar.gz",help="Create the given tarball with all results")
    parser.add_option("--psdfiles",action="store",default=None,type="string",metavar="H1,L1,V1",help="comma separater list of ASCII files with PSDs, one per IFO")
    parser.add_option("--kdecredibleregions",action="store_true",default=False,help="If given, will use 2-step KDE trees to estimate 1-d credible regions [default false: use greedy binning]")
    (opts,args)=parser.parse_args()

    datafiles=[]
    if args:
      datafiles=datafiles+args
    if opts.data:
      datafiles=datafiles + opts.data
    
    if opts.fixedBurnin:
      fixedBurnins = [int(fixedBurnin) for fixedBurnin in opts.fixedBurnin]
    else:
      fixedBurnins = None

    #List of parameters to plot/bin . Need to match (converted) column names.
    massParams=['m1','m2','chirpmass','mchirp','mc','eta','q','massratio','asym_massratio','mtotal']
    distParams=['distance','distMPC','dist']
    incParams=['iota','inclination','cosiota']
    polParams=['psi','polarisation','polarization']
    skyParams=['ra','rightascension','declination','dec']
    timeParams=['time','time_mean']
    spinParams=['spin1','spin2','a1','a2','phi1','theta1','phi2','theta2','costilt1','costilt2','chi','effectivespin','costheta_jn','cosbeta','tilt1','tilt2','phi_jl','theta_jn','phi12']
    phaseParams=['phase', 'phi0','phase_maxl']
    endTimeParams=['l1_end_time','h1_end_time','v1_end_time']
    ppEParams=['ppEalpha','ppElowera','ppEupperA','ppEbeta','ppElowerb','ppEupperB','alphaPPE','aPPE','betaPPE','bPPE']
    tigerParams=['dphi%i'%(i) for i in range(7)] + ['dphi%il'%(i) for i in [5,6] ]
    bransDickeParams=['omegaBD','ScalarCharge1','ScalarCharge2']
    massiveGravitonParams=['lambdaG']
    tidalParams=['lambda1','lambda2','lam_tilde','dlam_tilde','lambdat','dlambdat']
    statsParams=['logprior','logl','deltalogl','deltaloglh1','deltalogll1','deltaloglv1','deltaloglh2','deltaloglg1','flow']
    snrParams=['optimal_snr','matched_filter_snr']
    oneDMenu=massParams + distParams + incParams + polParams + skyParams + timeParams + spinParams + phaseParams + endTimeParams + ppEParams + tigerParams + bransDickeParams + massiveGravitonParams + tidalParams + statsParams + snrParams
    # ['mtotal','m1','m2','chirpmass','mchirp','mc','distance','distMPC','dist','iota','inclination','psi','eta','massratio','ra','rightascension','declination','dec','time','a1','a2','phi1','theta1','phi2','theta2','costilt1','costilt2','chi','effectivespin','phase','l1_end_time','h1_end_time','v1_end_time']
    ifos_menu=['h1','l1','v1']
    from itertools import combinations
    for ifo1,ifo2 in combinations(ifos_menu,2):
      oneDMenu.append(ifo1+ifo2+'_delay')
    #oneDMenu=[]
    twoDGreedyMenu=[]
    if not opts.no2D:
        for mp1,mp2 in combinations(massParams,2):
          twoDGreedyMenu.append([mp1, mp2])
        for mp in massParams:
            for d in distParams:
                twoDGreedyMenu.append([mp,d])
        for mp in massParams:
            for sp in spinParams:
                twoDGreedyMenu.append([mp,sp])
        for dp in distParams:
            for sp in snrParams:
                twoDGreedyMenu.append([dp,sp])
        for dp in distParams:
            for ip in incParams:
                twoDGreedyMenu.append([dp,ip])
        for dp in distParams:
            for sp in skyParams:
                twoDGreedyMenu.append([dp,sp])
        for dp in distParams:
            for sp in spinParams:
                twoDGreedyMenu.append([dp,sp])
        for ip in incParams:
            for sp in skyParams:
                twoDGreedyMenu.append([ip,sp])
        for ip in incParams:
            for sp in spinParams:
                twoDGreedyMenu.append([ip,sp])
            for phip in phaseParams:
                twoDGreedyMenu.append([ip,phip])
            for psip in polParams:
                twoDGreedyMenu.append([ip,psip])
        for sp1 in skyParams:
            for sp2 in skyParams:
                if not (sp1 == sp2):
                    twoDGreedyMenu.append([sp1, sp2])
        for sp1,sp2 in combinations(spinParams,2):
          twoDGreedyMenu.append([sp1, sp2])
        for mp in massParams:
             for tp in tidalParams:
                 if not (mp == tp):
                     twoDGreedyMenu.append([mp, tp])
        twoDGreedyMenu.append(['lambda1','lambda2'])
        twoDGreedyMenu.append(['lam_tilde','dlam_tilde'])
        twoDGreedyMenu.append(['lambdat','dlambdat'])
        twoDGreedyMenu.append(['optimal_snr','matched_filter_snr'])
        for psip in polParams:
            for phip in phaseParams:
                twoDGreedyMenu.append([psip,phip])
            for sp in skyParams:
                twoDGreedyMenu.append([psip,sp])
            for sp in spinParams:
                twoDGreedyMenu.append([psip,sp])

    #twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['ra','dec']]
    #Bin size/resolution for binning. Need to match (converted) column names.
    greedyBinSizes={'mc':0.025,'m1':0.1,'m2':0.1,'mass1':0.1,'mass2':0.1,'mtotal':0.1,'eta':0.001,'q':0.01,'asym_massratio':0.01,'iota':0.01,'cosiota':0.02,'time':1e-4,'time_mean':1e-4,'distance':1.0,'dist':1.0,'mchirp':0.025,'chirpmass':0.025,'spin1':0.04,'spin2':0.04,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'ra':0.05,'dec':0.05,'chi':0.05,'costilt1':0.02,'costilt2':0.02,'thatas':0.05,'costheta_jn':0.02,'beta':0.05,'omega':0.05,'cosbeta':0.02,'ppealpha':1.0,'ppebeta':1.0,'ppelowera':0.01,'ppelowerb':0.01,'ppeuppera':0.01,'ppeupperb':0.01,'polarisation':0.04,'rightascension':0.05,'declination':0.05,'massratio':0.001,'inclination':0.01,'phase':0.05,'tilt1':0.05,'tilt2':0.05,'phi_jl':0.05,'theta_jn':0.05,'phi12':0.05,'flow':1.0,'phase_maxl':0.05,'optimal_snr':0.05,'matched_filter_snr':0.05}
    for derived_time in ['h1_end_time','l1_end_time','v1_end_time','h1l1_delay','l1v1_delay','h1v1_delay']:
        greedyBinSizes[derived_time]=greedyBinSizes['time']
    if not opts.no2D:
        for dt1,dt2 in combinations(['h1_end_time','l1_end_time','v1_end_time'],2):
          twoDGreedyMenu.append([dt1,dt2])
        for dt1,dt2 in combinations( ['h1l1_delay','l1v1_delay','h1v1_delay'],2):
          twoDGreedyMenu.append([dt1,dt2])
    for param in tigerParams + bransDickeParams + massiveGravitonParams:
        greedyBinSizes[param]=0.01
    for param in tidalParams:
        greedyBinSizes[param]=2.5
    #Confidence levels
    for loglname in ['logl','deltalogl','deltaloglh1','deltaloglv1','deltalogll1','logll1','loglh1','loglv1']:
        greedyBinSizes[loglname]=0.1
    confidenceLevels=[0.67,0.9,0.95,0.99]
    #2D plots list
    #twoDplots=[['mc','eta'],['mchirp','eta'],['mc', 'time'],['mchirp', 'time'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['ra', 'dec'],['m1','dist'],['m2','dist'],['mc', 'dist'],['psi','iota'],['psi','distance'],['psi','dist'],['psi','phi0'], ['a1', 'a2'], ['a1', 'iota'], ['a2', 'iota'],['eta','time'],['ra','iota'],['dec','iota'],['chi','iota'],['chi','mchirp'],['chi','eta'],['chi','distance'],['chi','ra'],['chi','dec'],['chi','psi']]
    twoDplots=twoDGreedyMenu
    cbcBayesPostProc(
                        opts.outpath,datafiles,oneDMenu,twoDGreedyMenu,
                        greedyBinSizes,confidenceLevels,twoDplots,
                        #optional
                        injfile=opts.injfile,eventnum=opts.eventnum,
                        trigfile=opts.trigfile,trignum=opts.trignum,
                        skyres=opts.skyres,
                        # direct integration evidence
                        dievidence=opts.dievidence,boxing=opts.boxing,difactor=opts.difactor,
                        # Ellipitical evidence
                        ellevidence=opts.ellevidence,
                        #manual bayes factor entry
                        bayesfactornoise=opts.bsn,bayesfactorcoherent=opts.bci,
                        #manual input for SNR in the IFOs, optional.
                        snrfactor=opts.snr,
                        #nested sampling options
                        ns_flag=opts.ns,ns_Nlive=opts.Nlive,
                        #spinspiral/mcmc options
                        ss_flag=opts.ss,ss_spin_flag=opts.spin,
                        #LALInferenceMCMC options
                        li_flag=opts.lalinfmcmc,deltaLogL=opts.deltaLogL,fixedBurnins=fixedBurnins,nDownsample=opts.downsample,oldMassConvention=opts.oldMassConvention,
                        #followupMCMC options
                        fm_flag=opts.fm,
                        #injected spin frame
                        inj_spin_frame=opts.inj_spin_frame,
                        # Turn of ACF?
                        noacf=opts.noacf,
                        #Turn on 2D kdes
                        twodkdeplots=opts.twodkdeplots,
                        #Turn on R convergence tests
                        RconvergenceTests=opts.RconvergenceTests,
                        # Also save PDFs?
                        savepdfs=opts.nopdfs,
                        #List of covariance matrix csv files used as analytic likelihood
                        covarianceMatrices=opts.covarianceMatrices,
                        #List of meanVector csv files used, one csv file for each covariance matrix
                        meanVectors=opts.meanVectors,
                        #header file for parameter names in posterior samples
                        header=opts.header,
                        # ascii files (one per IFO) containing  freq - PSD columns
                        psd_files=opts.psdfiles,
                        greedy=not(opts.kdecredibleregions)
                    )

    if opts.archive is not None:
        import subprocess
        subprocess.call(["tar","cvzf",opts.archive,opts.outpath])

    # Send an email, useful for keeping track of dozens of jobs!
    # Will only work if the local host runs a mail daemon
    # that can send mail to the internet
    if opts.email:
        try:
            email_notify(opts.email,opts.outpath)
        except:
            print 'Unable to send notification email'
#
