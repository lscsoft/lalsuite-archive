#!/usr/bin/env python

# Load samples from the posterior PDF of an MCMC or nested sampling code
# and produce sky localisation plots and size estimates.

import sys
import os

from math import ceil,floor

from optparse import OptionParser
from ConfigParser import ConfigParser
from time import strftime

import numpy as np
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from pylal import SimInspiralUtils
from pylal import bayespputils as bppu

import cPickle as pickle

def posMode(pos_samps):
    pos_vals=pos_samps[:,-1]
    max_i=0
    max_pos=pos_vals[0]
    for i in range(len(pos_samps)):
        if pos_vals[i] > max_pos:
            max_pos=pos_vals[i]
            max_i=i
    return max_pos,max_i


def pickle_to_file(obj,fname):
    filed=open(fname,'w')
    pickle.dump(obj,filed)
    filed.close()
#

def oneD_dict_to_file(dict,fname):
    filed=open(fname,'w')
    for key,value in dict.items():
        filed.write("%s %s\n"%(str(key),str(value)) )
#

def loadDataFile(filename):
    print filename
    infile=open(filename,'r')
    formatstr=infile.readline().lstrip()
    formatstr=formatstr.replace('#','')
    header=formatstr.split()

    llines=[]
    import re
    dec=re.compile(r'[^Ee+\d.-]+')
    line_count=0
    for line in infile:
        sline=line.split()
        proceed=True
        if len(sline)<1:
            print 'Ignoring empty line in input file: %s'%(sline)
            proceed=False
        for s in sline:
            if dec.search(s) is not None:
                print 'Warning! Ignoring non-numeric data after the header: %s'%(sline)
                proceed=False
        if proceed:
            llines.append(array(map(float,sline)))
    flines=array(llines)
    for i in range(0,len(header)):
        if header[i].lower().find('log')!=-1 and header[i].lower()!='logl':
            print 'exponentiating %s'%(header[i])
            flines[:,i]=exp(flines[:,i])
            header[i]=header[i].replace('log','')
        if header[i].lower().find('sin')!=-1:
            print 'asining %s'%(header[i])
            flines[:,i]=arcsin(flines[:,i])
            header[i]=header[i].replace('sin','')
        if header[i].lower().find('cos')!=-1:
            print 'acosing %s'%(header[i])
            flines[:,i]=arccos(flines[:,i])
            header[i]=header[i].replace('cos','')
        header[i]=header[i].replace('(','')
        header[i]=header[i].replace(')','')
    print 'Read columns %s'%(str(header))
    return header,flines

#

def getinjpar(paramnames,inj,parnum):
    
    if paramnames[parnum]=='mchirp' or paramnames[parnum]=='mc': return inj.mchirp
    if paramnames[parnum]=='mass1' or paramnames[parnum]=='m1':
        (m1,m2)=bppu.mc2ms(inj.mchirp,inj.eta)
        return m1
    if paramnames[parnum]=='mass2' or paramnames[parnum]=='m2':
        (m1,m2)=bppu.mc2ms(inj.mchirp,inj.eta)
        return m2
    if paramnames[parnum]=='eta': return inj.eta
    if paramnames[parnum]=='time': return inj.get_end()
    if paramnames[parnum]=='phi0': return inj.phi0
    if paramnames[parnum]=='dist' or paramnames[parnum]=='distance': return inj.distance
    if paramnames[parnum]=='RA' or paramnames[parnum]=='long': return inj.longitude
    if paramnames[parnum]=='dec' or paramnames[parnum]=='lat': return inj.latitude
    if paramnames[parnum]=='psi': return inj.polarization
    if paramnames[parnum]=='iota': return inj.inclination
    return None
#


def cbcBayesSkyRes(outdir,data,oneDMenu,twoDGreedyMenu,GreedyRes,confidence_levels,twoDplots,injfile=None,eventnum=None,skyres=None,bayesfactornoise=None,bayesfactorcoherent=None):

    if eventnum is not None and injfile is None:
        print "You specified an event number but no injection file. Ignoring!"

    if data is None:
        print 'You must specify an input data file'
        exit(1)
    #
    if outdir is None:
        print "You must specify an output directory."
        exit(1)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #

    summary_fo=open(os.path.join(outdir,'summary.ini'),'w')

    summary_file=ConfigParser()
    summary_file.add_section('metadata')
    summary_file.set('metadata','group_id','X')
    if eventnum:
        summary_file.set('metadata','event_id',str(eventnum))
    summary_file.add_section('Confidence levels')
    summary_file.set('Confidence levels','confidence levels',str(confidence_levels))

    # Load in the main data
    paramnames, pos=loadDataFile(data[0])

    #Generate any required derived parameters
    if "m1" not in paramnames and "m2" not in paramnames and "mchirp" in paramnames and "eta" in paramnames:
        (m1,m2)=bppu.mc2ms(pos[:,paramnames.index('mchirp')],pos[:,paramnames.index('eta')])

        pos=np.column_stack((pos,m1,m2))
        paramnames.append("m1")
        paramnames.append("m2")
    #
    Nd=len(paramnames)
    print "Number of posterior samples: " + str(size(pos,0))
    # Calculate means
    means = mean(pos,axis=0)
    meanStr=map(str,means)
    out=reduce(lambda a,b:a+'||'+b,meanStr)
    print 'Means:'
    print '||'+out+'||'

    RAdim=paramnames.index('RA')
    decdim=paramnames.index('dec')

    injection=None

    # Select injections using tc +/- 0.1s if it exists or eventnum from the injection file
    if injfile:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])
        if(eventnum is not None):
            if(len(injections)<eventnum):
                print "Error: You asked for event %d, but %s contains only %d injections" %(eventnum,injfile,len(injections))
                sys.exit(1)
            else:
                injection=injections[eventnum]
        else:
            if(len(injections)<1):
                print 'Warning: Cannot find injection with end time %f' %(means[2])
            else:
                injection = itertools.ifilter(lambda a: abs(a.get_end() - means[2]) < 0.1, injections).next()

    #If injection parameter passed load object representation of injection
    #table entries.
    if injection:
        injpoint=map(lambda a: getinjpar(paramnames,injection,a),range(len(paramnames)))
        injvals=map(str,injpoint)
        out=reduce(lambda a,b:a+'||'+b,injvals)
        print 'Injected values:'
        print out

        #Add injection values to output file
        summary_file.add_section('Injection values')

        for parnum in range(len(paramnames)):
            summary_file.set('Injection values',paramnames[parnum],getinjpar(paramnames,injection,parnum))

    #

    twoDGreedyCL={}
    twoDGreedyInj={}

    #If sky resolution parameter has been specified try and create sky map.
    skyreses=None
    if skyres is not None:

        RAvec=array(pos)[:,paramnames.index('RA')]
        decvec=array(pos)[:,paramnames.index('dec')]
        skypos=column_stack([RAvec,decvec])
        injvalues=None
        if injection:
            injvalues=(injpoint[RAdim],injpoint[decdim])

        skyreses,toppoints,skyinjectionconfidence,min_sky_area_containing_injection=bppu.plotSkyMap(skypos,skyres,injvalues,confidence_levels,outdir)
        
        if skyinjectionconfidence:
            twoDGreedyInj['ra_sb,dec_sb']={}
            twoDGreedyInj['ra_sb,dec_sb']['confidence']=min_sky_area_containing_injection
            if min_sky_area_containing_injection:
                twoDGreedyInj['ra_sb,dec_sb']['area']=min_sky_area_containing_injection
            
            
    # Add bayes factor information to summary file
    summary_file.add_section('bayesfactor')
    if bayesfactornoise is not None:
        bfile=open(bayesfactornoise,'r')
        BSN=bfile.read()
        bfile.close()
        summary_file.set('bayesfactor','BSN',BSN)
    if bayesfactorcoherent is not None:
        bfile=open(bayesfactorcoherent,'r')
        BCI=bfile.read()
        bfile.close()
        summary_file.set('bayesfactor','BCI',BCI)

    #Loop over parameter pairs in twoDGreedyMenu and bin the sample pairs
    #using a greedy algorithm . The ranked pixels (toppoints) are used
    #to plot 2D histograms and evaluate Bayesian confidence intervals.

    summary_file.add_section('2D greedy cl')
    summary_file.add_section('2D greedy cl inj')

    ncon=len(confidence_levels)
    pos_array=np.array(pos)

    

    for par1_name,par2_name in twoDGreedyMenu:
        print "Binning %s-%s to determine confidence levels ..."%(par1_name,par2_name)

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

        #Get posterior samples
        try:
            par1_index=paramnames.index(par1_name)
        except ValueError:
            print "No input chain for %s, skipping %s/%s binning."%(par1_name,par1_name,par2_name)
            continue
        try:
            par2_index=paramnames.index(par2_name)
        except ValueError:
            print "No input chain for %s, skipping %s/%s binning."%(par2_name,par1_name,par2_name)
            continue

        pars_name=("%s,%s"%(par1_name,par2_name)).lower()

        par1pos=pos_array[:,par1_index]
        par2pos=pos_array[:,par2_index]

        injection_=None
        if injection:
            par1_injvalue=np.array(injpoint)[par1_index]
            par2_injvalue=np.array(injpoint)[par2_index]
            injection_=(par1_injvalue,par2_injvalue)

        posterior_array=column_stack([par1pos,par2pos])

        toppoints,injectionconfidence,twoDGreedyCL[pars_name],twoDGreedyInj[pars_name]=bppu.greedyBin2(posterior_array,(par1_bin,par2_bin),confidence_levels,par_names=(par1_name,par2_name),injection=injection_)
        
        #Plot 2D histograms of greedily binned points 
        if injection is not None and par1_injvalue is not None and par2_injvalue is not None:
            bppu.plot2Dbins(np.array(toppoints),(par1_bin,par2_bin),outdir,par_names=(par1_name,par2_name),injpoint=[par1_injvalue,par2_injvalue])
        else:
            bppu.plot2Dbins(np.array(toppoints),(par1_bin,par2_bin),outdir,par_names=(par1_name,par2_name))
    #
        summaryString='['
        for frac,area in twoDGreedyCL[pars_name].items():
            summaryString+=str(area)
        summary_file.set('2D greedy cl',pars_name,summaryString+']')
        summary_file.set('2D greedy cl inj',pars_name,str(injectionconfidence))

    if not os.path.exists(os.path.join(outdir,'pickle')):
        os.makedirs(os.path.join(outdir,'pickle'))


    pickle_to_file(twoDGreedyCL,os.path.join(outdir,'pickle','GreedyCL2.pickle'))
    pickle_to_file(twoDGreedyInj,os.path.join(outdir,'pickle','GreedyInj2.pickle'))

    for par in twoDGreedyCL.keys():
        oneD_dict_to_file(twoDGreedyCL[par],os.path.join(outdir,str(par)+'_greedy2.dat'))

    #1D binning
    summary_file.add_section('1D mean')
    summary_file.add_section('1D median')
    summary_file.add_section('1D mode')
    summary_file.add_section('1D contigious cl')
    summary_file.add_section('1D greedy cl')
    summary_file.add_section('1D stacc')
    
    oneDStats={}
    oneDGreedyCL={}
    oneDContCL={}
    oneDGreedyInj={}
    oneDContInj={}

    max_pos,max_i=posMode(pos_array)
    
    #Loop over each parameter and determine the contigious and greedy
    #confidence levels and some statistics.

    for par_name in oneDMenu:
        print "Binning %s to determine confidence levels ..."%par_name
        try:
            par_index=paramnames.index(par_name)
        except ValueError:
            print "No input chain for %s, skipping binning."%par_name
            continue
        try:
            par_bin=GreedyRes[par_name]
        except KeyError:
            print "Bin size is not set for %s, skipping binning."%par_name
            continue

        oneDGreedyCL[par_name]={}
        oneDStats[par_name]={}
        oneDContCL[par_name]={}

        oneDGreedyInj[par_name]={}
        oneDContInj[par_name]={}
        
        par_samps=pos_array[:,par_index]

        summary_file.set('1D mode',par_name,str(par_samps[max_i]))
        summary_file.set("1D mean",par_name,str(np.mean(par_samps)))
        summary_file.set("1D median",par_name,str(np.median(par_samps)))

        oneDStats[par_name]['mode']=par_samps[max_i]
        oneDStats[par_name]['mean']=np.mean(par_samps)
        oneDStats[par_name]['median']=np.median(par_samps)

        par_injvalue_=None
        if injection:
            par_injvalue_=np.array(injpoint)[par_index]

        oneDGreedyCL[par_name],oneDGreedyInj[par_name],toppoints,injectionconfidence = bppu.greedyBin1(par_samps,par_bin,confidence_levels,par_injvalue=par_injvalue_)

        oneDContCL[par_name],oneDContInj[par_name]=bppu.contCL1(par_samps,par_bin,confidence_levels,par_injvalue=par_injvalue_)
    
        #Ilya's standard accuracy statistic
        if injection:
            injvalue=np.array(injpoint)[par_index]
            if injvalue:
                stacc=bppu.stacc_stat(par_samps,injvalue)
                summary_file.set('1D stacc',par_name,str(stacc))
                oneDStats[par_name]['stacc']=stacc
    
    pickle_to_file(oneDGreedyCL,os.path.join(outdir,'pickle','GreedyCL1.pickle'))
    pickle_to_file(oneDContCL,os.path.join(outdir,'pickle','ContCL1.pickle'))
    pickle_to_file(oneDStats,os.path.join(outdir,'pickle','Stats1.pickle'))
    pickle_to_file(oneDContInj,os.path.join(outdir,'pickle','ContInj1.pickle'))
    pickle_to_file(oneDGreedyInj,os.path.join(outdir,'pickle','GreedyInj1.pickle'))
    
    for par in oneDGreedyCL.keys():
        oneD_dict_to_file(oneDGreedyCL[par],os.path.join(outdir,str(par)+'_greedy1.dat'))
    #

    for par in oneDGreedyCL.keys():
        oneD_dict_to_file(oneDContCL[par],os.path.join(outdir,str(par)+'_cont.dat'))
    #

    for par in oneDStats.keys():
        oneD_dict_to_file(oneDStats[par],os.path.join(outdir,str(par)+'_stats.dat'))
    #

    for par in oneDStats.keys():
        oneD_dict_to_file(oneDContInj[par],os.path.join(outdir,str(par)+'_cont_inj.dat'))
    #


    #####Generate 2D kde plots and webpage########
    margdir=os.path.join(outdir,'2D')
    if not os.path.isdir(margdir):
        os.makedirs(margdir)

    twoDKdePaths=[]
    
    for par1,par2 in twoDplots:

        try:
            i=paramnames.index(par1)
        except ValueError:
            print "No input chain for %s, skipping 2D plot of %s-%s."%(par1,par1,par2)
            continue
        try:
            j=paramnames.index(par2)
        except ValueError:
            print "No input chain for %s, skipping 2D plot of %s-%s."%(par2,par1,par2)
            continue

        print 'Generating %s-%s plot'%(paramnames[i],paramnames[j])

        if (size(np.unique(pos[:,i]))<2 or size(np.unique(pos[:,j]))<2):
            continue

        par_injvalues_=None
        if injection and reduce (lambda a,b: a and b, map(lambda idx: getinjpar(paramnames,injection,idx)>min(pos[:,idx]) and getinjpar(paramnames,injection,idx)<max(pos[:,idx]),[i,j])) :
            if getinjpar(paramnames,injection,i) is not None and getinjpar(paramnames,injection,j) is not None:
                par_injvalues_=( getinjpar(paramnames,injection,i) , getinjpar(paramnames,injection,j) )

        myfig=bppu.plot2Dkernel(pos[:,i],pos[:,j],50,50,par_names=(par1,par2),par_injvalues=par_injvalues_)

        twoDKdePath=os.path.join(margdir,paramnames[i]+'-'+paramnames[j]+'_2Dkernel.png')
        twoDKdePaths.append(twoDKdePath)

        myfig.savefig(twoDKdePath)

    htmlfile=open(os.path.join(outdir,'posplots.html'),'w')
    htmlfile.write('<HTML><HEAD><TITLE>Posterior PDFs</TITLE></HEAD><BODY><h3>'+str(means[2])+' Posterior PDFs</h3>')
    if(skyres is not None):
        htmlfile.write('<table border=1><tr><td>Confidence region<td>size (sq. deg)</tr>')
        for (frac,skysize) in skyreses:
            htmlfile.write('<tr><td>%f<td>%f</tr>'%(frac,skysize))
        htmlfile.write('</table>')
    if bayesfactornoise is not None:
        htmlfile.write('<p>log Bayes factor ( coherent vs gaussian noise) = %s, Bayes factor=%f</p>'%(BSN,exp(float(BSN))))
    if bayesfactorcoherent is not None:
        htmlfile.write('<p>log Bayes factor ( coherent vs incoherent OR noise ) = %s, Bayes factor=%f</p>'%(BCI,exp(float(BCI))))
    htmlfile.write('Produced from '+str(size(pos,0))+' posterior samples.<br>')
    htmlfile.write('Samples read from %s<br>'%(data[0]))
    htmlfile.write('<h4>Mean parameter estimates</h4>')
    htmlfile.write('<table border=1><tr>')
    paramline=reduce(lambda a,b:a+'<td>'+b,paramnames)
    htmlfile.write('<td>'+paramline+'<td>logLmax</tr><tr>')
    meanline=reduce(lambda a,b:a+'<td>'+b,meanStr)
    htmlfile.write('<td>'+meanline+'</tr>')
    if injection:
        htmlfile.write('<tr><th colspan='+str(len(paramnames))+'>Injected values</tr>')
        injline=reduce(lambda a,b:a+'<td>'+b,injvals)
        htmlfile.write('<td>'+injline+'<td></tr>')
    htmlfile.write('</table>')
    if injection:
        if skyinjectionconfidence:
            htmlfile.write('<p>Injection found at confidence interval %f in sky location</p>'%(skyinjectionconfidence))
        else:
            htmlfile.write('<p>Injection not found in posterior bins in sky location!</p>')
    htmlfile.write('<h5>2D Marginal PDFs</h5><br>')
    htmlfile.write('<table border=1 width=100%><tr>')
    #htmlfile.write('<td width=30%><img width=100% src="m1m2.png"></td>')
    #htmlfile.write('<td width=30%><img width=100% src="RAdec.png"></td>')
    #htmlfile.write('<td width=30%><img width=100% src="Meta.png"></td>')
    #htmlfile.write('</tr><tr><td width=30%><img width=100% src="2D/Mchirp (Msun)-geocenter time ISCO_2Dkernel.png"</td>')
    if skyres is not None:
        htmlfile.write('<td width=30%><img width=100% src="skymap.png"></td>')
    else:
        htmlfile.write('<td width=30%><img width=100% src="m1dist.png:></td>')
    #htmlfile.write('<td width=30%><img width=100% src="m2dist.png"></td>')

    row_switch=1
    for par1,par2 in twoDplots:

        if row_switch==3:
            row_switch=0

        plot_path=None
        if os.path.isfile(os.path.join(outdir,'2D',par1+'-'+par2+'_2Dkernel.png')):
            plot_path='2D/'+par1+'-'+par2+'_2Dkernel.png'

        elif os.path.isfile(os.path.join(outdir,'2D',par2+'-'+par1+'_2Dkernel.png')):
            plot_path='2D/'+par2+'-'+par1+'_2Dkernel.png'

        if plot_path:
            if row_switch==0:
                htmlfile.write('<tr>')

            htmlfile.write('<td width=30%><img width=100% src="'+plot_path+'"></td>')
            if row_switch==2:
                htmlfile.write('</tr>')
            row_switch+=1
    #
    if row_switch==2:
        htmlfile.write('<td></td></tr>')
    elif row_switch==1:
        htmlfile.write('<td></td><td></td></tr>')


    htmlfile.write('</table>')
    htmlfile.write('<br><a href="2D/">All 2D Marginal PDFs</a><hr><h5>1D marginal posterior PDFs</h5><br>')

    summary_file.add_section('1D ranking kde')
    summary_file.add_section('1D ranking bins')
    
    oneDplotPaths=[]
    for param in oneDMenu:

        try:
            par_index=paramnames.index(param)
            i=par_index
        except ValueError:
            print "No input chain for %s, skipping 1D plot."%param
            continue

        pos_samps=pos[:,i]

        injpar_=None
        if injection:
            injpar_=getinjpar(paramnames,injection,i)

        print "Generating 1D plot for %s."%param
        rbins,plotFig=bppu.plot1DPDF(pos_samps,param,injpar=injpar_)
        figname=param+'.png'
        oneDplotPath=os.path.join(outdir,figname)
        
        plotFig.savefig(os.path.join(outdir,param+'.png'))
        if rbins:
            print "r of injected value of %s (bins) = %f"%(param, rbins)

        ##Produce plot of raw samples
        myfig=plt.figure(figsize=(4,3.5),dpi=80)
        plt.plot(pos_samps,'.',figure=myfig)
        if injpar_:
            if min(pos_samps)<injpar_ and max(pos_samps)>injpar_:
                plt.plot([0,len(pos_samps)],[injpar_,injpar_],'r-.')
        myfig.savefig(os.path.join(outdir,figname.replace('.png','_samps.png')))
    
        #summary_file.set('1D ranking kde',param,rkde)
        summary_file.set('1D ranking bins',param,rbins)

        oneDplotPaths.append(figname)
    htmlfile.write('<table><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th>')
    for plotPath in oneDplotPaths:
        htmlfile.write('<tr><td><img src="'+plotPath+'"></td><td><img src="'+plotPath.replace('.png','_samps.png')+'"></td>')
    htmlfile.write('</table>')
    htmlfile.write('<hr><br />Produced using cbcBayesSkyRes.py at '+strftime("%Y-%m-%d %H:%M:%S"))
    htmlfile.write('</BODY></HTML>')
    htmlfile.close()


    # Save posterior samples too...

    posfilename=os.path.join(outdir,'posterior_samples.dat')
    posfile=open(posfilename,'w')
    for row in pos:
        for i in row:
            posfile.write('%f\t'%(i))
        posfile.write('\n')
    #

    #Close files
    posfile.close()
    summary_file.write(summary_fo)

if __name__=='__main__':

    parser=OptionParser()
    parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
    parser.add_option("-d","--data",dest="data",action="append",help="datafile")
    parser.add_option("-i","--inj",dest="injfile",help="SimInsipral injection file",metavar="INJ.XML",default=None)
    parser.add_option("--skyres",dest="skyres",help="Sky resolution to use to calculate sky box size",default=None)
    parser.add_option("--eventnum",dest="eventnum",action="store",default=None,help="event number in SimInspiral file of this signal",type="int",metavar="NUM")
    parser.add_option("--bsn",action="store",default=None,help="Optional file containing the bayes factor signal against noise",type="string")
    parser.add_option("--bci",action="store",default=None,help="Optional file containing the bayes factor coherent against incoherent models",type="string")

    (opts,args)=parser.parse_args()

    #List of parameters to plot/bin . Need to match (converted) column names.
    oneDMenu=['mtotal','m1','m2','mchirp','mc','distance','distMPC','dist','iota','psi','eta','RA','dec','a1','a2','phi1','theta1','phi2','theta2']
    #List of parameter pairs to bin . Need to match (converted) column names.
    twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['RA','dec']]
    #Bin size/resolution for binning. Need to match (converted) column names.
    greedyRes={'mc':0.025,'m1':0.1,'m2':0.1,'mtotal':0.1,'eta':0.001,'iota':0.01,'time':1e-4,'distance':1.0,'dist':1.0,'mchirp':0.025,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'RA':0.05,'dec':0.05}
    #Confidence levels
    confidenceLevels=[0.67,0.9,0.95,0.99]
    #2D plots list
    twoDplots=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['m1','dist'],['m2','dist'],['psi','iota'],['psi','distance'],['psi','dist'],['psi','phi0']]

    
    cbcBayesSkyRes(opts.outpath,opts.data,oneDMenu,twoDGreedyMenu,greedyRes,confidenceLevels,twoDplots,injfile=opts.injfile,eventnum=opts.eventnum,skyres=opts.skyres,bayesfactornoise=opts.bsn,bayesfactorcoherent=opts.bci)
#
