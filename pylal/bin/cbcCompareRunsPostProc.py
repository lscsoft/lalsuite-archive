#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcCompareRunsPostProc.py
#
#       Copyright 2011
#       Salvatore Vitale <salvatore.vitale@ligo.org>
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

#standard library imports
import sys
import os

from math import ceil,floor
import cPickle as pickle

from time import strftime

#related third party imports

import numpy as np
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack,cov,unique,hsplit,correlate,log,dot,power
from scipy import *
import matplotlib 
matplotlib.use("Agg")
from pylab import *
try:
    from xml.etree.cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser,XML
except ImportError:
    #Python < 2.5
    from cElementTree import Element, SubElement, ElementTree, Comment, tostring, XMLParser,XML

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

__author__="Salvatore Vitale <salvatore.vitale@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

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
        
def remove_ifo_from_string(name):
    ifos=['H1','L1','V1','Network']
    for ifo in ifos:
        if name.find(ifo)!=-1:
            return [name[0:len(ifo)],str(name[len(ifo)+1:])]

def linear_space( start, end, count ):
    delta = (end-start) / float(count)
    return [start,] + \
    map( lambda x:delta*x + start, range( 1, count ) ) + [end, ]
     
def linkImage(imageurl,width,height):
        return '<a href="'+ str(imageurl)+'" target="_blank"><img src="'+str(imageurl)+'" width="'+str(width) +'" height="'+str(height)+'" /></a>'
def link(url,text):
        return '<a href="'+ str(url)+'">'+str(text)+'</a>'
def checkDir(dir):
        if not os.access(dir,os.F_OK):
                try:
                        print 'Creating %s\n'%(dir)
                        os.makedirs(dir)
                except:
                        print 'Error: Unable to create directory %s'%(dir)
                        sys.exit(1)
        if not os.access(dir,os.W_OK):
                print 'Unable to write to %s'%(dir)
                sys.exit(1)
        return True
        
def std_dev(array,mean):
    return sqrt(sum((x-mean) ** 2 for x in array) / len(array))
def mean(array):
    return np.mean(array)
def skewness(array,mean,std_dev):
    return (sum((x-mean) ** 3 for x in array) / (len(array)*std_dev **3) )
def kurtosis(array,mean,std_dev):
    return (sum((x-mean) ** 4 for x in array) / (len(array)*std_dev **4) )
def read_snr(snrs,run,time,IFOs):
    snr_values=""
    snr_header=[]
    net_snr2=0.0
    for IFO in IFOs:
        path_to_file=os.path.join(snrs[run],'snr_'+IFO+'_'+str(time)+'.0.dat')
        path_to_multifile=os.path.join(snrs[run],'snr_'+''.join(IFO for IFO in IFOs)+'_'+str(time)+'.0.dat')
        if os.path.isfile(path_to_file):
            snr_file=open(path_to_file,'r')
            line=snr_file.readline()[:-1]
            snr_header.append('SNR_'+remove_ifo_from_string(line)[0])
            snr_values+=str(remove_ifo_from_string(line)[1])+' '
            net_snr2+=remove_ifo_from_string(line)[1]**2
            snr_file.close()
    if snr_header!=[]:
        snr_header.append('SNR_Network')
        snr_values+=str(sqrt(net_snr2))
    if os.path.isfile(path_to_multifile) and snr_header==[]:
        snr_file=open(path_to_multifile,'r')
        for line in snr_file.readlines()[:-1]:
            snr_header.append('SNR_'+remove_ifo_from_string(line)[0])
            snr_values +=remove_ifo_from_string(line)[1][:-1]
        snr_file.close()
    return [snr_values,snr_header]
def histN(mat,N):
    Nd=size(N)
    histo=zeros(N)
    scale=array(map(lambda a,b:a/b,map(lambda a,b:(1*a)-b,map(max,mat),map(min,mat)),N)) ### takes mat and calculates (max_omega - min_omega)/N, scale is a Nd-dims array
    axes=array(map(lambda a,N:linspace(min(a),max(a),N),mat,N)) ### takes man and N and returns the axes
    print shape(axes)
    bins=floor(map(lambda a,b:a/b , map(lambda a,b:a-b, mat, map(min,mat) ),scale*1.01))
    hbins=reshape(map(int,bins.flat),bins.shape)
    print str(bins.flat)
    for co in transpose(hbins):
        t=tuple(co)
        histo[t[::-1]]=histo[t[::-1]]+1
    return (axes,histo)
def MakeErrorPlots(time,outdir,in_data_path,run,f_0,f_up,IFOs):
    run=str(run)
    path_plots=os.path.join(outdir,run,'ErrorPlots')
    checkDir(path_plots)
    data={}
    for IFO in IFOs:
        path_to_data=os.path.join(in_data_path[int(run)-1],'calerr_'+IFO+'_'+str(time)+'.0.dat')
        data[IFO]=np.loadtxt(path_to_data)
    a=0
    for i in range(len(data[IFOs[0]][:,0])):
        if fabs(data[IFOs[0]][i,0]-f_0)<0.0001:
            a=i
            continue 
    if a==0:
        print "Could not fix f_low. Exiting...\n"
        exit(-1)
    b=0
    for i in range(len(data[IFOs[0]][a:,0])+a):
        if fabs(data[IFOs[0]][i,0]-f_up)<0.0001:
            b=i
            continue
    if b==0:
        print "Could not fix f_up. Exiting\n"
        exit(-1)
    
    myfig=figure(1,figsize=(10,8),dpi=80)
    for (IFO,color) in zip(IFOs,['r','b','k']):
        plot(data[IFO][a:b,0],data[IFO][a:b,1],color,label=IFO)
    xlabel('f[Hz]')
    ylabel('Amp_cal/Amp_uncal')
    grid()
    legend()
    myfig.savefig(os.path.join(path_plots,'amp_'+str(time)+'.png'))
    myfig.clear()
    phase={}
    phase_normalized={}

    def normalizer(phase,phase_normalized):
        for pha in phase:
            if pha<-1.5*pi:
                phase_normalized.append(pha+2*pi)
            elif pha>pi:
                phase_normalized.append(pha-2*pi)
            else:
                phase_normalized.append(pha)

    for IFO in IFOs:
        phase[IFO]=data[IFO][:,2]
        phase_normalized[IFO]=[]
        normalizer(phase[IFO],phase_normalized[IFO])

    myfig=figure(1,figsize=(10,8),dpi=80)
    for (IFO,color) in zip(IFOs,['r','b','k']):
        plot(data[IFO][a:b,0],phase_normalized[IFO][a:b],color,label=IFO)
    xlabel('f[Hz]')
    ylabel('Pha_cal-Pha_uncal [Rads]')
    legend()
    grid()
    myfig.savefig(os.path.join(path_plots,'pha_'+str(time)+'.png'))
    myfig.clear()


def RunsCompare(outdir,inputs,inj,raw_events,IFOs,snrs=None,calerr=None,path_to_result_pages=None):
    
    from pylal import SimInspiralUtils
    
    checkDir(outdir)
    number_of_inits=len(inputs)
    injTable=SimInspiralUtils.ReadSimInspiralFromFiles([inj])
    flow=20.0
    fup=500.0
    time_event=None
    if inj and raw_events:
        time_event={}
        events=[]
        times=[]
        raw_events=raw_events.replace('[','').replace(']','').split(',')
        if 'all' is raw_events or 'all' in raw_events:
            events=range(len(injTable))

        else:
            for raw_event in raw_events:
                if ':' in raw_event:
                    limits=raw_event.split(':')
                    if len(limits) != 2:
                        print "Error: in event config option; ':' must separate two numbers."
                        exit(0)
                    low=int(limits[0])
                    high=int(limits[1])+1   # The plus one is to include the rigthmost extremum
                    if low>high:
                        events.extend(range(int(high),int(low)))
                    elif high>low:
                        events.extend(range(int(low),int(high)))
                else:
                    events.append(int(raw_event))

        for event in events:
            times.append(injTable[event].get_end())
            time_event[times[-1]]=event

        starttime=min(times)-1.0
        endtime=max(times)+1.0
    else:
        print "Error, you must give --events and --inj option \n"

    BSN=[path for path in inputs]
    Combine=[path for path in inputs]
    snrs=[path for path in snrs]

    ## Remove the times for which posterior file is not present in either of the init ## TBD I only need to compare couple of ctrl-cali, not all of them
    for run in range(len(Combine)):
        for time in times:
            path_to_file=os.path.join(Combine[run],'posterior_samples_'+str(time)+'.000')
            if not os.path.isfile(path_to_file):
                times.remove(time)
                continue

 
    ## prepare files with means and other useful data ###
    for run in range(len(Combine)):
        if int(run)==0:
            summary=open(os.path.join(outdir,'summary_ctrl.dat'),'w')
        else:
            summary=open(os.path.join(outdir,'summary_'+str(run)+'.dat'),'w')
        header=open(os.path.join(outdir,'headers_'+str(run)+'.dat'),'w')
        header_l=[]
        header_l.append('injTime \t')
        for time in times:
            path_to_file=os.path.join(Combine[run],'posterior_samples_'+str(time)+'.000')
	    posterior_file=open(path_to_file,'r')
            peparser=bppu.PEOutputParser('common')
            commonResultsObj=peparser.parse(posterior_file)
            pos = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injTable[times.index(time)])
            posterior_file.close()
            parameters=pos.names
            
            parameters.remove('logl')
            summary.write(str(time)+'\t')
            for parameter in parameters:
                summary.write(repr(pos[parameter].mean) + '\t'+ repr(pos[parameter].stdev) +'\t')
                if time==times[0]:
                    header_l.append('mean_'+parameter)
                    header_l.append('stdev_'+parameter)
            if BSN is not None:
                path_to_file=os.path.join(BSN[run],'bayesfactor_'+str(time)+'.000.txt')
                bfile=open(path_to_file,'r')
                bsn=bfile.readline()[:-1] ## remove the newline tag
                bfile.close()
                summary.write(str(bsn)+'\t')
                if time==times[0]:
                    header_l.append('BSN')
            if snrs is not None:
                val,hea=read_snr(snrs,run,time,IFOs)
                summary.write(str(val)+'\t')
                if time==times[0]:
                    for he in hea:
                        header_l.append(he)
            summary.write('\n')
        header.write('\t'.join(str(n) for n in header_l)+'\n')
        summary.close()
        header.close()
    
    path_uncal=os.path.join(outdir,'summary_ctrl.dat')    
    for run in range(1,len(Combine)):
        run=str(run)
        path_cal=os.path.join(outdir,'summary_'+run+'.dat')
        MakePlots(outdir,path_cal,path_uncal,run,parameters)
        if snrs is not None:
            MakeSNRPlots(outdir,snrs,path_cal,path_uncal,run,parameters,header_l,IFOs)
        if BSN is not None and snrs is not None:
            MakeBSNPlots(outdir,path_cal,path_uncal,run,header_l)
        if calerr is not None:
            MakeErrorPlots(times[0],outdir,calerr,run,flow,fup,IFOs)
        WritePlotPage(outdir,run,parameters,times[0])
        
        WriteSummaryPage(outdir,run,path_to_result_pages[int(run)],path_to_result_pages[0],header_l,times,IFOs)



def MakePlots(outdir,path_cal,path_uncal,run,parameters):
    nbins=20
    data_cal=np.loadtxt(path_cal)
    data_uncal=np.loadtxt(path_uncal)
    path_plots=os.path.join(outdir,run,'ParametersPlots')
    checkDir(path_plots)
    
    for parameter in parameters:
        x=parameters.index(parameter)*2+1
        x_delta=(data_cal[:,x]-data_uncal[:,x])
        x_points=(map(lambda t:(min(x_delta) + (t/max(x_delta))*(max(x_delta)-min(x_delta))),x_delta)).sort
        x_points2=linear_space(min(x_delta),max(x_delta),len(x_delta))
        myfig=figure(2,figsize=(10,10),dpi=80)
        plot(x_delta,data_uncal[:,x+1],'.r',label='stdev')
        axvline(x=0, ymin=0, ymax=1,linewidth=2, color='b')
        plot(x_points2,fabs(x_points2)*2,'-k',label='$0.5 \sigma$')
        plot(x_points2,fabs(x_points2),'-y',label='$\sigma$')
        plot(x_points2,fabs(x_points2)/2,'-c',label='2$\sigma$')
        xlabel("delta_"+str(parameter))
        ylabel("sigma_"+str(parameter))
        grid()
        #legend(loc="2")
        myfig.savefig(os.path.join(path_plots,"delta_sigma_"+str(parameter)+'.png'))
        myfig.clear()
        mean_delta_omega=mean(x_delta)
        std_delta_omega=std_dev(x_delta,mean_delta_omega)
        skewness_delta_omega= skewness(x_delta, mean_delta_omega, std_delta_omega)
        kurtosis_delta_omega= kurtosis(x_delta, mean_delta_omega, std_delta_omega)
        #print "Mean : %e\n" % mean_delta_omega
        #print "Standard Deviation %e\n" % std_delta_omega
        #print "Skewness %e\n" % skewness_delta_omega
        #print "Kurtosis %e\n" % kurtosis_delta_omega
        
        bins=linear_space(x_delta.min(),x_delta.max(),nbins)
        myfig=figure(figsize=(4,3.5),dpi=80)
        xlabel(u"\u0394"+ parameter)
        hist(x_delta, bins=bins, normed="true")
        axvline(x=0, ymin=0, ymax=1,linewidth=2, color='r')
        grid()
        myfig.savefig(os.path.join(path_plots,'delta_'+parameter+'.png'))
        myfig.clear()
        ##### This part calculates effect size
        effect_size =(data_cal[:,x]-data_uncal[:,x])/data_uncal[:,x+1]  
        mean_effect_size=mean(effect_size)
        std_effect_size=std_dev(effect_size,mean_effect_size)
        skewness_effect_size= skewness(effect_size, mean_effect_size, std_effect_size)
        kurtosis_effect_size= kurtosis(effect_size, mean_effect_size, std_effect_size)
        #print "With effect size correction\n"
        #print "Mean : %e\n" % mean_effect_size
        #print "Standard Deviation %e\n" % std_effect_size
        #print "Skewness %e\n" % skewness_effect_size
        #print "Kurtosis %e\n" % kurtosis_effect_size
        bins=linear_space(effect_size.min(),effect_size.max(),nbins)
        myfig2=figure(figsize=(4,3.5),dpi=80)
        hist(effect_size, bins=bins, normed="true",color='r',fill=False, hatch='//', linewidth='2')
        xlabel('effect_'+parameter)
        axvline(x=0, ymin=0, ymax=1,linewidth=2, color='b')
        grid()
        legend()
        myfig2.savefig(os.path.join(path_plots,'effect_'+parameter+'.png'))
        myfig2.clear()
        
def MakeSNRPlots(outdir,snrs,path_cal,path_uncal,run,parameters,header_l,IFOs):
    
    data_cal=np.loadtxt(path_cal)
    data_ctrl=np.loadtxt(path_uncal)
    path_plots=os.path.join(outdir,run,'SNRPlots')
    checkDir(path_plots)
    network_snrs=data_cal[:,header_l.index('SNR_Network')]
    network_snrs_ctrl=data_ctrl[:,header_l.index('SNR_Network')]
    for parameter in parameters:
        i=parameters.index(parameter)*2+1
        y_delta=(data_cal[:,i]-data_ctrl[:,i])
        myfig=figure(2,figsize=(10,10),dpi=80)
        ax=myfig.add_subplot(111)
        ax.plot(network_snrs,y_delta,'ro',label='DeltavsSNR')
        ax.set_xlabel('Network SNR cal')
        ax.set_ylabel('delta_%s'%parameter)
        locs, labels = (ax.get_xticks(),ax.get_xticklabels)
        ax2=ax.twiny()
        grid()
        ax2.set_xticklabels(np.linspace(min(network_snrs_ctrl),max(network_snrs_ctrl),len(locs)))
        ax.set_xlabel('Network SNR ctrl')
        #legend()
        myfig.savefig(os.path.join(path_plots,'SNR_vs_'+parameter+'.png'))
        myfig.clear()
    return

def MakeBSNPlots(outdir,path_cal,path_uncal,run,header_l):
    data_cal=np.loadtxt(path_cal)
    data_uncal=np.loadtxt(path_uncal)
    path_plots=os.path.join(outdir,run,'BSNPlots')
    checkDir(path_plots)
    network_snrs=data_cal[:,header_l.index('SNR_Network')]
    bsns=data_cal[:,header_l.index('BSN')]
    myfig=figure(2,figsize=(10,10),dpi=80)
    plot(network_snrs,bsns,'ro',label='BSNvsSNR')
    title("$\mathrm{log\,B}_{coherent,gaussian}$ vs SNR")
    myfig.savefig(os.path.join(path_plots,'BSN_vs_SNR.png'))
    grid()
    legend()
    myfig.clear()

def WritePlotPage(outdir,run,parameters,time):
    run=str(run)
    time=str(time)
    wd=300
    hg=250
    abs_page_path=os.path.join(outdir,run)
    page_path="./"
    path_plots=os.path.join(page_path,'ParametersPlots')
    error_path_plots=os.path.join(page_path,'ErrorPlots')
    snr_plots=os.path.join(page_path,'SNRPlots')
    bsn_plots=os.path.join(page_path,'BSNPlots')

    html=bppu.htmlPage('CalibrationErrors',css=bppu.__default_css_string)
    html_err=html.add_section('Errors fit')
    html_err_st='<table><tr>'
    for plot in ['amp_','pha_']:
        html_err_st+='<td>'
        if os.path.isfile(os.path.join(abs_page_path,'ErrorPlots',plot +time+'.png')):
            html_err_st+=linkImage(os.path.join(error_path_plots,plot +time+'.png'),1.5*wd,1.5*hg)
        else:
            html_err_st+='<p> No calibration error curves found in ' + error_path_plots +'</p>'
        html_err_st+='</td>'
    html_err_st+='</tr></table>'
    html_err.write(html_err_st)

    html_plots=html.add_section('Summary plots')
    html_plots.write(link(os.path.join(page_path,'summary.html'),"Go to the summary table"))
    html_plots_st='<table>'
    for parameter in parameters:
        html_plots_st+='<tr>'
        for plot in ['delta_','effect_','delta_sigma_']:
            html_plots_st+='<td>'
            html_plots_st+=linkImage(os.path.join(path_plots,plot +parameter+'.png'),wd,hg)
            html_plots_st+='</td>'
        for plot in ['SNR_vs_']:
            html_plots_st+='<td>'
            html_plots_st+=linkImage(os.path.join(snr_plots,plot +parameter+'.png'),wd,hg)
            html_plots_st+='</td>'
        html_plots_st+='</tr>'
    html_plots_st+='<tr><td colspan="4">'
    html_plots_st+=linkImage(os.path.join(bsn_plots,'BSN_vs_SNR.png'),2*wd,2*hg)
    html_plots_st+='</td>'
    html_plots_st+='</tr>'
    html_plots_st+='</table>'
    html_plots.write(html_plots_st) 
    #Save results page
    plotpage=open(os.path.join(abs_page_path,'posposplots.html'),'w')
    plotpage.write(str(html))
    return

def go_home(path):
    current=os.getcwd()
    upo=''
    os.chdir(path)
    while os.getcwd()!=os.environ['HOME']:
        os.chdir(os.path.join(os.getcwd(),'../'))
        upo+='../'
    os.chdir(current)
    return upo

def WriteSummaryPage(outdir,run,path_to_result_pages,path_to_ctrl_result_pages,header_l,times,IFOs):
    
    run=str(run)
    abs_page_path=os.path.join(outdir,run)
    page_path="./"
    path_to_result_pages_from_LSC=path_to_result_pages[path_to_result_pages.find('LSC'):]
    path_to_ctrl_result_pages_from_LSC=path_to_ctrl_result_pages[path_to_ctrl_result_pages.find('LSC'):]
    snr_index={}
    ifos=['SNR_'+ifo for ifo in IFOs]
    ifos.append('SNR_Network')

    for ifo in ifos:
        snr_index[ifo]=header_l.index(ifo)  

    bsn_index=header_l.index('BSN')
    ctrl_data=np.loadtxt(os.path.join(outdir,'summary_ctrl.dat'))
    cal_data=np.loadtxt(os.path.join(outdir,'summary_'+run+'.dat'))

    html=bppu.htmlPage('SummaryPage',css=bppu.__default_css_string)
    html_table=html.add_section('Links to postprocessing pages')
    html_table.write(link(os.path.join(page_path,'posposplots.html'),"Go back to the plots page"))
    html_table_st='<table>'
    html_table_st+='<tr><th align="center" colspan="6"> Control Runs </th><th colspan="6" align="center"> Calibration Runs </th></tr>'
    html_table_st+='<tr>'
    html_table_st+=2*('<th> TriggerTime</th><th> BSN </th><th>'+ifos[0]+'</th><th>'+ifos[1]+'</th><th>'+ifos[2]+'</th><th>'+ifos[3]+'</th>')
    html_table_st+='</tr>'

    for time in times:
        time=str(time)
        html_table_st+='<tr>'
        ctrl_page=os.path.join(go_home(outdir),path_to_ctrl_result_pages_from_LSC,time,'posplots.html')
        cal_page=os.path.join(go_home(outdir),path_to_result_pages_from_LSC,time,'posplots.html')
        html_table_st+='<td>'+link(ctrl_page,time)+'</td>'
        html_table_st+='<td>'+'%4.2f'%(ctrl_data[times.index(time),bsn_index])+'</td>'
        for IFO in ifos:
            html_table_st+='<td>'+'%4.2f'%(ctrl_data[times.index(time),snr_index[IFO]])+'</td>'
        html_table_st+='<td>'+link(cal_page,time)+'</td>'
        html_table_st+='<td>'+'%4.2f'%(cal_data[times.index(time),bsn_index])+'</td>'
        for IFO in ifos:
            html_table_st+='<td>'+'%4.2f'%(cal_data[times.index(time),snr_index[IFO]])+'</td>'
        html_table_st+='</tr>'
    html_table_st+='</table>'
    html_table.write(html_table_st)
    #Save results page
    plotpage=open(os.path.join(abs_page_path,'summary.html'),'w')
    plotpage.write(str(html))

def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
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
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)        

if __name__=='__main__':
    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-o","--outpath", dest="outpath",type="string",help="make page and plots in DIR", metavar="DIR")
    parser.add_option("-d","--data", dest="indata",action="callback", callback=vararg_callback,help="The folders containing the posteriors and the BSN files of the runs", metavar="PathToPosterior1 PathToPosterior2 etc")
    parser.add_option("-s","--snr", dest="snr",action="callback", callback=vararg_callback,help="The folders containing the snrs of the runs",metavar="pathToSnr1 pathToSnr2 etc")
    parser.add_option("-r","--result_pages_path",default=None,dest="rp",action="callback",callback=vararg_callback,help="Paths to the folder containing the postplots pages (this folder must contain the timebins folders)",metavar="r")
    parser.add_option("-i","--inj",dest="inj",action="store",type="string",default=None,help="Injection xml table",metavar="injection.xml")
    parser.add_option("-e","--calerr",dest="calerr",action="callback",callback=vararg_callback,default=None,help="path to calibration errors path",metavar="/pathToCalerr1 /pathToCalerr2 etc")
    parser.add_option("-E","--events",dest="raw_events",action="store",type="string",default=None,metavar="\[0:50\]")
    parser.add_option("-I","--IFOS",dest="IFOs",action="callback", callback=vararg_callback,help="The IFOs used in the analysis", metavar="H1 L1 V1")
    (opts,args)=parser.parse_args()
    #if opts.num_of_init==1 and opts.uncal_path==None:
    #        print "Error, if -n is 1 it means that only jobs with calibration errors are running, then you must provide the path to the posterior of the (already run) corresponding jobs without calibration errors using the option -u /pathToUncalPosteriors \n"
    RunsCompare(opts.outpath,opts.indata,opts.inj,opts.raw_events,opts.IFOs,snrs=opts.snr,calerr=opts.calerr,path_to_result_pages=opts.rp)
