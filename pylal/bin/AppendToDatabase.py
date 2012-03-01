#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       AppendToDatabase.py
#
#       Copyright 2012
#      Salvatore Vitale <salvatore.vitale@ligo.org>
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
import getpass

import numpy as np
from numpy import log,loadtxt,vstack,array,exp,size,argsort
from numpy.random import rand

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version
from optparse import OptionParser
#from combine_evidence import combine_evidence

__author__="Salvatore  Vitale <salvatore.vitale@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date


### Functions borrowed from combine_evidence.py, as I did not find how to import them from here ###
def loaddata(datalist):
    out = list(map(np.loadtxt,datalist))
    Bfiles = list(map(getBfile,datalist))
    return out,Bfiles
def getBfile(datname):
    Bfile=datname+'_B.txt'
    print 'Looking for '+Bfile
    if os.access(Bfile,os.R_OK):
        outstat = np.loadtxt(Bfile)
        return outstat
    else:
        return None
def weightsamp(d,Nlive):
    #creates an array of vectors containing weights for every sample in every segment
    total_weight=[]
    for outfile in d:
        N=len(outfile[:,0])
        N_weighted = N - Nlive
        segment_weight=[]
        for i in range(1,N_weighted + 1):
            logw = -(i)/Nlive
            segment_weight.append(logw)
        for i in range(N_weighted + 1,N + 1):
            logw = -N_weighted / Nlive
            segment_weight.append(logw)
        total_weight += segment_weight
    return total_weight

def nest2pos(samps,weights):
    randoms=rand(size(samps,0))
    wt=weights+samps[:,-1]
    maxwt=max(wt)
    #posidx=find(wt>maxwt+log(randoms))
    posidx=[i for i in range(0,size(weights)) if wt[i]>maxwt+log(randoms[i]) ]
    pos=samps[posidx,:]
    return pos
    
def prodnoise(B):
    """
    Calculates sum (logZnoise[i] for i!=j) for each j to get logZ to add to each of the input files
    """
    N=len(B)
    totalnoise=[]
    for i in range(0,N):
        tn=0
    for j in range(0,N):
        if i!=j:
            tn+=B[j,2]
        totalnoise.append(tn)
    return totalnoise
    
def combine_evidence(data,xflag,Nlive):

    nfiles=len(data)

    #load in seperate data files#
    #datalist = makelist(path)
    (d,Bfiles) = loaddata(data)
    Barray = reduce(lambda x,y: vstack([x,y]), Bfiles)

    #Calculate total Bayes factor#
    if len(data)>1:
        ZnoiseTotal=sum(Barray[:,2])-log(len(data))
        totalBayes= reduce(logadd,Barray[:,0])
        totalBayes= float(totalBayes) - log(len(data)) #divide by 60 because we used 60 priors
    else:
        totalBayes=Barray[0]
        ZnoiseTotal=Barray[2]
    print "Total Bayes Factor= %f" %totalBayes

    #Scale likelihoods for entire sample#
    #Make list of sum(noise evidence), excepting current noise evidence)
    if len(data)>1:
        totalnoise = prodnoise(Barray)
    else: totalnoise=array([0])
    # Add logZnoise for other files to likelihoods for each sample
    if not None in Bfiles:
        for (outfile,noise) in zip(d,totalnoise):
            outfile[:,-1]+=noise

    #Remapping Parameters#
    #for outfile in d:
        #outfile[:,0]=exp(outfile[:,0])
        #outfile[:,4]=exp(outfile[:,4])
        #if xflag:
            #outfile[:,8]=x2iota(outfile[:,8])

    #Posterior Samples
    weights=weightsamp(d,Nlive)
    d_all = reduce(lambda x,y: vstack([x,y]), d)
    pos=nest2pos(d_all,weights)

    d_idx=argsort(d_all[:,-1])
    d_all=d_all[d_idx,:]

    return pos,d_all,totalBayes,ZnoiseTotal

###                                         ###

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

def my_sort(a,b):
    if a=='GR':
        return -1
    elif b=='GR':
        return 1
    elif len(a)<len(b):
        return -1
    elif len(a)>len(b):
        return 1
    elif str(a)<str(b):
        return -1
    elif str(a)>str(b):
        return 1
    else:
        print "Could not order %s and %s"%(a,b)
    

def AppendToDatabase(   path,
                        subhyp,
                        remote_script,
                        remote_database,
                        time,
                        testP,
                        seed,
                        local_pc,
                        Nlive=1000,
                        injection=0 #injnum
                        #misc. optional
                        #injfile=None,eventnum=None,skyres=None,
                       
                    ):
    
    ###### To do: I can actually make the safe_append.sh  to call a python to append. That will me able to check e.g. whether the line to be appended is already present 
    
    unwanted_pars=['psi','phi','logl']
    
    header=" "
    string_to_write=""
    posteriors_string=""
    string_to_write+=str(seed)+" "
    posteriors_string+=str(seed)+" "
    header+="InspinjSeed "
    header+="User "
    
    string_to_write+=str(getpass.getuser())+" "
    string_to_write+=str(local_pc)+" "
    posteriors_string+=str(getpass.getuser())+" "
    posteriors_string+=str(local_pc)+" "
    header+=str("cluster  time" )+ " "
    string_to_write+=str(time)+" "
    posteriors_string+=str(time)+" "
    found=0
    found_posteriors=0
    failed={}

    for hyp in sorted(subhyp,my_sort):
        path_to_snr=os.path.join(path,hyp,'SNR',"snr_H1L1V1_"+str(time)+".0.dat")
        path_to_Bfile=os.path.join(path,hyp,'nest',"outfile_"+str(time)+".000000_H1L1V1.dat_B.txt")
        path_to_outfile=os.path.join(path,hyp,'nest',"outfile_"+str(time)+".000000_H1L1V1.dat")
        path_to_headerfile=os.path.join(path,hyp,'nest',"outfile_"+str(time)+".000000_H1L1V1.dat_params.txt")
        if os.path.isfile(path_to_Bfile):
            bfile=np.loadtxt(path_to_Bfile)
            logB=bfile[0]
            string_to_write+="%.3f"%logB+" "
            header+=str(hyp)+" "
            found+=1
        else:
            print "WARNING: The Bfile %s was not found"%str(path_to_Bfile)
            failed[hyp]=str(path_to_Bfile)
            string_to_write+=str(path_to_Bfile)+" "
        if (os.path.isfile(path_to_outfile) and os.path.isfile(path_to_headerfile) and os.path.isfile(path_to_Bfile)):
            found_posteriors+=1
            data=path_to_outfile
            pos,d_all,totalBayes,ZnoiseTotal=combine_evidence([data],0,Nlive)
            headerfile=open(path_to_headerfile,'r')
            headerstr=headerfile.readline()
            headerfile.close()
            
            if  os.access(os.getenv("HOME"), os.W_OK):
                posfile_path=os.path.join(os.getenv("HOME"),'tmp.txt')
            else:
                print "Cannot write temporary files into %s!\n"%posfile_path
                ### I still need to catch a possible error here

            posfile=open(posfile_path,'w')
            posfile.write(headerstr+'\n')
            for row in pos:
                for i in row:
                    posfile.write('%f\t' %(i))
                posfile.write('\n')
            
            posfile.close()

            peparser=bppu.PEOutputParser('common')
            commonResultsObj=peparser.parse(open(posfile_path,'r'))
            injection=None
            pos = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injection)
            for name in pos.names:
                if not (name in unwanted_pars):
                    print "mean %s was %f std %f\n"%(name,pos[str(name)].mean,pos[str(name)].stdev)
                    posteriors_string+="%f %f "%(pos[str(name)].median,pos[str(name)].stdev)
            os.remove(posfile_path)

               
    if found < len(subhyp):
        print "ERROR, some of the B files were not found. This is probably due to some failed inspnest run(s)"
        remote_database=remote_database+"_failed"
 
    if os.path.isfile(path_to_snr):
        snrfile=np.loadtxt(path_to_snr,skiprows=3,usecols=(1,1))
        snr=snrfile[0]
        string_to_write+="%.2f"%snr+" "
    else:
        string_to_write+="SNRnotFound"
        ### If the SNR file it not present it just writes that in the string
    header+=str("NetSNR ")


    string_to_write+=str(testP)+" "
    header+=str("testParameter_ShiftPc ")
    remote_server=remote_script[0:remote_script.index(":")]
    path_remote_script=remote_script[remote_script.index(":")+1:]
    ### Writes inspinj seed, user cluster time  bayes and SNRs
    ### If any of the subhypotheses failed, it writes in a different file
    os.system("ssh "+remote_server + ' "'+ path_remote_script + " '" +string_to_write +"' "+ remote_database+'"') 
     
    if (not found < len(subhyp)) and (not found_posteriors< len(subhyp)) :
        ### If all the subhypotheses are ok (B files exist) AND the posterior files are present, write the parameters
        remote_database=remote_database+"_parameters"
        os.system("ssh "+remote_server + ' "'+ path_remote_script + " '" +posteriors_string +"' "+ remote_database+'"')  



if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-p","--path", dest="path",type="string",help="This dir contain all the sub-hypotheses dirs", metavar="DIR")
    parser.add_option("-s","--subhyp",dest="subhyp",action="callback",callback=vararg_callback,help="A space-separed list of the subhypothesis (comprised the GR)",metavar="GR dphi1 dphi2 dphi1dphi2")
    parser.add_option("-R","--remote-script",type="string",action="store",help="The path to the remote script which writes on the datafile", metavar="svitale@login.nikhef.nl:/project/gravwav/safe_append.sh")
    parser.add_option("-r","--remote-database",type="string",action="store",help="The remote database on which to append", metavar="datafile.txt")
    parser.add_option("-T","--testParam-ShiftPc",type="string",action="store",help="The shifted parameter underscore the percent value of the shift, or GR", metavar="dphi6_1pc")
    parser.add_option("-Q","--inspinj-seed",default=None,action="store",type="string",help="The unique value of the inspinj seed", metavar="700000")
    parser.add_option("-t","--event-time",type="string",action="store",help="The time of injection", metavar="939936910")
    parser.add_option("-c","--cluster",type="string",action="store",help="The cluster from which the pipeline is being run", metavar="Atlas1")
    parser.add_option("-N","--Nlive",default=1000,action="store",type="int",help="The number of live points used in the inspnest rum", metavar="1000")
    parser.add_option("-E","--eventnum",default=0,action="store",type="int",help="The eventnum in the injfile", metavar="0")
    parser.add_option("-I","--injfile",default=None,action="store",type="string",help="The injection file", metavar="My_injections.xml")
    (opts,args)=parser.parse_args()


    AppendToDatabase( 
                        opts.path,
                        opts.subhyp,
                        opts.remote_script,
                        opts.remote_database,
                        opts.event_time,
                        opts.testParam_ShiftPc,
                        opts.inspinj_seed,
                        opts.cluster,
                        Nlive=opts.Nlive,
                        injection=  opts.eventnum
                    )
#
