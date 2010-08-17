#!/usr/bin/env python

# Load samples from the posterior PDF of an MCMC or nested sampling code
# and produce sky localisation plots and size estimates
# Format for MCMC files should be:
# log(mc) eta tc phase log(dist) RA dec psi iota

#List of parameters to plot/bin . Need to match (converted) column names.
oneDMenu=['mtotal','m1','m2','mchirp','mc','distance','distMPC','dist','iota','eta','RA','dec','a1','a2','phi1','theta1','phi2','theta2']
#List of parameter pairs to bin . Need to match (converted) column names.
twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['RA','dec']]
#Bin size/resolution for binning. Need to match (converted) column names.
GreedyRes={'mc':0.025,'m1':0.1,'m2':0.1,'mtotal':0.1,'eta':0.001,'iota':0.01,'time':1e-6,'distance':0.01,'dist':0.1,'mchirp':0.025,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'RA':0.005,'dec':0.005}
#Confidence levels
confidence_levels=[0.67,0.9,0.95,0.99]
#2D plots list
twoDplots=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['m1','dist'],['m2','dist']]

debug=False

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

from scipy import stats

from pylal import SimInspiralUtils
from pylal import bayespputils

parser=OptionParser()
parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
parser.add_option("-d","--data",dest="data",action="append",help="datafile")
parser.add_option("-i","--inj",dest="injfile",help="SimInsipral injection file",metavar="INJ.XML",default=None)
parser.add_option("--skyres",dest="skyres",help="Sky resolution to use to calculate sky box size",default=None)
parser.add_option("--eventnum",dest="eventnum",action="store",default=None,help="event number in SimInspiral file of this signal",type="int",metavar="NUM")

(opts,args)=parser.parse_args()

if opts.eventnum is not None and opts.injfile is None:
    print "You specified an event number but no injection file. Ignoring!"

if opts.data is None:
    print 'You must specify an input data file'
    exit(1)

def mc2ms(mc,eta):
    root = sqrt(0.25-eta)
    fraction = (0.5+root) / (0.5-root)
    invfraction = 1/fraction

    m1= mc * pow((1+fraction),0.2) / pow(fraction,0.6)

    m2= mc* pow(1+invfraction,0.2) / pow(invfraction,0.6)
    return (m1,m2)

def histN(mat,N):
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

outdir=opts.outpath

def ang_dist(long1,lat1,long2,lat2):
# Find the angular separation of (long1,lat1) and (long2,lat2)
# which are specified in radians
    x1=cos(lat1)*cos(long1)
    y1=cos(lat1)*sin(long1)
    z1=sin(lat1)
    x2=cos(lat2)*cos(long2)
    y2=cos(lat2)*sin(long2)
    z2=sin(lat2)
    sep=math.acos(x1*x2+y1*y2+z1*z2)
    return(sep)

def pol2cart(long,lat):
    x=cos(lat)*cos(long)
    y=cos(lat)*sin(long)
    z=sin(lat)
    return array([x,y,z])

def sky_hist(skypoints,samples):
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
        print 'Done %d/%d iterations, minsep=%f degrees'%(j,len(samples),minsep*(180.0/3.1415926))
    return (skypoints,bins)

def skyhist_cart(skycarts,sky_samples):
    """
    Histogram the list of samples into bins defined by Cartesian vectors in skycarts
    """
    dot=np.dot
    N=len(skycarts)
    print 'operating on %d sky points'%(N)
    bins=np.zeros(N)
    for RAsample,decsample in sky_samples:
        sampcart=pol2cart(RAsample,decsample)
        maxdx=-1
        maxvalue=-1
        for i in xrange(0,N):
            dx=dot(sampcart,skycarts[i])
            if dx>maxvalue:
                    maxdx=i
                    maxvalue=dx
            #if debug:
                #print "sky co : "+str(skycarts[i])+"long : "+str(sample[RAdim])+" "+str(sample[decdim])
                #print "sample co: "+str(sampcart)
        if debug:
            print "maxdx : %i , maxvalue: %f sample0 : %f sample1: %f"%(maxdx,maxvalue,sample[RAdim],sample[decdim])
        bins[maxdx]+=1
    return bins

def loadDataFile(filename):
    print filename
    infile=open(filename,'r')
    formatstr=infile.readline().lstrip()
    formatstr=formatstr.replace('#','')
    header=formatstr.split()
  
    llines=[]
    import re
    dec=re.compile(r'[^\d.-]+')
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

summary_fo=open('summary.ini','w')

summary_file=ConfigParser()

summary_file.add_section('metadata')
summary_file.set('metadata','group_id','X')
if opts.eventnum:
    summary_file.set('metadata','event_id',str(opts.eventnum))
summary_file.add_section('Confidence levels')
summary_file.set('Confidence levels','confidence levels',str(confidence_levels))
 
# Load in the main data
paramnames, pos=loadDataFile(opts.data[0])

#Generate any required derived parameters
if "m1" not in paramnames and "m2" not in paramnames and "mchirp" in paramnames and "eta" in paramnames:
    (m1,m2)=mc2ms(pos[:,paramnames.index('mchirp')],pos[:,paramnames.index('eta')])
    
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
if(opts.injfile is not None):
    import itertools
    injections = SimInspiralUtils.ReadSimInspiralFromFiles([opts.injfile])
    if(opts.eventnum is not None):
        if(len(injections)<opts.eventnum):
            print "Error: You asked for event %d, but %s contains only %d injections" %(opts.eventnum,opts.opts.injfile,len(injections))
            sys.exit(1)
        else:
            injection=injections[opts.eventnum]
    else:
        if(len(injections)<1):
            print 'Warning: Cannot find injection with end time %f' %(means[2])
        else:
            injection = itertools.ifilter(lambda a: abs(a.get_end() - means[2]) < 0.1, injections).next()

def getinjpar(inj,parnum):
    if paramnames[parnum]=='mchirp' or paramnames[parnum]=='mc': return inj.mchirp
    if paramnames[parnum]=='mass1' or paramnames[parnum]=='m1':
        (m1,m2)=mc2ms(inj.mchirp,inj.eta)
        return m1
    if paramnames[parnum]=='mass2' or paramnames[parnum]=='m2':
        (m1,m2)=mc2ms(inj.mchirp,inj.eta)
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

if injection:
    injpoint=map(lambda a: getinjpar(injection,a),range(len(paramnames)))
    injvals=map(str,injpoint)
    out=reduce(lambda a,b:a+'||'+b,injvals)
    print 'Injected values:'
    print out

    #Add injection values to output file
    summary_file.add_section('Injection values')
    
    for parnum in range(len(paramnames)):
        summary_file.set('Injection values',paramnames[parnum],getinjpar(injection,parnum))
    
#

def calculateSkyConfidence_slow(shist,skypoints,injbin,skyres_,confidence_levels,lenpos):
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
            #print 'Nbins=%d, thisnum=%d, idx=%d, total=%d, cumul=%f\n'%(Nbins,maxbin,maxpos,len(pos),frac)
        print '%f confidence region: %f square degrees' % (frac,Nbins*float(skyres_)*float(skyres_))
        skyreses.append((frac,Nbins*float(skyres_)*float(skyres_)))
        toppoints=toppoints[:Nbins]
    return injectionconfidence,toppoints,skyreses
#
def plot2Dbins(toppoints,cl,par1_name,par1_bin,par2_name,par2_bin,injpoint):

    #Work out good bin size
    xbins=int(ceil((max(toppoints[:,0])-min(toppoints[:,0]))/par1_bin))
    ybins=int(ceil((max(toppoints[:,1])-min(toppoints[:,1]))/par2_bin))

    _dpi=120
    xsize_in_inches=6.
    xsize_points = xsize_in_inches * _dpi
    points_per_bin_width=xsize_points/xbins

    ysize_points=ybins*points_per_bin_width
    ysize_in_inches=ysize_points/_dpi
    #

    myfig=plt.figure(1,figsize=(xsize_in_inches+2,ysize_in_inches+2),dpi=_dpi,autoscale_on=True)

    cnlevel=[1-tp for tp in toppoints[:,3]]
    #

    myfig.gca().scatter(toppoints[:,0],toppoints[:,1],s=int(points_per_bin_width*1.5),faceted=False,marker='s',c=cnlevel,cmap=matplotlib.cm.jet)
    plt.colorbar()

    #Determine limits based on injection point (if any) and min/max values

    min_xlim=min(toppoints[:,0])
    max_xlim=max(toppoints[:,0])

    min_ylim=min(toppoints[:,1])
    max_ylim=max(toppoints[:,1])

    if injpoint is not None:
        myfig.gca().plot([injpoint[0]],[injpoint[1]],'r*',ms=20.)

        if injpoint[0] < min(toppoints[:,0]):
            min_xlim=injpoint[0]
        elif injpoint[0] > max(toppoints[:,0]):
            max_xlim=injpoint[0]

        if injpoint[1] < min(toppoints[:,1]):
            min_ylim=injpoint[1]
        elif injpoint[1] > max(toppoints[:,1]):
            max_ylim=injpoint[1]
#
    #Set limits on axes determined above
    myfig.gca().set_xlim(min_xlim,max_xlim)
    myfig.gca().set_ylim(min_ylim,max_ylim)

    #Reset figure size (above probably had no effect apart from to give correct bin sizes)
    myfig.set_figheight(6)
    myfig.set_figwidth(6)
    plt.title("%s-%s histogram (greedy binning)"%(par1_name,par2_name)) # add a title

    if not os.path.isdir(os.path.join(outdir,'2Dbins')):
        os.makedirs(os.path.join(outdir,'2Dbins'))
    myfig.savefig(os.path.join(outdir,'2Dbins',par1_name+'-'+par2_name+'.png'),dpi=_dpi)
    myfig.clear()
    return
#
#calculateSkyConfidence=calculateSkyConfidence_slow

def plotSkyMap(skypos,skyres,sky_injpoint):
    from pylal import skylocutils
    from mpl_toolkits.basemap import Basemap

    skypoints=array(skylocutils.gridsky(float(skyres)))
    skycarts=map(lambda s: pol2cart(s[1],s[0]),skypoints)
    skyinjectionconfidence=None
    
    shist=bayespputils.skyhist_cart(array(skycarts),skypos)

    #shist=skyhist_cart(skycarts,list(pos))
    bins=skycarts

    # Find the bin of the injection if available
    injbin=None
    if sky_injpoint:
        injhist=skyhist_cart(skycarts,array([sky_injpoint]))
        injbin=injhist.tolist().index(1)
        print 'Found injection in bin %d with co-ordinates %f,%f .'%(injbin,skypoints[injbin,0],skypoints[injbin,1])

    (skyinjectionconfidence,toppoints,skyreses)=bayespputils.calculateConfidenceLevels(shist,skypoints,injbin,float(opts.skyres),confidence_levels,len(pos))
    
    if injbin and skyinjectionconfidence:
        i=list(np.nonzero(np.asarray(toppoints)[:,2]==injbin))[0]
        
        min_sky_area_containing_injection=float(opts.skyres)*float(opts.skyres)*i
        print 'Minimum sky area containing injection point = %f square degrees'%min_sky_area_containing_injection

    myfig=plt.figure()
    plt.clf()
    m=Basemap(projection='moll',lon_0=180.0,lat_0=0.0)
    plx,ply=m(np.asarray(toppoints)[::-1,1]*57.296,np.asarray(toppoints)[::-1,0]*57.296)
    cnlevel=[1-tp for tp in np.asarray(toppoints)[::-1,3]]
    plt.scatter(plx,ply,s=5,c=cnlevel,faceted=False,cmap=matplotlib.cm.jet)
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
    np.savetxt('ranked_sky_pixels',column_stack([np.asarray(toppoints)[:,0:1],np.asarray(toppoints)[:,1],np.asarray(toppoints)[:,3]]))
    
    return skyreses,skyinjectionconfidence
#
skyreses=None
if opts.skyres is not None:

    RAvec=array(pos)[:,paramnames.index('RA')]
    decvec=array(pos)[:,paramnames.index('dec')]
    skypos=column_stack([RAvec,decvec])
    if injection:
        skyreses,skyinjectionconfidence=plotSkyMap(skypos,opts.skyres,(injpoint[RAdim],injpoint[decdim]))
    else:
        skyreses,skyinjectionconfidence=plotSkyMap(skypos,opts.skyres,None)
#

myfig=plt.figure(1,figsize=(6,4),dpi=80)
plt.clf()

#Bin 2D marginal pdfs and determine confidence intervals

summary_file.add_section('2D greedy cl')

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

    par1pos=pos_array[:,par1_index]
    par2pos=pos_array[:,par2_index]
    #Create 2D bin array
    par1pos_min=min(par1pos)
    par2pos_min=min(par2pos)

    par1pos_max=max(par1pos)
    par2pos_max=max(par2pos)

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

    injbin=None
    #if injection point given find which bin its in
    if injection:
        par1_injvalue=np.array(injpoint)[par1_index]
        par2_injvalue=np.array(injpoint)[par2_index]

        if par1_injvalue is not None and par2_injvalue is not None:

            par1_binNumber=floor((par1_injvalue-par1pos_min)/par1_bin)
            par2_binNumber=floor((par2_injvalue-par2pos_min)/par2_bin)

            injbin=int(par1_binNumber+par2_binNumber*par1pos_Nbins)
        elif par1_injvalue is None and par2_injvalue is not None:
            print "Injection value not found for %s!"%par1_name

        elif par1_injvalue is not None and par2_injvalue is None:
            print "Injection value not found for %s!"%par2_name

    #Bin posterior samples
    for par1_samp,par2_samp in zip(par1pos,par2pos):
        par1_binNumber=floor((par1_samp-par1pos_min)/par1_bin)
        par2_binNumber=floor((par2_samp-par2pos_min)/par2_bin)
        greedyHist[par1_binNumber+par2_binNumber*par1pos_Nbins]+=1

    #Now call usual confidence level function

    #print greedyHist,greedyPoints,injbin,sqrt(par1_bin*par2_bin),confidence_levels,len(par1pos)
    (injectionconfidence,toppoints,reses)=bayespputils.calculateConfidenceLevels(greedyHist,greedyPoints,injbin,float(sqrt(par1_bin*par2_bin)),confidence_levels,int(len(par1pos)))

    #Print confidence levels to file
    areastr=''
    for (frac,area) in reses:
        areastr+='%s,'%str(area)
    areastr=areastr.rstrip(',')
    summary_file.set('2D greedy cl',par1_name+','+par2_name,str(areastr))
    if injection is not None and injectionconfidence is not None:
        summary_file.set('2D greedy cl',par1_name+','+par2_name+'_inj',str(injectionconfidence))
    #Plot 2D histogram
    if injection is not None and par1_injvalue is not None and par2_injvalue is not None:
        plot2Dbins(np.array(toppoints),confidence_levels,par1_name,par1_bin,par2_name,par2_bin,[par1_injvalue,par2_injvalue])
    else:
        plot2Dbins(np.array(toppoints),confidence_levels,par1_name,par1_bin,par2_name,par2_bin,None)
#


#1D binning
summary_file.add_section('1D mean')
summary_file.add_section('1D median')
summary_file.add_section('1D mode')
summary_file.add_section('1D contigious cl')
summary_file.add_section('1D greedy cl')
summary_file.add_section('1D stacc')
max_i=0
max_pos=pos_array[0,-1] #???
par_samps=pos_array[:,-1]
for i in range(len(pos_array)):

    if par_samps[i] > max_pos:
        max_pos=par_samps[i]
        max_i=i

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

    par_samps=pos_array[:,par_index]
    summary_file.set('1D mode',par_name,str(par_samps[max_i]))
    summary_file.set("1D mean",par_name,str(np.mean(par_samps)))
    summary_file.set("1D median",par_name,str(np.median(par_samps)))

    parpos_min=min(par_samps)

    parpos_max=max(par_samps)
    par_point=parpos_min
    parpos_Nbins= int(ceil((parpos_max - parpos_min)/par_bin))+1

    greedyPoints=np.zeros((parpos_Nbins,2)) #2D so it can be put through same confidence level function
    greedyHist=np.zeros(parpos_Nbins,dtype='i8')

    #Bin up
    for i in range(parpos_Nbins):

        greedyPoints[i,0]=par_point
        greedyPoints[i,1]=par_point
        par_point+=par_bin

    for par_samp in par_samps:
        par_binNumber=int(floor((par_samp-parpos_min)/par_bin))
        try:
            greedyHist[par_binNumber]+=1
        except IndexError:
            print "IndexError: bin number: %i total bins: %i parsamp: %f bin: %f - %f"%(par_binNumber,parpos_Nbins,par_samp,greedyPoints[par_binNumber-1,0],greedyPoints[par_binNumber-1,0]+par_bin)
    injbin=None
    #Find injection bin
    if injection:
        #print par_index,numpy.array(injpoint)
        #print paramnames

        par_injvalue=np.array(injpoint)[par_index]
        if par_injvalue is not None:
            par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
            injbin=par_binNumber

    j=0
    print "Calculating contigious confidence intervals for %s..."%par_name
    len_par_samps=len(par_samps)

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
                if frac>confidence_level:
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
            if max_frac is not None:
                break

        if max_frac is None:
            print "Cant determine intervals at %f confidence!"%confidence_level
        else:
            summary_file.set('1D contigious cl',par_name,'['+str(max_left*par_bin)+','+str(max_right*par_bin)+','+str((max_right-max_left)*par_bin)+']')

            k=j
            while k+1<len(confidence_levels) :
                if confidence_levels[k+1]<max_frac:
                    j+=1
                k+=1
        j+=1


    #Determine confidence levels

    (injectionconfidence,toppoints,reses)=bayespputils.calculateConfidenceLevels(greedyHist,greedyPoints,injbin,sqrt(par_bin),confidence_levels,len(par_samps))

    areastr=''
    for (frac,area) in reses:
        areastr+='%s,'%str(area)
    areastr=areastr.rstrip(',')
    summary_file.set('1D greedy cl',par_name,'['+str(areastr)+']')

    if injection is not None and injectionconfidence is not None:
        summary_file.set('1D greedy cl',par_name+'_inj',str(injectionconfidence))

    #Ilya's standard accuracy statistic
    if injection:
        injvalue=np.array(injpoint)[par_index]
        if injvalue is not None:
            stacc=sqrt(np.var(par_samps)+pow((np.mean(par_samps)-injvalue),2) )

            summary_file.set('1D stacc',par_name,str(stacc))

#



def plot2Dkernel(xdat,ydat,Nx,Ny):
    xax=np.linspace(min(xdat),max(xdat),Nx)
    yax=np.linspace(min(ydat),max(ydat),Ny)
    x,y=np.meshgrid(xax,yax)
    samp=array([xdat,ydat])
    kde=stats.kde.gaussian_kde(samp)
    grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)
    
    z = kde(grid_coords.T)
    z = z.reshape(Nx,Ny)
    asp=xax.ptp()/yax.ptp()
#    if(asp<0.8 or asp > 1.6): asp=1.4
    plt.imshow(z,extent=(xax[0],xax[-1],yax[0],yax[-1]),aspect=asp,origin='lower')
    plt.colorbar()
#


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
    
    plot2Dkernel(pos[:,i],pos[:,j],50,50)
    
    if injection and reduce (lambda a,b: a and b, map(lambda idx: getinjpar(injection,idx)>min(pos[:,idx]) and getinjpar(injection,idx)<max(pos[:,idx]),[i,j])) :
        if getinjpar(injection,i) is not None and getinjpar(injection,j) is not None:
            plt.plot([getinjpar(injection,i)],[getinjpar(injection,j)],'go',scalex=False,scaley=False)
    
    plt.xlabel(paramnames[i])
    plt.ylabel(paramnames[j])
    plt.grid()
    margdir=os.path.join(outdir,'2D')
    
    if not os.path.isdir(margdir+'/'): 
        os.mkdir(margdir)
    
    myfig.savefig(os.path.join(margdir,paramnames[i]+'-'+paramnames[j]+'_2Dkernel.png'))
    myfig.clear()

htmlfile=open(os.path.join(outdir,'posplots.html'),'w')
htmlfile.write('<HTML><HEAD><TITLE>Posterior PDFs</TITLE></HEAD><BODY><h3>'+str(means[2])+' Posterior PDFs</h3>')
if(opts.skyres is not None):
    htmlfile.write('<table border=1><tr><td>Confidence region<td>size (sq. deg)</tr>')
    for (frac,skysize) in skyreses:
        htmlfile.write('<tr><td>%f<td>%f</tr>'%(frac,skysize))
    htmlfile.write('</table>')
htmlfile.write('Produced from '+str(size(pos,0))+' posterior samples.<br>')
htmlfile.write('Samples read from %s<br>'%(opts.data[0]))
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
htmlfile.write('<table border=1><tr>')
#htmlfile.write('<td width=30%><img width=100% src="m1m2.png"></td>')
#htmlfile.write('<td width=30%><img width=100% src="RAdec.png"></td>')
#htmlfile.write('<td width=30%><img width=100% src="Meta.png"></td>')
#htmlfile.write('</tr><tr><td width=30%><img width=100% src="2D/Mchirp (Msun)-geocenter time ISCO_2Dkernel.png"</td>')
if opts.skyres is not None:
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

for param in oneDMenu:

    try:
        par_index=paramnames.index(param)
        i=par_index
    except ValueError:
        print "No input chain for %s, skipping 1D plot."%param
        continue

    myfig=plt.figure(figsize=(4,3.5),dpi=80)

    histbins=50

    #histbinSize=max(pos[:,i])-min(pos[:,i])/float(histbins)
    
    (n, bins, patches)=plt.hist(pos[:,i],histbins,normed='true')
    histbinSize=bins[1]-bins[0]

    try:
        gkde=stats.kde.gaussian_kde(pos[:,i])
        pass
    except np.linalg.linalg.LinAlgError:
        print "Error occured generating plot for parameter %s: %s ! Trying next parameter."%(paramnames[i],'LinAlgError')
        continue

    print "Generating 1D plot for %s."%param
    ind=np.linspace(min(pos[:,i]),max(pos[:,i]),101)
    kdepdf=gkde.evaluate(ind)
    plt.plot(ind,kdepdf,label='density estimate')
    if injection and min(pos[:,i])<getinjpar(injection,i) and max(pos[:,i])>getinjpar(injection,i):
        plt.plot([getinjpar(injection,i),getinjpar(injection,i)],[0,max(kdepdf)],'r-.',scalex=False,scaley=False)
        
        rkde=gkde.integrate_box_1d(min(pos[:,i]),getinjpar(injection,i))
        print "r of injected value of %s (kde) = %f"%(param,rkde)

        #Find which bin the true value is in
        bins_to_inj=(getinjpar(injection,i)-bins[0])/histbinSize
        injbinh=int(floor(bins_to_inj))
        injbin_frac=bins_to_inj-float(injbinh)
        
        #Integrate over the bins
        rbins=(sum(n[0:injbinh-1])+injbin_frac*n[injbinh])*histbinSize
        
        print "r of injected value of %s (bins) = %f"%(param, rbins)

        summary_file.set('1D ranking kde',param,rkde)
        summary_file.set('1D ranking bins',param,rbins)
#
    plt.grid()
    plt.xlabel(paramnames[i])
    plt.ylabel('Probability Density')
    myfig.savefig(os.path.join(outdir,paramnames[i]+ '.png'))
    myfig=plt.figure(figsize=(4,3.5),dpi=80)
    plt.plot(pos[:,i],'.')
    if injection and min(pos[:,i])<getinjpar(injection,i) and max(pos[:,i])>getinjpar(injection,i):
        plt.plot([0,len(pos)],[getinjpar(injection,i),getinjpar(injection,i)],'r-.')
    myfig.savefig(os.path.join(outdir,paramnames[i]+'_samps.png'))
    htmlfile.write('<img src="'+paramnames[i]+'.png"><img src="'+paramnames[i]+'_samps.png"><br>')

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

