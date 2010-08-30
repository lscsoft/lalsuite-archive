#!/usr/bin/env python
#
#       bayespputils.py
#       
#       Copyright 2010 Benjamin Aylott <ben@star.sr.bham.ac.uk>
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
import os
from math import ceil,floor,sqrt

import numpy as np


from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm

import pylal

#Import from C extension module
from _bayespputils import skyhist_cart,calculateConfidenceLevels

def greedyBin2(posterior_array,par_bins,confidence_levels,par_names=None,injection=None):

    if par_names:
        par1_name,par2_name=par_names
    else:
        par1_name="Parameter 1"
        par2_name="Parameter 2"
    
    par1pos=posterior_array[:,0]
    par2pos=posterior_array[:,1]

    par1_bin,par2_bin=par_bins

    if injection:
        par1_injvalue,par2_injvalue=injection

    twoDGreedyCL={}
    twoDGreedyInj={}
    
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
    (injectionconfidence,toppoints,reses)=calculateConfidenceLevels(greedyHist,greedyPoints,injbin,float(sqrt(par1_bin*par2_bin)),confidence_levels,int(len(par1pos)))

    #Print confidence levels to file
    areastr=''
    for (frac,area) in reses:
        areastr+='%s,'%str(area)
        twoDGreedyCL[str(frac)]=area
    areastr=areastr.rstrip(',')
    

    if injection is not None and injectionconfidence is not None:
        twoDGreedyInj['confidence']=injectionconfidence
        #Recover area contained within injection point interval
        areasize=0
        while True:
            if injectionconfidence<np.asarray(toppoints)[areasize,3]:
                break
            areasize+=1
        areasize=areasize*par1_bin*par2_bin
        twoDGreedyInj['area']=areasize
    
    return toppoints,injectionconfidence,twoDGreedyCL,twoDGreedyInj
#

def pol2cart(long,lat):
    x=np.cos(lat)*np.cos(long)
    y=np.cos(lat)*np.sin(long)
    z=np.sin(lat)
    return np.array([x,y,z])
#

def skyhist_cart_slow(skycarts,sky_samples):
    """
    Histogram the list of samples into bins defined by Cartesian vectors in skycarts
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

def plotSkyMap(skypos,skyres,sky_injpoint,confidence_levels,outdir):

    from mpl_toolkits.basemap import Basemap
    from pylal import skylocutils

    np.seterr(under='ignore')

    skypoints=np.array(skylocutils.gridsky(float(skyres)))
    skycarts=map(lambda s: pol2cart(s[1],s[0]),skypoints)
    skyinjectionconfidence=None

    shist=skyhist_cart(np.array(skycarts),skypos)

    #shist=skyhist_cart(skycarts,list(pos))
    bins=skycarts

    # Find the bin of the injection if available
    injbin=None
    if sky_injpoint:
        injhist=skyhist_cart_slow(skycarts,np.array([sky_injpoint]))
        injbin=injhist.tolist().index(1)
        print 'Found injection in bin %d with co-ordinates %f,%f .'%(injbin,skypoints[injbin,0],skypoints[injbin,1])

    (skyinjectionconfidence,toppoints,skyreses)=calculateConfidenceLevels(shist,skypoints,injbin,float(skyres),confidence_levels,len(skypos))

    if injbin and skyinjectionconfidence:
        i=list(np.nonzero(np.asarray(toppoints)[:,2]==injbin))[0]

        min_sky_area_containing_injection=float(skyres)*float(skyres)*i
        print 'Minimum sky area containing injection point = %f square degrees'%min_sky_area_containing_injection

    myfig=plt.figure()
    plt.clf()
    m=Basemap(projection='moll',lon_0=180.0,lat_0=0.0)
    plx,ply=m(np.asarray(toppoints)[::-1,1]*57.296,np.asarray(toppoints)[::-1,0]*57.296)
    cnlevel=[1-tp for tp in np.asarray(toppoints)[::-1,3]]
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
    np.savetxt('ranked_sky_pixels',np.column_stack([np.asarray(toppoints)[:,0:1],np.asarray(toppoints)[:,1],np.asarray(toppoints)[:,3]]))

    return skyreses,skyinjectionconfidence
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


def plot2Dbins(toppoints,par_bins,outdir,par_names=None,injpoint=None):

    par1name,par2name=par_names

    par1_bin,par2_bin=par_bins
    #Work out good bin size
    xbins=int(ceil((max(toppoints[:,0])-min(toppoints[:,0]))/par1_bin))
    ybins=int(ceil((max(toppoints[:,1])-min(toppoints[:,1]))/par2_bin))

    if xbins==0:
        xbins=1
    if ybins==0:
        ybins=1

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

    coll=myfig.gca().scatter(toppoints[:,0],toppoints[:,1],s=int(points_per_bin_width*1.5),faceted=False,marker='s',c=cnlevel,cmap=mpl_cm.jet)
    plt.colorbar(mappable=coll,ax=myfig.gca(),cax=myfig.gca())

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
    plt.title("%s-%s histogram (greedy binning)"%(par1name,par2name)) # add a title

    if not os.path.isdir(os.path.join(outdir,'2Dbins')):
        os.makedirs(os.path.join(outdir,'2Dbins'))
    myfig.savefig(os.path.join(outdir,'2Dbins',par1name+'-'+par2name+'.png'),dpi=_dpi)
    myfig.clear()
    return
#

def mc2ms(mc,eta):
    
    root = np.sqrt(0.25-eta)
    fraction = (0.5+root) / (0.5-root)
    invfraction = 1/fraction

    m1= mc * np.power((1+fraction),0.2) / np.power(fraction,0.6)

    m2= mc* np.power(1+invfraction,0.2) / np.power(invfraction,0.6)
    return (m1,m2)
#
#
def ang_dist(long1,lat1,long2,lat2):
# Find the angular separation of (long1,lat1) and (long2,lat2)
# which are specified in radians
    
    x1=np.cos(lat1)*np.cos(long1)
    y1=np.cos(lat1)*np.sin(long1)
    z1=np.sin(lat1)
    x2=np.cos(lat2)*np.cos(long2)
    y2=np.cos(lat2)*np.sin(long2)
    z2=np.sin(lat2)
    sep=math.acos(x1*x2+y1*y2+z1*z2)
    return(sep)

#

def plot1DPDF(pos_samps,param,injpar=None,histbins=50):
    from scipy import stats
    from scipy import seterr as sp_seterr

    myfig=plt.figure(figsize=(4,3.5),dpi=80)

    (n, bins, patches)=plt.hist(pos_samps,histbins,normed='true')
    histbinSize=bins[1]-bins[0]

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    try:
        gkde=stats.kde.gaussian_kde(pos_samps)
    except np.linalg.linalg.LinAlgError:
        print "Error occured generating plot for parameter %s: %s ! Trying next parameter."%(param,'LinAlgError')
        return

    ind=np.linspace(min(pos_samps),max(pos_samps),101)
    kdepdf=gkde.evaluate(ind)
    plt.plot(ind,kdepdf,label='density estimate')

    rbins=None
    
    if injpar:
        if min(pos_samps)<injpar and max(pos_samps)>injpar:
            plt.plot([injpar,injpar],[0,max(kdepdf)],'r-.',scalex=False,scaley=False)

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
    
    return rbins,myfig#,rkde
#

def plot2Dkernel(xdat,ydat,Nx,Ny,par_names=None,par_injvalues=None):

    from scipy import seterr as sp_seterr
    from scipy import stats

    from matplotlib import pyplot as plt
    
    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    myfig=plt.figure(1,figsize=(6,4),dpi=80)
    plt.clf()

    xax=np.linspace(min(xdat),max(xdat),Nx)
    yax=np.linspace(min(ydat),max(ydat),Ny)
    x,y=np.meshgrid(xax,yax)
    samp=np.array([xdat,ydat])
    kde=stats.kde.gaussian_kde(samp)
    grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)

    z = kde(grid_coords.T)
    z = z.reshape(Nx,Ny)
    asp=xax.ptp()/yax.ptp()
#    if(asp<0.8 or asp > 1.6): asp=1.4
    plt.imshow(z,extent=(xax[0],xax[-1],yax[0],yax[-1]),aspect=asp,origin='lower')
    plt.colorbar()

    if par_injvalues:
        par_injvalue1,par_injvalue2=par_injvalues
        plt.plot([par_injvalue1],[par_injvalue2],'go',scalex=False,scaley=False)

    if par_names:
        par_name1,par_name2=par_names
        plt.xlabel(par_name1)
        plt.ylabel(par_name2)
    plt.grid()

    return myfig
#

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
#

def stacc_stat(par_samps,injvalue):
    return sqrt(np.var(par_samps)+pow((np.mean(par_samps)-injvalue),2) )
#


def greedyBin1(par_samps,par_bin,confidence_levels,par_injvalue=None):

    oneDGreedyCL={}
    oneDGreedyInj={}
    
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
            print "IndexError: bin number: %i total bins: %i parsamp: %f "%(par_binNumber,parpos_Nbins,par_samp)

    injbin=None
    #Find injection bin
    if par_injvalue:
        par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
        injbin=par_binNumber

    #Determine confidence levels
    (injectionconfidence,toppoints,reses)=calculateConfidenceLevels(greedyHist,greedyPoints,injbin,sqrt(par_bin),confidence_levels,len(par_samps))

    #areastr=''
    for (frac,area) in reses:
        oneDGreedyCL[frac]=area
    
    if par_injvalue and injectionconfidence:
        
        oneDGreedyInj['confidence']=injectionconfidence
        #Recover interval containing injection point
        interval=0
        while True:
            if injectionconfidence<np.asarray(toppoints)[interval,3]:
                break
            interval+=1
        interval=interval*par_bin
        oneDGreedyInj['interval']=interval

    return oneDGreedyCL,oneDGreedyInj,toppoints,injectionconfidence
    

#
def contCL1(par_samps,par_bin,confidence_levels,par_injvalue=None):

    oneDContCL={}
    oneDContInj={}

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
            print "IndexError: bin number: %i total bins: %i parsamp: %f bin: %f - %f"%(par_binNumber,parpos_Nbins,par_samp,greedyPoints[par_binNumber-1,0],greedyPoints[par_binNumber-1,0]+par_bin)

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
            #summary_file.set('1D contigious cl',par_name,'['+str(max_left*par_bin)+','+str(max_right*par_bin)+','+str((max_right-max_left)*par_bin)+']')
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
#
