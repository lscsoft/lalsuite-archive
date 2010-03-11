#!/usr/bin/python

#from numpy import *
import scipy
import matplotlib 
#matplotlib.use("Agg")
import sys
import math
from pylab import *
from optparse import OptionParser
import os
import numpy
from time import strftime
from scipy import stats
from pylal import SimInspiralUtils

parser=OptionParser()
parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
parser.add_option("-N","--Nlive",dest="Nlive",help="number of live points for each of the files")
parser.add_option("-d","--data",dest="data",action="append",help="datafile")
parser.add_option("-i","--inj",dest="injfile",help="SimInsipral injection file",metavar="INJ.XML",default=None)
parser.add_option("--inco0",dest="inco0",action="append",help="single-ifo runs for 0th ifo")
parser.add_option("--inco1",dest="inco1",action="append",help="single-ifo runs for 1th ifo")
parser.add_option("--inco2",dest="inco2",action="append",help="single-ifo runs for 2th ifo")
parser.add_option("--inco3",dest="inco3",action="append",help="single-ifo runs for 3th ifo")
parser.add_option("--skyres",dest="skyres",help="Sky resolution to use to calculate sky box size",default=None)
parser.add_option("--eventnum",dest="eventnum",action="store",default=None,help="event number in SimInspiral file of this signal",type="int",metavar="NUM")

(opts,args)=parser.parse_args()

if opts.eventnum is not None and opts.injfile is None:
    print "You specified an event number but no injection file. Ignoring!"

def logadd(a,b):
    if(a>b): (a,b)=(b,a)
    return (b+log(1+exp(a-b)))

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

def nest2pos_par(samps,weights):
    randoms=rand(size(samps,0))
    wt=weights+samps[:,-1]
    maxwt=max(wt)
    posidx=find(wt>maxwt+log(randoms))
    pos=samps[posidx,:]
    return pos

def nest2pos(samps,Nlive):
    weight = -linspace(1,len/Nlive,len)
    weight = weight + samps[:,-1]
    maxwt = max(weight)
    randoms = rand(len)
    pos = zeros(size(samps,1))
    posidx = find(weight>maxwt+log(randoms))
    pos=samps[posidx,:]
    return pos

def nestZ(d,Nlive):
    logw = log(1 - exp(-1.0/Nlive))
    logZ = logw + d[0,-1]
    logw = logw - 1.0/Nlive
    len=size(d,0)
    H=0
    for i in linspace(1,len-2,len):
        logZnew=logadd(logZ,logw+d[i,-1])
        H = exp(logw + d[i,-1] -logZnew)*d[i,-1] \
            + exp(logZ-logZnew)*(H+logZ) - logZnew
        logw = logw - 1.0/Nlive
        logZ=logZnew
    return (logZ,H)

def nestPar(d,Nlive):
    """
    Do nested sampling on parallel runs, taking into account
    different run lengths
    """
    maxes=[]
    for set in d:
        maxes.append(set[-1,-1]) # Max L point for each chain
    maxes=array(maxes)
    maxes.sort()
    N=len(d)
    print 'N chains = '+str(N)
    logw = log(1.0-exp(-1.0/(N*Nlive)))
    H=0
    alldat=reduce(lambda x,y: hstack([x,y]) , map(lambda x: x[:,-1],d))
    sidx=argsort(alldat)
    alldat=alldat[sidx]
    logZ=logw + alldat[0]
    logw-=1.0/(N*Nlive)
    weights = zeros(size(alldat,0))
    weights[0]=logw
    j=0
    numsamp=size(alldat,0)
    for samp in alldat[1:]:
        if samp>maxes[0]:
            maxes=maxes[1:]
            N-=1                
            print str(N)+' Parallel chains left at %d/%d'%(j,numsamp)
        logZnew = logadd(logZ,logw+samp)
        H = exp(logw + samp -logZnew)*samp \
            + exp(logZ-logZnew)*(H+logZ) - logZnew
        logw = logw -1.0/(N*Nlive)
        j+=1
        weights[j]=logw
        logZ=logZnew
    bigdata=reduce(lambda x,y: vstack([x,y]), d)
    bigdata=bigdata[sidx]
    return (logZ,H,bigdata,weights)

incoflag=0
outdir=opts.outpath
Nlive=int(opts.Nlive)

inco=[]
iBfiles=[]

def getBfile(datname):
    Bfile=datname+'_B.txt'
    print 'Looking for '+Bfile
    if os.access(Bfile,os.R_OK):
        outstat = loadtxt(Bfile)
        return outstat
    else:
        return None

def loaddata(datalist):
    out = list(map(loadtxt,datalist))
    Bfiles = list(map(getBfile,datalist))
    if not None in Bfiles: # Subtract off the noise evidence
        for (outfile,Bfile) in zip(out,Bfiles):
            outfile[:,-1]-=Bfile[2]
    return out,Bfiles

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
		seps=map(lambda s: ang_dist(sample[5],sample[6],s[1],s[0]),skypoints)
		minsep=math.pi
		for i in range(0,N):
			if seps[i]<minsep:
				minsep=seps[i]
				mindx=i
		bins[mindx]=bins[mindx]+1
		j=j+1
		print 'Done %d/%d iterations, minsep=%f degrees'%(j,len(samples),minsep*(180.0/3.1415926))
	return (skypoints,bins)

def skyhist_cart(skycarts,samples):
	N=len(skypoints)
	print 'operating on %d sky points'%(N)
	bins=zeros(N)
	j=0
	for sample in samples:
		sampcart=pol2cart(sample[5],sample[6])
		dots=map(lambda s: numpy.dot(sampcart,s),skycarts)
		maxdot=0
		for i in range(0,N):
			if dots[i]>maxdot:
				maxdot=dots[i]
				mindx=i
		bins[mindx]=bins[mindx]+1
		j=j+1
	#	print 'Done %d/%d iterations, minsep=%f degrees'%(j,len(samples),math.acos(maxdot)*(180.0/3.14159))
	return (skypoints,bins)

# Load in the main data
(d,Bfiles)=loaddata(opts.data)
if not None in Bfiles:
    Bflag=1
else:
    Bflag=0
# Load in the single-IFo data if specified
if opts.inco0 is not None:
    incoflag=1
    (a,iBfile)=loaddata(opts.inco0)
    inco.append(a)
    iBfiles.append(iBfile)
if opts.inco1 is not None:
    incoflag=2
    (a,iBfile)=loaddata(opts.inco1)
    inco.append(a)
    iBfiles.append(iBfile)
if opts.inco2 is not None:
    incoflag=3
    (a,iBfile)=loaddata(opts.inco2)
    inco.append(a)
    iBfiles.append(iBfile)
if opts.inco3 is not None:
    incoflag=4
    (a,iBfile)=loaddata(opts.inco3)
    inco.append(a)
    iBfiles.append(iBfile)

#len=size(d,0)
Nd=size(d[0],1)
#sidx=argsort(d[:,9])
#d=d[sidx,:]
#d[:,0]=exp(d[:,0])
print 'Exponentiated mc'
#maxL = max(d[-1,-1])

Zinco=0
for incox in inco:
    (Zincox,Hinco,d_inco_s,d_weights)=nestPar(incox,Nlive)
    Zinco=Zinco+Zincox

print 'Applying parallelised nested sampling algorithm to samples'

(logZ,H,d_sorted,d_weights)=nestPar(d,Nlive)

d_sorted[:,0]=exp(d_sorted[:,0])
d_sorted[:,4]=exp(d_sorted[:,4])
maxL=d_sorted[-1,-1]
print 'maxL = ' + str(maxL)
# Maximum likelihood point
print 'Max Likelihood point:'
maxPt= map(str,d_sorted[-1,:])
out=reduce(lambda a,b: a + ' || ' + b,maxPt)
print '|| ' + out + ' ||'

pos = nest2pos_par(d_sorted,d_weights)

print "Number of posterior samples: " + str(size(pos,0))
# Calculate means
means = mean(pos,axis=0)
meanStr=map(str,means)
out=reduce(lambda a,b:a+'||'+b,meanStr)
print 'Means:'
print '||'+out+'||'


#pos[:,2]=pos[:,2]-means[2]
injection=None
# Select injections using tc +/- 0.1s if it exists
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
    if parnum==0: return inj.mchirp
    if parnum==1: return inj.eta
    if parnum==2: return inj.get_end()
    if parnum==3: return inj.phi0
    if parnum==4: return inj.distance
    if parnum==5: return inj.longitude
    if parnum==6: return inj.latitude
    if parnum==7: return inj.polarization
    if parnum==8: return inj.inclination
    return None

if injection:
    injvals=map(str,map(lambda a: getinjpar(injection,a),range(0,9)))
    out=reduce(lambda a,b:a+'||'+b,injvals)
    print 'Injected values:'
    print out

if(Bflag==1):
    BayesFactor = logZ
    print 'log B = '+str(BayesFactor)

skyreses=[]
if(opts.skyres is not None):
	from pylal import skylocutils
	skypoints=array(skylocutils.gridsky(float(opts.skyres)))
	skycarts=map(lambda s: pol2cart(s[1],s[0]),skypoints)
	(bins,shist)=skyhist_cart(skycarts,pos)
	#(bins,hist)=sky_hist(skypoints,pos)
	frac=0
	Nbins=0
	toppoints=[]
	while(frac<0.67):
		maxbin=0
		for i in range(0,len(bins)):
			if shist[i]>maxbin:
				maxbin=shist[i]
				maxpos=i
		shist[maxpos]=0
		frac=frac+(float(maxbin)/float(len(pos)))
		Nbins=Nbins+1
		toppoints.append((skypoints[maxpos,0],skypoints[maxpos,1],maxbin))
		#print 'Nbins=%d, thisnum=%d, idx=%d, total=%d, cumul=%f\n'%(Nbins,maxbin,maxpos,len(pos),frac)
	print '%f confidence region: %f square degrees' % (frac,Nbins*float(opts.skyres)*float(opts.skyres))
	skyreses.append((frac,Nbins*float(opts.skyres)*float(opts.skyres)))
	while(frac<0.9):
                maxbin=0
                for i in range(0,len(bins)):
                        if shist[i]>maxbin:
                                maxbin=shist[i]
                                maxpos=i
                shist[maxpos]=0
                frac=frac+(float(maxbin)/float(len(pos)))
                Nbins=Nbins+1
		toppoints.append((skypoints[maxpos,0],skypoints[maxpos,1],maxbin))
		#print 'Nbins=%d, thisnum=%d, idx=%d, total=%d, cumul=%f\n'%(Nbins,maxbin,maxpos,len(pos),frac)
        print '%f confidence region: %f square degrees' % (frac,Nbins*float(opts.skyres)*float(opts.skyres))
        skyreses.append((frac,Nbins*float(opts.skyres)*float(opts.skyres)))
	while(frac<0.95):
                maxbin=0
                for i in range(0,len(bins)):
                        if shist[i]>maxbin:
                                maxbin=shist[i]
                                maxpos=i
                shist[maxpos]=0
                frac=frac+(float(maxbin)/float(len(pos)))
                Nbins=Nbins+1
		toppoints.append((skypoints[maxpos,0],skypoints[maxpos,1],maxbin))
		#print 'Nbins=%d, thisnum=%d, idx=%d, total=%d, cumul=%f\n'%(Nbins,maxbin,maxpos,len(pos),frac)
        print '%f confidence region: %f square degrees' % (frac,Nbins*float(opts.skyres)*float(opts.skyres))
        skyreses.append((frac,Nbins*float(opts.skyres)*float(opts.skyres)))
    
myfig=figure(1,figsize=(6,4),dpi=80)

def plot2Dkernel(xdat,ydat,Nx,Ny):
    xax=linspace(min(xdat),max(xdat),Nx)
    yax=linspace(min(ydat),max(ydat),Ny)
    x,y=numpy.meshgrid(xax,yax)
    samp=array([xdat,ydat])
    kde=stats.kde.gaussian_kde(samp)
    grid_coords = numpy.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)
    z = kde(grid_coords.T)
    z = z.reshape(Nx,Ny)
    asp=xax.ptp()/yax.ptp()
#    if(asp<0.8 or asp > 1.6): asp=1.4
    imshow(z,extent=(xax[0],xax[-1],yax[0],yax[-1]),aspect=asp,origin='lower')
    colorbar()

plot2Dkernel(pos[:,0],pos[:,1],100,100)
if injection and getinjpar(injection,0)<max(pos[:,0]) and getinjpar(injection,0)>min(pos[:,0]) and getinjpar(injection,1)>min(pos[:,1]) and getinjpar(injection,1)<max(pos[:,1]):
        plot(getinjpar(injection,0),getinjpar(injection,1),'go',scalex=False,scaley=False)
xlabel('chirp mass (Msun)')
ylabel('eta')
grid()
myfig.savefig(outdir+'/Meta.png')

if size(unique(pos[:,5]))>1 and size(unique(pos[:,6]))>1:
    myfig.clear()
    plot2Dkernel(pos[:,5],pos[:,6],100,100)
    if injection and getinjpar(injection,5)<max(pos[:,5]) and getinjpar(injection,5)>min(pos[:,5]) and getinjpar(injection,6)>min(pos[:,6]) and getinjpar(injection,6)<max(pos[:,6]):
        plot(getinjpar(injection,5),getinjpar(injection,6),'go',scalex=False,scaley=False)	
    xlabel('RA')
    ylabel('dec')
    grid()
    myfig.savefig(outdir+'/RAdec.png')

myfig.clear()
plot2Dkernel(pos[:,7],pos[:,8],100,100)
if injection and getinjpar(injection,7)<max(pos[:,7]) and getinjpar(injection,7)>min(pos[:,7]) and getinjpar(injection,8)<max(pos[:,8]) and getinjpar(injection,8)>min(pos[:,8]): plot(getinjpar(injection,7),getinjpar(injection,8),'go',scalex=False,scaley=False)
xlabel('psi')
ylabel('iota')
grid()
myfig.savefig(outdir+'/psiiota.png')
myfig.clear()

(m1,m2)=mc2ms(pos[:,0],pos[:,1])
plot2Dkernel(m1,m2,100,100)
if injection and injection.mass1>min(m1) and injection.mass1 < max(m1) and injection.mass2>min(m2) and injection.mass2<max(m2):
    plot(injection.mass1,injection.mass2,'go',scalex=False,scaley=False)
xlabel('mass 1')
ylabel('mass 2')
grid()
myfig.savefig(outdir+'/m1m2.png')
myfig.clear()

plot2Dkernel(m1,pos[:,4],100,100)
if injection and injection.mass1<max(m1) and injection.mass1>min(m1) and getinjpar(injection,4)<max(pos[:,4]) and getinjpar(injection,4)>min(pos[:,4]):
    plot(injection.mass1,injection.distance,'go',scalex=False,scaley=False)
xlabel('m1')
ylabel('Distance (Mpc)')
grid()
myfig.savefig(outdir+'/m1dist.png')
myfig.clear()
plot2Dkernel(m2,pos[:,4],100,100)
xlabel('m2')
ylabel('Distance (Mpc)')
grid()
myfig.savefig(outdir+'/m2dist.png')
myfig.clear()

plot2Dkernel(pos[:,4],pos[:,8],100,100)
if injection and getinjpar(injection,4)>min(pos[:,4]) and getinjpar(injection,4)<max(pos[:,4]) and getinjpar(injection,8)<max(pos[:,8]) and getinjpar(injection,8)>min(pos[:,8]):
    plot(getinjpar(injection,4),getinjpar(injection,8),'go',scalex=False,scaley=False)
xlabel('distance')
ylabel('iota')
grid()
myfig.savefig(outdir+'/Diota.png')
myfig.clear()

paramnames=('Mchirp (Msun)','eta','geocenter time ISCO','phi_c','Distance (Mpc)','RA (rads)','declination (rads)','psi','iota')

for i in range(0,Nd-1):
    for j in range(i+1,Nd-1):
        #print str(i)+','+str(j)+': '+str(size(unique(pos[:,i]))) + ' '+ str(size(unique(pos[:,j])))
        if (size(unique(pos[:,i]))<2 or size(unique(pos[:,j]))<2):   continue
        plot2Dkernel(pos[:,i],pos[:,j],50,50)
        if injection and reduce (lambda a,b: a and b, map(lambda idx: getinjpar(injection,idx)>min(pos[:,idx]) and getinjpar(injection,idx)<max(pos[:,idx]),[i,j])) :
            plot(getinjpar(injection,i),getinjpar(injection,j),'go',scalex=False,scaley=False)
        xlabel(paramnames[i])
        ylabel(paramnames[j])
        grid()
        margdir=outdir+'/2D'
        if not os.path.isdir(margdir+'/'): os.mkdir(margdir)
        myfig.savefig(margdir+'/'+paramnames[i]+'-'+paramnames[j]+'_2Dkernel.png')
        myfig.clear()

htmlfile=open(outdir+'/posplots.html','w')
htmlfile.write('<HTML><HEAD><TITLE>Posterior PDFs</TITLE></HEAD><BODY><h3>'+str(means[2])+' inspnest Posterior PDFs</h3>')
if(Bflag==1): htmlfile.write('<h4>log Bayes Factor: '+str(BayesFactor)+'</h4><br>')
htmlfile.write('signal evidence: '+str(logZ)+'. Information: '+str(H*1.442)+' bits.<br>')
if(Bflag==1): htmlfile.write('deltaLogLmax: '+str(d_sorted[-1,-1])+'<br>')
if(incoflag!=0): htmlfile.write('Odds of coherent vs incoherent: '+str(exp(logZ-Zinco))+'<br>')
if(opts.skyres is not None):
	htmlfile.write('<table border=1><tr><td>Confidence region<td>size (sq. deg)</tr>')
	for (frac,skysize) in skyreses:
		htmlfile.write('<tr><td>%f<td>%f</tr>'%(frac,skysize))
	htmlfile.write('</table>')
htmlfile.write('Produced from '+str(size(pos,0))+' posterior samples, in '+str(size(opts.data,0))+' parallel runs. Taken from '+str(size(d_sorted,0))+' NS samples using '+str(size(opts.data,0)*Nlive)+' live points<br>')
htmlfile.write('<h4>Mean parameter estimates</h4>')
htmlfile.write('<table border=1><tr>')
paramline=reduce(lambda a,b:a+'<td>'+b,paramnames)
htmlfile.write('<td>'+paramline+'<td>logLmax</tr><tr>')
meanline=reduce(lambda a,b:a+'<td>'+b,meanStr)
htmlfile.write('<td>'+meanline+'</tr>')
if injection:
    htmlfile.write('<tr><th colspan=9>Injected values</tr>')
    injline=reduce(lambda a,b:a+'<td>'+b,injvals)
    htmlfile.write('<td>'+injline+'<td></tr>')
htmlfile.write('</table>')
htmlfile.write('<h5>2D Marginal PDFs</h5><br>')
htmlfile.write('<table border=1><tr>')
htmlfile.write('<td width=30%><img width=100% src="m1m2.png"></td>')
htmlfile.write('<td width=30%><img width=100% src="RAdec.png"></td>')
htmlfile.write('<td width=30%><img width=100% src="Meta.png"></td>')
htmlfile.write('</tr><tr><td width=30%><img width=100% src="2D/Mchirp (Msun)-geocenter time ISCO_2Dkernel.png"</td>')
htmlfile.write('<td width=30%><img width=100% src="m1dist.png"></td>')
htmlfile.write('<td width=30%><img width=100% src="m2dist.png"></td>')
htmlfile.write('</table>')
htmlfile.write('<br><a href="2D/">All 2D Marginal PDFs</a><hr><h5>1D marginal posterior PDFs</h5><br>')

for i in [0,1,2,3,4,5,6,7,8]:
    myfig=figure(figsize=(4,3.5),dpi=80)
    hist(pos[:,i],50,normed='true')
    gkde=stats.gaussian_kde(pos[:,i])
    ind=linspace(min(pos[:,i]),max(pos[:,i]),101)
    kdepdf=gkde.evaluate(ind)
    plot(ind,kdepdf,label='density estimate')
    if injection and min(pos[:,i])<getinjpar(injection,i) and max(pos[:,i])>getinjpar(injection,i):   
        plot([getinjpar(injection,i),getinjpar(injection,i)],[0,max(kdepdf)],'r-.',scalex=False,scaley=False)
        print 'i=%i, %f' % (i,getinjpar(injection,i))
    grid()
    xlabel(paramnames[i])
    ylabel('Probability Density')
    myfig.savefig(outdir+'/'+paramnames[i]+ '.png')
    myfig=figure(figsize=(4,3.5),dpi=80)
    plot(pos[:,i],'.')
    if injection and min(pos[:,i])<getinjpar(injection,i) and max(pos[:,i])>getinjpar(injection,i):
	plot([0,len(pos)],[getinjpar(injection,i),getinjpar(injection,i)],'r-.')
    myfig.savefig(outdir+'/'+paramnames[i]+'_samps.png')
    htmlfile.write('<img src="'+paramnames[i]+'.png"><img src="'+paramnames[i]+'_samps.png><br>')

htmlfile.write('<hr><br>Produced using lalapps_inspnest and OddsPostProc.py at '+strftime("%Y-%m-%d %H:%M:%S"))
htmlfile.write('</BODY></HTML>')
htmlfile.close()


# Save posterior samples too...

posfilename=outdir+'/posterior_samples.dat'
posfile=open(posfilename,'w')
for row in pos:
	for i in row:
		posfile.write('%f\t'%(i))
	posfile.write('\n')

posfile.close()

