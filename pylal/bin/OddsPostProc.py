#!/usr/bin/env python

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

parser=OptionParser()
parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
parser.add_option("-N","--Nlive",dest="Nlive",help="number of live points for each of the files")
parser.add_option("-d","--data",dest="data",action="append",help="datafile")
parser.add_option("--inco0",dest="inco0",action="append",help="single-ifo runs for 0th ifo")
parser.add_option("--inco1",dest="inco1",action="append",help="single-ifo runs for 1th ifo")
parser.add_option("--inco2",dest="inco2",action="append",help="single-ifo runs for 2th ifo")
parser.add_option("--inco3",dest="inco3",action="append",help="single-ifo runs for 3th ifo")

(opts,args)=parser.parse_args()

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
print 'Loading ' + opts.data[0]

inco=[]
def loaddata(datalist):
    out=[]
    a=loadtxt(datalist[0])
    for infile in datalist[1:]:
        print 'Loading ' + infile
        tmp=loadtxt(infile)
        out.append(a)
    return out

d=loaddata(opts.data)

if opts.inco0 is not None:
    incoflag=1
    inco.append(loaddata(opts.inco0))
if opts.inco1 is not None:
    incoflag=2
    inco.append(loaddata(opts.inco1))
if opts.inco2 is not None:
    incoflag=3
    inco.append(loaddata(opts.inco2))
if opts.inco3 is not None:
    incoflag=4
    inco.append(loaddata(opts.inco3))

Bfile = opts.data[0]+'_B.txt'
print 'Looking for '+Bfile
if os.access(Bfile,os.R_OK):
    outstat = loadtxt(Bfile)
    NoiseZ = outstat[2]
    Bflag=1
else: Bflag=0

#len=size(d,0)
Nd=size(d[0],1)
#sidx=argsort(d[:,9])
#d=d[sidx,:]
#d[:,0]=exp(d[:,0])
print 'Exponentiated mc'
#maxL = max(d[-1,-1])

Zinco=0
for incox in inco:
    leninc=size(incox,0)
    sidx=argsort(incox[:,-1])
    incox=incox[sidx,:]
    (Zincox,Hinco)=nestZ(incox,Nlive)
    Zinco=Zinco+Zincox

print 'Applying parallelised nested sampling algorithm to samples'

(logZ,H,d_sorted,d_weights)=nestPar(d,Nlive)

d_sorted[:,0]=exp(d_sorted[:,0])
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


if(Bflag==1):
    BayesFactor = logZ - NoiseZ
    print 'log B = '+str(BayesFactor)
    
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
xlabel('chirp mass (Msun)')
ylabel('eta')
grid()
myfig.savefig(outdir+'/Meta.png')

myfig.clear()
plot2Dkernel(pos[:,5],pos[:,6],100,100)
xlabel('RA')
ylabel('dec')
grid()
myfig.savefig(outdir+'/RAdec.png')

myfig.clear()
plot2Dkernel(pos[:,7],pos[:,8],100,100)
xlabel('psi')
ylabel('iota')
grid()
myfig.savefig(outdir+'/psiiota.png')
myfig.clear()

(m1,m2)=mc2ms(pos[:,0],pos[:,1])
plot2Dkernel(m1,m2,100,100)
xlabel('mass 1')
ylabel('mass 2')
grid()
myfig.savefig(outdir+'/m1m2.png')
myfig.clear()

plot2Dkernel(m1,pos[:,4],100,100)
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
xlabel('distance')
ylabel('iota')
grid()
myfig.savefig(outdir+'/Diota.png')
myfig.clear()

paramnames=('Mchirp (Msun)','eta','geocenter time ISCO','phi_c','Distance (Mpc)','RA (rads)','declination (rads)','psi','iota')

for i in range(0,Nd-1):
    for j in range(i+1,Nd-1):
        plot2Dkernel(pos[:,i],pos[:,j],50,50)
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
if(Bflag==1): htmlfile.write('deltaLogLmax: '+str(d_sorted[-1,-1]-NoiseZ)+'<br>')
if(incoflag!=0): htmlfile.write('Odds of coherent vs incoherent: '+str(exp(logZ-Zinco))+'<br>')
htmlfile.write('Produced from '+str(size(pos,0))+' posterior samples, in '+str(size(opts.data,0))+' parallel runs. Taken from '+str(size(d_sorted,0))+' NS samples using '+str(size(opts.data,0)*Nlive)+' live points<br>')
htmlfile.write('<h4>Mean parameter estimates</h4>')
htmlfile.write('<table border=1><tr>')
paramline=reduce(lambda a,b:a+'<td>'+b,paramnames)
htmlfile.write('<td>'+paramline+'<td>logLmax</tr><tr>')
meanline=reduce(lambda a,b:a+'<td>'+b,meanStr)
htmlfile.write('<td>'+meanline+'</tr></table>')
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
    grid()
    xlabel(paramnames[i])
    ylabel('Probability Density')
    myfig.savefig(outdir+'/'+paramnames[i]+ '.png')
    htmlfile.write('<img src="'+paramnames[i]+'.png">')

htmlfile.write('<hr><br>Produced using lalapps_inspnest and OddsPostProc.py at '+strftime("%Y-%m-%d %H:%M:%S"))
htmlfile.write('</BODY></HTML>')
htmlfile.close()
