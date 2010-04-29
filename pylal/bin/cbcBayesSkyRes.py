#!/usr/bin/python

# Load samples from the posterior PDF of an MCMC or nested sampling code
# and produce sky localisation plots and size estimates
# Format for MCMC files should be:
# log(mc) eta tc phase log(dist) RA dec psi iota


#from numpy import *
import scipy
import matplotlib 
matplotlib.use("Agg")
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

def skyhist_cart(skycarts,samples):
	"""
	Histogram the list of samples into bins defined by Cartesian vectors in skycarts
	"""	
	dot=numpy.dot
	N=len(skycarts)
	print 'operating on %d sky points'%(N)
	bins=zeros(N)
	for sample in samples:
		sampcart=pol2cart(sample[RAdim],sample[decdim])
		maxdx=max(xrange(0,N),key=lambda i:dot(sampcart,skycarts[i]))
		bins[maxdx]+=1
	return (skycarts,bins)

def loadDataFile(filename):
	print filename
	infile=open(filename,'r')
	formatstr=infile.readline().lstrip()
	header=formatstr.split()
	llines=[]
	import re
	dec=re.compile(r'[^\d.-]+')
	for line in infile:
		sline=line.split()
		proceed=True
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
	
# Load in the main data
paramnames, pos=loadDataFile(opts.data[0])
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
    if paramnames[parnum]=='mchirp': return inj.mchirp
    if paramnames[parnum]=='eta': return inj.eta
    if paramnames[parnum]=='time': return inj.get_end()
    if paramnames[parnum]=='phi0': return inj.phi0
    if paramnames[parnum]=='dist': return inj.distance
    if paramnames[parnum]=='RA': return inj.longitude
    if paramnames[parnum]=='dec': return inj.latitude
    if paramnames[parnum]=='psi': return inj.polarization
    if paramnames[parnum]=='iota': return inj.inclination
    return None

if injection:
    injpoint=map(lambda a: getinjpar(injection,a),range(0,9))
    injvals=map(str,map(lambda a: getinjpar(injection,a),range(0,9)))
    out=reduce(lambda a,b:a+'||'+b,injvals)
    print 'Injected values:'
    print out

skyreses=[]
if(opts.skyres is not None):
	from pylal import skylocutils
	skypoints=array(skylocutils.gridsky(float(opts.skyres)))
	skycarts=map(lambda s: pol2cart(s[1],s[0]),skypoints)
	(bins,shist)=skyhist_cart(skycarts,pos)
	#(bins,hist)=sky_hist(skypoints,pos)
	# Find the bin of the injection if available
	if injection:
		(injbins,injhist)=skyhist_cart(skycarts,array([injpoint]))
		injbin=injhist.tolist().index(1)
		print 'Found injection in bin %d with co-ordinates %f,%f\n'%(injbin,skypoints[injbin,0],skypoints[injbin,1])
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
		if injection:
			if (injbin==maxpos):
				print 'Injection sky point found at confidence %f'%(frac)
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
		if injection:
                        if (injbin==maxpos):
                                print 'Injection sky point found at confidence %f'%(frac)
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
		if injection:
                        if (injbin==maxpos):
                                print 'Injection sky point found at confidence %f'%(frac)
		#print 'Nbins=%d, thisnum=%d, idx=%d, total=%d, cumul=%f\n'%(Nbins,maxbin,maxpos,len(pos),frac)
        print '%f confidence region: %f square degrees' % (frac,Nbins*float(opts.skyres)*float(opts.skyres))
        skyreses.append((frac,Nbins*float(opts.skyres)*float(opts.skyres)))
	from mpl_toolkits.basemap import Basemap
        myfig=figure()
        clf()
        m=Basemap(projection='moll',lon_0=180.0,lat_0=0.0)
	plx,ply=m(numpy.asarray(toppoints)[::-1,1]*57.296,numpy.asarray(toppoints)[::-1,0]*57.296)
        scatter(plx,ply,s=5,c=numpy.asarray(toppoints)[::-1,2],faceted=False,cmap=matplotlib.cm.jet)
        m.drawmapboundary()
        m.drawparallels(numpy.arange(-90.,120.,45.),labels=[1,0,0,0],labelstyle='+/-')
        # draw parallels
        m.drawmeridians(numpy.arange(0.,360.,90.),labels=[0,0,0,1],labelstyle='+/-')
        # draw meridians
        title("Skymap") # add a title
        colorbar()
        myfig.savefig(outdir+'/skymap.png')
	
myfig=figure(1,figsize=(6,4),dpi=80)
clf()

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

plot2Dkernel(pos[:,paramnames.index('mchirp')],pos[:,paramnames.index('eta')],100,100)
if injection and getinjpar(injection,paramnames.index('mchirp'))<max(pos[:,paramnames.index('mchirp')]) and getinjpar(injection,paramnames.index('mchirp'))>min(pos[:,paramnames.index('mchirp')]) and getinjpar(injection,paramnames.index('eta'))>min(pos[:,paramnames.index('eta')]) and getinjpar(injection,paramnames.index('eta'))<max(pos[:,paramnames.index('eta')]):
        plot(getinjpar(injection,0),getinjpar(injection,1),'go',scalex=False,scaley=False)
xlabel('chirp mass (Msun)')
ylabel('eta')
grid()
myfig.savefig(outdir+'/Meta.png')

if size(unique(pos[:,paramnames.index('RA')]))>1 and size(unique(pos[:,paramnames.index('dec')]))>1:
    myfig.clear()
    plot2Dkernel(pos[:,RAdim],pos[:,decdim],100,100)
    if injection and getinjpar(injection,RAdim)<max(pos[:,RAdim]) and getinjpar(injection,RAdim)>min(pos[:,RAdim]) and getinjpar(injection,decdim)>min(pos[:,decdim]) and getinjpar(injection,decdim)<max(pos[:,decdim]):
        plot(getinjpar(injection,RAdim),getinjpar(injection,decdim),'go',scalex=False,scaley=False)	
    xlabel('RA')
    ylabel('dec')
    grid()
    myfig.savefig(outdir+'/RAdec.png')

myfig.clear()
xdim=paramnames.index('psi')
ydim=paramnames.index('iota')
plot2Dkernel(pos[:,xdim],pos[:,ydim],100,100)
if injection and getinjpar(injection,xdim)<max(pos[:,xdim]) and getinjpar(injection,xdim)>min(pos[:,xdim]) and getinjpar(injection,ydim)<max(pos[:,ydim]) and getinjpar(injection,ydim)>min(pos[:,ydim]): plot(getinjpar(injection,xdim),getinjpar(injection,ydim),'go',scalex=False,scaley=False)
xlabel('psi')
ylabel('iota')
grid()
myfig.savefig(outdir+'/psiiota.png')
myfig.clear()

(m1,m2)=mc2ms(pos[:,paramnames.index('mchirp')],pos[:,paramnames.index('eta')])
plot2Dkernel(m1,m2,100,100)
if injection and injection.mass1>min(m1) and injection.mass1 < max(m1) and injection.mass2>min(m2) and injection.mass2<max(m2):
    plot(injection.mass1,injection.mass2,'go',scalex=False,scaley=False)
xlabel('mass 1')
ylabel('mass 2')
grid()
myfig.savefig(outdir+'/m1m2.png')
myfig.clear()

distdim=paramnames.index('dist')
plot2Dkernel(m1,pos[:,distdim],100,100)
if injection and injection.mass1<max(m1) and injection.mass1>min(m1) and getinjpar(injection,distdim)<max(pos[:,distdim]) and getinjpar(injection,distdim)>min(pos[:,distdim]):
    plot(injection.mass1,injection.distance,'go',scalex=False,scaley=False)
xlabel('m1')
ylabel('Distance (Mpc)')
grid()
myfig.savefig(outdir+'/m1dist.png')
myfig.clear()
plot2Dkernel(m2,pos[:,distdim],100,100)
xlabel('m2')
ylabel('Distance (Mpc)')
grid()
myfig.savefig(outdir+'/m2dist.png')
myfig.clear()

iotadim=paramnames.index('iota')
plot2Dkernel(pos[:,distdim],pos[:,iotadim],100,100)
if injection and getinjpar(injection,distdim)>min(pos[:,distdim]) and getinjpar(injection,distdim)<max(pos[:,distdim]) and getinjpar(injection,iotadim)<max(pos[:,iotadim]) and getinjpar(injection,iotadim)>min(pos[:,iotadim]):
    plot(getinjpar(injection,distdim),getinjpar(injection,iotadim),'go',scalex=False,scaley=False)
xlabel('distance')
ylabel('iota')
grid()
myfig.savefig(outdir+'/Diota.png')
myfig.clear()
#paramnames=('Mchirp (Msun)','eta','geocenter time ISCO','phi_c','Distance (Mpc)','RA (rads)','declination (rads)','psi','iota')

for i in range(0,Nd-1):
    for j in range(i+1,Nd-1):
        print 'Generating %s-%s plot'%(paramnames[i],paramnames[j])
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
if opts.skyres is not None:
        htmlfile.write('<td width=30%><img width=100% src="skymap.png"></td>')
else:
        htmlfile.write('<td width=30%><img width=100% src="m1dist.png:></td>')
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
    htmlfile.write('<img src="'+paramnames[i]+'.png"><img src="'+paramnames[i]+'_samps.png"><br>')

htmlfile.write('<hr><br>Produced using cbcBayesSkyRes.py at '+strftime("%Y-%m-%d %H:%M:%S"))
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

