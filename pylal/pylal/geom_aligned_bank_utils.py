from __future__ import division
import os,sys,math,copy
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab

# FIXME: Use lal/pylals own variables here
mtsun = 4.92529095E-6
LAL_GAMMA = 0.5772156649015328606065120900824024

# This function is taken from Stackoverflow:
# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def determine_eigen_directions(psd,order,f0,f_low,f_upper,delta_f,\
                               return_gs=False,verbose=False,elapsed_time=None):

  if verbose:
    print >>sys.stdout,"Beginning to calculate moments at %d." %(elapsed_time())
  
  moments = get_moments(psd,f0,f_low,f_upper,delta_f)

  if verbose:
    print >>sys.stdout,"Moments calculated, transforming to metric at %d." \
                       %(elapsed_time())

  Js = []
  for i in range(18):
    Js.append(moments['J%d'%(i)])
  Js.append(moments['J%d'%(-1)])

  logJs = []
  for i in range(18):
    logJs.append(moments['log%d'%(i)])
  logJs.append(moments['log%d'%(-1)])


  loglogJs = []
  for i in range(18):
    loglogJs.append(moments['loglog%d'%(i)])
  loglogJs.append(moments['loglog%d'%(-1)])


  logloglogJs = []
  for i in range(18):
    logloglogJs.append(moments['logloglog%d'%(i)])
  logloglogJs.append(moments['logloglog%d'%(-1)])


  loglogloglogJs = []
  for i in range(18):
    loglogloglogJs.append(moments['loglogloglog%d'%(i)])
  loglogloglogJs.append(moments['loglogloglog%d'%(-1)])

  mapping = {}

  if order == 'twoPN':
    maxLen = 4
    gs = numpy.matrix(numpy.zeros(shape=(maxLen,maxLen),dtype=float))
    mapping['Lambda0'] = 0
    mapping['Lambda2'] = 1
    mapping['Lambda3'] = 2
    mapping['Lambda4'] = 3
  elif order == 'threePointFivePN':
    maxLen = 8
    gs = numpy.matrix(numpy.zeros(shape=(maxLen,maxLen),dtype=float))
    mapping['Lambda0'] = 0
    mapping['Lambda2'] = 1
    mapping['Lambda3'] = 2
    mapping['Lambda4'] = 3
    mapping['LogLambda5'] = 4
    mapping['Lambda6'] = 5
    mapping['Lambda7'] = 6
    mapping['LogLambda6'] = 7
  elif order == 'taylorF4_45PN':
    maxLen = 12
    gs = numpy.matrix(numpy.zeros(shape=(maxLen,maxLen),dtype=float))
    mapping['Lambda0'] = 0
    mapping['Lambda2'] = 1
    mapping['Lambda3'] = 2
    mapping['Lambda4'] = 3
    mapping['LogLambda5'] = 4
    mapping['Lambda6'] = 5
    mapping['Lambda7'] = 6
    mapping['LogLambda6'] = 7
#    mapping['Lambda8'] = 8
    mapping['LogLambda8'] = 8
    mapping['LogLogLambda8'] = 9
    mapping['Lambda9'] = 10
    mapping['LogLambda9'] = 11
  else:
    raise BrokenError
 
  for i in range(16):
    for j in range(16):
      # Normal terms
      if mapping.has_key('Lambda%d'%i) and mapping.has_key('Lambda%d'%j):
        gs[mapping['Lambda%d'%i],mapping['Lambda%d'%j]] = 0.5 * (Js[17-i-j] - Js[12-i]*Js[12-j] - (Js[9-i] - Js[4]*Js[12-i]) * (Js[9-j] - Js[4] * Js[12-j])/(Js[1] - Js[4]*Js[4]))
      # Normal,log cross terms
      if mapping.has_key('Lambda%d'%i) and mapping.has_key('LogLambda%d'%j):
        gammaij = logJs[17-i-j] - logJs[12-j] * Js[12-i]
        gamma0i = (Js[9-i] - Js[4] * Js[12-i])
        gamma0j = logJs[9-j] - logJs[12-j] * Js[4]
        gs[mapping['Lambda%d'%i],mapping['LogLambda%d'%j]] = \
            gs[mapping['LogLambda%d'%j],mapping['Lambda%d'%i]] = \
            0.5 * (gammaij - gamma0i*gamma0j/(Js[1] - Js[4]*Js[4]))
      # Log,log terms
      if mapping.has_key('LogLambda%d'%i) and mapping.has_key('LogLambda%d'%j):
        gammaij = loglogJs[17-i-j] - logJs[12-j] * logJs[12-i]
        gamma0i = (logJs[9-i] - Js[4] * logJs[12-i])
        gamma0j = logJs[9-j] - logJs[12-j] * Js[4]
        gs[mapping['LogLambda%d'%i],mapping['LogLambda%d'%j]] = \
            0.5 * (gammaij - gamma0i*gamma0j/(Js[1] - Js[4]*Js[4]))
      # Normal,loglog cross terms
      if mapping.has_key('Lambda%d'%i) and mapping.has_key('LogLogLambda%d'%j):
        gammaij = loglogJs[17-i-j] - loglogJs[12-j] * Js[12-i]
        gamma0i = (Js[9-i] - Js[4] * Js[12-i])
        gamma0j = loglogJs[9-j] - loglogJs[12-j] * Js[4]
        gs[mapping['Lambda%d'%i],mapping['LogLogLambda%d'%j]] = \
            gs[mapping['LogLogLambda%d'%j],mapping['Lambda%d'%i]] = \
            0.5 * (gammaij - gamma0i*gamma0j/(Js[1] - Js[4]*Js[4]))
      # log,loglog cross terms
      if mapping.has_key('LogLambda%d'%i) and mapping.has_key('LogLogLambda%d'%j):
        gammaij = logloglogJs[17-i-j] - loglogJs[12-j] * logJs[12-i]
        gamma0i = (logJs[9-i] - Js[4] * logJs[12-i])
        gamma0j = loglogJs[9-j] - loglogJs[12-j] * Js[4]
        gs[mapping['LogLambda%d'%i],mapping['LogLogLambda%d'%j]] = \
            gs[mapping['LogLogLambda%d'%j],mapping['LogLambda%d'%i]] = \
            0.5 * (gammaij - gamma0i*gamma0j/(Js[1] - Js[4]*Js[4]))
      # Loglog,loglog terms
      if mapping.has_key('LogLogLambda%d'%i) and mapping.has_key('LogLogLambda%d'%j):
        gammaij = loglogloglogJs[17-i-j] - loglogJs[12-j] * loglogJs[12-i]
        gamma0i = (loglogJs[9-i] - Js[4] * loglogJs[12-i])
        gamma0j = loglogJs[9-j] - loglogJs[12-j] * Js[4]
        gs[mapping['LogLogLambda%d'%i],mapping['LogLogLambda%d'%j]] = \
            0.5 * (gammaij - gamma0i*gamma0j/(Js[1] - Js[4]*Js[4]))

  evals,evecs = numpy.linalg.eig(gs)
  gs = numpy.matrix(gs)

  if verbose:
    print >>sys.stdout,"Metric and eigenvalues calculated at %d." \
                       %(elapsed_time())

  for i in range(len(evals)):
    if evals[i] < 0:
      print "WARNING: Negative eigenvalue %e. Setting as positive." %(evals[i])
      evals[i] = -evals[i]
  if return_gs:
    return evals,evecs,gs

  return evals,evecs

def get_moments(psd_file,f0,f_low,f_high,deltaF):
  psd = numpy.loadtxt(psd_file)
  psd_f = psd[:,0]
  psd_amp = psd[:,1]
  psd_amp = psd_amp * psd_amp
  new_f,new_amp = interpolate_psd(psd_f,psd_amp,deltaF)

  # Need I7
  funct = lambda x: 1
  I7 = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)

  # Do all the J moments
  moments = {}
  for i in range(-1,18):
    funct = lambda x: x**((-i+7)/3.)
    moments['J%d' %(i)] = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)/I7

  # Do the logx multiplied by some power terms
  for i in range(-1,18):
    funct = lambda x: (numpy.log(x**(1./3.))) * x**((-i+7)/3.)
    moments['log%d' %(i)] = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)/I7

  # Do the loglog term
  for i in range(-1,18):
    funct = lambda x: (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * x**((-i+7)/3.)
    moments['loglog%d' %(i)] = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)/I7

  # Do the logloglog term
  for i in range(-1,18):
    funct = lambda x: (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * x**((-i+7)/3.)
    moments['logloglog%d' %(i)] = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)/I7

  # Do the logloglog term
  for i in range(-1,18):
    funct = lambda x: (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * (numpy.log(x**(1./3.))) * x**((-i+7)/3.)
    moments['loglogloglog%d' %(i)] = calculate_moment(new_f,new_amp,f_low,f_high,f0,funct)/I7

  return moments

def interpolate_psd(psd_f,psd_amp,deltaF):
  new_psd_f = []
  new_psd_amp = []
  fcurr = psd_f[0]

  for i in range(len(psd_f) - 1):
    f_low = psd_f[i]
    f_high = psd_f[i+1]
    amp_low = psd_amp[i]
    amp_high = psd_amp[i+1]
    while(1):
      if fcurr > f_high:
        break
      new_psd_f.append(fcurr)
      gradient = (amp_high - amp_low) / (f_high - f_low)
      fDiff = fcurr - f_low
      new_psd_amp.append(amp_low + fDiff * gradient)
      fcurr = fcurr + deltaF
  return numpy.asarray(new_psd_f),numpy.asarray(new_psd_amp)


def calculate_moment(psd_f,psd_amp,fmin,fmax,f0,funct):
  # Must ensure deltaF in psd_f is constant
  psd_x = psd_f / f0
  deltax = psd_x[1] - psd_x[0]

  comps = (psd_x)**(-7./3.) * funct(psd_x) * deltax/ psd_amp
  comps_red = comps[numpy.logical_and(psd_f > fmin, psd_f < fmax)]
  return comps_red.sum()

def estimate_mass_range_slimline(numiters,order,evals,evecs,maxmass1,minmass1,maxmass2,minmass2,maxspin,f0,covary=True,maxBHspin=None,evecsCV=None):
  out = []
  valsF = get_random_mass_slimline(numiters,minmass1,maxmass1,minmass2,maxmass2,maxspin,maxBHspin = maxBHspin,return_spins=True)
  valsF = numpy.array(valsF)
  mass = valsF[0]
  eta = valsF[1]
  beta = valsF[2]
  sigma = valsF[3]
  gamma = valsF[4]
  chis = 0.5*(valsF[5] + valsF[6])
  if covary:
    lambdas = get_cov_params(mass,eta,beta,sigma,gamma,chis,f0,evecs,evals,evecsCV,order)
  else:
    lambdas = get_conv_params(mass,eta,beta,sigma,gamma,chis,f0,evecs,evals,order)

  return numpy.array(lambdas)

def get_random_mass_slimline(N,minmass1,maxmass1,minmass2,maxmass2,maxspin,maxBHspin = None,return_spins=False,qm_scalar_fac=1):
  # WARNING: We expect mass1 > mass2 ALWAYS
  minmass = minmass1 + minmass2
  maxmass = maxmass1 + maxmass2
  mincompmass = minmass2
  maxcompmass = maxmass1

  mass = numpy.random.random(N) * (minmass**(-5./3.)-maxmass**(-5./3.)) + maxmass**(-5./3.)
  mass = mass**(-3./5.)
  maxmass2 = numpy.minimum(mass/2.,maxmass2)
  minmass1 = numpy.maximum(minmass1,mass/2.)
  mineta = numpy.maximum(mincompmass * (mass-mincompmass)/(mass*mass), maxcompmass*(mass-maxcompmass)/(mass*mass))
  maxeta = numpy.minimum(0.25,maxmass2 * (mass - maxmass2) / (mass*mass))
  maxeta = numpy.minimum(maxeta,minmass1 * (mass - minmass1) / (mass*mass))
  if (maxeta < mineta).any():
    print "WARNING: Max eta is smaller than min eta!!"
  eta = numpy.random.random(N) * (maxeta - mineta) + mineta
  diff = (mass*mass * (1-4*eta))**0.5
  mass1 = (mass + diff)/2.
  mass2 = (mass - diff)/2.
  if (mass1 > maxmass1).any() or (mass1 < minmass1).any():
    print "WARNING: Mass1 outside of mass range"
  if (mass2 > maxmass2).any() or (mass2 < minmass2).any():
    print "WARNING: Mass1 outside of mass range"
  if maxspin > 0:
    mspin = numpy.zeros(len(mass1))
    mspin += maxspin
    if maxBHspin:
      mspin[mass1 > 3] = maxBHspin
    spin1z = numpy.random.random(N) * mspin*2 - mspin
    mspin = numpy.zeros(len(mass2))
    mspin += maxspin
    if maxBHspin:
      mspin[mass2 > 3] = maxBHspin
    spin2z = numpy.random.random(N) * mspin*2 - mspin

    spinspin = spin1z*spin2z
  else:
    spinspin = 0
    spin1z = 0
    spin2z = 0

  chiS = 0.5 * (spin1z + spin2z)
  chiA = 0.5 * (spin1z - spin2z)
  delta = (mass1 - mass2) / (mass1 + mass2)
 
  beta = (113. / 12. - 19./3. * eta) * chiS
  beta += 113. / 12. * delta * chiA
  sigma = eta / 48. * (474 * spinspin)
  gamma = numpy.zeros(len(sigma))
  sigma += (1 - 2*eta) * (81./16. * (chiS*chiS + chiA*chiA))
  sigma += delta * (81. / 8. * (chiS*chiA))
  gamma = (732985./2268. - 24260./81.*eta - 340./9.*eta*eta)*chiS
  gamma += (732985. / 2268. + 140./9. * eta) * delta * chiA

  if return_spins:
    return mass,eta,beta,sigma,gamma,spin1z,spin2z
  else:
    return mass,eta,beta,sigma,gamma

def rotate_vector(evecs,old_vector,rescale_factor,index,length):
  temp = 0
  for i in range(length):
    temp += evecs[i,index] * old_vector[i]
  temp *= rescale_factor
  return temp

def rotate_vector_inv(evecs,old_vector,rescale_factor,index,length):
  temp = 0
  for i in range(length):
    temp += evecs[index,i] * old_vector[i]
  temp *= rescale_factor
  return temp

def get_conv_params(totmass,eta,beta,sigma,gamma,chis,f0,evecs,evals,order):

  lambdas = get_chirp_params(totmass,eta,beta,sigma,gamma,chis,f0,order)

  if order == 'twoPN':
    lam0 = rotate_vector(evecs,lambdas,math.sqrt(evals[0]),0,4)
    lam2 = rotate_vector(evecs,lambdas,math.sqrt(evals[1]),1,4)
    lam3 = rotate_vector(evecs,lambdas,math.sqrt(evals[2]),2,4)
    lam4 = rotate_vector(evecs,lambdas,math.sqrt(evals[3]),3,4)
    return lam0,lam2,lam3,lam4
  elif order[0:16] == 'threePointFivePN' or order[0:13] == 'taylorF4_35PN':
    lams = []
    for i in range(8):
      lams.append(rotate_vector(evecs,lambdas,math.sqrt(evals[i]),i,8))
    return lams
  elif order[0:13] == 'taylorF4_45PN':
    lams = []
    for i in range(12):
      lams.append(rotate_vector(evecs,lambdas,math.sqrt(evals[i]),i,12))
    return lams
  else:
    raise BrokenError

def get_cov_params(totmass,eta,beta,sigma,gamma,chis,f0,evecs,evals,evecsCV,order):
  mus = get_conv_params(totmass,eta,beta,sigma,gamma,chis,f0,evecs,evals,order)
  xis = get_covaried_params(mus,evecsCV)
  return xis

def get_covaried_params(lambdas,evecsCV):
  length = len(evecsCV)
  lams = []
  for i in range(length):
    lams.append(rotate_vector(evecsCV,lambdas,1.,i,length))
  return lams

def get_chirp_params(totmass,eta,beta,sigma,gamma,chis,f0,order):
  # Convert mass to seconds
  totmass = totmass * mtsun
  pi = numpy.pi
  lambda0 = 3 / (128 * eta * (pi * totmass * f0)**(5/3))
  lambda2 = 5 / (96 * pi * eta * totmass * f0) * (743/336 + 11/4 * eta)
  lambda3 = (-3 * pi**(1/3))/(8 * eta * (totmass*f0)**(2/3)) * (1 - beta/ (4 * pi))
  lambda4 = 15 / (64 * eta * (pi * totmass * f0)**(1/3)) * (3058673/1016064 + 5429/1008 * eta + 617/144 * eta**2 - sigma)
  if order == 'twoPN':
    return lambda0,lambda2,lambda3,lambda4
  elif order[0:16] == 'threePointFivePN' or order[0:8] == 'taylorF4':
    lambda5 = 3. * (38645.*pi/756. - 65.*pi*eta/9. - gamma)
    lambda5 = lambda5 * (3./(128.*eta))
    lambda6 = 11583231236531./4694215680. - (640.*pi*pi)/3. - (6848.*LAL_GAMMA)/21.
    lambda6 -= (6848./21.)  * numpy.log(4 * (pi * totmass * f0)**(1./3.))
    lambda6 += (-15737765635/3048192. + 2255.*pi*pi/12.)*eta
    lambda6 += (76055.*eta*eta)/1728. - (127825.*eta*eta*eta)/1296.;
    lambda6 = lambda6
    lambda6 = lambda6 * 3./(128.*eta) * (pi * totmass * f0)**(1/3.)
    lambda7 = (77096675.*pi)/254016. + (378515.*pi*eta)/1512. - (74045.*pi*eta*eta)/756.
    lambda7 = lambda7
    lambda7 = lambda7 * 3./(128.*eta) * (pi * totmass * f0)**(2/3.)
    lambda6log =  -( 6848./21)
    lambda6log = lambda6log * 3./(128.*eta) * (pi * totmass * f0)**(1/3.)
    if order[0:16] == 'threePointFivePN':
      return lambda0,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda6log
    elif order[0:13] == 'taylorF4_45PN' or order[0:13] == 'taylorF4_35PN':
      # Add 3PN spin corr term
      lambda6spin = 502.6548245743669 * beta + 88.45238095238095 * sigma + \
                   (110. * eta * sigma) - 20. * beta * beta
      lambda6spin = lambda6spin * 3./(128.*eta) * (pi * totmass * f0)**(1/3.)
      lambda6 += lambda6spin
      # Add 3.5PN spin corr term
      lambda7spin = -510.0603994773098*beta - 368.01525846326734*beta*eta + \
          1944.363555525455*chis*eta - 502.6548245743669*sigma + \
          40.*beta*sigma + 241.47615535889872*beta*eta*eta + \
          2961.654024441635*chis*eta*eta + 676.0619469026549*chis*eta*eta*eta
      lambda7spin = lambda7spin * 3./(128.*eta) * (pi * totmass * f0)**(2/3.)
      lambda7 += lambda7spin
      # Add 4PN non spin term (this has log and loglog terms)
      lambda8 = 342.6916926002232 + 2869.024558661873*eta - \
          3773.659169914512*eta*eta + 172.0609149438239*eta*eta*eta - \
          24.284336419753085*eta*eta*eta*eta
      lambda8log = -1028.0750778006693 - 8607.073675985623*eta + \
          11320.977509743536*eta*eta - 516.1827448314717*eta*eta*eta + \
          72.85300925925927*eta*eta*eta*eta
      lambda8loglog = 480.7316704459562 + 597.8412698412699*eta
      # And the 4PN spin terms
      lambda8spin = 936.7471880419974*beta - 311.03929625364435*beta*eta - \
          2455.4171568883194*chis*eta + 195.39588893571195*beta*chis*eta + \
          48.491201779065534*sigma + 101.92901234567901*eta*sigma - \
          58.81315259633844*beta*beta + \
          8.918387413962636*eta*beta*beta - 686.5167663065837*chis*eta*eta \
          + 54.631268436578175*beta*chis*eta*eta + \
          71.69753086419753*sigma*eta*eta - \
          4.444444444444445*sigma*sigma
      lambda8logspin = -2810.241564125992*beta + 933.117888760933*beta*eta + \
          7366.251470664957*chis*eta - 586.1876668071359*beta*chis*eta - \
          145.4736053371966*sigma - 305.78703703703707*eta*sigma \
          + 176.4394577890153*beta*beta - \
          26.755162241887906*eta*beta*beta + \
          2059.5502989197507*chis*eta*eta - \
          163.89380530973452*beta*chis*eta*eta - \
          215.0925925925926*sigma*eta*eta + \
          13.333333333333334*sigma*sigma
      # Construct the combined 4PN terms
      lambda8 = (lambda8 + lambda8spin) * 3./(128.*eta) * \
          (pi * totmass * f0)**(3./3.)
      lambda8log = (lambda8log + lambda8logspin) * 3./(128.*eta) * \
          (pi * totmass * f0)**(3./3.)
      lambda8loglog = lambda8loglog * 3./(128.*eta) * \
          (pi * totmass * f0)**(3./3.)
      # And now the 4.5PN non-spin term
      lambda9 = 20021.24571514093 - 42141.161261993766*eta - \
          4047.211701119762*eta*eta - 2683.4848475303043*eta*eta*eta
      lambda9log = -4097.833617482457
      # And now the 4.5PN spin term
      lambda9spin = 2105.9471080635244*beta + 3909.271818583914*beta*eta - \
          2398.354686411564*chis*eta + 1278.4225104920606*sigma - \
          198.6688790560472*beta*sigma + 589.0486225480862*eta*sigma - \
          62.43362831858406*beta*eta*sigma + \
          439.6407501053519*chis*eta*sigma - 376.99111843077515*beta*beta + \
          10.*beta*beta*beta - 202.1451909795383*beta*eta*eta - \
          5711.929102446965*chis*eta*eta + \
          122.9203539823009*chis*sigma*eta*eta - \
          493.00738145963066*beta*eta*eta*eta - \
          4955.659484448894*chis*eta*eta*eta - \
          991.4721607669617*chis*eta*eta*eta*eta
      lambda9logspin = 326.0952380952381*beta
      lambda9 = (lambda9 + lambda9spin) * 3./(128.*eta) * \
          (pi * totmass * f0)**(4./3.)
      lambda9log = (lambda9log + lambda9logspin) * 3./(128.*eta) * \
          (pi * totmass * f0)**(4./3.)
      if order[0:13] == 'taylorF4_45PN':
        return lambda0,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda6log,lambda8log,lambda8loglog,lambda9,lambda9log
      elif order[0:13] == 'taylorF4_35PN':
        return lambda0,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda6log
    else:
      raise BrokenError
  else:
    raise BrokenError

def make_plots(a,b,c,d,aname,bname,cname,dname,paper_plots=False):
  if paper_plots:
    paper_plot()
  if not os.path.isdir('plots'):
    os.makedirs('plots')
  vals = [a,b,c,d]
  names = [aname,bname,cname,dname]
  for i in range(4):
    for j in range(i+1,4):
      pylab.plot(vals[i],vals[j],'b.')
      pylab.xlabel(names[i])
      pylab.ylabel(names[j])
#      if i == 0:
#       pylab.xlim([-700,0])
#       if j == 1:
#         pylab.ylim([-3.5,0.5])
#       if j == 2:
#         pylab.ylim([-2,1])
#       if j == 3:
#         pylab.ylim([0.,0.1])
      pylab.savefig('plots/%s_vs_%s.png' % (names[i],names[j]))
      pylab.clf()
 
  # 3D plots
#  for i in range(4):
#    for j in range(i+1,4):
#      for k in range(j+1,4):
#        fig = pylab.figure()
#        ax = Axes3D(fig)
#        ax.plot(vals[i],vals[j],'b.',zs=vals[k])
#        Axes3D.set_xlabel(ax,names[i])
#        Axes3D.set_ylabel(ax,names[j])
#        Axes3D.set_zlabel(ax,names[k])
#        pylab.savefig('plots/%s_vs_%s_vs_%s.png' % (names[i],names[j],names[k]))
#        pylab.clf()

def generate_hexagonal_lattice(maxv1,minv1,maxv2,minv2,mindist):
  # Place first point
  v1s = []
  v2s = []
  initPoint = [minv1,minv2]
  # Place first line
  initLine = [initPoint]
  tmpv1 = minv1
  while (tmpv1 < maxv1):
    tmpv1 = tmpv1 + (3 * mindist)**(0.5)
    initLine.append([tmpv1,minv2])
    v1s.append(tmpv1)
    v2s.append(minv2)
  initLine = numpy.array(initLine)
  initLine2 = copy.deepcopy(initLine)
  initLine2[:,0] += 0.5 * (3*mindist)**0.5
  initLine2[:,1] += 1.5 * (mindist)**0.5
  for i in xrange(len(initLine2)):
    v1s.append(initLine2[i,0])
    v2s.append(initLine2[i,1])
  tmpv2_1 = initLine[0,1]
  tmpv2_2 = initLine2[0,1]
  while tmpv2_1 < maxv2 and tmpv2_2 < maxv2:
    tmpv2_1 = tmpv2_1 + 3.0 * (mindist)**0.5
    tmpv2_2 = tmpv2_2 + 3.0 * (mindist)**0.5 
    initLine[:,1] = tmpv2_1
    initLine2[:,1] = tmpv2_2
    for i in xrange(len(initLine)):
      v1s.append(initLine[i,0])
      v2s.append(initLine[i,1])
    for i in xrange(len(initLine2)):
      v1s.append(initLine2[i,0])
      v2s.append(initLine2[i,1])
  v1s = numpy.array(v1s)
  v2s = numpy.array(v2s)
  return v1s,v2s

def get_physical_covaried_masses(xis,bestMasses,bestXis,f0,temp_number,req_match,order,evecs,evals,evecsCV,maxmass1,minmass1,maxmass2,minmass2,maxNSspin,maxBHspin,return_smaller=False,nsbh_flag=False):
  # TUNABLE PARAMETERS GO HERE!
  origScaleFactor = 1

  # Set up
  xi_size = len(xis)
  scaleFactor = origScaleFactor
  bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)
  count = 0
  unFixedCount = 0
  currDist = 100000000000000000
  while(1):
#    if not (count % 100):
#      print '\rTemplate %d  Distance %e   Count %d  Scale Factor %d       ' %(temp_number,currDist,count,scaleFactor),
  # If we are a long way away we use larger jumps
    if count:
      if currDist > 1 and scaleFactor == origScaleFactor:
        scaleFactor = origScaleFactor*10
    chirpmass,totmass,eta,spin1z,spin2z,diff,mass1,mass2,beta,sigma,gamma,chis,new_xis = get_mass_distribution(bestChirpmass,bestMasses[1],bestMasses[2],bestMasses[3],scaleFactor,order,evecs,evals,evecsCV,maxmass1,minmass1,minmass2,maxmass2,maxNSspin,maxBHspin,f0,nsbh_flag = nsbh_flag)
    cDist = (new_xis[0] - xis[0])**2
    for j in range(1,xi_size):
      cDist += (new_xis[j] - xis[j])**2
    if (cDist.min() < req_match):
      idx = cDist.argmin()
      scaleFactor = origScaleFactor
      return mass1[idx],mass2[idx],spin1z[idx],spin2z[idx],count,cDist.min(),new_xis[0][idx],new_xis[1][idx],new_xis[2][idx],new_xis[3][idx]
    if (cDist.min() < currDist):
      idx = cDist.argmin()
      bestMasses[0] = totmass[idx]
      bestMasses[1] = eta[idx]
      bestMasses[2] = spin1z[idx]
      bestMasses[3] = spin2z[idx]
      bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)
      currDist = cDist.min()
      unFixedCount = 0
      scaleFactor = origScaleFactor
    count += 1
    unFixedCount += 1
    if unFixedCount > 5000:
      if return_smaller:
        diff = (bestMasses[0]*bestMasses[0] * (1-4*bestMasses[1]))**0.5
        mass1 = (bestMasses[0] + diff)/2.
        mass2 = (bestMasses[0] - diff)/2.
        return mass1,mass2,bestMasses[2],bestMasses[3],count,currDist,new_xis[0][0],new_xis[1][0],new_xis[2][0],new_xis[3][0]
      # Give up
      else:
        raise BrokenError
    if not unFixedCount % 100:
      scaleFactor *= 2
    if scaleFactor > 64:
      scaleFactor = 1
  # Shouldn't be here!
  raise BrokenError

def get_mass_distribution(bestChirpmass,bestEta,bestSpin1z,bestSpin2z,scaleFactor,order,evecs,evals,evecsCV,maxmass1,minmass1,minmass2,maxmass2,maxNSspin,maxBHspin,f0,nsbh_flag = False):
  chirpmass = bestChirpmass * ( 1 - (numpy.random.random(100) - 0.5) * 0.0001 * scaleFactor)
  chirpmass[0] = bestChirpmass
  eta = bestEta * ( 1 - (numpy.random.random(100) - 0.5) * 0.01 * scaleFactor)
  eta[0] = bestEta
#   eta = bestMasses[2] + ( (numpy.random.random(100) - 0.5) * 0.01 * scaleFactor)
  eta[eta > 0.25] = 0.25
  eta[eta < 0.0001] = 0.0001
  totmass = chirpmass / (eta**(3./5.))
  spin1z = bestSpin1z + ( (numpy.random.random(100) - 0.5) * 0.01 * scaleFactor)
  spin1z[0] = bestSpin1z
  spin2z = bestSpin2z + ( (numpy.random.random(100) - 0.5) * 0.01 * scaleFactor)
  spin2z[0] = bestSpin2z
  beta,sigma,gamma,chis = get_beta_sigma_from_aligned_spins(totmass,eta,spin1z,spin2z)

  diff = (totmass*totmass * (1-4*eta))**0.5
  mass1 = (totmass + diff)/2.
  mass2 = (totmass - diff)/2.

  numploga1 = numpy.logical_and(abs(spin1z) < maxBHspin,mass1 > 2.99)
  if nsbh_flag:
    numploga = numploga1
  else:
    numploga2 = numpy.logical_and(abs(spin1z) < maxNSspin,mass1 < 3.01)
    numploga = numpy.logical_or(numploga1,numploga2)
  numplogb2 = numpy.logical_and(abs(spin2z) < maxNSspin,mass2 < 3.01)
  if nsbh_flag:
    numplogb = numplogb2
  else:
    numplogb1 = numpy.logical_and(abs(spin2z) < maxBHspin,mass2 > 2.99)
    numplogb = numpy.logical_or(numplogb1,numplogb2)
  numplog1 = numpy.logical_and(numploga,numplogb)
  numplog = numpy.logical_not(numplog1)
  beta[numplog] = 0
  sigma[numplog] = 0
  gamma[numplog] = 0
  chis[numplog] = 0
  spin1z[numplog] = 0
  spin2z[numplog] = 0

  totmass[mass1 < minmass1] = 0.0001
  totmass[mass1 > maxmass1] = 0.0001
  totmass[mass2 < minmass2] = 0.0001
  totmass[mass2 > maxmass2] = 0.0001

  new_xis = get_cov_params(totmass,eta,beta,sigma,gamma,chis,f0,evecs,evals,evecsCV,order)
  return chirpmass,totmass,eta,spin1z,spin2z,diff,mass1,mass2,beta,sigma,gamma,chis,new_xis

def stack_xi_direction_brute(xis,bestMasses,bestXis,f0,temp_number,direction_num,req_match,order,evecs,evals,evecsCV,maxmass1,minmass1,maxmass2,minmass2,maxNSspin,maxBHspin,nsbh_flag=False):
  # TUNING PARAMETERS GO HERE
  origScaleFactor = 0.8
  numTestPoints = 3000

  # Setup
  xi_size = len(xis)
  origMasses = copy.deepcopy(bestMasses)
  bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)
  count = 0
  unFixedCount = 0
  xi3min = 10000000000
  xi3max = -100000000000

  for i in range(numTestPoints):
    # Evaluate upper extent of xi3
    scaleFactor = origScaleFactor
    chirpmass,totmass,eta,spin1z,spin2z,diff,mass1,mass2,beta,sigma,gamma,chis,new_xis = get_mass_distribution(bestChirpmass,bestMasses[1],bestMasses[2],bestMasses[3],scaleFactor,order,evecs,evals,evecsCV,maxmass1,minmass1,minmass2,maxmass2,maxNSspin,maxBHspin,f0,nsbh_flag = nsbh_flag)
    cDist = (new_xis[0] - xis[0])**2
    for j in range(1,xi_size):
      cDist += (new_xis[j] - xis[j])**2
    redCDist = cDist[cDist < req_match]
    redXis = (new_xis[direction_num])[cDist < req_match]

    if len(redCDist):
      new_xis[direction_num][cDist > req_match] = -10000000
      maxXi3 = (new_xis[direction_num]).max()
      idx = (new_xis[direction_num]).argmax()
      if maxXi3 > xi3max:
        xi3max = maxXi3
        bestMasses[0] = totmass[idx]
        bestMasses[1] = eta[idx]
        bestMasses[2] = spin1z[idx]
        bestMasses[3] = spin2z[idx]
        m1 = mass1[idx]
        m2 = mass2[idx]
        bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)

  bestMasses = copy.deepcopy(origMasses)
  bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)

  for i in range(numTestPoints):
    # Evaluate lower extent of xi3
    scaleFactor = origScaleFactor
    chirpmass,totmass,eta,spin1z,spin2z,diff,mass1,mass2,beta,sigma,gamma,chis,new_xis = get_mass_distribution(bestChirpmass,bestMasses[1],bestMasses[2],bestMasses[3],scaleFactor,order,evecs,evals,evecsCV,maxmass1,minmass1,minmass2,maxmass2,maxNSspin,maxBHspin,f0,nsbh_flag = nsbh_flag)

    cDist = (new_xis[0] - xis[0])**2
    for j in range(1,xi_size):
      cDist += (new_xis[j] - xis[j])**2
    redCDist = cDist[cDist < req_match]
    redXis = (new_xis[direction_num])[cDist < req_match]

    if len(redCDist):
      new_xis[direction_num][cDist > req_match] = 10000000
      maxXi3 = (new_xis[direction_num]).min()
      idx = (new_xis[direction_num]).argmin()
      if maxXi3 < xi3min:
        xi3min = maxXi3
        bestMasses[0] = totmass[idx]
        bestMasses[1] = eta[idx]
        bestMasses[2] = spin1z[idx]
        bestMasses[3] = spin2z[idx]
        m1 = mass1[idx]
        m2 = mass2[idx]
        bestChirpmass = bestMasses[0] * (bestMasses[1])**(3./5.)

  return xi3min,xi3max,xi3max-xi3min

def get_beta_sigma_from_aligned_spins(mass,eta,spin1z,spin2z):
  diff = (mass*mass * (1-4*eta))**0.5
  mass1 = (mass + diff)/2.
  mass2 = (mass - diff)/2.
  chiS = 0.5 * (spin1z + spin2z)
  chiA = 0.5 * (spin1z - spin2z)
  delta = (mass1 - mass2) / (mass1 + mass2)
  spinspin = spin1z*spin2z

  beta = (113. / 12. - 19./3. * eta) * chiS
  beta += 113. / 12. * delta * chiA
  sigma = eta / 48. * (474 * spinspin)
  sigma += (1 - 2*eta) * (81./16. * (chiS*chiS + chiA*chiA))
  sigma += delta * (81. / 8. * (chiS*chiA))
  gamma = (732985./2268. - 24260./81.*eta - 340./9.*eta*eta)*chiS
  gamma += (732985. / 2268. + 140./9. * eta) * delta * chiA
  return beta,sigma,gamma,chiS

