#!/usr/bin/python

import git_version

__author__ = "Larry Price <larry.price@ligo.org> and Patrick Brady <patrick.brady@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import sys
import os
import operator
from math import sqrt, sin, cos
import gzip 

from pylal import date
from pylal import CoincInspiralUtils, SnglInspiralUtils, SimInspiralUtils
from pylal.xlal import tools, inject
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS 

import glue.iterutils
from glue.ligolw import utils, table as tab, lsctables




##############################################################################
#
#          global variables
#
##############################################################################

pi = 3.14159265358979324

#set the detector locations
detector_locations = {}
detector_locations["L1"] =tools.cached_detector["LLO_4k"].location   
detector_locations["H1"] =tools.cached_detector["LHO_4k"].location 
detector_locations["V1"] =tools.cached_detector["VIRGO"].location

#set the detector responses
detector_responses = {}
detector_responses["L1"] = tools.cached_detector["LLO_4k"].response
detector_responses["H1"] = tools.cached_detector["LHO_4k"].response
detector_responses["V1"] = tools.cached_detector["VIRGO"].response

##############################################################################
#
#          function definitions
#
##############################################################################

def get_delta_t_rss(pt,coinc,reference_frequency=None):
  """
  returns the rss timing error for a particular location in
  the sky (longitude,latitude)
  """
  latitude,longitude = pt
  earth_center = (0.0,0.0,0.0)
  tref = {}
  tgeo={}
  for ifo in coinc.ifo_list:
    
    if reference_frequency:
      tFromRefFreq = get_signal_duration(ifo,coinc,reference_frequency)
      tref[ifo] = LIGOTimeGPS(int(tFromRefFreq), 1.e9*(tFromRefFreq-int(tFromRefFreq)))
    else:
      tref[ifo] = 0.0   
         
    #compute the geocentric time from each trigger
    tgeo[ifo] = coinc.gps[ifo] - tref[ifo] - \
                LIGOTimeGPS(0,1.0e9*date.XLALArrivalTimeDiff(detector_locations[ifo],\
                earth_center,longitude,latitude,coinc.gps[ifo]))  
      
  #compute differences in these geocentric times
  time={}
  delta_t_rms = 0.0
  for ifos in coinc.ifo_coincs:
    time[ifos[0]+ifos[1]] = 1.0e-9*date.XLALGPSToINT8NS(tgeo[ifos[0]] - tgeo[ifos[1]])
    delta_t_rms += time[ifos[0]+ifos[1]] * time[ifos[0]+ifos[1]]
        
  return sqrt(delta_t_rms)
  
def get_signal_duration(ifo,coinc,frequency):
  """
  determine the amount by which the coalescence time must be translated to get the reference time
  """
  M = coinc.mass1[ifo]+coinc.mass2[ifo]
  mu = coinc.mass1[ifo]*coinc.mass2[ifo]/M
  eta = mu/M
  chirpM = pow(mu*mu*mu*M*M,1./5.)
  M = M*4.92549095e-6
  mu = mu*4.92549095e-6
  chirpM = chirpM*4.92549095e-6
  tau0 = 5./256. * pow(chirpM,-5./3.) * pow(pi*frequency,-8./3.)
  tau1 = 5./(192.*mu*pi*pi*frequency*frequency) * (743./336. + 11./4.*eta)
  tau1p5 = 1./(8.*mu) * pow(M/(pi*pi*pow(frequency,5.)),1./3.)
  tau2 = 5./(128.*mu) * pow(M/(pi*pi*frequency*frequency),2./3.)\
         *(3058673./1016064. + 5429./1008.*eta + 617./144.*eta*eta)        
  duration = tau0 + tau1 - tau1p5 + tau2
    
  return duration
  
def get_delta_D_rss(pt,coinc):
  """
  compute the rms difference in the ratio of the difference of the squares of Deff to
  the sum of the squares of Deff between the measured values and a "marginalized" effective
  distance this is just the squared Deff integrated over inclination and polarization which
  is proportional to (F+^2 + Fx^2)^(-1)
  """
  latitude,longitude = pt
  gmst = {}
  D_marg_sq = {}
  F_plus = {}
  F_cross = {}
  for ifo in coinc.ifo_list:
    gmst[ifo] = date.XLALGreenwichMeanSiderealTime(coinc.gps[ifo])
    F_plus[ifo], F_cross[ifo] = inject.XLALComputeDetAMResponse(detector_responses[ifo],\
                                longitude,latitude,0,gmst[ifo])
    D_marg_sq[ifo] = 1/(F_plus[ifo]*F_plus[ifo]+F_cross[ifo]*F_cross[ifo])

  delta_D = {}
  effD_diff = 0.0
  effD_sum = 0.0
  Dmarg_diff = 0.0
  Dmarg_sum = 0.0
  delta_D_rms = 0.0
  for ifos in coinc.ifo_coincs:
    effD_diff = coinc.eff_distances[ifos[0]] * coinc.eff_distances[ifos[0]]\
                - coinc.eff_distances[ifos[1]] * coinc.eff_distances[ifos[1]]
    effD_sum = coinc.eff_distances[ifos[0]] * coinc.eff_distances[ifos[0]]\
               + coinc.eff_distances[ifos[1]] * coinc.eff_distances[ifos[1]]
    Dmarg_diff = D_marg_sq[ifos[0]] - D_marg_sq[ifos[1]]
    Dmarg_sum = D_marg_sq[ifos[0]] + D_marg_sq[ifos[1]]
    delta_D[ifos[0]+ifos[1]] = (effD_diff/effD_sum) - (Dmarg_diff/Dmarg_sum)
    delta_D_rms += delta_D[ifos[0]+ifos[1]]*delta_D[ifos[0]+ifos[1]]

  return sqrt(delta_D_rms)

def gridsky(resolution,shifted=False):
  """
  grid the sky up into roughly square regions
  resolution is the length of a side 
  the points get placed at the center of the squares and to 
  first order each square has an area of resolution^2
  if shifted is true, the grids are reported with latitudes
  in (0,pi).  otherwise (default) they lie in (-pi/2,pi/2)
  """
  latitude = 0.0
  longitude = 0.0
  ds = pi*resolution/180.0
  if shifted:
    dlat = 0.0
  else:
    dlat = 0.5*pi
  while latitude <= pi:
    latitude += ds
    longitude = 0.0
    points.append((latitude-dlat, longitude))
    while longitude <= 2.0*pi:
      longitude += ds / abs(sin(latitude))
      points.append((latitude-dlat, longitude))
  #add a point at the south pole
  points.append((0.0-dlat,0.0))
  #there's some slop so get rid of it and only focus on points on the sphere
  sphpts = []
  if shifted:
    latmin = 0.0
    latmax = pi
  else:
    latmin = -pi/2
    latmax = pi/2
  for pt in points:
    if pt[0] > latmax or pt[0] < latmin or pt[1] > 2*pi or pt[1] < 0.0:
      pass
    else:
      sphpts.append(pt)
  return sphpts

def map_grids(coarsegrid,finegrid,coarseres=4.0):
  """
  takes the two grids (lists of lat/lon tuples) and returns a dictionary
  where the points in the coarse grid are the keys and lists of tuples of
  points in the fine grid are the values
  """
  fgtemp = finegrid[:]
  coarsedict = {}
  
  ds = coarseres*pi/180.0
  for cpt in coarsegrid:
    flist = []
    for fpt in fgtemp:
      if (cpt[0]-fpt[0])*(cpt[0]-fpt[0]) <= ds*ds/4.0 and \
         (cpt[1]-fpt[1])*(cpt[1]-fpt[1])*sin(cpt[0])*sin(cpt[0]) <=  ds*ds/4.0:
        flist.append(fpt)
    coarsedict[cpt] = flist
    for rpt in flist:
      fgtemp.remove(rpt)

  return coarsedict, fgtemp

##############################################################################
#
#          class definitions
#
##############################################################################

class SkyPoints(list):
  """
  useful class for doing sky localization.
  assumes that it's a list of (latitude,longitude,L) tuples
  and contains: 
    * a method for sorting those lists
    * a method for writing itself to disk
  """
  def nsort(self,n):
    """
    in place sort of (latitude,longitude,dt,dD...) tuples 
    according to the values in the nth column
    """
    super(SkyPoints,self).sort(key=operator.itemgetter(n))

  def write(self,fname,comment=None,gz=True):
    """
    write the grid to a text file
    """ 
    grid = '#  ra' + '\t' + 'dec' + '\t' + 'dt' + '\t' + 'dD' + '\n'
    self.nsort(3)
    self.nsort(2)
    for pt in self:
      grid += str(pt[1]) + '\t' + str(pt[0]) + '\t' + str(pt[2]) + '\t' + str(pt[3]) + '\n'
    if comment:
      grid += '# ' + comment
    if gz:
      f = gzip.open(fname, 'w')
    else:
      f = open(fname, 'w')
    f.write(grid)
    f.close()  

class CoincData(object):
  """
  simple container for the information needed to run the sky localization code
  """
  def __init__(self):
    """
    here are all the things we need
    """
    #start with data needed for every coinc
    self.ifo_list = []
    self.ifo_coincs = []
    
    self.snr = {}
    self.gps = {}
    self.eff_distances = {}
    self.mass1 = {}
    self.mass2 = {}
    
    self.time = None
    
    #this stuff is only needed for injections
    self.is_injection = False
    self.latitude_inj = None
    self.longitude_inj = None
    self.mass1_inj = None 
    self.mass2_inj = None
    self.eff_distances_inj = {}

  
  def set_ifos(self,ifolist):
    """
    set the ifo_list ans ifo_coincs from the list of ifos involved
    in the coincidence
    """
    self.ifo_list = ifolist
    self.ifo_coincs.extend(list(glue.iterutils.choices(self.ifo_list,2)))

  def set_snr(self,snrdict):
    self.snr = snrdict
 
  def set_gps(self,gpsdict):
    self.gps = gpsdict
    self.time = min(t for t in self.gps.values())
 
  def set_effDs(self,effDdict):
    self.eff_distances = effDdict

  def set_masses(self,m1,m2):
    self.mass1 = m1
    self.mass2 = m2

  def set_inj_params(self,lat,lon,m1,m2,effDs):
    """
    set all of the injection parameters at once
    """
    self.latitude_inj = lat
    self.longitude_inj = lon
    self.mass1_inj = m1
    self.mass2_inj = m2
    self.eff_distances_inj = effDs

class Coincidences(list):
  """
  takes input in either the form of coire files or coinc tables (xml format)
  and produces a list of CoincData objects
  """
  
  def __init__(self,files,filetype='coinctable'):
    """
    files is assumend to be a list of filenames
    """
    if filetype == 'coinctable':
      self.get_coincs_from_coinctable(files)
    elif filetype == 'coire':
      self.get_coincs_from_coire(files)
    else:
      print >>sys.stdout, 'Unknown input file type.'
      sys.exit(0)
   
  def get_coincs_from_coinctable(self,files):
    """
    read data from coinc tables (xml format)
    
    FIXME: currently assumes one coinc per file!!!
    """
    for file in files:
      coinc = CoincData()
      xmldoc = utils.load_filename(file)
      sngltab = tab.get_table(xmldoc,lsctables.SnglInspiralTable.tableName)
      coinc.set_snr(dict((row.ifo, row.snr) for row in sngltab))
      coinc.set_gps(dict((row.ifo, LIGOTimeGPS(row.get_end())) for row in sngltab))
      coinc.set_effDs(dict((row.ifo,row.eff_distance) for row in sngltab))
      coinc.set_masses(dict((row.ifo, row.mass1) for row in sngltab), \
                       dict((row.ifo, row.mass2) for row in sngltab))
      ctab = tab.get_table(xmldoc,lsctables.CoincInspiralTable.tableName)
      coinc.set_ifos(ctab[0].get_ifos())
      
      try:
        simtab = tab.get_table(xmldoc,lsctables.SimInspiralTable.tableName)
        row = siminsptab[0]
        effDs_inj = {}
        for ifo in coinc.ifo_list:
          if ifo == 'H1':
            effDs_inj[ifo] = row.eff_dist_h
          elif ifo == 'L1':
            effDs_inj[ifo] = row.eff_dist_l
          elif ifo == 'V1':
            effDs_inj[ifo] = row.eff_dist_v
        coinc.set_inj_params(row.latitude,row.longitude,row.mass1,row.mass2, \
                             effDs_inj)
        coinc.is_injection = True
      #FIXME: name the exception!
      except:
        pass

      self.append(coinc)
  
  def get_coincs_from_coire(self,files,stat='snr'):
    """
    uses CoincInspiralUtils to get data from old-style (coire'd) coincs
    """
    coincTrigs = CoincInspiralUtils.coincInspiralTable()
    inspTrigs = SnglInspiralUtils.ReadSnglInspiralFromFiles(files, \
                                  mangle_event_id = True,verbose=None)
    statistic = CoincInspiralUtils.coincStatistic(stat,None,None)
    coincTrigs = CoincInspiralUtils.coincInspiralTable(inspTrigs,statistic)

    try:
      inspInj = SimInspiralUtils.ReadSimInspiralFromFiles(files)
      coincTrigs.add_sim_inspirals(inspInj)
    #FIXME: name the exception!
    except:
      pass

    #now extract the relevant information into CoincData objects
    for ctrig in coincTrigs:
      coinc = CoincData()
      coinc.set_ifos(ctrig.get_ifos()[1])
      coinc.set_gps(dict((trig.ifo,LIGOTimeGPS(trig.get_end())) for trig in ctrig))
      coinc.set_snr(dict((trig.ifo,getattr(ctrig,trig.ifo).snr) for trig in ctrig))
      coinc.set_effDs(dict((trig.ifo,getattr(ctrig,trig.ifo).eff_distance) for trig in ctrig))
      coinc.set_masses(dict((trig.ifo,getattr(ctrig,trig.ifo).mass1) for trig in ctrig), \
                       dict((trig.ifo,getattr(ctrig,trig.ifo).mass2) for trig in ctrig))
      
      try:
        effDs_inj = {}
        for ifo in coinc.ifo_list:
          if ifo == 'H1':
            effDs_inj[ifo] = getattr(ctrig,'sim').eff_dist_h
          elif ifo == 'L1':
            effDs_inj[ifo] = getattr(ctrig,'sim').eff_dist_l
          elif ifo == 'V1':
            effDs_inj[ifo] = getattr(ctrig,'sim').eff_dist_v
        coinc.set_inj_params(getattr(ctrig,'sim').latitude,getattr(ctrig,'sim').longitude, \
                             getattr(ctrig,'sim').mass1,getattr(ctrig,'sim').mass2, effDs_inj)
        coinc.is_injection = True
        #FIXME: name the exception!
      except:
        pass
      
      self.append(coinc)

##############################################################################
#
#          table definitions and functions for populating them
#
##############################################################################

class SkyLocTable(tab.Table):
  tableName = "SkyLoc:table"
  validcolumns = {
    "end_time": "int_4s",
    "comb_snr": "real_4",
    "ifos": "lstring",
    "ra": "real_4",
    "dec": "real_4",
    "a60dt": "real_4",
    "a90dt": "real_4",
    "a60dt60dD": "real_4",
    "a90dt90dD": "real_4",
    "min_eff_distance": "real_4",
    "skymap": "lstring",
    "grid": "lstring"
    }
    
class SkyLocRow(object):
  __slots__ = SkyLocTable.validcolumns.keys()
  
  def get_ifos(self):
    """
    Return a set of the instruments for this row.
    """
    return lsctables.instrument_set_from_ifos(self.ifos)

  def set_ifos(self, instruments):
    """
    Serialize a sequence of instruments into the ifos
    attribute.  The instrument names must not contain the ","
    character.
    """
    self.ifos = lsctables.ifos_from_instrument_set(instruments)

SkyLocTable.RowType = SkyLocRow

class SkyLocInjTable(tab.Table):
  tableName = "SkyLocInj:table"
  validcolumns = {
    "end_time": "int_4s",
    "ifos": "lstring",
    "comb_snr": "real_4",
    "h1_snr": "real_4",
    "l1_snr": "real_4",
    "v1_snr": "real_4",
    "ra": "real_4",
    "dec": "real_4",
    "dt_area": "real_4",
    "rank_area": "real_4",
    "delta_t_rss": "real_8",
    "delta_D_rss": "real_8",
    "h1_eff_distance": "real_4",
    "l1_eff_distance": "real_4",
    "v1_eff_distance": "real_4",
    "mass1": "real_4",
    "mass2": "real_4"
    }
    
class SkyLocInjRow(object):
  __slots__ = SkyLocInjTable.validcolumns.keys()

  def get_ifos(self):
    """
    Return a set of the instruments for this row.
    """
    return lsctables.instrument_set_from_ifos(self.ifos)

  def set_ifos(self, instruments):
    """
    Serialize a sequence of instruments into the ifos
    attribute.  The instrument names must not contain the ","
    character.
    """
    self.ifos = lsctables.ifos_from_instrument_set(instruments)


SkyLocInjTable.RowType = SkyLocInjRow

class GalaxyTable(tab.Table):
  tableName = "Galaxy:table"
  validcolumns = {
    "end_time": "int_4s",
    "name": "lstring",
    "ra": "real_8",
    "dec": "real_8",
    "distance_kpc": "real_8",
    "distance_error": "real_8",
    "luminosity_mwe": "real_8",
    "metal_correction": "real_8",
    "magnitude_error": "real_8"
    }

class GalaxyRow(object):
  __slots__ = GalaxyTable.validcolumns.keys()

GalaxyTable.RowType = GalaxyRow

def populate_SkyLocTable(skyloctable,coinc,adt60,adt90,adt60dD60,adt90dD90,\
                         pt,grid_fname,skymap_fname=None):
  """
  populate a row in a skyloctable
  """
  row = skyloctable.RowType()
  
  row.end_time = coinc.time
  row.set_ifos(coinc.ifo_list)
  rhosquared = 0.0
  for ifo in coinc.ifo_list:
    rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
  row.comb_snr = sqrt(rhosquared)
  row.dec,row.ra  = pt[0],pt[1]
  row.a60dt = adt60
  row.a90dt = adt90
  row.a60dt60dD = adt60dD60
  row.a90dt90dD = adt90dD90
  row.min_eff_distance = min(effD for effD in coinc.eff_distances.values())
  if skymap_fname:
    row.skymap = os.path.basename(str(skymap_fname))
  else:
    row.skymap = skymap_fname
  row.grid = os.path.basename(str(grid_fname))

  skyloctable.append(row)
  
def populate_SkyLocInjTable(skylocinjtable,coinc,dt_area,rank_area, \
                            dtrss_inj,dDrss_inj):
  """
  record injection data in a skylocinjtable
  """
  row = skylocinjtable.RowType()

  row.end_time = coinc.time
  row.set_ifos(coinc.ifo_list)
  rhosquared = 0.0
  for ifo in coinc.ifo_list:
    rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
  row.comb_snr = sqrt(rhosquared)
  try:  
    row.h1_snr = coinc.snr['H1']
  except:
    row.h1_snr = None
  try:  
    row.l1_snr = coinc.snr['L1']
  except:
    row.l1_snr = None
  try:  
    row.v1_snr = coinc.snr['V1']
  except:
    row.v1_snr = None
  row.ra = coinc.longitude_inj
  row.dec = coinc.latitude_inj
  row.dt_area = dt_area
  row.rank_area = rank_area
  row.delta_t_rss = dtrss_inj
  row.delta_D_rss = dDrss_inj
  try:
    row.h1_eff_distance = coinc.eff_distances_inj['H1']
  except:
    row.h1_eff_distance = None
  try:
    row.l1_eff_distance = coinc.eff_distances_inj['L1']
  except:
    row.l1_eff_distance = None
  try:
    row.v1_eff_distance = coinc.eff_distances_inj['V1']
  except:
    row.v1_eff_distance = None
  row.mass1 = coinc.mass1_inj
  row.mass2 = coinc.mass2_inj

  skylocinjtable.append(row)

def populate_GalaxyTable(galaxytable,coinc,galaxy):
  """
  record galaxy data in a galaxytable
  """
  row = galaxytable.RowType()

  row.end_time = coinc.time
  row.name = galaxy.name
  row.ra = galaxy.ra
  row.dec = galaxy.dec
  row.distance_kpc = galaxy.distance_kpc
  row.luminosity_mwe = galaxy.luminosity_mwe
  row.magnitude_error = galaxy.magnitude_error
  row.distance_error = galaxy.distance_error
  row.metal_correction = galaxy.metal_correction

  galaxytable.append(row)





