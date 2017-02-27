#!/usr/bin/env python

"""
This script provides a way to inject very specific sources.
It generates an xml file containing the sources specified in
a textfile that is provided as argument.
"""

__author__ = "Jeroen Meidam"
__credits__ = ["Jeroen Meidam"]
__maintainer__ = "Jeroen Meidam"
__email__ = "jeroen.meidam@ligo.org"
__status__ = ""

usage="""  tom_xml_injections.py [options]
  Will generate an xml file containing sources with
  parameters as given in sourcefile.
  Example sourcefile (all fields assumed to be present in this file, except nonGRparams):
    event_n name mass1 mass2 phi1 phi2 theta1 theta2 a1 a2 distance inclination polarization ra dec injtime
    0 gw150914 40.835 33.263 0.0 0.0 0.0 0.0 0.0 0.0 555.597 2.837 1.429 -1.262 1.950 966391866
    1 gw151226 17.285 7.477 0.0 0.0 0.0 0.0 0.0 0.0 513.201 2.654 2.897 0.710 0.348 966391866
    2 lvt151012 22.467 19.101 0.0 0.0 0.0 0.0 0.0 0.0 598.574 2.161 2.009 0.667 0.131 966391866
"""

###############################################################################
#
# LOAD LIBRARIES
#
###############################################################################

from optparse import OptionParser
from os import path,system


# import modules
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

import lal
from numpy import cos,sin,sqrt,genfromtxt

LLO_index = lal.LALDetectorIndexLLODIFF
LHO_index = lal.LALDetectorIndexLHODIFF
Virgo_index = lal.LALDetectorIndexVIRGODIFF

# define a content handler
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)


###############################################################################
#
# DEFINITIONS
#
###############################################################################

nonGRparamnames = ["dchi0","dchi1","dchi2","dchi3","dchi4","dchi5","dchi5l","dchi6","dchi6l","dchi7"]
nonGRparamnames+= ["dbeta1","dbeta2","dbeta3"]
nonGRparamnames+= ["dalpha1","dalpha2","dalpha3","dalpha4","dalpha5"]
nonGRparamnames+= ["dsigma1","dsigma2","dsigma3","dsigma4"]
nonGRparamnames+= ["dxi1","dxi2","dxi3","dxi4","dxi5","dxi6"]

class source():
  def __init__(self,label,mass1=1.4,mass2=1.4,spin1=[0.,0.,0.],spin2=[0.,0.,0.],
               dist=1.0,incl=0.,pol=0.,ra=0.,dec=0.,epoch=0,slide_h=0,slide_l=0,slide_v=0,
               nonGRparams={}):
    """
    Distance in Mpc
    Mass in Msol
    """
    self.label = label
    self.mass1 = mass1
    self.mass2 = mass2
    self.spin1x = spin1[0]
    self.spin1y = spin1[1]
    self.spin1z = spin1[2]
    self.spin2x = spin2[0]
    self.spin2y = spin2[1]
    self.spin2z = spin2[2]
    self.distance = dist
    self.inclination = incl
    self.polarization = pol
    self.ra = ra
    self.dec = dec
    self.geocent_end_time = epoch
    self.nonGRparams=nonGRparams


def get_distances(distance,inclination,ra,dec,psi,gmst):
  """
  compute the effective distance of the signal
  from SimInspiralUtils.c (e.g. LALInspiralSiteTimeAndDist)
  """

  #initialize distances with real distance and compute splus and scross
  eff_dist_h = eff_dist_l = eff_dist_v = 2.0*distance
  cosiota = cos(inclination)
  splus = -(1.0 + cosiota*cosiota)
  scross = -2.0 * cosiota

  #compute the response of the detectors
  rp_lho, rc_lho = lal.ComputeDetAMResponse(lal.CachedDetectors[LHO_index].response,
                                        ra, dec, psi, gmst)
  rp_llo, rc_llo = lal.ComputeDetAMResponse(lal.CachedDetectors[LLO_index].response,
                                        ra, dec, psi, gmst)
  rp_vir, rc_vir = lal.ComputeDetAMResponse(lal.CachedDetectors[Virgo_index].response,
                                        ra, dec, psi, gmst)

  #compute the effective distance for LHO
  eff_dist_h = eff_dist_h/sqrt(splus*splus*rp_lho*rp_lho + scross*scross*rc_lho*rc_lho )

  #compute the effective distance for LLO
  eff_dist_l = eff_dist_l/sqrt(splus*splus*rp_llo*rp_llo + scross*scross*rc_llo*rc_llo )

  #compute the effective distance for Virgo
  eff_dist_v = eff_dist_v/sqrt(splus*splus*rp_vir*rp_vir + scross*scross*rc_vir*rc_vir )

  return eff_dist_h,eff_dist_l,eff_dist_v


def get_times(geocent_end_time,ra,dec,gmst):
  """
  compute gps end times for each detector
  from SimInspiralUtils.c (e.g. LALInspiralSiteTimeAndDist)
  """
  #initialize end times with geocentric value
  h_end_time = l_end_time = v_end_time = geocent_end_time;

  time_diff_lho = lal.TimeDelayFromEarthCenter(lal.CachedDetectors[LHO_index].location,
                                               ra,dec,gmst)
  time_diff_llo = lal.TimeDelayFromEarthCenter(lal.CachedDetectors[LLO_index].location,
                                               ra,dec,gmst)
  time_diff_vir = lal.TimeDelayFromEarthCenter(lal.CachedDetectors[Virgo_index].location,
                                               ra,dec,gmst)

  h_end_time += time_diff_lho
  l_end_time += time_diff_llo
  v_end_time += time_diff_vir

  return h_end_time,l_end_time,v_end_time

def get_component_spins(phi,theta,a):
  spinz = a*cos(theta)
  spinx = a*sin(theta)*cos(phi)
  spiny = a*sin(theta)*sin(phi)
  return [spinx,spiny,spinz]


def run_inspinj(approximant,gps_start_time,distance,flow,longitude,latitude,polarization,inclination,output,inspinj="lalapps_inspinj"):
  """Create an xml document containing only one row with specified parameters
  """

  xmltmpname = path.join(path.dirname(output),"tmp.xml")
  command = "%s "%inspinj
  command += "--output %s "%xmltmpname
  command += "--seed 1000 "
  command += "--waveform %s "%approximant
  command += "--gps-start-time %d "%gps_start_time
  command += "--gps-end-time %d "%(gps_start_time+1000)
  command += "--time-step %d "%1000
  command += "--min-distance %f "%(distance*1000.0)
  command += "--max-distance %f "%(distance*1000.0)
  command += "--d-distr volume "
  command += "--m-distr componentMass "
  command += "--min-mass1 5.0 "
  command += "--max-mass1 40.0 "
  command += "--min-mass2 5.0 "
  command += "--max-mass2 40.0 "
  command += "--f-lower %f "%flow
  command += "--l-distr fixed "
  command += "--longitude %f "%longitude
  command += "--latitude %f "%latitude
  command += "--i-distr fixed "
  command += "--polarization %f "%polarization
  command += "--fixed-inc %f "%inclination
  command += "--t-distr fixed "
  command += "--enable-spin "
  command += "--min-spin1 0.01 "
  command += "--max-spin1 0.9 "
  command += "--min-spin2 0.01 "
  command += "--max-spin2 0.9"

  system(command)
  xmldoc = ligolw_utils.load_filename(xmltmpname, contenthandler = LIGOLWContentHandler, verbose = False)
  system("rm %s"%xmltmpname)

  return xmldoc




def main(testarguments=""):

  #############################################################################
  #
  # ARGUMENT PARSING
  #
  #############################################################################

  """
  For testing, main can be called from within e.g. ipython with commandline provided as:
  ['--sourcefile', 'testfile', '--flow', '10.0']
  """

  parser=OptionParser(usage)
  parser.add_option("-f","--flow",default=20.0,help="Lower frequency cutoff used for all the waveforms")#,dest="flow"
  parser.add_option("-a","--approximant",default="IMRPhenomPv2threePointFivePN",help="Approximant name including PN order")
  parser.add_option("-s","--sourcefile",help="File containing sources with parameters")
  parser.add_option("-o","--output",help="output xml file")
  parser.add_option("-i","--inspinj",default='lalapps_inspinj',help="lalapps_inspinj executable location")
  if testarguments=="":
    (opts,args) = parser.parse_args()
  else:
    (opts,args) = parser.parse_args(testarguments)

  if opts.sourcefile:
    sourcefile = opts.sourcefile
  else:
    print "must provide sourcefile with --sourcefile"
    return 1

  if opts.output:
    output = opts.output
  else:
    output = "./custom.xml"

  flow = float(opts.flow)
  inspinj = opts.inspinj

  approximant=opts.approximant


  data = genfromtxt(sourcefile,converters={1: lambda s: str(s)})
  data = data[1:] #skip header
  N_sources = len(data)

  with open(sourcefile,"r") as fp:
    header = fp.readline().strip().split()

  index = {}
  for p,i in zip(header,range(len(header))):
    index[p] = i

  #Create a temporary xml file that contains enough sources (rows) to function as template
  xmltmpname = path.join(path.dirname(output),"tmp_template.xml")
  xmlout = output
  command = "%s "%inspinj
  command += "--output %s "%xmltmpname
  command += "--seed 1000 "
  command += "--waveform %s "%approximant
  command += "--gps-start-time %d "%966384015
  command += "--gps-end-time %d "%(966384015+N_sources*100)
  command += "--time-step 100 "
  command += "--min-distance 100000 "
  command += "--max-distance 150000 "
  command += "--d-distr volume "
  command += "--m-distr componentMass "
  command += "--min-mass1 5.0 "
  command += "--max-mass1 40.0 "
  command += "--min-mass2 5.0 "
  command += "--max-mass2 40.0 "
  command += "--f-lower %f "%flow
  command += "--l-distr random "
  command += "--i-distr uniform "
  command += "--t-distr fixed "
  command += "--enable-spin "
  command += "--min-spin1 0.01 "
  command += "--max-spin1 0.9 "
  command += "--min-spin2 0.01 "
  command += "--max-spin2 0.9"
  #p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
  system(command)
  xmldoc_template = ligolw_utils.load_filename(xmltmpname, contenthandler = LIGOLWContentHandler, verbose = False)

  # retrieve the sim_inspiral table.  these are list-like
  # objects of rows.  the row objects' attributes are the column names
  sim_inspiral_table_template = lsctables.SimInspiralTable.get_table(xmldoc_template)

  sources = []
  for s in data:
    print "adding",s[index["name"]]
    spin1 = get_component_spins(s[index['phi1']],s[index['theta1']],s[index['a1']])
    spin2 = get_component_spins(s[index['phi2']],s[index['theta2']],s[index['a2']])
    nonGRparams = {}
    for n in nonGRparamnames:
      if index.has_key(n):
        nonGRparams[n]=s[index[n]]
    sources.append(source(s[index["name"]],mass1=s[index["mass1"]],mass2=s[index["mass2"]],
                  spin1=spin1,spin2=spin2,
                  dist=s[index["distance"]],incl=s[index["inclination"]],pol=s[index["polarization"]],
                  ra=s[index["ra"]],dec=s[index["dec"]],epoch=s[index["injtime"]],
                  nonGRparams=nonGRparams))

  #update the rows in the xml with the source parameters
  print "using lalapps_inspinj to calculate gps times from detector orientations/locations and ra,dec..."
  for row,s in zip(sim_inspiral_table_template,sources):

    #run inspinj for each row to properly calculate the gps times
    xmldoc = run_inspinj(approximant,s.geocent_end_time,s.distance,flow,s.ra,s.dec,s.polarization,s.inclination,path.join(path.dirname(output),"tmp.xml"),inspinj=inspinj)
    sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)
    params = sim_inspiral_table[0]

    row.geocent_end_time    = params.geocent_end_time
    row.geocent_end_time_ns = params.geocent_end_time_ns
    row.h_end_time    = params.h_end_time
    row.h_end_time_ns = params.h_end_time_ns
    row.l_end_time    = params.l_end_time
    row.l_end_time_ns = params.l_end_time_ns
    row.v_end_time    = params.v_end_time
    row.v_end_time_ns = params.v_end_time_ns
    row.end_time_gmst = params.end_time_gmst

    row.mass1 = s.mass1
    row.mass2 = s.mass2
    M = s.mass1+s.mass2
    eta = (s.mass1*s.mass2)/(M*M)
    mchirp = M*pow(eta,3./5.)
    row.mchirp = mchirp
    row.eta = eta

    row.spin1x = s.spin1x
    row.spin1y = s.spin1y
    row.spin1z = s.spin1z
    row.spin2x = s.spin2x
    row.spin2y = s.spin2y
    row.spin2z = s.spin2z

    row.longitude = params.longitude
    row.latitude = params.latitude
    row.inclination = params.inclination
    row.coa_phase = 0.0
    row.polarization = params.polarization

    row.distance = params.distance
    row.eff_dist_h = params.eff_dist_h
    row.eff_dist_l = params.eff_dist_l
    row.eff_dist_v = params.eff_dist_v

    k = s.nonGRparams.keys()
    if len(k) > 0:
      if "dchi0" in k: row.dchi0 = s.nonGRparams["dchi0"]
      if "dchi1" in k: row.dchi1 = s.nonGRparams["dchi1"]
      if "dchi2" in k: row.dchi2 = s.nonGRparams["dchi2"]
      if "dchi3" in k: row.dchi3 = s.nonGRparams["dchi3"]
      if "dchi4" in k: row.dchi4 = s.nonGRparams["dchi4"]
      if "dchi5" in k: row.dchi5 = s.nonGRparams["dchi5"]
      if "dchi5l" in k: row.dchi5l = s.nonGRparams["dchi5l"]
      if "dchi6" in k: row.dchi6 = s.nonGRparams["dchi6"]
      if "dchi6l" in k: row.dchi6l = s.nonGRparams["dchi6l"]
      if "dchi7" in k: row.dchi7 = s.nonGRparams["dchi7"]
      if "dalpha1" in k: row.dalpha1 = s.nonGRparams["dalpha1"]
      if "dalpha2" in k: row.dalpha2 = s.nonGRparams["dalpha2"]
      if "dalpha3" in k: row.dalpha3 = s.nonGRparams["dalpha3"]
      if "dalpha4" in k: row.dalpha4 = s.nonGRparams["dalpha4"]
      if "dalpha5" in k: row.dalpha5 = s.nonGRparams["dalpha5"]
      if "dbeta1" in k: row.dbeta1 = s.nonGRparams["dbeta1"]
      if "dbeta2" in k: row.dbeta2 = s.nonGRparams["dbeta2"]
      if "dbeta3" in k: row.dbeta3 = s.nonGRparams["dbeta3"]
      if "dsigma1" in k: row.dsigma1 = s.nonGRparams["dsigma1"]
      if "dsigma2" in k: row.dsigma2 = s.nonGRparams["dsigma2"]
      if "dsigma3" in k: row.dsigma3 = s.nonGRparams["dsigma3"]
      if "dsigma4" in k: row.dsigma4 = s.nonGRparams["dsigma4"]
      if "dxi1" in k: row.dxi1 = s.nonGRparams["dxi1"]
      if "dxi2" in k: row.dxi2 = s.nonGRparams["dxi2"]
      if "dxi3" in k: row.dxi3 = s.nonGRparams["dxi3"]
      if "dxi4" in k: row.dxi4 = s.nonGRparams["dxi4"]
      if "dxi5" in k: row.dxi5 = s.nonGRparams["dxi5"]
      if "dxi6" in k: row.dxi6 = s.nonGRparams["dxi6"]


  # Save the new xml and remove the temp file
  ligolw_utils.write_filename(xmldoc_template,xmlout)
  print "saved xml as %s"%xmlout
  system("rm %s"%xmltmpname)




###############################################################################
#
# START THE MAIN FUNCTION
#
###############################################################################

if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS
	exit(main())




