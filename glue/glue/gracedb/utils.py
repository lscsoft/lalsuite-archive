#!/usr/bin/python

from pylal import Fr
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables

##############################################################################
#
#          useful variables
#
##############################################################################

#Need these for each search: inspiral, burst, etc.
InspiralCoincDef = lsctables.CoincDef(search = u"MBTA_inspiral", \
                                      search_coinc_type = 0, \
                                      description = \
                                      u"sngl_inspiral<-->sngl_inspiral \
                                      coincidences")
InspiralCoincId = 'coinc_inspiral:coinc_event_id:0'

#list of detectors participating in the coinc
#MBTA only sends triples to gracedb at the time being so this list is
#simply for convenience.  burst and future inspiral searches should
#construct this list on the fly
all_detectors = ['H1', 'L1', 'V1']

#this is the subset of SnglInspiralTable.validcolumn.keys() that
#are assigned from MBTA coinc triggers
MBTA_set_keys = ['ifo', 'search', 'end_time', 'end_time_ns', 'mass1', 'mass2',\
               'mchirp', 'mtotal', 'eta', 'snr', 'eff_distance', 'event_id',\
               'process_id', 'channel']

#this dictionary is the simplest way to assign event_id's
#collisions are are taken care of in the process of conversion to sqlite
insp_event_id_dict = {'H1': 'sngl_inspiral:event_id:0',\
                 'L1': 'sngl_inspiral:event_id:1',\
                 'V1': 'sngl_inspiral:event_id:2'}

##############################################################################
#
#          convenience functions
#
##############################################################################

def compute_mchirp_eta(m1,m2):
  """
  compute and return mchirp and eta for a given pair of masses 
  """
  
  mtot = m1 + m2
  mu = m1*m2/mtot
  eta = mu/mtot
  mchirp = pow(mu,3/5)*mtot
  
  return float(mchirp), float(eta)

def write_output_files(root_dir, xmldoc, log_content, \
                       xml_fname = 'coinc.xml', log_fname = 'event.log'):
  """
  write the xml-format coinc tables and log file
  """

  f = open(root_dir+'/'+xml_fname,'w')
  xmldoc.write(f)
  f.close()

  f = open(root_dir+'/'+log_fname,'w')
  f.write(log_content)
  f.close()

##############################################################################
#
#          table populators
#
##############################################################################

def populate_inspiral_tables(MBTA_frame, set_keys = MBTA_set_keys, \
                             process_id = 'process:process_id:0', \
                             event_id_dict = insp_event_id_dict, \
                             coinc_event_id='coinc_inspiral:coinc_event_id:0', \
                             detectors=all_detectors):
  """
  create xml file and populate the SnglInspiral and CoincInspiral tables from a
  coinc .gwf file from MBTA
  xmldoc: xml file to append the tables to
  MBTA_frame: frame file to get info about triggers from
  set_keys: columns in the SnglInspiral Table to set
  process_id: process_id
  event_id_dict: {ifo:event_id} dictionary to assign event_id's
  coinc_event_id: coinc_event_id
  detectors: detectors participating in the coinc

  returns xmldoc and contents of the comment field
  """
  #initialize xml document
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  detectors.sort()
  #dictionaries to store about individual triggers
  end_time_s = {}
  end_time_ns = {}
  snr = {}
  mass1 = {}
  mass2 = {}
  Deff = {}
  mchirp = {}
  eta = {}

  #extract the information from the frame file
  events = Fr.frgetevent(MBTA_frame)
  for event in events:
    if 'MbtaHLV' in event['name']:
      log_data = event['comment']
      far = [line.split(':')[1].split()[0] for line in log_data.splitlines() if \
             'False Alarm Rate' in line][0]
      for ifo in detectors:
        end_time_s[ifo], end_time_ns[ifo] = str(event[ifo+':end_time']).split('.')
        snr[ifo] = float(event[ifo+':SNR'])
        mass1[ifo] = float(event[ifo+':mass1'])
        mass2[ifo] = float(event[ifo+':mass2'])
        mchirp[ifo], eta[ifo] = compute_mchirp_eta(mass1[ifo],mass2[ifo])
        Deff[ifo] = float(event[ifo+':eff_distance'])

  #fill the SnglInspiralTable
  sin_table = lsctables.New(lsctables.SnglInspiralTable)
  xmldoc.childNodes[0].appendChild(sin_table)
  for ifo in detectors:
    row = sin_table.RowType()
    row.ifo = ifo
    row.search = 'MBTA'
    row.end_time = int(end_time_s[ifo])
    row.end_time_ns = int(end_time_ns[ifo])
    row.mass1 = mass1[ifo]
    row.mass2 = mass2[ifo]
    row.mchirp = mchirp[ifo]
    row.mtotal = mass1[ifo] + mass2[ifo]
    row.eta = eta[ifo]
    row.snr = snr[ifo]
    row.eff_distance = Deff[ifo]
    row.event_id = event_id_dict[ifo]
    row.process_id = process_id
    row.channel = ''
    #zero out the rest of the columns
    #should work in chi2 and chi2cut 
    for key in sin_table.validcolumns.keys():
      if key not in set_keys:
        setattr(row,key,None)
    sin_table.append(row)

  #CoincInspiralTable
  #using the conventions found in:
  #https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/S6Plan/ 
  #090505160219S6PlanningNotebookCoinc_and_Experiment_Tables_ihope_implementation?
  #highlight=%28coinc%29|%28table%29

  cin_table = lsctables.New(lsctables.CoincInspiralTable)
  xmldoc.childNodes[0].appendChild(cin_table)
  row = cin_table.RowType()
  row.set_ifos(detectors)
  row.coinc_event_id = coinc_event_id
  row.end_time = int(end_time_s['H1'])
  row.end_time_ns = int(end_time_ns['H1'])
  row.mass = (sum(mass1.values()) + sum(mass2.values()))/3
  row.mchirp = sum(mchirp.values())/3
  #the snr here is really the snr NOT effective snr
  row.snr = pow(sum([x*x for x in snr.values()]),0.5)
  #far is triggers/day
  row.false_alarm_rate = float(far)
  row.combined_far = 0
  cin_table.append(row)

  return xmldoc, log_data
    
def populate_coinc_tables(xmldoc, coinc_event_id, event_id_dict, CoincDef, \
                          process_id = 'process:process_id:0', \
                          coinc_def_id ='coinc_definer:coinc_def_id:0', \
                          time_slide_id = None, likelihood = None, nevents = 3, \
                          detectors = all_detectors):
  """
  populate a set of coinc tables
  xmldoc:  xml file to append the tables to
  CoincDef: pre-initialized CoincDef table row
  detectors: detectors participating in the coinc
  """
  #CoincTable
  coinc_table = lsctables.New(lsctables.CoincTable)
  xmldoc.childNodes[0].appendChild(coinc_table)
  row = coinc_table.RowType()
  row.process_id = process_id
  row.coinc_event_id =  coinc_event_id #'coinc_inspiral:coinc_event_id:0'
  row.coinc_def_id = coinc_def_id
  row.time_slide_id = time_slide_id
  row.set_instruments(detectors)
  row.nevents = nevents
  row.likelihood = likelihood
  coinc_table.append(row)

  #CoincMapTable
  coinc_map_table = lsctables.New(lsctables.CoincMapTable)
  xmldoc.childNodes[0].appendChild(coinc_map_table)
  for ifo in detectors:
    row = coinc_map_table.RowType()
    row.coinc_event_id = coinc_event_id
    row.table_name = lsctables.SnglInspiralTable.tableName.split(':')[0]
    row.event_id = event_id_dict[ifo]
    coinc_map_table.append(row)

  #CoincDefTable
  coinc_def_table = lsctables.New(lsctables.CoincDefTable)
  xmldoc.childNodes[0].appendChild(coinc_def_table)
  row = coinc_def_table.RowType()
  row.coinc_def_id = coinc_def_id
  row.search = CoincDef.search
  row.search_coinc_type = CoincDef.search_coinc_type
  row.description = InspiralCoincDef.description
  coinc_def_table.append(row)
  

##############################################################################
#
#          usage example
#
##############################################################################

#populate the tables
#xmldoc, log_data = populate_inspiral_tables("MbtaFake-930909680-16.gwf")
#populate_coinc_tables(xmldoc,InspiralCoincId,insp_event_id_dict,InspiralCoincDef)
#write the output
#write_output_files('.', xmldoc, log_data)

