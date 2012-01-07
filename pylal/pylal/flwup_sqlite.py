# Copyright (C) 2012 Cristina Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
######################################################################

"""
This is the followup pipeline (flwup) library for manipulating sqlite
files and other related database/xml data manipulation.  Any followup
method that needs to access external metadata should rely on this
methods to do so.
"""

__author__ = "Cristina Valeria Torres <cristina.torres@ligo.org>,\
Sarah Caudill <sarah.caudill@ligo.org>"


def create_is_playground_func(connection, playground_segs):
    """
    Create a temporary sqlite function which can be used to evaluate
    gps timestamps as either in or outside of the specified playground
    segments.  Invoke this method to create a function which can be
    used by other sqlite queries during the same CONNECTION
    """
    connection.create_function("is_playground",
                               2,
                               lambda seconds, nanoseconds:
    LIGOTimeGPS(seconds, nanoseconds) in playground_segs) 


class DB_summary(object):	
	def __init__(self, connection, live_time_program, file_name, veto_segments_name = None, verbose = False, base=None):
		"""
		Compute and record some summary information about the
		database.
		"""

		self.base = base
		self.connection = connection
		xmldoc = dbtables.get_xml(connection)
		self.file_name = filename
		self.sim_file_name = None
		self.verbose = verbose

		cursor = connection.cursor()

		# find the tables
		try:
			self.sngl_inspiral_table = table.get_table(xmldoc, dbtables.lsctables.SnglInspiralTable.tableName)
		except ValueError:
			self.sngl_inspiral_table = None
		try:
			self.sim_inspiral_table = table.get_table(xmldoc, dbtables.lsctables.SimInspiralTable.tableName)
			# write out the injection file to use in later inspiral jobs
			newxmldoc = ligolw.Document()
			newxmldoc.appendChild(ligolw.LIGO_LW())
			newxmldoc.childNodes[-1].appendChild(self.sim_inspiral_table)
			self.sim_file_name = "sim_"+os.path.split(filename)[1].replace('sqlite','xml')
			utils.write_filename(newxmldoc, self.sim_file_name, gz=False, verbose=verbose)

		except ValueError:
			self.sim_inspiral_table = None
		try:
			self.coinc_def_table = table.get_table(xmldoc, dbtables.lsctables.CoincDefTable.tableName)
			self.coinc_table = table.get_table(xmldoc, dbtables.lsctables.CoincTable.tableName)
			self.time_slide_table = table.get_table(xmldoc, dbtables.lsctables.TimeSlideTable.tableName)
		except ValueError:
			self.coinc_def_table = None
			self.coinc_table = None
			self.time_slide_table = None
		try:
			self.coinc_inspiral_table = table.get_table(xmldoc, dbtables.lsctables.CoincInspiralTable.tableName)
		except ValueError:
			self.coinc_inspiral_table = None
		try:
			self.experiment_summary_table = table.get_table(xmldoc, dbtables.lsctables.ExperimentSummaryTable.tableName)
		except ValueError:
			self.experiment_summary_table = None

		# determine a few coinc_definer IDs
		# FIXME:  don't hard-code the numbers
		if self.coinc_def_table is not None:
			try:
				self.ii_definer_id = self.coinc_def_table.get_coinc_def_id("inspiral", 0, create_new = False)
			except KeyError:
				self.ii_definer_id = None
			try:
				self.si_definer_id = self.coinc_def_table.get_coinc_def_id("inspiral", 1, create_new = False)
			except KeyError:
				self.si_definer_id = None
			try:
				self.sc_definer_id = self.coinc_def_table.get_coinc_def_id("inspiral", 2, create_new = False)
			except KeyError:
				self.sc_definer_id = None
		else:
			self.ii_definer_id = None
			self.si_definer_id = None
			self.sc_definer_id = None

		# retrieve the distinct on and participating instruments
		self.on_instruments_combos = [frozenset(dbtables.lsctables.instrument_set_from_ifos(x)) for x, in cursor.execute("SELECT DISTINCT(instruments) FROM coinc_event WHERE coinc_def_id == ?", (self.ii_definer_id,))]

		# get the segment lists
		self.seglists = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program)
		self.playground_segs = segmentsUtils.S2playground(self.seglists.extent_all())
		self.instruments = set(self.seglists)
		if veto_segments_name is not None:
			self.veto_segments = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
		else:
			self.veto_segments = segments.segmentlistdict()
		self.seglists -= self.veto_segments

class Sngl(object):
	def __init__(self, row, contents,ignore_proc_params=False):
		self.data_tuple = row
		self.contents = contents
        	self.ifo = row[0]
		self.snr = float(row[1])
		self.chisq = float(row[2])
		self.mass1 = float(row[3])
		self.mass2 = float(row[4])
		self.mchirp = float(row[5])
		self.time = float(row[6])
		self.process_id = row[7]
		self.row = contents.sngl_inspiral_table.row_from_cols(row[8:])
		self.inj_file_name = contents.sim_file_name

		self.process_params = []

		if not ignore_proc_params:
			for val in self.contents.connection.cursor().execute("""
SELECT
                program, process_id, param, type, value
        FROM
                process_params
        WHERE
                process_id == ?
			""", (self.process_id,) ):
				self.process_params.append(ProcParam(val[0], val[1], val[2], val[3], val[4]))	


	def get_gps_start_time(self):
		#FIXME probably really slow
		for row in self.process_params:
			if row.param == '--gps-start-time': 
				return int(row.value)

	def get_gps_end_time(self):
		#FIXME probably really slow
		for row in self.process_params:
			if row.param == '--gps-end-time': 
				return int(row.value)
	
	def get_sample_rate(self):
		#FIXME probably really slow
		for row in self.process_params:
			if row.param == '--sample-rate': 
				return int(row.value)

	def get_proc_param(self, param):
		for row in self.process_params:
			if row.param.strip('-') == param.strip('-'):
				return row.value
		
	def switch_ifo(self,ifo):
		self.ifo = ifo
		self.process_params = []

		for val in self.contents.connection.cursor().execute("""
SELECT program, process_id, param, type, value FROM process_params AS proc_params WHERE proc_params.process_id = (SELECT p1.process_id FROM process_params AS p1 JOIN search_summary ON search_summary.process_id == p1.process_id WHERE search_summary.out_start_time <= ? AND search_summary.out_end_time > ? and search_summary.ifos = ? LIMIT 1)
		""", (self.time, self.time, self.ifo) ):
			self.process_params.append(ProcParam(val[0], val[1], val[2], val[3], val[4]))
		#change necessary values in process params table
		channel = None
		for row in self.process_params:
			if row.param == '--channel-name':  
				type, channel = stfu_pipe.figure_out_type(self.time,ifo)
				row.value = channel
				break
		if not channel: 
			print sys.stderr, "Coudn't find process params for trigger aborting"
			sys.exit(1)
		#change necessary values in sngl_inspiral table
		self.row.ifo = ifo	
		self.row.channel = channel.split(':')[1]
		
class flwupRingCoinc(dict):
    """
    Base class responsible for creating an appropriately populated
    object for following up a ringdown event with the followup
    pipeline.
    """
    def __init__(self,connection=None,inputPickle=None):
        """
        This method creates a flwupRingCoinc object, either by
        using a connected sqlite database or importing a previously
        created pickle object of this type.
        """
        if (type(connection) == type(None) \
                and
            type(inputPickle) == type(None)):
            raise Exception("invalid NoneType arguments encountered!")
        if (type(connection) != type(None) \
                and
            type(inputPickle) != type(None)):
            raise Exception("both input arguments not NoneType!")
        self.sngl=list()
        dict.__init__(self)
        ### Need to add method to parse pickle file infor
        ### Also add method that uses sqlite connection seeking info

#class flwupInspiralCoinc(object):

class Coinc(object):
	def __init__(self, row, contents):
		self.contents = contents
		self.sngl_inspiral = {}
		self.sngl_inspiral_coh = {}
		self.coinc_event_id = ilwd.get_ilwdchar(row[0])
		self.combined_far = float(row[1])
		self.snr = float(row[2])
		self.time = float(row[3])
		self.mass = float(row[4])
		self.mchirp = float(row[5])
		self.ifos = row[6]
		self.ifos_set = frozenset(dbtables.lsctables.instrument_set_from_ifos(self.ifos))
		self.instruments = row[7]
		self.instruments_set = frozenset(dbtables.lsctables.instrument_set_from_ifos(self.instruments))
		for val in self.contents.connection.cursor().execute("""
SELECT
                sngl_inspiral.ifo, sngl_inspiral.snr, sngl_inspiral.chisq, sngl_inspiral.mass1, sngl_inspiral.mass2, sngl_inspiral.mchirp, sngl_inspiral.end_time + sngl_inspiral.end_time_ns * 1.0e-9, sngl_inspiral.process_id, sngl_inspiral.*
        FROM
                sngl_inspiral
                JOIN coinc_event_map ON (
                        coinc_event_map.coinc_event_id == ?
                )
        WHERE
                sngl_inspiral.event_id == coinc_event_map.event_id
		""", (str(self.coinc_event_id),) ): 
			self.sngl_inspiral.setdefault(val[0],None)
			self.sngl_inspiral[val[0]] = Sngl(val, contents)

		# add instruments that were on but didn't participate in coinc
		for ifo in self.instruments_set:
			if ifo not in self.ifos_set:
				self.sngl_inspiral[ifo] = self.make_sngl_from_max_ifo(ifo)

		# make a list of sngls appropriate for the coherent code (use max template)
		for ifo in self.instruments_set:
			self.sngl_inspiral_coh[ifo] = self.make_sngl_from_max_ifo(ifo)
		# FIX ME: CURRENTLY THIS WOULD RUN FOR HOURS
		#self.extract_sim()
		self.sim = None

	def extract_sim(self):
		#FIXME FINISH!
		if self.contents.sim_inspiral_table:
			self.sim = self.contents.connection.cursor().execute("""
SELECT
		sim_inspiral.* 
	FROM 
		sim_inspiral
	JOIN 
		coinc_event_map AS mapA ON (mapA.event_id == sim_inspiral.simulation_id)
	JOIN 
		coinc_event_map AS mapB ON (mapA.coinc_event_id == mapB.coinc_event_id)
	JOIN 
		coinc_event ON (mapB.event_id == coinc_event.coinc_event_id)
	WHERE 
		coinc_event.coinc_event_id == ?
			""", (str(self.coinc_event_id),) ).fetchone()

		else: self.sim = None
		
	def get_sample_rate(self):
		#FIXME NEEDS TO CHECK ALL SAMPLE RATES ARE THE SAME!
		for sngl in self.sngl_inspiral.values():
			return sngl.get_sample_rate()

	def max_trigger_ifo(self):
		snr_tuple = [(sngl.snr, sngl.ifo) for sngl in self.sngl_inspiral.values()]
		return max(snr_tuple)[1]

	def make_sngl_from_max_ifo(self, ifo):
		max_sngl = self.sngl_inspiral[self.max_trigger_ifo()]
		sngl = Sngl(max_sngl.data_tuple, max_sngl.contents, ignore_proc_params=True)
		sngl.switch_ifo(ifo)
		return sngl
		
class FUTriggers(object):
	def __init__(self, num=10):
		self.num = num
		self.playground_candidates = []
		self.candidates = []
		self.background_candidates = []

	def add_grb_contents(self, contents, stat='combined_far'):
		for (l, v, n) in [(self.candidates, (self.num,), 'non playground')]:
			if contents.verbose: print >>sys.stderr, "querying %s" % (n,)
			for i,row in enumerate(contents.connection.cursor().execute("""
SELECT
	DISTINCT(coinc_inspiral.coinc_event_id),
	coinc_inspiral.combined_far,
	coinc_inspiral.snr,
	coinc_inspiral.end_time + coinc_inspiral.end_time_ns * 1.0e-9,
	coinc_inspiral.mass,
	coinc_inspiral.mchirp,
	coinc_inspiral.ifos,
	coinc_event.instruments
FROM
	coinc_inspiral
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
	)
ORDER BY
	combined_far
LIMIT ?
		""",v)):
				l.append(Coinc(row, contents))

	def add_contents(self, contents, stat='combined_far'):
		for (l, v, n) in [(self.candidates, (0,0,0,self.num), 'non playground'),(self.playground_candidates, (1,0,0,self.num), 'playground'),(self.background_candidates, (0,1,0,self.num), 'time slides'),(self.candidates, (1,0,0,self.num), 'HACK incorrect playground: NOT playground candidates listed as INSIDE PLAYGROUND SEGMENTS(Check SegDB Versions)'),(self.playground_candidates, (0,0,0,self.num), 'HACK incorrect playground: ARE playground candidates listed OUT OF PLAYGROUND SEGMENTS(Check SegDB Versions)')]:
			if contents.verbose: print >>sys.stderr, "querying %s" % (n,)
			for i,row in enumerate(contents.connection.cursor().execute("""
SELECT
	DISTINCT(coinc_inspiral.coinc_event_id),
	coinc_inspiral.combined_far,
	coinc_inspiral.snr,
	coinc_inspiral.end_time + coinc_inspiral.end_time_ns * 1.0e-9,
	coinc_inspiral.mass,
	coinc_inspiral.mchirp,
	coinc_inspiral.ifos,
	coinc_event.instruments,
	is_playground(coinc_inspiral.end_time, coinc_inspiral.end_time_ns) AS is_play,
	EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) AS is_slide,
	experiment_summary.datatype = "simulation" AS is_sim
FROM
	coinc_inspiral
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
	)
JOIN 
	experiment_map ON (
		experiment_map.coinc_event_id == coinc_event.coinc_event_id 
	)
JOIN experiment_summary ON (
		experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id 
	)
WHERE is_play = ? AND is_slide = ? AND is_sim = ?
ORDER BY
	combined_far
LIMIT ?
		""",v)):
				l.append(Coinc(row, contents))

	def topN(self):
		if len(self.candidates) < self.num: self.num = len(self.candidates)
		trigs = [(t.combined_far, t) for t in self.candidates]
		trigs.sort()
		self.candidates = [t[1] for t in trigs[0:self.num] ]

                if len(self.playground_candidates) < self.num: self.num = len(self.playground_candidates)
		trigs = [(t.combined_far, t) for t in self.playground_candidates]
		trigs.sort()
		self.playground_candidates = [t[1] for t in trigs[0:self.num] ]

                if len(self.background_candidates) < self.num: self.num = len(self.background_candidates)
		trigs = [(t.combined_far, t) for t in self.background_candidates]
		trigs.sort()
		self.background_candidates = [t[1] for t in trigs[0:self.num] ]


def exportCoincToDisk(coinc, dir, cnt, dag, filename=None):
	"""
	This method takes an input from the Coinc Class and writes
	that object to disk as ascii.  The output file created by this
	method is intended to serve as input to new makechecklistwiki
	script. This method returns the base filename in order to know
	which file was just created with the data.
	"""
	#If directy coincHeadings not there make the directory if
	#filename is None
	headingDir= dir + "/coinc_info"
	instruments = "".join(coinc.instruments.split(','))
	ifos = "".join(coinc.ifos.split(','))

	idKeyStr="%s_%s" % (str(coinc.time), instruments)
	if filename==None:
		filename="coincEvent.info"
		stfu_pipe.mkdir(headingDir)
		filename=os.path.normpath(headingDir+'/'+idKeyStr+'_'+ filename)
	fp=open(filename,'w')
	fp.write("#DIR\t\t\tRANK\tFAR\t\tSNR\tIFOS\tINSTRUMENTS\tTIME\t\tMASS\tMCHIRP\t\n")
	fp.write("%-16s\t%d\t%0.2e\t%.2f\t%s\t%s\t\t%.3f\t%.2f\t%.2f\n" % (dir, cnt, float(coinc.combined_far), float(coinc.snr), ifos, instruments, float(coinc.time), float(coinc.mass), float(coinc.mchirp)) )
	fp.write("#DIR\t\t\tIFO\tTIME\t\tSNR\tCHISQ\tMASS1\tMASS2\tMCHIRP\n")
	rowString="%-16s\t%s\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
	content=list()
	for ifo,snglEvent in coinc.sngl_inspiral.items():
		content.append(rowString%(dir, 
					  snglEvent.ifo,
					  float(snglEvent.time),  
					  float(snglEvent.snr),
					  float(snglEvent.chisq),
					  float(snglEvent.mass1),
					  float(snglEvent.mass2),
					  float(snglEvent.mchirp)))
	cache = lal.CacheEntry(instruments, "COINC_INFO_"+dir.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+os.path.abspath(filename))
	dag.output_cache.append(cache)
	fp.writelines(content)
	fp.close()
	return os.path.split(filename)[1]

def exportGPSEventToDisk(tevent, dir, cnt, dag, filename=None):
	"""
	"""
	#If directy coincHeadings not there make the directory if
	#filename is None
	headingDir= dir + "/coinc_info"
	ifos = tevent.instruments
	instruments =tevent.instruments
	time = tevent.time

	idKeyStr="%s_%s" % (str(time), instruments)
	if filename==None:
		filename="coincEvent.info"
		stfu_pipe.mkdir(headingDir)
		filename=os.path.normpath(headingDir+'/'+idKeyStr+'_'+ filename)
	fp=open(filename,'w')
	fp.write("#DIR\t\t\tRANK\tFAR\t\tSNR\tIFOS\tINSTRUMENTS\tTIME\t\tMASS\tMCHIRP\t\n")
	fp.write("%-16s\t%d\t%0.2e\t%.2f\t%s\t%s\t\t%.3f\t%.2f\t%.2f\n" % (dir, cnt, 0, 0, ifos, instruments, float(time), 0, 0) )
	fp.write("#DIR\t\t\tIFO\tTIME\t\tSNR\tCHISQ\tMASS1\tMASS2\tMCHIRP\n")
	rowString="%-16s\t%s\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
	content=list()
	for ifo in tevent.ifos_list:
		content.append(rowString%(dir, 
					  ifo,
					  float(time),  
					  float(0),
					  float(0),
					  float(0),
					  float(0),
					  float(0)))
	cache = lal.CacheEntry(instruments, "COINC_INFO_"+dir.upper(), segments.segment(float(time), float(time)), "file://localhost/"+os.path.abspath(filename))
	dag.output_cache.append(cache)
	fp.writelines(content)
	fp.close()
	return os.path.split(filename)[1]

