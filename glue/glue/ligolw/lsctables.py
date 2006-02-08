"""
LSC Table definitions.  These have been painstakingly copied from
support/include/LIGOLwXMLHeaders.h.  Yes, I'm sure there are typos.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import re
from xml import sax

from glue import lal
from glue import segments
import ligolw
import metaio


#
# =============================================================================
#
#                            Convenience Functions
#
# =============================================================================
#

def New(Type, columns = None):
	"""
	Convenience function for constructing pre-defined LSC tables.  The
	optional columns argument is a list of the names of the columns the
	table should be constructed with.  If columns = None, then the
	table is constructed with all valid columns included (pass columns
	= [] to create a table with no columns).

	Example:
		import lsctables

		table = lsctables.New(lsctables.ProcessTable)
	"""
	table = Type(sax.xmlreader.AttributesImpl({u"Name": Type.tableName}))
	for key, value in table.validcolumns.items():
		if (columns == None) or (key in columns):
			table.appendChild(metaio.Column(sax.xmlreader.AttributesImpl({u"Name": ":".join(Type.tableName.split(":")[:-1]) + ":" + key, u"Type": value})))
	table.appendChild(metaio.TableStream(sax.xmlreader.AttributesImpl({u"Name": Type.tableName})))
	return table


def IsTableElement(Type, elem):
	"""
	Convenience function to check that an element is a Table of type
	Type.
	"""
	if elem.tagName != ligolw.Table.tagName:
		return False
	return metaio.CompareTableNames(elem.getAttribute("Name"), Type.tableName) == 0


def IsTableProperties(Type, tagname, attrs):
	"""
	Convenience function to check that the given tag name and
	attributes match those of a Table of type Type.
	"""
	if tagname != ligolw.Table.tagName:
		return False
	return metaio.CompareTableNames(attrs["Name"], Type.tableName) == 0


#
# =============================================================================
#
#                 Table and Row Class With a Mapping Protocol
#
# =============================================================================
#

class LSCTableDict(object):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows.
	"""
	def __init__(self, rows):
		"""
		Initialize the mapping for the list of rows rows.
		"""
		self.rows = rows

	def __getitem__(self, key):
		"""
		Return the row matching key.
		"""
		for row in self.rows:
			if row._has_key(key):
				return row
		raise KeyError, repr(key)

	def __setitem__(self, key, value):
		"""
		If a row has key equal to key, replace it with value,
		otherwise append value as a new row.  Note: the row key
		carried by value need not equal key, but this behaviour may
		change in the future.
		"""
		for i in range(len(self.rows)):
			if self.rows[i]._has_key(key):
				self.rows[i] = value
				return
		# FIXME: should we call _set_key() on value to force it to have
		# the key that was searched for?
		self.append(value)

	def __delitem__(self, key):
		"""
		Delete all rows having the given key.
		"""
		for i in range(len(self.rows)):
			if self.rows[i]._has_get(key):
				del self.rows[i]
				return
		raise KeyError, repr(key)

	def __contains__(self, key):
		"""
		Return True if a row has key equal to key, otherwise return
		False.
		"""
		for row in self.rows:
			if row._has_key(key):
				return True
		return False

	def keys(self):
		return [row._get_key() for row in self.rows]

	def get_idmap(self):
		"""
		Return the key --> row object mapping for this table.
		"""
		map = {}
		for row in self.rows:
			key = row._get_key()
			if key in map:
				raise ligolw.ElementError, "duplicate key %s" % repr(key)
			map[key] = row
		return map



class LSCTable(metaio.Table):
	def __init__(self, attrs):
		metaio.Table.__init__(self, attrs)
		self.dict = LSCTableDict(self.rows)


# We don't subclass metaio.TableRow because that defeats the __slots__
# feature.
class LSCTableRow(object):
	__slots__ = []

	# Prefix with underscores to avoid collision with column names
	def _get_key(self):
		"""
		Get the unique ID for this row.
		"""
		raise KeyError, "row object does not define a key column"

	def _set_key(self, key):
		"""
		Set the unique ID for this row.
		"""
		raise KeyError, "row object does not define a key column"

	def _has_key(self, key):
		"""
		Check if this row's unique ID is equal to key.
		"""
		raise KeyError, "row object does not define a key column"


#
# =============================================================================
#
#                                process:table
#
# =============================================================================
#

class ProcessTable(LSCTable):
	tableName = "process:table"
	validcolumns = {
		"program": "lstring",
		"version": "lstring",
		"cvs_repository": "lstring",
		"cvs_entry_time": "int_4s",
		"comment": "lstring",
		"is_online": "int_4s",
		"node": "lstring",
		"username": "lstring",
		"unix_procid": "int_4s",
		"start_time": "int_4s",
		"end_time": "int_4s",
		"jobid": "int_4s",
		"domain": "lstring",
		"ifos": "lstring",
		"process_id": "ilwd:char"
	}

class Process(LSCTableRow):
	__slots__ = ProcessTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return key == self.process_id

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in ProcessTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

ProcessTable.RowType = Process

class ProcessIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "process:process_id", n)

#
# =============================================================================
#
#                                lfn:table
#
# =============================================================================
#

class LfnTable(LSCTable):
	tableName = "lfn:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"lfn_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"start_time": "int_4s",
		"end_time": "int_4s",
	}

class Lfn(LSCTableRow):
	__slots__ = LfnTable.validcolumns.keys()

	def _get_key(self):
		return self.lfn_id

	def _set_key(self, key):
		self.lfn_id = key

	def _has_key(self, key):
		return key == self.lfn_id

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in LfnTable.validcolumns.keys():
			if key == "lfn_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

LfnTable.RowType = Lfn

class LfnIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "lfn:lfn_id", n)


#
# =============================================================================
#
#                             process_params:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class ProcessParamsTable(metaio.Table):
	tableName = "process_params:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"param": "lstring",
		"type": "lstring",
		"value": "lstring"
	}

	def append(self, row):
		if row.type not in metaio.Types:
			raise ligolw.ElementError, "ProcessParamsTable.append():  unrecognized type \"%s\"" % row.type
		metaio.Table.append(self, row)

	def get_program(self, key):
		"""
		Return the name of the program associated with process ID
		key.
		"""
		for row in self:
			if row.process_id == key:
				return row.program
		raise KeyError, repr(key)

	def set_program(self, key, value):
		"""
		Set the program for all entries with process ID key to
		value.
		"""
		for row in self:
			if row.process_id == key:
				row.program = value

	def __getitem__(self, key):
		"""
		Return a list of rows matching the process ID key.
		"""
		params = []
		for row in self:
			if row.process_id == key:
				params.append(row)
		if not len(params):
			raise KeyError, repr(key)
		# sort by process ID, then parameter name (all rows should
		# be unique by this measure).
		params.sort(lambda a, b: cmp((a.process_id, a.param), (b.process_id, b.param)))
		return params

	def __setitem__(self, key, params):
		"""
		Replace the rows having process ID key with the
		ProcessParams in the list params, appending the list to the
		table if there are no rows with that ID.
		"""
		del self[key]
		map(self.append, params)

	def __delitem__(self, key):
		"""
		Delete all rows having process ID key.
		"""
		self.filterRows(lambda row: row.process_id != key)

	def __contains__(self, key):
		"""
		Return True if a row has process ID equal to key, otherwise
		return False.
		"""
		for row in self:
			if row.process_id == key:
				return True
		return False

class ProcessParams(LSCTableRow):
	__slots__ = ProcessParamsTable.validcolumns.keys()

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in ProcessParamsTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

ProcessParamsTable.RowType = ProcessParams


#
# =============================================================================
#
#                             search_summary:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class SearchSummaryTable(metaio.Table):
	tableName = "search_summary:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"shared_object": "lstring",
		"lalwrapper_cvs_tag": "lstring",
		"lal_cvs_tag": "lstring",
		"comment": "lstring",
		"ifos": "lstring",
		"in_start_time": "int_4s",
		"in_start_time_ns": "int_4s",
		"in_end_time": "int_4s",
		"in_end_time_ns": "int_4s",
		"out_start_time": "int_4s",
		"out_start_time_ns": "int_4s",
		"out_end_time": "int_4s",
		"out_end_time_ns": "int_4s",
		"nevents": "int_4s",
		"nnodes": "int_4s"
	}

	def __getitem__(self, key):
		"""
		Return a list of rows matching the process ID key.
		"""
		summaries = []
		for row in self:
			if row.process_id == key:
				summaries.append(row)
		if not len(summaries):
			raise KeyError, repr(key)
		# sort by process ID, then output segment (all rows should
		# be unique by this measure).
		summaries.sort(lambda a, b: cmp((a.process_id, a.get_out()), (b.process_id, b.get_out())))
		return summaries

	def __setitem__(self, key, params):
		"""
		Replace the rows having process ID key with the
		ProcessParams in the list params, appending the list to the
		table if there are no rows with that ID.
		"""
		del self[key]
		map(self.append, params)

	def __delitem__(self, key):
		"""
		Delete all rows having process ID key.
		"""
		self.filterRows(lambda row: row.process_id != key)

	def __contains__(self, key):
		"""
		Return True if a row has process ID equal to key, otherwise
		return False.
		"""
		for row in self:
			if row.process_id == key:
				return True
		return False

	def get_inlist(self):
		return segments.segmentlist([row.get_in() for row in self])

	def get_outlist(self):
		return segments.segmentlist([row.get_out() for row in self])

class SearchSummary(LSCTableRow):
	__slots__ = SearchSummaryTable.validcolumns.keys()

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in SearchSummaryTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

	def get_in(self):
		"""
		Return the input segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.in_start_time, self.in_start_time_ns), lal.LIGOTimeGPS(self.in_end_time, self.in_end_time_ns))

	def set_in(self, seg):
		"""
		Set the input segment.
		"""
		self.in_start_time, self.in_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.in_end_time, self.in_end_time_ns = seg[1].seconds, seg[1].nanoseconds

	def get_out(self):
		"""
		Get the output segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.out_start_time, self.out_start_time_ns), lal.LIGOTimeGPS(self.out_end_time, self.out_end_time_ns))

	def set_out(self, seg):
		"""
		Set the output segment.
		"""
		self.out_start_time, self.out_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.out_end_time, self.out_end_time_ns = seg[1].seconds, seg[1].nanoseconds

SearchSummaryTable.RowType = SearchSummary


#
# =============================================================================
#
#                            search_summvars:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class SearchSummVarsTable(metaio.Table):
	tableName = "search_summvars:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"name": "lstring",
		"string": "lstring",
		"value": "real_8"
	}

class SearchSummVars(LSCTableRow):
	__slots__ = SearchSummVarsTable.validcolumns.keys()

SearchSummVarsTable.RowType = SearchSummVars


#
# =============================================================================
#
#                               sngl_burst:table
#
# =============================================================================
#

class SnglBurstTable(LSCTable):
	tableName = "sngl_burst:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"peak_time": "int_4s",
		"peak_time_ns": "int_4s",
		"duration": "real_4",
		"central_freq": "real_4",
		"bandwidth": "real_4",
		"amplitude": "real_4",
		"snr": "real_4",
		"confidence": "real_4",
		"clusterT": "real_4",
		"peak_dof": "real_4",
		"event_id": "int_8s"
	}

class SnglBurst(LSCTableRow):
	__slots__ = SnglBurstTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

	def set_peak(self, gps):
		self.peak_time, self.peak_time_ns = gps.seconds, gps.nanoseconds

	def get_period(self):
		start = lal.LIGOTimeGPS(self.start_time, self.start_time_ns)
		return segments.segment(start, start + self.duration)

	def set_period(self, period):
		# FIXME: should duration be forced to type float?
		self.start_time, self.start_time_ns = period[0].seconds, period[0].nanoseconds
		self.duration = float(period.duration())

	def get_band(self):
		return segments.segment(self.central_freq - self.bandwidth/2.0, self.central_freq + self.bandwidth/2.0)

	def set_band(self, band):
		self.central_freq = (band[0] + band[1])/2.0
		self.bandwidth = band.duration()

SnglBurstTable.RowType = SnglBurst

class SnglBurstIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sngl_burst:event_id", n)


#
# =============================================================================
#
#                             sngl_inspiral:table
#
# =============================================================================
#

class SnglInspiralTable(LSCTable):
	tableName = "sngl_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"template_duration": "real_8",
		"event_duration": "real_8",
		"amplitude": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"mtotal": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"f_final": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "int_4s",
		"sigmasq": "real_8",
		"rsqveto_duration": "real_4",
		"event_id": "int_8s"
	}

class SnglInspiral(LSCTableRow):
	__slots__ = SnglInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

SnglInspiralTable.RowType = SnglInspiral

class SnglInspiralIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sngl_inspiral:event_id", n)


#
# =============================================================================
#
#                             sngl_ringdown:table
#
# =============================================================================
#

class SnglRingDownTable(LSCTable):
	tableName = "sngl_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"frequency": "real_4",
		"quality": "real_4",
		"mass": "real_4",
		"spin": "real_4",
		"snr": "real_4",
		"eff_distance": "real_4",
		"sigma_sq": "real_8",
		"event_id": "int_8s"
	}

class SnglRingDown(LSCTableRow):
	__slots__ = SnglRingDownTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

SnglRingDownTable.RowType = SnglRingDown

class SnglRingDownIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sngl_ringdown:event_id", n)


#
# =============================================================================
#
#                             multi_inspiral:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class MultiInspiralTable(metaio.Table):
	tableName = "multi_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifos": "lstring",
		"search": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"amplitude": "real_4",
		"ifo1_eff_distance": "real_4",
		"ifo2_eff_distance": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"ifo1_snr": "real_4",
		"ifo2_snr": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "real_4",
		"sigmasq": "real_4",
		"ligo_axis_ra": "real_4",
		"ligo_axis_dec": "real_4",
		"ligo_angle": "real_4",
		"ligo_angle_sig": "real_4",
		"inclination": "real_4",
		"polarization": "real_4"
	}

class MultiInspiral(LSCTableRow):
	__slots__ = MultiInspiralTable.validcolumns.keys()

MultiInspiralTable.RowType = MultiInspiral


#
# =============================================================================
#
#                              sim_inspiral:table
#
# =============================================================================
#

class SimInspiralTable(LSCTable):
	tableName = "sim_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_end_time": "int_4s",
		"geocent_end_time_ns": "int_4s",
		"h_end_time": "int_4s",
		"h_end_time_ns": "int_4s",
		"l_end_time": "int_4s",
		"l_end_time_ns": "int_4s",
		"g_end_time": "int_4s",
		"g_end_time_ns": "int_4s",
		"t_end_time": "int_4s",
		"t_end_time_ns": "int_4s",
		"v_end_time": "int_4s",
		"v_end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"source": "lstring",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"distance": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"inclination": "real_4",
		"coa_phase": "real_4",
		"polarization": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"spin1x": "real_4",
		"spin1y": "real_4",
		"spin1z": "real_4",
		"spin2x": "real_4",
		"spin2y": "real_4",
		"spin2z": "real_4",
		"theta0": "real_4",
		"phi0": "real_4",
		"f_lower": "real_4",
		"f_final": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"eff_dist_g": "real_4",
		"eff_dist_t": "real_4",
		"eff_dist_v": "real_4",
		"simulation_id": "ilwd:char"
	}

class SimInspiral(LSCTableRow):
	__slots__ = SimInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

SimInspiralTable.RowType = SimInspiral

class SimInspiralIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sim_inspiral:simulation_id", n)


#
# =============================================================================
#
#                               sim_burst:table
#
# =============================================================================
#

class SimBurstTable(LSCTable):
	tableName = "sim_burst:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_peak_time": "int_4s",
		"geocent_peak_time_ns": "int_4s",
		"h_peak_time": "int_4s",
		"h_peak_time_ns": "int_4s",
		"l_peak_time": "int_4s",
		"l_peak_time_ns": "int_4s",
		"peak_time_gmst": "real_8",
		"dtminus": "real_4",
		"dtplus": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"coordinates": "lstring",
		"polarization": "real_4",
		"hrss": "real_4",
		"hpeak": "real_4",
		"distance": "real_4",
		"freq": "real_4",
		"tau": "real_4",
		"zm_number": "int_4s",
		"simulation_id": "ilwd:char"
	}

class SimBurst(LSCTableRow):
	__slots__ = SimBurstTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

SimBurstTable.RowType = SimBurst

class SimBurstIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sim_burst:simulation_id", n)


#
# =============================================================================
#
#                              sim_ringdown:table
#
# =============================================================================
#

class SimRingDownTable(LSCTable):
	tableName = "sim_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"coordinates": "lstring",
		"geocent_start_time": "int_4s",
		"geocent_start_time_ns": "int_4s",
		"h_start_time": "int_4s",
		"h_start_time_ns": "int_4s",
		"l_start_time": "int_4s",
		"l_start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"mass": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"distance": "real_4",
		"inclination": "real_4",
		"polarization": "real_4",
		"epsilon": "real_4",
		"spin": "real_4",
		"frequency": "real_4",
		"quality": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"h0": "real_4",
		"hrss": "real_4",
		"hrss_h": "real_4",
		"hrss_l": "real_4",
		"simulation_id": "ilwd:char"
	}

class SimRingDown(LSCTableRow):
	__slots__ = SimRingDownTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

SimRingDownTable.RowType = SimRingDown

class SimRingDownIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sim_ringdown:simulation_id", n)


#
# =============================================================================
#
#                               summ_value:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class SummValueTable(metaio.Table):
	tableName = "summ_value:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"ifo": "lstring",
		"name": "lstring",
		"value": "real_4",
		"comment": "lstring"
	}

class SummValue(LSCTableRow):
	__slots__ = SummValueTable.validcolumns.keys()

SummValueTable.RowType = SummValue


#
# =============================================================================
#
#                            sim_inst_params:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class SimInstParamsTable(metaio.Table):
	tableName = "sim_inst_params:table"
	validcolumns = {
		"simulation_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"value": "real_8"
	}

class SimInstParams(LSCTableRow):
	__slots__ = SimInstParamsTable.validcolumns.keys()

SimInstParamsTable.RowType = SimInstParams

class SimInstParamsIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "sim_inst_params:simulation_id", n)


#
# =============================================================================
#
#                               stochastic:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class StochasticTable(metaio.Table):
	tableName = "stochastic:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"duration": "int_4s",
		"duration_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"cc_stat": "real_8",
		"cc_sigma": "real_8"
	}

class Stochastic(LSCTableRow):
	__slots__ = StochasticTable.validcolumns.keys()

StochasticTable.RowType = Stochastic


#
# =============================================================================
#
#                               stochsumm:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class StochSummTable(metaio.Table):
	tableName = "stochsumm:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"y_opt": "real_8",
		"error": "real_8"
	}

class StochSumm(LSCTableRow):
	__slots__ = StochSummTable.validcolumns.keys()

StochSummTable.RowType = StochSumm


#
# =============================================================================
#
#                            external_trigger:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class ExtTriggersTable(metaio.Table):
	tableName = "external_trigger:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"det_alts": "lstring",
		"det_band": "lstring",
		"det_fluence": "lstring",
		"det_fluence_int": "lstring",
		"det_name": "lstring",
		"det_peak": "lstring",
		"det_peak_int": "lstring",
		"det_snr": "lstring",
		"email_time": "int_4s",
		"event_dec": "real_4",
		"event_dec_err": "real_4",
		"event_epoch": "lstring",
		"event_err_type": "lstring",
		"event_ra": "real_4",
		"event_ra_err": "real_4",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"event_type": "lstring",
		"event_z": "real_4",
		"event_z_err": "real_4",
		"notice_comments": "lstring",
		"notice_id": "lstring",
		"notice_sequence": "lstring",
		"notice_time": "int_4s",
		"notice_type": "lstring",
		"notice_url": "lstring",
		"obs_fov_dec": "real_4",
		"obs_fov_dec_width": "real_4",
		"obs_fov_ra": "real_4",
		"obs_fov_ra_width": "real_4",
		"obs_loc_ele": "real_4",
		"obs_loc_lat": "real_4",
		"obs_loc_long": "real_4",
		"ligo_fave_lho": "real_4",
		"ligo_fave_llo": "real_4",
		"ligo_delay": "real_4",
		"event_number_gcn": "int_4s",
		"event_number_grb": "lstring",
		"event_status": "int_4s"
	}

class ExtTriggers(LSCTableRow):
	__slots__ = ExtTriggersTable.validcolumns.keys()

ExtTriggersTable.RowType = ExtTriggers


#
# =============================================================================
#
#                                 filter:table
#
# =============================================================================
#

# FIXME: how to subclass LSCTable?
class FilterTable(metaio.Table):
	tableName = "filter:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"program": "lstring",
		"start_time": "int_4s",
		"filter_name": "lstring",
		"comment": "lstring"
	}

class Filter(LSCTableRow):
	__slots__ = FilterTable.validcolumns.keys()

FilterTable.RowType = Filter


#
# =============================================================================
#
#                                segment:table
#
# =============================================================================
#

class SegmentTable(LSCTable):
	tableName = "segment:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"active": "int_4s",
		"segnum": "int_4s",
		"insertion_time": "int_4s"
	}

	def segmentlist(self, key, active = 1):
		"""
		Return the segment list having process_id equal to key, and
		with the activity flag having the same sign as active.  The
		list represents exactly the rows in the table, in order;
		in other words it has not been coalesced.
		"""
		if not active:
			raise ValueError, "SegmentTable.segmentlist(): activity flag must != 0."
		return segments.segmentlist([row.get_segment() for row in self if (row.process_id == key) and (row.active * active > 0)])

class Segment(LSCTableRow):
	__slots__ = SegmentTable.validcolumns.keys()

	def _get_key(self):
		return self.segment_id

	def _set_key(self, key):
		self.segment_id = key

	def _has_key(self, key):
		return self.segment_id == key

	def get_segment(self):
		"""
		Return the segment described by this row.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.start_time, self.start_time_ns), lal.LIGOTimeGPS(self.end_time, self.end_time_ns))

	def set_segment(self, segment):
		"""
		Set the segment described by this row.
		"""
		self.start_time, self.start_time_ns = segment[0].seconds, segment[0].nanoseconds
		self.end_time, self.end_time_ns = segment[1].seconds, segment[1].nanoseconds

SegmentTable.RowType = Segment

class SegmentIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "segment:segment_id", n)


#
# =============================================================================
#
#                            segment_def_map:table
#
# =============================================================================
#

class SegmentDefMapTable(LSCTable):
	tableName = "segment_def_map:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"seg_def_map_id": "ilwd:char",
		"segment_cdb": "int_4s",
		"segment_id": "ilwd:char",
		"segment_def_cdb": "int_4s",
		"segment_def_id": "ilwd:char",
		"state_vec_map": "int_4s",
		"insertion_time": "int_4s"
	}

class SegmentDefMap(LSCTableRow):
	__slots__ = SegmentDefMapTable.validcolumns.keys()

	def _get_key(self):
		return self.seg_def_map_id

	def _set_key(self, key):
		self.seg_def_map_id = key

	def _has_key(self, key):
		return self.seg_def_map_id == key

SegmentDefMapTable.RowType = SegmentDefMap

class SegmentDefMapIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "segment_def_map:segment_def_id", n)


#
# =============================================================================
#
#                            segment_definer:table
#
# =============================================================================
#

class SegmentDefTable(LSCTable):
	tableName = "segment_definer:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_def_id": "ilwd:char",
		"run": "lstring",
		"ifos": "lstring",
		"name": "lstring",
		"version": "int_4s",
		"comment": "lstring",
		"state_vec_major": "int_4s",
		"state_vec_minor": "int_4s",
		"insertion_time": "int_4s"
	}

class SegmentDef(LSCTableRow):
	__slots__ = SegmentDefTable.validcolumns.keys()

	def _get_key(self):
		return self.segment_def_id

	def _set_key(self, key):
		self.segment_def_id = key

	def _has_key(self, key):
		return self.segment_def_id == key

SegmentDefTable.RowType = SegmentDef


#
# =============================================================================
#
#                                 slide:table
#
# =============================================================================
#

class TimeSlideTable(metaio.Table):
	tableName = "time_slide:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"ifo": "lstring",
		"offset": "int_4s",
		"offset_ns": "int_4s"
	}

class TimeSlide(LSCTableRow):
	__slots__ = TimeSlideTable.validcolumns.keys()

	def _get_key(self):
		return self.slide_id

	def _set_key(self, key):
		self.slide_id = key

	def _has_key(self, key):
		return self.slide_id == key

TimeSlideTable.RowType = TimeSlide

class TimeSlideIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "time_slide:time_slide_id", n)


#
# =============================================================================
#
#                                 coinc:table
#
# =============================================================================
#

class CoincTable(LSCTable):
	tableName = "coinc:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"coinc_id": "ilwd:char",
		"time_slide_id": "ilwd:char"
	}

class Coinc(LSCTableRow):
	__slots__ = CoincTable.validcolumns.keys()

	def _get_key(self):
		return self.coinc_id

	def _set_key(self, key):
		self.coinc_id = key

	def _has_key(self, key):
		return self.coinc_id == key

CoincTable.RowType = Coinc

class CoincIDs(metaio.ILWD):
	def __init__(self, n = 0):
		metaio.ILWD.__init__(self, "coinc:coinc_id", n)


#
# =============================================================================
#
#                            coinc_event_map:table
#
# =============================================================================
#

class CoincMapTable(LSCTable):
	tableName = "coinc_event_map:table"
	validcolumns = {
		"coinc_id": "ilwd:char",
		"event_id": "ilwd:char"
	}

class CoincMap(LSCTableRow):
	__slots__ = CoincMapTable.validcolumns.keys()

CoincMapTable.RowType = CoincMap


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#

# Table name ---> table type mapping.

TableByName = {
	ProcessTable.tableName: ProcessTable,
	LfnTable.tableName: LfnTable,
	ProcessParamsTable.tableName: ProcessParamsTable,
	SearchSummaryTable.tableName: SearchSummaryTable,
	SearchSummVarsTable.tableName: SearchSummVarsTable,
	SnglBurstTable.tableName: SnglBurstTable,
	SnglInspiralTable.tableName: SnglInspiralTable,
	SnglRingDownTable.tableName: SnglRingDownTable,
	MultiInspiralTable.tableName: MultiInspiralTable,
	SimInspiralTable.tableName: SimInspiralTable,
	SimBurstTable.tableName: SimBurstTable,
	SimRingDownTable.tableName: SimRingDownTable,
	SummValueTable.tableName: SummValueTable,
	SimInstParamsTable.tableName: SimInstParamsTable,
	StochasticTable.tableName: StochasticTable,
	StochSummTable.tableName: StochSummTable,
	ExtTriggersTable.tableName: ExtTriggersTable,
	FilterTable.tableName: FilterTable,
	SegmentTable.tableName: SegmentTable,
	SegmentDefMapTable.tableName: SegmentDefMapTable,
	SegmentDefTable.tableName: SegmentDefTable,
	TimeSlideTable.tableName: TimeSlideTable,
	CoincTable.tableName: CoincTable,
	CoincMapTable.tableName: CoincMapTable
}


class LIGOLWContentHandler(metaio.LIGOLWContentHandler):
	"""
	ContentHandler that redirects Table elements with known structure
	to the definitions in this module, using the Table element in
	metaio as a fall-back for unrecognized tables.
	"""
	def startTable(self, attrs):
		try:
			return TableByName[metaio.StripTableName(attrs["Name"])](attrs)
		except KeyError:
			return metaio.Table(attrs)
