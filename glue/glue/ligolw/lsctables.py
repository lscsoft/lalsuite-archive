"""
LSC Table definitions.  These have been painstakingly copied from
support/include/LIGOLwXMLHeaders.h.  Yes, I'm sure there are typos.
"""

from xml import sax

import metaio
from glue import lal
from glue import segments
import ligolw


#
# =============================================================================
#
#                            Convenience Functions
#
# =============================================================================
#


def New(Type):
	"""
	Convenience function for constructing pre-defined LSC tables.

	Example:
		import lsctables

		table = lsctables.New(lsctables.ProcessTable)
	"""
	table = Type(sax.xmlreader.AttributesImpl({u"Name": Type.tableName}))
	for key, value in table.validcolumns.items():
		table.appendChild(metaio.Column(sax.xmlreader.AttributesImpl({u"Name": ":".join(Type.tableName.split(":")[:-1]) + ":" + key, u"Type": value})))
	table.appendChild(metaio.Stream(sax.xmlreader.AttributesImpl({u"Name": Type.tableName})))
	return table


def Is(Type, elem):
	"""
	Convenience function to check that elem is a Table of type Type.
	"""
	if elem.tagName != ligolw.Table.tagName:
		return False
	return elem.getAttribute("Name") == Type.tableName


#
# =============================================================================
#
#                          processgroup:process:table
#
# =============================================================================
#

class ProcessTable(metaio.Table):
	tableName = "processgroup:process:table"
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

	def _appendRow(self, row):
		metaio.Table.appendRow(self, row)

	def appendRow(self, row):
		if row.process_id in self.keys():
			raise ligolw.ElementError, "duplicate process ID %s" % row.process_id
		self._appendRow(row)

	def __getitem__(self, key):
		"""
		Return the row having process ID equal to key.
		"""
		for row in self:
			if row.process_id == key:
				return row
		raise KeyError, "process ID %s not found" % key

	def __setitem__(self, key, value):
		"""
		If a row has proces ID equal to key, replace it with value,
		otherwise append value as a new row.  Note:
		value.process_id need not equal key.
		"""
		for i in range(len(self)):
			if self.rows[i].process_id == key:
				self.rows[i] = value
				return
		self._appendRow(value)

	def __delitem__(self, key):
		"""
		Delete all rows having process ID key.
		"""
		for i in range(len(self)):
			if self.rows[i].process_id == key:
				del self.rows[i]
				return
		raise KeyError, "process ID %s not found" % key

	def __contains__(self, key):
		"""
		Return True if a row has process ID equal to key, otherwise
		return False.
		"""
		for row in self:
			if row.process_id == key:
				return True
		return False

	def keys(self):
		return [row.process_id for row in self]

class Process(metaio.TableRow):
	__slots__ = ProcessTable.validcolumns.keys()

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


#
# =============================================================================
#
#                   process_paramsgroup:process_params:table
#
# =============================================================================
#

class ProcessParamsTable(metaio.Table):
	tableName = "process_paramsgroup:process_params:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"param": "lstring",
		"type": "lstring",
		"value": "lstring"
	}

	def _appendRow(self, row):
		metaio.Table.appendRow(self, row)

	def appendRow(self, row):
		if (row.process_id, row.param) in [(r.process_id, r.param) for r in self]:
			raise ligolw.ElementError, "duplicate parameter %s for process ID %s" % (row.param, row.process_id)
		if row.type not in metaio.Types:
			raise ligolw.ElementError, "unrecognized Type attribute %s" % row.type
		self._appendRow(row)

	def get_program(self, key):
		"""
		Return the name of the program associated with process ID
		key.
		"""
		for row in self:
			if row.process_id == key:
				return row.program
		raise KeyError, "process ID %s not found" % key

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
			raise KeyError, "process ID %s not found" % key
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
		map(self.appendRow, params)

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

class ProcessParams(metaio.TableRow):
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
#                   search_summarygroup:search_summary:table
#
# =============================================================================
#

class SearchSummaryTable(metaio.Table):
	tableName = "search_summarygroup:search_summary:table"
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
			raise KeyError, "process ID %s not found" % key
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
		map(self.appendRow, params)

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

class SearchSummary(metaio.TableRow):
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
#                  search_summvarsgroup:search_summvars:table
#
# =============================================================================
#

class SearchSummVarsTable(metaio.Table):
	tableName = "search_summvarsgroup:search_summvars:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"name": "lstring",
		"string": "lstring",
		"value": "real_8"
	}

class SearchSummVars(metaio.TableRow):
	__slots__ = SearchSummVarsTable.validcolumns.keys()

SearchSummVarsTable.RowType = SearchSummVars


#
# =============================================================================
#
#                       sngl_burstgroup:sngl_burst:table
#
# =============================================================================
#

class SnglBurstTable(metaio.Table):
	tableName = "sngl_burstgroup:sngl_burst:table"
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

class SnglBurst(metaio.TableRow):
	__slots__ = SnglBurstTable.validcolumns.keys()

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


#
# =============================================================================
#
#                    sngl_inspiralgroup:sngl_inspiral:table
#
# =============================================================================
#

class SnglInspiralTable(metaio.Table):
	tableName = "sngl_inspiralgroup:sngl_inspiral:table"
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

class SnglInspiral(metaio.TableRow):
	__slots__ = SnglInspiralTable.validcolumns.keys()

SnglInspiralTable.RowType = SnglInspiral


#
# =============================================================================
#
#                    sngl_ringdowngroup:sngl_ringdown:table
#
# =============================================================================
#

class SnglRingDownTable(metaio.Table):
	tableName = "sngl_ringdowngroup:sngl_ringdown:table"
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

class SnglRingDown(metaio.TableRow):
	__slots__ = SnglRingDownTable.validcolumns.keys()

SnglRingDownTable.RowType = SnglRingDown


#
# =============================================================================
#
#                   multi_inspiralgroup:multi_inspiral:table
#
# =============================================================================
#

class MultiInspiralTable(metaio.Table):
	tableName = "multi_inspiralgroup:multi_inspiral:table"
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

class MultiInspiral(metaio.TableRow):
	__slots__ = MultiInspiralTable.validcolumns.keys()

MultiInspiralTable.RowType = MultiInspiral


#
# =============================================================================
#
#                     sim_inspiralgroup:sim_inspiral:table
#
# =============================================================================
#

class SimInspiralTable(metaio.Table):
	tableName = "sim_inspiralgroup:sim_inspiral:table"
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
		"simulation_id": "real_4"
	}

class SimInspiral(metaio.TableRow):
	__slots__ = SimInspiralTable.validcolumns.keys()

SimInspiralTable.RowType = SimInspiral


#
# =============================================================================
#
#                        sim_burstgroup:sim_burst:table
#
# =============================================================================
#

class SimBurstTable(metaio.Table):
	tableName = "sim_burstgroup:sim_burst:table"
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

class SimBurst(metaio.TableRow):
	__slots__ = SimBurstTable.validcolumns.keys()

SimBurstTable.RowType = SimBurst


#
# =============================================================================
#
#                     sim_ringdowngroup:sim_ringdown:table
#
# =============================================================================
#

class SimRingDownTable(metaio.Table):
	tableName = "sim_ringdowngroup:sim_ringdown:table"
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

class SimRingDown(metaio.TableRow):
	__slots__ = SimRingDownTable.validcolumns.keys()

SimRingDownTable.RowType = SimRingDown


#
# =============================================================================
#
#                       summ_valuegroup:summ_value:table
#
# =============================================================================
#

class SummValueTable(metaio.Table):
	tableName = "summ_valuegroup:summ_value:table"
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

class SummValue(metaio.TableRow):
	__slots__ = SummValueTable.validcolumns.keys()

SummValueTable.RowType = SummValue


#
# =============================================================================
#
#                  sim_inst_paramsgroup:sim_inst_params:table
#
# =============================================================================
#

class SimInstParamsTable(metaio.Table):
	tableName = "sim_inst_paramsgroup:sim_inst_params:table"
	validcolumns = {
		"simulation_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"value": "real_8"
	}

class SimInstParams(metaio.TableRow):
	__slots__ = SimInstParamsTable.validcolumns.keys()

SimInstParamsTable.RowType = SimInstParams


#
# =============================================================================
#
#                       stochasticgroup:stochastic:table
#
# =============================================================================
#

class StochasticTable(metaio.Table):
	tableName = "stochasticgroup:stochastic:table"
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

class Stochastic(metaio.TableRow):
	__slots__ = StochasticTable.validcolumns.keys()

StochasticTable.RowType = Stochastic


#
# =============================================================================
#
#                        stochsummgroup:stochsumm:table
#
# =============================================================================
#

class StochSummTable(metaio.Table):
	tableName = "stochsummgroup:stochsumm:table"
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

class StochSumm(metaio.TableRow):
	__slots__ = StochSummTable.validcolumns.keys()

StochSummTable.RowType = StochSumm


#
# =============================================================================
#
#                            external_trigger:table
#
# =============================================================================
#

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

class ExtTriggers(metaio.TableRow):
	__slots__ = ExtTriggersTable.validcolumns.keys()

ExtTriggersTable.RowType = ExtTriggers


#
# =============================================================================
#
#                           filtergroup:filter:table
#
# =============================================================================
#

class FilterTable(metaio.Table):
	tableName = "filtergroup:filter:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"program": "lstring",
		"start_time": "int_4s",
		"filter_name": "lstring",
		"comment": "lstring"
	}

class Filter(metaio.TableRow):
	__slots__ = FilterTable.validcolumns.keys()

FilterTable.RowType = Filter


#
# =============================================================================
#
#                          segmentgroup:segment:table
#
# =============================================================================
#

class SegmentTable(metaio.Table):
	tableName = "segmentgroup:segment:table"
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

class Segment(metaio.TableRow):
	__slots__ = SegmentTable.validcolumns.keys()

SegmentTable.RowType = Segment


#
# =============================================================================
#
#                  segment_def_mapgroup:segment_def_map:table
#
# =============================================================================
#

class SegmentDefMapTable(metaio.Table):
	tableName = "segment_def_mapgroup:segment_def_map:table"
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

class SegmentDefMap(metaio.TableRow):
	__slots__ = SegmentDefMapTable.validcolumns.keys()

SegmentDefMapTable.RowType = SegmentDefMap


#
# =============================================================================
#
#                  segment_definergroup:segment_definer:table
#
# =============================================================================
#

class SegmentDefTable(metaio.Table):
	tableName = "segment_definergroup:segment_definer:table"
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

class SegmentDef(metaio.TableRow):
	__slots__ = SegmentDefTable.validcolumns.keys()

SegmentDefTable.RowType = SegmentDef


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
	SegmentDefTable.tableName: SegmentDefTable
}


class LIGOLWContentHandler(metaio.LIGOLWContentHandler):
	"""
	ContentHandler that redirects Table elements with known structure
	to the definitions in this module, using the Table element in
	metaio as a fall-back for unrecognized tables.
	"""
	def startTable(self, attrs):
		try:
			return TableByName[attrs["Name"]](attrs)
		except KeyError:
			return metaio.Table(attrs)
