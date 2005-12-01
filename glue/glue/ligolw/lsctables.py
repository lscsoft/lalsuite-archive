"""
Table definitions in use by LAL.
"""

#
# These have been painstakingly copied from
# support/include/LIGOLwXMLHeaders.h.  Yes, I'm sure there are typos.
#

from xml import sax

import ligolw
import metaio
from glue import lal
from glue import segments


def New(Type):
	"""
	Convenience functions for constructing pre-defined LAL tables.

	Example:
		import laltables

		table = laltables.New(laltables.ProcessTable)
	"""
	attrs = sax.xmlreader.AttributesImpl({u"Name": Type.tableName})
	table = Type(attrs)
	for name, type in table.validcolumns.items():
		table.appendChild(ligolw.Column(sax.xmlreader.AttributesImpl({u"Name": ":".join(Type.tableName.split(":")[:-1]) + ":" + name, u"Type": type})))
	table.appendChild(metaio.Stream(attrs))
	return table


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

class Process(object):
	__slots__ = ProcessTable.validcolumns.keys()

ProcessTable.RowType = Process


class ProcessParamsTable(metaio.Table):
	tableName = "process_paramsgroup:process_params:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"param": "lstring",
		"type": "lstring",
		"value": "lstring"
	}

class ProcessParams(object):
	__slots__ = ProcessParamsTable.validcolumns.keys()

ProcessParamsTable.RowType = ProcessParams


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

class SearchSummary(object):
	__slots__ = SearchSummaryTable.validcolumns.keys()

	def get_in(self):
		return segments.segment(lal.LIGOTimeGPS(self.in_start_time, self.in_start_time_ns), lal.LIGOTimeGPS(self.in_end_time, self.in_end_time_ns))

	def set_in(self, seg):
		self.in_start_time, self.in_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.in_end_time, self.in_end_time_ns = seg[1].seconds, seg[1].nanoseconds

	def get_out(self):
		return segments.segment(lal.LIGOTimeGPS(self.out_start_time, self.out_start_time_ns), lal.LIGOTimeGPS(self.out_end_time, self.out_end_time_ns))

	def set_out(self, seg):
		self.out_start_time, self.out_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.out_end_time, self.out_end_time_ns = seg[1].seconds, seg[1].nanoseconds

SearchSummaryTable.RowType = SearchSummary


class SearchSummVarsTable(metaio.Table):
	tableName = "search_summvarsgroup:search_summvars:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"name": "lstring",
		"string": "lstring",
		"value": "real_8"
	}

class SearchSummVars(object):
	__slots__ = SearchSummVarsTable.validcolumns.keys()

SearchSummVarsTable.RowType = SearchSummVars


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

	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

	def set_peak(self, gps):
		self.peak_time, self.peak_time_ns = gps.seconds, gps.nanoseconds

class SnglBurst(object):
	__slots__ = SnglBurstTable.validcolumns.keys()

SnglBurstTable.RowType = SnglBurst


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
		"inpulse_time_ns": "int_4s",
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
		"chisq_dof": "real_4",
		"sigmasq": "real_4",
		"rsqveto_duration": "real_4",
		"event_id": "int_8s"
	}

class SnglInspiral(object):
	__slots__ = SnglInspiralTable.validcolumns.keys()

SnglInspiralTable.RowType = SnglInspiral


class SnglRingDownTable(metaio.Table):
	tableName = "sngl_inspiralgroup:sngl_ringdown:table"
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

class SnglRingDown(object):
	__slots__ = SnglRingDownTable.validcolumns.keys()

SnglRingDownTable.RowType = SnglRingDown


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

class MultiInspiral(object):
	__slots__ = MultiInspiralTable.validcolumns.keys()

MultiInspiralTable.RowType = MultiInspiral


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

class SimInspiral(object):
	__slots__ = SimInspiralTable.validcolumns.keys()

SimInspiralTable.RowType = SimInspiral


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
		"simulation_id": "ilwd_char"
	}

class SimBurst(object):
	__slots__ = SimBurstTable.validcolumns.keys()

SimBurstTable.RowType = SimBurst


class SimRingDownTable(metaio.Table):
	tableName = "sim_inspiralgroup:sim_ringdown:table"
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
		"simulation_id": "ilwd_char"
	}

class SimRingDown(object):
	__slots__ = SimRingDownTable.validcolumns.keys()

SimRingDownTable.RowType = SimRingDown


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

class SummValue(object):
	__slots__ = SummValueTable.validcolumns.keys()

SummValueTable.RowType = SummValue


class SimInstParamsTable(metaio.Table):
	tableName = "sim_inst_paramsgroup:sim_inst_params:table"
	validcolumns = {
		"simulation_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"value": "real_8"
	}

class SimInstParams(object):
	__slots__ = SimInstParamsTable.validcolumns.keys()

SimInstParamsTable.RowType = SimInstParams


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

class Stochastic(object):
	__slots__ = StochasticTable.validcolumns.keys()

StochasticTable.RowType = Stochastic


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

class StochSumm(object):
	__slots__ = StochSummTable.validcolumns.keys()

StochSummTable.RowType = StochSumm


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

class ExtTriggers(object):
	__slots__ = ExtTriggersTable.validcolumns.keys()

ExtTriggersTable.RowType = ExtTriggers


class FilterTable(metaio.Table):
	tableName = "filtergroup:filter:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"program": "lstring",
		"start_time": "int_4s",
		"filter_name": "lstring",
		"comment": "lstring"
	}

class Filter(object):
	__slots__ = FilterTable.validcolumns.keys()

FilterTable.RowType = Filter


TableByName = {
	ProcessTable.tableName: ProcessTable,
	ProcessParamsTable.tableName: ProcessParamsTable,
	SearchSummaryTable.tableName: SearchSummaryTable,
	SearchSummVarsTable.tableName: SearchSummVarsTable,
	SnglBurstTable.tableName: SnglBurstTable,
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
	FilterTable.tableName: FilterTable
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
