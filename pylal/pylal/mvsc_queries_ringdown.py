try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables
from glue.ligolw import table
from glue.ligolw import ilwd
from glue import segments
from pylal import SnglInspiralUtils
from pylal import db_thinca_rings
from pylal import git_version
from time import clock,time
from optparse import *
import glob
import sys
import random
import math

usage="""
this is a module for use in mvsc_get_doubles
"""

__author__ = "Kari Hodge <khodge@ligo.caltech.edu>"

class CandidateEventQuery:
	# this is the list of parameters that will describe each event in the training and testing sets:
	parameters = "ds_sq delta_t df dQ a_snr b_snr coinc_snr snr_ratio"  #null_stat eff_coh_snr
	# these are the sqlite queries used to extract these parameters (the dimensions to be considered in the multivariate statitical classification algorithm)
	select_dimensions="""
		SELECT
			coinc_ringdown.coinc_event_id,
			snglA.*,
			snglB.*,
			insp_coinc_event.time_slide_id,
			calc_delta_t(snglA.ifo, snglA.start_time, snglA.start_time_ns, snglB.ifo, snglB.start_time, snglB.start_time_ns, insp_coinc_event.time_slide_id),
			abs(snglA.frequency - snglB.frequency),
			abs(snglA.Quality - snglB.Quality),
			snglA.snr,
			snglB.snr,
			coinc_ringdown.snr,
			max(snglA.snr/snglB.snr,snglB.snr/snglA.snr)"""
#			coinc_ringdown.null_stat,
#			coinc_ringdown.eff_coh_snr"""
	add_join_injections="""
		FROM
			coinc_ringdown
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN sngl_ringdown AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_ringdown AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event_map AS mapC ON (mapC.event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapD ON (mapD.coinc_event_id == mapC.coinc_event_id)
			JOIN sim_inspiral ON (sim_inspiral.simulation_id == mapD.event_id)
			JOIN coinc_event AS sim_coinc_event ON (sim_coinc_event.coinc_event_id == mapD.coinc_event_id)
			JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == mapA.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == sim_coinc_event.coinc_def_id)
		WHERE
			( coinc_definer.search == 'ring' OR coinc_definer.search =='ringdown' )
			AND coinc_definer.search_coinc_type == 2
			AND mapA.table_name == 'sngl_ringdown'
			AND mapB.table_name == 'sngl_ringdown'
			AND mapC.table_name == 'coinc_event'
			AND mapD.table_name == 'sim_inspiral'
			AND snglA.ifo == ?
			AND snglB.ifo == ?"""
	add_join_fulldata="""
		, experiment_summary.datatype
		FROM
			coinc_ringdown
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN sngl_ringdown AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_ringdown AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event AS insp_coinc_event ON (mapA.coinc_event_id == insp_coinc_event.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == insp_coinc_event.coinc_def_id)
			JOIN experiment_map ON (experiment_map.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN experiment_summary ON (experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
		WHERE
			coinc_definer.search == 'ring'
			AND coinc_definer.search_coinc_type == 0
			AND mapA.table_name == 'sngl_ringdown'
			AND mapB.table_name == 'sngl_ringdown'
			AND snglA.ifo == ?
			AND snglB.ifo == ?"""
