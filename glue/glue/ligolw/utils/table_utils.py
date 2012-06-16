# Copyright (C) 2012  Matthew West
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


#
# =============================================================================
#
#					Preamble
#
# =============================================================================
#

"""

"""

import itertools
import math
from optparse import OptionParser
import re
import sys
import os
try:
	any
	all
except NameError:
	# Python < 2.5
	from glue.iterutils import any, all

from glue import iterutils
from glue import lal
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import git_version

__author__ = "Matt West <matthew.west@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date



#
# =============================================================================
#
#		     Depopulate ligolw_xml tables 
#
# =============================================================================
#

def depopulate_sngl_inspiral(xmldoc, verbose = False):
	"""
	This function takes the lists of event_ids from the sngl_inspiral and coinc_event_map 
	tables and determine the difference, such that the newlist contains only  non-coinc 
	single-ifo triggers. Then it remove these non-coinc triggers from the sngl_inspiral table.
	"""
	sngls_tbl = lsctables.table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
	sngls_tbl_eid = list(sngls_tbl.getColumnByName("event_id"))

	coinc_map_tbl = lsctables.table.get_table(xmldoc, lsctables.CoincMapTable.tableName)

	if len(coinc_map_tbl) == 0:
		del sngls_tbl[:]
		if verbose:
			 print >> sys.stderr, "This file lacks any coincident events. All %i single-ifo inspiral triggers have been removed." %( len(sngls_tbl_eid) )
	else:
		coinc_map_tbl_eid = list(coinc_map_tbl.getColumnByName("event_id"))
		non_coincs = list(set(sngls_tbl_eid) - set(coinc_map_tbl_eid))

		coinc_sngls_tbl = xmldoc.childNodes[0].insertBefore( lsctables.New(lsctables.SnglInspiralTable), sngls_tbl)
		for (idx,event_id) in enumerate(coinc_map_tbl_eid):
			coinc_sngls_tbl.insert(idx, sngls_tbl[sngls_tbl_eid.index(event_id)] )
		xmldoc.childNodes[0].removeChild(sngls_tbl)

		if verbose:
			print >> sys.stderr, "%i single-ifo inspiral triggers not associated with a coincident event have been removed." %( len(non_coincs) )

	return xmldoc


def depopulate_experiment_tables(xmldoc, verbose = False):
	"""
	Removes entries from the experiment tables that do not have events
	or durations in them. In other words, if none of the rows in the
	experiment_summary table that are assoicated with a single experiment_id
	have neither an event in them nor any duration then all of the rows in the
	experiment_summary table associated with that experiment_id are deleted,
	as well as the corresponding row in the experiment table. If, however, just
	one of the rows in the experiment_summary table associated with an experiment_id
	have at least a 1 in their nevents column or at least 1 second in the
	durations column then nothing associated with that experiment_id are deleted.
	(In other words, if just one time slide in an experiment has just 1 event or
	is just 1 second long, then none of the time slides in that experiment are deleted,
	even if all of the other time slides in that experiment have nothing in them.)
	"""

	if verbose:
		print >>sys.stderr, "Depopulating the experiment tables...",

	#
	# find the experiment and experiment summary table
	#

	try:
		experiment_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
	except ValueError:
		# no table --> no-op
		if verbose:
			print >>sys.stderr, "Cannot find the experiment table"
		return

	try:
		experiment_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)
	except ValueError:
		# no table --> no-op
		if verbose:
			print >>sys.stderr, "Cannot find the experiment_summary table"
		return

	del_eid_indices = []
	del_esid_indices = []

	for mm, erow in enumerate(experiment_table):
		this_eid = erow.experiment_id
		es_index_list = []
		for nn, esrow in enumerate(experiment_summ_table):
			if esrow.experiment_id == this_eid and (esrow.duration or esrow.nevents):
				# something in this experiment, go on to next experiment_id
				break
			if esrow.experiment_id == this_eid:
				es_index_list.append(nn)
			if nn == len(experiment_summ_table) - 1:
				# if get to here, nothing in that experiment, mark id and indices
				# for removal
				del_eid_indices.append(mm)
				del_esid_indices += es_index_list

	# delte all experiments who's eids fall in del_eid_indices
	del_eid_indices.sort(reverse = True)
	for eid_index in del_eid_indices:
		del experiment_table[eid_index]
	# delete all experiment_summaries whose esids fall in del_esid_indices
	del_esid_indices.sort(reverse = True)
	for esid_index in del_esid_indices:
		del experiment_summ_table[esid_index]

	if verbose:
		print >> sys.stderr, "removed %i empty experiment(s) from the experiment table and %i associated time slides from the experiment_summary table." %( len(del_eid_indices), len(del_esid_indices) )


#
# =============================================================================
#
#		experiment and experiment_summ tables 
#
# =============================================================================
#

def get_experiment_times(xmldoc):
	"""
	Use the start & end-times stored in the segment_summary table to define 
	the experiment times.  This presumes that the vetoes file has been added
	to the file being analyzed and that the program used to make said vetoes file
	is ligolw_segments_from_cats.
	"""
	if ".executable/ligolw_segments_from_cats" in process_tbl.getColumnByName("program")
		# get the segment_summary table
		segsum_tbl = lsctables.table.get_table(xmldoc, lsctables.SegmentSummaryTable.tableName)
		expr_start_time = []
		expr_end_time = []
		for row segment_summary_tbl:
			start_time.append(row.start_time)
			end_time.append(row.end_time)
		expr_start_time = min(expr_start_time)
		expr_end_time = max(expr_end_time)
	else:
		# if the segments tables are not in the file, set these times to None
		expr_start_time = None
		expr_end_time = None

	return expr_start_time, expr_end_time

def populate_experiment_table(
	xmldoc,
	search_group,
	trigger_program,
	lars_id,
	instruments,
	comments = None,
	add_inst_subsets = False,
	verbose = False
):
	"""
	Populate the experiment table using the given entries. If
	add_inst_subsets is set to True, will write additional
	entries for every possible sub-combination of the given instrument
	set. Returns a dictionary of experiment_ids keyed by the instrument
	set.

	@xmldoc: xmldoc to get/write table to
	@lars_id: lars_id of the experiment
	@search_group: lsc group that performed the experiment (e.g., cbc)
	@trigger_program: name of the program that performed the analysis
		(e.g., inspiral, ringdown, etc.)
	@comments: any desired comments
	@add_inst_subsets: will write an entry for every possible subset
		of @instruments
	@verbose: be verbose
	"""

	if verbose:
		print >> sys.stderr, "\tPopulating the Experiment table..."

	# find the experiment table or create one if needed
	try:
		expr_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
	except ValueError:
		expr_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentTable))

	# determine experiment start and end times
	expr_start_time, expr_end_time = get_experiment_times(xmldoc)

	# write entry to the experiment table for the given instruments if it doesn't already exist
	experiment_ids = {}

	experiment_ids[frozenset(instruments)] =  expr_table.write_new_expr_id(
		search_group,
		trigger_program,
		lars_id,
		instruments,
		expr_start_time,
		expr_end_time,
		comments = comments
	)

	#  add every possible sub-combination of the instrument set if
	# they're not already in the table
	if add_inst_subsets:
		for nn in range(2, len(instruments) ):
			for sub_combo in iterutils.choices( list(instruments), nn ):
				if frozenset(sub_combo) not in experiment_ids:
					experiment_ids[frozenset(sub_combo)] = expr_table.write_new_expr_id(
						search_group,
						trigger_program,
						lars_id,
						sub_combo,
						expr_start_time,
						expr_end_time,
						comments = comments
					)
	return experiment_ids

def get_experiment_type(xmldoc, time_slide_dict):
	"""
	Determine the which experiment type(s) the coincident triggers in this
	file belong to.  It uses information from the inspiral files stored in
	the process params table to decide if the triggers come from playground
	time or are from an injection run.  If the time_slide_dict has more than
	one entry, then the triggers are from a slide run.
	"""
	# get the param column of the process_params table
	process_params_tbl = lsctables.table.get_table(xmldoc, lsctables.ProcessParamsTable.tableName)
	pp_value = set(process_params_tbl.getColumnByName("value"))

	process_tbl = lsctables.table.get_table(xmldoc, lsctables.ProcessTable.tableName)
	program = set(process_tbl.getColumnByName("program"))

	# determine experiment type(s)
	if len(time_slide_dict) == 1:
		if 'inspinj' in program:
			datatypes = ['simulation']
			sim_proc_id = str(process_tbl.get_ids_by_program("inspinj").pop())
		elif 'rinj' in program:
			datatypes = ['simulation']
			sim_proc_id = str(process_tbl.get_ids_by_program("rinj").pop())
		else:
			datatypes = ['all_data']
			sim_proc_id = None
			if 'PLAYGROUND' in pp_value:
				datatypes += ['playground']
			else:
				datatypes += ['exclude_play']
	elif len(time_slide_dict) > 1:
		if 'PLAYGROUND' in pp_value:
			datatypes = ['play_slide']
		else:
			datatypes = ['slide']
		sim_proc_id = None

	return datatypes, sim_proc_id


def populate_experiment_summ_table(
	xmldoc,
	experiment_id,
	time_slide_dict,
	veto_def_name,
	return_dict = False,
	verbose = False
):
	"""
	Populate the experiment_summ_table using an experiment_id, a
	veto_def_name, and a list of time_slide ids.

	@xmldoc: xmldoc to get/write table to
	@experiment_id: experiment_id to be added to the table.
	@veto_def_name: veto_def_name to be added to the table.
	@time_slide_dict: time_slide table as dictionary; used to set time_slide_id
		column and figure out whether or not is zero-lag. Can either be the result
		of lsctables.time_slide_table.as_dict or any dictionary having same format.
	@return_dict: will return the experiment_summary table as an id_dict
	"""

	if verbose:
		print >> sys.stderr, "\tPopulating the Experiment Summary table..."

	# find the experiment_summary table or create one if needed
	try:
		expr_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)
	except ValueError:
		expr_summ_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentSummaryTable))

	# populate the experiment_summary table
	datatypes, sim_proc_id = get_experiment_type(xmldoc, time_slide_dict)

	for type in datatypes:
		for slide_id in time_slide_dict:
			expr_summ_table.write_experiment_summ(
				experiment_id,
				slide_id,
				veto_def_name,
				type,
				sim_proc_id = sim_proc_id
			)


def get_on_instruments(xmldoc, trigger_program):
	process_tbl = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
	instruments = set([])
	for row in process_tbl:
		if row.program == trigger_program:
			instruments.add(row.ifos)
	return instruments


def generate_experiment_tables(xmldoc, **cmdline_opts):
	"""
	Create or adds entries to the experiment table and experiment_summ
	table using instruments pulled from the search summary table and
	offsets pulled from the time_slide table.
	"""

	if cmdline_opts["verbose"]:
		print >> sys.stderr, "Populating the experiment and experiment_summary tables using search_summary and time_slide tables..."

	# Get the instruments that were on
	instruments = get_on_instruments(xmldoc, cmdline_opts["trigger_program"])

	# Populate the experiment table
	experiment_ids = populate_experiment_table(
		xmldoc,
		cmdline_opts["search_group"],
		cmdline_opts["trigger_program"],
		cmdline_opts["lars_id"],
		instruments,
		comments = cmdline_opts["comment"],
		add_inst_subsets = True,
		verbose = cmdline_opts["verbose"]
	)

	# Get the time_slide table as dict
	time_slide_dict = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict()

	# Populate the experiment_summary table
	for instruments in experiment_ids:
		populate_experiment_summ_table(
			xmldoc,
			experiment_ids[instruments],
			time_slide_dict,
			cmdline_opts["vetoes_name"],
			return_dict = False,
			verbose = cmdline_opts["verbose"]
		)


def populate_experiment_map(xmldoc, veto_def_name, verbose = False):
	from glue.pipeline import s2play as is_in_playground

	#
	# find the experiment_map table or create one if needed
	#
	if verbose:
		print >> sys.stderr, "\tMapping coinc events to experiment_summary table..."

	try:
		expr_map_table = table.get_table(xmldoc, lsctables.ExperimentMapTable.tableName)
	except ValueError:
		expr_map_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentMapTable))

	#
	# find the coinc_event table
	#

	coinc_event_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)

	#
	# Index the coinc_inspiral table as a dictionary
	#

	coinc_index = dict((row.coinc_event_id, row) for row in table.get_table(xmldoc, lsctables.CoincInspiralTable.tableName))

	#
	# Get the time_slide_table as dict
	#

	time_slide_dict = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict()

	#
	# find the experiment & experiment summary tables
	#

	expr_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
	expr_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)

	#
	# determine what experiment datatype this file belongs to
	#

	datatypes, sim_proc_id = get_experiment_type(xmldoc, time_slide_dict)

	#
	# cycle through the coincs in the coinc_inspiral table
	#
	for coinc in coinc_event_table:

		#
		# get the experiment and experiment_summ_id for this coinc
		#

		for expr in expr_table:
			if expr.instruments == coinc.instruments:
				expr_id =  expr.experiment_id

		for type in datatypes:
			# map the coinc to an experiment
			expr_map = lsctables.ExperimentMap()
			expr_map.coinc_event_id = coinc.coinc_event_id

			expr_map.experiment_summ_id = expr_summ_table.get_expr_summ_id(
				expr_id,
				coinc.time_slide_id,
				veto_def_name,
				type,
				sim_proc_id = sim_proc_id
			)
			if not expr_map.experiment_summ_id:
				raise ValueError, "%s experiment_summ_id could not be found with %s" \
				%( type, ','.join([ str(expr_id), str(coinc.time_slide_id), veto_def_name ]))

			# map the experiment
			expr_map_table.append(expr_map)
			# Increment number of events in nevents column by 1
			expr_summ_table.add_nevents( expr_map.experiment_summ_id, 1 )

def make_experiment_tables(xmldoc, **cmdline_opts):
	"""
	This collection of functions creates and populates the experiment,
	experiment_summary, and experiment_map tables.	
	"""

	generate_experiment_tables(xmldoc, **cmdline_opts)

	populate_experiment_map(xmldoc, cmdline_opts["vetoes_name"], verbose = cmdline_opts["verbose"])

	depopulate_experiment_tables(xmldoc, verbose = cmdline_opts["verbose"])

