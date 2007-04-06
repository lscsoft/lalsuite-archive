# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
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
#                                   Preamble
#
# =============================================================================
#


import sys


from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_bucluster"


def append_process(doc, **kwargs):
	process = llwapp.append_process(doc, program = process_program_name, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	llwapp.append_process_params(doc, process, [("--cluster-algorithm", "lstring", kwargs["cluster_algorithm"])])

	return process


#
# =============================================================================
#
#                        Add "Most Significant" Columns
#
# =============================================================================
#


#
# FIXME:  these columns should be generated by the C code, but that would
# change the sngl_burst table definition and interfere with the string
# search.  Something to sort out later.
#


def add_ms_columns(sngl_burst_table):
	# add columns if required
	added = False
	for colname in ("peak_frequency", "ms_start_time", "ms_start_time_ns", "ms_duration", "ms_flow", "ms_bandwidth", "ms_hrss", "ms_snr", "ms_confidence"):
		try:
			sngl_burst_table.getColumnByName(colname)
		except KeyError:
			sngl_burst_table.appendColumn(colname)
			added = True
	if not added:
		# didn't add any columns, so don't muck their contents
		return

	# at least one column was added, intialize them all
	for row in sngl_burst_table:
		row.peak_frequency = row.central_freq
		row.set_ms_period(row.get_period())
		row.set_ms_band(row.get_band())
		row.ms_hrss = row.amplitude
		row.ms_snr = row.snr
		row.ms_confidence = row.confidence


#
# =============================================================================
#
#                            Clustering Algorithms
#
# =============================================================================
#


#
# "excess power" clustering algorithm
#


def ExcessPowerPreFunc(sngl_burst_table):
	"""
	For speed, convert peak times to floats relative to epoch.
	"""
	if not len(sngl_burst_table):
		return
	offset = sngl_burst_table[0].get_peak()
	for row in sngl_burst_table:
		row.peak_time = float(row.get_peak() - offset)
	return offset



def ExcessPowerPostFunc(sngl_burst_table, offset):
	"""
	Restore peak times to absolute LIGOTimeGPS values.
	"""
	if not len(sngl_burst_table):
		return
	for row in sngl_burst_table:
		row.set_peak(offset + row.peak_time)


def ExcessPowerBailoutFunc(a, b):
	"""
	Orders a and b by ifo, then by channel, then by search, then by
	time interval.  Returns 0 if a and b are from the same channel of
	the same instrument and their time intervals are not disjoint.
	"""
	return cmp(a.ifo, b.ifo) or cmp(a.channel, b.channel) or cmp(a.search, b.search) or llwapp.cmp_seg_intervals(a.get_period(), b.get_period())


def ExcessPowerTestFunc(a, b):
	"""
	Orders a and b by ifo, then by channel, then by search, then time
	interval, then by frequency band.  Returns 0 if a and b are from
	the same channel of the same instrument, and their time-frequency
	tiles are not disjoint.
	"""
	return cmp(a.ifo, b.ifo) or cmp(a.channel, b.channel) or cmp(a.search, b.search) or llwapp.cmp_seg_intervals(a.get_period(), b.get_period()) or llwapp.cmp_seg_intervals(a.get_band(), b.get_band())


def ExcessPowerClusterFunc(a, b):
	"""
	Replace a with a cluster constructed from a and b.  The cluster's
	time-frequency tile is the smallest tile that contains the original
	two tiles, and the "most signficiant" contributor for the cluster
	is the most confident of the two input tiles most significant
	contributors.  The event a is returned.
	"""
	#
	# Save the properties of the most significant contributor
	#

	if b.ms_confidence > a.ms_confidence:
		a.set_ms_period(b.get_ms_period())
		a.set_ms_band(b.get_ms_band())
		a.ms_hrss = b.ms_hrss
		a.ms_snr = b.ms_snr
		a.ms_confidence = b.ms_confidence

	#
	# Compute the SNR-weighted peak time and frequency (recall that the
	# peak times have been converted to floats relative to epoch, and
	# stored in the peak_time column).
	#

	a.peak_time = (a.snr * a.peak_time + b.snr * b.peak_time) / (a.snr + b.snr)
	a.peak_frequency = (a.snr * a.peak_frequency + b.snr * b.peak_frequency) / (a.snr + b.snr)

	#
	# Compute the combined hrss and snr by summing the original ones.
	# Note that no accounting of the overlap of the events is made, so
	# these parameters are being horribly overcounted, but the SNR in
	# particular must be summed like this in order to carry the
	# information needed to continue computing the SNR-weighted peak
	# time and frequencies.
	#

	a.amplitude += b.amplitude
	a.snr += b.snr

	#
	# The confidence is the confidence of the most significant tile.
	# FIXME:  correctly, this should be computed from some sort of
	# joint distribution, but it's not currently used anywhere so
	# there's no point in obsessing over it right now.
	#

	a.confidence = a.ms_confidence

	#
	# The cluster's frequency band is the smallest band containing the
	# bands of the two original events
	#

	a.set_band(llwapp.smallest_enclosing_seg(a.get_band(), b.get_band()))

	#
	# The cluster's time interval is the smallest interval containing
	# the intervals of the two original events
	#

	a.set_period(llwapp.smallest_enclosing_seg(a.get_period(), b.get_period()))

	#
	# Success
	#

	return a


#
# =============================================================================
#
#                               Clustering Loop
#
# =============================================================================
#


def ClusterSnglBurstTable(sngl_burst_table, testfunc, clusterfunc, bailoutfunc = None):
	"""
	Cluster the candidates in the sngl_burst table.  testfunc will be
	passed a pair in random order, and must return 0 (or False) if they
	should be clustered.  clusterfunc will be passed a pair of
	candidates in random order, and must modify the contents of the
	first so as to be a "cluster" of the two.

	If bailoutfunc is not None, the candidates will be sorted into
	"increasing" order using testfunc as a comparison operator, and
	then only pairs of candidates for which bailoutfunc returns 0 (or
	False) will be considered for clustering.  When used this way,
	testfunc must return a numeric result indicating the sort order of
	the two candidates it has been passed:  >0 if the first is
	"greater" than the second, <0 if the first is "less" than the
	second, and 0 if the order does not matter (like a subtraction
	operator).

	The return value is True if the sngl_burst table was modified, and
	False if it was not.
	"""
	# enter loop
	table_changed = False
	did_cluster = True
	while did_cluster:
		did_cluster = False

		if bailoutfunc is not None:
			sngl_burst_table.sort(testfunc)

		# loop over the candidates from the end of the list
		# backwards;  both loops are done in reverse to reduce the
		# memory copies needed when candidates are removed from the
		# list
		for i in xrange(len(sngl_burst_table) - 2, -1, -1):
			# determine the range of candidates to be compared
			# to the current
			# FIXME: make sure all the corner cases are right
			if bailoutfunc is not None:
				for end in xrange(i + 1, len(sngl_burst_table)):
					if bailoutfunc(sngl_burst_table[i], sngl_burst_table[end]):
						break
			else:
				end = len(sngl_burst_table)
			# loop through the comparisons in reverse
			for j in xrange(end - 1, i, -1):
				if not testfunc(sngl_burst_table[i], sngl_burst_table[j]):
					clusterfunc(sngl_burst_table[i], sngl_burst_table.pop(j))
					table_changed = did_cluster = True
	return table_changed


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_bucluster(doc, **kwargs):
	"""
	Run the clustering algorithm on the list of burst candidates.  The
	return value is the tuple (doc, changed), where doc is the input
	document, and changed is a boolean that is True if the contents of
	the sngl_burst table were altered, and False if the triggers were
	not modified by the clustering process.

	If the document does not contain a sngl_burst table, then the
	document is not modified (including no modifications to the process
	metadata tables).
	"""

	#
	# Extract live time segment and sngl_burst table
	#

	try:
		sngl_burst_table = table.get_table(doc, lsctables.SnglBurstTable.tableName)
	except ValueError:
		# no-op:  document does not contain a sngl_burst table
		if kwargs["verbose"]:
			print >>sys.stderr, "document does not contain a sngl_burst table, skipping ..."
		return doc, False
	seg = llwapp.segmentlistdict_fromsearchsummary(doc, program = kwargs["program"]).extent_all()

	# FIXME:  don't do this:  fix lalapps_power's output
	add_ms_columns(sngl_burst_table)

	#
	# Add process information
	#

	process = append_process(doc, **kwargs)

	#
	# Preprocess candidates
	#

	if kwargs["verbose"]:
		print >>sys.stderr, "pre-processing ..."
	preprocess_output = kwargs["prefunc"](sngl_burst_table)

	#
	# Cluster
	#

	if kwargs["verbose"]:
		print >>sys.stderr, "clustering ..."

	table_changed = ClusterSnglBurstTable(sngl_burst_table, kwargs["testfunc"], kwargs["clusterfunc"], kwargs["bailoutfunc"])

	#
	# Postprocess candidates
	#

	if kwargs["verbose"]:
		print >>sys.stderr, "post-processing ..."
	kwargs["postfunc"](sngl_burst_table, preprocess_output)

	#
	# Add search summary information
	#

	llwapp.append_search_summary(doc, process, inseg = seg, outseg = seg, nevents = len(sngl_burst_table))

	#
	# Finish process information
	#

	llwapp.set_process_end_time(process)

	#
	# Done
	#

	return doc, table_changed
