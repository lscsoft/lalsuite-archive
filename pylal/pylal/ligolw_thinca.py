# $Id$
#
# Copyright (C) 2008  Kipp C. Cannon
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


import bisect
import math
import sys


from glue import iterutils
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from pylal import llwapp
from pylal import snglcoinc
from pylal.date import LIGOTimeGPS
from pylal import tools
from pylal.xlal import tools as xlaltools
try:
	all
except NameError:
	# Python < 2.5.x
	from glue.iterutils import all as all


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


#
# Use C row type for coinc_event_map table
#


lsctables.CoincMapTable.RowType = lsctables.CoincMap = xlaltools.CoincMap


#
# Construct a subclass of the C sngl_inspiral row class with the methods
# that are needed
#


class SnglInspiral(xlaltools.SnglInspiralTable):
	__slots__ = ()

	def get_end(self):
		return LIGOTimeGPS(self.end_time, self.end_time_ns)

	def set_end(self, gps):
		self.end_time, self.end_time_ns = gps.seconds, gps.nanoseconds

	def get_effective_snr(self, fac):
		return self.snr/ (1 + self.snr**2/fac)**(0.25)/(self.chisq/(2*self.chisq_dof - 2) )**(0.25)

	def __eq__(self, other):
		return not (
			cmp(self.ifo, other.ifo) or
			cmp(self.end_time, other.end_time) or
			cmp(self.end_time_ns, other.end_time_ns) or
			cmp(self.mass1, other.mass1) or
			cmp(self.mass2, other.mass2) or
			cmp(self.search, other.search)
		)

	def __cmp__(self, other):
		# compare self's end time to the LIGOTimeGPS instance
		# other.  allows bisection searches by GPS time to find
		# ranges of triggers quickly
		return cmp(self.end_time, other.seconds) or cmp(self.end_time_ns, other.nanoseconds)


#
# Use C LIGOTimeGPS type
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


#
# Use C segments module
#


def use___segments(modulename):
	from glue import __segments
	modulename.segments.infinity = __segments.infinity
	modulename.segments.NegInfinity = __segments.NegInfinity
	modulename.segments.PosInfinity = __segments.PosInfinity
	modulename.segments.segment = __segments.segment
	modulename.segments.segmentlist = __segments.segmentlist
use___segments(llwapp)
use___segments(lsctables)


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_thinca"


def append_process(xmldoc, comment = None, force = None, program = None, e_thinca_parameter = None, verbose = None):
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [
		(u"--program", u"lstring", program),
		(u"--e-thinca-parameter", u"real_8", e_thinca_parameter)
	]
	if comment is not None:
		params += [(u"--comment", u"lstring", comment)]
	if force is not None:
		params += [(u"--force", None, None)]
	if verbose is not None:
		params += [(u"--verbose", None, None)]

	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# The sngl_inspiral <--> sngl_inspiral coinc type.
#


InspiralCoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 0, description = u"sngl_inspiral<-->sngl_inspiral coincidences")


#
# Custom snglcoinc.CoincTables subclass.
#


class InspiralCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc):
		snglcoinc.CoincTables.__init__(self, xmldoc)

		#
		# find the coinc_inspiral table or create one if not found
		#

		try:
			self.coinc_inspiral_table = lsctables.table.get_table(xmldoc, lsctables.CoincInspiralTable.tableName)
		except ValueError:
			self.coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
			xmldoc.childNodes[0].appendChild(self.coinc_inspiral_table)

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events, effective_snr_factor):
		#
		# populate the coinc_event and coinc_event_map tables
		#

		coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, coinc_def_id, events)

		#
		# populate the coinc_inspiral table:
		#
		# - end_time is the end time of the first trigger in
		#   alphabetical order by instrument (!?) time-shifted
		#   according to the coinc's offset vector
		# - mass is average of total masses
		# - mchirp is average of mchirps
		# - snr is root-sum-square of SNRs
		# - false-alarm rates are blank
		#

		events = sorted(events, lambda a, b: cmp(a.ifo, b.ifo))

		coinc_inspiral = self.coinc_inspiral_table.RowType()
		coinc_inspiral.coinc_event_id = coinc.coinc_event_id
		coinc_inspiral.mass = sum(event.mass1 + event.mass2 for event in events) / len(events)
		coinc_inspiral.mchirp = sum(event.mchirp for event in events) / len(events)
		if all(event.chisq for event in events):
			coinc_inspiral.snr = math.sqrt(sum(event.get_effective_snr(fac = effective_snr_factor)**2 for event in events))
		else:
			# would get divide-by-zero without a \chi^{2} value
			coinc_inspiral.snr = None
		coinc_inspiral.false_alarm_rate = None
		coinc_inspiral.combined_far = None
		coinc_inspiral.set_end(events[0].get_end())
		coinc_inspiral.set_ifos(event.ifo for event in events)
		self.coinc_inspiral_table.append(coinc_inspiral)

		return coinc


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#


class InspiralEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the inspiral
	search.
	"""
	def make_index(self):
		"""
		Sort events by end time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglInspiral row class having
		previously been set to compare the event's end time to a
		LIGOTimeGPS.
		"""
		self.sort(lambda a, b: cmp(a.end_time, b.end_time) or cmp(a.end_time_ns, b.end_time_ns))

	def set_dt(self, dt):
		"""
		If an event's end time differs by more than this many
		seconds from the end time of another event then it is
		*impossible* for them to be coincident.
		"""
		# add 1% for safety, and pre-convert to LIGOTimeGPS to
		# avoid doing type conversion in loops
		self.dt = LIGOTimeGPS(dt * 1.01)

	def _add_offset(self, delta):
		"""
		Add an amount to the end time of each event.
		"""
		for event in self:
			event.set_end(event.get_end() + delta)

	def get_coincs(self, event_a, e_thinca_parameter, comparefunc):
		#
		# event_a's end time
		#

		end = event_a.get_end()

		#
		# extract the subset of events from this list that pass
		# coincidence with event_a (use bisection searches for the
		# minimum and maximum allowed end times to quickly identify
		# a subset of the full list)
		#

		return [event_b for event_b in self[bisect.bisect_left(self, end - self.dt) : bisect.bisect_right(self, end + self.dt)] if not comparefunc(event_a, event_b, e_thinca_parameter)]


#
# =============================================================================
#
#                              Coincidence Tests
#
# =============================================================================
#


def inspiral_max_dt(events, e_thinca_parameter):
	"""
	Given an e-thinca parameter and a list of sngl_inspiral events,
	return the greatest \Delta t that can separate two events and they
	still be considered coincident.
	"""
	# for each instrument present in the event list, compute the
	# largest \Delta t interval for the events from that instrument,
	# and return the sum of the largest two such \Delta t's.
	return sum(sorted(max(xlaltools.XLALSnglInspiralTimeError(event, e_thinca_parameter) for event in events if event.ifo == instrument) for instrument in set(event.ifo for event in events))[-2:])


def inspiral_coinc_compare(a, b, e_thinca_parameter):
	"""
	Returns False (a & b are coincident) if they pass the ellipsoidal
	thinca test.
	"""
	try:
		# FIXME:  should it be ">" or ">="?
		return xlaltools.XLALCalculateEThincaParameter(a, b) > e_thinca_parameter
	except ValueError:
		# ethinca test failed to converge == events are not
		# coincident
		return True


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def replicate_threshold(e_thinca_parameter, instruments):
	"""
	From a single threshold and a list of instruments, return a
	dictionary whose keys are every instrument pair (both orders), and
	whose values are all the same single threshold.

	Example:

	>>> replicate_threshold(6, ["H1", "H2"])
	{("H1", "H2"): 6, ("H2", "H1"): 6}
	"""
	instruments = sorted(instruments)
	thresholds = dict((pair, e_thinca_parameter) for pair in iterutils.choices(instruments, 2))
	instruments.reverse()
	thresholds.update(dict((pair, e_thinca_parameter) for pair in iterutils.choices(instruments, 2)))
	return thresholds


def ligolw_thinca(
	xmldoc,
	program,
	process_id,
	EventListType,
	CoincTables,
	coinc_definer_row,
	event_comparefunc,
	thresholds,
	ntuple_comparefunc = lambda events: False,
	get_max_segment_gap = lambda xmldoc, thresholds: float("inf"),
	effective_snr_factor = 250.0,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."
	coinc_tables = CoincTables(xmldoc)
	coinc_def_id = llwapp.get_coinc_def_id(xmldoc, coinc_definer_row.search, coinc_definer_row.search_coinc_type, create_new = True, description = coinc_definer_row.description)

	#
	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence
	#

	eventlists = snglcoinc.make_eventlists(xmldoc, EventListType, lsctables.SnglInspiralTable.tableName, get_max_segment_gap(xmldoc, thresholds), program)

	#
	# set the \Delta t parameter on all the event lists
	#

	max_dt = inspiral_max_dt(lsctables.table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName), thresholds)
	if verbose:
		print >>sys.stderr, "event bisection search window will be %.16g s" % max_dt
	for eventlist in eventlists.values():
		eventlist.set_dt(max_dt)

	#
	# replicate the ethinca parameter for every possible instrument
	# pair
	#

	avail_instruments = set(eventlists.keys())
	thresholds = replicate_threshold(thresholds, avail_instruments)

	#
	# iterate over time slides
	#

	time_slides = coinc_tables.get_ordered_time_slides()
	for n, (time_slide_id, offsetdict) in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "time slide %d/%d: %s" % (n + 1, len(time_slides), ", ".join(("%s = %+.16g s" % x) for x in sorted(offsetdict.items())))

		#
		# can we do it?
		#

		offset_instruments = set(offsetdict.keys())
		if len(offset_instruments) < 2:
			if verbose:
				print >>sys.stderr, "\tsingle-instrument time slide: skipped"
			continue
		if not offset_instruments.issubset(avail_instruments):
			if verbose:
				print >>sys.stderr, "\twarning: do not have data for instrument(s) %s: skipping" % ", ".join(offset_instruments - avail_instruments)
			continue

		#
		# apply offsets to events
		#

		if verbose:
			print >>sys.stderr, "\tapplying time offsets ..."
		eventlists.set_offsetdict(offsetdict)

		#
		# search for and record coincidences
		#

		if verbose:
			print >>sys.stderr, "\tsearching ..."
		for ntuple in snglcoinc.CoincidentNTuples(eventlists, event_comparefunc, offset_instruments, thresholds, verbose = verbose):
			if not ntuple_comparefunc(ntuple):
				coinc_tables.append_coinc(process_id, time_slide_id, coinc_def_id, ntuple, effective_snr_factor)

	#
	# remove time offsets from events
	#

	eventlists.remove_offsetdict()

	#
	# done
	#

	return xmldoc
