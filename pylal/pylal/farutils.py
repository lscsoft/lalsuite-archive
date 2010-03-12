# Copyright (C) 2010 Chad Hanna
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

from glue import iterutils
from glue import segments
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from pylal import SnglBurstUtils
from pylal import llwapp
from pylal import db_thinca_rings

# get choices from a set (useful for on/off ifos)
def detector_combos( ifo_set ):
	out = []
	for i in range(1, len(ifo_set)):
		X = list(iterutils.choices(ifo_set, i))
		Y = list(iterutils.choices(ifo_set, len(ifo_set) - i))
		Y.reverse()
		out.extend(zip(X, Y))
	ifo_set = list(ifo_set)
	ifo_set.sort()
	ifo_set = tuple(ifo_set)
	out.append((ifo_set, ()))
	return out

def add_livetime_nonring(connection, live_time_program, verbose = False):
	# get the segment lists and live time
	# FIXME veto segments not handled yet
	xmldoc = dbtables.get_xml(connection)
	zero_lag_time_slides, background_time_slides = SnglBurstUtils.get_time_slides(connection)
	seglists = llwapp.segmentlistdict_fromsearchsummary(xmldoc, live_time_program).coalesce()
	instruments = frozenset(seglists.keys())
	cached_livetime = {}
	for on_inst, off_inst in detector_combos(list(instruments)):
		on_inst = frozenset(on_inst)
		off_inst = frozenset(off_inst)
		for time_slide in background_time_slides.values():
			seglists.offsets.update(time_slide)
			segs=seglists.intersection(list(on_inst))-seglists.union(list(off_inst))
			key = lsctables.ifos_from_instrument_set(on_inst)
			cached_livetime.setdefault(key, 0)
			cached_livetime[key] += float(abs(segs))
	return cached_livetime

def add_livetime_ring(connection, live_time_program, veto_segments_name=None, verbose = False):
	cached_livetime = {}
	if veto_segments_name is not None:
		if verbose:
			print >>sys.stderr, "\tretrieving veto segments \"%s\" ..." % veto_segments_name
			veto_segments = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
		else:
			veto_segments = segments.segmentlistdict()
		if verbose:
			print >>sys.stderr, "\tcomputing livetimes:",
		for on_instruments, livetimes in db_thinca_rings.get_thinca_livetimes(db_thinca_rings.get_thinca_rings_by_available_instruments(connection, program_name = live_time_program), veto_segments, db_thinca_rings.get_background_offset_vectors(connection), verbose = verbose).items():
			on_instruments = lsctables.ifos_from_instrument_set(on_instruments)
			try:
				cached_livetime[on_instruments] += sum(livetimes)
			except KeyError:
				cached_livetime[on_instruments] = sum(livetimes)
		if verbose:
			print >>sys.stderr
	return cached_livetime
