# $Id$
#
# Copyright (C) 2007  Kipp C. Cannon
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


import math
import numpy
from scipy.interpolate import interpolate
import sys


from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import itertools
from pylal import ligolw_burca_tailor
from pylal import llwapp
from pylal.date import LIGOTimeGPS


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                  Likelihood
#
# =============================================================================
#


#
# starting from Bayes' theorem:
#
# P(coinc is a g.w. | its parameters)
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#                                P(parameters)
#
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * P(coinc is not g.w.) + P(inj params) * P(coinc is g.w.)
#
#                       P(inj params) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * [1 - P(coinc is g.w.)] + P(inj params) * P(coinc is g.w.)
#
#                        P(inj params) * P(coinc is g.w.)
#   = ----------------------------------------------------------------------
#     P(noise params) + [P(inj params) - P(noise params)] * P(coinc is g.w.)
#


def make_interps(pair, x, y):
	# extrapolate x and y arrays by one element at each end.  this has
	# to be done because the Rate class in pylal.rate returns the x
	# co-ordinates as the bin centres, which is correct, but it means
	# that an event can have a set of parameter values that lie beyond
	# the end of the x co-ordinate array (the parameters are still in
	# the bin, but in the outer half), meanwhile the scipy interpolator
	# insists on only interpolating (go figger) so it throws an error
	# when such an event is encountered.
	x = numpy.hstack((x[0] + (x[0] - x[1]), x, x[-1] + (x[-1] - x[-2])))
	y = numpy.hstack((y[0] + (y[0] - y[1]), y, y[-1] + (y[-1] - y[-2])))

	# construct interpolator
	interp_dict = {}
	interp_dict[pair] = interpolate.interp1d(x, y)

	# construct reversed interpolator (for when the instrument pair is
	# reversed)
	interp_dict[(pair[1], pair[0])] = interpolate.interp1d(-x[::-1], y[::-1])

	# done
	return interp_dict


class Likelihood(object):
	def __init__(self):
		self.dt_bak = {}
		self.dt_inj = {}
		self.df_bak = {}
		self.df_inj = {}
		self.dh_bak = {}
		self.dh_inj = {}
		self.P_gw = None

	def set_dt_from_xml(self, pair, xml):
		self.dt_bak.update(make_interps(pair, xml.array[0], xml.array[2]))
		self.dt_inj.update(make_interps(pair, xml.array[0], xml.array[1]))

	def set_df_from_xml(self, pair, xml):
		self.df_bak.update(make_interps(pair, xml.array[0], xml.array[2]))
		self.df_inj.update(make_interps(pair, xml.array[0], xml.array[1]))

	def set_dh_from_xml(self, pair, xml):
		self.dh_bak.update(make_interps(pair, xml.array[0], xml.array[2]))
		self.dh_inj.update(make_interps(pair, xml.array[0], xml.array[1]))

	def set_P_gw(self, P):
		self.P_gw = P

	def P(self, events, offsets):
		P_inj = 1.0
		P_bak = 1.0
		for event1, event2 in itertools.choices(events, 2):
			ifos = (event1.ifo, event2.ifo)
			if ifos not in self.dt_inj:
				# no statistics for this instrument pair
				continue

			t1 = LIGOTimeGPS(event1.peak_time, event1.peak_time_ns) + offsets[event1.ifo]
			t2 = LIGOTimeGPS(event2.peak_time, event2.peak_time_ns) + offsets[event2.ifo]
			dt = float(t1 - t2) / ((event1.ms_duration + event2.ms_duration) / 2)

			df = (event1.peak_frequency - event2.peak_frequency) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)

			dh = (event1.ms_hrss - event2.ms_hrss) / ((event1.ms_hrss + event2.ms_hrss) / 2)

			P_inj *= self.dt_inj[ifos](dt)[0] * self.df_inj[ifos](df)[0] * self.dh_inj[ifos](dh)[0]
			P_bak *= self.dt_bak[ifos](dt)[0] * self.df_bak[ifos](df)[0] * self.dh_bak[ifos](dh)[0]
		return P_bak, P_inj

	def __call__(self, events, offsets):
		"""
		Compute the likelihood that the coincident n-tuple of
		events are the result of a gravitational wave:  the
		probability that the hypothesis "the events are a
		gravitational wave" is correct, in the context of the
		measured background and foreground distributions, and the
		intrinsic rate of gravitational wave coincidences.  offsets
		is a dictionary of instrument --> offset mappings to be
		used to time shift the events before comparison.
		"""
		P_bak, P_inj = self.P(events, offsets)
		return (P_inj * self.P_gw) / (P_bak + (P_inj - P_bak) * self.P_gw)


class Confidence(Likelihood):
	def __call__(self, events, offsets):
		"""
		Compute the confidence that the list of events are the
		result of a gravitational wave:  -ln[1 - P(gw)], where
		P(gw) is the likelihood returned by the Likelihood class.
		A set of events very much like gravitational waves will
		have a likelihood of being a gravitational wave very close
		to 1, so 1 - P is a small positive number, and so -ln of
		that is a large positive number.
		"""
		P_bak, P_inj = self.P(events, offsets)
		return  math.log(P_bak + (P_inj - P_bak) * self.P_gw) - math.log(P_inj) - math.log(self.P_gw)


#
# =============================================================================
#
#                         Load Distribution Functions
#
# =============================================================================
#


#
# Load likelihood distribution functions.  This must be kept synchronized
# with the document generator in pylal.ligolw_burca_tailor.
#


def load_likelihood_control(filename, verbose = False):
	#
	# Load the XML document
	#

	xmldoc = utils.load_filename(filename, verbose = verbose, gz = filename[-3:] == ".gz")

	#
	# Find the first LIGO_LW element with the name
	# "ligolw_burca_tailor"
	#

	subtree = None
	for node in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
		try:
			if node.getAttribute("Name") == ligolw_burca_tailor.process_program_name:
				subtree = node
				break
		except KeyError:
			# doesn't have a Name attribute, so can't be it
			pass
	if subtree is None:
		raise ValueError, "%s does not contain likelihood information" % filename

	#
	# Extract the pickled thresholds object
	#

	thresholds = llwapp.pickle_from_param(subtree, "thresholds")

	#
	# Extract the likelihood arrays
	#

	likelihood = Likelihood()
	for pair in thresholds.keys():
		rev = (pair[1], pair[0])

		try:
			a = array.get_array(subtree, "%s_%s_dt" % pair)
			likelihood.set_dt_from_xml(pair, a)
		except ValueError:
			try:
				a = array.get_array(subtree, "%s_%s_dt" % rev)
				likelihood.set_dt_from_xml(rev, a)
			except ValueError:
				raise ValueError, "no delta t distribution information for instrument pair %s in %s" % (pair, filename)

		try:
			a = array.get_array(subtree, "%s_%s_df" % pair)
			likelihood.set_df_from_xml(pair, a)
		except ValueError:
			try:
				a = array.get_array(subtree, "%s_%s_df" % rev)
				likelihood.set_df_from_xml(rev, a)
			except ValueError:
				raise ValueError, "no delta f distribution information for instrument pair %s in %s" % (pair, filename)

		try:
			a = array.get_array(subtree, "%s_%s_dh" % pair)
			likelihood.set_dh_from_xml(pair, a)
		except ValueError:
			try:
				a = array.get_array(subtree, "%s_%s_dh" % rev)
				likelihood.set_dh_from_xml(rev, a)
			except ValueError:
				raise ValueError, "no delta h_rss distribution information for instrument pair %s in %s" % (pair, filename)
	
	return likelihood


#
# =============================================================================
#
#                              Library Interface
#
# =============================================================================
#


def ligolw_burca2(xmldoc, likelihood, verbose = False):
	"""
	Assigns likelihood values to excess power coincidences.  xmldoc is
	an XML document tree to process, and likelihood is a Likelihood
	class instance initialized from a likelihood control file (for
	example as obtained by calling load_likelihood_control()).
	"""
	#
	# Find document parts.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."

	definer_ids = [llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName], create_new = False)]
	try:
		definer_ids.append(llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName], create_new = False))
	except KeyError:
		# guess there are no injections in this file?
		pass
	time_slides = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).get_offsets()
	zero_lag_tisi_ids = llwapp.get_zero_lag_time_slides(xmldoc).keys()
	coinc_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)

	#
	# Index the sngl_burst table.
	#

	index = {}
	for row in table.get_table(xmldoc, lsctables.SnglBurstTable.tableName):
		index[row.event_id] = row

	#
	# Index the coinc_event_map table.
	#

	for row in table.get_table(xmldoc, lsctables.CoincMapTable.tableName):
		if row.table_name == table.StripTableName(lsctables.SnglBurstTable.tableName):
			if row.coinc_event_id not in index:
				index[row.coinc_event_id] = []
			index[row.coinc_event_id].append(index[row.event_id])

	#
	# Iterate over coincs, assigning likelihoods to burst+burst coincs,
	# and sim+burst coincs if the document contains them.
	#

	if verbose:
		print >>sys.stderr, "computing likelihoods ..."
		n_coincs = len(coinc_table)

	for n, coinc in enumerate(coinc_table):
		if verbose and not n % 30:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / n_coincs),
		if coinc.coinc_def_id in definer_ids:
			coinc.likelihood = likelihood(index[coinc.coinc_event_id], time_slides[coinc.time_slide_id])
	if verbose:
		print >>sys.stderr, "\t100.0%"

	#
	# Done
	#

	return xmldoc
