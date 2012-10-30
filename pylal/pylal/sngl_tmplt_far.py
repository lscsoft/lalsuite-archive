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
#				 Preamble
#
# =============================================================================
#

"""
A collection of utilities that allow one to estimate a false alarm rate from 
triggers of a single template without the need for doing time-slides.
"""

import sqlite3
import numpy
import math

from glue.ligolw import dbtables
from glue import segments

from pylal import ligolw_sqlutils
from pylal import ligolw_compute_durations as compute_dur



#
# =============================================================================
#
#			       Library API 
#
# =============================================================================
#

def get_newsnr(snr, chisq, chisq_dof):
	fac = 6.0
	rchisq = chisq/(2*chisq_dof-2)
	newsnr = snr/(0.5*(1+rchisq**(fac/2.0)))**(1.0/fac)
	if rchisq < 1:
		newsnr = snr
	return newsnr

def get_effsnr(snr, chisq, chisq_dof):
	fac = 250.0
	rchisq = chisq/(2*chisq_dof-2)
	effsnr = snr/((1 + snr**2/fac)*rchisq)**(1./4)
	return effsnr

def get_snr(snr, chisq, chisq_dof):
	return snr

def get_snr_over_chi(snr, chisq, chisq_dof):
	return snr/chisq**(1./2)

def set_getsnr_function(connection, statistic):
	if statistic == "snr":
		connection.create_function('get_snr', 3, get_snr)
	elif statistic == "newsnr":
		connection.create_function('get_snr', 3, get_newsnr)
	elif statistic == "effsnr":
		connection.create_function('get_snr', 3, get_effsnr)
	elif statistic == "snroverchi":
		connection.create_function('get_snr', 3, get_snr_over_chi)


def end_time_w_ns(end_time, end_time_ns):
	time = end_time + 1e-9*end_time_ns
	return time



def get_sngl_snrs(
	connection,
	min_snr,
	sngls_width,
	ifo,
	tmplt,
	snr_stat = "newsnr",
	usertag = "FULL_DATA",
	datatype = None,
	sngls_bins = None):

	# create function for the desired snr statistic
	set_getsnr_function(connection, snr_stat)

	connection.create_function('end_time_w_ns', 2, end_time_w_ns)

	# split tmplt tuple into its component parameters
	mchirp = tmplt[0]
	eta = tmplt[1]

	# SQLite query to get a list of (snr, gps-time) tuples for inspiral triggers
	sqlquery = """
	SELECT DISTINCT
		get_snr(snr, chisq, chisq_dof) as snr_stat,
		end_time_w_ns(end_time, end_time_ns)
	FROM sngl_inspiral
		JOIN process_params ON (
			process_params.process_id == sngl_inspiral.process_id)
	"""
	# if a datatype is given, get only the inspiral triggers from coincs of that type
	if datatype:
		sqlquery += """
		JOIN coinc_event_map, experiment_map, experiment_summary ON (
			coinc_event_map.event_id == sngl_inspiral.event_id
			AND experiment_map.coinc_event_id == coinc_event_map.coinc_event_id
			AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
		"""
	sqlquery += ''.join(["""
	WHERE
		sngl_inspiral.ifo = \"""", ifo, """\" 
		AND snr_stat >= """, str(min_snr), """
		AND sngl_inspiral.mchirp = """, str(mchirp), """ 
		AND sngl_inspiral.eta = """, str(eta), """
		AND process_params.value = \"""", usertag, """\" """])
	if datatype:
		sqlquery += ''.join(["""
		AND experiment_summary.datatype = \"""", datatype, """\" """])

	# execute query
	trig_list = connection.cursor().execute( sqlquery ).fetchall()

	# get dq-veto segments
	xmldoc = dbtables.get_xml(connection)
	veto_segments = compute_dur.get_veto_segments(xmldoc, False)
	veto_segments = veto_segments[ veto_segments.keys()[0] ]

	snrlist = []
	# apple vetoes to the list of trigger times
	for snr, trig_time in trig_list:
		trig_segment = segments.segment(trig_time, trig_time)
		if not veto_segments[ifo].intersects_segment( trig_segment ):
			snrlist.append( snr )

	if not sngls_bins:
		sngls_bins = numpy.arange(min_snr, max(snrlist) + sngls_width, sngls_width)
	
	# make the binned snr histogram
	sngls_hist, junk = numpy.histogram(snrlist, bins=sngls_bins)

	return sngls_hist, sngls_bins


def get_coinc_snrs(
	connection,
	sngls_threshold,
	tmplt,
	datatype = None,
	little_dog = True,
	combined_bins = None,
	snr_stat = "newsnr"):

	# create function for the desired snr statistic
	set_getsnr_function(connection, snr_stat)

	# split tmplt tuple into its component parameters
	mchirp = tmplt[0]
	eta = tmplt[1]

	sqlquery = ''.join(["""
	SELECT coinc_inspiral.snr
	FROM coinc_inspiral
		JOIN sngl_inspiral, coinc_event_map, experiment_map, experiment_summary ON (
			coinc_inspiral.coinc_event_id == coinc_event_map.coinc_event_id
			AND coinc_event_map.event_id == sngl_inspiral.event_id
			AND experiment_map.coinc_event_id == coinc_inspiral.coinc_event_id
			AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
	WHERE
		experiment_summary.datatype = \"""", datatype, """\"
		AND get_snr(sngl_inspiral.snr, chisq, chisq_dof) >= """, str(sngls_threshold), """
		AND sngl_inspiral.mchirp = """, str(mchirp), """ 
		AND sngl_inspiral.eta = """, str(eta)])

	if not little_dog:
		sqlquery += """
		AND sngl_inspiral.event_id NOT IN (
			SELECT event_id
			FROM coinc_event_map
				JOIN experiment_map, experiment_summary ON (
					coinc_event_map.coinc_event_id == experiment_map.coinc_event_id
					AND experiment_map.experiment_summ_id == experiment_summary.experiment_summ_id)
			WHERE experiment_summary.datatype == "all_data")
		"""

	# execute query
	snrlist = connection.cursor().execute( sqlquery ).fetchall()

	# make histogram of coinc snr
	binned_snr, junk = numpy.histogram(snrlist, bins=combined_bins)

	return binned_snr


def combined_snr_hist(counts, snrs, combined_bins):
	minval = min(combined_bins)
	binWidth = combined_bins[1] - combined_bins[0]

	combined_counts = zeros(len(combined_bins)-1)
	ifos = counts.keys()
	if set(ifos) != set(snrs.keys()):
		raise ValueError, "The histogram and snr-bin dictionary have different sets of keys (ifos)"

	for n0, snr0 in zip(counts[ifos[0]], snrs[ifos[0]]):
		for n1, snr1 in zip(counts[ifos[1]], snrs[ifos[1]]):
			combined_snr = hypot(snr0, snr1)
			index =  int( math.floor((combined_snr - minval)/binWidth) )
			combined_counts[index] += n0 * n1

	return combined_counts


def all_possible_coincs(
	sngls_connection,
	sngls_width,
	min_snr,
	coinc_connection,
	coinc_width,
	ifos,
	tmplt,
	little_dog = True):

	# current code can only handle FAR estimation for doubles
	if len(ifos) > 2:
		raise ValueError, "Can only estimate FARs for doubles"

	# set up needed dictionaries
	mid_bins = {}
	sngls_hist = {}

	for ifo in ifos:
		# get all single ifo triggers for the given tmplt & ifo
		sngls_hist[ifo], sngls_bins = get_sngl_snrs(
			sngls_connection,
			min_snr,
			sngls_width,
			ifo,
			tmplt)
		mid_bins[ifo] = 0.5*( sngls_bins[1:] + sngls_bins[:-1] )

		# to avoid "little dogs", remove all sngls that make up a zerolag coinc
		if not little_dog:
			zerolag_sngls, junk = get_sngl_snrs(
				coinc_connection,
				min_snr,
				sngls_width,
				ifo,
				tmplt,
				datatype='all_data',
				sngls_bins = sngls_bins)
			for idx, N in enumerate(zerolag_sngls):
				sngls_hist[ifo][idx] -= N

	# define bins for combined-snr
	max_comb_snr = numpy.sum([max(snr_bins)**2.0 for snr_bins in mid_bins.values()])**(1./2)
	combined_bins = numpy.arange(2**(1./2)*min_snr, max_comb_snr + 2*coinc_width, coinc_width)

	coinc_hist = combined_snr_hist(sngls_hist, mid_bins, combined_bins)

	# remove only zerolag coincs from the coinc_hist
	if little_dog:
		zerolag_coincs = get_coinc_snrs(
			coinc_connection,
			min_snr,
			tmplt,
			datatype = 'all_data',
			combined_bins = combined_bins,
			little_dog = True)
		for idx, N in enumerate( zerolag_coincs ):
			coinc_hist[idx] -= N

	return coinc_hist, combined_bins

def get_coinc_time(connection, type):
	sqlquery = """
	SELECT duration
	FROM experiment_summary
	WHERE datatype = ?
	"""
	livetime = numpy.sum( connection.cursor().execute( sqlquery, (type,) ).fetchall() )

	return livetime

def get_singles_times( connection, verbose = False ):
	# get single-ifo filtered segments
	ifo_segments = compute_dur.get_single_ifo_segments(
		connection, program_name = "inspiral", usertag = "FULL_DATA")

	# get all veto segments
	xmldoc = dbtables.get_xml(connection)
	veto_segments = compute_dur.get_veto_segments(xmldoc, verbose)

	sngls_durations = {}

	for veto_def_name, veto_seg_dict in veto_segments.items():
		post_vetoes_ifosegs = ifo_segments - veto_seg_dict
		for ifo, segs_list in post_vetoes_ifosegs.items():
			sngls_durations[ifo] = float( abs(segs_list) )

	return sngls_durations


