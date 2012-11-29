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

from glue.ligolw import dbtables
from glue import segments

from pylal import ligolw_sqlutils
from pylal import ligolw_compute_durations as compute_dur



#
# =============================================================================
#
#				   Library API 
#
# =============================================================================
#

def get_newsnr(snr, chisq, chisq_dof):
	fac = 6.0
	rchisq = chisq/(2*chisq_dof-2)
	newsnr = snr/(0.5*(1+rchisq**(fac/2.0)))**(1.0/fac)
	if rchisq < 1.0:
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
	"""
	Set the get_snr SQL function. Available options are:
	-- "rawsnr", "newsnr", "effsnr", "snroverchi"
	"""
	if statistic == "rawsnr":
		connection.create_function('get_snr', 3, get_snr)
	elif statistic == "newsnr":
		connection.create_function('get_snr', 3, get_newsnr)
	elif statistic == "effsnr":
		connection.create_function('get_snr', 3, get_effsnr)
	elif statistic == "snroverchi":
		connection.create_function('get_snr', 3, get_snr_over_chi)
	else:
		raise ValueError, "%s is not a valid snr statistic" % statistic

def end_time_w_ns(end_time, end_time_ns):
	time = end_time + 1e-9*end_time_ns
	return time

def quadrature_sum(tuple):
	snrsq_array = numpy.array(tuple)**2.
	comb_snr = numpy.sum( snrsq_array )**(1/2.)
	return comb_snr

def get_sngl_snrs(
	connection,
	ifo,
	mchirp,
	eta,
	min_snr,
	sngls_width = None,
	usertag = "FULL_DATA",
	datatype = None,
	sngls_bins = None):
	"""
	Creates a histogram of sngl_inspiral triggers and returns a list of counts
	and the associated snr bins.

	@connection: connection to a SQLite database with lsctables
	@ifo: the instrument one desires triggers from
	@mchirp: the chirp mass from the desired template
	@eta: the symmetric mass ratio from the desired template
	@min_snr: a lower threshold on the value of the snr_stat
	@sngls_width: the bin width for the histogram
	@usertag: the usertag for the triggers. The default is "FULL_DATA".
	@datatype: the datatype (all_data, slide, ...) if single-ifo triggers from
		coincident events is desired. The default is to collect all triggers.
	@sngls_bins: a list of bin edges for the snr-histogram
	"""

	# create function for the desired snr statistic
	sqlquery = """
	SELECT value
	FROM process_params
	WHERE
		program == 'ligolw_thinca'
		AND param == '--weighted-snr'
	"""
	snr_stat = connection.cursor().execute( sqlquery ).fetchone()[0]
	set_getsnr_function(connection, snr_stat)

	connection.create_function('end_time_w_ns', 2, end_time_w_ns)

	# set SQL statement parameters
	sql_params_dict = {
		"mchirp": mchirp, "eta": eta,
		"min_snr": min_snr, "ifo": ifo,
		"type": datatype}

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
	sqlquery += """
	WHERE
		sngl_inspiral.ifo == :ifo 
		AND snr_stat >= :min_snr
		AND sngl_inspiral.mchirp == :mchirp
		AND sngl_inspiral.eta == :eta
		AND process_params.value == :usertag
	"""
	if datatype:
		sqlquery += """
		AND experiment_summary.datatype == :type
		"""
		sql_params_dict["type"] = datatype

	# execute query
	trig_list = connection.cursor().execute( sqlquery, sql_params_dict ).fetchall()

	# get dq-veto segments
	xmldoc = dbtables.get_xml(connection)
	veto_segments = compute_dur.get_veto_segments(xmldoc, False)
	veto_segments = veto_segments[ veto_segments.keys()[0] ]

	snr_array = array([])
	# apply vetoes to the list of trigger times
	for snr, trig_time in trig_list:
		trig_segment = segments.segment(trig_time, trig_time)
		if not veto_segments[ifo].intersects_segment( trig_segment ):
			numpy.append( snr_array, snr )

	if not sngls_bins:
		sngls_bins = numpy.arange(min_snr, max(snrlist) + sngls_width, sngls_width)
	
	# make the binned snr histogram
	sngls_hist, junk = numpy.histogram(snr_array, bins=sngls_bins)

	return sngls_hist, sngls_bins

def all_sngls_hist(
	connection,
	min_snr,
	sngls_width,
	mchirp,
	eta,
	all_ifos
	little_dog = True):
	"""
	Creates a pair of dictionaries containing single-ifo snr histograms and
	their associated snr bins. The keys are the instruments filtered in the
	analysis.

	@connection: connection to a SQLite database with lsctables
	@min_snr: a lower threshold on the value of the snr_stat
	@sngls_width: the bin width for the histogram
	@mchirp: the chirp mass from the desired template
	@eta: the symmetric mass ratio from the desired template
	@all_ifos: a list containing the instruments used in the analysis 
	@little_dog: if argument is False, all coincs with a single-ifo trigger that 
		also constitutes part of a zerolag coinc are NOT included. The 
		default value is True.
	"""

	sngl_ifo_hist = {}
	sngl_ifo_midbins = {}
	# loop over the analyzed ifos
	for ifo in all_ifos:
		sngl_ifo_hist[ifo], bins = sngl_tmpltget_sngl_snrs(
			connection,
			ifo,
			mchirp, eta
			min_snr,
			sngls_width = sngls_width)
		# define the midpoint of each snr bin
		sngl_ifo_midbins[ifo] = 0.5*( bins[1:] + bins[:-1] )

		if not little_dog:
			# if one does not want "little dogs"
			zerolag_hist, junk = get_sngl_snrs(
				connection,
				ifo,
				mchirp, eta,
				min_snr
				datatype = "all_data",
				sngls_bins = bins)
			# remove zerolag-coinc constituents from singles histogram
			for idx, N in enumerate(zerolag_hist):
				sngl_ifo_hist[ifo][idx] -= N
	
	return sngl_ifo_hist, sngl_ifo_midbins

def get_coinc_snrs(
	connection,
	sngls_threshold,
	ifos,
	mchirp,
	eta,
	datatype = None,
	little_dog = True,
	combined_bins = None):

	"""
	Creates a histogram of coinc_inspiral triggers and returns a list of counts
	in each of the snr bins.

	@connection: connection to a SQLite database with lsctables
	@sngls_threshold: a lower threshold on the value of the snr_stat
	@ifos: a tuple containing the ifos a coinc must come from 
	@mchirp: the chirp mass from the desired template
	@eta: the symmetric mass ratio from the desired template
	@datatype: the datatype (all_data, slide, ...) if single-ifo triggers from
		coincident events is desired. The default is to collect all triggers.
	@little_dog: if argument is False, all coincs with a single-ifo trigger that 
		also constitutes part of a zerolag coinc are NOT included. The 
		default value is True.
	@combined_bins: a list of bin edges for the snr-histogram
	"""
	# create function for the desired snr statistic
	sqlquery = """
	SELECT value
	FROM process_params
	WHERE
		program == 'ligolw_thinca'
		AND param == '--weighted-snr'
	"""
	snr_stat = connection.cursor().execute( sqlquery ).fetchone()[0]
	set_getsnr_function(connection, snr_stat)

	# set SQL statement parameters
	query_params_dict = {
		"ifo1": ifos[0], "ifo2": ifos[1],
		"type": datatype,
		"min_snr": sngls_threshold,
		"mchirp": mchirp, "eta": eta }

	# get a list of the combined snrs from the coinc_inspiral table
	sqlquery = """
	SELECT
		get_snr(si_ifo1.snr, si_ifo1.chisq, si_ifo1.chisq_dof),
		get_snr(si_ifo2.snr, si_ifo2.chisq, si_ifo2.chisq_dof)"""
	if len(ifos) > 2:
		sqlquery += """,
		get_snr(si_ifo3.snr, si_ifo3.chisq, si_ifo3.chisq_dof),
		"""
		query_params_dict["ifo3"] = ifos[2]

	sqlquery += """
	FROM coinc_inspiral
		JOIN sngl_inspiral AS si_ifo1, coinc_event_map AS cem_ifo1 ON (
			coinc_inspiral.coinc_event_id == cem_ifo1.coinc_event_id
			AND cem_ifo1.event_id == si_ifo1.event_id
			AND si_ifo1.ifo == :ifo1)
		JOIN sngl_inspiral AS si_ifo2, coinc_event_map AS cem_ifo2 ON (
			coinc_inspiral.coinc_event_id == cem_ifo2.coinc_event_id
			AND cem_ifo2.event_id == si_ifo2.event_id
			AND si_ifo2.ifo == :ifo2)"""
	if len(ifos) > 2:
		sqlquery += """
		JOIN sngl_inspiral AS si_ifo3, coinc_event_map AS cem_ifo3 ON (
			coinc_inspiral.coinc_event_id == cem_ifo3.coinc_event_id
			AND cem_ifo3.event_id == si_ifo3.event_id
			AND si_ifo3.ifo == :ifo3)"""

	sqlquery += """
		JOIN experiment_map, experiment_summary ON (
			experiment_map.coinc_event_id == coinc_inspiral.coinc_event_id
			AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
	WHERE
		experiment_summary.datatype == :type
		AND si_ifo1.mchirp == :mchirp
		AND si_ifo1.eta == :eta
		AND get_snr(si_ifo1.snr, si_ifo1.chisq, si_ifo1.chisq_dof) >= :min_snr
		AND get_snr(si_ifo2.snr, si_ifo2.chisq, si_ifo2.chisq_dof) >= :min_snr"""
	if len(ifos) > 2:
		sqlquery += """
		AND get_snr(si_ifo3.snr, si_ifo3.chisq, si_ifo3.chisq_dof) >= :min_snr"""

	if not little_dog:
		zerolag_eids_script = """
		CREATE TEMP TABLE zerolag_eids AS
			SELECT event_id
			FROM coinc_event_map
				JOIN experiment_map, experiment_summary ON (
					coinc_event_map.coinc_event_id == experiment_map.coinc_event_id
					AND experiment_map.experiment_summ_id == experiment_summary.experiment_summ_id)
			WHERE experiment_summary.datatype == "all_data"
		"""
		connection.cursor().execute( zerolag_eids_script )

		sqlquery += """
		AND si_ifo1.event_id NOT IN (SELECT event_id FROM zerolag_eids)
		AND si_ifo2.event_id NOT IN (SELECT event_id FROM zerolag_eids)"""
		if len(ifos) > 2:
			sqlquery += """
			AND si_ifo3.event_id NOT IN (SELECT event_id FROM zerolag_eids)
			"""

	# execute query
	snr_array = numpy.array([quadrature_sum(snrs) for snrs in connection.execute( sqlquery, query_params_dict )])

	if not little_dog:
		connection.cursor().execute('DROP TABLE zerolag_eids')

	# make histogram of coinc snr
	binned_snr, junk = numpy.histogram(snr_array, bins=combined_bins)

	return map(float, binned_snr)


def combined_snr_hist(counts, snrs, combined_bins):
	minval = min(combined_bins)
	binWidth = combined_bins[1] - combined_bins[0]

	combined_counts = numpy.zeros(len(combined_bins)-1)
	ifos = counts.keys()
	if set(ifos) != set(snrs.keys()):
		raise ValueError, "The histogram and snr-bin dictionary have different sets of keys (ifos)"

	for n0, snr0 in zip(counts[ifos[0]], snrs[ifos[0]]):
		for n1, snr1 in zip(counts[ifos[1]], snrs[ifos[1]]):
			combined_snr = quadrature_sum( (snr0, snr1) )
			index =  int( numpy.floor((combined_snr - minval)/binWidth) )
			combined_counts[index] += n0 * n1

	return combined_counts


def all_possible_coincs(
	sngl_ifo_hist,
	sngl_ifo_midbins,
	combined_bins,
	zerolag_coinc_hist,
	ifos
	little_dog = True):
	"""
	Creates a histogram of all possible coincident events and returns a list of counts
	in each of the snr bins. This is made using the single-ifo snr histograms.

	@sngl_ifo_hist: a dictionary containing an snr histogram for each IFO
	@sngl_ifo_midbins: a dictionary the snr mid-bin values for sngl_ifo_hist
	@combined_bins: a list of bin edges for the combined-snr histogram
	@zerolag_coinc_hist: the snr histogram for zerolag coincident events
	@ifos: a list of analyzed instruments to generate a coinc
	@little_dog: if argument is False, all coincs with a single-ifo trigger that 
		also constitutes part of a zerolag coinc are NOT included. The 
		default value is True.
	"""

	if len(ifos) > 3:
		raise ValueError, "Can only estimate FARs for doubles & triples"

	binWidth = combined_bins[1] - combined_bins[0]
	combined_counts = numpy.zeros(len(combined_bins)-1)

	# for computing doubles rate: N01_ij = n0_i*n1_j
	N01 = numpy.outer(sngl_ifo_hist[ifo[0]], sngl_ifo_hist[ifo[1]])
	len0, len1 = N_01.shape
	for idx0, snr0 in enumerate(sngl_ifo_midbins[ifos[0]]:
		for idx1, snr1 in enumerate(sngl_ifo_midbins[ifos[1]]:
			if len(ifos) == 2:
				combined_snr = quadrature_sum( (snr0, snr1) )
				index =  int( numpy.floor((combined_snr - min(combined_bins))/binWidth) )
				combined_counts[index] += N01[idx0,idx1]
			else:
				# for computing triples rate: N012_ijk = n0_i*n1_j*n2_k
				N012 = numpy.outer(N01, sngl_ifo_hist[ifo[2]])
				N012 = N_012.reshape(len0, len1, N012.size/(len0*len1))
				for idx2, snr2 in enumerate(sngl_ifo_midbins[ifos[2]]:
					combined_snr = quadrature_sum( (snr0, snr1, snr2) )
					index =  int( numpy.floor((combined_snr - min(combined_bins))/binWidth) )
					combined_counts[index] += N012[idx0,idx1,idx2]

	# if keeping "little-dogs", remove only zerolag coincs from the coinc_hist
	if little_dog:
		for idx, N in enumerate( zerolag_coinc_hist ):
			combined_counts[idx] -= N

	return combined_counts


def inclusive_coinc_time(connection, type, inc_ifo_list):
	sqlquery = """
	SELECT duration, instruments
	FROM experiment_summary
		JOIN experiment ON (
			experiment.experiment_id == experiment_summary.experiment_id)
	WHERE datatype = ?
	"""
	def is_ifosubset(inclusive_list, list):
		if len(set(inclusive_list) & set(list)) == len(inclusive_list):
			return True
		else:
			return False

	livetime = numpy.sum([T for T,ifos in connection.execute(sqlquery, (type,)) if is_ifosubset(inc_ifo_list,ifos)])

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


def get_coinc_window(connection, ifos):
	sqlquery = """
	SELECT DISTINCT si_ifo1.end_time, si_ifo2.end_time, si_ifo1.end_time_ns, si_ifo2.end_time_ns
	FROM coinc_inspiral
		JOIN sngl_inspiral AS si_ifo1, coinc_event_map AS cem_ifo1 ON (
			coinc_inspiral.coinc_event_id == cem_ifo1.coinc_event_id 
			AND cem_ifo1.event_id == si_ifo1.event_id
			AND si_ifo1.ifo == ?)
		JOIN sngl_inspiral AS si_ifo2, coinc_event_map AS cem_ifo2 ON (
			coinc_inspiral.coinc_event_id == cem_ifo2.coinc_event_id 
			AND cem_ifo2.event_id == si_ifo2.event_id
			AND si_ifo2.ifo == ?)
	"""
	coincs = connection.cursor().execute( sqlquery, tuple(ifos) ).fetchall()
	
	shift = connection.cursor().execute('SELECT MIN(ABS(offset)) FROM time_slide WHERE offset != 0.0 ').fetchone()[0]
	toa_diff = []
	for coinc in coincs:
		dT_sec = (coinc[0]-coinc[1]) + (coinc[2]-coinc[3])*1e-9
		num_shifts = round(dT_sec/shift)
		toa_diff.append( dT_sec - num_shifts*shift )

	return max(toa_diff) - min(toa_diff)

def compute_cumrate(hist, T_bkgd):
	cum_hist = numpy.cumsum( hist )[::-1]
	rate = {}
	# moments of rate pdf based on the assumption that bkgd coincs are a Poisson process
	rate["mode"] = cum_hist * secINyr/T_bkgd
	rate["mean"] = (cum_hist +`1) * secINyr/T_bkgd
	rate["stdev"] = (cum_hist +`1)**(1/2.) * secINyr/T_bkgd

	return rate
