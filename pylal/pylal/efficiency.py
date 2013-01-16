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
#                                 Preamble
#
# =============================================================================
#

"""
Collection of functions to compute the efficiency and effective 4-volume
"""

import sqlite3
import math
from operator import itemgetter
import numpy
from scipy import special

from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import dbtables

from pylal import antenna
from pylal import ligolw_sqlutils as sqlutils
from pylal import ligolw_compute_durations as compute_dur


#
# =============================================================================
#
#                                 
#
# =============================================================================
#


def chirp_dist(distance, mchirp):
        mchirp_DNS = (1.4+1.4) * (1./4)**(3.0/5.0)

        return distance * (mchirp_DNS/mchirp)**(5.0/6.0)

def decisive_dist(
    h_dist, l_dist, v_dist, 
    mchirp, weight_dist, ifos):
    
    dist_list = []
    if 'H1' in ifos or 'H2' in ifos:
        dist_list.append(h_dist)
    if 'L1' in ifos:
        dist_list.append(l_dist)
    if 'V1' in ifos:
        dist_list.append(v_dist) 

    if weight_dist:
        return chirp_dist(sorted(dist_list)[1], mchirp) 
    else:
        return sorted(dist_list)[1]

def end_time_with_ns(end_time, end_time_ns):
    time = end_time + 1e-9*end_time_ns
    return time

def get_livetime(connection, veto_cat, on_ifos, datatype):
    sqlquery = """
    SELECT duration
    FROM experiment_summary
        JOIN experiment ON (
            experiment_summary.experiment_id == experiment.experiment_id)
    WHERE
        datatype = :0
        AND veto_def_name = :1
        AND instruments = :2 """

    # total livetime in seconds 
    total_dur = numpy.sum(connection.execute(sqlquery, (datatype, veto_cat, on_ifos)).fetchall() )

    return total_dur

#
# =============================================================================
#
#                         Injections Functions
#
# =============================================================================
#

def inj_dist_range(d_min, d_max, dist_scale = "linear", step = 4.0):

    if dist_scale == "linear":
        dist_bin_edges = numpy.arange(d_min-step, d_max+step, step)
    elif dist_scale == "log":
        log_limits = numpy.log10([d_min, d_max])/numpy.log10(step)
        dist_bin_edges = numpy.power(
            step,
            numpy.arange(log_limits[0]-1, log_limits[1]+1)
        )

    return dist_bin_edges


def successful_injections(
    connection,
    tag,
    on_ifos,
    veto_cat,
    dist_type = "distance",
    weight_dist = False,
    verbose = False):

    """
    My attempt to get a list of the simulations that actually made
    it into some level of coincident time
    """

    xmldoc = dbtables.get_xml(connection)
    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    # Get the veto segments as dictionaries, keyed by veto category
    veto_segments = compute_dur.get_veto_segments(xmldoc, verbose)

    # ------------------------ Get List of Injections ------------------------ #
    sql_params_dict = {}
    sqlquery = """
        SELECT DISTINCT
            simulation_id,
            end_time_with_ns(geocent_end_time, geocent_end_time_ns),"""
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 2, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp)
        FROM sim_inspiral """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['ifos'] = on_ifos
        sql_params_dict['weight_dist'] = weight_dist
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :weight_dist, :ifos)
        FROM sim_inspiral """

    if tag != 'ALL_INJ':
        # if a specific injection set is wanted
        sqlquery += """
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
        WHERE process_params.value = :usertag) """
        sql_params_dict["usertag"] = tag
    else:
        # for all injections
        tag = 'FULL_DATA'

    # Get segments that define which time was filtered
    ifo_segments = compute_dur.get_single_ifo_segments(
        connection,
        program_name = "inspiral",
        usertag = tag)

    zero_lag_dict = dict([(ifo, 0.0) for ifo in ifo_segments])

    successful_inj = []
    # determine coincident segments for that veto category 
    coinc_segs = compute_dur.get_coinc_segments(
        ifo_segments - veto_segments[veto_cat],
        zero_lag_dict)

    # Apply vetoes to single-ifo filter segments
    for injection in connection.execute(sqlquery, sql_params_dict):
        inj_segment = segments.segment(injection[1], injection[1])
        if coinc_segs[on_ifos].intersects_segment( inj_segment ):
            successful_inj.append( injection )

    return successful_inj


def found_injections(
    connection,
    tag,
    on_ifos,
    dist_type = "distance",
    weight_dist = False,
    verbose = False):

    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    sql_params_dict = {}
    sqlquery = """
    SELECT DISTINCT
        sim_inspiral.simulation_id,
        end_time_with_ns(geocent_end_time, geocent_end_time_ns), """
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 2, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp), """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['weight_dist'] = weight_dist
        sql_params_dict['ifos'] = on_ifos
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :weight_dist, :ifos), """

    sqlquery += """
        false_alarm_rate,
        coinc_inspiral.snr
    FROM
        coinc_event_map AS coincs
        JOIN coinc_event_map AS sims, coinc_inspiral, coinc_event, sim_inspiral ON (
            coincs.coinc_event_id == sims.coinc_event_id
            AND coinc_event.coinc_event_id == coincs.event_id
            AND coinc_inspiral.coinc_event_id == coincs.event_id
            AND sim_inspiral.simulation_id == sims.event_id)
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
    WHERE
        coincs.table_name = "coinc_event"
        AND sims.table_name = "sim_inspiral" """

    if tag != 'ALL_INJ':
        sqlquery += """
        AND coinc_event.instruments = :ifos
        AND process_params.value = :usertag
        """
        sql_params_dict["ifos"] = on_ifos
        sql_params_dict["tag"] = tag

    injections = connection.execute(sqlquery, sql_params_dict).fetchall()
    injections.sort( key=itemgetter(3), reverse=True)

    found_inj = [inj[0:3] for inj in injections]
    inj_fars = [inj[3] for inj in injections]
    inj_snrs = [inj[4] for inj in injections]

    return found_inj, inj_fars, inj_snrs


def binomial_confidence(K, N, eff, confidence):
    # Calculate pdf peak and credible interval for p(eff|k,n)
    # posterior generated with binomial p(k|eff,n) and a uniform p(eff)

    values = {
        'peak': numpy.float_(K)/N,
        'low': numpy.zeros(len(K)),
        'high': numpy.zeros(len(K))
        }
    min_range = 1
    for I in range(0,len(K)):
        for a in eff:
            B_a = special.betainc(K[I]+1, N[I]-K[I]+1, a)
            if B_a > 1 - confidence:
                break
            b = special.betaincinv(K[I]+1, N[I]-K[I]+1, confidence + B_a)
            if b - a < min_range:
                min_range = b - a

        values['low'][I] = a
        values['high'][I] = b

    return values

def detection_efficiency(
    successful_inj,
    found_inj,
    found_fars,
    far_list,
    r,
    confidence):

    # catching any edge cases were the injection end_time is nearly on a second boundary
    successful_inj = set(successful_inj) | set(found_inj)
    # histogram of successful injections into coincidence time post vetoes
    successful_dist = [inj[2] for inj in successful_inj]
    N, _ = numpy.histogram(successful_dist, bins = r)

    significant_dist = [inj[2] for inj in found_inj]

    eff = {}
    eff_min = numpy.linspace(0, 1, 1e3+1)
    for threshold in far_list:
        for idx, far in enumerate(found_fars):
            if far <= threshold:
                new_start = idx
                break
        # Histogram found injections with FAR < threshold
        K, _ = numpy.histogram(significant_dist[new_start:], bins = r)

        # binomial confidence returns an array where each element is (peak, eff_low, eff_high)
        eff[threshold] = binomial_confidence(K, N, eff_min, confidence)

    return eff


def rescale_dist(on_ifos, distbins, old_distbins, dist_type, weight_dist):
    N_signals = int(1e6)
    trigTime = 0.0

    # if decisive distance is desired, get the antenna responses for each signal
    if dist_type == 'decisive_distance':
        # sky position (right ascension & declination)
        ra = 360 * numpy.random.rand(N_signals)
        dec = 180 * numpy.random.rand(N_signals) - 90
        # additional angles
        inclination = 180 * numpy.random.rand(N_signals)
        polarization = 360 * numpy.random.rand(N_signals)
   
        f_q = {} 
        for ifo in on_ifos:
            f_q[ifo] = numpy.zeros(N_signals)
            for index in range(N_signals):
                _, _, _, f_q[ifo][index] = antenna.response(
                   trigTime,
                   ra[index], dec[index],
                   inclination[index], polarization[index],
                   'degree', ifo )
    
    prob_d_d = {}
    for j in range(len(distbins)-1):
        # for this physical distance range, create signals that are uniform in volume
        volume = 4*numpy.pi/3 * numpy.random.uniform(
            low = distbins[j]**3.0,
            high = distbins[j+1]**3.0,
            size = N_signals)
        dist = numpy.power(volume*(3/(4*numpy.pi)), 1./3)

        # create decisive distance (if desired)
        if dist_type == 'decisive_distance':
            dist_eff = {}
            for ifo in on_ifos:
                dist_eff[ifo] = dist / f_q[ifo]
            dist_dec = numpy.sort(dist_eff.values(), 0)[1]

        # weight distance measure by chirp mass (if desired)
        if weight_dist:
            # Component masses are Gaussian distributed around the Chandrasekar mass
            mass1, mass2 = 0.13 * numpy.random.randn(2, N_signals) + 1.40
            mchirp = numpy.power(mass1+mass2, -1./5) * numpy.power(mass1*mass2, 3./5)
            if dist_type == 'decisive_distance':
                dist_chirp = chirp_dist(dist_dec, mchirp)
            if dist_type == 'distance':
                dist_chirp = chirp_dist(dist, mchirp)
            N_d, _ = numpy.histogram(dist_chirp, bins=old_distbins)
        else:
            N_d, _ = numpy.histogram(dist_dec, bins=old_distbins)
    
        prob_d_d[distbins[j+1]] = numpy.float_(N_d)/N_signals

    return prob_d_d

def eff_vs_dist(measured_eff, prob_d_d):
    eff_dist = {}
    for far, values in measured_eff.items():
        eff_dist[far] = {
            'peak': numpy.zeros(len(prob_d_d)),
            'var_plus': numpy.zeros(len(prob_d_d)),
            'var_minus': numpy.zeros(len(prob_d_d))
            }
        for idx, p in enumerate(prob_d_d.values()):
            eff_dist[far]['peak'][idx] = numpy.sum(values['peak'] * p)
            eff_dist[far]['var_plus'][idx] = numpy.sum( (values['high'] - values['peak'])**2 * p**2 )
            eff_dist[far]['var_minus'][idx] = numpy.sum( (values['low'] - values['peak'])**2 * p**2 )


def volume_efficiency(measured_eff, r, old_distbins, prob_d_d):

    vol_shell = 4./3 * numpy.pi * (r[1:]**3 - r[:-1]**3)
    
    prob_d = [0]*len(old_distbins[1:])
    for i in range( len(old_distbins[1:]) ):
        B = numpy.zeros( len(prob_d_d) )
        for j, p in enumerate(prob_d_d.values()):
            B[j] = p[i]*vol_shell[j]
        prob_d[i] = numpy.cumsum(B)/numpy.cumsum(vol_shell)
    
    vol_eff = {}
    for far, values in measured_eff.items():
        vol_eff[far] = {}
        vol_eff[far]['peak'] = numpy.sum([ values['peak'][i] * prob_d[i] 
            for i in range(len(prob_d))], 0 )
        vol_eff[far]['var_plus'] = numpy.sum([ (values['high'][i]-values['peak'])**2 * prob_d[i]**2
            for i in range(len(prob_d))], 0 )
        vol_eff[far]['var_minus'] = numpy.sum([ (values['low'][i]-values['peak'])**2 * prob_d[i]**2
            for i in range(len(prob_d))], 0 )


