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

from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import dbtables

from pylal import ligolw_sqlutils as sqlutils
from pylal import ligolw_compute_durations as compute_dur


#
# =============================================================================
#
#                                 
#
# =============================================================================
#


def chirp_dist(distance, mchirp, source_ref):
    if not source_ref:
        return distance
    else:
        if source_ref == "NSNS":
            mchirp_ref = (1.4+1.4) * (1./4)**(3.0/5.0)
        elif source_ref == "NSBH":
            mchirp_ref = (1.4+10.0) * (14.0/11.4**2)**(3.0/5.0)
        elif source_ref == "BHBH":
            mchirp_ref = (10.0+10.0) * (1./4)**(3.0/5.0)

        return distance * (mchirp_ref/mchirp)**(5.0/6.0)

def decisive_dist(
    h_dist, l_dist, v_dist, 
    mchirp, source_ref, ifos):
    
    dist_list = []
    if 'H1' in ifos or 'H2' in ifos:
        dist_list.append(h_dist)
    if 'L1' in ifos:
        dist_list.append(l_dist)
    if 'V1' in ifos:
        dist_list.append(v_dist) 

    return chirp_dist(sorted(dist_list)[-2], mchirp, source_ref) 
 

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

def inj_dist_range(dist_array, dist_scale = "linear", step = 4.0):

    d_min = numpy.min(dist_array)
    d_max = numpy.max(dist_array)

    if dist_scale == "linear":
        dist_bin_edges = numpy.arange(d_min, d_max+2*step, step)
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
    source_ref = None,
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
    sql_params_dict = {'source_ref': source_ref}
    sqlquery = """
        SELECT DISTINCT
            simulation_id,
            end_time_with_ns(geocent_end_time, geocent_end_time_ns),"""
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 3, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp, :source_ref)
        FROM sim_inspiral """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['ifos'] = on_ifos
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :source_ref, :ifos)
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
    source_ref = None,
    verbose = False):

    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    ifos_str = ','.join(on_ifos)

    sql_params_dict = {'source_ref': source_ref}
    sqlquery = """
    SELECT DISTINCT
        sim_inspiral.simulation_id,
        end_time_with_ns(geocent_end_time, geocent_end_time_ns), """
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 3, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp, :source_ref), """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['ifos'] = on_ifos
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :source_ref, :ifos), """

    sqlquery += """
        false_alarm_rate
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
        sql_params_dict["ifos"] = ifos_str
        sql_params_dict["tag"] = tag

    injections = connection.execute(sqlquery, sql_params_dict).fetchall()
    injections.sort( key=itemgetter(3), reverse=True)

    found_inj = [inj[0:3] for inj in injections]
    inj_fars = [inj[3] for inj in injections]

    return found_inj, inj_fars


def detection_efficiency(
    successful_inj,
    found_inj,
    found_fars,
    far_list,
    r):

    # catching any edge cases were the injection end_time is nearly on a second boundary
    successful_inj = set(successful_inj) | set(found_inj)
    # histogram injections that went into Cat-N time
    successful_dist = [inj[2] for inj in successful_inj]
    N_success, junk = numpy.histogram(successful_dist, bins = r)

    significant_dist = [inj[2] for inj in found_inj]

    eff = {}
    eff_stdev = {}
    for threshold in far_list:
        for idx, far in enumerate(found_fars):
            if far <= threshold:
                new_start = idx
                break
        # Histogram found injections with FAR < threshold
        N_significant, junk = numpy.histogram(significant_dist[new_start:], bins = r)
        eff[threshold] = numpy.float_(N_significant)/N_success
        eff_stdev[threshold] = numpy.sqrt(eff[threshold]*(1-eff[threshold]) / N_success)

    return eff, eff_stdev

def get_four_volume(eff, eff_stdev, r, T_fgd):
    # calculate 3 volume in each shell
    V = 4./3 * numpy.pi * (r[1:]**3. - r[:-1]**3.)

    VT = {}
    VT_var = {}
    # for a given FAR threshold
    for threshold, eff_array in eff.items():
        # make the cumulative VxT
        VT[threshold] = numpy.cumsum(numpy.nan_to_num(eff_array * V * T_fgd))
        # compute the variance for the cumulative VxT
        VT_var[threshold] = numpy.cumsum(
            numpy.power(numpy.nan_to_num(eff_stdev[threshold] * V * T_fgd), 2.0)
        )

    return VT, VT_var

