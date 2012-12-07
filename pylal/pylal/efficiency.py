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


def chirp_distance(distance, mchirp, source_ref):
    if source_ref == "NSNS":
        mchirp_ref = (1.4+1.4) * (1./4)**(3.0/5.0)
    elif source_ref == "NSBH":
        mchirp_ref = (1.4+10.0) * (14.0/11.4**2)**(3.0/5.0)
    elif source_ref == "BHBH":
        mchirp_ref = (10.0+10.0) * (1./4)**(3.0/5.0)
    chirp_dist = distance * (mchirp_ref/mchirp)**(5.0/6.0)

    return chirp_dist


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

def inj_dist_range(
    connection,
    tag,
    dist_type = "distance",
    dist_scale = "linear",
    source_ref = None,
    num_bins = 10):

    sql_params_dict = {}
    # SQL query to get a list of the distances for desired injections
    if dist_type == "chirp_dist":
        connection.create_function('chirp_distance', 3, chirp_distance)
        sqlquery += """
        SELECT DISTINCT
            chirp_distance(distance, mchirp, :source_ref)
        FROM sim_inspiral """
        sql_params_dict['source_ref'] = source_ref
    elif dist_type == "distance":
        sqlquery += """
        SELECT DISTINCT
            distance
        FROM sim_inspiral """

    if tag != 'ALL_INJ':
        # if a specific injection set is wanted
        sqlquery += """
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
        WHERE process_params.value = :usertag) """
        sql_params_dict["usertag"] = tag

    distances = numpy.array([ dist[0]
        for dist in connection.execute(sqlquery, sql_params_dict) ])
    min_dist = numpy.min(distances)
    max_dist = numpy.max(distances)

    if dist_scale == "linear":
       d = numpy.linspace(start=min_dist, stop=max_dist, num=num_bins)
    elif dist_scale == "log":
       d = numpy.logspace(
           start=numpy.log10(min_dist),
           stop=numpy.log10(max_dist),
           num=num_bins)

    return d


def successful_injections(
    connection,
    tag,
    on_ifos,
    veto_cat,
    distance_type = "distance",
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
    sql_params_dict = {}
    sqlquery = """
        SELECT
            simulation_id,
            end_time_with_ns(geocent_end_time, geocent_end_time_ns),"""
    # add the desired distance measure to the SQL query
    if distance_type == "chirp_dist":
        connection.create_function('chirp_distance', 3, chirp_distance)
        sql_params_dict['source_ref'] = source_ref
        sqlquery += """
            chirp_distance(distance, mchirp, :source_ref)
        FROM sim_inspiral """
    elif distance_type == "distance":
        sqlquery += """
            distance
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
        ifo_segments - veto_seg_dict,
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
    veto_cat,
    verbose = False):

    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    ifos_str = ','.join(on_ifos)

    sql_params_dict = {}
    sqlquery = """
    SELECT DISTINCT
        sim_inspiral.simulation_id,
        end_time_with_ns(geocent_end_time, geocent_end_time_ns),
        false_alarm_rate,"""
    # add the desired distance measure to the SQL query
    if distance_type == "chirp_dist":
        connection.create_function('chirp_distance', 3, chirp_distance)
        sqlquery += "chirp_distance(distance, mchirp, :source_ref)"
        sql_params_dict['source_ref'] = source_ref
    elif distance_type == "distance":
        sqlquery += "distance"

    sqlquery += """
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
    far_list):

    # catching any edge cases were the injection end_time is nearly on a second boundary
    successful_inj = set(successful_inj) & set(found_inj)
    # histogram injections that went into Cat-N time
    successful_dist = [inj[2] for inj in successful_inj]
    N_success, junk = numpy.histogram(successful_dist, bins = r)

    significant_dist = [inj[2] for inj in found_inj]

    eff = {}
    for threshold in far_list:
        for idx, far in enumerate(found_fars):
            if far <= threshold:
                new_start = idx
                break
        # Histogram found injections with FAR < threshold
        N_significant, junk = numpy.histogram(significant_dist[new_start:], bins = r)
        eff[threshold] = numpy.float_(N_significant)/N_success

    return eff

def get_four_volume(eff, r, T_z):
    # calculate 3 volume in each shell
    V = 4./3 * numpy.pi * (r[:-1]**3. - r[1:]**3.)

    VT = {}
    for threshold in eff:
        VT[threshold] = eff[threshold] * V * T_z
        for index, vt in VT[threshold]:
            if numpy.isnan(vt):
                VT[threshold][index] = 0.0

    return VT

