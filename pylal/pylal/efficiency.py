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


def chirp_distance(eff_dist_h, eff_dist_l, mchirp):
    mchirp_dns = (1.4+1.4)*pow(0.25,3.0/5.0)
    dec_chirp_dist = min(eff_dist_h, eff_dist_l)*pow(mchirp_dns/mchirp, 5.0/6.0)

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

def successful_injections(connection, tag, verbose = False):
    """
    My attempt to get a list of the simulations that actually made
    it into some level of coincident time
    """

    xmldoc = dbtables.get_xml(connection)
    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    # Get the veto segments as dictionaries, keyed by veto category
    veto_segments = compute_dur.get_veto_segments(xmldoc, verbose)

    # sqlite query to get a list of the simulations for a given tag
    sqlquery = ''.join(["""
    SELECT
        simulation_id,
        end_time_with_ns(geocent_end_time, geocent_end_time_ns),
        distance
    FROM sim_inspiral
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
    WHERE process_params.value = \"""", tag, """\" """])

    if tag == 'ALL_INJ':
        sqlquery = '\n'.join( sqlquery.split('\n')[0:6] )
        # redefine tag for getting single-ifo segments
        tag = 'FULL_DATA'

    all_inj = connection.execute( sqlquery ).fetchall()

    # Get segments that define which time was filtered
    livetime_program = "inspiral"
    ifo_segments = compute_dur.get_single_ifo_segments(
        connection,
        program_name = livetime_program,
        usertag = tag)

    zero_lag_dict = {}
    for ifo in ifo_segments:
        zero_lag_dict[ifo] = 0.0

    successful_inj = {}
    post_veto_ifosegs = segments.segmentlistdict()
    # Apply vetoes to single-ifo filter segments
    for category, veto_seg_dict in veto_segments.items():
        successful_inj[category] = {}

        # determine coincident segments for that veto category 
        coinc_segs = compute_dur.get_coinc_segments(
            ifo_segments - veto_seg_dict,
            zero_lag_dict)

        for on_ifos in coinc_segs:
            successful_inj[category][on_ifos] = []
            for idx, injection in enumerate(all_inj):
                inj_segment = segments.segment(injection[1], injection[1])
                if coinc_segs[on_ifos].intersects_segment( inj_segment ):
                    successful_inj[category][on_ifos].append( injection )

    return successful_inj

def found_injections(connection, tag, verbose):
    connection.create_function('end_time_with_ns', 2, end_time_with_ns)
    xmldoc = dbtables.get_xml(connection)
    veto_def_names = set( table.get_table(xmldoc, lsctables.SegmentDefTable.tableName).getColumnByName('name') )

    ifos = [ifo[0] for ifo in connection.execute('SELECT DISTINCT ifos FROM process WHERE program = "inspiral"')]

    found_inj = {}
    inj_fars = {}
    for category in veto_def_names:
        found_inj[category] = {}
        inj_fars[category] = {}
        for on_ifos in compute_dur.get_allifo_combos(ifos, 2)[0]:
            sql_params_dict = {}
            sqlquery = """
            SELECT DISTINCT
                sim_inspiral.simulation_id,
                end_time_with_ns(geocent_end_time, geocent_end_time_ns),
                distance,
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
                sql_params_dict["ifos"] = on_ifos
                sql_params_dict["tag"] = tag

            injections = connection.execute(sqlquery, sql_params_dict).fetchall()
            injections.sort( key=itemgetter(3), reverse=True)

            found_inj[category][on_ifos] = [inj[0:3] for inj in injections]
            inj_fars[category][on_ifos] = [inj[3] for inj in injections]

    return found_inj, inj_fars


def get_four_volume(
    successful_inj,
    found_inj,
    found_fars,
    r,
    connection,
    on_ifos,
    veto_cat):

    # catching any edge cases were the injection end_time is nearly on a second boundary
    successful_inj = set(successful_inj) & set(found_inj)
    # histogram injections that went into Cat-N time
    successful_dist = [inj[2] for inj in successful_inj]
    N_success, junk = numpy.histogram(successful_dist, bins = r)

    significant_dist = [inj[2] for inj in found_inj]
    # Determine the minimum false_alarm_rate one can estimate from slides
    slide_time = get_livetime(connection, veto_cat, on_ifos, 'slide')
    minFAR = (3600.0*24.0*365.25) / slide_time
    far_list = minFAR * numpy.arange(1000,-1,-1)
    # Calculate the foreground search time in years
    zerolag_time = get_livetime(connection, veto_cat, on_ifos, 'all_data')
    T = zerolag_time / (3600.0*24.0*365.25)
    eff = {}
    VT = {}
    for threshold in far_list:
        for idx, far in enumerate(found_fars):
            if far <= threshold:
                new_start = idx
                break
        # Histogram found injections with FAR < threshold
        N_significant, junk = numpy.histogram(significant_dist[new_start:], bins = r)
        eff[threshold] = numpy.float_(N_significant)/N_success
        VT[threshold] = 0.0
        for idx, e in enumerate(eff[threshold]):
            if not math.isnan(e) and e > 0.0:
                VT[threshold] += 4.0/3.0*numpy.pi*e*T*(pow(r[idx+1]/2.0,3) - pow(r[idx]/2.0,3))

    return eff, VT

