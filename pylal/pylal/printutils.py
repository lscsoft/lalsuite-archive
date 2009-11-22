#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

'''
A collection of utilities to assist in printing out information from an xmldoc.
'''

import sys
import time
import datetime

from glue.ligolw.utils import print_tables
from glue.ligolw import ligolw
from glue.ligolw import table

from pylal.xlal.date import XLALGPSToUTC
try:
    from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
except ImportError:
    # s6 code
    from pylal.xlal.date import LIGOTimeGPS


__author__ = "Collin Capano <cdcapano@physics.syr.edu>"
__date__ = "$Date$" 
__version__ = "$Revision$"



# =============================================================================
#
#                           Utilities
#
# =============================================================================

def generic_get_pyvalue(obj):
    if obj.value is None:
        return None
    return ligolwtypes.ToPyType[obj.type or "lstring"](obj.value)


def get_columns_to_print(xmldoc, tableName):
    """
    Retrieves canonical columns to print for the given tableName.
    Returns a columnList, row_span and rspan_break lists.
    """
    tableName = tableName.endswith(":table") and tableName or tableName+":table"
    summTable = table.get_table(xmldoc, tableName )
    # get rankname
    rankname = [col.getAttribute("Name").split(":")[-1]
        for col in summTable.getElementsByTagName(u'Column') if "rank" in col.getAttribute("Name")][0]

    if tableName == "loudest_events:table":
        durname = [col.getAttribute("Name").split(":")[-1]
            for col in summTable.getElementsByTagName(u'Column') if "duration" in col.getAttribute("Name")][0]
        columnList = [
            rankname,
            'combined_far',
            'fap',
            'fap_1yr',
            'snr',
            'end_time',
            'end_time_utc__Px_click_for_daily_ihope_xP_',
            'ifos__Px_click_for_elog_xP_',
            'mass',
            'mchirp',
            'mini_followup',
            durname]
        row_span_columns = rspan_break_columns = [durname]
    elif tableName == "selected_found_injections:table":
        durname = [col.getAttribute("Name").split(":")[-1]
            for col in summTable.getElementsByTagName(u'Column') if "duration" in col.getAttribute("Name")][0]
        columnList = [
            rankname,
            'injected_end_time',
            'injected_end_time_utc__Px_click_for_daily_ihope_xP_',
            'elogs',
            'mini_followup',
            'injected_eff_dist_h',
            'injected_eff_dist_l',
            'injected_eff_dist_v',
            'injected_mchirp',
            'injected_mass1',
            'injected_mass2',
            'recovered_match_rank',
            'recovered_ifos',
            'recovered_combined_far',
            'recovered_fap',
            'recovered_fap_1yr',
            'recovered_snr',
            'recovered_end_time',
            'recovered_mchirp',
            'recovered_mass']
        row_span_columns = rspan_break_columns = [
            rankname,
            'injected_end_time',
            'injected_end_time_utc__Px_click_for_daily_ihope_xP_',
            'elogs',
            'mini_followup',
            'injected_eff_dist_h',
            'injected_eff_dist_l',
            'injected_eff_dist_v',
            'injected_mchirp',
            'injected_mass1',
            'injected_mass2']
    elif tableName == "close_missed_injections:table":
        columnList = [
            'rank',
            'decisive_distance',
            'end_time',
            'end_time_utc__Px_click_for_daily_ihope_xP_',
            'elogs',
            'mchirp',
            'mass1',
            'mass2',
            'eff_dist_h',
            'eff_dist_l',
            'eff_dist_v',
            'mini_followup'
            ]
        row_span_columns = rspan_break_columns = []
    else:
        # unrecognized table, just return all the columns in the table
        columnList = [col.getAttribute("Name").split(":")[-1] for col in summTable.getElementsByTagName(u'Column')]
        row_span_columns = rspan_break_columns = []
        
    return columnList, row_span_columns, rspan_break_columns


#
#   Some helper functions for manipulating times
#

def get_dst_start_end(ifo, year):
    """
    Figures out what dates daylight savings time starts and ends at a given site on a given year.
    """
    # in the United Stats, prior to 2007, DST began on the first Sunday in April
    # and ended on the last Sunday in October (http://aa.usno.navy.mil/faq/docs/daylight_time.php)
    if ("H" in ifo  or "L" in ifo) and year < 2007:
        for ii in range(1,28):
            dst_start = datetime.datetime(year, 4, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday':
                break
        for ii in range(31,0,-1):
            dst_end = datetime.datetime(year, 10, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    # in the US, starting in 2007, DST begins on the second Sunday in March and ends on the first
    # Sunday in November
    elif ("H" in ifo  or "L" in ifo) and year >= 2007:
        nn = 1
        for ii in range(1,31):
            dst_start = datetime.datetime(year, 3, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday' and nn == 2:
                break
            elif dst_start.strftime('%A') == 'Sunday':
                nn += 1
        for ii in range(1,28):
            dst_end = datetime.datetime(year, 11, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    # in Europe, DST begins on the last Sunday of March and ends on the last Sunday of October
    # source: http://www.timeanddate.com/news/time/europe-dst-starts-march-29-2009.html
    elif ("V" in ifo or "G" in ifo):
        for ii in range(31,0,-1):
            dst_start = datetime.datetime(year, 3, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday':
                break
        for ii in range(31,0,-1):
            dst_end = datetime.datetime(year, 10, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    else:
        raise ValueError, "unrecognized ifo %s" % ifo
    
    return dst_start, dst_end
        

def get_sitelocaltime_from_gps(ifo, gpstime):
    # get the utc time in datetime.datetime format
    utctime = XLALGPSToUTC(LIGOTimeGPS(gpstime, 0))
    utctime = datetime.datetime(utctime[0],utctime[1],utctime[2],utctime[3],utctime[4],utctime[5],utctime[6])
    # figure out if daylight savings time was on or not
    dst_start, dst_end = get_dst_start_end(ifo, utctime.year)
    # figure out the appropriate time offset
    if "H" in ifo:
        toffset = datetime.timedelta(hours=-7)
    elif "L" in ifo:
        toffset = datetime.timedelta(hours=-5)
    elif ("V" in ifo or "G" in ifo):
        toffset = datetime.timedelta(hours=+2)
    # apply the dst time offset to see if daylight savings was on; if not, adjust the toffset
    if not (utctime + toffset >= dst_start and utctime + toffset < dst_end):
        toffset = toffset + datetime.timedelta(hours=-1)

    return utctime + toffset


def format_end_time_in_utc(gps_sec):
    return time.strftime("%a %d %b %Y %H:%M:%S", XLALGPSToUTC(LIGOTimeGPS(gps_sec, 0)))


def get_elog_page(ifo, gpstime):
    # set site_address
    if "H" in ifo:
        site_address = "http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?group=detector"
    elif "L" in ifo:
        site_address = "http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?group=detector"
    elif "V" in ifo:
        #FIXME: What's the site address and format for Virgo log book?
        site_address = "https://pub3.ego-gw.it/logbook/"
    # get local time at the site
    site_localtime = get_sitelocaltime_from_gps(ifo, gpstime)
    # set the address
    if "H" in ifo or "L" in ifo:
        site_address = "%s&date_to_view=%s" % ( site_address, site_localtime.strftime("%m/%d/%Y") )

    return site_address

def get_daily_ihope_page(gpstime, pages_location = "https://ldas-jobs.ligo.caltech.edu/~cbc/ihope_daily"):
    utctime = XLALGPSToUTC(LIGOTimeGPS(gpstime, 0))
    return "%s/%s/%s/" %(pages_location, time.strftime("%Y%m", utctime), time.strftime("%Y%m%d", utctime))


def create_hyperlink(address, link):
    return '<a href="%s">%s</a>' % (address, link)


def create_filter( connection, tableName, param_name = None, param_ranges = None, 
    exclude_coincs = None, include_only_coincs = None, sim_tag = 'ALLINJ'):
    """
    Strings together param_name, param_ranges, exclude/include_only_coincs, and
    sim_tag options into a filter string that can be stuck in a sqlite WHERE clause.
    """
    in_this_filter = ''
    
    # Get param and param-ranges if specified
    if param_name is not None:
        param_name = sqlutils.validate_option(param_name)
        param_filters = sqlutils.parse_param_ranges( tableName, param_name, 
            param_ranges, verbose = verbose ).get_param_filters()
        # since want triggers that fall within all the parameters, concatenate
        # all param ranges
        param_filters = '\n\t\tOR '.join( param_filters )
        in_this_filter = ''.join([ in_this_filter, '\n\tAND (\n\t\t', param_filters, '\n\t)' ])
    
    # Get exclude_coincs list if specified
    if exclude_coincs is not None:
        exclude_coinc_filters = sqlutils.parse_coinc_options( exclude_coincs, 
            verbose = verbose ).get_coinc_filters()
        # concatenate exclude_coinc_filters
        exclude_coinc_filters = '\n\t\tOR '.join( exclude_coinc_filters )
        # add to in_this_filter
        in_this_filter = ''.join([ in_this_filter, '\n\tAND NOT (\n\t\t', exclude_coinc_filters, '\n\t)' ]) 
    
    # Get include_only_coincs list if specified
    if include_only_coincs is not None:
        include_coinc_filters = sqlutils.parse_coinc_options( include_only_coincs, 
            verbose = verbose ).get_coinc_filters()
        # concatenate include_coinc_filters
        include_coinc_filters = '\n\t\tOR '.join( include_coinc_filters )
        # add to in_this_filter
        in_this_filter = ''.join([ in_this_filter, '\n\tAND (\n\t\t', include_coinc_filters, '\n\t)' ])
    
    # if sim-tag specified add the sim-tag to the filter
    if sim_tag != 'ALLINJ':
        # create a map between sim_proc_id and sim-tag
        sim_map = sqlutils.sim_tag_proc_id_mapper( connection )
        # check that sim_tag is in the the map
        sim_tag = sqlutils.validate_option(sim_tag, lower = False).upper()
        if sim_tag not in sim_map.tag_id_map.keys():
            raise ValueError, "sim-tag %s not found in database" % sim_tag
        # create the filter
        connection.create_function( 'get_sim_tag', 1, sim_map.get_sim_tag )
        sim_filter = ''.join(['get_sim_tag(experiment_summary.sim_proc_id) == "', sim_tag, '"' ])
        # add to in_this_filter
        in_this_filter = ''.join([ in_this_filter, '\n\tAND ', sim_filter ])

    return in_this_filter


# =============================================================================
#
#                           Library API
#
# =============================================================================

def printsims(connection, recovery_table, simulation_table, ranking_stat, rank_by, comparison_datatype,
    param_name = None, param_ranges = None, exclude_coincs = None, include_only_coincs = None,
    sim_tag = 'ALLINJ', rank_range = None, convert_durations = 's',
    daily_ihope_pages_location = 'https://ldas-jobs.ligo.caltech.edu/~cbc/ihope_daily' verbose = False):

    # check and format options appropriately
    simulation_table = sqlutils.validate_option(simulation_table)
    recovery_table = sqlutitls.validate_option(recovery_table)
    ranking_stat = sqlutils.validate_option(ranking_stat)
    comparison_datatype = sqlutils.validate_option(comparison_datatype)
    convert_durations = sqlutils.validate_option(convert_durations)

    if not ranking_stat.startswith(recovery_table):
        ranking_stat = '.'.join([recovery_table, ranking_stat])

    #
    #   Set up sim_rec_map table
    #
    sqlutils.create_sim_rec_map_table(connection, simulation_table, recovery_table, ranking_stat)
    
    #
    #   Set recovery table filters
    #
    in_this_filter = create_filter(connection, recovery_table, param_name = param_name, param_ranges = param_ranges,
        exclude_coincs = exclude_coincs, include_only_coincs = include_only_coincs, sim_tag = sim_tag)
    
    
    #
    #   Initialize ranking. Statistics for ranking are collected from non-injections
    #   in the recovery table.
    #
    if verbose:
        print >> sys.stdout, "Getting statistics for ranking..."
    ranker = sqlutils.rank_stats(recovery_table, ranking_stat, rank_by)
    # add requirement that stats not be found in the sim_rec_table to in_this_filter
    rank_filter = ''.join([
            recovery_table, '''.coinc_event_id NOT IN (
                SELECT
                    rec_id
                FROM
                    sim_rec_map)
            AND
                experiment_summary.datatype == "''', comparison_datatype, '"'])
    
    if in_this_filter == '':
        rank_filter = '\n\tAND '.join([ in_this_filter, rank_filter ])
    else:
        rank_filter = '\n\t'.join([ sqlutils.join_experiment_tables_to_coinc_table(recovery_table), 'WHERE', rank_filter ])
    
    ranker.populate_stats_list(connection, limit = None, filter = rank_filter)
    connection.create_function( 'rank', 1, ranker.get_rank )
    
    
    # Now apply the filter to the sim_rec_map table: this will delete all sim/rec maps where the simulation id is
    # mapped to a recovered event that falls outside the filter, even if that particular sim/rec map is in the
    # filter. For example, if the filter is recovery_table.combined_far != 0., and there are two entries in the
    # sim_rec_map table sharing the same sim_id:
    #   sim_id:0 | rec_id:0 | rec_id's combined_far = 0.
    #   sim_id:0 | rec_id:1 | rec_id's combined_far = 1.2
    # both entries will get deleted even though the second entry's combined_far is not 0.
    if in_this_filter != '':
        # join the needed tables to in_this_filter
        in_this_filter = ''.join([ sqlutils.join_experiment_tables_to_coinc_table(recovery_table), "\n    WHERE\n\t", re.sub('\n\tAND', '', in_this_filter, 1) ])
        sqlscript = ''.join([ """
            CREATE TEMP TABLE del_ids AS
                SELECT
                    sim_id AS del_id
                FROM
                    sim_rec_map
                WHERE
                    rec_id NOT IN (
                    SELECT
                        """, recovery_table, """.coinc_event_id
                    FROM
                        """, recovery_table, """
                    """, in_this_filter, """
                        AND experiment_summary.datatype == "simulation"
                    );
        
            DELETE FROM
                sim_rec_map
            WHERE
                sim_id IN (
                    SELECT
                        del_id
                    FROM
                        del_ids );
        
            DROP TABLE del_ids;""" ])
        connection.cursor().executescript(sqlscript)
    
    #
    #   Set other needed functions
    #
    
    # establish what units will be converting duration to
    def convert_duration( duration ):
        return sqlutils.convert_duration( duration, convert_durations )
    connection.create_function( 'convert_duration', 1, convert_duration )
    
    # Get range ranks
    if rank_range is not None:
        rank_range_parser = sqlutils.parse_param_ranges( 'rank(sim_rec_map', 'ranking_stat)',
            rank_range, verbose = verbose )
    
    #
    #   Create and prepare the SelectedFoundTable to store summary information
    #
    
    # Get recovery table and simulation table column names from database
    simulation_table_columns = sqlutils.get_column_names_from_table( connection, simulation_table )
    recovery_table_columns = sqlutils.get_column_names_from_table( connection, recovery_table )
    
    # Get list of column name for injected parameters
    injected_cols = []
    for col in simulation_table_columns:
        if col == 'simulation_id':
            injected_cols.append('simulation_id')
        else:
            injected_cols.append('injected_'+col)
    injected_cols.extend(['injected_end_time', 'injected_end_time_ns', 'injected_end_time_utc__Px_click_for_daily_ihope_xP_'])
    
    # Get list of column names from the recovery table
    recovered_cols = []
    for col in recovery_table_columns:
        if col == 'coinc_event_id':
            recovered_cols.append(u'coinc_event_id')
        else:
            recovered_cols.append(u'recovered_'+col)
    # generate list of column names for the summary table
    column_names = injected_cols + recovered_cols
    # add instruments on, duration, mini_followups
    rankname = 'rank_in_' + comparison_datatype.strip().lower() + '_using_' + ranking_stat.split('.')[-1]
    durname = ''.join([ 'simulation', u'_duration__Px_', convert_durations, '_xP_' ])
    column_names.extend([u'instruments_on', durname, u'mini_followup' ])
    column_names.extend( [ rankname, 'recovered_match_rank', 'instruments_on', 'elogs', durname, 'mini_followup' ] )
    
    #
    # define needed tables
    #
    class tmpSimTable(table.Table):
        tableName = "tmp_sim_table:table"
        validcolumns = dict([ [col, sqlutils.get_col_type(simulation_table, col)] for col in simulation_table_columns ])
    class tmpSim(object):
        __slots__ = tmpSimTable.validcolumns.keys()
        def get_pyvalue(self):
            return generic_get_pyvalue(self)
    
    class tmpRecTable(table.Table):
        tableName = "tmp_rec_table:table"
        validcolumns = dict([ [col, sqlutils.get_col_type(recovery_table, col)] for col in recovery_table_columns ])
    class tmpRec(object):
        __slots__ = tmpRecTable.validcolumns.keys()
        def get_pyvalue(self):
            return generic_get_pyvalue(self)
    
    class SelectedFoundTable(table.Table):
        tableName = "selected_found_injections:table"
        validcolumns = {}
        for col_name in column_names:
            if 'rank_in_' in col_name:
                validcolumns[col_name] = "int_4u"
            elif '_duration_' in col_name:
                validcolumns[col_name] = "real_8"
            elif 'instruments_on' == col_name:
                validcolumns[col_name] = lsctables.ExperimentTable.validcolumns['instruments']
            elif col_name == 'injected_end_time' or col_name == 'injected_end_time_ns':
                validcolumns[col_name] = "int_4s"
            elif 'injected_' in col_name or col_name == 'simulation_id':
                validcolumns[col_name] = sqlutils.get_col_type(simulation_table, re.sub('injected_', '', col_name))
            elif 'recovered_' in col_name or col_name == 'coinc_event_id':
                validcolumns[col_name] = sqlutils.get_col_type(recovery_table, re.sub('recovered_', '', col_name))
            # if custom columns exist in the database, just set them to lstrings
            else:
                validcolumns[col_name] = "lstring"
        # add FAP columns
        validcolumns['recovered_fap'] = "real_8"
        validcolumns['recovered_fap_1yr'] = "real_8"
    
    class SelectedFound(object):
        __slots__ = SelectedFoundTable.validcolumns.keys()
    
        def get_pyvalue(self):
            return generic_get_pyvalue(self)
    
    # connect the rows to the tables
    tmpSimTable.RowType = tmpSim
    tmpRecTable.RowType = tmpRec
    SelectedFoundTable.RowType = SelectedFound
    
    #
    #   Get the Data
    #
    tmp_sim_table = lsctables.New(tmpSimTable) 
    tmp_rec_table = lsctables.New(tmpRecTable) 
    tmp_sftable = lsctables.New(SelectedFoundTable)
    prior_sim_id = ''
    group_rank = None
    current_match_rank = 1
    sqlquery = ''.join(["""
        SELECT
            """, simulation_table, """.*,
            """, recovery_table, """.*,
            rank(sim_rec_map.ranking_stat),
            NULL AS match_rank,
            experiment.instruments,
            convert_duration(experiment_summary.duration),
            NULL AS mini_followup
        FROM
            sim_rec_map
        JOIN
            """, ', '.join([simulation_table, recovery_table]), """, experiment, experiment_summary, experiment_map ON (
            sim_rec_map.sim_id == """, simulation_table, """.simulation_id AND
            sim_rec_map.rec_id == """, recovery_table, """.coinc_event_id AND
            sim_rec_map.rec_id == experiment_map.coinc_event_id AND
            experiment_map.experiment_summ_id == experiment_summary.experiment_summ_id AND
            experiment_summary.experiment_id == experiment.experiment_id)
        ORDER BY
            sim_rec_map.sim_id, sim_rec_map.ranking_stat """, rank_by])
    
    if verbose:
        print >> sys.stdout, "Getting coincs..."
        print >> sys.stdout, "SQLite query used is:"
        print >> sys.stdout, sqlquery
    
    for values in connection.cursor().execute( sqlquery ).fetchall():
        # sort the data
        tmp_sim_row = tmpSim()
        tmp_rec_row = tmpRec()
        [ setattr(tmp_sim_row, column, values[ii]) for ii, column in enumerate(simulation_table_columns) ]
        [ setattr(tmp_rec_row, column, values[ii+1+jj]) for jj, column in enumerate(recovery_table_columns) ]
        # figure out the rank
        this_inj_rank = values[-5]
        this_sim_id = tmp_sim_row.simulation_id
        this_ranking_stat = getattr(tmp_rec_row, ranking_stat.split('.')[-1])
        if this_sim_id == prior_sim_id and this_ranking_stat != prior_ranking_stat:
            current_match_rank += 1
        elif this_sim_id != prior_sim_id:
            current_match_rank = 1
            group_rank = this_inj_rank
        prior_sim_id = this_sim_id
        prior_ranking_stat = this_ranking_stat 
        # only store data if the group_rank falls in the desired rank ranges
        if rank_range and rank_range_parser.group_by_param_range(group_rank) is None:
            continue
        on_instruments = lsctables.instrument_set_from_ifos(values[-3])
        duration = values[-2]
        # now that have all the information from this row, create a row for the selected found table
        sfrow = SelectedFound()
        # set the ranks
        setattr(sfrow, rankname, group_rank)
        sfrow.recovered_match_rank = current_match_rank
        # set the injected parameters
        use_this_site = sorted(on_instruments)[0][0].lower()
        sfrow.injected_end_time = getattr( tmp_sim_row, use_this_site+'_end_time' )
        sfrow.injected_end_time_ns = getattr( tmp_sim_row, use_this_site+'_end_time_ns' )
        for col in simulation_table_columns:
            if col == "simulation_id":
                sfrow.simulation_id = tmp_sim_row.simulation_id 
            else:
                setattr(sfrow, 'injected_'+col, getattr( tmp_sim_row, col) )
        # set the recovered parameters
        for col in recovery_table_columns:
            if col == "coinc_event_id":
                sfrow.coinc_event_id = tmp_rec_row.coinc_event_id
            else:
                setattr(sfrow, 'recovered_'+col, getattr( tmp_rec_row, col) )
        # calculate and add faps
        if sfrow.recovered_combined_far is not None:
            t_in_s = float(values[-2]) / sqlutils.convert_duration(1, convert_durations)
            sfrow.recovered_fap = 1 - math.exp(-sqlutils.convert_duration(t_in_s, 'yr') * sfrow.recovered_combined_far)
            sfrow.recovered_fap_1yr = 1 - math.exp(-sfrow.recovered_combined_far)
        else:
            sfrow.recovered_fap = None
            sfrow.recovered_fap_1yr = None
        # set elog page
        elog_pages = [(ifo, get_elog_page(ifo, sfrow.injected_end_time)) for ifo in on_instruments]
        sfrow.elogs = ','.join([ create_hyperlink(elog[1], elog[0]) for elog in sorted(elog_pages) ])
        # set daily_ihope page
        end_time_utc = format_end_time_in_utc( sfrow.injected_end_time ) 
        daily_ihope_address = get_daily_ihope_page(sfrow.injected_end_time, pages_location = daily_ihope_pages_location)
        sfrow.injected_end_time_utc__Px_click_for_daily_ihope_xP_ = create_hyperlink( daily_ihope_address, end_time_utc ) 
        # set any other info
        sfrow.instruments_on = ','.join(sorted(on_instruments))
        sfrow.mini_followup = None
        setattr(sfrow, durname, duration)
    
        # add the row
        tmp_sftable.append(sfrow)
    
    # Re-sort the sftable by rank, recovered_match_rank
    sftable = lsctables.New(SelectedFoundTable)
    for sfrow in sorted([ row for row in tmp_sftable ], key = lambda row: getattr(row, rankname) ):
        if sfrow.simulation_id not in [row.simulation_id for row in sftable]:
            sftable.append(sfrow)
            sftable.extend(sub_row for sub_row in sorted([row for row in tmp_sftable
                if row.simulation_id == sfrow.simulation_id
                and row.coinc_event_id != sfrow.coinc_event_id],
                key = lambda row: row.recovered_match_rank))
    
    return sftable
