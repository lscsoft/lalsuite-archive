#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

'''
A collection of utilities to assist in carrying out operations on a SQLite
database containing lsctables.
'''

try:
    import sqlite3
except ImportError:
    # pre 2.5.x
    from pysqlite2 import dbapi2 as sqlite3

import sys
import re
import os
import bisect

from glue.ligolw import dbtables

__author__ = "Collin Capano <cdcapano@physics.syr.edu>"
__date__ = "$Date$" 
__version__ = "$Revision$"



# =============================================================================
#
#                           Generic Utilities
#
# =============================================================================

# Following utilities can be used with any coinc_table

class parse_param_ranges:
    
    param = None
    param_ranges = []

    def __init__( self, table_name, table_param, param_ranges_opt, verbose = False ):
        """
        Parse --param-ranges option. Creates self.param which is the table_name and
        the table_param appended together (with a '.') and self.param_ranges, which is
        a list of tuples that give the lower parameter value, whether it is an open or
        closed boundary, and the same for the upper parameter. For example, if 
        table_name is coinc_inspiral, table_param is mchirp and param_ranges_opt 
        is '[2,8);[8,17]' will get:
        self.param = 'coinc_inspiral.mchirp'
        self.param_ranges = 
          [ ( ('>=',2.0), ('<',8.0) ),
            ( ('>=',8.0), ('<=', 17.0) ) ]

        @table_name: Name of coinc_table in which the desired parameter is a column.
        @table_param: Parameter in the table on which to separate rows.
        @param_ranges_opt: string from the --param-ranges option. Param-ranges must
         follow these format rules:
            * A '(' or ')' implies an open boundary, a '[' or ']' a closed boundary.
            * To specify multiple ranges, separate each range by a ';'.
        @verbose: be verbose
        """
        if verbose:
            print >> sys.stderr, "Parsing param-ranges..."
        
        # check that table_name and table_param have no spaces in them
        if len( table_name.split(' ') ) > 1:
            raise ValueError, "table_name cannot have spaces in it."
        if len ( table_param.split(' ') ) > 1:
            raise ValueError, "table_param cannot have spaces in it."

        # make param unique by appending table_name to the param_name
        self.param = '.'.join([ table_name, table_param ])

        ranges = param_ranges_opt.split(';')

        for this_range in ranges:

            # get lower-bound
            lowerparam = this_range.split(',')[0].strip()
            # check if lower boundary open or closed
            if lowerparam.find('[') != -1:
                lowerbndry = '>='
                lowerparam = float( lowerparam.lstrip('[') )
            elif lowerparam.find('(') != -1:
                lowerbndry = '>'
                lowerparam = float( lowerparam.lstrip('(') )
            else:
                raise ValueError, "Parameter range %s not formatted correctly" % this_range
  
            # get upper-bound (similar to lower bound method)
            upperparam = this_range.split(',')[1].strip()
            if upperparam.find(']') != -1:
                upperbndry = '<='
                upperparam = float( upperparam.rstrip(']') )
            elif upperparam.find(')') != -1:
                upperbndry = '<'
                upperparam = float( upperparam.rstrip(')') )
            else:
                raise ValueError, "Parameter range %s not formatted correctly" % this_range

            # add param to filters
            self.param_ranges.append( ( (lowerbndry, lowerparam), (upperbndry, upperparam) ))

        if verbose:
            print >> sys.stderr, "done."
    

    def get_param_name( self ):
        return self.param


    def get_param_ranges( self ):
        return self.param_ranges


    def get_param_filters( self ):
        """
        Converts param_ranges into a list of strings that can be used in 
        a SQLite WHERE clause. For example, if table_name is coinc_inspiral, 
        table_param is mchirp and param_ranges_opt is '[2,8);[8,17]' the 
        elements in the returned list will be:
        ['coinc_inspiral.mchirp >= 2.0 AND coinc_inspiral.mchirp < 8.0',
         'coinc_inspiral.mchirp >= 8.0 AND coinc_inspiral.mchirp <= 17.0']
        """
        self.param_filters = []
        # construct paramfilter for SQL statement
        for range in self.param_ranges:
            lowerbndry = range[0][0]
            lowerparam = str( range[0][1] )
            upperbndry = range[1][0]
            upperparam = str( range[1][1] )
            self.param_filters.append( ' '.join([ '(', self.param, lowerbndry, lowerparam, 
              'AND', self.param, upperbndry, upperparam, ')' ]) )

        return self.param_filters


    def group_by_param_range( self, param_value ):
        """
        Takes in a value and returns a number corresponding to
        which value param_range it falls in.
        """
        for n, range in enumerate(self.param_ranges):
            # set boundry conditions and parameters
            lowerbndry = range[0][0]
            lowerparam = range[0][1]
            upperbndry = range[1][0]
            upperparam = range[1][1]
            # the following works by checking what the boundaries are
            # and then checking if the param value is within those boundaries:
            # if [a,b]
            if ((lowerbndry, upperbndry) == ('>=', '<=')) and \
               (param_value >= lowerparam and param_value <= upperparam):
                return n
            # if (a,b]
            if ((lowerbndry, upperbndry) == ('>', '<=')) and \
               (param_value > lowerparam and param_value <= upperparam):
                return n
            # if [a,b)
            if ((lowerbndry, upperbndry) == ('>=', '<')) and \
               (param_value >= lowerparam and param_value < upperparam):
                return n
            # if (a,b)
            if ((lowerbndry, upperbndry) == ('>', '<')) and \
               (param_value > lowerparam and param_value < upperparam):
                return n

        # if get to here, param_value falls outside all the ranges; 
        # just return None
        return None


class parse_coinc_options:

    coinc_types = {}
    
    def __init__( self, coincs_opt, verbose = False):
        """
        Parses --exclude-coincs and --include-coincs options. The class doesn't
        care whether it's --include or exclude; it just takes in the input and 
        creates self.coinc_types, which is a dictionary of coinc types in which the
        keys are the type of time and the values are the coincidence type. For example,
        if either --include-coincs or --exclude coincs is set to "[h2,l1 in h1,h2,l1]"
        self.coinc_types will be:
            coinc_types[ frozenset(H1,H2,L1) ] = set(H2,L1)
            
        
        @coincs_opt: the input from either --exclude-coincs or --include-coincs.
         This input must follow these format rules (the following can be copied
         into the help message for --input/exclude-coincs opts):
            * Coinc-types and detector time must be separated by 
            an ' in '. When specifying a coinc_type or detector    
            time, detectors and/or ifos must be separated by 
            commas, e.g. 'H1,L1' not 'H1L1'.                     
            * To specify multiple coinc-types in one type of time,
            separate each coinc-type by a '+', e.g., 
            '[H1,H2 + H2,L1 in H1,H2,L1]'.                        
            * To specify all the coincs in a detector time 
            or a specific coinc-type in all times, use 'ALL'. E.g.,
            to exclude/include all H1,H2 triggers, use '[H1,H2 in ALL]' 
            or to exclude/include all H2,L1 time use '[ALL in H2,L1]'.   
            * To specify multiple exclusions, separate each 
            bracket by a ';'.                              
            * Order of the instruments nor case of the letters 
            matter. So if your pinky is broken and you're      
            dyslexic you can type '[h2,h1 in all]' without a 
            problem.
        @verbose: be verbose.
        """

        if verbose:
            print >> sys.stderr, "Parsing coinc options..."

        for rule in coincs_opt.split(';'):
            rule = rule.strip().lstrip('[').rstrip(']').upper()

            # get coinc_instruments, instruments_on 
            [ coinc_instruments, instruments_on ] = rule.split(' IN ')
            instruments_on = instruments_on.strip()
            coinc_instruments = coinc_instruments.strip()

            # Parse instruments_on
            if instruments_on != 'ALL':
                instruments_on = frozenset(sorted([instrument.strip() for instrument in instruments_on.split(',')]))
                if instruments_on not in self.coinc_types:
                    # sanity check
                    if len(instruments_on) <= 1:
                        raise ValueError, "Must delimit instruments by commas."
                    # following is to try to protect against injection attacks
                    for instrument in instruments_on:
                        if len(instrument.split(' ')) > 1:
                            raise ValueError, "Instrument names cannot have spaces in them."
                    self.coinc_types[ instruments_on ] = []
            elif 'ALL' not in self.coinc_types:
                self.coinc_types[ 'ALL' ] = []

            # Parse coinc_instruments
            if coinc_instruments != 'ALL':
                for coinc_instset in coinc_instruments.split('+'):
                    coinc_instset = set(sorted( instrument.strip() for instrument in coinc_instset.split(',') ))
                    if coinc_instset not in self.coinc_types[ instruments_on ]:
                        # sainity check
                        if len(coinc_instset) <= 1:
                            raise ValueError, "Must delimit instruments by commas."
                        for instrument in coinc_instset:
                            if len(instrument.split(' ')) > 1:
                                raise ValueError, "Instrument names cannot have spaces in them."
                        # add instset to coinc_types
                        self.coinc_types[ instruments_on ].append( coinc_instset )
            else:
                self.coinc_types[ instruments_on ] = ['ALL']


    def get_coinc_types( self ):
        return self.coinc_types


    def get_coinc_filters( self ):
        """
        Converts self.coinc_types to a list of strings that can be used
        in a SQLite WHERE clause to filter coincs by coinc_type,
        by coinc_instruments (which is stored in the coinc_inspiral table)
        in instruments_on (which is stored in the experiment table).
        """
        self.coinc_filters = []
        # import ifos_from_instrument_set in lsctables for converting
        # instrument sets in self.coinc_types to strings
        from glue.ligolw.lsctables import ifos_from_instrument_set
        
        # cycle through instruments_on in coinc_types
        for instruments_on in self.coinc_types:
            this_coincfilter = ''
            if instruments_on != 'ALL':
                this_coincfilter = ''.join([
                    'experiment.instruments == "', ifos_from_instrument_set(instruments_on), '"' ])
                # now cycle through coinc_instruments in self.coinc_types[ instruments_on ],
                # concatenate each coinc_instruments set with an OR;
                # append the concatenated string this_coincfilter with an AND
                if 'ALL' not in self.coinc_types[ instruments_on ]:
                    this_coincfilter = ' '.join([ this_coincfilter, 'AND (' ])
                    for coinc_instruments in self.coinc_types[ instruments_on ]:
                        this_coincfilter = ''.join([ this_coincfilter,
                            ' coinc_inspiral.ifos == "', ifos_from_instrument_set(coinc_instruments), '"', ' OR' ])
                    # strip the last 'OR' and replace with a ')' to close out the coinc_instruments
                    this_coincfilter = this_coincfilter.rstrip('OR') + ')'
            # if instruments_on is 'ALL', just add what coincs to filter
            elif instruments_on == 'ALL' and 'ALL' not in self.coinc_types[ instruments_on ]:
                for coinc_instruments in self.coinc_types[ instruments_on ]:
                    this_coincfilter = ''.join([ this_coincfilter,
                        ' coinc_inspiral.ifos == "', ifos_from_instrument_set(coinc_instruments), '"', ' OR' ])
                # strip the last 'OR'
                this_coincfilter = this_coincfilter.rstrip('OR')

            self.coinc_filters.append( ''.join([ '(', this_coincfilter, ')' ]) )

        return self.coinc_filters
            


def del_rows_from_table( connection, del_table, del_table_id, join_conditions, del_filters = None, save_filters = None, verbose = False ):
    """
    Deletes triggers from any specified table in the del_table option.
    @connection: DBTables connection to a database
    del_table: Any coinc_table (coinc_inspiral, sngl_inspiral, coinc_event,
    etc.) from which to delete triggers.
    @del_table_id: name of ID column in the del_table on which will be deleting
    triggers.
    @join_conditions: SQLite string that draws connections between different
     coinc_tables. Must be of format 'JOIN table1 ON table1-link-to-other-table 
     JOIN table2 ON table2-link', etc.
    @del_filter: List of filters. Triggers that fall within will be deleted.
    @save_filter: List of filters. Triggers that fall within will NOT be deleted.

    NOTE: Save filters will override del_filters if they overlap. For example,
    say del filter species H1,H2 triggers in H1,H2,L1 time and save filters are
    for triggers with mchirp between 2 and 8. Then all triggers with chirp mass
    between 2 and 8 will be saved, even if they are H1,H2. All other H1,H2
    triggers will be deleted. What this means is if you want to do a global
    delete -- say you wan to delete all H1,H2 triggers, do not specify a
    save_filter that overlaps with it.
    """
    # append table name to table_id to ensure uniqueness
    del_table_id = '.'.join([ del_table, del_table_id])
  
    # set where clause to be used in delete statement based on del and save
    # filters
    where_clause = ''
    if del_filters:
        del_filters = [ ''.join([ '(', filter, ')' ]) for filter in del_filters ]
        del_filters = ' OR '.join( del_filters )
        where_clause = ' '.join([ 'WHERE (', del_filters, ')' ])
    if save_filters:
        save_filters = [ ''.join([ '(', filter, ')' ]) for filter in save_filters ]
        save_filters = ' OR '.join( save_filters )
        if not del_filters:
            where_clause = ' '.join([ 'WHERE NOT (', save_filters, ')' ])
        else:
            where_clause = ' '.join([ where_clause, 'AND NOT (', save_filters, ')' ])
    # if no filters, warn user
    if not where_clause:
        print >> sys.stderr, '''WARNING: No filters specified in delete statement.
          Deleting all rows from %s''' % del_table
    elif verbose:
        print >> sys.stderr, "Deleting rows from %s table %s..." % (del_table, where_clause)
  
    sqlquery = ' '.join([
          'DELETE FROM', del_table,
          'WHERE', del_table_id, 'IN (', 
              'SELECT', del_table_id,
              'FROM', del_table, join_conditions,
                 where_clause, ')' ])
    connection.cursor().execute( sqlquery )
    connection.commit()

    if verbose:
        print >> sys.stderr, "done."


def get_column_names_from_table( connection, table_name ):
    """
    Gets the column names from a table and returns them as a list.
    """
    sqlquery = ''.join(['PRAGMA table_info(', table_name, ')' ])
    column_names = [ name[1] for name in connection.cursor().execute( sqlquery).fetchall() ]
    
    return column_names


def convert_duration( duration, convert_to ):
    """
    Converts durations stored in the experiment_summary_table from seconds 
    to other units of time.

    @duration: duration to convert; assumed to be a float or long in seconds
    @convert_to: the unit to convert to. Options are:
        's': to seconds - will just divide by 1.
            This can be useful if need to convert
            the duration from a long int to a float.
        'min': to minutes - will divide by 60.
        'hr': to hours - will divide by 3600.
        'd': to days - will divide by 86400.
        'yr': to years - will divide by 31557600. 
            This is the Julian year, which is the
            accepted astronomical year
    """
    if not duration:
        return 0. 
    if convert_to == 's':
        return duration / 1.
    elif convert_to == 'min':
        return duration / 60.
    elif convert_to == 'hr':
        return duration / 3600.
    elif convert_to == 'd':
        return duration / 86400.
    elif convert_to == 'yr':
        return duration / 31557600.
    else:
        raise ValueError, "Unrecognized unit."

def get_next_id(connection, table, id_column):
    """
    Gets the next available id in the specified id_column in the specified table.
    """
    sqlquery = ' '.join(['SELECT', id_column, 'FROM', table ])
    ids = [id[0] for id in connection.cursor().execute(sqlquery)]
    idnums = [int(id.split(':')[2]) for id in ids]
    new_idnum = max(idnums) + 1
    new_id = ':'.join([ids[0].split(':')[0], ids[1].split(':')[1], str(new_idnum)])

    return new_id
    

class Summaries:
    """
    This class stores information about the foreground and background in a 
    database for making calculation of uncombined fars and combined fars quick 
    and efficient.
    
    bkg_stats groups triggers by experiment_id, ifos, and param_group 
    (param_group is an arbitrary integer representing the param bin, e.g., 
    mchirp [3.48,7.4), to which a trigger belongs; if no binning is done, then
    it is 0 for all triggers). It stores ALL the triggers in all the time 
    slides (except zero-lag) within that group.

    sngl_slide_stats groups triggers by experiment_id, experiment_summ_id, ifos, and 
    param_group. It therefore groups all triggers within each time slide 
    separately. It is used to subtract triggers within the same slide when calculating
    uncombined fars for the background. Therefore, it only stores slide triggers;
    for any zero-lag datatype sngl_slide_stats is just an empty list.

    frg_durs stores the duration for each experiment_summ_id. It's keys are 
    [experiment_id][experimen_summ_id].

    bkg_durs stores the background duration for each time-slide and zero-lag, 
    i.e., for each experiment_summ_id. This is the sum of all other slide
    datatypes sharing the same experiment_id except for the given slide. 

    max_bkg_fars stores the maximum background fars of all the categories 
    within each time slide. It's keys are (experiment_summ_id, ifo_group). 
    The maximum background far is just the total number of triggers within a 
    category divided by the background duration for that time slide.
    If opts.combine_fars is set to across_all a category is defined by the 
    param bin in which a trigger exists and the ifos that took part in the 
    trigger. So, if there are three param bins and we've excluded H2,L1 triggers 
    in H1,H2,L1 time, then there are 6 categories for H1,H2,L1 time: three param
    bins each for H1,L1 and H1,H2,L1 coincident triggrs. Thus, ifo_group will 
    be set to "ALL_IFOS" and there will be 6 max_bkg_fars stored for each 
    experiment_summ_id in triple time.
    If opts.combine_fars is set to across_param_only, then a category is 
    defined only by the param bins; ifo coincidences are treated as 
    separate experiments. Thus, ifo_group will be set to whatever
    coinc. trigger we are considering and there will only be 3 max_bkg_fars 
    stored for that entry.

    zero_lag_ids stores the esid and datatype (all_data, playground, or exclude_play)
    of the zero-lag slide for an experiment:
        zero_lag_ids[ experiment_id ][ datatype ] = experiment_summ_id
    """
    def __init__(self):
        self.bkg_stats = {}
        self.sngl_slide_stats = {}
        self.frg_durs = {}
        self.bkg_durs = {}
        self.max_bkg_fars = {}
        self.zero_lag_ids = {}

    def add_to_bkg_stats(self, experiment_id, experiment_summ_id, ifos, param_group, stat):
        """
        Adds a stat to bkg_stats and sngl_slide_stats. What stat is added is determined on the command
        line by the ranking-stat option.
        """
        if experiment_summ_id in self.zero_lag_ids[experiment_id].values():
            self.sngl_slide_stats[(experiment_id, experiment_summ_id, ifos, param_group)] = []
        else:
            if (experiment_id, ifos, param_group) not in self.bkg_stats:
                self.bkg_stats[(experiment_id, ifos, param_group)] = []
            if (experiment_id, experiment_summ_id, ifos, param_group) not in self.sngl_slide_stats:
                self.sngl_slide_stats[(experiment_id, experiment_summ_id, ifos, param_group)] = []
            self.bkg_stats[(experiment_id, ifos, param_group)].append( stat )
            self.sngl_slide_stats[(experiment_id, experiment_summ_id, ifos, param_group)].append(stat)

    def sort_bkg_stats(self):
        """
        Sorts each list in bkg_stats and sngl_slide_stats from smallest to largest value.
        """
        for thislist in self.bkg_stats.values():
            thislist.sort()
        for thislist in self.sngl_slide_stats.values():
            thislist.sort()

    def append_zero_lag_id(self, experiment_id, zero_lag_esid, datatype):
        """
        Adds a zero_lag_id to the zero_lag_ids dictionary.
        """
        self.zero_lag_ids[experiment_id] = dict({ datatype: zero_lag_esid })

    def append_duration(self, experiment_id, experiment_summ_id, duration):
        """
        Adds a duration to frg_durs.
        """
        if experiment_id not in self.frg_durs:
            self.frg_durs[experiment_id] = {}
        self.frg_durs[experiment_id][experiment_summ_id] = duration

    def calc_bkg_durs(self):
        """
        Sums the background durs for each time-slide (experiment_summ_id).
        """
        for eid in self.frg_durs:
            for this_esid in self.frg_durs[eid]:
                self.bkg_durs[this_esid] = sum([self.frg_durs[eid][bkg_esid] for bkg_esid in self.frg_durs[eid].keys() if bkg_esid != this_esid and bkg_esid not in self.zero_lag_ids[eid].values()])

    def append_max_bkg_far(self, experiment_summ_id, ifo_group, max_bkg_far):
        """
        Adds a max_bkg_far to the appropiate list; lists are grouped by 
        experiment_summ_id and ifo_group. If one wants to combined fars across 
        param_bins and coincident_ifos (as was done in the low-mass S51yr and 
        12-18 month analyses), ifo_group should be set to "ALL_IFOS".
        """
        if (experiment_summ_id, ifo_group) not in self.max_bkg_fars:
            self.max_bkg_fars[(experiment_summ_id, ifo_group)] = []
        self.max_bkg_fars[(experiment_summ_id, ifo_group)].append(max_bkg_far)

    def sort_max_bkg_fars(self):
        """
        Sorts the max_bkg_fars lists from smallest to highest values.
        """
        for thislist in self.max_bkg_fars.values():
            thislist.sort()

    def calc_ufar_by_max(self, eid, esid, ifos, param_group, stat):
        """
        Calculates the uncombined false alarm rate for a trigger by counting 
        the number of background triggers in the same category as it that have
        a stat value greater than or equal to the trigger's stat value and 
        dividing by the background duration for that slide.
        To do this quickly, bisect.bisect_left is used (see python 
        documentation for more info) on the bkg_stats list. Since bkg_stats 
        contains all the triggers in all the slides for some experiment_id,
        this will result in counting the triggers that are in the same slide
        (given by the esid) as the trigger we are considering (except for zero-lag).
        To correct for this, the trigger's place in it's sngl_slide_stats list is
        subtracted from this value. The "background" considered for some trigger is
        therefore all the triggers sharing the same experiment_id, excluding
        zero-lag triggers and triggers in the same time-slide as the trigger. This
        means that uncombined far for non-zero-lag triggers will use one less time
        slide than zero-lag triggers.
        """
        return (\
            ( len(self.bkg_stats[(eid, ifos, param_group)]) - bisect.bisect_left(self.bkg_stats[(eid, ifos, param_group)], stat) ) \
            - \
            ( len(self.sngl_slide_stats[(eid, esid, ifos, param_group)]) - bisect.bisect_left(self.sngl_slide_stats[(eid, esid, ifos, param_group)], stat) ) \
            ) / self.bkg_durs[esid]

    def calc_ufar_by_min(self, eid, esid, ifos, param_group, stat):
        """
        Same as calc_ufar_by_max, except that the uncombined far is calculated
        by counting background triggers that have a stat value less than or 
        equal to the given stat. (Done by using bisect.bisect_right as opposed to 
        len(list) - bisect.bisect_left).
        Note: if stat is 0, will just return 0. This is because a 0 when caclulating
        FARs by minimum value is equivalent to inf. when caclulating FARs by maximum
        value.
        """
        if stat == 0.:
            return stat

        return ( \
            bisect.bisect_right(self.bkg_stats[(eid, ifos, param_group)], stat) \
            - \
            bisect.bisect_right(self.sngl_slide_stats[(eid, esid, ifos, param_group)], stat) \
            ) / self.bkg_durs[esid]

    def calc_cfar( self, esid, ifo_group, ufar ):
        """
        Calculates the combined far for the given uncombined far (ufar). This 
        is defined as the ufar times the number of categories that are active 
        at that point plus the sum of the max_bkg_fars of all the categories
        that are inactive. Whether or not a category is "active" is determined 
        by it's max_bkg_far. If the given ufar is greater than some max_bkg_far, 
        then the category which that max_bkg_far represents is considered 
        inactive. If the given ufar is less than some max_bkg_far, then 
        the category is considered active.
        """
        return \
            (len( self.max_bkg_fars[(esid, ifo_group)] ) - bisect.bisect_left( self.max_bkg_fars[(esid,ifo_group)], ufar ))*ufar \
            + sum([self.max_bkg_fars[(esid,ifo_group)][ii] for ii in range(bisect.bisect_left( self.max_bkg_fars[(esid,ifo_group)], ufar))])


# =============================================================================
#
#                          ExperimentSummary Utilities
#
# =============================================================================

# Following utilities are specific to the experiment_summary table

def update_experiment_summ_nevents( connection, verbose = False ):
    """
    Updates the number of events in the num_events column of the
    experiment_summary table. This should be used whenever coincs
    are deleted from the experiment_map table or when new files
    are added to a database.
    """
    if verbose:
        print >> sys.stderr, "Updating nevents column in experiment_summary table..."

    sqlquery = """
        UPDATE experiment_summary
        SET nevents = (
            SELECT COUNT(*)
            FROM experiment_map
            WHERE experiment_map.experiment_summ_id == experiment_summary.experiment_summ_id )
        """
    connection.cursor().execute(sqlquery)
    if verbose:
        print >> sys.stderr, "done."



# =============================================================================
#
#                             CoincInspiral Utilities
#
# =============================================================================

# Following utilities are specific to the coinc_inspiral table

def join_experiment_tables_to_coinc_inspiral():
    """
    Writes JOIN string to join the experiment, experiment_summary,
    and experiment_map table to the coinc_inspiral table. 
    NOTE: Should only use when querying the coinc_inspiral table (i.e.,
    the only table listed in the FROM statement is the coinc_inspiral).
    """

    return """ 
    JOIN
        experiment, experiment_summary, experiment_map 
    ON ( 
        experiment.experiment_id == experiment_summary.experiment_id
        AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id
        AND experiment_map.coinc_event_id == coinc_inspiral.coinc_event_id )"""


def apply_inclusion_rules_to_coinc_inspiral( connection, exclude_coincs = None, include_coincs = None, 
        param_filters = None, verbose = False ):
    """
    Clears coinc_inspiral table of coinc triggers falling outside of the
    desired ranges, as specified by parse_param_ranges and parse_coinc_opts.

    @connection: connection to a SQLite database with lsctables
    @param_filters: output of parse_param_ranges(...).get_param_filters()
    @include_coincs: output of parse_coinc_opts(...).get_coinc_filters().
       The coincs that are specified in this list will be SAVED.
    @exclude_coincs: output of parse_coinc_opts(...).get_coinc_filters().
        The coincs that are specified in this list will be DELETED.
    Note: exclude_coincs is applied first, so anything falling in it will 
    be deleted, regardless of wether or not the same falls in include_coincs.
    To avoid confusion, it is best to only specify one or the other, not both.
    """
    if verbose:
        print >> sys.stderr, "Removing coincs from coinc_inspiral table that " + \
            "fall outside of desired ranges and coinc-types..."

    join_conditions = join_experiment_tables_to_coinc_inspiral()

    if exclude_coincs:
        del_rows_from_table( connection, 'coinc_inspiral', 'coinc_event_id', 
            join_conditions,
            del_filters = exclude_coincs, verbose = verbose )
    if include_coincs:
        del_rows_from_table( connection, 'coinc_inspiral', 'coinc_event_id',
            join_conditions,
            save_filters = include_coincs, verbose = verbose )
    if param_filters:
        del_rows_from_table( connection, 'coinc_inspiral', 'coinc_event_id',
            join_conditions,
            save_filters = param_filters, verbose = verbose )



def clean_inspiral_tables( connection, verbose = False ):
    """
    Clears experiment_map, coinc_event, coinc_event_map, and all tables pointing to the
    coinc_event_map of triggers that are no longer in the coinc_inspiral table.
    Note that the experiment_summary, experiment, and time_slide_tables are left alone.
    This is because even if no events are left in an experiment, we still want info. about
    the experiment that was performed.
    """

    # Delete from experiment_map
    if verbose:
        print >> sys.stderr, '''Cleaning the experiment_map table...'''
    sqlquery = ' '.join([ 
            'DELETE',
            'FROM experiment_map',
            'WHERE coinc_event_id NOT IN (',
                'SELECT coinc_event_id',
                'FROM coinc_inspiral )' ])
    connection.cursor().execute( sqlquery )
    connection.commit()

    # Delete from coinc_event
    if verbose:
        print >> sys.stderr, '''Cleaning the coinc_event table...'''
    sqlquery = ' '.join([ 
            'DELETE',
            'FROM coinc_event',
            'WHERE coinc_event_id NOT IN (',
                'SELECT coinc_event_id',
                'FROM coinc_inspiral )' ])
    connection.cursor().execute( sqlquery )
    connection.commit()
  
    # Delete from coinc_definer
    if verbose:
        print >> sys.stderr, '''Cleaning the coinc_definer table...'''
    sqlquery = ' '.join([
            'DELETE',
            'FROM coinc_definer',
            'WHERE coinc_def_id NOT IN (',
                'SELECT coinc_def_id',
                'FROM coinc_event )' ])

    # Find tables listed in coinc_event_map
    sqlquery = 'SELECT DISTINCT table_name FROM coinc_event_map'
    table_names = connection.cursor().execute( sqlquery ).fetchall()

    # Delete from coinc_event_map
    if verbose:
        print >> sys.stderr, '''Cleaning the coinc_event_map table...'''
    sqlquery = ' '.join([
            'DELETE',
            'FROM coinc_event_map',
            'WHERE coinc_event_id NOT IN (',
                'SELECT coinc_event_id',
                'FROM coinc_event )' ])
    connection.cursor().execute( sqlquery )
    connection.commit()

    # Delete events from tables that were listed in the coinc_event_map
    for table in table_names:
        table = table[0]
        if verbose:
            print >> sys.stderr, '''Cleaning the %s table...''' % table
        sqlquery = ' '.join([
            'DELETE',
            'FROM', table,
            'WHERE event_id NOT IN (',
                'SELECT event_id',
                'FROM coinc_event_map )' ])
        connection.cursor().execute( sqlquery )
        connection.commit()

    if verbose:
        print >> sys.stderr, "done."



# =============================================================================
#
#                             TimeSlide Utilities
#
# =============================================================================

# Following utilities are specific to the time_slide table


def get_zero_lag_time_slide_ids( connection ):
    """
    Gets zero-lag time_slide_id's from the time_slide_table.
    """
    sqlquery = 'SELECT time_slide_id, offset FROM time_slide GROUP BY time_slide_id'
    slide_ids = connection.cursor().execute( sqlquery )
    zero_lag_ids = [slide_id[0] for slide_id in slide_ids if slide_id[1] == 0.]

    return zero_lag_ids


def get_zero_lag_instrument_sets( connection ):
    """
    Gets instrument sets from time slide table by using the ids of the zero-lag
    time-slides (Assumption is there is a zero-lag row in the time-slide table).
    """
    zero_lag_ids = get_zero_lag_time_slide_ids( connection )
  
    # sanity check
    if not zero_lag_ids:
        raise ValueError, "No zero-lag ids in time slide table, cannot get instrument set."
  
    zero_lag_instrument_sets = {}
    for id in zero_lag_ids:
        id = ''.join([ '"', id, '"' ])
        sqlquery = ' '.join(['SELECT instrument',
            'FROM time_slide',
            'WHERE time_slide_id ==', id ])
        instruments =  sorted(instrument[0] for instrument in connection.cursor().execute( sqlquery ).fetchall()) 
        if instruments not in zero_lag_instrument_sets:
            zero_lag_instrument_sets[ instruments ] = [ id ]
        else:
            zero_lag_instrument_sets[ instruments ].append( id )

    return zero_lag_instrument_sets


def get_instrument_sets_and_time_slide_ids( connection ):
    """
    Gets all instrument sets available in the time slide table and gets all
    time-slide ids associated with that instrument set. Since this only uses the
    time-slide table, will get everything even if there were no coincident
    events during a time-slide.
    """
    # get zero_lag ids and instrument set
    zero_lag_instrument_sets = get_zero_lag_instrument_sets( connection )
    # will save all ids to this new dictionary
    instrument_set_time_slide_ids = {}

    for instrument_set in zero_lag_instrument_sets:
        instrument_clause = ' AND '.join([ ''.join([ 
            'instrument == ', '"', instrument, '"' ]) for instrument in instrument_set ])
        sqlquery = ' '.join([ 
            'SELECT time_slide_id',
            'FROM time_slide',
            'WHERE', instrument_set ])
        instrument_set_time_slide_ids[ instrument_set ] = [ id[0] for id in
        connection.cursor().execute(sqlquery) ]

    return instrument_set_time_slide_ids

