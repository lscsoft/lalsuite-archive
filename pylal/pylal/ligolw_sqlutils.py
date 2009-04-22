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

from glue.ligolw import dbtables

__author__ = "Collin Capano <cdcapano@physics.syr.edu>"
__date__ = "$Date$" 
__version__ = "$Revision$"



# =============================================================================
#
#                             Generic coinc_table Utilities
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
      @param_ranges_opt: string from the --param-ranges option.
      """
      if verbose:
        print >> sys.stderr, "Parsing param-ranges..."

      # make param unique by appending 'coinc_inspiral' to the param name
      self.param = '.'.join([ table_name, table_param ])

      ranges = param_ranges_opt.split(';')
      self.param_ranges = []

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
          raise ValueError, "Parameter range %s not formatted correctly!" % this_range
  
        # get upper-bound (similar to lower bound method)
        upperparam = this_range.split(',')[1].strip()
        if upperparam.find(']') != -1:
          upperbndry = '<='
          upperparam = float( upperparam.rstrip(']') )
        elif upperparam.find(')') != -1:
          upperbndry = '<'
          upperparam = float( upperparam.rstrip(')') )
        else:
          raise ValueError, "Parameter range %s not formatted correctly!" % this_range

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
            self.param_filters.append( ' '.join([ self.param, lowerbndry, lowerparam, 
              'AND', self.param, upperbndry, upperparam ]) )

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


def parse_exclude_coincs( exclude_coincs_opt, verbose = False ):
  """
  Parse --exclude-coincs option. Returns a list of strings that can be used in a SQLite
  WHERE clause to filter out coinc types.
  """

  if verbose:
    print >> sys.stderr, "Parsing exclude-coincs..."

  exclude_coinc_filters = []

  for rule in exclude_coincs_opt.split(';'):
    exclude_filter = ''
    rule = rule.strip().lstrip('[').rstrip(']').upper()

    # get instruments_on, coinc_types 
    [ coinc_types, instruments_on ] = rule.split(' IN ')
    
    if instruments_on.strip() != 'ALL':
      instruments_on = instruments_on.split(',')
      # sanity check
      if len(instruments_on) <= 1:
        raise ValueError, "Must delimit instruments in --exclude-coincs by commas."
      instruments_on = ','.join( sorted( instrument.strip() for instrument in instruments_on ))
      instruments_on = ''.join([ '"', instruments_on, '"' ])
      exclude_filter = ' '.join([ 'experiment.instruments ==', instruments_on ])
    
    # if coinc_types are ALL, means to filter out all coincs in specified
    # instruments_on: just go on to next rule. Otherwise, evaluate coinc type
    coinc_filters = ''
    if coinc_types.strip() != 'ALL':
      for coinc_type in coinc_types.split('+'):
        coinc_type = ','.join( sorted( instrument.strip() for instrument in coinc_type.split(',') ))
        coinc_type = ''.join([ '"', coinc_type, '"' ])
        if not coinc_filters:
          coinc_filters = ' '.join([ 'coinc_inspiral.ifos ==', coinc_type ])
        else:
          coinc_filters = ' '.join([ coinc_filters, 'OR', 'coinc_inspiral.ifos ==', coinc_type ])
      # join coinc_filters to exclude_filter
      if not exclude_filter:
        exclude_filter = coinc_filters
      else:
        exclude_filter = ' '.join([ exclude_filter, 'AND', '(', coinc_filters, ')' ])
    
    # add filter to coinc_filters
    exclude_coinc_filters.append( exclude_filter )
  
  if verbose:
    print >> sys.stderr, "done."

  return exclude_coinc_filters


def del_rows_from_table( connection, del_table, del_table_id, join_conditions, del_filters = None, save_filters = None, verbose = False ):
  """
  Deletes triggers from any specified table in the del_table option.
  @connection: DBTables connection to a database
  @del_table: Any coinc_table (coinc_inspiral, sngl_inspiral, coinc_event,
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

  join_conditions = ' '.join([ 
        'JOIN experiment, experiment_summary, experiment_map',
        'ON experiment.experiment_id == experiment_summary.experiment_id',
          'AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id',
          'AND experiment_map.coinc_event_id == coinc_inspiral.coinc_event_id' ])

  return join_conditions



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
    instruments =  sorted(instrument[0]
      for instrument in connection.cursor().execute( sqlquery ).fetchall()) 
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

