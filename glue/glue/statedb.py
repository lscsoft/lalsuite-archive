"""
The statedb module is the interface between the segment publishing script
and the metadata database.

$Id$
"""
import os
import socket
import pwd
import sys
import re
import time
import random
import exceptions
from glue import gpstime
import mx.ODBC.DB2



class StateSegmentDatabaseException(exceptions.Exception):
  """
  Exceptions raised by the classes and methods in this module
  will be instances of this class.
  """
  def __init__(self, args=None):
    """
    Create an instance of this class, ie. a StateSegmentDatabaseException
    exception.

    @param args: 

    @return: Instance of class StateSegmentDatabaseException
    """
    self.args = args



class StateSegmentDatabaseSegmentExistsException(exceptions.Exception):
  """
  Exceptions raised by the classes and methods in this module
  will be instances of this class.
  """
  def __init__(self, args=None):
    """
    Create an instance of this class, ie. a StateSegmentDatabaseException
    exception.

    @param args: 

    @return: Instance of class StateSegmentDatabaseSegmentExistsException
    """
    self.args = args



class StateSegmentDatabaseLFNExistsException(exceptions.Exception):
  """
  Exceptions raised by the classes and methods in this module
  will be instances of this class.
  """
  def __init__(self, args=None):
    """
    Create an instance of this class exception.

    @param args: 

    @return: Instance of class StateSegmentDatabaseLFNExistsException
    """
    self.args = args



class StateSegmentDatabase:
  """
  Class that represents an instance of a state segment database
  """
  def __init__(self, dbname, dbuser = '', dbpasswd = '', debug = False):
    """
    Open a connection to the state segment database.

    @param dbname: the name of the database containing the segments
    @type dbname: string
  
    @param username: the username which has permission to write to 
      the state segment database
    @type username: string

    @param passwd: the password of the user who has permission to write
      to the state segment database (optional)
    @type username: string
    """
    self.debug = debug
    self.db = None
    self.cursor = None
    self.state_vec = {}
    self.state_vec['H1'] = {}
    self.state_vec['H2'] = {}
    self.state_vec['L1'] = {}
    self.lfn_id = None
    self.framereg = re.compile(r'^([A-Za-z]+)\-(\w+)\-(\d+)\-(\d+)\.gwf$')
    self.process_id = None

    # seed the random number generator with the current time
    random.seed()

    # connect to the database
    try:
      self.db = mx.ODBC.DB2.Connect( dbname )
      self.db.setconnectoption( mx.ODBC.DB2.SQL.ATTR_AUTOCOMMIT, 
        mx.ODBC.DB2.SQL.AUTOCOMMIT_OFF )
      self.db.setconnectoption( mx.ODBC.DB2.SQL.ATTR_TXN_ISOLATION, 
        mx.ODBC.DB2.SQL.TXN_READ_UNCOMMITTED )
      self.cursor = self.db.cursor()
    except Exception, e:
      msg = "Error connecting to database: %s" % e
      raise StateSegmentDatabaseException, e

    # generate a process table for this instance of the publisher
    try:
      cvs_date = time.strptime( 
        '$Date$'[7:-2], '%Y/%m/%d %H:%M:%S' )

      sql = "VALUES GENERATE_UNIQUE()"
      self.cursor.execute(sql)
      self.db.commit()
      self.process_id = self.cursor.fetchone()[0]

      process_tuple = ( os.path.basename(sys.argv[0]), 
        '$Revision$'[11:-2], 
        '$Source$' [9:-2],
        gpstime.GpsSecondsFromPyUTC( time.mktime(cvs_date) ),
        1,
        socket.gethostname(),
        pwd.getpwuid(os.geteuid())[0],
        os.getpid(),
        gpstime.GpsSecondsFromPyUTC(time.time()),
        self.process_id )
    
      sql = ' '.join(["INSERT INTO process (program,version,cvs_repository,",
      "cvs_entry_time,is_online,node,username,unix_procid,start_time,",
      "process_id) VALUES (", ','.join(['?' for x in process_tuple]), ")" ])

      self.cursor.execute(sql,process_tuple)
      self.db.commit()

    except Exception, e:
      msg = "Unable to initialize process table: %s" % e
      raise StateSegmentDatabaseException, e

    # build a dictionary of the state vector types
    sql = "SELECT ifos, state_vec_major, state_vec_minor, segment_def_id "
    sql += "FROM segment_definer WHERE state_vec_major IS NOT NULL "
    sql += "AND state_vec_minor IS NOT NULL"

    try:
      self.cursor.execute(sql)
      result = self.cursor.fetchall()
    except Exception, e:
      msg = "Unable to get existing segment_def_id from database: %s" % e
      raise StateSegmentDatabaseException, e

    for r in result:
      self.state_vec[r[0].strip()][tuple([r[1],r[2]])] = r[3]

    if self.debug:
      print "DEBUG: current state known state vec types are: "
      print self.state_vec
      

  def __del__(self):
    try:
      self.close()
    except:
      pass
    try:
      del self.cursor
    except:
      pass
    try:
      del self.db
    except:
      pass


  def close(self):
    """
    Close the connection to the database.
    """
    try:
      now = gpstime.GpsSecondsFromPyUTC(time.time())
      sql = "UPDATE process SET (end_time) = (?) WHERE process_id = (?)" 
      self.cursor.execute(sql,(now, self.process_id))
      self.db.commit()
    except Exception, e:
      msg = "Unable to update end_time in database: %s" % e
      raise StateSegmentDatabaseException, msg

    try:
      self.cursor.close()
      self.db.close()
    except Exception, e:
      msg = "Error closing connection to database: %s" % e
      raise StateSegmentDatabaseException, msg
      

  def register_lfn(self,lfn,start=None,end=None):
    """
    Start publishing state information for a new logical file name

    @param lfn: logical file name of object containing state data
    @type lfn: string
    """
    self.lfn_id = None

    r = self.framereg.search(lfn)
    if not r:
      if not start or not end:
        msg = "File %s not a frame file. Start and end times must be given" % \
          lfn
        raise StateSegmentDatabaseException, msg
    else:
      start = long(self.framereg.search(lfn).group(3))
      end = long(self.framereg.search(lfn).group(3)) + \
        long(self.framereg.search(lfn).group(4))
    
    # try and get an lfn_id for this file
    sql = "SELECT lfn_id from lfn WHERE lfn = '%s'" % lfn
    try:
      self.cursor.execute(sql)
      result = self.cursor.fetchone()
    except Exception, e:
      msg = "Unable to query database for lfn_id (%s): %s" % (lfn,e)
      raise StateSegmentDatabaseException, msg

    if result:
      # use the lfn_id returned by the database and raise an exception
      self.lfn_id = result[0]
      raise StateSegmentDatabaseLFNExistsException

    else:
      # generate an lfn_id for this file
      sql = "VALUES GENERATE_UNIQUE()"
      try:
        self.cursor.execute(sql)
        self.lfn_id = self.cursor.fetchone()[0]
      except Exception, e:
        msg = "Unable to generate a new lfn_id (%s): %s" % (lfn,e)
        raise StateSegmentDatabaseException, msg

      # create the entry in the lfn table for this file
      sql = "INSERT INTO lfn (process_id,lfn_id,lfn,start_time,end_time) "
      sql += "values (?,?,?,?,?)"
      try:
        self.cursor.execute(sql,(self.process_id,self.lfn_id,lfn,start,end))
        self.db.commit()

      except (mx.ODBC.DB2.InterfaceError, mx.ODBC.DB2.InternalError), e:
        self.db.rollback()
        if ( (e[1] == -803 and int(e[0]) == 23505) or 
             (e[1] == -911 and int(e[0]) == 40001) ):
          # someone may have just beaten us to inserting this file name
          # try another select to get the lfn_id and abort if it fails
          e_prev = e
          sql = "SELECT lfn_id from lfn WHERE lfn = '%s'" % lfn
          try:
            self.cursor.execute(sql)
            self.lfn_id = self.cursor.fetchone()[0]
          except Exception, e:
            msg = "Unable to resolve race on lfn_id (%s) caused by %s: %s" \
              % (lfn,e_prev,e)
            raise StateSegmentDatabaseException, msg
        else:
          # some other database error we cannot handle: just give up
          msg = "Unable to insert entry into lfn table (%s): %s" % (lfn,e)
          raise StateSegmentDatabaseException, msg

      except Exception, e:
        # give up on all other errors
        self.db.rollback()
        msg = "Unable to insert entry into lfn table (%s): %s" % (lfn,e)
        raise StateSegmentDatabaseException, msg

     
  def publish_state(self, 
      ifo, start_time, start_time_ns, end_time, end_time_ns, ver, val ):
    """
    Publish a state segment for a state vector in the database
    """

    # check that we have an lfn registered
    if not self.lfn_id:
      msg = "No LFN registered to publish state information"
      raise StateSegmentDatabaseException, msg

    # see if we need a new state val or if we know it already
    if (ver, val) not in self.state_vec[ifo]:

      # generate a unique segment_def_id for this segment_definer
      sql = "VALUES GENERATE_UNIQUE()"
      try:
        self.cursor.execute(sql)
        self.state_vec[ifo][(ver,val)] = self.cursor.fetchone()[0]
      except Exception, e:
        msg = "Could not generate a unique segment_def_id: %s" % e
        raise StateSegmentDatabaseException, msg

      # insert a segment_definer row for this state vector value
      sql = "INSERT INTO segment_definer (process_id,segment_def_id,ifos,"
      sql += "name,version,comment,state_vec_major,state_vec_minor) VALUES "
      sql += "(?,?,?,?,?,?,?,?)"

      try:
        self.cursor.execute(sql,(self.process_id, self.state_vec[ifo][(ver,val)], 
          ifo, 'STATEVEC.%d.%d' % (ver, val), 0, 
          'Created automatically by StateSegmentDatabase', ver, val))
        self.db.commit()

      except (mx.ODBC.DB2.InterfaceError, mx.ODBC.DB2.InternalError), e:
        self.db.rollback()
        if ( (e[1] == -803 and int(e[0]) == 23505) or 
             (e[1] == -911 and int(e[0]) == 40001) ):
          # someone may have just beaten us to inserting this segment_definer
          # try another select to get the segment_def_id and abort if it fails
          e_prev = e
          sql = "SELECT segment_def_id FROM segment_definer WHERE "
          sql += "ifos = '" + ifo + "' AND "
          sql += "state_vec_major = %d AND state_vec_minor = %d" % (ver,val)
          try:
            self.cursor.execute(sql)
            self.state_vec[ifo][(ver,val)] = self.cursor.fetchone()[0]
          except Exception, e:
            msg = "Unable to resolve race on segment_def_id caused by %s: %s" \
              % (e_prev,e)
            raise StateSegmentDatabaseException, e
        else:
          # some other database error we cannot handle: just give up
          msg = "Unable to insert entry into segment_definer table: %s" % e
          raise StateSegmentDatabaseException, msg

      except Exception, e:
        # give up on all other errors
        msg = "Unable to insert entry into segment_definer table: %s" % e
        raise StateSegmentDatabaseException, e

      if self.debug:
        print ("DEBUG: using new state vec type (%s,%d,%d), id = " % \
          (ifo,ver,val)), 
        print tuple([self.state_vec[ifo][(ver,val)]])

    # save the state segment_def_id for use below
    sv_id = self.state_vec[ifo][(ver,val)]

    # generate a unique id for the segment we are going to insert
    sql = "VALUES GENERATE_UNIQUE()"
    try:
      self.cursor.execute(sql)
      segment_id = self.cursor.fetchone()[0]
    except Exception, e:
      msg = "Unable to generate unique segment_id for segment: %s" % e
      raise StateSegmentDatabaseException, e

    # insert the segment
    sql = "INSERT INTO state_segment (process_id,segment_id,segment_def_id,"
    sql += "start_time,start_time_ns,end_time,end_time_ns,lfn_id)"
    sql += "VALUES (?,?,?,?,?,?,?,?)"

    for attempt in range(3):
      try:
        self.cursor.execute(sql, (self.process_id,segment_id,sv_id,
          start_time,start_time_ns,end_time,end_time_ns,self.lfn_id))
        self.db.commit()
        break

      except mx.ODBC.DB2.InternalError, e:
        self.db.rollback()
        if e[1] == -911 and int(e[0]) == 40001:
          # retry up to three times on a deadlock or timeout
          if attempt < 3:
            time.sleep(random.randrange(0,5,1))
          else:
            msg = "Unable to insert segment after 3 deadlocks: %s" % e
            raise StateSegmentDatabaseException, msg
        else:
          # give up on all other internal errors
          msg = "Unable to insert segment: %s" % e
          raise StateSegmentDatabaseException, msg

      except mx.ODBC.DB2.InterfaceError, e:
        self.db.rollback()
        if e[1] == -438 and int(e[0]) == 70001:
          # if this is a duplicate entry, throw the exists exception
          raise StateSegmentDatabaseSegmentExistsException
        else:
          # give up on all other interface errors
          msg = "Unable to insert segment: %s" % e
          raise StateSegmentDatabaseException, msg

      except Exception, e:
        # catch all other exceptions and give up
        msg = "Unable to insert segment: %s" % e
        raise StateSegmentDatabaseException, msg

    if self.debug:
      print "DEBUG: inserted with segment_id "
      print tuple([segment_id])
