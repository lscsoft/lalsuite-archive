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
import exceptions
from glue import gpstime
import mx.ODBC.DB2 as mxdb
from mx.ODBC.DB2 import SQL

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

    @return: Instance of class StateSegmentDatabaseException
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

    # connect to the database
    try:
      self.db = mxdb.Connect(dbname)
      self.db.setconnectoption(SQL.AUTOCOMMIT, SQL.AUTOCOMMIT_OFF)
      self.cursor = self.db.cursor()
    except Exception, e:
      msg = "Error connecting to database: %s" % e
      raise StateSegmentDatabaseException, e

    # generate a process table for this instance of the publisher
    try:
      cvs_date = time.strptime( '$Date$'[7:-2], '%Y/%m/%d %H:%M:%S' )

      sql = "VALUES GENERATE_UNIQUE()"
      self.cursor.execute(sql)
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
      msg = "Could not initialize process table: %s" % e
      raise StateSegmentDatabaseException, e

    # build a dictionary of the state vector types
    sql = "SELECT ifos, state_vec_major, state_vec_minor, segment_def_id "
    sql += "FROM segment_definer WHERE state_vec_major IS NOT NULL "
    sql += "AND state_vec_minor IS NOT NULL"
    try:
      self.cursor.execute(sql)
      result = self.cursor.fetchall()
    except Exception, e:
      msg = "Error fetching state vector values from database: %s" % e
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
      sql = "UPDATE process SET (end_time) = (?) WHERE process_id = '%s'" % self.process_id
      self.cursor.execute(sql,tuple([now]))
      self.db.commit()
    except Exception, e:
      msg = "Error inserting end_time into database: %s" % e
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
    
    try:
      # get a unique id for this lfn
      sql = "VALUES GENERATE_UNIQUE()"
      self.cursor.execute(sql)
      self.lfn_id = self.cursor.fetchone()[0]
     
      sql = "INSERT INTO lfn (process_id,lfn_id,lfn,start_time,end_time) "
      sql += "values (?,?,?,?,?)"
      self.cursor.execute(sql,(self.process_id,self.lfn_id,lfn,start,end))
      self.db.commit()

    except mxdb.InterfaceError, e:
      if e[1] == -803:
        sql = "SELECT lfn_id from lfn WHERE lfn = '%s'" % lfn
        self.cursor.execute(sql)
        self.lfn_id = self.cursor.fetchone()[0]
      else:
        msg = "Unable to obtain unique id for LFN : %s : %s" % (lfn,e)
        raise StateSegmentDatabaseException, msg

    except Exception, e:
      msg = "Unable to create entry for LFN : %s : %s" % (lfn,e)
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
    if (ver, val) not in self.state_vec:
      try:
        sql = "VALUES GENERATE_UNIQUE()"
        self.cursor.execute(sql)
        self.state_vec[ifo][(ver,val)] = self.cursor.fetchone()[0]

        sql = "INSERT INTO segment_definer (process_id,segment_def_id,ifos,"
        sql += "name,version,comment,state_vec_major,state_vec_minor) VALUES "
        sql += "(?,?,?,?,?,?,?,?)"

        self.cursor.execute(sql,(self.process_id, self.state_vec[(ver,val)], ifo,
          'STATEVEC.%d.%d' % (ver, val), 0, 
          'Created automatically by StateSegmentDatabase', ver, val))
      except:
        self.db.rollback()
        msg = "Error inserting new state vector type into database : %s" % e
        raise StateSegmentDatabaseException, e
      if self.debug:
        print ("DEBUG: created a new state vec type (%d,%d), id = " % \
          (ver,val)), 
        print tuple([self.state_vec[(ver,val)]])

    sv_id = self.state_vec[(ver,val)]

    # insert the state segment 
    sql = "VALUES GENERATE_UNIQUE()"
    self.cursor.execute(sql)
    segment_id = self.cursor.fetchone()[0]

    sql = "INSERT INTO segment (process_id, segment_id,"
    sql += "start_time,end_time,active) VALUES (?,?,?,?,?)"

    try:
      self.cursor.execute(sql,
        (self.process_id,segment_id,start_time,end_time,1))
    except Exception, e:
      self.db.rollback()
      msg = "error inserting segment information : %s" % e
      raise StateSegmentDatabaseException, msg

    sql = "INSERT INTO segment_lfn_map (process_id,segment_id,lfn_id) "
    sql += "VALUES (?,?,?)"
    try:
      self.cursor.execute(sql,(self.process_id, segment_id, self.lfn_id))
    except Exception, e:
      self.db.rollback()
      msg = "error inserting segment information : %s" % e
      raise StateSegmentDatabaseException, msg

    sql = "INSERT INTO segment_def_map (process_id,segment_id,segment_def_id) "
    sql += "VALUES (?,?,?)"
    try:
      self.cursor.execute(sql,(self.process_id, segment_id, sv_id))
    except Exception, e:
      self.db.rollback()
      msg = "error inserting segment information : %s" % e
      raise StateSegmentDatabaseException, msg

    # all the transactions went through ok, so commit everything
    self.db.commit()

    if self.debug:
      print "DEBUG: inserted with segment_id "
      print tuple([segment_id])

