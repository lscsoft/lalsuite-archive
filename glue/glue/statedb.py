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
  def __init__(self, dbname, dbuser, dbpasswd = ''):
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
    self.debug = True
    self.db = None
    self.cursor = None
    self.state_vec = {}
    self.lfn_id = None
    self.framereg = re.compile(r'^([A-Za-z]+)\-(\w+)\-(\d+)\-(\d+)\.gwf$')
    self.process_id = None
    self.version = '$Revision$'[11:-2]
    self.cvs_repository = '$Source$' [9,-2]
    self.cvs_date = time.strptime( '$Date$'[7,-2], '%Y/%m/%d %T' )
    self.cvs_entry_time = str( gpstime.GpsSecondsFromPyUTC( self.cvs_date ) )
    self.node = socket.gethostname()
    self.username = pwd.getpwuid(os.geteuid())[0]
    self.unix_procid = os.getpid()
    self.start_time = str( gpstime.GpsSecondsFromPyUTC(time.time()) )

    # connect to the database
    try:
      self.db = mxdb.Connect(dbname)
      self.cursor = self.db.cursor()
    except Exception, e:
      msg = "Error connecting to database: %s" % e
      raise StateSegmentDatabaseException, e

    # generate a process table for this instance of the publisher
    try:
      sql = "VALUES GENERATE_UNIQUE()"
      self.cursor.execute(sql)
      self.process_id = self.cursor.fetchone()
      sql = "INSERT INTO process (program,version,cvs_repository,"
      sql += "cvs_entry_time,is_online,node,username,unix_procid,start_time,"
      sql += "process_id) VALUES ('StateSegmentDatabase', '" + self.version + 
      sql += "', '" + self.cvs_repository + "', '" + self.cvs_entry_time
      sql += "', 1, '" + self.node + "', '" + self.username + "', " 
      sql +=  self.unix_procid + ", " + self.start_time + ", " 
      sql += self.process_id + ")"

    

    # build a dictionary of the state vector types
    sql = "SELECT version, value, state_vec_id from state_vec"
    try:
      self.cursor.execute(sql)
      result = self.cursor.fetchall()
    except Exception, e:
      msg = "Error fetching state vector values from database: %s" % e
      raise StateSegmentDatabaseException, e
    for r in result:
      if r[0] and r[1]:
        self.state_vec[tuple([r[0],r[1]])] = r[2]

    if self.debug:
      print "DEBUG: current state known state vec types are: "
      print self.state_vec
     
      
  def close(self):
    """
    Close the connection to the database.
    """
    try:
      self.cursor.close()
      self.db.close()
    except Exception, e:
      msg = "Error disconnecting from database: %s" % e
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
      
      
      sql = "INSERT INTO lfn (lfn,start_time,end_time) values (%s,%s,%s)"
      self.cursor.execute(sql,(lfn,start,end))
      sql = "SELECT LAST_INSERT_ID()"
      self.cursor.execute(sql)
    except _mysql_exceptions.IntegrityError, e:
      sql = "SELECT lfn_id from lfn WHERE lfn = '%s'" % lfn
      self.cursor.execute(sql)
    except Exception, e:
      msg = "Unable to obtain unique id for LFN : %s : %s" % (lfn,e)
      raise StateSegmentDatabaseException, msg
      
    self.lfn_id = self.cursor.fetchall()[0][0]


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
        sql = "INSERT INTO state_vec (version,value) VALUES (%s,%s)"
        self.cursor.execute(sql,(ver,val))
        sql = "SELECT LAST_INSERT_ID()"
        self.cursor.execute(sql)
        self.state_vec[(ver,val)] = self.cursor.fetchall()[0][0]
      except:
        msg = "Error inserting new state vector type into database : %s" % e
        raise StateSegmentDatabaseException, e
      if self.debug:
        print "DEBUG: create a new state vec type (%d,%d) -> %d" % \
          (ver,val,self.state_vec[(ver,val)])

    sv_id = self.state_vec[(ver,val)]

    # insert the state segment 
    sql = "INSERT INTO state_segment (ifo,start_time,start_time_ns,"
    sql += "end_time,end_time_ns,state_vec_id,lfn_id) VALUES "
    sql += "(%s,%s,%s,%s,%s,%s,%s)"

    try:
      self.cursor.execute(sql,(ifo,start_time,start_time_ns,
            end_time,end_time_ns,sv_id,self.lfn_id))
    except _mysql_exceptions.IntegrityError, e:
      msg = "error : this state segment already exists in the database"
      raise StateSegmentDatabaseSegmentExistsException, msg
    except Exception, e:
      msg = "error inserting segment information : %s" % e
      raise StateSegmentDatabaseException, msg

    if self.debug:
      sql = "SELECT LAST_INSERT_ID()"
      self.cursor.execute(sql)
      id = self.cursor.fetchall()[0][0]
      print "DEBUG: inserted with id", id


  def set_state_name(self,ver,val,name):
    try:
      sql = "INSERT INTO state_vec (version,value,state) VALUES (%s,%s,%s)"
      self.cursor.execute(sql,(ver,val,name))
    except _mysql_exceptions.IntegrityError, e:
      sql = "UPDATE state_vec SET state = '%s'" % name
      sql += "WHERE version = %s AND value = %s" % (ver,val)
      self.cursor.execute(sql)
    except Exception, e:
      msg = "Unable to state state name for state %s,%s : %s" % (ver,val,e)
      raise StateSegmentDatabaseException, msg


  def get_state_name(self,ver,val):
    """
    Get the name of a state by version and value
    """
    sql = "SELECT state from state_vec WHERE version = %s AND value = %s" % \
      (ver,val)
    try:
      self.cursor.execute(sql)
      r = str(self.cursor.fetchall()[0][0])
    except Excpetion, e:
      msg = "error getting state name : %s" % e
      raise StateSegmentDatabaseException, msg
