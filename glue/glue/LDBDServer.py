"""
The LDBDServer module provides an API for responding to request from the
LDBDClient by connecting to the DB2 database.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.

$Id: LDBDServer.py,v 1.54 2009/02/10 16:34:05 duncan Exp $

This file is part of the Grid LSC User Environment (GLUE)

GLUE is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__version__ = '$Revision: 1.54 $'[11:-2]

import os
import sys
import re
import types
import pyRXP
import exceptions
import socket
import SocketServer
import cPickle
from glue import ldbd
try:
  import rlsClient
except:
  pass


def dtd_uri_callback(uri):
  if uri == 'http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt':
    # if the XML file contants a http pointer to the ligolw DTD at CIT then
    # return a local copy to avoid any network problems
    return 'file://localhost' + os.path.join( os.environ["GLUE_LOCATION"],
      'etc/ligolw_dtd.txt' )
  else:
    # otherwise just use the uri in the file
    return uri

def initialize(configuration,log):
  # define the global variables used by the server
  global logger, max_bytes, xmlparser, dbobj, rls
  global dmt_proc_dict, dmt_seg_def_dict, creator_db
  global ldbd_com, db2_com, port

  # initialize the logger
  logger = log
  logger.info("Initializing server module %s" % __name__ )
  
  # initialize the database hash table
  dbobj = ldbd.LIGOMetadataDatabase(configuration['dbname'])
  max_bytes = configuration['max_client_byte_string']

  # initialize the ldbd commands and db2 command restrictions
  port = configuration['port']
  ldbd_com = configuration['ldbd_com'].split(',')
  db2_com = configuration['db2_com'].split(',')

  # create the xml parser
  xmlparser = pyRXP.Parser()

  # use a local copy of the DTD, if one is available
  try:
    GLUE_LOCATION = os.environ["GLUE_LOCATION"]
    xmlparser.eoCB = dtd_uri_callback
    logger.info("Using local DTD in " + 
      'file://localhost' + os.path.join( GLUE_LOCATION, 'etc') )
  except KeyError:
    logger.warning('GLUE_LOCATION not set, unable to use local DTD') 

  # open a connection to the rls server
  rls_server = configuration['rls']
  cert = configuration['certfile']
  key = configuration['keyfile']
  try:
    rls = rlsClient.RlsClient(rls_server,cert,key)
  except:
    rls = None

  # initialize dictionaries for the dmt processes and segments definers
  dmt_proc_dict = {}
  dmt_seg_def_dict = {}
  creator_db = None

def shutdown():
  global logger, max_bytes, xmlparser, dbobj, rls
  global dmt_proc_dict, dmt_seg_def_dict
  logger.info("Shutting down server module %s" % __name__ )
  if rls:
    del rls
  del xmlparser
  del dbobj
  del dmt_proc_dict
  del dmt_seg_def_dict

class ServerHandlerException(exceptions.Exception):
  """
  Class representing exceptions within the ServerHandler class.
  """
  def __init__(self, args=None):
    """
    Initialize an instance.

    @param args: 

    @return: Instance of class ServerHandlerException
    """
    self.args = args
        
class ServerHandler(SocketServer.BaseRequestHandler):
  """
  An instance of this class is created to service each request of the server.
  """
  def handle(self):
    """
    This method does all the work of servicing a request to the server. See
    the documentation for the standard module SocketServer.

    The input from the socket is parsed for the method with the remaining
    strings stripped of null bytes passed to the method in a list.

    There are no parameters. When the instance of the class is created to
    process the request all necessary information is made attributes of the
    class instance.

    @return: None
    """
    global logger
    global max_bytes

    logger.debug("handle method of %s class called" % __name__)

    # mapping of ldbdd RPC protocol names to methods of this class
    methodDict = {
      'PING' : self.ping,
      'QUERY' : self.query,
      'INSERT' : self.insert,
      'INSERTMAP' : self.insertmap,
      'INSERTDMT' : self.insertdmt
    }

    if True:
      # from the socket object create a file object
      self.sfile = self.request.makefile("rw")
      f = self.sfile

      # read all of the input up to limited number of bytes
      input = f.read(size=max_bytes,waitForBytes=2)

      # try 10 more times if we don't have a null byte at the end
      while input[-1] != '\0':
        input += f.read(size=max_bytes,waitForBytes=2)

      # the format should be a method string, followed by a null byte
      # followed by the arguments to the method encoded as null
      # terminated strings

      # check if the last byte is a null byte
      if input[-1] != '\0':
        logger.error("Bad input on socket: %s" % input)
        raise ServerHandlerException, "Last byte of input is not null byte"
    #except Exception, e:
    #  logger.error("Error reading input on socket: %s" %  e)
    #  return

    try:
      # parse out the method and arguments 
      stringList = input.split('\0')
      methodString = stringList[0]
      argStringList = stringList[1:-1]
                        
    except Exception, e:
      logger.error("Error parsing method and argument string: %s" % e)

      msg = "ERROR LDBDServer Error: " + \
        "Error parsing method and argument string: %s" % e
      self.__reply__(1, msg)
      return

    # Set read and read/write access to different ports according to their configuration file
    try:
      # ldbd_com lists the allowed methodString for ldbd server based on the port number it is running on
      if methodString in ldbd_com:
        pass
      else:
        msg = "\nOnly authorized users can %s\n" % methodString
        msg += "To %s, authorized users please explicitly append port number 30020 to " % methodString
        msg += "your --server agument at command line" 
        logmsg = "ldbd server on port %d DO NOT support %s" % (port, methodString)
        logger.error(logmsg) 
        self.__reply__(1,msg)
        raise ServerHandlerException, msg 
      # list allowed sql commands when methodString is QUERY
      if methodString=='QUERY':
        # get the sql command, for example "SELECT, UPDATE, INSERT"
        # see if the command is in the .ini file "db2_com" variable 
        dbcommand =  argStringList[0].split(' ')[0].upper()
        if dbcommand in db2_com:
          pass
        else:
           msg = 'ldbd server on port %d DO NOT support "%s"' % (port, dbcommand)
           logger.error(msg) 
           self.__reply__(1,msg)
           raise ServerHandlerException, msg
    except Exception, e:
      logger.error("Error filtering allowed commands: %s" % e) 
      return

                
    try:
      # look up method in dictionary
      method = methodDict[methodString]
    except Exception, e:
      msg = "Error converting method string %s to method call: %s" % \
        (methodString, e)
      logger.error(msg)
                        
      self.__reply__(1, msg)
      return

    try:
      # call the method requested with the rest of strings as input
      result = method(argStringList) 
      self.__reply__( result[0], result[1] )
    except Exception, e:
      logger.error("Error while calling method %s: %s" % (methodString, e))

    return
        
  def __reply__(self, code, msg):
    """
    Format and send a reply back down the socket to the client. The file
    representing the socket is closed at the end of this method.

    @param code: integer representing the error code with 0 for success
                
    @param msg: object to be passed back to the client, either a string
    or a list of items that can be represented by strings
                        
    @return: None
    """
    f = self.sfile
    reply = "%d\0%s\0" % (code, msg)
    f.write(reply)

    # close the file associated with the socket
    f.close()

  def ping(self, arg):
    """
    Bounce back alive statment. Corresponds to the PING method in the
    ldbdd RPC protocol.

    @param arg: list (perhaps empty) of strings representing message sent
      by client to server

    @return: None
    """
    global logger

    logger.debug("Method ping called")
    try:
      hostname = socket.getfqdn()
      msg = "%s at %s is alive" % (__name__, hostname)
    except Exception, e:
      msg = "%s is alive" % __name__

    return (0, msg)

  def query(self, arg):
    """
    Execute an SQL query on the database and return the result as LIGO_LW XML

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """
    global logger
    global xmlparser, dbobj

    # get the query string and log it
    querystr = arg[0]
    logger.debug("Method query called with %s" % querystr)

    # assume failure
    code = 1

    try:
      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # execute the query
      rowcount = ligomd.select(querystr)

      # convert the result to xml
      result = ligomd.xml()

      logger.debug("Method query: %d rows returned" % rowcount)
      code = 0
    except Exception, e:
      result = ("Error querying metadata database: %s" % e)
      logger.error(result)

    try:
      del ligomd
      del lwtparser
    except Exception, e:
      logger.error(
        "Error deleting metadata object in method query: %s" % e)

    return (code,result)

  def insert(self, arg):
    """
    Insert some LIGO_LW xml data in the metadata database

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """
    global logger
    global xmlparser, dbobj

    logger.debug("Method insert called")

    # assume failure
    code = 1

    try:
      # capture the remote users DN for insertion into the database
      cred = self.request.get_delegated_credential()
      remote_dn = cred.inquire_cred()[1].display()

      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # parse the input string into a metadata object
      ligomd.parse(arg[0])

      # add a gridcert table to this request containing the users dn
      ligomd.set_dn(remote_dn)

      # insert the metadata into the database
      result = str(ligomd.insert())

      logger.info("Method insert: %s rows affected by insert" % result)
      code = 0
    except Exception, e:
      result = ("Error inserting metadata into database: %s" % e)
      logger.error(result)

    try:
      del ligomd
      del lwtparser
    except Exception, e:
      logger.error(
        "Error deleting metadata object in method insert: %s" % e)

    return (code,result)

  def insertmap(self, arg):
    """
    Insert some LIGO_LW xml data in the metadata database with an LFN to
    PFN mapping inserted into the RLS database.

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """
    msg = "server is not initialized for RLS connections"
    logger.error(msg)
    return (1, msg)

  def insertdmt(self, arg):
    """
    Insert LIGO_LW xml data from the DMT in the metadata database. For
    DMT inserts, we need to check for existing process_id and
    segment_definer_id rows and change the contents of the table to be
    inserted accordingly. We must also update the end_time of any 
    existing entries in the process table.

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """
    global logger
    global xmlparser, dbobj
    global dmt_proc_dict, dmt_seg_def_dict, creator_db
    proc_key = {}
    known_proc = {}
    seg_def_key = {}

    logger.debug( "Method dmtinsert called." )
    logger.debug( "Known processes %s, " % str(dmt_proc_dict) )
    logger.debug( "Known segment_definers %s" % str(dmt_seg_def_dict) )

    # assume failure
    code = 1

    try:
      # capture the remote users DN for insertion into the database
      cred = self.request.get_delegated_credential()
      remote_dn = cred.inquire_cred()[1].display().strip()

      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # parse the input string into a metadata object
      logger.debug("parsing xml data")
      ligomd.parse(arg[0])

      # store the users dn in the process table
      ligomd.set_dn(remote_dn)

      # determine the local creator_db number
      if creator_db is None:
        sql = "SELECT DEFAULT FROM SYSCAT.COLUMNS WHERE "
        sql += "TABNAME = 'PROCESS' AND COLNAME = 'CREATOR_DB'"
        ligomd.curs.execute(sql)
        creator_db = ligomd.curs.fetchone()[0]

      # determine the locations of columns we need in the process table
      process_cols = ligomd.table['process']['orderedcol']
      node_col = process_cols.index('node')
      prog_col = process_cols.index('program')
      upid_col = process_cols.index('unix_procid')
      start_col = process_cols.index('start_time')
      end_col = process_cols.index('end_time')
      pid_col = process_cols.index('process_id')

      # determine and remove known entries from the process table
      rmv_idx = []
      for row_idx,row in enumerate(ligomd.table['process']['stream']):
        uniq_proc = (row[node_col],row[prog_col],row[upid_col],row[start_col])
        try:
          proc_key[str(row[pid_col])] = dmt_proc_dict[uniq_proc]
          known_proc[str(dmt_proc_dict[uniq_proc])] = row[end_col]
          rmv_idx.append(row_idx)
        except KeyError:
          # we know nothing about this process, so query the database
          sql = "SELECT BLOB(process_id) FROM process WHERE "
          sql += "creator_db = " + str(creator_db) + " AND "
          sql += "node = '" + row[node_col] + "' AND "
          sql += "program = '" + row[prog_col] + "' AND "
          sql += "unix_procid = " + str(row[upid_col]) + " AND "
          sql += "start_time = " + str(row[start_col])
          ligomd.curs.execute(sql)
          db_proc_ids = ligomd.curs.fetchall()
          if len(db_proc_ids) == 0:
            # this is a new process with no existing entry
            dmt_proc_dict[uniq_proc] = row[pid_col]
          elif len(db_proc_ids) == 1:
            # the process_id exists in the database so use that insted
            dmt_proc_dict[uniq_proc] = db_proc_ids[0][0]
            proc_key[str(row[pid_col])] = dmt_proc_dict[uniq_proc]
            known_proc[str(dmt_proc_dict[uniq_proc])] = row[end_col]
            rmv_idx.append(row_idx)
          else:
            # multiple entries for this process, needs human assistance
            raise ServerHandlerException, "multiple entries for dmt process"

      # delete the duplicate processs rows and clear the table if necessary
      newstream = []
      for row_idx,row in enumerate(ligomd.table['process']['stream']):
        try:
          rmv_idx.index(row_idx)
        except ValueError:
          newstream.append(row)
      ligomd.table['process']['stream'] = newstream
      if len(ligomd.table['process']['stream']) == 0:
        del ligomd.table['process']

      # turn the known process_id binary for this insert into ascii
      for pid in known_proc.keys():
        pid_str = "x'"
        for ch in pid:
          pid_str += "%02x" % ord(ch)
        pid_str += "'"
        known_proc[pid] = (pid_str, known_proc[pid])

      # determine the locations of columns we need in the segment_definer table
      seg_def_cols = ligomd.table['segment_definer']['orderedcol']
      ifos_col = seg_def_cols.index('ifos')
      name_col = seg_def_cols.index('name')
      vers_col = seg_def_cols.index('version')
      sdid_col = seg_def_cols.index('segment_def_id')

      # determine and remove known entries in the segment_definer table
      rmv_idx = []
      for row_idx,row in enumerate(ligomd.table['segment_definer']['stream']):
        uniq_def = (row[ifos_col],row[name_col],row[vers_col])
        try:
          seg_def_key[str(row[sdid_col])] = dmt_seg_def_dict[uniq_def]
          rmv_idx.append(row_idx)
        except KeyError:
          # we know nothing about this segment_definer, so query the database
          sql = "SELECT BLOB(segment_def_id) FROM segment_definer WHERE "
          sql += "creator_db = " + str(creator_db) + " AND "
          sql += "ifos = '" + row[ifos_col] + "' AND "
          sql += "name = '" + row[name_col] + "' AND "
          sql += "version = " + str(row[vers_col])
          ligomd.curs.execute(sql)
          db_seg_def_id = ligomd.curs.fetchall()
          if len(db_seg_def_id) == 0:
            # this is a new segment_defintion with no existing entry
            dmt_seg_def_dict[uniq_def] = row[sdid_col]
          else:
            dmt_seg_def_dict[uniq_def] = db_seg_def_id[0][0]
            seg_def_key[str(row[sdid_col])] = dmt_seg_def_dict[uniq_def]
            rmv_idx.append(row_idx)

      # delete the necessary rows. if the table is empty, delete it
      newstream = []
      for row_idx,row in enumerate(ligomd.table['segment_definer']['stream']):
        try:
          rmv_idx.index(row_idx)
        except ValueError:
          newstream.append(row)
      ligomd.table['segment_definer']['stream'] = newstream
      if len(ligomd.table['segment_definer']['stream']) == 0:
        del ligomd.table['segment_definer']

      # now update the values in the xml with the values we know about
      for tabname in ligomd.table.keys():
        table = ligomd.table[tabname]
        if tabname == 'process':
          # we do nothing to the process table
          pass
        elif tabname == 'segment' or tabname == 'segment_summary':
          # we need to update the process_id and the segment_def_id columns
          pid_col = table['orderedcol'].index('process_id')
          sdid_col = table['orderedcol'].index('segment_def_id')
          row_idx = 0
          for row in table['stream']:
            try:
              repl_pid = proc_key[str(row[pid_col])]
            except KeyError:
              repl_pid = row[pid_col]
            try:
              repl_sdid = seg_def_key[str(row[sdid_col])]
            except KeyError:
              repl_sdid = row[sdid_col]
            row = list(row)
            row[pid_col] = repl_pid
            row[sdid_col] = repl_sdid
            table['stream'][row_idx] = tuple(row)
            row_idx += 1
        else:
          # we just need to update the process_id column
          pid_col = table['orderedcol'].index('process_id')
          row_idx = 0
          for row in table['stream']:
            try:
              repl_pid = proc_key[str(row[pid_col])]
              row = list(row)
              row[pid_col] = repl_pid
              table['stream'][row_idx] = tuple(row)
            except KeyError:
              pass
            row_idx += 1

      # insert the metadata into the database
      logger.debug("inserting xml data")
      result = str(ligomd.insert())
      logger.debug("insertion complete")

      # update the end time of known processes in the process table
      for pid in known_proc.keys():
        # first check to see if we are backfilling missing segments
        sql = "SELECT end_time,domain FROM process "
        sql += " WHERE process_id = " + known_proc[pid][0]
        ligomd.curs.execute(sql)
        last_end_time = ligomd.curs.fetchone()

        # check the dn in the row we are about to update matches the users dn
        dn = last_end_time[1].strip()
        if remote_dn != dn:
          msg = "%s does not have permission to update row entries" % remote_dn
          msg += " created by %s (process_id %s)" % (dn, known_proc[pid][0])
          raise ServerHandlerException, msg
        else:
          logger.debug('"%s" updating process_id %s' % (dn, known_proc[pid][0]))

        if int(known_proc[pid][1]) <= int(last_end_time[0]):
          logger.debug("Backfilling missing segments for process_id " +
            known_proc[pid][0] + " not updating end_time")
        else:
          # if we are not backfilling, update the end_time of the process
          sql = "UPDATE process SET end_time = " + str(known_proc[pid][1])
          sql += " WHERE process_id = " + known_proc[pid][0]
          sql += " AND end_time < " + str(known_proc[pid][1])
          ligomd.curs.execute(sql)
      ligomd.dbcon.commit()

      logger.info("Method insert: %s rows affected by insert" % result)
      code = 0
    except Exception, e:
      result = ("Error inserting metadata into database: %s" % e)
      logger.error(result)

    try:
      del ligomd
      del lwtparser
      del known_proc
      del seg_def_key
      del proc_key
    except Exception, e:
      logger.error(
        "Error deleting metadata object in method insertdmt: %s" % e)

    return (code,result)

