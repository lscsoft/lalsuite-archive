"""
The LSCsegFindServer module provides an API for responding to request from the
LSCsegFindClient by connecting to the DB2 database.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.

$Id$
"""

__version__ = '$Revision$'[11:-2]

import os
import sys
import re
import imp
import types
import copy
import exceptions
import socket
import SocketServer
import cPickle
import mx.ODBC.DB2
from mx.ODBC.DB2 import SQL

def initialize(configuration,log):
  # define the global variables used by the server
  global logger, max_bytes, db
  
  # initialize the logger
  logger = log
  logger.info("Initializing server module %s" % __name__ )
  
  # initialize the database hash table
  db = mx.ODBC.DB2.Connect(configuration['dbname'])
  max_bytes = configuration['max_client_byte_string']

def shutdown():
  global logger, max_bytes, db
  logger.info("Shutting down server module %s" % __name__ )
  try:
    db.close()
    del db
  except:
    pass


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

    logger.debug("handle method of %s class called" % __name__)

    # mapping of LSCsegFindServer RPC protocol names to methods of this class
    methodDict = {
      'PING' : self.ping,
      'DISTINCT' : self.distinctAttribute,
      'METASEGS' : self.segmentFindWithMetadata_v1,
      'METASEGSVX' : self.segmentFindWithMetadata_vx
    }

    try:
      # from the socket object create a file object
      self.sfile = self.request.makefile("rw")
      f = self.sfile

      # read all of the input up to limited number of bytes
      input = f.read(size=max_bytes,waitForBytes=2)
      if input[-1] != '\0':
        input += f.read(size=max_bytes,waitForBytes=2)

      # the format should be a method string, followed by a null byte
      # followed by the arguments to the method encoded as null
      # terminated strings

      # check if the last byte is a null byte
      if input[-1] != '\0':
        raise ServerHandlerException, "Last byte of input is not null byte"
    except Exception, e:
      logger.error("Error reading input on socket: %s" %  e)
      return

    logger.debug("Input on socket: %s" % input[0:-1])

    try:
      # parse out the method and arguments 
      stringList = input.split('\0')
      methodString = stringList[0]
      argStringList = stringList[1:-1]
                        
    except Exception, e:
      logger.error("Error parsing method and argument string: %s" % e)

      msg = "ERROR LSCsegFindServer Error: " + \
        "Error parsing method and argument string: %s" % e
      self.__reply__(1, msg)
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
      if result is not None:
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
    pmsg = cPickle.dumps(msg)

    f = self.sfile
    reply = "%d\0%s\0" % (code, pmsg)
    f.write(reply)

    # close the file associated with the socket
    f.close()

  def ping(self, arg):
    """
    Bounce back alive statment. Corresponds to the PING method in the
    LSCsegFindServer RPC protocol.

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

  def distinctAttribute(self, arg):
    """
    Find distinct values for an attribute. Corresponds to the DISTINCT method
    in the LSCsegFindServer RPC protocol.

    @param arg: list of strings representing messages sent by client; first
      one should be the attribute to be examined

    @return: None
    """
    global logger, db

    logger.debug("Method distinctAttribute called")

    # assume failure
    code = 1

    try:
      attribute = arg[0]

      sql = "SELECT DISTINCT "
      if attribute == "interferometers":
        sql += "ifos FROM segment_definer"
      elif attribute == "state":
        sql += "name FROM segment_definer WHERE state_vec_major IS NOT NULL "
        sql += "AND state_vec_minor IS NOT NULL"
      elif attribute == "explainstate":
        sql += "name,comment FROM segment_definer WHERE state_vec_major IS NOT NULL "
        sql += "AND state_vec_minor IS NOT NULL"
      else:
        msg = "Unknown select distinct method " + str(attribute)
        raise ServerHandlerException, msg

      c = db.cursor()
      c.execute(sql)
      res = c.fetchall()
      c.close()

      result = ""
      for x in res:
        if isinstance(x, types.StringTypes):
          result += x.strip() + '\n'
        else:
          result += str(x) + '\n'
          
      result = result.rstrip()    
      
      logger.debug("Method distinctAttribute: %d results found" % 
        len(result))
      code = 0
    except Exception, e:
      result = "Error querying metadata for distinct values of attribute " + \
        "%s: %s" % (attribute, e)
      logger.error(result)

    return (code, result)

  def segmentFindWithMetadata_v1(self, arg):
    """
    Given a SQL WHERE-type clause find segment(s) that have the matching
    attributes. Corresponds to the METASEGS method in the LSCsegFindServer RPC
    protocol.

    @param arg: list of strings representing messages sent by the client; see
    the RPC protocol 

    @return: None
    """
    global logger, db

    logger.debug("Method segmentFindWithMetadata_v1 called")
    time_rx = \
      re.compile(r"state_segment\.start_time\s*BETWEEN\s*(\d*)\s*AND\s*(\d*)")
    ifo_rx = re.compile(r"state_segment\.ifo\s*=\s*'([A-Za-z0-9]*)'")
    state_rx = re.compile(r"state_vec\.state = '([A-Za-z0-9_.,]+)'")
    try:
      start_time = time_rx.search(arg[0]).group(1)
      end_time = time_rx.search(arg[0]).group(2)
      ifo = ifo_rx.search(arg[0]).group(1)
      states = ','.join(state_rx.findall(arg[0]))
      vx_query = ['1',start_time,end_time,ifo,states]

      self.segmentFindWithMetadata_vx(vx_query)

    except Exception, e:
      msg = "Error parsing version 1 protocol string %s : %s" % (arg[0],e)
      raise ServerHandlerException, msg

    return None

  def segmentFindWithMetadata_vx(self, arg):
    """
    Given a list of attributes to query on, the first of which is the segfind
    client/server communication protocol version.  Corresponds to the
    METASEGSVX method in the LSCsegFindServer RPC protocol.

    @param arg: list of strings representing messages sent by the client; see
    the RPC protocol 

    @return: None
    """

    logger.debug("Method segmentFindWithMetadata_vx called")
    code = 0

    # get the protocol version
    try:
      protocol = int(arg[0])
    except Exception, e:
      result = "Error parsing protocol version %s: %s" % (str(arg), e)
      logger.error(result)
      code = 1
    if protocol > 3:
      result = "Unknown protocol version %d" % (protocol)
      logger.error(result)
      code = 1

    if protocol < 3:
      # load revision 1.7 of the segments module
      logger.debug("Importing segments verson 1.7 from glue")
      m_fp, m_path, m_desc = imp.find_module('segments',
        [os.path.join(
          os.environ["GLUE_LOCATION"],
          'lib/python/glue/segfindserver/segments_1_7')
        ])
      segments = imp.load_module('glue.segments',m_fp,m_path,m_desc)
      m_fp.close()
    
    else:
      # load the current revision of the segments module
      logger.debug("Importing segments from glue")
      segments = __import__('glue.segments', globals(), locals(), ['segments'])

    # parse out query information
    try:
      start = arg[1]
      end = arg[2]
      interferometer = arg[3]
      type = arg[4]
    except Exception, e:
      result = "Error parsing query arguments from list %s: %s" % (str(arg), e)
      logger.error(result)
      code = 1

    # split the interferometers and types into lists
    try:
      ifoList = interferometer.split(',')
      typeList = type.split(',')
    except:
      result = "Unable to parse ifo or type string %s : %s" % (str(arg), e)
      logger.error(result)
      code = 1
    
    # only executle a query if we managed to construct one
    if code is not 1:

      logger.debug("Query parameters [%s,%s); %s; %s" % 
        (start,end,str(ifoList),str(typeList)))

      # assume failure
      code = 1

      ifoSegList = []

      try:
        c = db.cursor()

        for ifo in ifoList:
          typeList_tmp = copy.deepcopy(typeList)
          sql = "SELECT segment.start_time, segment.end_time FROM "
          sql += "segment,segment_def_map,segment_definer WHERE "
          sql += "segment.segment_id = segment_def_map.segment_id AND "
          sql += "segment_def_map.segment_def_id = segment_definer.segment_def_id AND"
          sql += "((segment.start_time BETWEEN %s AND %s) OR " % \
            (start, end)
          sql += "(segment.end_time BETWEEN %s AND %s)) " % \
            (start, end)
          sql += "AND segment_definer.ifos = '%s' AND (segment_definer.name = '%s'" % \
            (ifo, typeList_tmp.pop(0))
          for x in typeList_tmp:
            sql += " OR segment_definer.name = '%s'" % (x)
          sql += ")"
          del typeList_tmp

          c.execute(sql)
          result = c.fetchall()
          logger.debug("Method : segmentFindWithMetadata_vx %d results found" % 
            len(result))
          c.close()

          ifoSegs = segments.segmentlist()
          try:
            for r in result:
              ifoSegs.append(segments.segment(r[0],r[1]))
            ifoSegs.coalesce()
          except:
            pass
          ifoSegList.append(ifoSegs)

      except Exception, e:
        result = "Error querying metadata for segments " + \
          "%s: %s" % (arg[0], e)
        logger.error(result)

      # take the intersection if there is more than one ifo queried
      result = ifoSegList.pop(0)
      try:
        for iseg in ifoSegList: result &= iseg
      except:
        pass

      code = 0
    else:
      logger.error( "segmentFindWithMetadata_vx skipped query" )

    self.__reply__(code, result)

    del result

    # unload the segments module
    logger.debug("Unloading segments module")
    del segments
    del sys.modules["glue.segments"]

    return None
