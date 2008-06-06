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
from glue import segments
import cPickle
import mx.ODBC.DB2

def initialize(configuration,log):
  # define the global variables used by the server
  global logger, max_bytes, db, dbname
  
  # initialize the logger
  logger = log
  logger.info("Initializing server module %s" % __name__ )
  
  # initialize the database hash table
  dbname = configuration['dbname']
  db = mx.ODBC.DB2.Connect(dbname)
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
    global logger, db, dbname

    logger.debug("Method distinctAttribute called")

    # assume failure
    code = 1

    try:
      attribute = arg[0]

      sql = "SELECT DISTINCT "
      if attribute == "interferometers":
        sql += "ifos AS x FROM segment_definer "
      elif attribute == "state":
        sql += "name AS x FROM segment_definer "
      elif attribute == "explainstate":
        sql += "run AS x, ifos, name, version, comment FROM segment_definer "
      else:
        msg = "Unknown select distinct method " + str(attribute)
        raise ServerHandlerException, msg
      sql += "ORDER BY x ASC"

      try:
        c = db.cursor()
      except mx.ODBC.DB2.OperationalError, e:
        logger.debug( "OperationalError: %s" % str(e) )
        if ( int(e[0]) == 8003 and e[1] == -99999 ):
          logger.info("Reconnecting to database due to error %s" % str(e))
          db = mx.ODBC.DB2.Connect(dbname)
          c = db.cursor()
        else:
          raise
      except:
        logger.debug( "Unhandled Error: %s" % str(e) )
        raise

      c.execute(sql)
      res = c.fetchall()
      c.close()

      result = ""
      for x in res:
        if len(x) == 1:
          result += x[0].strip() + '\n'
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
    global logger, db, dbname

    logger.debug("Method segmentFindWithMetadata_v1 called")
    logger.debug("arg="+repr(arg))


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
    global logger, db, dbname

    logger.debug("Method segmentFindWithMetadata_vx called")
    code = 0

    logger.error(repr(arg))

    # get the protocol version
    try:
      protocol = int(arg[0])
    except Exception, e:
      result = "Error parsing protocol version %s: %s" % (str(arg), e)
      logger.error(result)
      code = 1


    if protocol > 4:
      result = "Unknown protocol version %d" % (protocol)
      logger.error(result)
      code = 1

    if protocol < 3:
      # return an error message
      result = "Protocol versions less than 3 are not supported by this server"
      logger.error(result)
      code = 1
    
    # parse out query information
    try:
      if(protocol==3):
        start = arg[1]
        end = arg[2]
        interferometer = arg[3]
        type = arg[4]
      elif(protocol==4):
        run=arg[1]
        interferometer = arg[2]
        type = arg[3]
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
      if(protocol==3):
        logger.debug("Query parameters [%s,%s); %s; %s" % 
                     (start,end,str(ifoList),str(typeList)))
      elif(protocol==4):
        logger.debug("Query parameters %s %s; %s" % 
                     (run,str(ifoList),str(typeList)))

      # assume failure
      code = 1

      try:
        try:
          c = db.cursor()
        except mx.ODBC.DB2.OperationalError, e:
          logger.debug( "OperationalError: %s" % str(e) )
          if ( int(e[0]) == 8003 and e[1] == -99999 ):
            logger.info("Reconnecting to database due to error %s" % str(e))
            db = mx.ODBC.DB2.Connect(dbname)
            c = db.cursor()
          else:
            raise
        except:
          logger.debug( "Unhandled Error: %s" % str(e) )
          raise

        if(protocol==4):
          logger.debug("Protocol 4")
          result=self.segmentFindWithMetadata_vx_p4(c, run, ifoList, typeList)
        elif(protocol==3):
          logger.debug("Protocol 3")
          if(len(typeList)==1 and typeList[0]=='All'):
            result=self.segmentFindWithMetadata_vx_p3_all(c, start, end, ifoList)
          else:
            result=self.segmentFindWithMetadata_vx_p3(c, start, end, ifoList, typeList)
            
        # close the database cursor
        c.close()
      except Exception, e:
        result = "Error querying metadata for segments " + \
          "%s: %s" % (arg[0], e)
        logger.error(result)

      # take the intersection of the query result with the requested range
      if(protocol==3):
        query_range = segments.segmentlist()
        query_range.append(segments.segment(int(start),int(end)))
        try:
          result &= query_range
        except:
          pass
        
      code = 0
    else:
      logger.error( "segmentFindWithMetadata_vx skipped query" )

    self.__reply__(code, result)
    del result
    return None


  def segmentFindWithMetadata_vx_p4(self, c, run, ifoList, typeList):
    global logger
    logger.debug("segmentFindWithMetadata_vx_p4 called")
    
    iSegs=[]

    for ifo in ifoList:
      typeSegs=segments.segmentlist()
      for tt in typeList:          
        only=0
        if(tt.find(":")!=-1):
          try:
            [t,v,only]=tt.split(":")
            only=int(only)
            v=int(v)
          except:
            [t,v]=tt.split(":")
            v=int(v)
          maxversion=int(v)
        else:
          t=tt
          sql = "SELECT max(version) FROM segment_definer where run='%s' and ifos= '%s' and name = '%s'" % (run,ifo,t)
          logger.debug(sql)
          c.execute(sql)
          result=c.fetchall()
          maxversion=result[0][0]
        logger.debug("maxversion="+str(maxversion)+ " for ifo="+ifo+" and type="+t)
            
        sql = "SELECT segment.start_time, segment.end_time FROM "
        sql += "segment,segment_def_map,segment_definer A WHERE "
        sql += "segment.active = 1 AND "
        sql += "segment.segment_id = segment_def_map.segment_id AND "
        sql += "segment.creator_db = segment_def_map.segment_cdb AND "
        sql += "segment_def_map.segment_def_id = A.segment_def_id AND "
        sql += "segment_def_map.segment_def_cdb = A.creator_db AND "
        sql += "A.run = '%s' AND A.ifos = '%s' AND A.name = '%s' AND A.version=%d" % (run, ifo, t, maxversion)

        logger.debug(sql)
        c.execute(sql)
        result=c.fetchall()

        mergedSegs=segments.segmentlist()

        try:
          for r in result:
            mergedSegs.append(segments.segment(r[0],r[1]))
          mergedSegs.coalesce()
          logger.debug("len(mergedSegs)="+repr(len(mergedSegs)))
        except:
          logger.debug("Why did mergeSegs try failed?")

        if(only==0):
          try:
            sql = "SELECT min(segment.start_time), max(segment.end_time) FROM "
            sql += "segment,segment_def_map,segment_definer A WHERE "
            sql += "segment.active = 1 AND "
            sql += "segment.segment_id = segment_def_map.segment_id AND "
            sql += "segment.creator_db = segment_def_map.segment_cdb AND "
            sql += "segment_def_map.segment_def_id = A.segment_def_id AND "
            sql += "segment_def_map.segment_def_cdb = A.creator_db AND "
            sql += "A.run = '%s' AND A.ifos = '%s' AND A.name = '%s'" % (run, ifo, t)
            logger.debug(sql)
            c.execute(sql)
            result=c.fetchall()
            runMaxRange=segments.segment(result[0][0],result[0][1])
            logger.debug("runMaxRange="+repr(runMaxRange))
            maxVRange=mergedSegs.extent()
            logger.debug("maxVRange="+str(maxVRange))
            remainderRange1=segments.segmentlist([runMaxRange])-segments.segmentlist([maxVRange])
            logger.debug("remainderRange1="+str(remainderRange1))
            if(len(remainderRange1)>0):
              remainderRange=remainderRange1[-1]
            else:
              remainderRange=segments.segment(-1,-1)
          except:
            logger.debug("Failed to find runMaxRange")

          newSegs=segments.segmentlist()            
          if(abs(remainderRange)>0):
            sql = "SELECT segment.start_time, segment.end_time FROM "
            sql += "segment,segment_def_map,segment_definer WHERE "
            sql += "segment.active = 1 AND "
            sql += "segment.segment_id = segment_def_map.segment_id AND "
            sql += "segment.creator_db = segment_def_map.segment_cdb AND "
            sql += "segment_def_map.segment_def_id = segment_definer.segment_def_id AND "
            sql += "segment_def_map.segment_def_cdb = segment_definer.creator_db AND "
            sql += "segment.start_time <= %d AND " % remainderRange[1]
            sql += "segment.end_time >= %d " % remainderRange[0]
            sql += "AND segment_definer.ifos = '%s' AND segment_definer.name = '%s'" % \
                   (ifo, t)
            
            logger.debug(sql)
            c.execute(sql)
            result=c.fetchall()
              
            try:
              for r in result:
                newSegs.append(segments.segment(r[0],r[1]))
              newSegs.coalesce()
              logger.debug("len(newSegs)="+repr(len(newSegs)))
            except:
              pass
        elif(only==1):
          newSegs=segments.segmentlist([])
        segs=mergedSegs | newSegs  
        typeSegs|=segs
            
      typeSegs.coalesce()
      logger.debug("len(typeSegs)="+repr(len(typeSegs)))
      iSegs.append(typeSegs)
      
    logger.debug("len(iSegs)="+repr(len(iSegs)))
    result=reduce(lambda x,y: x & y,iSegs)
    result.coalesce()
    logger.debug("len(result)="+repr(len(result)))
    return result

  def segmentFindWithMetadata_vx_p3_all(self, c, start, end, ifoList):
    global logger
    logger.debug("segmentFindWithMetadata_vx_p3_all called")
    
    tvseList={}
    logger.debug("All segment types and versions are requested")
    try:
      sql = "SELECT A.ifos, A.name, A.version, B.start_time, B.end_time FROM "
      sql += "segment B,segment_def_map,segment_definer A WHERE "
      sql += "B.active = 1 AND "
      sql += "B.segment_id = segment_def_map.segment_id AND "
      sql += "B.creator_db = segment_def_map.segment_cdb AND "
      sql += "segment_def_map.segment_def_id = A.segment_def_id AND "
      sql += "segment_def_map.segment_def_cdb = A.creator_db AND "
      sql += "B.start_time <= %s AND " % end
      sql += "B.end_time >= %s " % start
      sql += "order by A.ifos, A.name, A.version, B.start_time"

      logger.debug(sql)
      c.execute(sql)
      result=c.fetchall()
          
      for r in result:
        logger.debug("r="+repr(r))
        try:
          tvseList[(r[0].strip(),r[1],r[2])].append(segments.segment(r[3],r[4]))
        except:
          tvseList[(r[0].strip(),r[1],r[2])]=segments.segmentlist([segments.segment(r[3],r[4])])
              
      kkk=tvseList.keys()
      kkk.sort()
      logger.debug("kkk="+repr(kkk))
      print kkk
      for kkkk in kkk:
        tvseList[kkkk].coalesce()
      result=tvseList
      code=0
    except:
      raise
    return result
                                                                                                                                           
  def segmentFindWithMetadata_vx_p3(self, c, start, end, ifoList, typeList):
    global logger
    logger.debug("segmentFindWithMetadata_vx_p3 called")    
    iSegs=[]

    for ifo in ifoList:
      typeSegs=segments.segmentlist()
      for tt in typeList:
        only=0
        if(tt.find(":")!=-1):
          try:
            [t,v,only]=tt.split(":")
            only=int(only)
            v=int(v)
          except:
            [t,v]=tt.split(":")
            v=int(v)
          maxversion=int(v)
        else:
          t=tt
          sql = "SELECT max(version) FROM segment_definer where ifos= '%s' and name = '%s'" % (ifo,t)
          logger.debug(sql)
          c.execute(sql)
          result=c.fetchall()
          maxversion=result[0][0]
        logger.debug("maxversion="+str(maxversion)+ " for ifo="+ifo+" and type="+t)
        
        sql = "SELECT segment.start_time, segment.end_time FROM "
        sql += "segment,segment_def_map,segment_definer A WHERE "
        sql += "segment.active = 1 AND "
        sql += "segment.segment_id = segment_def_map.segment_id AND "
        sql += "segment.creator_db = segment_def_map.segment_cdb AND "
        sql += "segment_def_map.segment_def_id = A.segment_def_id AND "
        sql += "segment_def_map.segment_def_cdb = A.creator_db AND "
        sql += "segment.start_time <= %s AND " % end
        sql += "segment.end_time >= %s " % start
        sql += "AND A.ifos = '%s' AND A.name = '%s'" % \
               (ifo, t)
        sql += "AND A.version=%d" % maxversion

        logger.debug(sql)
        c.execute(sql)
        result=c.fetchall()

        mergedSegs=segments.segmentlist()
        try:
          for r in result:
            mergedSegs.append(segments.segment(r[0],r[1]))
          mergedSegs.coalesce()
          logger.debug("len(mergedSegs)="+repr(len(mergedSegs)))
        except:
          logger.debug("megedSegs failed")

        if(only==0):
          try:
            maxVRange=mergedSegs.extent()
            logger.debug("maxVRange="+str(maxVRange))
            remainderRange1=segments.segmentlist([segments.segment(int(start),int(end))])-segments.segmentlist([maxVRange])
            logger.debug("remainderRange1="+str(remainderRange1))            
            remainderRange=remainderRange1[-1]
          except:
            remainderRange=segments.segment(int(start),int(end))
            logger.debug("remainderRange="+str(remainderRange))

          newSegs=segments.segmentlist()            
          if(remainderRange[0]>int(start) and remainderRange[0]<int(end)):
            sql = "SELECT segment.start_time, segment.end_time FROM "
            sql += "segment,segment_def_map,segment_definer WHERE "
            sql += "segment.active = 1 AND "
            sql += "segment.segment_id = segment_def_map.segment_id AND "
            sql += "segment.creator_db = segment_def_map.segment_cdb AND "
            sql += "segment_def_map.segment_def_id = segment_definer.segment_def_id AND "
            sql += "segment_def_map.segment_def_cdb = segment_definer.creator_db AND "
            sql += "segment.start_time <= %d AND " % remainderRange[1]
            sql += "segment.end_time >= %d " % remainderRange[0]
            sql += "AND segment_definer.ifos = '%s' AND segment_definer.name = '%s'" % \
                   (ifo, t)

            logger.debug(sql)
            c.execute(sql)
            result=c.fetchall()
              
            try:
              for r in result:
                newSegs.append(segments.segment(r[0],r[1]))
              newSegs.coalesce()
              logger.debug("len(newSegs)="+repr(len(newSegs)))
            except:
              pass
        elif(only==1):
          newSegs=segments.segmentlist()
        segs=mergedSegs | newSegs
        logger.debug("len(segs)="+repr(len(segs)))
        typeSegs|=segs
            
      typeSegs.coalesce()
      logger.debug("len(typeSegs)="+repr(len(typeSegs)))
    iSegs.append(typeSegs)
    logger.debug("len(iSegs)="+repr(len(iSegs)))
                   
    result=reduce(lambda x,y: x & y,iSegs)
    result.coalesce()

    logger.debug("len(result)="+repr(len(result)))
    
    return result
  

