"""
The LDBDServer module provides an API for responding to request from the
LDBDClient by connecting to the DB2 database.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.

$Id$
"""

__version__ = '$Revision$'[11:-2]

import re
import types
import pyRXP
import exceptions
import SocketServer
import cPickle
from glue import ldbd
import rlsClient

def initialize(configuration,log):
  # define the global variables used by the server
  global logger, max_bytes, xmlparser, dbobj, xmlparser, lwtparser, rls
  
  # initialize the logger
  logger = log
  logger.info("Initializing server module %s" % __name__ )
  
  # initialize the database hash table
  dbobj = ldbd.LIGOMetadataDatabase(configuration['dbname'])
  max_bytes = configuration['max_client_byte_string']

  # create the xml and ligolw parsers
  xmlparser = pyRXP.Parser()
  lwtparser = ldbd.LIGOLwParser()

  # open a connection to the rls server
  rls_server = configuration['rls']
  cert = configuration['certfile']
  key = configuration['keyfile']
  rls = rlsClient.RlsClient(rls_server,cert,key)

def shutdown():
  global logger, max_bytes, xmlparser, dbobj, xmlparser, lwtparser, rls
  logger.info("Shutting down server module %s" % __name__ )
  del rls
  del lwtparser
  del xmlparser
  del dbobj

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
    global logger, max_bytes

    logger.debug("handle method of LDBDServer class called")

    # mapping of ldbdd RPC protocol names to methods of this class
    methodDict = {
      'PING' : self.ping,
      'QUERY' : self.query,
      'INSERT' : self.insert,
      'INSERTMAP' : self.insertmap
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
        raise ServerHandlerException, \
          "Last byte of input is not null byte"
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

      msg = "ERROR ldbdd Error: " + \
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

    logger.debug("Method ping called")
    try:
      hostname = socket.getfqdn()
      msg = "ldbdd at %s is alive" % hostname
    except Exception, e:
      msg = "ldbdd is alive"

    return (0, msg)


  def query(self, arg):
    """
    Execute an SQL query on the database and return the result as LIGO_LW XML

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """

    # get the query string and log it
    querystr = arg[0]
    logger.debug("Method query called with %s" % querystr)

    # assume failure
    code = 1

    try:
      # create a ligo metadata object
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

    logger.debug("Method insert called")

    # assume failure
    code = 1

    try:
      # capture the remote users DN for insertion into the database
      cred = self.request.get_delegated_credential()
      remote_dn = cred.inquire_cred()[1].display()

      # create a ligo metadata object
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
    except Exception, e:
      logger.error(
        "Error deleting metadata object in method insert: %s" % e)

    return (code,result)

  def insertmap(self, arg):
    """
    Insert some LIGO_LW xml data in the metadata database with an LFN to
    PFN mapping inserted into the RLS database

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """

    logger.debug("Method insertmap called")

    # assume failure
    code = 1

    try:
      # unpickle the PFN/LFN mappings from the client
      lfnpfn_dict = cPickle.loads(arg[1])
      if not isinstance(lfnpfn_dict, dict):
        raise ServerHandlerException, \
          "LFN/PFN mapping from client is not dictionary"

      # capture the remote users DN for insertion into the database
      cred = self.request.get_delegated_credential()
      remote_dn = cred.inquire_cred()[1].display()

      # create a ligo metadata object
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # parse the input string into a metadata object
      ligomd.parse(arg[0])

      # add a gridcert table to this request containing the users dn
      ligomd.set_dn(remote_dn)

      # add the lfns to the metadata insert to populate the lfn table
      for lfn in lfnpfn_dict.keys():
        ligomd.add_lfn(lfn)

      # insert the metadata into the database
      result = str(ligomd.insert())
      logger.info("Method insert: %s rows affected by insert" % result)

      # insert the PFN/LFN mappings into the RLS
      for lfn in lfnpfn_dict.keys():
        pfns = lfnpfn_dict[lfn]
        if not isinstance( pfns, types.ListType ):
          raise ServerHandlerException, \
            "PFN must be a single string or a list of PFNs"
        rls.lrc_create_lfn( lfn, pfns[0] )
        for pfn in pfns[1:len(pfns)]:
          rls.lrc_add( lfn, pfn )
          
      logger.info("Method insertmap: insert LFN mappings for %s" % 
        str(lfnpfn_dict.keys()))
      code = 0

    except Exception, e:
      result = ("Error inserting LFN/PFN mapping into RLS: %s" % e)
      logger.error(result)
      return (code,result)

    return (code,result)
