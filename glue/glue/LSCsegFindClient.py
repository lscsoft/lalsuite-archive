"""
The LSCsegFindClient module provides an API for connecting to
and making requests of a LSCsegFindServer.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.

$Id$
"""

__version__ = '$Revision$'[11:-2]

import sys
import os
import exceptions
import types
import cPickle

from pyGlobus import io
from pyGlobus import security

from glue import segments


def version():
  return __version__


class LSCsegFindException(exceptions.Exception):
  """
  Exceptions raised by the classes and methods in this client
  will be instances of this class.
  """
  def __init__(self, args=None):
    """
    Create an instance of this class, ie. an LSCSegFindException.

    @param args:

    @return: Instance of class LSCsegFindException
    """
    self.args = args


class LSCsegFindClientException(exceptions.Exception):
  """
  Exceptions raised by the classes and methods in this module
  will be instances of this class.
  """
  def __init__(self, args=None):
    """
    Create an instance of this class, ie. an LSCsegFindClient
    exception.

    @param args: 

    @return: Instance of class LSCsegFindClientException
    """
    self.args = args

class LSCsegFindClient(object):
  """
  Class that represents a client interacting with a LSCsegFindServer. It is
  expected that clients for interacting with a LSCsegFindServer will inherit
  from this base class.
  """
  def __init__(self, host, port):
    """
    Open a connection to a LSCsegFindServer and return an instance of
    class LSCsegFindClient. One of the public methods can then be 
    called to send a request to the server.

    @param host: the host on which the LSCsegFindServer runs
    @type host: string

    @param port: port on which the LSCsegFindServer listens
    @type port: integer


    @return: Instance of LSCsegFindClient
    """
    # check arguments
    if type(host) != types.StringType:
      msg = "Argument 'host' must be a string"
      raise LSCsegFindClientException, msg

    if type(port) != types.IntType:
      msg = "Argument 'port' must be a positive integer"
      raise LSCsegFindClientException, msg

    if port <= 0:
      msg = "Argument 'port' must be a positive integer"
      raise LSCsegFindClientException, msg
                
    # try to connect and if there are any exceptions catch them
    # and turn them into LSCsegFindClientException
    try:
      self.__connect__(host, port)
    except Exception, e:
      raise LSCsegFindClientException, e


  def __del__(self):
    """
    Disconnect from the LSCsegFindServer.

    @return: None
    """
    self.__disconnect__()


  def __connect__(self, host, port):
    """
    Attempt to open a connection to the LSCsegFindServer
    using the 'host' and 'port' and expecting the server
    to identify itself with a corresponding host certificate.

    A IOException is raised if the connection cannot be made,
    but this is caught by the __init__ method above and 
    turned into a LSCdataFindClient exception.
        
    @param host: the host on which the LSCsegFindServer runs
    @type host: string

    @param port: port on which the LSCsegFindServer listens
    @type port: integer

    @return: None
    """

    self.host = host
    self.port = port

    # redirect stdout and stderror for now
    try:
      f = open("/dev/null", "w")
      sys.stdout = f
      sys.stderr = f
    except:
      pass

    try:
      # create TCPIOAttr instance and set authentication mode to be GSSAPI 
      myAttr = io.TCPIOAttr()
      myAttr.set_authentication_mode(
        io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_GSSAPI)

      # create AuthData instance
      authData = io.AuthData()
      authData.set_identity("/DC=org/DC=doegrids/OU=Services/CN=lscsegfind/%s" 
        % (host))

      # set authorization, channel, and delegation modes
      myAttr.set_authorization_mode(
        io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_IDENTITY, authData)
      myAttr.set_channel_mode(io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
      myAttr.set_delegation_mode(
        io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_LIMITED_PROXY)

      # create socket instance and attempt to connect
      s = io.GSITCPSocket()
      s.connect(host, port, myAttr)
      self.socket = s
      self.sfile = s.makefile("rw")

    finally:
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__
      f.close()

  def __disconnect__(self):
    """
    Disconnect from the LSCsegFindServer.

    @return: None
    """
    try:
      self.socket.shutdown(2)
    except:
      pass

  def __response__(self):
    """
    Read the response sent back by the LSCsegFindServer. Parse out the
    return code with 0 for success and non-zero for error, and then
    the list of strings representing the returned result(s).

    @return: tuple containing the integer error code and the list of 
    strings representing the output from the server
    """
    f = self.sfile
       
    response = ""

    # read in 512 byte chunks until there is nothing left to read
    # this blocks until the socket is ready for reading
    while 1: 
      input = f.read(512)
      if input == "": break
      response += input
    
    # the response from the server must always end in a null byte
    if response[-1] != '\0':
      msg = "Bad format for response from server"
      raise LSCsegFindClientException, msg

    # delete the last \0 before splitting into strings
    response = response[0:-1]

    try:
      stringList = response.split('\0')
      if len(stringList) is not 2:
        msg = "Malformatted response from server"
        raise LSCsegFindClientException, msg
      code = int(stringList[0])
      output = cPickle.loads(stringList[1])
    except Exception, e:
      msg = "Error parsing response from server : %s" % e
      try:
        f.close()
      except:
        pass
      raise LSCsegFindClientException, msg

    f.close()

    return code, output

  def ping(self):
    """
    Ping the LSCsegFindServer and return any message received back as
    a string.

    @return: message received (may be empty) from LSCsegFindServer as a string
    """

    msg = "PING\0"
    self.sfile.write(msg)

    ret, output = self.__response__()

    if ret:
      msg = "Error pinging server %d:%s" % (ret, reply)
      raise LSCsegFindClientException, msg

    return output


  def distinctAttrValues(self, attr):
    """
    Query LSCsegFindServer metadata tables for the distince values of an
    attribute and return the values as a list.

    Note that the values will always be returned as strings. Any conversion
    should be done by the client.

    @param attr: name of attribute for which to find distinct values

    @return: list of strings representing the distinct values that the attribute
    has in the metadata tables
    """

    # check argument
    if type(attr) != types.StringType:
      msg = "Argument 'attr' must be a string"
      raise LSCsegFindClientException, msg
        
    msg = "DISTINCT\0%s\0" % str(attr)
    self.sfile.write(msg)

    ret, output = self.__response__()

    if ret:
      msg = \
        "Error querying LSCsegFindServer for distinct values of attributes: %s"\
        % str(output)
      raise LSCsegFindClientException, msg

    return output

  def segmentQueryWithMetadata(self, queryList):
    """
    Query LSCsegFindServer to find the segment(s) with the appropriate
    metadata values.

    @param queryList: list of instances of the query parameters, each
      describing a query to be done
                        
    @return: list of strings representing the union of all segments(s) found
      from the queries
    """

    # check arguments
    if not isinstance(queryList, list):
      msg = "Argument must be a list of query strings"
      raise LSCsegFindClientException, msg
    for query in queryList:
      if not isinstance(query, str):
        msg = "Argument must be a query string"
        raise LDRdataFindClientException, msg

    # prepare the messange to send down the socket
    msg = "METASEGS\0"
    for q in queryList:
      msg += "%s" % q
    msg += "\0"

    self.sfile.write(msg)

    ret, output = self.__response__()

    if ret:
      msg = \
        "Error querying LSCsegFindServer for segments with metadata %s : %s" \
        % (str(queryList), str(output[0]))
      raise LSCsegFindClientException, msg

    return output


class LSCsegFind(LSCsegFindClient):
  """
  Class that represents this client interacting with a LSCsegFindServer in
  order to find state segments.
  """
  def __init__(self, host, port):
    """
    Open a connection to a LSCsegFindServer and return an instance of
    class LSCsegFind. One of the public methods can then be called to send a
    request to the server.

    @param host: the host on which the LSCsegFindServer runs
    @type host: string

    @param port: port on which the LSCsegFindServer listens
    @type port: integer


    @return: Instance of LSCsegFindClient
    """
    LSCsegFindClient.__init__(self, host, port)

  def __check_gps(self, gpsString):
    """
    Minimal checking on GPS time strings. Raises a LSCdataFindClientException if
    the GPS time string is not 9 digits long.

    @param gpsString: The string representing the 9 digit GPS time.

    @returns: None
    """
    if len(gpsString) != 9:
      msg = "GPS times must be 9 digits"
      raise LSCsegFindException, msg

    try:
      a = int(gpsString)
    except Exception, e:
      msg = "GPS times must be 9 digits"
      raise LSCsegFindException, msg


  def ping(self, argDict):
    """
    Ping the LSCsegFindServer and print any response sent back.

    @param argDict: Dictionary of arguments passed to all methods.

    @return: None
    """
    response = LSCsegFindClient.ping(self)
    return str(response)


  def showInterferometers(self, argDict):
    """
    Query LSCsegFindServer for the distinct values for the 'state_segment.ifo'
    attribute in the metadata table.

    @param argDict: Dictionary of arguments passed to all methods.

    @return: None
    """

    response = LSCsegFindClient.distinctAttrValues(self, "interferometers")
    return str(response)


  def showTypes(self, argDict):
    """
    Query LSCsegFindServer for the distinct values for the 'state_vec.state'
    attribute in the metadata table.

    @param argDict: Dictionary of arguments passed to all methods.

    @return: None
    """
    response = LSCsegFindClient.distinctAttrValues(self, "state")
    return str(response)
                

  def findStateSegments(self, argDict):
    """
    Query the LDRdataFindServer for state segments from a particular
    interferometer, with a particular state type, for a particular range of
    GPS times.

    @param argDict: Dictionary of arguments passed to all methods.

    @return: None
    """
    interferometer = argDict['interferometer']
    type = argDict['type']
    start = argDict['start']
    end = argDict['end']
    type = argDict['type']
    strict = argDict['strict']
    format = argDict['format']

    # check that combination of command-line arguments is sound
    if (not interferometer) or (not type) or (not start) or (not end):
      msg = """\
Bad combination of command line arguments:
--interferometer --type --gps-start-time --gps-end-time must all
be present when searching for groups of segments 
"""
      raise LSCsegFindException, msg

    self.__check_gps(start)
    self.__check_gps(end)

    try:
      typeList = type.split(',')
    except:
      msg = "Unable to parse segment type string %s : %s" % (typeString, e)
      raise LSCsegFindException, msg

    query = "((state_segment.start_time BETWEEN %s AND %s) OR " % \
      (start, end)
    query += "(state_segment.end_time BETWEEN %s AND %s)) " % \
      (start, end)
    query += "AND state_segment.ifo = '%s' AND (state_vec.state = '%s'" % \
      (interferometer, typeList.pop(0))
    for x in typeList:
      query += " OR state_vec.state = '%s'" % (x)
    query += ")"
    
    seglist = LSCsegFindClient.segmentQueryWithMetadata(self,[query])

    if strict:
      range = segments.segmentlist([segments.segment(long(start),long(end))])
      seglist &= range

    return seglist
