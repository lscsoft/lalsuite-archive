"""
The LDRdataFindClient module provides an API for connecting to
and making requests of a LDRdataFindServer.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.

$Id$
"""

__version__ = '$Revision$'[11:-2]

import sys
import os
import exceptions
import types

from pyGlobus import io
from pyGlobus import security

def version():
        return __version__

class LDRMetadataQuery(object):
        """
        """
        def __init__(self, attrQuery=None):
                """
                """
                self.attrQuery = attrQuery
                self.regex = None
                self.offset = 0
                self.limit = None
                self.sortAttr = None
                self.sortOrder = None

        def set_query(self, query):
                self.attrQuery = query

        def set_regex_filter(self, filter):
                self.regex = filter

        def set_offset(self, offset):
                self.offset = offset

        def set_limit(self, limit):
                self.limit = limit

        def set_sort_attribute(self, attribute):
                self.sortAttr = attribute

        def set_sort_order(self, sortOrder):
                self.sortOrder = sortOrder

        def __verify(self):
                # we should insure here that
                # - attrQuery is a string and well behaved
                # - regex is a string and well behaved
                # - offset is an integer >= 0
                # - limit is None or integer > 0
                # - sortAttr is None or string
                # - sortOrder is one of "ASC" or "DESC"
                pass

        def __str__(self):
                self._LDRMetadataQuery__verify()

                if not self.limit: 
                        limit = "-1"
                else:
                        limit = str(self.limit)

                if not self.sortAttr:
                        sortAttr = "NONE"
                else:
                        sortAttr = self.sortAttr

                if not self.regex:
                        regex = "NONE"
                else:
                        regex = self.regex

                if not self.sortOrder:
                        sortOrder = "ASC"
                else:
                        sortOrder = self.sortOrder
                        

                s = "%s\0%s\0%d\0%s\0%s\0%s\0" % (
                        self.attrQuery,
                        regex,
                        self.offset,
                        limit,
                        sortAttr,
                        sortOrder
                        )

                return s


class LDRdataFindClientException(exceptions.Exception):
        """
        Exceptions raised by the classes and methods in this module
        will be instances of this class.
        """
        def __init__(self, args=None):
                """
                Create an instance of this class, ie. an LDRdataFindClient
                exception.

                @param args: 

                @return: Instance of class LDRdataFindClientException
                """
                self.args = args

class LDRdataFindClient(object):
        """
        Class that represents a client interacting with a LDRdataFindServer. It is expected
        that clients for interacting with a LDRdataFindServer will inherit from this base
        class.
        """
        def __init__(self, host, port):
                """
                Open a connection to a LDRdataFindServer and return an instance of
                class LDRdataFindClient. One of the public methods can then be 
                called to send a request to the server.

                @param host: the host on which the LDRdataFindServer runs
                @type host: string

                @param port: port on which the LDRdataFindServer listens
                @type port: integer


                @return: Instance of LDRdataFindClient
                """
                # check arguments
                if type(host) != types.StringType:
                        msg = "Argument 'host' must be a string"
                        raise LSCdataFindClientException, msg

                if type(port) != types.IntType:
                        msg = "Argument 'port' must be a positive integer"
                        raise LSCdataFindClientException, msg

                if port <= 0:
                        msg = "Argument 'port' must be a positive integer"
                        raise LSCdataFindClientException, msg
                 
                
                # try to connect and if there are any exceptions catch them
                # and turn them into LDRdataFindClientException
                try:
                        self.__connect__(host, port)
                except Exception, e:
                        raise LDRdataFindClientException, e


        def __del__(self):
                """
                Disconnect from the LDRdataFindServer.

                @return: None
                """
                self.__disconnect__()


        def __connect__(self, host, port):
                """
                Attempt to open a connection to the LDRdataFindServer
                using the 'host' and 'port' and expecting the server
                to identify itself with a corresponding host certificate.

                A IOException is raised if the connection cannot be made,
                but this is caught by the __init__ method above and 
                turned into a LSCdataFindClient exception.
        
                @param host: the host on which the LDRdataFindServer runs
                @type host: string

                @param port: port on which the LDRdataFindServer listens
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
                        myAttr.set_authentication_mode(io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_GSSAPI)

                        # create AuthData instance
                        authData = io.AuthData()

                        # set authorization, channel, and delegation modes
                        myAttr.set_authorization_mode(io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_HOST, authData)
                        myAttr.set_channel_mode(io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
                        myAttr.set_delegation_mode(io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_LIMITED_PROXY)

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
                Disconnect from the LDRdataFindServer.

                @return: None
                """
                try:
                        self.socket.shutdown(2)
                except:
                        pass

        def __response__(self):
                """
                Read the response sent back by the LDRdataFindServer. Parse out the
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
                        raise LDRdataFindClientException, msg

                # delete the last \0 before splitting into strings
                response = response[0:-1]

                try:
                        stringList = response.split('\0')
                        code = int(stringList[0])
                        output = stringList[1:]
                except Exception, e:
                        msg = "Error parsing response from server : %s" % e
                        try:
                                f.close()
                        except:
                                pass
                        raise LDRdataFindClientException, msg

                f.close()

                return code, output

        def ping(self):
                """
                Ping the LDRdataFindServer and return any message received back as
                a string.

                @return: message received (may be empty) from LDRdataFindServer as a string
                """

                msg = "PING\0"
                self.sfile.write(msg)

                ret, output = self.__response__()

                reply = str(output[0])

                if ret:
                        msg = "Error pinging server %d:%s" % (ret, reply)
                        raise LDRdataFindClientException, msg

                return reply


        def distinctAttrValues(self, attr):
                """
                Query LDRdataFindServer metadata tables for the distince values of an attribute
                and return the values as a list.

                Note that the values will always be returned as strings. Any conversion
                should be done by the client.

                @param attr: name of attribute for which to find distinct values

                @return: list of strings representing the distinct values that the attribute
                        has in the metadata tables
                
                """

                # check argument
                if type(attr) != types.StringType:
                        msg = "Argument 'attr' must be a string"
                        raise LDRdataFindClientException, msg
        
                msg = "DISTINCT\0%s\0" % str(attr)
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for distinct values of attributes: %s" % str(output)
                        raise LDRdataFindClientException, msg

                return output

        def pfnQuery(self, lfn):
                """
                Query LDRdataFindServer to find the PFN(s) associated with a LFN and return them
                as a list. More then one PFN may be returned.

                @param lfn: the LFN with which to query the LDRdataFindServer

                @return: list of strings representing the PFN(s) for the LFN
                """

                # check argument
                if type(lfn) != types.StringType:
                        msg = "Argument 'lfn' must be a string"
                        raise LDRdataFindClientException, msg

                msg = "LFNPFN\0%s\0" % str(lfn)
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFN with LFN %s : %s" % (str(lfn), str(output))
                        raise LDRdataFindClientException, msg

                return output

        def lfnQueryWithMetadata(self, queryList):
                """
                Query LDRdataFindServer to find the LFN(s) with the appropriate metadata values.

                @param queryList: list of instances of the LDRMetadataQuery class, each describing
                                  a query to be done
                        
                @return: list of strings representing the union of all LFN(s) found from the queries
                """

                # check arguments
                if not isinstance(queryList, list):
                        msg = "Argument must be a list of instances of LDRMetadataQuery"
                        raise LDRdataFindClientException, msg
                for query in queryList:
                        if not isinstance(query, LDRMetadataQuery):
                                msg = "Argument must be an instance of LDRMetadataQuery"
                                raise LDRdataFindClientException, msg

                # prepare the messange to send down the socket
                msg = "METALFN\0"
                for q in queryList:
                        msg += "%s" % str(q)

                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for LFNs with metadata query %s : %s" % (sql, str(output[0]))
                        raise LDRdataFindClientException, msg


                return output
                
        def pfnQueryWithMetadata(self, queryList):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values.

                @param queryList: list of instances of the LDRMetadataQuery class, each describing
                                  a query to be done

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if not isinstance(queryList, list):
                        msg = "Argument must be a list of instances of LDRMetadataQuery"
                        raise LDRdataFindClientException, msg
                for query in queryList:
                        if not isinstance(query, LDRMetadataQuery):
                                msg = "Argument must be an instance of LDRMetadataQuery"
                                raise LDRdataFindClientException, msg

                msg = "METAPFN\0" 
                for q in queryList:
                        msg += "%s" % str(q)

                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs with metadata query %s : %s" % (sql, str(output[0]))
                        raise LDRdataFindClientException, msg

                return output

        def pfnQueryWithMetadataRegExp(self, sql, rexp, offset=None, number=None):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values and return those that match a regular expression.

                @param sql: clause that will be part of a SQL query done to find the LFN(s)
                        with the appropriate metadata values. A an example would be
                        
                        gpsStart <= '777777777' AND gpsEnd >= '666666666' AND instrument = 'H' AND runTag = 'S2'

                @param rexp: regular expression against which to match found PFNs

                @param offset: the offset into the list of matching PFNs at which to begin 
                        returning results

                @param number: the total number of PFNs to return

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if type(sql) != types.StringType:
                        msg = "Argument 'sql' must be a string"
                        raise LDRdataFindClientException, msg

                if type(rexp) != types.StringType:
                        msg = "Argument 'rexp' must be a string"
                        raise LDRdataFindClientException, msg

                if offset:
                    if type(offset) != types.IntType:
                            msg = "Argument 'offset' must be a positive integer or zero"
                            raise LDRdataFindClientException, msg
                
                    if offset < 0:
                            msg = "Argument 'offset' must be a positive integer or zero"
                            raise LDRdataFindClientException, msg
                
                if number:
                    if type(number) != types.IntType:
                            msg = "Argument 'number' must be a positive integer"
                            raise LDRdataFindClientException, msg
                
                    if number <= 0:
                            msg = "Argument 'number' must be a positive integer"
                            raise LDRdataFindClientException, msg

                msg = "METAREPFN\0%s\0%s\0" % (str(sql), str(rexp))
                if offset:
                        msg += "%d\0" % offset
                else:
                        msg += "0\0"
                if number:
                        msg += "%d\0" % number

                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs with metadata query %s : %s" % (sql, str(output[0]))
                        raise LDRdataFindClientException, msg

                

                return output

        def pfnQueryWithMetadataUnion(self, sql1, sql2):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values found by the union of two queries.

                @param sql1: clause that will be part of a SQL query done to find the LFN(s)
                        with the appropriate metadata values. A an example would be
                        
                        gpsStart >= '777777777' AND gpsStart <= '888888888' AND instrument = 'H' AND runTag = 'S2'

                @param sql2: second clause with similar form as above. An example would be

                        gpsEnd >= '777777777' AND gpsEnd <= '888888888' AND instrument = 'H' AND runTag = 'S2'

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if type(sql1) != types.StringType:
                        msg = "Argument 'sql1' must be a string"
                        raise LDRdataFindClientException, msg

                if type(sql2) != types.StringType:
                        msg = "Argument 'sql2' must be a string"
                        raise LDRdataFindClientException, msg

                msg = "METAPFNUNION\0%s\0%s\0" % (str(sql1), str(sql2))

                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs with metadata queryies %s AND %s : %s" % (sql1, sql2, str(output[0]))
                        raise LDRdataFindClientException, msg

                

                return output
