"""
The LSCsegFind.py module contains Classes e.g., used by LSCsegFind
"""

__author__ = "Greg Mendell: module that contains Classes e.g., used by LSCsegFind"
__date__ = '$Date$'
__version__ = '$Revision$'[0:0]

import sys
import os
import exceptions
import types

#from pyGlobus import io
#from pyGlobus import security

try:
        from LSCsegments import *
except ImportError:
        print >>sys.stderr, "\nPython LSCsegments.py library not found; please add its location to the PYTHONPATH environmental variable.\n"
        print >>sys.stderr, "For example, from the bash or csh shell respectively, run:\n"
        print >>sys.stderr, "    export PYTHONPATH=$PYTHONPATH:$LSCSOFT/src/glue/lib"
        print >>sys.stderr, "or "
        print >>sys.stderr, "    setenv PYTHONPATH ${PYTHONPATH}:$LSCSOFT/src/glue/lib"
        print >>sys.stderr, "\nwhere $LSCSOFT needs to be replaced with the local directory with software from the lscsoft cvs repository."
        sys.exit(1)

class LSCsegFindException(exceptions.Exception):
        """
        Exceptions raised by the classes and methods in this client
        will be instances of this class.
        """
        def __init__(self, args=None):
                """
                Create an instance of this class, ie. an LSCsegFindClientException.

                @param args:

                @return: Instance of class LSCsegFindClientException
                """
                print 'An exception was raised in Class LSCsegFind:'
                self.args = args
                

class LSCsegFind(object):
        """
        Class that finds GPS segment data.
        """
        
        # Class data:
        
        # tuple of allowed sources of GPS segments
        segSourceTypes = ('none','url')
        
        # tuple of allowed interferometers
        ifos = ('H1', 'H2', 'L1')

        # tuple of allowed segment types; this is currently unused
        # segTypes = ('none', 'SCIENCE_MODE')
        
        def __init__(self, segSourceType=None, ifo=None, startTime=None, endTime=None,\
        segURL=None, segType=None, cfgFile=None, minLength=None, addStart=None, reduceEnd=None,\
        outputFormat=None, strict=None, saveFiles=None):
                """
                  See help for LSCsegFind for definitions of parameters.
                """
                                
                # get defaults from local config file
                self.urlDefaultsExist = False
                self.urlDefaultName = 'urlDefault_%s' % ifo
                if cfgFile != None:
                   try:
                       exec ( "from %s import %s" % (cfgFile, self.urlDefaultName) )
                       self.urlDefaultsExist = True
                   except ImportError:
                       # config could be a local file; try execfile:
                       try:
                           execfile(cfgFile)
                           self.urlDefaultsExist = True                           
                       except:
                           # Will handle this error below only when values from config file are needed.
                           pass
                
                # Check that required arguments are set
                if segSourceType == None:
                      msg = "\nMissing or bad input data: 'segSourceType' not found"
                      raise LSCsegFindException, msg
                
                if ifo == None:
                      msg = "\nMissing or bad input data: 'ifo' not found"
                      raise LSCsegFindException, msg
                         
                if startTime == None:
                      msg = "\nMissing or bad input data: 'startTime' not found"
                      raise LSCsegFindException, msg
                        
                if endTime == None:
                      msg = "\nMissing or bad input data: 'endTime' not found"
                      raise LSCsegFindException, msg
                
                if segURL == None:
                      if self.urlDefaultsExist:
                           try:
                                exec ( "segURL = %s" % self.urlDefaultName)
                           except:
                                msg = "\nMissing or bad input data: a default for 'segURL' was not found in config file %s for ifo %s" % (cfgFile, ifo)
                                raise LSCsegFindException, msg                           
                      else:
                           if cfgFile != None:
                                msg = "\nMissing or bad input data: 'segURL' not found; the config file %s does not exist or is corrupt" % cfgFile
                                raise LSCsegFindException, msg
                           else:
                                msg = "\nMissing or bad input data: 'segURL' not found"
                                raise LSCsegFindException, msg

                if minLength == None:
                      msg = "\nMissing or bad input data: 'minLength' not found"
                      raise LSCsegFindException, msg
                
                if addStart == None:
                      msg = "\nMissing or bad input data: 'addStart' not found"
                      raise LSCsegFindException, msg
                
                if reduceEnd == None:
                      msg = "\nMissing or bad input data: 'reduceEnd' not found"
                      raise LSCsegFindException, msg
                
                if outputFormat == None:
                      msg = "\nMissing or bad input data: 'outputFormat' not found"
                      raise LSCsegFindException, msg
                      
                if strict == None:
                      msg = "\nMissing or bad input data: 'strict' not found"
                      raise LSCsegFindException, msg
                
                if saveFiles == None:
                      msg = "\nMissing or bad input data: 'saveFiles' not found"
                      raise LSCsegFindException, msg
                      
                # Check arguments are arguments are of the required type and have valid values
                            
                if type(segSourceType) != types.StringType:
                        msg = "\nMissing or bad input data: 'segSourceType' must be a string"
                        raise LSCsegFindException, msg

                if not (segSourceType in LSCsegFind.segSourceTypes):
                        msg = "\nMissing or bad input data: 'segSourceType' = '%s' is invalid; current valid values are %s" % (segSourceType, LSCsegFind.segSourceTypes)
                        raise LSCsegFindException, msg

                if type(ifo) != types.StringType:
                        msg = "\nMissing or bad input data: 'ifo' must be a string"
                        raise LSCsegFindException, msg                        

                if not (ifo in LSCsegFind.ifos):
                        msg = "\nMissing or bad input data: 'ifo' = '%s' is invalid; current valid values are %s" % (ifo, LSCsegFind.ifos)
                        raise LSCsegFindException, msg
                                                                              
                try:
                        self.__check_gps(str(startTime))
                except Exception, e:
                        msg = "\nMissing or bad input data: GPS 'startTime' must be positive integer at least 9 digits long"
                        raise LSCsegFindException, msg

                try:
                        self.__check_gps(str(endTime))
                except Exception, e:
                        msg = "\nMissing or bad input data: GPS 'endTime' must be positive integer at least 9 digits long"
                        raise LSCsegFindException, msg
                        
                if long(startTime) >= long(endTime):
                        msg = "\nMissing or bad input data: 'startTime' must be less than argument 'endTime'"
                        raise LSCsegFindException, msg
                        

                if type(segURL) != types.StringType:
                        msg = "\nMissing or bad input data: 'segURL' must be a string"
                        raise LSCsegFindException, msg
                
                # TO DO: test URL
                
                if type(minLength) != types.LongType:
                        msg = "\nMissing or bad input data: 'minLength' must be an integer "
                        raise LSCsegFindException, msg
                
                if type(minLength) < 0L:
                        msg = "\nMissing or bad input data: 'minLength' must be a positive integer "
                        raise LSCsegFindException, msg

                if type(addStart) != types.LongType:
                        msg = "\nMissing or bad input data: 'addStart' must be an integer "
                        raise LSCsegFindException, msg
                
                if type(addStart) < 0L:
                        msg = "\nMissing or bad input data: 'addStart' must be a positive integer "
                        raise LSCsegFindException, msg
                        
                if type(reduceEnd) != types.LongType:
                        msg = "\nMissing or bad input data: 'reduceEnd' must be an integer "
                        raise LSCsegFindException, msg
                
                if type(reduceEnd) < 0L:
                        msg = "\nMissing or bad input data: 'reduceEnd' must be a positive integer "
                        raise LSCsegFindException, msg
                
                if type(outputFormat) != types.StringType:
                        msg = "\nMissing or bad input data: 'outputFormat' must be a string"
                        raise LSCsegFindException, msg
                                    
                try:
                        if (strict != True) and (strict != False):
                           msg = "\nMissing or bad input data: 'strict' must be True or False"
                           raise LSCsegFindException, msg
                except Exception, e:
                        msg = "\nMissing or bad input data: 'strict' must be True or False"
                        raise LSCsegFindException, msg

                try:
                        if (saveFiles != True) and (saveFiles != False):
                           msg = "\nMissing or bad input data: 'saveFiles' must be True or False"
                           raise LSCsegFindException, msg
                except Exception, e:
                        msg = "\nMissing or bad input data: 'saveFiles' must be True or False"
                        raise LSCsegFindException, msg
                        

                self.segSourceType = segSourceType
                self.ifo = ifo
                self.startTime = long(startTime)
                self.endTime = long(endTime)
                self.segURL = segURL
                self.segType = segType
                self.minLength = long(minLength)
                self.addStart = long(addStart)
                self.reduceEnd = long(reduceEnd)
                self.outputFormat = outputFormat
                self.strict = strict
                self.saveFiles = saveFiles

                self.myScienceData = ScienceData()
                self.mySegList = []
                
        def __check_gps(self, gpsString):
                """
                Minimal checking on GPS time strings. Raises a LSCsegFindException if
                the GPS time string is not at least 9 digits long.

                @param gpsString: The string representing the 9+ digit GPS time.

                @returns: None
                """
                if len(gpsString) < 9:
                        msg = "GPS time must be at least 9 digits"
                        raise LSCsegFindException, msg

                try:
                        a = long(gpsString)
                except Exception, e:
                        msg = "GPS time must be an integer"
                        raise LSCsegFindException, msg
                
                if a < 0:
                        msg = "GPS time must be a positive integer"
                        raise LSCsegFindException, msg

        def GetSegments(self):
                
                fileName = 'tmpLSCsegFindWebPage.txt'

                try:
                    webPageString = os.system('/usr/bin/curl %s 1> %s 2>/dev/null' % (self.segURL, fileName))
                except:
                    msg = "Could not retrieve web page '%s', make sure /usr/bin/curl and web page exists." % self.segURL
                    raise LSCsegFindException, msg
  
                if webPageString:
                    msg = "Error retrieving web page '%s', make sure web page exists." % self.segURL
                    raise LSCsegFindException, msg
                
                try: 
                    self.myScienceData.read(fileName,self.minLength)
                except:
                    msg = "Could not read file %s" % fileName
                    raise LSCsegFindException, msg
                                                    
                formatList = self.outputFormat.split()
                for i in range (0,self.myScienceData.__len__()):
                    segment = self.myScienceData.__getitem__(i)
                    segID = segment.id()
                    segStart = long(segment.start())
                    segEnd = long(segment.end())
                    segDuration = segment.dur()
                    
                    if  (segStart >= self.endTime) or (segEnd <= self.startTime):
                         # This segment is not in the interval [self.startTime,self.endTime)
                         continue
                    else:
                         if self.strict and (segEnd > self.endTime):
                            # Truncate this segment at the end of the requested interval
                            segEnd = self.endTime
                         if self.strict and (segStart < self.startTime):
                            # Truncate this segment at the beginning of the requested interval
                            segStart = self.startTime
                    if self.addStart > 0L:
                            segStart += self.addStart
                    if self.reduceEnd > 0L:                            
                            segEnd -= self.reduceEnd
                    if (segEnd - segStart) < self.minLength:
                         # This segment is too short
                         continue
                    
                    # Print out the segment
                    
                    if (self.outputFormat == 'PYTHON_LIST') or (self.outputFormat == 'TCL_LIST'):
                         # Just build up the list for now and print it below
                         self.mySegList.append([segStart, segEnd])                    
                    elif (self.outputFormat == 'SCIENCE_SEGMENTS'):
                         # prints segments as objects of the LSC python Class ScienceSegment
                         print ScienceSegment(tuple([segID,segStart,segEnd,segDuration]))                    
                    else:
                         # use format to print                         
                         outputDictionary = {'%i':segID, '%s':segStart, '%e':segEnd, '%d':segDuration}
                         for formatChar in formatList:
                             try:
                                 print outputDictionary[formatChar],
                             except:
                                 print formatChar,
                         print ''                       
                # End for i in range (0,self.myScienceData.__len__())
                
                if (self.outputFormat == 'PYTHON_LIST'):                
                     print self.mySegList
                elif (self.outputFormat == 'TCL_LIST'):                
                     print "{",
                     for thisSegment in self.mySegList:
                         print "{ %s %s }" % (thisSegment[0],thisSegment[1]),
                     print "}"                      
                
                # clean up                
                # TODO: handle saving and reusing of files; for now just keep the temporary file if saveFiles is true
                if self.saveFiles:
                     savedFileName = 'LSCsegFindWebPage.txt'                
                     try:
                         rmOut = os.system('/bin/mv %s %s 1>/dev/null 2>/dev/null' % (fileName, savedFileName))
                     except:
                         # Just leave the file on disk for now; will not hurt anything.
                         pass
                else:

                     try:
                         rmOut = os.system('/bin/rm -f %s 1>/dev/null 2>/dev/null' % fileName)
                     except:
                         # Just leave the file on disk for now; will not hurt anything.
                         pass
                    
        # End def GetSegments(self)
        
# End class LSCsegFind(object)

#def checkCredentials():
#        """
#        Check to make sure that the proper Grid Credentials (a proxy certificate) is
#        available in order to authenticate to the remote LDRsegFindServer.
#        """
#        # verify that we have access to credentials
#        try:
#                proxyText = security.grid_proxy_info()
#        except Exception, e:
#                print >>sys.stderr, "Error verifying credentials: %s" % e
#                print >>sys.stderr, "Run 'grid-proxy-init' to generate a proxy certificate"
#                sys.exit(1)
#
#        pat = re.compile(r'timeleft : (\d{1,3}):(\d\d):(\d\d)')
#
#        try:
#                if isinstance(proxyText, str):
#                        m = pat.search(proxyText)
#                elif isinstance(proxyText, tuple):
#                        m = pat.search(proxyText[0])
#                else:
#                        raise RuntimeError, "bad format for proxyText in checkCredentials"
#                        
#                hours, minutes, seconds = map(int, m.groups())
#        except Exception, e:
#                print >>sys.stderr, "Error parsing proxy information: %s" % e
#                sys.exit(1)
#
#        timeleft = seconds + 60 * minutes + 3600 * hours
#
#        if timeleft < 300:
#                print >>sys.stderr, "Less than 5 minutes left for proxy certificate."
#                print >>sys.stderr, "Run 'grid-proxy-init' to generate a new proxy certificate"
#                sys.exit(1)
