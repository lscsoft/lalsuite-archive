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
                if cfgFile != None:
                   self.cfgFile = self.__get_STR('cfgFile',cfgFile)
                   try:
                       #exec ( "from %s import %s" % (cfgFile, self.urlDefaultName) )
                       exec ( "from %s import urlDefaults" % self.cfgFile )
                       self.urlDefaultsExist = True
                   except ImportError:
                       # config could be a local file; try execfile:
                       try:
                           execfile(self.cfgFile)
                           self.urlDefaultsExist = True
                       except:
                           # Will handle this error below only when values from config file are needed.
                           pass
                else:
                   self.cfgFile = cfgFile
                   
                # Check that required parameters are set
                
                self.segSourceType = self.__get_STR('segSourceType',segSourceType)
                if not (segSourceType in LSCsegFind.segSourceTypes):
                        msg = "\nMissing or bad input data: 'segSourceType' = '%s' is invalid; current valid values are %s" % (segSourceType, LSCsegFind.segSourceTypes)
                        raise LSCsegFindException, msg
                
                self.ifo = self.__get_STR('ifo',ifo)
                
                self.startTime = self.__get_GPSLONG('startTime',startTime)
                self.endTime = self.__get_GPSLONG('endTime',endTime)
                if self.startTime >= self.endTime:
                      msg = "\nMissing or bad input data: 'startTime' must be less than argument 'endTime'"
                      raise LSCsegFindException, msg
                                         
                if segURL == None:
                      # if no segURL input, then try to get this from the urlDefaults dictionary set in the config file.
                      if self.urlDefaultsExist:
                           try:
                                segURL = urlDefaults[self.ifo]
                           except:
                                try:
                                     msg = "\nMissing or bad input data: a default URL for ifo '%s' was not found in config file '%s'; defaults exist for these ifos: %s" % (ifo, self.cfgFile, urlDefaults.keys())
                                except:
                                     msg = "\nMissing or bad input data: a default URL for ifo '%s' was not found in config file '%s'; no defaults were found." % (ifo, self.cfgFile)
                                raise LSCsegFindException, msg                           
                      else:
                           if self.cfgFile != None:
                                msg = "\nMissing or bad input data: 'segURL' not found; the config file %s does not exist or is corrupt" % self.cfgFile
                                raise LSCsegFindException, msg
                           else:
                                msg = "\nMissing or bad input data: 'segURL' not found"
                                raise LSCsegFindException, msg
                self.segURL = self.__get_STR('segURL',segURL)
                
                # segType is currently unused, so no testing:
                self.segType = segType
                
                self.minLength = self.__get_PLONG('minLength',minLength)
                self.addStart = self.__get_ULONG('addStart',addStart)
                self.reduceEnd = self.__get_ULONG('reduceEnd',reduceEnd)
                self.outputFormat = self.__get_STR('outputFormat',outputFormat)
                self.strict = self.__get_BOOLEAN('strict',strict)
                self.saveFiles = self.__get_BOOLEAN('saveFiles',saveFiles)
                
                # data not from parameters
                self.myScienceData = ScienceData()
                self.mySegList = []
                           
        def __check_exists(self, inputName, inputValue):
                """
                   Check that a value exists, else raise an exception
                """        
                if inputValue == None:
                      msg = "\nMissing or bad input data: '%s' not found" % inputName
                      raise LSCsegFindException, msg
                                                           
        def __get_STR(self, inputName, inputValue):
                """
                   Convert inputValue to type str; raise exceptions if value does not exist or conversion fails
                """        
                self.__check_exists(inputName, inputValue)
                try:
                   testValue = str(inputValue)
                   return testValue
                except: 
                   msg = "\nMissing or bad input data: '%s' must be a string" % inputName
                   raise LSCsegFindException, msg
                        
        def __get_ULONG(self, inputName, inputValue):
                """
                   Convert inputValue to nonnegative long integer; raise exceptions if value does not exist or conversion fails
                """        
                self.__check_exists(inputName, inputValue)
                msg = "\nMissing or bad input data: '%s' must be integer >= 0; value given was %s" % ( str(inputName), str(inputValue) )
                try:
                   testValue = long(inputValue)
                   if testValue < 0L:
                      raise LSCsegFindException, msg                   
                   else:
                       return testValue
                except:
                   raise LSCsegFindException, msg
        
        def __get_PLONG(self, inputName, inputValue):
                """
                   Convert inputValue to nonnegative long integer; raise exceptions if value does not exist or conversion fails
                """        
                self.__check_exists(inputName, inputValue)
                msg = "\nMissing or bad input data: '%s' must be integer > 0; value given was %s" % ( str(inputName), str(inputValue) )
                try:
                   testValue = long(inputValue)
                   if testValue <= 0L:
                      raise LSCsegFindException, msg                   
                   else:
                       return testValue
                except:
                   raise LSCsegFindException, msg
        
        def __get_GPSLONG(self, inputName, inputValue):
                """
                   Convert inputValue to positive long integer with at least 9 digits.
                   Raise exceptions if the value does not exist or conversion fails                
                """
                self.__check_exists(inputName, inputValue)
                msg = "\nMissing or bad input data: '%s' must be a valid GPS time with at least 9 digits; value given was %s" % ( str(inputName), str(inputValue) )
                gpsString = str(inputValue)
                if len(gpsString) < 9:
                        raise LSCsegFindException, msg
                try:
                        testValue = long(inputValue)
                except Exception, e:
                        raise LSCsegFindException, msg
                
                if testValue < 0L:
                        raise LSCsegFindException, msg
                else:
                        return testValue
                
        def __get_BOOLEAN(self, inputName, inputValue):
                """
                   Convert inputValue to nonnegative long integer; raise exceptions if value does not exist or conversion fails
                """        
                self.__check_exists(inputName, inputValue)
                msg = "\nMissing or bad input data: '%s' must be True or False; value given was %s" % ( str(inputName), str(inputValue) )
                try:
                        if (inputValue != True) and (inputValue != False):
                           raise LSCsegFindException, msg
                        else:
                           return inputValue
                except:
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
                    segDuration = long(segment.dur())
                    
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
                    
                    # update value and check against minLength
                    segDuration = segEnd - segStart                  
                    if segDuration < self.minLength:
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
