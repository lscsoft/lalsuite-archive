"""
The LSCsegFind.py module contains Classes e.g., used by LSCsegFind
"""

__author__ = "Greg Mendell: module that contains Classes e.g., used by LSCsegFind"
__date__ = '$Date$'
__version__ = '$Revision$'[0:0]

import sys
import os
import exceptions

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
                print 'An exception was raised and handled in Class LSCsegFind. See last message:'
                self.args = args
                

class LSCsegFind(object):
        """
        Class that finds GPS segment data.
        """
        
        # Class data:
          
        
        def __init__(self, server=None, startTime=None, endTime=None, segType=None,\
        minLength=None, addStart=None, reduceEnd=None, outputFormat=None, strict=None,\
        saveFiles=False, cfgFile=None, localSegFile=None, coalesce=False, showTypes=False):
                """
                  See help for LSCsegFind for definitions of parameters.
                """

                # Define defaults for local temporary files and local files to save if saveFiles is True.
                self.tmpCfgDefaultsFileName  = 'tmpCfgDefaults_LSCsegFind.py'
                self.tmpCfgDefaultsFileNameC = 'tmpCfgDefaults_LSCsegFind.pyc'
                self.savCfgDefaultsFileName  = 'savCfgDefaults_LSCsegFind.py'
                self.tmpSegDataFileName = 'tmpSegData_LSCsegFind.txt'
                self.savSegDataFileName = 'savSegData_LSCsegFind.txt'
                self.saveFiles = saveFiles
                self.saveFiles = self.__get_BOOLEAN('saveFiles',saveFiles)                
                self.showTypes = self.__get_BOOLEAN('showTypes',showTypes)
                                
                if localSegFile == None:
                   # Usual case; use server or config file to find segments.
                   # get defaults from server or local config file
                   if cfgFile == None:
                      # Usual case: get the default config file from a web server; save as temporary file
                      # Set up the server host name 
                      self.server = self.__get_STR('server',server)
                      if self.server.find(':') < 0:
                          # no port specified
                          self.host = self.server
                          self.port = 0L
                      else:
                          # server and port specified
                          hostString, portString = self.server.split(':')
                          self.host = self.__get_STR('host',hostString)
                          self.port = self.__get_PLONG('minLength',portString)

                      # set up serverURL (note that the file name on the server is hard-coded here.)                                      
                      self.serverFileName = "lscsegfind/LSCsegFindCfgDefaults.txt"
                      self.serverURL = "%s/%s" % (self.host, self.serverFileName)
                   
                      self.cfgFile = self.tmpCfgDefaultsFileName                
                      self.__getWebPage(self.serverURL,self.cfgFile)
                   else: 
                      self.cfgFile = self.__get_STR('cfgFile',cfgFile)
                
                   self.cfgDefaultsExist = False
                   try:
                       #exec ( "from %s import %s" % (cfgFile, self.urlDefaultName) )
                       exec ( "from %s import cfgDefaults" % self.cfgFile )
                       self.cfgDefaultsExist = True
                   except ImportError:
                       # config file could be a local file; try execfile:
                       try:
                          execfile(self.cfgFile)
                          self.cfgDefaultsExist = True
                       except:
                          # Will handle this error below only when values from config file are needed.
                          pass

                   if self.showTypes:
                       # just display the default segTypes and return
                       print "\n Available types : description "
                       setTypeList = cfgDefaults.keys()
                       setTypeList.sort()
                       for segType in setTypeList:
                               print "%16s : %s  " % (segType,  cfgDefaults[segType]['desc'])
                       return


                   # Try to get segURL from the cfgDefaults dictionary set in the config file.
                   self.segType = self.__get_STR('type',segType)                   
                   segURL = None
                   if self.cfgDefaultsExist:
                        try:
                            segURL = cfgDefaults[self.segType]['url']
                            coalesce = cfgDefaults[self.segType]['coalesce']                                
                        except:
                            try:
                               msg = "\nMissing or bad input data: the type '%s' was not found; found these types: %s" % (self.segType, cfgDefaults.keys())
                            except:
                               msg = "\nMissing or bad input data: the type '%s' was not found; failed to find any defined types" % self.segType
                               raise LSCsegFindException, msg                           
                   else:
                        if self.cfgFile != None:
                            msg = "\nMissing or bad input data: the type '%s' was not found." % self.segType
                            raise LSCsegFindException, msg
                        else:
                            msg = "\nMissing or bad input data: invalid server or config file"
                            raise LSCsegFindException, msg
                   # Will get the segment file from a URL; set the local file to the tmp file name
                   self.segURL = self.__get_STR('segURL',segURL)
                   self.localSegFile = self.tmpSegDataFileName                   
                else:   
                   # Will use the localSegFile parameter
                   self.segURL = None                   
                   self.localSegFile = self.__get_STR('localSegFile',localSegFile)                   
                   if self.showTypes:
                      msg = "\nMissing or bad input data:: cannot show types when localSegFile is given"
                      raise LSCsegFindException, msg
                # End if localSegFile == None
                
                # Set up remaining parameters
                self.coalesce = self.__get_BOOLEAN('coalesce',coalesce)
                self.startTime = self.__get_GPSLONG('startTime',startTime)
                self.endTime = self.__get_GPSLONG('endTime',endTime)
                if self.startTime >= self.endTime:
                      msg = "\nMissing or bad input data: 'startTime' must be less than argument 'endTime'"
                      raise LSCsegFindException, msg
                self.minLength = self.__get_PLONG('minLength',minLength)
                self.addStart = self.__get_ULONG('addStart',addStart)
                self.reduceEnd = self.__get_ULONG('reduceEnd',reduceEnd)
                self.outputFormat = self.__get_STR('outputFormat',outputFormat)
                self.strict = self.__get_BOOLEAN('strict',strict)
                self.saveFiles = self.__get_BOOLEAN('saveFiles',saveFiles)
                
                # data not from parameters
                self.myScienceData = ScienceData()
                self.mySegList = []

        def __del__(self):
                """
                   Clean up and remove temporary files or save if requested
                """ 
                
                # Always remove the tmp .pyc cfg file
                try:
                    rmOut = os.system('/bin/rm -f %s 1>/dev/null 2>/dev/null' % self.tmpCfgDefaultsFileNameC)

                except:
                    # Just leave the file on disk for now; will not hurt anything.
                    pass                         
                
                # Remove tmp files or mv tmp files to sav files
                if self.saveFiles:

                     # Append the segtype name to the beginning the name of a segment file to save.
                     try:
                         tmpName = '%s_%s' % (self.segType, self.savSegDataFileName)
                         self.savSegDataFileName = tmpName
                     except:
                         pass
                     
                     try:
                         mvOut = os.system('/bin/mv %s %s 1>/dev/null 2>/dev/null' % (self.tmpCfgDefaultsFileName, self.savCfgDefaultsFileName))
                         mvOut = os.system('/bin/mv %s %s 1>/dev/null 2>/dev/null' % (self.tmpSegDataFileName, self.savSegDataFileName))
                     except:
                         # Just leave the file on disk for now; will not hurt anything.
                         pass
                else:
                     try:
                         rmOut = os.system('/bin/rm -f %s 1>/dev/null 2>/dev/null' % self.tmpCfgDefaultsFileName)
                         rmOut = os.system('/bin/rm -f %s 1>/dev/null 2>/dev/null' % self.tmpSegDataFileName)                         
                     except:
                         # Just leave the file on disk for now; will not hurt anything.
                         pass                         

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
                   Convert inputValue to positive long integer; raise exceptions if value does not exist or conversion fails
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

        def __getWebPage(self,urlString,localFileName):
                """
                   returns a url to a local file.
                """
                try:
                    for i in range (0,2):
                        curlExit = os.system('/usr/bin/curl --fail --connect-timeout 100 %s 1> %s 2>/dev/null' % (urlString, localFileName))
                        if long(curlExit) == 0L:
                             curlFailed = False
                             break
                    if curlFailed:
                        msg = "Error retrieving web page '%s', make sure web page exists; error was %s" % (urlString, str(curlExit))
                        raise LSCsegFindException, msg
                except:
                    msg = "Could not retrieve web page '%s', make sure /usr/bin/curl and web page exist." % urlString
                    raise LSCsegFindException, msg
        
        def GetSegments(self):
                """
                   Get segments and print those that match LSCsegFind input parameters to stdout
                """
                
                if self.segURL != None:
                   # Get the segments from a web page; save them to a local file.                
                   self.__getWebPage(self.segURL,self.tmpSegDataFileName)

                # Read in the segments in the ScienceData Class
                try: 
                    if self.coalesce:
                         # Get all the segments; will coalesce segments below
                         self.myScienceData.read(self.localSegFile,1)
                    else:
                         self.myScienceData.read(self.localSegFile,self.minLength)
                except:
                    msg = "Could not read file %s" % self.localSegFile
                    raise LSCsegFindException, msg

                if self.coalesce:
                    try:                 
                        numSciSegs = self.myScienceData.coalesce()
                    except:
                        msg = "An error occured trying to coalesce the data." 
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
