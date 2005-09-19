"""
LSCfileAdd.py module contains the classes for the LSCfileAdd utility
"""

__author__="Ben Johnson: module contains the classes for the LSCfileAdd utility"
__date__='$Date$'
__version__='$Revision$'[0:0]

import os
import sys
import getopt
import exceptions
import md5
import re
import urlparse
from types import *
#from pyGlobus import security

## LDR specific(?) modules
#try:
#        ldrHome = os.environ["LDR_LOCATION"]
#except:
#        sys.stderr.write("LDR_LOCATION environment variable undefined\n")
#        sys.exit(1)

#sys.path.append(os.path.join(ldrHome, "ldr/lib"))

import LDRUtil
import RLS
import LDRMetadataCatalog

class LSCfileAddMetadataAttr(object):
        """
        This class sets up a dictionary containing the LDR metadata fields
        available to this routine, and defines how they are used in the
        context of this script.
        
        The dictionary is layed out as follows:
                FieldName:DescriptionDict
                
                FieldName is the metadata field, e.g. "size" or "md5"
                
                DescriptionList contains the following:
                       0) Value    # field value which is put into database
                       1) Type     # data type, must be acceptable by database, e.g. int, string, float
                       2) Default  # default if any
                       3) UserSet  # If this can be set via the program's user, True or False
                       4) Cli_arg_long # long version of cli argument with NO DASHES, e.g. "gpsStart"
                       5) Cli_arg_short # short version of cli argument with NO DASHES, e.g. "s"
                       6) Description # description of this field
                       7) Test_method # Function pointers, baby!
                                      # method used to test the validity of this parameter
                                      # format is 'Test_method':getattr(self,"FUNCTION_NAME")
                                      
        
        """
        attr = {}
        
        def gps_processor(self):
                """
                Performs a consistancy check on start and end gps times
                along with the duration parameter. 
                
                If the three fields are not None, makes sure they are consistent.

                If two are specified, makes sure they are sane (e.g. start time < end time),
                and if so, calculates third field.
                
                If only one field is specified, returns an error?
                
                Returns 0 if there are no errors, otherwise a string containing the error message.
                """
                # set some "booleans"
                if self.attr['gpsStart']['Value']:
                        start = 1
                else:
                        start = 0
                if self.attr['gpsEnd']['Value']:
                        end = 1
                else:
                        end = 0
                if self.attr['duration']['Value']:
                        duration = 1
                else:
                        duration = 0
                
                ## Now, according to what's been set perform various tests
                
                # make sure at least two parameters are set
                if start + end + duration < 2:
                        msg = "No fields specified!"
                        if start:
                                msg = "Only know the start time."
                        if end:
                                msg = "Only know the end time."
                        if duration:
                                msg = "Only know the data duration."
                        return "Not enough information to fill gps start, end, and duration fields. %s" % (msg,)
                # now perform basic "is this a valid gps time at all?" checks
                if start:
                        inputValue = self.attr['gpsStart']['Value']
                        gpsString = str(inputValue)
                        
                        try:
                                testValue = long(inputValue)
                        except Exception, e:
                                return "GPS start time must be an integer. Received \"%s\"" % (str(inputValue),)
                        if testValue < 0L:
                                return "GPS start time must be positive. Received \"%s\"" % (str(inputValue),)
                        if testValue < 100000000:
                                return "GPS start time must be at least nine digits long. Received \"%s\"" % (str(inputValue),)
                if end:
                        inputValue = self.attr['gpsEnd']['Value']
                        gpsString = str(inputValue)
                        
                        try:
                                testValue = long(inputValue)
                        except Exception, e:
                                return "GPS end time must be an integer. Received \"%s\"" % (str(inputValue),)
                        if testValue < 0L:
                                return "GPS end time must be positive. Received \"%s\"" % (str(inputValue),)
                        if testValue < 100000000:
                                return "GPS end time must be at least nine digits long. Received \"%s\"" % (str(inputValue),)
                # see if duration is an integer etc.
                if duration:
                        inputValue = self.attr['duration']['Value']
                        durString = str(inputValue)
                        try:
                                testValue = long(inputValue)
                        except Exception, e:
                                return "Duration time must be an integer. Received \"%s\"" % (str(inputValue),)
                # now perform consistancy checks, and missing parameter calculations
                starttime = self.attr['gpsStart']['Value']
                endtime = self.attr['gpsEnd']['Value']
                durtime = self.attr['duration']['Value']
                if start and duration:
                        self.attr['gpsEnd']['Value'] = starttime + durtime
                        return 0
                if end and duration:
                        self.attr['gpsStart']['Value'] = endtime - durtime
                        return 0
                if start and end:
                        if endtime <= starttime:
                                return "GPS end time must be greater than start time. Received gps-start = %s, gps-end = %s" % (str(starttime), str(endtime))
                        if duration:
                                durtime = self.attr['duration']['Value']
                                if durtime != (endtime - starttime):
                                        return "Duration parameter does not match GPS start and end times. \
                                        Recevied gps-start = %s, gps-end = %s, duration = %s" (str(starttime),str(endtime),str(durtime))
                        else:
                                self.attr['duration']['Value'] = endtime - starttime
                                return 0
        ## END gps_processor()
        
        def ifo_processor(self):
                """
                Checks given interferometer against list of accepted interferometers.
                returns 0 if found, an error string if not.
                """
                ifo = self.attr['interferometer']['Value']
                accepted_ifos = self.attr['interferometer']['Accepted_values']
                if not ifo:
                        return "Interferomter(s) must be specified. Accepted values are %s or any combination thereof." % (str(accepted_ifos),)
                # now split ifos into \w\d+ groups, though I don't know how to do this with regular expressions
                templist = list(ifo)
                ifolist = []
                for item in templist:
                        if item.isalpha():
                                ifolist.append(item)
                                continue
                        if item.isdigit():
                                ifolist[len(ifolist) - 1] += item
                for myifo in ifolist:
                        if not accepted_ifos.count(myifo):
                                return "Interferometer \"%s\" was not found in the list of accepted interferometers. Valid ifos are any combination of the following: %s." % (ifo, str(accepted_ifos))
                        if accepted_ifos.count(myifo) > 1:
                                return "An interferometer value can only be specified once. At least two of ifo \"%s\" were specified." % (ifo,)
                # all is well, then
                return 0
        ## END ifo_processor()
        
        def site_processor(self):
                """
                Check given site against list of accepted sites.
                returns 0 if found, an error string if not
                """
                site = self.attr['site']['Value']
                if not site:
                        return site
                # split site list
                sitelist = list(site)
                accepted_sites = self.attr['site']['Accepted_values']
                for mysite in sitelist:
                        if not accepted_sites.count(site):
                                return "Site \"%s\" was not found in the list of accepted sites. Valid sites are %s" % (site,str(accepted_sites))
                        if sitelist.count(site) > 1:
                                return "Site \"%s\" can only be specified once." % (mysite,)
                return 0
        ## END site_processor()
        
        def filetype_processor(self):
                """
                Makes sure that the file type (a.k.a. file extension), matches one that is allowed in
                the database.
                """
                thistype = self.attr['fileType']['Value']
                accepted_types = self.attr['fileType']['Accepted_values']
                if not accepted_types.count(thistype):
                        return "File type \"%s\" is not accepted. Acceptable values %s" % (thistype,str(accepted_types))
                else:
                        return 0
        ##END filetype_processor()
        
        def group_processor(self):
                """
                Makes sure the given group paramter is one of the allowed values.
                """
                group = self.attr['group']['Value']
                accepted_groups = self.attr['group']['Accepted_values']
                if not accepted_groups.count(group):
                        return "Group \"%s\" is not among the accepted list. Acceptable values %s" % (group,str(accepted_groups))
                else:
                        return 0
        ##END group_processor()
        
        
        def __init__(self):
                self.attr = {
"size":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Size of file",
        'Test_method':None
       },
"md5":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"md5 sum of file",
        'Test_method':None
        },
"interferometer":{
        'Value':None,
        'Default':"Null",
        'Type':"string",
        'UserSet':True,
        'Cli_arg_long':"interferometer",
        'Cli_arg_short':"i",
        'Description':"Interferometer site+ifonumber, e.g. H1,H2,L1,H1H2",
        'Accepted_values':['H1','H2','L1','G1'],
        'Test_method':getattr(self,"ifo_processor")
        },
"site":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'Type':"string",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Interferometer site. e.g. H, L, GHLV",
        'Accepted_values':['H','L','G','V'],
        'Test_method':getattr(self,"site_processor")
        },
"fileType":{
        'Value':None,
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':"type",
        'Cli_arg_short':"t",
        'Description':"Type of file",
        'Accepted_values':['gwf','sft','xml'],
        'Test_method':getattr(self,"filetype_processor")
        },
"frameType":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"frame type, e.g. R, RDS_R_L1, RDS_R_L1, h_of_t",
        'Test_method':None
        },
"gpsStart":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"gps-start-time",
        'Cli_arg_short':"s",
        'Description':"GPS start time of file data",
        'Test_method':getattr(self,"gps_processor")
        },
"gpsEnd":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"gps-end-time",
        'Cli_arg_short':"e",
        'Description':"GPS end time of file data",
        'Test_method':getattr(self,"gps_processor")
        },
"duration":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Time duration of file",
        'Test_method':getattr(self,"gps_processor")
        },
"locked":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"IFO locked bitmask",
        'Test_method':None
        },
"scienceMode":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"science mode bitmask",
        'Test_method':None
        },
"playground":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Bit mask for playground data",
        'Test_method':None
        },
"runTag":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"LDR run tag",
        'Test_method':None
        },
"group":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"group",
        'Cli_arg_short':"g",
        'Accepted_values':['pulsar','burst','cw','other'],
        'Description':"Analysis group to which this data is relevant. e.g. pulsar, burst.",
        'Test_method':getattr(self,"group_processor")
        },
"publisher":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"publisher",
        'Cli_arg_short':"p",
        'Description':"Name of person publishing this file.",
        'Test_method':None
        },
"author":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"author",
        'Cli_arg_short':"a",
        'Description':"Name of file's creator. That is the person that ran the code which generated this data.",
        'Test_method':None
        },
"comment":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"comment",
        'Cli_arg_short':"c",
        'Description':"An arbitrary comment describing this file. e.g. \"generated Big Bang Search v0.1a. with fftw v16.2.4 alpha\". (quotes on command line are necessary)",
        'Test_method':None
        }
                             } # end of attr dict.
## END class LSCfileAddMetadataAttr(object)
                
class LSCFile(LSCfileAddMetadataAttr):
        """
        This class contains routines for publishing and unpublishing files
        into and out of an LDR instance.
        """
        
        def __init__(self):
                LSCfileAddMetadataAttr.__init__(self)
                
                # flags for metadata/rls removal
                self.RMV_PFN      = 0x0001
                self.RMV_RLS      = 0x0002
                self.RMV_METADATA = 0x0004
                self.RMV_COMPLETE = self.RMV_RLS | self.RMV_METADATA
                
                # Initialize LDR stuff
                self.config = LDRUtil.getConfig("local", "DEFAULT")
                self.metadata = LDRMetadataCatalog.LDRMetadataCatalog(self.config)
                self.rliurl = self.config["rli"]
                self.lrc = RLS.getCatalog(self.rliurl)
                self.gsiurl = self.config.get("Storage", "hostname")
                
                # list containing files successfully published
                self.successes = successList()
                # list of failure tuples
                self.failures = failureList()
        ## END __init__(self)
        
        def __attr_self_test_methods(self,filename = "Filename not supplied"):
                for field,vals in self.attr.iteritems():
                                if vals['Test_method'] is not None:
                                        result = vals['Test_method']()
                                        if result:
                                                self.failures.append((filename,result))
                                                raise LSCfileAddException, "Error, skipping file: %s" % (result,)
        ## END __attr_self_test_methods(self,filename = "Filename not supplied")
        
        def __create_lfn_pfn_strings(self,filename = "NO_NAME",urlType = ['file']):
                """
                Creates appropriately formatted lfn and pfn strings.
                """
                if urlType is "file":
                        pfn = "file://localhost" + filename # need to switch on --url-type (e.g. file, gsiftp, etc.)
                else:
                        msg = "Error creating lfn pfn pair for %s . urlType %s not yet defined." % (filename,urlType)
                        raise LSCfileAddException, msg
                
                lfn = os.path.basename(filename)
                
                return (lfn,pfn) # MUST BE IN EL FN to PEA FN ORDER!!!!
        
        def __get_attribs_from_filename(self,filename,lfn):
                """
                Sets the following metadata attributes from the given filename.
                (filename must be an lfn)
                
                site
                frameType
                gpsStart
                duration
                (gpsEnd generated from the above two later)
                size
                """
                filepat = re.compile(r'^(\w+)-([\w\d]+)\-(\d+)\-(\d+)\..+')
                if filepat.search(lfn):
                        parsedfilename = filepat.search(lfn).groups()
                else:
                        msg = "Invalid filename format \"%s\"" % (lfn,)
                        print >>sys.stderr, "%s Skipping file %s" % (msg,filename)
                        self.failures.append((filename,msg))
                        raise LSCfileAddException
                      
                # Must have 4 parts to name field
                if len(parsedfilename) is not 4:
                        msg = "Invalid filename format \"%s\"" % (lfn,)
                        print >>sys.stderr, "%s Skipping file %s" % (msg,filename)
                        self.failures.append((filename,msg))
                        raise LSCfileAddException
                        
                # set the name fields in attribute dictionary
                #  NEEDS TO BE RECONCILED WITH --gps-start-time and --gps-end-time OPTIONS!!!!!
                self.attr['site']['Value'] = parsedfilename[0]
                self.attr['frameType']['Value'] = parsedfilename[1]
                self.attr['gpsStart']['Value'] = parsedfilename[2]
                self.attr['duration']['Value'] = parsedfilename[3]
                
                # fill in appropriate fields
                self.attr['size']['Value'] = os.path.getsize(filename)
        ## END __get_attribs_from_filename(self,filename,lfn)
                
        def publish(self,attributes = {}, filelist = [], urlType = "file", host = "", port = 0):
                """
                Adds a lfn <-> pfn mapping. After checking for existance
                of previous mapping, and calculating md5s and any file
                format specific checksums?
                """

                # import use specified attributes
                self.attr = attributes # dumb importation at the moment
                # attempt to publish each specified file
                for filename in filelist:
                        # see if the physical file exists?
                        filename = os.path.abspath(filename)
                        if not os.path.isfile(filename):
                                # print error and _skip_ files which do not exist (or are not files)
                                msg = "Filename %s does not exist (or is a directory)." % (filename,)
                                print >>sys.stderr,  "%s Skipping." % msg
                                failures.append((filename,msg))
                                continue
                        # create lfn<->pair
                        lfn, pfn = self.__create_lfn_pfn_strings(filename,urlType)
                        
                        # see if it already exists in database (respect --replace???)
                        #   Check also for LDR version, for the S4/S5 LDR, no metadatadeletion
                        #   will be supported.
                        #metaexists = self.metadata.exists(lfn)
                        ## DEBUG
                        metaexists = 0
                        ## END DEBUG
                        #if metaexists:
                        #        failures.append((filename,"Metadata for this lfn already exists."))
                        if not metaexists:
                                # Get extension here, may want to process files
                                #   in an extension dependant manner in the future
                                dummy, extension = os.path.splitext(lfn)
                                extension = extension.strip(".")
                                self.attr['fileType']['Value'] = extension
                                
                                try:
                                        self.__get_attribs_from_filename(filename,lfn)
                                
                                except LSCfileAddException:
                                        # skip this file
                                        continue

                                # perform any other consistancy checks
                                try:
                                        self.__attr_self_test_methods(filename)
                                        
                                except LSCfileAddException, e:
                                        print >>sys.stderr, e
                                        continue
                                
                                # calc md5sum, and other checksums
                                # switch on fileType?, perform data format specific checksums?
                                self.attr['md5']['Value'] = self.computeMD5(filename)
                                # enter metadata into database
                                self.addmetadata(lfn)
                                
                        # create lfn, pfn pair in LRC....
                        #if self.rli_lfn_exists(lfn):
                        #        self.lrc_add(lfn,pfn)
                                print "Will create lfn<->pfn mapping for\n%s <-> %s" % (lfn,pfn)
                        #else:
                        #        self.lrc_create_lfn(lfn,pfn)
                        # if all DB additions and checks were successful, then
                        self.successes.append(filename)
                        
                        # END loop over filelist 
        ## END def publish(self)
                
        def remove(self,dalist,flags):
                """
                Removes files described by dalist. dalist can be either
                1) a string containing an LFN or PFN
                2) a plain list of PFNs or LFNs
                   the PFN^LFN type is determined by
                   a urlpase of the first element.
                   If its a URL, assume all elements are PFNs
                   else all elements are assumed to be LFNs
                3) a tuple (a,[b]):
                        a) descriptor string, can be one of "SQL","LFN","PFN"
                        b) list, depends on the value of a
                           if a is "SQL"
                                an SQL query
                           if a is "LFN"
                                list of LFNs to remove
                           if b is "PFN"
                                list of PFNs to remove
                        
                """
                
                
                # convert list into appropriate tupletypes
                if type(dalist) is StringType:
                        if dalist == os.path.basename(dalist):
                                dalist = ("LFN",[dalist])
                        else:
                                dalist = ("PFN",[dalist])
                elif type(dalist) is ListType:
                        # see if first element is a url
                        if not urlparse.urlsplit(dalist[0])[0]:
                                # first element is url
                                # assume all others are urls as well
                                # --> it's a PFN list
                                dalist = ("PFN",dalist)
                        else:
                                dalist = ("LFN",dalist)
                
                if not type(dalist) is TupleType:
                        msg = "remove: Description of files (first arg) is wrong type. Please see doc string."
                        raise LSCfileAddException, msg
                
                # Finally, if dalist is an SQL query,
                # make query and build list from there
                if dalist[0] == "SQL":
                        try:
                                templist = files.query_for_filelist(dalist[1])
                        except LSCfileAddException, e:
                                msg = "remove: could not query for filelist with query \"%s\". reason: %s" % (dalist[1],e)
                                raise LSCfileAddException, msg
                                
                        # see if first element of result is a url
                        if not urlparse.urlsplit(templist[0])[0]:
                                # first element is url
                                # assume all others are urls as well
                                # --> it's a PFN list
                                dalist = ("PFN",templist)
                        else:
                                dalist = ("LFN",templist)
                        
                        del templist
                        
                !!!!# CHECK FOR EMPTY FILE LISTS!!!!!!!!!
                
                # now all lists should be in decent shape
                # metadata must be removed first, as it is keyed on LFNs
                # remove all metadata associated with dalist
                
                if flags | self.RMV_METADATA:
                        # must be done BEFORE rls unmapping,
                        #  at least if dalist[1] is a pfn list
                        # if given a list of PFNs
                        if dalist[0] == "PFN":
                                for pfn in dalist[1]:
                                        try:
                                                lfn = lrc.get_lfn(pfn)
                                        except rlsClient.RlsClientException:
                                                ######### FAILURE #########
                                                continue
                                
                                        try:
                                                self.metadata.delete(lfn)
                                        except:
                                                ######### FAILURE #########
                                                continue
                        elif dalist[0] == "LFN":
                                for lfn in dalist[1]:
                                        try:
                                                self.metadata.delete(lfn)
                                        except:
                                                ########## FAILURE ########
                                                continue
                                
                if flags | self.RMV_PFN:
                        # remove all LFNs associated with dalist in RLS database
                        if dalist[0] == "PFN":
                                for pfn in dalist[1]:
                                        try:
                                                lfn = lrc.get_lfn(pfn)
                                        except rlsClient.RlsClientException:
                                                ######### FAILURE #########
                                                continue
                                        try:
                                                lrc.delete(lfn,pfn)
                                        except rlsClient.RlsClientException:
                                                ######### FAILURE #########
                                                continue
                        else:
                                ########### FAILURE ###########
                                ######### CANNOT SIMPLY DELETE PFN IF ONLY HAVE LFNs ####
                                raise LSCfileAddException
                                
                if flags | self.RMV_RLS:
                        if dalist[0] == "PFN":
                                #for pfn in dalist[1]:
                                try:
                                        templist = lrc.get_lfn_bulk(dalist[1])
                                except rlsClient.RlsClientException:
                                        ######### FAILURE #########
                                        continue
                                try:
                                        lrc.delete_bulk(templist)
                                except rlsClient.RlsClientException:
                                        ######### FAILURE #########
                                        continue
                                del templist
                        elif dalist[0] == "LFN":
                                #### MAY WANT TO CHANGE THE METHODS IN THE FUTURE !!!!!!!!!!!!!!!!!!!!!!!!!
                                try:
                                        templist = lrc.get_pfn_bulk(dalist[1])
                                except rlsClient.RlsClientException:
                                        ############## FAILURE ###########
                                        continue
                                try:
                                        lrc.delete_bulk(templist)
                                except rlsClient.RlsClientException:
                                        ######### FAILURE #######
                                        continue
                                del templist
                                        
                if flags | self.RMV_COMPLETE:
                        # remove all RLS and metadata entries associated with dalist
                        # should already been completed by now due to "flags | " matching
                        pass
                
        ## END def remove(self,dalist,flags)
        
        def mv(self,source,destination):
                """
                Used like UNIX /bin/mv on RLS entries.
                
                usage:
                        foo.mv_rls("<Original PFN>","<Destination PFN>")
                
                Works by 
                1) Checking for existance
                2) Adding <Destination PFN>
                3) Deleting <Original PFN>
                4) Verifying that <Destination PFN> is mapped by it's LFN
                """
                lrc.attr_modify(lfn,attr)
                try:
                        if not lrc.pfn_exists(source):
                                msg = "Source PFN, %s, does not exist. Nothing done." % (source,)
                                raise LSCfileAddException, msg
                except rlsClient.RlsClientException, e:
                       msg = "Caught exception running lrc.pfn_exists(%s)" % (source,)
                       msg = msg + str(e)
                       raise LSCfileAddException, msg
                else:
                        try:
                                lfn = lrc.get_lfn(source)
                        except rlsClient.RlsClientException, e:
                                msg = "Caught exception running lrc.get_lfn(%s)" % (source,)
                                msg = msg + str(e)
                                raise LSCfileAddException, msg
                        
                        try:
                                lrc.add(lfn,destination)
                        except:
                                msg = "Caught exception running lrc.add(%s,%s)" % (lfn,destination)
                                msg = msg + str(e)
                                raise LSCfileAddException, msg
                        
                        try:
                                lrc.delete(lfn,destination)
                        except rlsClient.RlsClientException, e:
                                msg = "Caught exception running lrc.delete(%s,%s)" % (lfn,destination)
                                msg = msg + str(e)
        ## END mv(self,source,destination)
        
        def addmetadata(self,lfn):
                self.metadata.add(lfn)
                for field, val in self.attr.iteritems():
                        self.metadata.set_attr(lfn,field,val['Value'])
        ## END def addmetadata(self,lfn)
        
        def computeMD5(self,filename):
                """
                Compute md5sum of a file. Must be given a filesystem
                path.
                """
                m = md5.new()
                f = open(filename,"r")
                aMeg = 1048576  # 1024 *1024
                line = f.read(aMeg)
                while line:
                        m.update(line)
                        line = f.read(aMeg)
                return m.hexdigest()
        ## END computeMD5(self,filename)
        
        
        
## END class LSCFile(LSCfileAddMetadataAttr)


class CLIUtil(LSCfileAddMetadataAttr):
        """
        Contains methods etc. for handling the command line.
        Some of these methods set up the metadata field dictionaries
        """
        
        def __init__(self):
                """
                Sets up appropriate strings and dictionaries.
                The parameters here should reflect available database fields.
                """
                LSCfileAddMetadataAttr.__init__(self)
		# Some class attributes
		self.shortop = ""
		self.longop = []
                self.urlType = None
		self.hostPortString = None
                self.host = None
		self.port = None
		# defaults
		self.default_port = 30100
		self.default_urlType = ['file']
		# maps command line options to their appropriate fields
		self.cli_short_name = {}
		self.cli_long_name = {}
		# contains results for non-metadata cli parameters
		self.nonmetaparam = {}
		# list of parameter tuples
		self.params = []
		self.filelist = []
                #Initializes some shorthand variables from the attr dictionary.
                for field, vals in self.attr.iteritems():
                        if vals['UserSet']:
                                if vals['Cli_arg_short']:
                                        self.cli_short_name[vals['Cli_arg_short']] = field
                                if vals['Cli_arg_long']:
                                        self.cli_long_name[vals['Cli_arg_long']] = field
                # Non-metadata specific fields, e.g. url-type etc.
                self.shortop = "u:s:hv"
                self.longop = [
                        "help",
                        "url-type=",
                        "server=",
                        "verbose"
                        ]
                # more defaults
                self.port = self.default_port
     
                # now collect non-metadata specific args with LSCfileAddMetadataAttr args
                #    for use in getopts
                for field, vals in self.cli_short_name.iteritems():
                        self.shortop = self.shortop + field + ":"
                for field, vals in self.cli_long_name.iteritems():
                        field = field + "="
                        self.longop.append(field)
        ## END __init__
        
        
        def class_sanity_check(self):
                """
                Meant to be run by the programmer to make sure that this class is sane.
                For example, this checks to make sure the attribute dictionary is self consistant
                (should be in LSCfileAddMetadataAttr...), but also compares with non-metadata
                specific CLI args to make sure nothing gets clobbered inappropriately.
                """
		### PRESNENTLY NOT WORKING EITHER
                exit = "NO"
                for op in self.shortop.split(":"):
                        for field,vals in self.attr.iteritems():
                                if vals['Cli_arg_short'] == op:
                                        exit = "YES"
                                        print >>sys.stderr, "Option collision shortop is \"%s\"\nField is %s" % (op,field)
                for op in self.longop:
                        for field,vals in self.attr.iteritems():
                                if vals['Cli_arg_short'] == op:
                                        exit = "YES"
                                        print >>sys.stderr, "Option collision shortop is \"%s\"\nField is %s" % (op,field)
                
                if exit is "YES":
                        print >>sys.stderr, "Sanity Check Failed."
                        sys.exit(123)
        ## END class_sanity_check(self)
        
        def __put_opts_in_place(self,opts):
                """
                Populates various data structures as per command line.
                """
                for o, a in opts:
                        # strip leading "-"'s
                        o = o.lstrip("-")
                        
                        if o == "h" or o == "help":
                                self.print_usage()
                                sys.exit(0)
                        elif o == "u" or o == "url-type":
                                self.nonmetaparam['url-type'] = a
                        elif o == "s" or o == "server":
                                self.nonmetaparam['server'] = a
                        elif o == "v" or o == ['verbose']:
                                self.nonmetaparam['verbose'] = True
                        elif not self.cli_short_name.has_key(o) and not self.cli_long_name.has_key(o):
                                # invalid parameter
                                print >>sys.stderr, "Bad option, %s" % (o,)
                                print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                                sys.exit(123)
                        elif self.cli_short_name.has_key(o):
                                self.attr[self.cli_short_name[o]]['Value'] = a
                        elif self.cli_long_name.has_key(o):
                                self.attr[self.cli_long_name[o]]['Value'] = a
        ## END __put_opts_in_place(self,opts,args)
        
        def __process_nonmetadata_opts(self,args):
                """
                processes nonmeteadata options. e.g. host and port URL parts.
                """
                ## environment variables override defaults but not
                ## command line options
                # Configure serverl url
                hostPortString = None
                try:
                        hostPortString = os.environ['LSC_FILEADD_SERVER']
                except:
                        pass
                try:
                        hostPortString = self.nonmetaparam['server']
                except:
                        pass
                
                # URL type of pfns to publish for this session
                if not self.nonmetaparam.has_key('urlType'):
                        try:   
                                self.nonmetaparam['urlType'] = os.environ['LSC_FILEADD_URL_TYPE'] 
                        except:
                                self.nonmetaparam['urlType'] = self.default_urlType;
                self.urlType = self.nonmetaparam['urlType']
                
                # determine server and port
                if not hostPortString:
                        print >>sys.stderr, "No LDRfileAddServer specified"
                        print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                        sys.exit(1)

                if hostPortString.find(':') < 0:
                        # no port specified
                        host = hostPortString
                else:
                        # server and port specified
                        host, portString = hostPortString.split(':')
                        self.host = host
                        self.port = int(portString)
                
                # See if any files were specified
                if not args: # empty file list
                        print >>sys.stderr, "You must specify at least one filename."
                        print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                        sys.exit(10)
                else:
                        self.filelist = args
        ## __process_nonmetadata_opts(self)
        
        def get_user_parameters(self):
                """
                Grabs data from command line, user environment, etc.
                and sets the appropriate variables to be used later.
                """
                # Get options and args from command line
                try:
                        opts, args = getopt.getopt(sys.argv[1:], self.shortop, self.longop)
                except getopt.GetoptError:
                        print >>sys.stderr, "Error parsing command line"
                        print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                        sys.exit(1)
                        
                # Process options and arguments
                self.__put_opts_in_place(opts)
                
                ## Handle non-metadata options
                self.__process_nonmetadata_opts(args)
                
        ## END def get_user_parameters
                        
        def print_usage(self):
                """
                Prints a usage message to stderr.
                Usage will be partially generated dynamically from the
                LSCfileAddMetadataAttr dictionary.
                """
                msg = """\
NAME
        LSCfileAdd: Publishes a file with appropriate metadata to an LDR 
        database.

SYNOPSIS
        LSCfileAdd --help
        
        LSCfileAdd <options> file1 file2 ...
        
\
"""
                msg += "DESCRIPTION\n"
                for field, vals in self.attr.iteritems():
                        if vals['UserSet']:
                                msg += "\t-%s, --%s\n" % (vals['Cli_arg_short'],vals['Cli_arg_long'])
                                # format long description lines properly for standard terminal
                                description = vals['Description']
                                length = len(description)
                                newdes = ""
                                oldidx = 0
                                idx = 55
                                while idx < length:
                                        while not description[idx].isspace() and idx > 0:
                                                idx -= 1
                                        newdes += "%s\n\t\t" % (description[oldidx:idx],)
                                        oldidx = idx+1
                                        idx += 55
                                if newdes:
                                        description = "%s%s" % (newdes,description[oldidx:length])
                                msg += "\t\t%s\n\n" % (description,)
                
                msg += """\
ENVIRONMENT

        LSC_FILEADD_SERVER defines the database server where metadata and 
                           lfn<->pfn mappings will go, overridden by the 
                           --server option.

        LSC_FILEADD_URLTYPE defines the url type for lfn<->pfn mappings. 
                            Same arguments as --url-type. This is also 
                            overridden by that option.

        ....

EXAMPLE
$ LSCfileAdd .....

\
"""

                print >>sys.stderr, msg
        ## END def print_usage():

## END class CLIUtil(LSCfileAddMetadataAttr)


class GridUtil:
        """
        Class containing grid-specific utilities.
        """
        def checkCredentials(self):
                """
                Check to make sure that the proper Grid Credentials (a proxy certificate) is
                available in order to authenticate to the remote LDRdataFindServer.
                """
                # verify that we have access to credentials
                try:
                        proxyText = security.grid_proxy_info()
                except Exception, e:
                        print >>sys.stderr, "Error verifying credentials: %s" % e
                        print >>sys.stderr, "Run 'grid-proxy-init' to generate a proxy certificate"
                        sys.exit(1)

                pat = re.compile(r'timeleft : (\d{1,3}):(\d\d):(\d\d)')

                try:
                        if isinstance(proxyText, str):
                                m = pat.search(proxyText)
                        elif isinstance(proxyText, tuple):
                                m = pat.search(proxyText[0])
                        else:
                                raise RuntimeError, "bad format for proxyText in checkCredentials"
                                
                        hours, minutes, seconds = map(int, m.groups())
                except Exception, e:
                        print >>sys.stderr, "Error parsing proxy information: %s" % e
                        sys.exit(1)

                timeleft = seconds + 60 * minutes + 3600 * hours

                if timeleft < 300:
                        print >>sys.stderr, "Less than 5 minutes left for proxy certificate."
                        print >>sys.stderr, "Run 'grid-proxy-init' to generate a new proxy certificate"
                        sys.exit(1)
            ## END checkCredentials()
                        
## END GridUtil(object)


class LSCfileAddException(exceptions.Exception):
        """
        Exceptions raised by the classes and methods in this client
        will be instances of this class.
        """
        def __init__(self, args=None):
                """
                Create an instance of this class, ie. an LSCfileAddException

                @param args:

                @return: Instance of class LSCfileAddException
                """
                self.args = args
                
## END LSCaddFileClientException(exceptions.Exception)


class LSCfileAddClient:
        """
        Class that interacts with the LDRfileAddServer.
        """
        pass
## END LSCaddFileClient(object)


