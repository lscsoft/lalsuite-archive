"""
LSCfileAdd.py module contains the classes for the LSCfileAdd utility
"""

__author__="Ben Johnson: module contains the classes for the LSCfileAdd utility"
__date__='$Date$'
__version__='$Revision$'[0:0]

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
                       7) Requirements # list of other fields required to be set if this one is set by user,
                                       # if Requirements contains the field name
                                       # then it is required to be set on the command line.
                                       # e.g. self.attr["runTag"]["Requirements"] contains a "runTag" element.
                                       
        
        """
        attr = {}
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
        'Requirements':['size']
       },
"md5":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"md5 sum of file",
        'Requirements':['md5']
        },
"interferometer":{
        'Value':None,
        'Default':"Null",
        'Type':"string",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Interferometer site+ifonumber, e.g. H1,H2,L1",
        'Requirements':None
        },
"site":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'Type':"string",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Interferometer site+. e.g. H, L, GHLV",
        'Requirements':None
        },
"fileType":{
        'Value':None,
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':"type",
        'Cli_arg_short':"t",
        'Description':"Type of file",
        'Requirements':['filetype']
        },
"frameType":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"frame type, e.g. R, RDS_R_L1, RDS_R_L1, h_of_t",
        'Requirements':None
        },
"gpsStart":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"gps-start-time",
        'Cli_arg_short':"s",
        'Description':"GPS start time of file data",
        'Requirements':['gpsStart']
        },
"gpsEnd":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"gps-end-time",
        'Cli_arg_short':"e",
        'Description':"GPS end time of file data",
        'Requirements':['gpsEnd']
        },
"duration":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Time duration of file",
        'Requirements':['duration']
        },
"locked":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"IFO locked bitmask",
        'Requirements':None
        },
"scienceMode":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"science mode bitmask",
        'Requirements':None
        },
"playground":{
        'Value':None,
        'Type':"int",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"Bit mask for playground data",
        'Requirements':None
        },
"runTag":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':False,
        'Cli_arg_long':None,
        'Cli_arg_short':None,
        'Description':"LDR run tag",
        'Requirements':None
        },
"group":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"group",
        'Cli_arg_short':"g",
        'Description':"Analysis group to which this data is relevant. e.g. pulsar, burst",
        'Requirements':None
        },
"publisher":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"publisher",
        'Cli_arg_short':"p",
        'Description':"Name of person publishing this file.",
        'Requirements':None
        },
"author":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UserSet':True,
        'Cli_arg_long':"author",
        'Cli_arg_long':"a",
        'Description':"Name of file's creator. That is the person ran the code which generated the data."
        'Requirements':None
        },
"comment":{
        'Value':None,
        'Type':"string",
        'Default':"Null",
        'UsetSet':True,
        'Cli_arg_long':"comment"
        'Cli_arg_short':"c",
        'Description':"An arbitrary comment describing this file. e.g. \"generated Big Bang Search v0.1a. with fftw v12.2.4\"",
        'Requirements':None
        }
                             } # end of attr dict.
## END class LSCfileAddMetadataAttr(dict)
                
class LSCFile(LSCfileAddMetadataAttr):
        """
        This class contains routines for publishing and unpublishing files
        into and out of an LDR instance.
        """
        
        def __init__(self):
                LSCfileAddMetadataAttr.__init__(self)
                # File is essentially an lfn,pfn pair
                self.lfn = None
                self.pfn = None
        ## END def __init__(self)
        
        
        def publish(self, filelist = []):
                """
                Adds a lfn <-> pfn mapping. After checking for existance
                of previous mapping, and calculating md5s and any file
                format specific checksums?
                """
                # authentication stuff
                #blah
                for filename in filelist:
                        # Set up necessary attributes???
                        #blah?
                        # see if the physical file exists?
                        #blah?
                        # see if it already exists in database (respect --replace???)
                        # Check also for LDR version, for the S4/S5 LDR, no metadatadeletion
                        # will be supported.
                        metaexists = self.metadata.exists(lfn)
                        if metaexists:
                                pass
                                if not metaexists:
                                        # calc md5sum, and other checksums
                                        # switch on fileType?, perform data format specific checksums?
                                        attr.md5['Value'] = computeMD5(self,filename)
                                        # perform any other consistancy checks
                                        #blah
                                        # enter metadata into database
                                        self.addmetadata(fields)
                        # create lfn, pfn pair in LRC....
                        if self.rli_lfn_exists(lfn):
                                self.lrc_add(lfn,pfn)
                        else:
                                self.lrc_create_lfn(lfn,pfn)
        ## END def publish(self)
                
        def remove_all(self):
                """
                Removes ALL metadata and lfn<->pfn mappings associated
                with a particular lfn.
                """
                print >>sys.stderr "Function not yet implemented."
                # authentication stuff
                #blah
                # see if it already exists in database, delete if so
                #if metadata.exists(lfn):
                #       metadata.delete(lfn)
                        # NEED TO MAKE EXCEPTION ABOUT NON EXISTANCE!!!
                # remove ALL lfn, pfn pairs
                #pfn_list = lrc.lrc_get_pfn(lfn,0,0)
                #for mypfn in pfn_list:
                #        self.lrc_delete(lfn,mypfn)
        ## END def remove_all(self)
        
        def addmetadata(self):
                # Fields to be defined
                # the following is from Publisher.py, certainly doesn't work now
                self.metadata.add(lfn)
                for field, val in self.attr.iteritems():
                        self.metadata.set_attr(lfn,field,val["Value"])
        ## END def addmetadata(self)
        
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
        # Some class attributes
        shortop = ""
        longop = []
        hostPortString = None
        port = None
        clientMethodArgDict = {}
        clientMethod = ''
        # maps command line options to their appropriate fields
        cli_short_name = {}
        cli_long_name = {}
        # list of parameter tuples
        params = []
        
        def __init__(self):
                """
                Sets up appropriate strings and dictionaries.
                The parameters here should reflect available database fields.
                """
                LSCfileAddMetadataAttr.__init__
                #Initializes some shorthand variables from the attr dictionary.
                for field, vals in attr.iteritems():
                        if vals['UserSet']:
                                if vals['Cli_arg_short']:
                                        cli_short_name[vals['Cli_arg_short']] = field
                                if vals['Cli_arg_long']:
                                        cli_long_name[vals['Cli_arg_long']] = field
                # Non-metadata specific fields, e.g. url-type etc.
                shortop = "u:s"
                longop = [
                        "url-type=",
                        "server="
                        ]
                # more defaults
                port = 30010
                clientMethodArgDict = {
                        'urlType': None
                }
                # default method OBVIOUSLY NEEDS TO BE CHANGED FOR THIS CLASS
                clientMethod = 'findFrameURLs'
                # now collect non-metadata specific args with LSCfileAddMetadataAttr args
                #    for use in getopts
                for field, vals in cli_short_name.iteritems():
                        shortop = shortop + "field"
                for field, vals in cli_long_name.iteritems():
                        field = field + "="
                        longop.append(field)
        ## END __init__
        
        
        def class_sanity_check(self):
                """
                Meant to be run by the programmer to make sure that this class is sane.
                For example, this checks to make sure the attribute dictionary is self consistant
                (should be in LSCfileAddMetadataAttr...), but also compares with non-metadata
                specific CLI args to make sure nothing gets clobbered inappropriately.
                """
                pass
        ## END class_sanity_check(self)
        
        
        def get_user_parameters(self):
                """
                Grabs data from command line, user environment, etc.
                and sets the appropriate variables to be used later.
                """
                try:
                        opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
                except getopt.GetoptError:
                        print >>sys.stderr, "Error parsing command line"
                        print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                        sys.exit(1)
               
                for o, a in opts:
                        if not cli_short_name[o] and not cli_long_name[o]:
                                # invalid parameter
                                print >>sys.stderr, "Bad option, %s" % (o,)
                                print >>sys.stderr, "Enter 'LSCfileAdd --help' for usage"
                                sys.exit(123)
                        elif cli_short_name[o]:
                                # Still need to check and convert these to appropriate types etc.
                                self.attr[cli_short_name[o]] = a
                                
                # environment variables override defaults but not
                # command line options
                # NEEDS TO BE CUSTOMIZED FOR THE LSCfileAdd SCRIPT
                #try:
                #        hostPortString = os.environ['LSC_FILEADD_SERVER']
                #except:
                #        pass

                #try:   
                #        clientMethodArgDict['urlType'] = os.environ['LSC_FILEADD_URL_TYPE'] 
                #        clientMethod = 'addFileURLsFilter'
                #except:
                #        pass

                #try:
                #        clientMethodArgDict['match'] = os.environ['LSC_FILEADD_MATCH']
                #        clientMethod = 'addFileURLsFilter'
                #except:
                #        pass
                # Now actually process command line arguments
                #if not clientMethod:
                #        print >>sys.stderr, "Bad combination or missing options"
                #        print >>sys.stderr, "Enter 'LSCdataFind --help' for usage"
                #        sys.exit(1)

                # determine server and port
                #if not hostPortString:
                #        print >>sys.stderr, "No LDRdataFindServer specified"
                #        print >>sys.stderr, "Enter 'LSCdataFind --help' for usage"
                #        sys.exit(1)

                #if hostPortString.find(':') < 0:
                        # no port specified
                #        host = hostPortString
                #else:
                        # server and port specified
                #        host, portString = hostPortString.split(':')
                #        port = int(portString)
        ## END def get_user_parameters
                        
        def print_usage(self):
                """
                Prints a usage message to stderr.
                Usage will be partially generated dynamically from the
                LSCfileAddMetadataAttr dictionary.
                """
                msg = """\
NAME
        LSCfileAdd

SYNOPSIS
        LSCfileAdd <options> file1 file2...

DESCRIPTION
        TBD (to be done)

ENVIRONMENT

        LSC_ADDFILES_SERVER ...

        ....

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


