"""
Main classes for publishing S4 data 
and possibly <reverb>beyond...</reverb>

Now being used as LSCfileAdd.py

"""

_author__="Ben Johnson: module contains the classes for the LSCfileAdd utility"
__date__='$Date$'
__version__='$Revision$'[0:0]

import exceptions
import os
import md5
import time
import sys
import re
from types import *



# pulled from LSCfileAdd
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

# LDR-specific modules
# make sure cert/key are set up in environment
if not os.environ.has_key('X509_USER_CERT'):
        msg = 'Must point X509_USER_CERT env variable to cert file!'
        raise LSCfileAddException
if not os.environ.has_key('X509_USER_KEY'):
        msg = 'Must point X509_USER_KEY env variable to key file!'
        raise LSCfileAddException
import rlsClient
from LDRMetadata import metadata
from Metadata import Metadata
from Metadata import MetadataException

class Publisher(object):
        """
        Class for publishing Engineering and/or science run data at
        at the LIGO observatories.
        """
        def __init__(self,mysqlURL,rlsURL):
                #LSCrunMetadataAttr.__init__(self)
                # set some default flags
                self.PRESERVE_METADATA  = 0x01
                self.OVERWRITE_METADATA = 0x02
                
                # store urls for error reporting purposes
                self.mdURL = mysqlURL
                self.rlsURL = rlsURL
                
                # Initialize LDR stuff
                ### NEED TO CATCH EXCEPTION HERE!!!
                self.md = metadata(mysqlURL)
                self.rls = rlsClient.RlsClient(rlsURL)
                
                # setup variables for periodic handle regeneration
                self.regenHandle_timeout =  300.0
                self.regenHandle_timer = time.time()
                
                # set default attrib generator handle
                self.attribute_generator_attach(self.default_attribute_generator)
                
                # list containing files successfully published
                self.successes = []
                # list of tuples containing files which failed to be published
                #    and the reasons they were not (hopefully).
                self.failures = []
        ## END __init__
        
        def checkHandle(self):
                """
                Sees if regenHandle_timer is older than now by,
                regenHandle_timeout. If so it updates the timer to the current time,
                and regenerates the metadata and rls handles.
                """
                now = time.time()
                if now - self.regenHandle_timer > self.regenHandle_timeout:
                        print "checkHandle branch taken."
                        sys.stdout.flush()
                        self.regenHandle_timer = now
                        self.regenHandle()
        ## END checkHandle(self)
        
        def regenHandle(self):
                """
                Regenerates handle to metadata and rls hosts
                """
                # "close" connections
                try:
                        self.md.close()
                        del self.rls
                except Exception,e:
                        msg = "Caught exception while closing metadata and/or RLS handles. Going ahead with regeneration attempt. Error was: %s " \
                        % (str(e),)
                        print >>sys.stderr, msg
                
                # reopen connections
                try:
                        self.md = metadata(self.mdURL)
                except Exception, e:
                        msg = "Caught exception while attempting to regenerate handle to metadata host, %s. Error was: %s" % (self.mdURL,str(e))
                        raise LSCfileAddException, msg
                try:
                        self.rls = rlsClient.RlsClient(self.rlsURL)
                except Exception, e:
                        msg = "caught exception while attempting to generate handle to RLS host, %s. Error was: %s" % (self.rlsURL,str(e))
                        raise LSCfileAddException, msg
                        
                print >>sys.stderr, "Regenerated the metadata and rls handles at %s" % time.asctime()
                sys.stderr.flush()
        ## END regenHandle(self)
        
        def default_attribute_generator(self,data = None):
                # simply return given object unmodified
                return data
        ## END default_attribute_generator(lfn)
        
        def attribute_generator_attach(self,myobject):
                self.attribute_generator = myobject
        ## END attribute_generator_attach(myobject)
        
        def publish(self, datalist, flags = 0x01):
                """
                Adds a lfn <-> pfn mapping, and adds metadata provided or
                calculated with the Publisher.attribute_generator() method.
                
                Grabs appropriate locked/sciencemode segments from segments dictionary.
                
                Paramters:
                        datalist == dictionary containing LFNs and associated
                                    PFNs and metadata to publish. See below
                                    for details.
                        
                        flags == Determines whether or not to overwrite
                                 Existing entries. Defaults to PRESERVE_METADATA.
                                 The other option is OVERWRITE_METADATA.
                
                datalist format:
                        datalist is a list of dictionaries, ultimately keyed on LFN
                        
                        datalist[n] --> ["name"] --> "<LFN>"
                                    --> ["urls"] --> ["<PFN1>","<PFN2>",...]
                                    --> ["pset"] --> "<name of publishing set>"
                                    --> ["metadata"] --> ["<attribute>"]["Value"]
                """
                
                ### Need to perform some basic checking of the
                ### data list structure...                
                
                # to be implemented
                
                # see if database handles need to be regenerated
                self.checkHandle()
                
                #### Attempt to publish metadata
                metadata_failure = True  # though, any metadata failure should 
                                         # most likely warrant raising an exception
                for data in datalist:
                        lfn = data["name"]
                       
                        # check for existance
                        if self.md.getAttribute(lfn,"md5"):
                                metaexists = True
                                print "meta md5 %s" % (self.md.getAttribute(lfn,"md5"),)
                        else:
                                metaexists = False

                        if metaexists & (flags & self.PRESERVE_METADATA):
                                metadata_failure = False
                                #self.failures.append((filename,"Metadata for this lfn already exists."))
                        else:
                                try:
                                        data = self.attribute_generator(data)
                                except LSCfileAddException, e:
                                        msg = "Caught exception while running the attribute_generator() on LFN \"%s\":\
                                        Message received was: %s" % (str(lfn),str(e))
                                        metadata_failure = True
                                        raise LSCfileAddException
                                        
                                # enter metadata into database
                                try:
                                        self.add_metadata(data["pset"],lfn,data["metadata"])
                                except Exception, e:
                                        metadata_failure = True
                                        msg = "Caught exception while publishing metadata self.add_metadata(%s,%s,%s)! Error was: %s" % (data["pset"],lfn,data["metadata"],str(e))
                                        raise LSCfileAddException, msg
                                
                                # everything successful? hope so ;)
                                metadata_failure = False
                #### Attempt to publish RLS entries
                rls_failure = True
                # list of mapdicts
                maps = []
                
                for data in datalist:
                        # build mapdicts of the form
                        # mapdict[lfn] = SINGLE_URL_STRING
                        # then make list of these mapdicts
                        if not len(data['urls']):
                                # no urls to publish
                                print "Publisher.publish(): No urls to publish for lfn \"%s\"" % (data['name'],)
                        else:
                                i = 0
                                lfn = data['name']
                                for url in data['urls']:
                                        if i == len(maps):
                                                maps.append({})
                                        maps[i][lfn] = url
                                        i = i + 1
                                
                for mapdict in maps:                              
                        if mapdict:
                                try:
                                        rls_failed = self.rls.lrc_create_lfn_bulk(mapdict)
                                except rlsClient.RlsClientException,e:
                                        msg = "Caught exception from self.rls.lrc_create_lfn_bulk(). Error was: %s" % (str(e),)
                                        msg += "\nRegenerating handle to LDR and RLS and trying again."
                                        print >>sys.stderr,msg
                                        self.regenHandle()
                                        try:
                                                rls_failed = self.rls.lrc_create_lfn_bulk(mapdict)
                                        except rlsClient.RlsClientException,e:
                                                msg = "Caught exception from self.rls.lrc_create_lfn_bulk(). Again!. Error was: %s" % (str(e),)
                                                print >>sys.stderr,msg
                                                raise LSCfileAddException,msg
                                        
                                if rls_failed:
                                        msg = "Some RLS mappings failed while running self.rls.lrc_create_lfn_bulk()."
                                        msg += "\nTrying the self.rls.lrc_add_bulk() method."
                                        print >>sys.stderr,msg
                                        try:
                                                temp_failed = self.rls.lrc_add_bulk(rls_failed)
                                        except rlsClient.RlsClientException,e:
                                                msg = "Caught exception while running self.rls.lrc_add_bulk(). Error was: %s" % (str(e),)
                                                msg += "\nRegenerating handle to LDR and RLS and trying again."
                                                print >>sys.stderr,msg
                                                try:
                                                        temp_failed = self.rls.lrc_add_bulk(rls_failed)
                                                except rlsClient.RlsClientException,e:
                                                        msg = "Caught exception while running self.rls.lrc_add_bulk(). Again!. Attempted to map, %s" % (rls_failed,)
                                                        msg += " Error was: %s" % (str(e),)
                                                        raise LSCfileAddException,msg
                                        if temp_failed:
                                                ## Check for existance
                                                temp_failed2 = self.rls.lrc_get_pfn_bulk(temp_failed.keys())
                                                for key in temp_failed2.keys():
                                                        if isinstance(temp_failed[key],ListType):
                                                                for item in temp_failed[key]:
                                                                        if not temp_failed2[key].count(item):
                                                                                msg = "For some reason, could not add PFN %s" % (item)
                                                                                raise LSCfileAddException,e
                                                        else: # assume StringType
                                                                if not temp_failed2[key].count(temp_failed[key]):
                                                                        msg = "For some reason, could not add PFN %s" % (item)
                                                                        raise LSCfileAddException,e
                                                # At this point this code decalres all is good,
                                                # these URLs were published at some earlier time
                                                msg = "Could not add, %s. These URLs already exist in the database." % (temp_failed,)
                                                rls_failure = False
                                        else:
                                                # declare rls success
                                                rls_failure = False
                                else:
                                        # declare rls success
                                        rls_failure = False
        
        
                                # if all DB additions and checks were successful, then
                                if not rls_failure and not metadata_failure:
                                        #self.successes.append(data)
                                        pass
                                else:
                                        #self.failures.append(data)
                                        temp = ""
                                        if metadata_failure:
                                                temp += "metadata failure, "
                                        if rls_failure:
                                                temp += "rls failure"
                                        msg = "There is an algorithmic failure with the publish() method. Exception wasn't handled properly? The failures were of type, %s" % (temp,)
                                        raise LSCfileAddException,msg

        ## END def publish(self)
                
        def completely_remove(self,lfnlist):
                """
                Removes ALL metadata and lfn<->pfn mappings associated
                with a particular lfn.
                """
                # see if database handles need to be regenerated
                self.checkHandle()
                
                for lfn in lfnlist:
                        try:
                                # remove ALL lfn, pfn pairs
                                pfnlist = self.lrc.get_pfn(lfn)
                                for mypfn in pfnlist:
                                        self.lrc.delete(lfn,mypfn)
                        except rlsClient.RlsClientException,e:
                                msg = "completely_remove(): Caught RlsClientException while attempting to delete \"%s\" Error was: %s" % (lfn,str(e))
                                raise LSCfileAddException, msg
                        try:
                                # see if it already exists in database, delete if so
                                if self.md.getAttribute(lfn,"md5"):
                                        self.md.delete(lfn)
                        except Exception,e:
                                msg = "completely_remove(): Caught an exception while calling either metadata.exists() or metadtaa.delete() on \"%s\". Error was: %s" % (lfn,str(e))
                                raise LSCfileAddException, msg
        ## END def completely_remove(self)
        
        def mv_url(self,list_to_move,moves_per_loop = 50):
                """
                Used like UNIX /bin/mv on RLS entries.
                
                usage:
                        foo.mv_url(list_of_lists,bulk_quantity)
                
                Where list_of_lists, is a list of sublists., each sublist containing
                a source,destination pair. e.g. list_of_lists = [[source1,dest1],[source2,dest2],...]
                
                bulk_quantity is an optional argument. It is the number of source,destination 
                pairs to publish at once. This parameter can also safely be greater than the 
                length of list_of_lists. It defaults to 50.
                
                
                Upon error, appends error message to failures[] and raises
                an LSCfileAddException. Otherwise adds a message to successes[].
                """
                # first, perform some basic type checking
                if not isinstance(list_to_move,types.ListType):
                        msg = "Must be given a list of lists to move. e.g. [[source1,dest1],[source2,dest2],...]"
                        raise LSCfileAddException,msg
                elif not isinstance(list_to_move[0],types.ListType):
                        msg = "Must be given a list of lists to move. e.g. [[source1,dest1],[source2,dest2],...]"
                        raise LSCfileAddException,msg
                
                # grab first "slice" of list_to_move
                i = 0
                j = moves_per_loop
                f = lambda x: x - len(list_to_move) < 0 and x or len(list_to_move)
                littlelist =  list_to_move[i:f(j)]
                
                while len(littlelist):
                        # see if database handles need to be regenerated
                        self.checkHandle()
                        try:
                                print "Attempting to move,"
                                for name in littlelist:
                                        print name
                                failed = self.rls.rename_pfn_bulk(littlelist)
                        except rlsClient.RlsClientException,e:
                                msg = "Caught exception while attempting a rename_pfn_bulk(). Error was: %s" % (str(e),)
                                raise LSCfileAddException
                        if failed:
                                print >>sys.stderr, "Could not move the following source+destination pair(s): %s" % (failed,)
                        i = j
                        j = i + moves_per_loop
                        littlelist = list_to_move[i:f(j)]
                                                
        ## END mv_url(self,list_to_move,moves_per_loop = 50)
        
        def rm_url(self,pfnlist,rms_per_loop = 50):
                """
                Used similar to UNIX /bin/rm.
                
                usage:
                        Multiple file removal: foo.rm_url(pfnlist)
                        Single file removal: foo.rm_url("single_pfn_string")
                        
                NOTE:
                        If the URL happens to be the only one to which the 
                        listing (LFN) is mapped, the entire listing will also be deleted.
                        This code presently does not alert the user to this, so beware!
                        
                Upon error, appends error message to failures[] and raises
                an LSCfileAddException. Otherwise adds a message to successes[].
                """
               
                
                if isinstance(pfnlist,StringType):
                        pfnlist = [pfnlist]
                elif not isinstance(pfnlist,ListType):
                        msg = "rm_url(var): var must be either a string or a list of strings."
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
                        
                # grab first "slice" of pfnlist
                i = 0
                j = rms_per_loop
                f = lambda x: x - len(pfnlist) < 0 and x or len(pfnlist)
                littlelist =  pfnlist[i:f(j)]
                while len(littlelist):
                        # see if database handles need to be regenerated
                        self.checkHandle()
                        try:
                                for pfn in pfnlist:
                                        if self.rls.lrc_exists(pfn,rlsClient.obj_lrc_pfn):
                                                lfn = self.rls.lrc_get_lfn(pfn)
                                                print "RLS lfn %s" % (lfn,)
                                                self.rls.lrc_delete(lfn[0],pfn)
                                        else:
                                                print 'PFN %s does not exist?' % (pfn,)
                                                
                                #mapdict = self.rls.lrc_get_lfn_bulk(pfnlist)
                                #failed = self.rls.lrc_delete_bulk(mapdict)
                                
                        except rlsClient.RlsClientException, e:
                                msg = "rm_url(%s): Caught RLS exception. Message was: %s" % (littlelist,str(e))
                                raise LSCfileAddException, msg
                        #if failed:
                               # print >>sys.stderr, "Could not delete the following urls %s" % (failed,)
                        i = j
                        j = i + rms_per_loop
                        littlelist = pfnlist[i:f(j)]
        ## END rm_url(self,pfnlist)
        
        def rm_listing(self,lfnlist,rms_per_loop=50):
                """
                Used similar to UNIX /bin/rm.
                
                usage:
                        Multiple LFN removal: foo.rm_listing(lfnlist)
                        Single LFN removal: foo.rm_listing("single_lfn_string")
                        
                Removes a an LFN -> PFN url listing. i.e. After successful
                usage of rm_listing(), clients should not be able to "see"
                that a particular LFN exists.
                        
                Upon error, appends error message to failures[] and raises
                an LSCfileAddException. Otherwise adds a message to successes[].
                """
                if isinstance(lfnlist,StringType):
                        lfnlist = [lfnlist]
                elif not isinstance(lfnlist,ListType):
                        msg = "rm_listing(var): var must be either a string or a list of strings."
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
                # grab first "slice" of pfnlist
                i = 0
                j = rms_per_loop
                f = lambda x: x - len(lfnlist) < 0 and x or len(lfnlist)
                littlelist =  lfnlist[i:f(j)]
                while len(littlelist):
                        # see if database handles need to be regenerated
                        self.checkHandle()
                        try:
                               mapdict = self.rls.pfn_get_bulk(lfnlist)
                               failed = self.rls.lrc_delete_bulk(mapdict)
                        except rlsClient.RlsClientException, e:
                                msg = "rm_listing(%s): Caught RLS exception. Message was: %s" % (name,str(e))
                                self.failures.append(msg)
                                raise LSCfileAddException, msg
                        if failed:
                                print >>sys.stderr, "Could not delete listings for %s" % (failed,)
                                self.failures.append(failed)
                        i = j
                        j = i + rms_per_loop
                        littlelist = lfnlist[i:f(j)]
        ## END rm_listing(self,lfnlist)
        
        def add_metadata(self,pset,lfn,attribs):
                # new method
                # add lfn
                try:
                        self.md.add(lfn,pset)
                except Exception,e:
                        msg = "Caught exception while calling self.md.add(%s). Message was: %s" % (lfn,str(e))   
                        raise LSCfileAddException, msg
                try:
                        for field, val in attribs.iteritems():
                                self.md.setAttribute(lfn,field,val['Value'])
                except Exception,e:
                        msg = "Caught exception while calling self.md.setAttribute() for lfn, \"%s\". Error was: %s" % (lfn,str(e))
                        raise LSCfileAddException, msg
                
        ## END def add_metadata(self)
        
        def get_metadata(self,lfnlist):
                """
                Returns dictionary containing metadata for each member of the lfn list.
                If metadata does not exist for a member, a that lfn entry will point to
                the string "LFN does not exist in database.".
                
                A dicionary element will be of the format,
                        foo_dict[<some_lfn>][<field>]["Value"] = <some_value>
                        
                Unless the entry does not exist in the database, then it will be
                        foo_dict[<missing_lfn>] = "LFN does not exist in database."
                Also, this will result in the lfn + string mentioned above to be appended
                to the failures[].
                        
                Other failures will be reported by adding to the failures[], and raising
                an LSCfileAddException.
                """
                if isinstance(lfnlist,StringType):
                        lfnlist = [lfnlist]
                elif not isinstance(lfnlist,ListType):
                        msg = "get_metadata(var): var must be either a string or a list of strings."
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
                ret = {}
                # see if database handles need to be regenerated
                self.checkHandle()
                try: # catch metadata exceptions.
                        for lfn in lfnlist:
                                ret[lfn] = {}
                                if not self.md.getAttribute(lfn,"md5"):
                                        msg = "LFN does not exist in database."
                                        self.failures.append(msg)
                                        ret[lfn] = msg
                                else:
                                        for item in self.md.getAttributes(lfn).iteritems():
                                                temp = {}
                                                temp["Value"] = item[1]
                                                ret[lfn][item[0]] = temp
                                        
                        # return lfn attribute dictionary
                        return ret
                except MetadataException,e:
                        msg = "get_metadata(%s): Caught MetadataException. Message was: %s" % (lfn,str(e))
                        raise LSCfileAddException, msg
        ## END get_metadata(self,lfnlist)
        
        def get_urls(self,lfnlist,gets_per_loop=50):
                """
                Given a list of LFNs, lfnlist, this will return
                a mapping dictionary of the form,
                
                retdict[LFN] = [PFN1,PFN2,...]
                
                If the LFN does not exist or does not have any PFNs
                associated with it, the PFN list will be of length 0,
                and an entry will be appended to failures[].
                """
                if isinstance(lfnlist,StringType):
                        lfnlist = [lfnlist]
                if not isinstance(lfnlist,ListType):
                        msg = "get_urls: Argument must be a single lfn string, or an lfn list."
                        raise LSCfileAddException, msg
                # now, query RLS for the PFN list, for each lfn in lfnlist
                ret = {}
                
                # grab first "slice" of pfnlist
                i = 0
                j = gets_per_loop
                f = lambda x: x - len(pfnlist) < 0 and x or len(pfnlist)
                littlelist =  lfnlist[i:f(j)]
                while len(littlelist):
                        # see if database handles need to be regenerated
                        self.checkHandle()
                        try:
                                ret.update(self.lrc.pfn_get_bulk(littlelist))
                                if len(ret) != len(littlelist):
                                        if not len(ret):
                                                msg = "get_urls: Warning, could not get any urls for the following, %s " % (littlelist,)
                                                #self.failures.append(msg)
                                                print >>sys.stderr, msg
                                        else:   # find which LFNs do not have pfns mapped to them
                                                temp = ""
                                                for item in littlelist:
                                                        if not ret.has_key(item):
                                                                temp += item
                                                msg = "get_urls: Warning, could not get any urls for the following, %s " % (littlelist,)
                                                #self.failures.append(msg)
                                                print >>sys.stderr, msg
                                                del temp
                                
                        except rlsClient.RlsClientException, e:
                                msg = "get_urls: Caught RlsClientException on \"%s\"" % (littlelist,)
                                self.failures.append(msg)
                                raise LSCfileAddException, msg
                        i = j
                        j = i + gets_per_loop
                        littlelist = pfnlist[i:f(j)]
                        
                return ret
        ## END get_urls(self,lfnlist)
        
        def print_urls(self,lfnlist):
                """
                Prints urls for each entry in lfnlist in the following format.
                
                ${LFN}:\t${PFN1}\n\t\t${PFN2}\n...\t\t${PFNn}\n\n
                
                lfnlist can be either a python list of lfns or a single-lfn string.
                """
                if isinstance(lfnlist,StringType):
                        lfnlist = [lfnlist]
                elif not isinstance(lfnlist,ListType):
                        msg = "get_metadata(var): var must be either a string or a list of strings."
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
                # print urls as they are found
                i = 0
                step = 32
                j = step
                f = lambda x: x - len(lfnlist) < 0 and x or len(lfnlist)
                templist = lfnlist[i:f(j)]
                # see if database handles need to be regenerated
                self.checkHandle()
                try:
                        while(templist):
                                for name in templist:
                                        urls = self.get_urls(name)[name]
                                        mystr = name + ':\n'
                                        # Handle strings given 
                                        #(most likely, lfn doesn't exist messages)
                                        if isinstance(urls,StringType):
                                                mystr += '\t' + urls + '\n'
                                                continue
                                        if len(urls) >= 1:
                                                for pfn in urls:
                                                        mystr += "\t\t%s\n" % (pfn,)
                                                print "%s" % (mystr,)
                                        else:
                                                print "%s\n" % mystr
                                # end for name in templist
                                sys.stdout.flush()
                                i = f(j)
                                j = i + step
                                templist =lfnlist[i:f(j)]
                except LSCfileAddException, e:
                        msg = "print_urls: Generated LSCfileAddException which processing lfn \"%s\". Message was: %s" % (name,str(e))
                        print >>sys.stderr, msg
                        self.failures.append(msg)
                        
        ## END print_urls(self,lfnlist)
        
        def print_metadata(self,lfnlist):
                """
                Prints metadata for each entry in lfnlist in the following format.
                
                ${LFN}:\t${PFN1}\n\t\t${PFN2}\n...\t\t${PFNn}\n\n
                
                lfnlist can be either a python list of lfns or a single-lfn string.
                """
                if isinstance(lfnlist,StringType):
                        lfnlist = [lfnlist]
                elif not isinstance(lfnlist,ListType):
                        msg = "get_metadata(var): var must be either a string or a list of strings."
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
                # print urls as they are found
                mydict = {}
                # Query and print 32 records at a time.
                #   though get_metadata() ultimately queries
                #   one at a time presently.
                i = 0
                step = 32
                j = step
                f = lambda x: x - len(lfnlist) < 0 and x or len(lfnlist)
                templist = lfnlist[i:f(j)]
                # see if database handles need to be regenerated
                self.checkHandle()
                try:
                        
                        while len(templist):
                                for name in templist:
                                        # Convert LFN's attributes
                                        # to a list of tuples????
                                        meta = self.get_metadata(name)[name]
                                        if isinstance(meta,StringType):
                                                mystr = name + ':\n' + '\t' + meta + '\n'
                                                print mystr
                                                continue
                                        myitems = meta.items()
                                        mystr = name + ':\n'
                                        if len(myitems):
                                                for temp in myitems:
                                                        field = str(temp[0])
                                                        value = str(temp[1]["Value"])
                                                        #value = value.rjust(10)
                                                        field = field.rjust(12)
                                                        mystr += "\t\t%s: %s\n" % (field,value)

                                                print "%s" % (mystr,)
                                        else:
                                                continue
                                # end for name in templist
                                sys.stdout.flush()
                                i = f(j)
                                j = i + step
                                templist = lfnlist[i:f(j)]
                        # end while len(templist)        
                except LSCfileAddException, e:
                        msg = "print_metadata: Generated LSCfileAddException which processing lfn \"%s\". Message was: %s" % (name,str(e))
                        self.failures.append(msg)
        ## END print_metadata(self,lfnlist)
        
        def new_pset(self,name,createnew = False):
                """
                Creates a new pset if it does not already exist.
                """
                if not isinstance(name,StringType):
                        msg = "new_pset(): name, %s, is not a string." % (name,)
                        raise LSCfileAddException
                        
                if not self.md.psets().count(name):
                        if createnew:
                                self.md.addPSet(name)
                        else:
                                msg = 'pset, \"%s\", does not exist. User has elected not to create new psets, so the requested pset was not created.' % (name,)
                                raise LSCfileAddException,msg
        ## END new_pset(self,name)
        
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
## END class Publisher(LSCrunMetadataAttr)

class successList(list):
        """
        Simple list class for recording and handling successful actions.
        """
        def __init__(self):
                list.__init__()
                # by default the handler is never called
                self.limit = 0
                self.count = 0
                self.Handler = getattr(self,"default_Handler")
        def append(self,item):
                self.data.append(item)
                self.count += 1
                if self.count - self.limit == 0:
                        self.count = 0
                        self.Handler()
        def default_Handler():
                pass
        def SetHandler(self,object):
                """
                Sets hander to be called when success limit is reached
                """
                self.Handler = object
        def SetLimit(self,newlimit):
                """
                Number of events before calling handler.
                If zero is specified, handler is _NEVER_ called
                """
                self.limit[level] = newlimit
        def clear(self):
                self.data = []
## END class successList(list)

class failureList(list):
        """
        List class for recording and reporting failed actions.
        """
        def __init__(self):
                """
                Inits class. Sets the three default failure levels.
                """
                list.__init__()
                # by default the handler is never called
                self.limit['NOTICE'] = self.count['NOTICE'] = NoticeLimit = 0
                self.limit['WARNING'] = self.count['WARNING'] = WarningLimit = 0
                self.limit['ERROR'] = self.count['ERROR'] = ErrorLimit = 0
                self.Handler['NOTICE'] = getattr(self,"default_Handler")
                self.Handler['WARNING'] = getattr(self,"default_Handler")
                self.Handler['ERROR'] = getattr(self,"default_Handler")
                
        def append(self,item):
                """
                Typical list append, but checks size of
                """
                self.data.append(item)
                # In case of syntax errors in script (i.e. bad levels).
                #   at least the error will be recorded somehow
                #   but this does result in possibly skipping a desired
                #   error handler.
                if not self.limit.has_key(item[0]):
                        self.NewFailLevel(item[0])
                
                self.count[item[0]] += 1
                
                if self.count[item[0]] - self.limit[item[0]] == 0:
                        # reset error counter
                        self.count[item[0]] = 0
                        # call handler
                        self.Handler[item[0]]()
                        
        def default_Handler():
                """
                simply passes
                """
                pass
                
        def NewFailLevel(self,level):
                """
                Creates a new (or redefines an old) failure level.
                """
                self.limit[level] = 0
                self.count[level] = 0
                self.Handler[level] = getattr(self,"default_Handler")
                
        def SetHandler(self,level,object):
                """
                Ties a handler to a given failure level. Creates/redefines new level as needed.
                1st arg is the failure level string
                2nd arg is the object that is to become the handler
                """
                if not self.limit.has_key(level):
                        self.NewFailLevel(level)
                self.Handler[level] = object
                
        def SetLimit(self,level,newlimit):
                """
                Sets number of failure events (at level of first arg) to second arg.
                If the limit is set to 0, handler is _NEVER_ called (unless limit is changed later).
                """
                if not self.limit.has_key(level):
                        self.NewFailLevel(level)
                self.limit[level] = newlimit
        def clear(self):
                self.data = []
## END class failureList(list)
