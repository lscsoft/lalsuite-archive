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
import urlparse
import md5
import time
import sys
import re
from types import *

# LDR-specific modules
import LDRUtil
import RLS
import LDRMetadataCatalog
import rlsClient

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


class Publisher(object):
        """
        Class for publishing Engineering and/or science run data at
        at the LIGO observatories.
        """
        def __init__(self):
                #LSCrunMetadataAttr.__init__(self)
                # set some default flags
                self.PRESERVE_METADATA  = 0x01
                self.OVERWRITE_METADATA = 0x02
                
                # Initialize LDR stuff
                self.config = LDRUtil.getConfig("local", "DEFAULT")
                self.metadata = LDRMetadataCatalog.LDRMetadataCatalog(self.config)
                self.rliurl = self.config["rli"]
                self.lrc = RLS.getCatalog(self.rliurl)
                self.gsiurl = self.config.get("Storage", "hostname")
                
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
                Regenerates handle to rls... in case of rlsClient.RlsClientException
                during publishing.
                """
                self.metadata = LDRMetadataCatalog.LDRMetadataCatalog(self.config)
                self.lrc = RLS.getCatalog(self.rliurl)
                print >>sys.stderr, "Regenerated the lrc and metadata handle at %s" % time.asctime()
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
                                    --> ["metadata"] --> ["<attribute>"]["Value"]
                """
                
                ### Need to perform some basic checking of the
                ### data list structure...                
                
                # to be implemented
                
                # see if database handles need to be regenerated
                self.checkHandle()
                
                #### Attempt to publish metadata
                metadata_failure = False # though, any metadata failure should 
                                         # most likely warrant raising an exception
                for data in datalist:
                        lfn = data["name"]
                        
                        metaexists = self.metadata.exists(lfn)

                        if metaexists & (flags & self.PRESERVE_METADATA):
                                pass
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
                                        self.add_metadata(lfn,data["metadata"])
                                except Exception, e:
                                        msg = "Caught exception while publishing metadata \"%s\"! Error was: %s" % (lfn,str(e))
                                        metadata_failure = True
                                        raise LSCfileAddException, msg
                                        
                #### Attempt to publish RLS entries
                rls_failure = False
                # Map each pfn one at a time
                for data in datalist:
                        lfn = data["name"]
                        
                        # Check for LFN existance
                        lfn_exists = False
                        try:
                                lfn_exists = self.lrc.lfn_exists(lfn) 
                        except rlsClient.RlsClientException, e:
                                msg = "Caught excfeption from self.lrc.lfn_exists(%s). Error message was: %s" % (lfn,str(e))
                                msg += "\nRegenerating handle to LDR and RLS and trying again."
                                print >>sys.stderr,msg
                                self.regenHandle()
                                try:
                                        lfn_exists = self.lrc.lfn_exists(lfn)
                                except rlsClient.RlsClientException, e:
                                        msg = "Caught exception from self.lrc.lfn_exists(%s) Again! Error message was: %s" % (lfn,str(e))
                                        print >>sys.stderr,msg
                                        raise LSCfileAddException, msg
                        
                        pfnlist = data["urls"]
                        for pfn in pfnlist:
                                # Now see if pfn exists
                                pfn_exists = False
                                if lfn_exists:
                                        try:
                                                pfn_exists = self.lrc.pfn_exists(pfn) 
                                        except rlsClient.RlsClientException, e:
                                                msg = "Caught excfeption from self.lrc.lfn_exists(%s). Error message was: %s" % (lfn,str(e))
                                                msg += "\nRegenerating handle to LDR and RLS and trying again."
                                                print >>sys.stderr,msg
                                                self.regenHandle()
                                                try:
                                                        pfn_exists = self.lrc.pfn_exists(pfn)
                                                except rlsClient.RlsClientException, e:
                                                        msg = "Caught exception from self.lrc.lfn_exists(%s) Again! Error message was: %s" % (lfn,str(e))
                                                        print >>sys.stderr,msg
                                                        raise LSCfileAddException, msg
                                if pfn_exists and lfn_exists:
                                        msg = "LFN \"%s\" and PFN \"%s\" already exists in database. This script will not perform the pfn insertion operation." % (lfn,pfn)
                                        print >>sys.stderr, msg
                                        self.failures.add((lfn,msg))
                                else:
                                        try:
                                                self.lrc.add(lfn,pfn)
                                        except rlsClient.RlsClientException, e:
                                                msg = "Caught Exception from self.lrc.add(%s,%s). \
                                                Error message was \"%s\"" % (lfn,pfn,e)
                                                msg += "\nRegenerating handle to LDR and RLS and trying again."
                                                print >>sys.stderr, msg
                                                self.regenHandle()
                                                try:
                                                        self.lrc.add(lfn,pfn)
                                                except rlsClient.RlsClientException, e:
                                                        msg = "Caught Exception from self.lrc.add(lfn,pfn[i]). \
                                                        Again!. Error message from rlsClient \"%s\"" % (e,)
                                                        print >>sys.stderr, msg
                                                        self.failures.append((lfn,msg))
                                                        rls_failure = True
                                                
                ### End RLS entry publishing
                
                # if all DB additions and checks were successful, then
                if not rls_failure and not metadata_failure:
                        self.successes.append(lfn)
                        
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
                                if self.metadata.exists(lfn):
                                       self.metadata.delete(lfn)
                        except Exception,e:
                                msg = "completely_remove(): Caught an exception while calling either metadata.exists() or metadtaa.delete() on \"%s\". Error was: %s" % (lfn,str(e))
                                raise LSCfileAddException, msg
       ## END def completely_remove(self)
        
        def mv_url(self,source,destination):
                """
                Used like UNIX /bin/mv on RLS entries.
                
                usage:
                        foo.mv_url("<Original PFN>","<Destination PFN>")
                        
                Upon error, appends error message to failures[] and raises
                an LSCfileAddException. Otherwise adds a message to successes[].
                """
                # see if database handles need to be regenerated
                self.checkHandle()
                try:
                        if not self.lrc.pfn_exists(source):
                                msg = "mv_url(%s,%s): Source PFN, %s, does not exist. Nothing done." % (source,destination,source)
                                self.failures.append(msg)
                                raise LSCfileAddException, msg
                        else:
                                lfn = self.lrc.get_lfn(source)[0]
                                self.lrc.add(lfn,destination)
                                # assume new entry exists, remove source pfn
                                self.lrc.delete(lfn,source)
                                msg = "Moved \"%s\" to \"%s\"." % (source,destination)
                                successes.append(msg)
                                return
                except rlsClient.RlsClientException, e:
                        msg = "mv_url(%s,%s): Caught RLS exception. Error message received was: %s" % (source,destination,str(e))
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
        ## END mv_url(self,source,destination)
        
        def rm_url(self,pfnlist):
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
                        
                # see if database handles need to be regenerated
                self.checkHandle()
                try:
                        for name in pfnlist:
                                if not self.lrc.pfn_exists(name):
                                        msg = "rm_url(): URL %s, does not exist in RLS database. Nothing done." % (name,)
                                        self.failures.append(msg)
                                        #raise LSCfileAddException, msg
                                else:
                                        lfn = self.lrc.get_lfn(name)[0]
                                        self.lrc.delete(lfn,name)
                                        msg = "Deleted URL %s. " % (name,)
                                        self.successes.append(msg)
                except rlsClient.RlsClientException, e:
                        msg = "rm_url(%s): Caught RLS exception. Message was: %s" % (name,str(e))
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
        ## END rm_url(self,pfnlist)
        
        def rm_listing(self,lfnlist):
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
                # see if database handles need to be regenerated
                self.checkHandle()
                try:
                        for name in lfnlist:
                                if not self.lrc.lfn_exists(name):
                                        msg = "rm_listing(%s,%s): URL %s, does not exist in RLS database. Nothing done." % (source,destination,source)
                                        self.failures.append(msg)
                                        raise LSCfileAddException, msg
                                else:
                                        pfns = self.lrc.get_pfn(name)
                                        for element in pfns:
                                                self.lrc.delete(name,element)
                                        # Check again, in case this code is buggy!
                                        if not self.lrc.lfn_exists(name):
                                                msg = "rm_listing: Deleted listing for LFN %s." % (name,)
                                                successes.append(msg)
                                                return
                                        else:
                                                msg = "Failed to remove listing for LFN %s. This code is buggy!" % (name,)
                                                self.failures.append(msg)
                                                raise LSCfileAddException, msg
                except rlsClient.RlsClientException, e:
                        msg = "rm_listing(%s): Caught RLS exception. Message was: %s" % (name,str(e))
                        self.failures.append(msg)
                        raise LSCfileAddException, msg
        ## END rm_listing(self,lfnlist)
        
        def add_metadata(self,lfn,attribs):
                try:
                        self.metadata.add(lfn)
                except LDRMetadataCatalogException,e:
                        msg = "Caught LDRMetadataCatalogException while calling self.metadata.add(%s). Message was: %s" \
                        % (lfn,str(e))   
                        raise LSCfileAddException
                try:
                        for field, val in attribs.iteritems():
                                self.metadata.set_attr(lfn,field,val['Value'])
                except LDRMetadataCatalogException,e:
                        msg = "Caught LDRMetadataCatalogException while calling self.metadata.set_attr(%s,%s,%s). Message was: %s" \
                        % (lfn,field,val['Value'],str(e))
                        raise LSCfileAddException
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
                                if not self.metadata.exists(lfn):
                                        msg = "LFN does not exist in database."
                                        self.failures.append(msg)
                                        ret[lfn] = msg
                                else:
                                        for item in self.metadata.get_attrs(lfn).iteritems():
                                                temp = {}
                                                temp["Value"] = item[1]
                                                ret[lfn][item[0]] = temp
                                        
                        # return lfn attribute dictionary
                        return ret
                except LDRMetadataCatalog.LDRMetadataCatalogException,e:
                        msg = "get_metadata(%s): Caught LDRMetadataCatalogException. Message was: %s" % (lfn,str(e))
                        raise LSCfileAddException, msg
        ## END get_metadata(self,lfnlist)
        
        def get_urls(self,lfnlist):
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
                # see if database handles need to be regenerated
                self.checkHandle()
                for name in lfnlist:
                        try:
                                pfns = self.lrc.get_pfn(name)
                                ret[name] = pfns
                                if len(pfns) == 0:
                                        msg = "get_urls: could not find any URLs that map to LFN \"%s\"." % (name,)
                                        self.failures.append(msg)
                                return ret
                        except rlsClient.RlsClientException, e:
                                msg = "get_urls: Caught RlsClientException on \"%s\"" % (name,)
                                self.failures.append(msg)
                                raise LSCfileAddException, msg
                
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
                                templist = lfnlist[i:f(j)]
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
