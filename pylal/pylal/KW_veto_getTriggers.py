#!/usr/bin/env python
"""
Tomoki Isogai (isogait@carleton.edu)
function save_dict is adapted from io.py by Nickolas Fotopolos and modified.

This program gets time and snr data for inspiral triggers from a file, filter 
triggers by the minimum snr and segment file, and saves the data in a file.
This code tries to deal with txt, dat, xml, mat and pickle file.
See note.txt for more detail about supported format.
To read/save mat file, scipy must be in the path

$Id$
"""
from __future__ import division
import cPickle
import optparse
import os.path as p
import sys, os
from numpy import *
from glue import segmentsUtils
from glue.segments import segment, segmentlist

__author__ = "Tomoki Isogai <isogait@carleton.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,\
                                 version="$Id: get_triggers.py,v 1.10 2008/7/1")
    parser.add_option("-t", "--trigger_file",\
                      help="file which contains triggers")
    parser.add_option("-s", "--segment_file",
                      help="file which contains segments by which "+\
                                                      "triggers are filtered")
    parser.add_option("-m", "--min_thresh", default=8, type="float",\
                      help="filter triggers below this minimum snr value")
    parser.add_option("-o", "--outname", default='output.pickle',\
                      help="output file name")
    parser.add_option("-v", "--verbose", action="store_true",\
                      default=False, help="run verbosely")
    
    opts, args = parser.parse_args()
    
    # check if necessary input exists
    for o in ("trigger_file","segment_file"):
        if getattr(opts,o) is None:
            print >>sys.stderr, "Error: --%s is a required parameter"%o
            sys.exit(1)
        if not os.path.isfile(getattr(opts,o)):
            print >>sys.stderr, "Error: %s not found"%o
            sys.exit(1)
        
        
    # show parameters
    if opts.verbose:
        print "running get_triggers..."
        print
        print "********************** PARAMETERS ******************************"
        print 'trigger file:'; print opts.trigger_file; 
        print 'segment file:'; print opts.segment_file;
        print 'minimum snr:'; print opts.min_thresh; 
        print 'output file name:'; print opts.outname;
        print
    
    return opts

def get_trigs_txt(trigger_file,min_thresh,verbose):
    """
    read triggers data from text file
    """
    trigs_file=open(trigger_file)
    # read the first line after header and find out the number of columns
    comment = True # avoid header
    while comment:
        line=trigs_file.readline().split()
        if line[0][0] != "#" and line[0][0] != "%":
            columnNum = len(line)
            comment = False
    if verbose: print "number of columns in the file:"; print columnNum
      
    # if the number of the columns is 2-4, read the first line as time and the 
    # second line as snr
    if columnNum > 1 and columnNum < 5:
        t_column=0; s_column=1
        
    # if the number of the columns is 5, read the first line as time and the 
    # 5th line as snr (virgo's trigger file)
    elif columnNum == 5:
        t_column=0; s_column=4;
        
    # otherwise unsupported
    else:
        print >>sys.stderr, "Error: unrecognized file format, please read "\
                            "note.txt and modify the file to a supported format"
        sys.exit(1)
        
    if verbose:
        print "column # used for times:"; print t_column+1
        print "column # used for snr:"; print s_column+1; print
        
    ## get times and snrs
    times=[]; snr=[]
    for line in trigs_file:
        trig=line.split()
        if trig[0][0]!="#" and trig[0][0]!="%": # exclude comment lines
            s = float(trig[s_column])
            if s >= min_thresh and s != float('inf'): # inf is bad
                times.append(float(trig[t_column]))
                snr.append(s)
    
    # check if first entry of time list is bigger than 800000000 to reduce the 
    # possibility to read the wrong column
    # (sometimes first colume is row number...)
    return check_time_list(times,snr) 

def get_trigs_xml(trigger_file,min_thresh,verbose):
    """
    read triggers data from xml file that has ligo_lw format.
    many lines in this function are adapted from ligolw_print by Kipp Cannon and 
    inspiral_veto_evaluator by Kipp Cannon, Chad Hanna and Jake Slutsky, and 
    modified to suit the purpose.
    """
    from glue.lal import CacheEntry
    from glue.ligolw import ligolw
    from glue.ligolw import table
    from glue.ligolw import lsctables
    from glue.ligolw import utils
    # speed hacks
    # replace Glue's pure Python LIGOTimeGPS class with pyLAL's C version
    #from pylal.xlal.date import LIGOTimeGPS
    #lsctables.LIGOTimeGPS = LIGOTimeGPS
    
    # Enable column interning to save memory
    table.RowBuilder = table.InterningRowBuilder
    
    # for now, hardcode the names
    urls = [trigger_file]
    tables = ['sngl_inspiral']
    columns = ['end_time','end_time_ns','snr']
    
    # don't parse other tables so as to improve parsing speed and reduce memory 
    # requirements.  Because we do this, we can assume later that we should 
    # print all the tables that can be found in the document.
    
    def ContentHandler(xmldoc):
        return ligolw.PartialLIGOLWContentHandler(xmldoc, lambda name, attrs:\
                                (name == ligolw.Table.tagName) and\
                                (table.StripTableName(attrs["Name"]) in tables))
	
    utils.ContentHandler = ContentHandler
    
    # loop over and get the triggers whose snr is above minimum snr specified
    snr = []; times = []
    for url in urls:
        # some xml file uses the old event_id int_8s
        # try the new version and if failed try int_8s for event_id
        for i in range(2):
            try:
                xmldoc = utils.load_url(url, verbose = verbose,\
                                          gz = (url or "stdin").endswith(".gz"))
                break
            except (ligolw.ElementError):
                print 'old version'
                lsctables.SnglInspiralTable.validcolumns["event_id"] = "int_8s"
        table.InterningRowBuilder.strings.clear()
        for table_elem in xmldoc.getElements(lambda e:\
                                           (e.tagName == ligolw.Table.tagName)):
            for n, row in enumerate(table_elem):
                s = row.snr
                if s >= min_thresh and s != float('inf'): # inf is bad
                    times.append(float(row.get_end()))
                    snr.append(s)
        xmldoc.unlink()
    return times,snr

def get_trigs_mat(trigger_file,min_thresh,verbose):
    """
    read triggers data from mat file that contain lists carrying the trigger
    information.
    ! Don't use it! 
    It's possible that code can grab the wrong list as times or snr and check is
    not robust.
    Only reason it's here is that the old code was matlab code and I have many
    triggers in .mat format.
    """
    # scipy is necessary to read mat file
    try:
        import scipy.io
    except (ImportError):
        print >>sys.stderr, "Error: need scipy to read a mat file."
        raise
    try:
        # load and see what's inside
        mat_contents=scipy.io.loadmat(trigger_file)
        # exclude __version__, __header__, __globals__, etc.
        list_names = [l for l in mat_contents.keys() if l[:2]!="__"]
        times=[]; snr=[]
        # case 1: only one list
        if len(list_names) == 1:
            if verbose:
                    print "list name used for getting triggers info:"
                    print list_names[0]; print
            trigs = mat_contents[list_names[0]].tolist()
            # assume there are more than three triggers in the file
            # issue an error otherwise
            # find the dimension of the matrix
            n = len(trigs); m = len(trigs[0])
            # case 1: 2 * n matrix
            if n == 2 and m > 2:
                if verbose: print "2 * n matrix found, parsing..."
                # figure out which list is time
                # (check the first entry of each and assume that the list with 
                # its entry > 800000000 is time list)
                t, s = check_time_list(trigs[0],trigs[1])
                # loop over and filter triggers below the minimum snr defined
                for i in range(m):
                    if s[i] >= min_thresh and s[i]!=float('inf'):
                        times.append(t[i])
                        snr.append(s[i])
            # case 2: n * 2 matrix
            elif n > 2 and m == 2:
                if verbose: print "n * 2 matrix found, parsing..."
                # transpose
                trigs = map(lambda *row: list(row), *trigs)
                # figure out which list is time
                t, s = check_time_list(trigs[0],trigs[1])
                # loop over and filter triggers below the minimum threshold
                for i in range(n):
                    if s[i] >= min_thresh and s[i]!=float('inf'):
                        times.append(t[i])
                        snr.append(s[i])
            # issue an error for other cases
            else:
                print >>sys.stderr, "Error: only 2 * n or n * 2 matrix is"\
                            "supported for mat file. Please read note.txt and"\
                            "modify the file to a supported format"
                sys.exit(1)
        # case 2: 2 lists    
        elif len(list_names) == 2:
            if verbose:
                    print "list names used for getting triggers info:"
                    print list_names[0], list_names[1]; print
            # find out which list is time list
            # sanity check for list containing another list is done in
            # check_time_list()
            times, snr = check_time_list(mat_contents[list_names[0]].tolist(),\
                                         mat_contents[list_names[1]].tolist())
        # case 3: more than two lists 
        # the program can't handle so ask user to modify the file
        else:
            print >>sys.stderr, "Error: unrecognized format: please read",\
                            "note.txt and modify the file to a supported format"
            sys.exit(1)
    except: # improve as errors are found
        raise
    
    return times, snr

def get_trigs_pickle(trigger_file,min_thresh,verbose):
    """
    Get triggers info from pickle file
    Assumption is that the info stored in the pickle is a dictionary whose key
    named times corresponds to times and key named snr corresponds to snr
    """
    trigs=cPickle.load(open(trigger_file))
    # check if it's dictionary
    if type(trigs)==type({}):
        t_name='times';s_name='snr'
        if verbose:
            print "dictionary key used for times:"; print t_name;
            print "dictionary key used for snr:"; print s_name; print
        try:
            # filter out low snr
            times=[]; snr=[]
            for trig in zip(trigs[t_name],trigs[s_name]):
                # in case entries are str, float() it
                s = float(trig[1]) # snr
                if s >= min_thresh and s!=float('inf'):
                    times.append(float(trig[0]))
                    snr.append(s)
        except (KeyError):
            print >>sys.stderr, """
            Error: dictionary must have keys 'times' and 'snr' for times and 
                   snr: please read note.txt and modify the file to a supported
                   format
            """
            raise
    # issue an error if it's not dictionary
    else:
        print >>sys.stderr, """
        Error: unrecognized file format (no dictionary found):
               please read note.txt and modify the file to a supported format
               """
        sys.exit(1)
    return times,snr

def check_time_list(list1,list2):
    """
    two functionality: 
    1. error checking if the program got the right time list
    2. figure out which of two list is time list
    return time list, snr list
    assumption: time is always > 800000000
    """
    # make sure list contains triggers
    if list1 == []:
        print >>sys.stderr, "Error: no triggers"
        sys.exit(1)
    # make sure element is float, exclude list of list
    if type(list1[0])!=type(float(1)) or type(list2[0])!=type(float(1)):
        print >>sys.stderr, """
        Error: unrecognized file format (possibly too many dimensions): 
               please read note.txt and modify the file to a supported format
        """
        sys.exit(1)
        
    # go on until one of the list element is > 800000000 but not the other
    # expectation: not so many SNRs are > 800000000
    for i in xrange(len(list1)):
        if list1[i]>800000000 and list2[i]<800000000: 
            return list1, list2
        elif list2[0]>800000000 and list1[i]<800000000:
            return list2, list1

def check_triggers(times, snr, segfile, verbose):
    """
    This function checks if all triggers falls in the segment list.
    If not, exclude those triggers outside of segment list and issue a message.
    """
    ## read in segment data
    seg_list =\
           segmentsUtils.fromsegwizard(open(segfile),coltype=float,strict=False)
    if seg_list==[]:
        print >> sys.stderr, """
        Error: file contains no segments or glue.segmentsUtils.fromsegwizard is
               not reading segments correctly. Please check the seg file. 
               (possibly comments in the file is causing this)
        """
        sys.exit(1)
        
    ## check if triggers are in the segments
    # trigs[0] is the time list and trigs[1] is the snr list of the triggers
    time_list=[]; snr_list=[]
    outside=False
    for i in xrange(len(times)):
        if times[i] in seg_list:
            time_list.append(times[i]); snr_list.append(snr[i])
        else:
            outside=True
    if outside:
        print >> sys.stderr, """
        Warning: Some of the triggers are outside of the segment list.
                 Unless intentional (using DQ flags etc.), make sure you
                 are using the right segment list.
                 Ignoring..."""
        
    # check at least one trigger remains
    if len(time_list)==0:
        print >> sys.stderr, """
        Error: No triggers. Check trigger file and segment file.
        """
        sys.exit(1)
        
    return time_list, snr_list

def save_dict(filename, data_dict, file_handle=None):
    """
    Adapted from io.py by Nickolas Fotopoulos, modified.
    Dump a dictionary to file.  Do something intelligent with the filename
    extension.  Supported file extensions are:
    * .pickle - Python pickle file (dictionary serialized unchanged)
    * .pickle.gz - gzipped pickle (dictionary serialized unchanged)
    * .mat - Matlab v4 file (dictionary keys become variable names; requires
             Scipy)
    * .txt - ASCII text 
    * .txt.gz - gzipped ASCII text
    For all formats but .mat, we can write to file_handle directly.
    """
    if filename == '':
        raise ValueError, "Empty filename"
    
    ext = p.splitext(filename)[-1]
    
    ## .mat is a special case, unfortunately
    if ext == '.mat':
        if file_handle is not None:
            print >>sys.stderr, "Warning: Cannot write to a given file_handle",\
                "with .mat format.  Attempting to ignore."
        import scipy.io
        scipy.io.savemat(filename, data_dict)
        return
    
    ## Set up file_handle
    if file_handle is None:
        file_handle = file(filename, 'wb')
    
    # For gz files, bind file_handle to a gzip file and find the new extension
    if ext == '.gz':
        import gzip
        file_handle = gzip.GzipFile(fileobj=file_handle, mode='wb')
        ext = p.splitext(filename[:-len(ext)])[-1]
    
    ## Prepare output
    if ext == '.pickle':
        import cPickle
        output = cPickle.dumps(data_dict, -1) # -1 means newest format
    elif ext == '.txt':
        line=[]
        line.append("#time"+" "*15+"#snr")
        for trig in zip(data_dict['times'],data_dict['snr']):
            # align
            line.append("%f"%trig[0]+" "*(20-len("%f"%trig[0]))+"%f"%trig[1])
        output = '\n'.join(line)
    else:
        raise ValueError, "Unrecognized file extension"
    
    ## Write output
    file_handle.write(output)

def main(trigger_file,min_thresh,segment_file,verbose):
    # find out file format and read in the data
    ext=os.path.splitext(trigger_file)[-1]
    if ext=='.pickle':
        if verbose: print "getting triggers from pickle file..."
        times,snr=get_trigs_pickle(trigger_file,min_thresh,verbose)
    elif ext == '.txt' or ext == '.dat':
        if verbose: print "getting triggers from txt/dat file..."
        times,snr=get_trigs_txt(trigger_file,min_thresh,verbose)
    elif ext == '.xml':
        if verbose: print "getting triggers from xml file..."
        times,snr=get_trigs_xml(trigger_file,min_thresh,verbose)
    elif ext == '.mat':
        if verbose: print "getting triggers from mat file..."
        times,snr=get_trigs_mat(trigger_file,min_thresh,verbose)
    else:
        print >>sys.stderr, """
        Error: unrecognized file format: please read note.txt and modify the 
               file to a supported format
        """
        sys.exit(1)
        
    if verbose: print "times and snr for triggers retrieved!"
        
    # filter triggers by segment list
    times, snr = check_triggers(times,snr,segment_file,verbose)
    
    ## sanity check
    assert len(times)==len(snr)
        
    return [times, snr]
    
if __name__=="__main__":
    # parse commandline
    opts = parse_commandline()
    trigs = main(opts.trigger_file, opts.min_thresh,\
                                        opts.segment_file, opts.verbose)
    
    # convert the data to dictionary to save
    trigs={"times": trigs[0], "snr": trigs[1]}
    
    #save
    if opts.verbose: print "saving in '%s'..." % opts.outname
    save_dict(opts.outname,trigs)
    if opts.verbose: print "get_triggers done!"

