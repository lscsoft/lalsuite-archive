#!/usr/bin/env python
"""
Tomoki Isogai (isogait@carleton.edu)
function save_dict is adapted from io.py by Nickolas Fotopolos and modified.

This program gets all the KW triggers that fall into a specified segment list
above a certain threshold for the channel specified.
Supported ifos are H1, H2, L1, and V1.
Only applicable for S5 so far since triggers location is hard coded.
Also, from the same reason, this needs to be ran at cit

Supported channels are listed in LIGO_channel_list.txt and 
VIRGO_channel_list.txt

This program supports to save in pickle, pickle.gz, txt, txt.gz and mat file.
Requires scipy if saving in mat file.

$Id$
"""

from __future__ import division
import optparse
import os.path as p
import os, sys
from glue.segments import segment, segmentlist
from glue import segmentsUtils
import cPickle

__author__ = "Tomoki Isogai <isogait@carleton.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,\
                                   version="$Id: get_KW,v 1.00 2008/6/28")
    parser.add_option("-i", "--ifo", help="ifo you want KW from")
    parser.add_option("-a", "--segment_file",\
                      help="file which contains analyzed time segments")
    parser.add_option("-c", "--channel_name", \
      help="channel names you want KW from. Use comma to specify more than one")
    parser.add_option("-m", "--min_thresh", help="minimum threshold for KW")
    parser.add_option("-s", "--save_option", action="store_true", \
        default=False, help="save it in file")
    parser.add_option("-o", "--outname", help="output name in case of saving")
    parser.add_option("-v", "--verbose", action="store_true", default=False, \
                      help="run verbosely")
    
    opts, args = parser.parse_args()
    
    # check if necessary input exists
    if opts.ifo is None:
        print >>sys.stderr, "--ifo is a required parameter"
        sys.exit(1)
        
    if opts.segment_file is None:
        print >>sys.stderr, "--segment_file is a required parameter"
        sys.exit(1)
    if not os.path.isfile(opts.segment_file):
        print >>sys.stderr, "segment_file not found"
        sys.exit(1)
    
    if opts.channel_name is None:
        print >>sys.stderr, "--channel_file is a required parameter"
        sys.exit(1)
    
    if opts.min_thresh is None:
        print >>sys.stderr, "--min_thresh is a required parameter"
        sys.exit(1)
        
    if opts.save_option is True and opts.outname is None:
        print >>sys.stderr, "--outname is a required parameter if saving"
        sys.exit(1)
        
    # show parameters
    if opts.verbose:
        print
        print "********************** PARAMETERS ******************************"
        print 'segment file:'; print opts.segment_file; 
        print 'channel:'; print opts.channel_name; 
        print 'minimum threshold:'; print opts.min_thresh;
        print 'save:'; print opts.save_option;
        if opts.save_option is True:
            print 'output file name:'; print opts.outname;
        print
        
    return opts

def get_trigs(ifo,channel,segfile,min_thresh,verbose):
    """
    get time and snr of KW triggers for a particular channel that occured  in 
    the specified segment list and above specified snr threshold
    ifo has to be either H1, H2, L1, V1
    """
    
    # read segment file
    seg_list =\
          segmentsUtils.fromsegwizard(open(segfile),coltype=float,strict=False)
    seg_list.coalesce()
    if seg_list==[]:
        print >> sys.stderr, """
        Error: file contains no segments or glue.segmentsUtils.fromsegwizard is
               not reading segments correctly. Please check the seg file. 
               (possibly comments in the file is causing this)
        """
        sys.exit(1)
        
    times=[]; snr=[]
    
    # deal ligo and virgo case separately, since they are stored in a different
    # way
    
    ## virgo case
    if ifo=="V1":
        data_loc="/archive/home/mabizoua/public_html/VSR1/KW/cleandata_DQ_cat1/"
        try:
            # check local first
            # if not there, use the Internet to get KW triggers
            if p.exists(data_loc):
                trigs_loc=os.path.join(data_loc,"%s.dat"%channel)
                if verbose: print "retreiving data from %s ..."%trigs_loc
                trig_file=open(trigs_loc)
            else:
                sys.exit(1)
        except (IOError):
            print >>sys.stderr, "Error: file for channel %s not found" % channel 
            sys.exit(1)
        for line in trig_file:
            trig=line.split()
            # exclude comment part
            if trig[0][0]!="#":
                # check if that KW event is in analyzed segment and also
                # check if its snr is bigger than specified minimum snr
                t = float(trig[0]) # time
                s = float(trig[7]) # snr
                if t in seg_list and s > min_thresh and s != float('inf'):
                    times.append(t)
                    snr.append(s)
    ## ligo case
    else:
        trigs_loc="/archive/home/lindy/public_html/triggers/s5/"
        ## figure out which days we need from segment file
        day=86400
        
        # express all the possible days in segments
        alldays=segmentlist([segment(startTime,startTime+day) for startTime\
            in range(815068800,875376000,day)])
        
        analyze_days=segmentlist() # this will be the list of days we need
        
        start=False; end=False
        for d in alldays:
            # if the very fist time in the seg list is in the day, that's the 
            # first day to be analyzed
            if seg_list[0][0] in d: start=True
            
            # if the very last time in the seg list is in the day, that's the
            # last day to be analyzed
            if start and not end: analyze_days.append(d)
            if seg_list[-1][1] in d: end=True
            
        # convert segments to string list
        analyze_days=segmentsUtils.to_range_strings(analyze_days)
        ## retreive the data
        for day in analyze_days:
            try:
                # check local first
                # if not there, use the Internet to get KW triggers
                if p.exists(trigs_loc):
                    data_loc=trigs_loc+"_".join(day.split(":"))+"/s5_"+channel+".trg"
                    if verbose: print "retreiving data from %s..."%data_loc
                    trig_file=open(data_loc)
                else:
                    sys.exit(1)
            except (IOError):
                # some channel are not in all the daily folder
                print >>sys.stderr, "channel %s not found in day %s, ignoring"\
                                     %(channel, "-".join(day.split(":")))
            
                        
                        
    ## sanity check
    assert len(times)==len(snr)
    
    return times, snr
    
def save_dict(filename, data_dict, file_handle=None):
    """
    Originally written by Nickolas Fotopoulos, modified.
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
            print >>sys.stderr, "Cannot write to a given file_handle with",\
                ".mat format.  Attempting to ignore."
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
        for chan in data_dict.keys():
            line.append("# "+chan)
            line.append("#time"+" "*15+"#snr")
            for trig in zip(data_dict[chan][0],data_dict[chan][1]):
                line.append("%f"%(trig[0])+'    '+"%f"%(trig[1]))
        output = '\n'.join(line)
    else:
        raise ValueError, "Unrecognized file extension"
    
    ## Write output
    file_handle.write(output)
    
def main(opts):
    # KW trigs is a dictionary whose key is the channel name and item is a list 
    # of two list: time and snr
    KWtrigs={}
    
    # make a list of channels from input
    # deal with white space
    channels=[chan.strip() for chan in opts.channel_name.split(",")]
    for c in channels:
        if opts.verbose: print "getting data for",c 
        KWtrigs[c]=get_trigs(opts.ifo,c,opts.segment_file,int(opts.min_thresh),\
            opts.verbose)
    
    # save the result
    if opts.verbose: print "saving the data in %s..." % opts.outname
    if opts.save_option==True:
        save_dict(opts.outname,KWtrigs)
    
if __name__=="__main__":
    # parse commandline
    opts = parse_commandline()
    # do the job
    main(opts)                
