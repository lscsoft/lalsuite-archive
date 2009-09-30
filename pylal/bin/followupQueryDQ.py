#!/usr/bin/env python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
This script queries the available DQ information given a gps time with
an asymetric windows specified at command line. This routine will
return a text string in MoinMoin as a table for inclusion in the
candidate checklist.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__    = 'followupQueryDQ.py'


import optparse
import sys
import os
from pylal import git_version
from pylal.fu_utils import followupDQV

sys.path.append('@PYTHONLIBDIR@')

usage = """usage: %prog [options]"""

parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
#Add all options to setup the query
parser.add_option("-X","--segment-url",action="store",type="string",\
                      metavar="URL",default=None,\
                      help="Using this argument specify a URL the LDBD \
server that you want to query DQ Veto segment information from for\
example ldbd://metaserver.phy.syr.edu:30015")
parser.add_option("-w","--window",action="store",type="string",\
                      metavar="frontWin,backWin",default="30,15",\
                      help="Using this argument you can specify a \
asymetric window around the trigger time for performing DQ queries. \
The two times should be postive values seperated by a comma. If only \
one time is specified then the window is assumed to be symetric. \
Example: --window='100,10'")
parser.add_option("-t","--trigger-time",action="store",type="string", \
                      metavar="GPStime", default=None,\
                      help="Using this argument you can specify the \
GPS time of the trigger to check the data quality flags on.")
parser.add_option("-p","--output-format",action="store",type="string",\
                      metavar="OUTPUT_TYPE",default="MOINMOIN", \
                      help="The valid output options here are \
LIST(python var), MOINMOIN, and HTML.")
parser.add_option("-o","--output-file",action="store",type="string",\
                  metavar="OUTPUTFILE.wiki",default=None,\
                  help="This sets the name of the file to save the DQ\
 result into.")
                  


######################################################################

(opts,args) = parser.parse_args()


server=opts.segment_url
triggerTime=opts.trigger_time
outputType=opts.output_format
outputFile=opts.output_file

if len(opts.window.split(',')) == 1:
    frontWindow=backWindow=opts.window
if len(opts.window.split(',')) == 2:
    frontWindow,backWindow=opts.window.split(',')
    if len(backWindow) == 0:
        backWindow=frontWindow

x=followupDQV(server)
x.fetchInformationDualWindow(triggerTime,frontWindow,backWindow)
result=""
if outputType.upper().strip() == "LIST":
    result=x.generateResultList()
if outputType.upper().strip() == "MOINMOIN":
    result=x.generateMOINMOINTable("DQ")
if outputType.upper().strip() == "HTML":
    result=x.generateHTMLTable("DQ")

if outputFile == None:
    sys.stdout.write("%s"%(result))
else:
    fp=file(outputFile,"w")
    fp.writelines("%s"%(result))
    fp.close()
