# $Id$
#
# Copyright (C) 2009  Larne Pekowsky
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Utilities to create process and process_param entries in a ligolw document
"""

from optparse import OptionParser

import sys
import os
import time
import socket

from glue.ligolw import lsctables

from glue import gpstime
from glue import LDBDClient
from glue import gsiserverutils



__author__ = "Larne Pekowsky <lppekows@physics.syr.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]



def add_process(doc, opts, version = None, cvs_date = None):
    """Takes a ligolw.Document, options returned from
    OptionsParser.parse_args(), and optional version and date as
    strings. Adds a process table to doc, with a row representing the
    current process.  Gets program name, pid, node and user from the
    environment, ifo from the options, and version and cvs date from
    the optional parameters """

    PROGRAM_NAME = sys.argv[0].replace('./','')
    PROGRAM_PID  = os.getpid()
    USER_NAME    = os.getlogin()

    # Do we already have a process table?
    proctable = None
    children  = doc.childNodes

    # This should be the LIGO_LW node
    if len(children) == 1:
        # Find the process table if its there to be found
        proctables = filter(lambda x: str(type(x[0])).find('lsctables.Process') != -1, children[0].childNodes)

        if len(proctables) > 0:
            proctable = proctables[0]

    if proctable == None:
        proctable = lsctables.New(lsctables.ProcessTable)
        doc.childNodes[0].appendChild(proctable)

    proc_id = proctable.get_next_id()

    process                = lsctables.Process()
    process.program        = PROGRAM_NAME
    process.cvs_repository = "lscsoft"

    if version:
        process.version        = version

    if cvs_date:
        process.cvs_entry_time = gpstime.GpsSecondsFromPyUTC(time.mktime(time.strptime(cvs_date, "%Y/%m/%d %H:%M:%S")))
        
    process.comment        = ""
    process.is_online      = 0
    process.node           = socket.gethostbyaddr(socket.gethostname())[0]
    process.username       = USER_NAME
    process.unix_procid    = PROGRAM_PID
    process.start_time     = gpstime.GpsSecondsFromPyUTC(time.time())
    process.end_time       = 0
    process.jobid          = 0
    process.domain         = ""
    process.process_id     = proc_id

    if opts:
        if 'ifo' in opts.__dict__:
            process.ifos = opts.ifo
        elif 'ifos' in opts.__dict__:
            process.ifos = opts.ifos
        else:
            process.ifos = ""
        
    proctable.append(process)

    return proc_id



def add_params(doc, proc_id, opts):
    """Takes a ligolw.Document, a process_id, and options returned
    from OptionsParser.parse_args().  Adds a process_param table to
    doc, populated with the parameters present in the opts"""

    PROGRAM_NAME = sys.argv[0].replace('./','')
    

    # Do we already have a processParams table?
    proc_param_table = None
    children         = doc.childNodes

    # This should be the LIGO_LW node
    if len(children) == 1:
        # Find the process params table if its there to be found
        proc_param_tables = filter(lambda x: str(type(x[0])).find('lsctables.ProcessParams') != -1, children[0].childNodes)

        if len(proc_param_tables) > 0:
            proc_param_table = proc_param_tables[0]

    if proc_param_table == None:
        proc_param_table = lsctables.New(lsctables.ProcessParamsTable)
        doc.childNodes[0].appendChild(proc_param_table)


    for key in opts.__dict__:
        value = opts.__dict__[key]

        if value is not None:
            param            = lsctables.ProcessParams()
            param.program    = PROGRAM_NAME
            param.process_id = proc_id
            param.param      = key
            param.value      = value
            param.type       = "lstring"

            proc_param_table.append(param)


def add_process_and_params(doc, opts = {}, version = None, cvs_date = None):
    proc_id = add_process(doc, opts, version, cvs_date)

    add_params(doc, proc_id, opts)

    return proc_id


