#!/usr/bin/env python
#
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

import os
import re

def get_all_files_in_range(dirname, starttime, endtime):
    ret = []

    first_four_start = starttime / 100000
    first_four_end   = endtime   / 100000

    for filename in os.listdir(dirname):
        if re.match('.*-[0-9]{4}$', filename):
            dirtime = int(filename[-4:])
            if dirtime >= first_four_start and dirtime <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime)
        elif re.match('.*-[0-9]*-[0-9]*\.xml', filename):
            file_time = int(filename.split('-')[-2])
            if file_time >= (starttime-16) and file_time <= (endtime+16):
                ret.append(os.path.join(dirname,filename))

    return ret



def setup_database(host_and_port):
    global PROGRAM_NAME

    port = 30020
    
    if host_and_port.find(':') < 0:
        host = host_and_port
    else:
        # server and port specified
        host, portString = host_and_port.split(':')
        port = int(portString)


    identity = "/DC=org/DC=doegrids/OU=Services/CN=ldbd/%s" % host

    # open connection to LDBD Server
    client = None

    try:
        client = LDBDClient.LDBDClient(host, port, identity)
    except Exception, e:
        print >>sys.stderr, \
              "Unable to connect to LDBD Server %s:%d" % (host, port)
        if gsiserverutils.checkCredentials():
            print >>sys.stderr, "Got the following error : " + str(e)
            print >>sys.stderr, "Enter '%s --help' for usage" % PROGRAM_NAME
        sys.exit(-1)

    return client


