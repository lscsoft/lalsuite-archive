#!/usr/bin/env python
#
# ldbd setup file
#
# $Id$
#
# Copyright (C) 2003 Duncan Brown
# 
# This file is part of the lightweight datapase dumper (ldbd)
#
# ldbd is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# ldbd is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# ldbd; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA
#

from distutils.core import setup

setup( name = "glue",
  version = "0.1",
  author = "Duncan Brown",
  author_email = "duncan@gravity.phys.uwm.edu",
  description = "Grid LSC User Engine",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/",
  license = 'See file LICENSE',
  package_dir = {'': 'lib'},
  py_modules=["glue","gpsTime"]
  )
