#!/usr/bin/env python
#
# $Id$
# 
# setup for glue

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
