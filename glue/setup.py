# $Id$
# 
# setup for glue

import os
from distutils.core import setup

setup( name = "glue",
  version = "0.1",
  author = "Duncan Brown",
  author_email = "dbrown@ligo.caltech.edu",
  description = "Grid LSC User Engine",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/",
  license = 'See file LICENSE',
  packages = [ 'glue' ],
  scripts = [ os.path.join('bin','LSCdataFind'),
    os.path.join('bin','LSCdataFind'),
    os.path.join('bin','LSCsegFind'),
    os.path.join('bin','WebsegFind'),
    os.path.join('bin','LSCfileAdd'),
    os.path.join('bin','ldbdc'),
    os.path.join('sbin','ldbdd'),
    os.path.join('sbin','LSCsegFindServer') ],
  data_files = [ ('etc',[
    os.path.join('etc','glue-user-env.sh'),
    os.path.join('etc','glue-user-env.csh'),
    os.path.join('etc','lscsegfindserver.ini')
    ] ) ]
  )
