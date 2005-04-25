# $Id$
# 
# setup for glue

import os
from distutils.core import setup
from setupext import build_metaio

ext_modules = []

packages = [ 'lgen' ]

build_metaio(ext_modules, packages)

setup( name = "lgen",
  version = "0.1",
  author = "Patrick Brady",
  author_email = "patrick@gravity.phys.uwm.edu",
  description = "LSC Graphics Toolkit",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/",
  license = 'See file LICENSE',
  packages = packages,
  ext_modules = ext_modules,
  scripts = [ os.path.join('bin','plotburst'),
              os.path.join('bin','plotsiminspiral'),
              os.path.join('bin','plotinspiral') ],
  data_files = [ ('etc',[
    os.path.join('etc','lgen-user-env.sh'),
    os.path.join('etc','lgen-user-env.csh')
    ] ) ]
  )
