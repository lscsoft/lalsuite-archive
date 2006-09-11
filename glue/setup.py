# $Id$
# 
# setup script for glue

import os
from distutils.core import setup, Extension
from sys import version_info

# Add py_modules to the argument list for setup() if this is at least Python
# 2.3.  Duplicating this function call is ugly but it solves the problem.
# When, finally, nobody we care about it using RH9, we can kill off this mess
# altogether.
if version_info >= (2, 3):
  setup( name = "glue",
    version = "1.8",
    author = "Duncan Brown",
    author_email = "dbrown@ligo.caltech.edu",
    description = "Grid LSC User Engine",
    url = "http://www.lsc-group.phys.uwm.edu/daswg/",
    license = 'See file LICENSE',
    packages = [ 'glue', 'glue.ligolw', 'glue.ligolw.utils' ],
    py_modules = [ 'glue.segfindserver.segments_1_7.segments' ],
    ext_modules = [
      Extension("glue.ligolw.tokenizer", ["glue/ligolw/tokenizer.c"]),
      #Extension("glue.__segments", ["src/segments/segments.c", "src/segments/infinity.c"], include_dirs = ["src/segments"])
    ],
    scripts = [
      os.path.join('bin','LSCdataFind'),
      os.path.join('bin','LSCsegFind'),
      os.path.join('bin','LSCfileAdd'),
      os.path.join('bin','ldbdc'),
      os.path.join('bin','dmtdq_seg_insert'),
      os.path.join('bin','ligolw_add'),
      os.path.join('bin','ligolw_cut'),
      os.path.join('bin','ligolw_burst2mon'),
      os.path.join('bin','ligolw_inspiral2mon'),
      os.path.join('bin','ligolw_print'),
      os.path.join('sbin','ldbdd'),
      os.path.join('sbin','segpagegen'),
      os.path.join('sbin','LSCdqInsert'),
      os.path.join('sbin','publishstatefromfile'),
      os.path.join('sbin','bulkpublishstate'), ],
    data_files = [ ('etc',[
      os.path.join('etc','vdsproperties'),
      os.path.join('etc','glue-user-env.sh'),
      os.path.join('etc','glue-user-env.csh'),
      os.path.join('etc','lscsegfindserver.ini'),
      os.path.join('etc','segpagegen.ini'),
      os.path.join('etc','ldbdserver.ini')
      ] ) ]
    )
else:
  setup( name = "glue",
    version = "1.8",
    author = "Duncan Brown",
    author_email = "dbrown@ligo.caltech.edu",
    description = "Grid LSC User Engine",
    url = "http://www.lsc-group.phys.uwm.edu/daswg/",
    license = 'See file LICENSE',
    packages = [ 'glue', 'glue.ligolw' ],
    ext_modules = [
      Extension("glue.ligolw.tokenizer", ["glue/ligolw/tokenizer.c"]),
      #Extension("glue.__segments", ["src/segments/segments.c", "src/segments/infinity.c"], include_dirs = ["src/segments"])
    ],
    scripts = [ os.path.join('bin','LSCdataFind'),
      os.path.join('bin','LSCdataFind'),
      os.path.join('bin','LSCsegFind'),
      os.path.join('bin','LSCfileAdd'),
      os.path.join('bin','ldbdc'),
      os.path.join('bin','dmtdq_seg_insert'),
      os.path.join('bin','ligolw_add'),
      os.path.join('sbin','ldbdd'),
      os.path.join('sbin','segpagegen'),
      os.path.join('sbin','LSCdqInsert'),
      os.path.join('sbin','publishstatefromfile'),
      os.path.join('sbin','bulkpublishstate'), ],
    data_files = [ ('etc',[
      os.path.join('etc','vdsproperties'),
      os.path.join('etc','glue-user-env.sh'),
      os.path.join('etc','glue-user-env.csh'),
      os.path.join('etc','lscsegfindserver.ini'),
      os.path.join('etc','segpagegen.ini'),
      os.path.join('etc','ldbdserver.ini')
      ] ) ]
    )
