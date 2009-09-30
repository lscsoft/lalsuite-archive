# $Id$
# 
# setup script for glue

import os, sys

try:
  from sys import version_info
except:
  print >> sys.stderr, "Unable to determine the python version"
  print >> sys.stderr, "Please check that your python version is >= 2.4"
  sys.exit(1)

if version_info < (2, 4):
  print >> sys.stderr, "Your python version " + str(version_info) + " appears to be less than 2.4"
  print >> sys.stderr, "Please check that your python version is >= 2.4"
  print >> sys.stderr, "Glue requires at least version 2.4"
  sys.exit(1)

from misc import determine_git_version
from distutils.core import setup, Extension
from distutils.command import build_py
from distutils.command import install
from distutils.command import sdist
from distutils.command import clean
from distutils import log

ver = "1.28.1"

def remove_root(path,root):
  if root:
    return os.path.normpath(path).replace(os.path.normpath(root),"")
  else:
    return os.path.normpath(path)
class glue_build_py(build_py.build_py):
  def run(self):
    # create the git_version module
    if determine_git_version.in_git_repository():
      try:
        log.info("generating glue/git_version.py")
        git_version_fileobj = open("glue/git_version.py", "w")
        determine_git_version.write_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()
    elif os.path.exists("glue/git_version.py"):
      # We're probably being built from a release tarball; don't overwrite
      log.info("not in git checkout; using existing glue/git_version.py")
    else:
      log.info("not in git checkout; writing empty glue/git_version.py")
      try:
        git_version_fileobj = open("glue/git_version.py", "w")
        determine_git_version.write_empty_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()

    # resume normal build procedure
    build_py.build_py.run(self)

class glue_install(install.install):
  def run(self):

    # create the user env scripts
    if self.install_purelib == self.install_platlib:
      glue_pythonpath = self.install_purelib
    else:
      glue_pythonpath = self.install_platlib + ":" + self.install_purelib

    glue_prefix = remove_root(self.prefix,self.root)
    glue_install_scripts = remove_root(self.install_scripts,self.root)
    glue_pythonpath = remove_root(glue_pythonpath,self.root)
    glue_install_platlib = remove_root(self.install_platlib,self.root)
    
    log.info("creating glue-user-env.sh script")
    env_file = open(os.path.join('etc','glue-user-env.sh'),'w')
    print >> env_file, "# Source this file to access GLUE"
    print >> env_file, "GLUE_PREFIX=" + glue_prefix
    print >> env_file, "export GLUE_PREFIX"
    print >> env_file, "PATH=" + glue_install_scripts + ":${PATH}"
    print >> env_file, "PYTHONPATH=" + glue_pythonpath + ":${PYTHONPATH}"
    print >> env_file, "LD_LIBRARY_PATH=" + glue_install_platlib + ":${LD_LIBRARY_PATH}"
    print >> env_file, "DYLD_LIBRARY_PATH=" + glue_install_platlib + ":${DYLD_LIBRARY_PATH}"
    print >> env_file, "export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH"
    env_file.close()

    log.info("creating glue-user-env.csh script")
    env_file = open(os.path.join('etc','glue-user-env.csh'),'w')
    print >> env_file, "# Source this file to access GLUE"
    print >> env_file, "setenv GLUE_PREFIX " + glue_prefix
    print >> env_file, "setenv PATH " + glue_install_scripts + ":${PATH}"
    print >> env_file, "if ( $?PYTHONPATH ) then"
    print >> env_file, "  setenv PYTHONPATH " + glue_pythonpath + ":${PYTHONPATH}"
    print >> env_file, "else"
    print >> env_file, "  setenv PYTHONPATH " + glue_pythonpath
    print >> env_file, "endif"
    print >> env_file, "if ( $?LD_LIBRARY_PATH ) then"
    print >> env_file, "  setenv LD_LIBRARY_PATH " + glue_install_platlib + ":${LD_LIBRARY_PATH}"
    print >> env_file, "else"
    print >> env_file, "  setenv LD_LIBRARY_PATH " + glue_install_platlib
    print >> env_file, "endif"
    print >> env_file, "if ( $?DYLD_LIBRARY_PATH ) then"
    print >> env_file, "  setenv DYLD_LIBRARY_PATH " + glue_install_platlib + ":${DYLD_LIBRARY_PATH}"
    print >> env_file, "else"
    print >> env_file, "  setenv DYLD_LIBRARY_PATH " + glue_install_platlib
    print >> env_file, "endif"
    env_file.close()

    # now run the installer
    install.install.run(self)

class glue_clean(clean.clean):
  def finalize_options (self):
    clean.clean.finalize_options(self)
    self.clean_files = [ 'misc/__init__.pyc', 'misc/determine_git_version.pyc' ]

  def run(self):
    clean.clean.run(self)
    for f in self.clean_files:
      self.announce('removing ' + f)
      try:
        os.unlink(f)
      except:
        log.warn("'%s' does not exist -- can't clean it" % f)

class glue_sdist(sdist.sdist):
  def run(self):
    # remove the automatically generated user env scripts
    for script in [ 'glue-user-env.sh', 'glue-user-env.csh' ]:
      log.info( 'removing ' + script )
      try:
        os.unlink(os.path.join('etc',script))
      except:
        pass

    # create the git_version module
    if determine_git_version.in_git_repository():
      log.info("generating glue/git_version.py")
      try:
        git_version_fileobj = open("glue/git_version.py", "w")
        determine_git_version.write_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()
    else:
      log.info("not in git checkout; writing empty glue/git_version.py")
      try:
        git_version_fileobj = open("glue/git_version.py", "w")
        determine_git_version.write_empty_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()

    # now run sdist
    sdist.sdist.run(self)

setup(
  name = "glue",
  version = ver,
  author = "Duncan Brown",
  author_email = "dbrown@ligo.caltech.edu",
  description = "Grid LSC User Engine",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/",
  license = 'See file LICENSE',
  packages = [ 'glue', 'glue.gracedb', 'glue.lars', 'glue.lars.cli', 'glue.lars.util', 'glue.ligolw', 'glue.ligolw.utils', 'glue.lvalert', 'glue.segmentdb' ],
  cmdclass = {
    'build_py' : glue_build_py,
    'install' : glue_install,
    'clean' : glue_clean,
    'sdist' : glue_sdist
  },
  ext_modules = [
    Extension(
      "glue.ligolw.tokenizer",
      [
        "glue/ligolw/tokenizer.c",
        "glue/ligolw/tokenizer.Tokenizer.c",
        "glue/ligolw/tokenizer.RowBuilder.c",
        "glue/ligolw/tokenizer.RowDumper.c"
      ],
      include_dirs = [ "glue/ligolw" ]
    ),
    Extension(
      "glue.ligolw.__ilwd",
      [
        "glue/ligolw/ilwd.c"
      ],
      include_dirs = [ "glue/ligolw" ]
    ),
    Extension(
      "glue.__segments",
      [
        "src/segments/segments.c",
        "src/segments/infinity.c",
        "src/segments/segment.c",
        "src/segments/segmentlist.c"
      ],
      include_dirs = [ "src/segments" ]
    )
  ],
  scripts = [
    os.path.join('bin','gracedb'),
    os.path.join('bin','LSCdataFind'),
    os.path.join('bin','LSCdataFindcheck'),
    os.path.join('bin','ligo_data_find'),
    os.path.join('bin','lars'),
    os.path.join('bin','lars_add'),
    os.path.join('bin','lars_search'),
    os.path.join('bin','ldbdc'),
    os.path.join('bin','ldg_submit_dax'),
    os.path.join('bin','dmtdq_seg_insert'),
    os.path.join('bin','ligolw_add'),
    os.path.join('bin','ligolw_cut'),
    os.path.join('bin','ligolw_burst2mon'),
    os.path.join('bin','ligolw_inspiral2mon'),
    os.path.join('bin','ligolw_print'),
    os.path.join('bin','ligolw_sqlite'),
    os.path.join('bin','ligolw_segments_from_cats'),
    os.path.join('bin','ligolw_cbc_glitch_page'),
    os.path.join('bin','ligolw_segment_insert'),
    os.path.join('bin','ligolw_segment_intersect'),
    os.path.join('bin','ligolw_segment_diff'),
    os.path.join('bin','ligolw_segment_union'),
    os.path.join('bin','ligolw_segment_query'),
    os.path.join('bin','ligolw_veto_sngl_trigger'),
    os.path.join('bin','ligolw_dq_query'),
    os.path.join('bin','ligolw_dq_active'),
    os.path.join('bin','ligolw_dq_active_cats'),
    os.path.join('bin','lvalert_admin'),
    os.path.join('bin','lvalert_send'),
    os.path.join('bin','lvalert_listen'),
       os.path.join('bin','ldbdd'),
    os.path.join('bin','ligolw_publish_dqxml'),
    os.path.join('bin','segdb_coalesce'),
    os.path.join('bin', 'ligolw_print_tables') ],
  data_files = [
    ( 'etc',
      [ os.path.join('etc','ldg-sites.xml'),
        os.path.join('etc','pegasus-properties.bundle'),
        os.path.join('etc','glue-user-env.sh'),
        os.path.join('etc','glue-user-env.csh'),
        os.path.join('etc','ldbdserver.ini'),
        os.path.join('etc','ldbduser.ini'),
        os.path.join('etc','ligolw.xsl'),
        os.path.join('etc','ligolw.js'),
        os.path.join('etc','LDBDWServer.wsgi'),
        os.path.join('etc','ligolw_dtd.txt') ]
    ),
    ( os.path.join( 'etc', 'httpd', 'conf.d' ),
      [
        os.path.join('etc', 'segdb.conf')
      ]
    ),
    ( os.path.join( 'var', 'php', 'seginsert' ),
      [
        os.path.join('src', 'php', 'seginsert','index.php'),
        os.path.join('src', 'php', 'seginsert','flagcheck.php'),
        os.path.join('src', 'php', 'seginsert','ligolw.xsl'),
        os.path.join('src', 'php', 'seginsert','listflags.php'),
        os.path.join('src', 'php', 'seginsert','submitflag.php')
      ]
    ),
    ( os.path.join( 'var', 'php', 'seginsert', 'img' ),
      [
        os.path.join('src', 'php', 'seginsert','img','LIGOLogo.gif'),
        os.path.join('src', 'php', 'seginsert','img','brace.gif'),
        os.path.join('src', 'php', 'seginsert','img','lsc.gif'),
        os.path.join('src', 'php', 'seginsert','img','plus.gif')
      ]
    ),
    ( os.path.join( 'var', 'php', 'seginsert', 'scripts' ),
      [
        os.path.join('src', 'php', 'seginsert','scripts','footer.php'),
        os.path.join('src', 'php', 'seginsert','scripts','form_day_list.php'),
        os.path.join('src', 'php', 'seginsert','scripts','form_month_list.php'),
        os.path.join('src', 'php', 'seginsert','scripts','form_year_list.php'),
        os.path.join('src', 'php', 'seginsert','scripts','header.php'),
        os.path.join('src', 'php', 'seginsert','scripts','style.css'),
        os.path.join('src', 'php', 'seginsert','scripts','styletitle.php'),
        os.path.join('src', 'php', 'seginsert','scripts','time_conv_functions.php')
      ]
    ),
    ( os.path.join( 'var', 'php', 'dq_report' ),
      [
        os.path.join('src', 'php', 'dq_report','index.php'),
        os.path.join('src', 'php', 'dq_report','get_report.php'),
        os.path.join('src', 'php', 'dq_report','header.php')
      ]
    )
  ]
)
