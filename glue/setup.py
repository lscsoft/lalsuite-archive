# $Id$
# 
# setup script for glue

import os, sys
import subprocess
import time

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

from distutils.core import setup, Extension
from distutils.command import build_py
from distutils.command import install
from distutils.command import sdist
from distutils.command import clean
from distutils import log

from misc import generate_vcs_info as gvcsi

ver = "1.33"

def remove_root(path,root):
  if root:
    return os.path.normpath(path).replace(os.path.normpath(root),"")
  else:
    return os.path.normpath(path)

def write_build_info():
  """
  Get VCS info from misc/generate_vcs_info.py and add build information.
  Substitute these into misc/git_version.py.in to produce glue/git_version.py.
  """
  vcs_info = gvcsi.generate_git_version_info()

  # determine current time and treat it as the build time
  build_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())

  # determine builder
  retcode, builder_name = gvcsi.call_out(('git', 'config', 'user.name'))
  if retcode:
    builder_name = "Unknown User"
  retcode, builder_email = gvcsi.call_out(('git', 'config', 'user.email'))
  if retcode:
    builder_email = ""
  builder = "%s <%s>" % (builder_name, builder_email)

  sed_cmd = ('sed',
             '-e', 's/@ID@/%s/' % vcs_info.id,
             '-e', 's/@DATE@/%s/' % vcs_info.date,
             '-e', 's/@BRANCH@/%s/' % vcs_info.branch,
             '-e', 's/@TAG@/%s/' % vcs_info.tag,
             '-e', 's/@AUTHOR@/%s/' % vcs_info.author,
             '-e', 's/@COMMITTER@/%s/' % vcs_info.committer,
             '-e', 's/@STATUS@/%s/' % vcs_info.status,
             '-e', 's/@BUILDER@/%s/' % builder,
             '-e', 's/@BUILD_DATE@/%s/' % build_date,
             'misc/git_version.py.in')

  # FIXME: subprocess.check_call becomes available in Python 2.5
  sed_retcode = subprocess.call(sed_cmd,
    stdout=open('glue/git_version.py', 'w'))
  if sed_retcode:
    raise gvcsi.GitInvocationError

class glue_build_py(build_py.build_py):
  def run(self):
    # create the git_version module
    log.info("Generating glue/git_version.py")
    try:
      write_build_info()
    except gvcsi.GitInvocationError:
      if os.path.exists("glue/git_version.py"):
        # We're probably being built from a release tarball; don't overwrite
        log.info("Not in git checkout or cannot find git executable; "\
            "using existing glue/git_version.py")
      else:
        log.error("Not in git checkout or cannot find git executable "\
            "and no glue/git_version.py. Exiting.")
        sys.exit(1)

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
    self.clean_files = [ 'misc/__init__.pyc', 'misc/generate_vcs_info.pyc' ]

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
    log.info("Generating glue/git_version.py")
    try:
      write_build_info()
    except gvcsi.GitInvocationError:
      log.error("Not in git checkout or cannot find git executable and no "\
        "glue/git_version.py. Exiting.")
      sys.exit(1)

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
  packages = [ 'glue', 'glue.gracedb', 'glue.lars', 'glue.lars.cli', 'glue.lars.util', 'glue.ligolw', 'glue.ligolw.utils', 'glue.lvalert', 'glue.segmentdb', 'glue.nmi' ],
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
    os.path.join('bin','ligolw_segments_from_cats_split'),
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
    os.path.join('bin','ligolw_dq_grapher'),
    os.path.join('bin','lvalert_admin'),
    os.path.join('bin','lvalert_send'),
    os.path.join('bin','lvalert_listen'),
       os.path.join('bin','ldbdd'),
    os.path.join('bin','ligolw_publish_dqxml'),
    os.path.join('bin','ligolw_diff'),
    os.path.join('bin','ligolw_geo_fr_to_dq'),
    os.path.join('bin','segdb_coalesce'),
    os.path.join('bin', 'glue_nmi_genbuild'),
    os.path.join('bin', 'ligolw_print_tables'),
    os.path.join('bin', 'ligolw_veto_def_check')],
  data_files = [
    ( 'etc',
      [ os.path.join('etc','ldg-sites.xml'),
        os.path.join('etc','cbcwebpage.css'),
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
    ( os.path.join( 'share','nmi' ),
      [ 
        os.path.join('src', 'nmi', 'builds', 'lalsuite-build.spec'),
        os.path.join('src', 'nmi', 'builds', 'lalsuite-build-scripts.location'),
      ]
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
