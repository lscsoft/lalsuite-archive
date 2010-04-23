#!/usr/bin/env python
#
# git_version.py - determine git version info
#
# Copyright (C) 2009,2010 Adam Mercer <adam.mercer@ligo.org>
#                         Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#
# preamble
#

# metadata
__author__ = 'Adam Mercer <adam.mercer@ligo.org>'

# import required system modules
import exceptions
import subprocess
import time
import os

#
# process management functions
#

class GitInvocationError(exceptions.LookupError):
  pass

# return output from running given command
def call_out(command):
  """
  Run the given command (with shell=False) and return a tuple of
  (int returncode, str output). Strip the output of enclosing whitespace.
  """
  # start external command process
  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # get outputs
  out, _ = p.communicate()

  return p.returncode, out.strip()

def check_call_out(command):
  """
  Run the given command (with shell=False) and return the output as a
  string. Strip the output of enclosing whitespace.
  If the return code is non-zero, throw GitInvocationError.
  """
  # start external command process
  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # get outputs
  out, _ = p.communicate()

  # throw exception if process failed
  if p.returncode != 0:
    raise GitInvocationError, 'failed to run "%s"' % " ".join(command)

  return out.strip()


#
# primary functions
#
def in_git_repository():
  """
  Return True if we are in a git repository and False if not.

  NB: Unfortunately there is a magic number without any documentation to back
  it up. It turns out that git status returns non-zero exit codes for all sorts
  of success conditions, but I cannot find any documentation of them. 128 was
  determined empirically. I sure hope that it's portable.
  """
  ret_code, git_path = call_out(('/usr/bin/which', 'git'))
  return (ret_code != 1) and not git_path.startswith('no') \
      and (call_out((git_path, 'status'))[0] != 128)

def write_git_version(fileobj):
  """
  Query git to determine current repository status and write a Python module
  with this information.

  Ex:
  >>> write_git_version(open("git_version.py", "w"))
  >>> import git_version
  >>> print git_version.id
  1b0549019e992d0e001f3c28e8488946f825e873
  """
  git_path = check_call_out(('/usr/bin/which', 'git'))

  # determine current time and treat it as the build time
  build_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())

  # determine builder
  git_builder_name = check_call_out((git_path, 'config', 'user.name'))
  git_builder_email = check_call_out((git_path, 'config', 'user.email'))
  git_builder = "%s <%s>" % (git_builder_name, git_builder_email)

  # determine git id
  git_id = check_call_out((git_path, 'log', '-1', '--pretty=format:%H'))

  # determine commit date, iso utc
  git_udate = float(check_call_out((git_path, 'log', '-1', '--pretty=format:%ct')))
  git_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime(git_udate))

  # determine branch
  branch_match = check_call_out((git_path, 'rev-parse', '--symbolic-full-name',
    'HEAD'))
  if branch_match == "HEAD":
    git_branch = None
  else:
    git_branch = os.path.basename(branch_match)

  # determine tag
  status, git_tag = call_out((git_path, 'describe', '--exact-match',
    '--tags', git_id))
  if status != 0:
    git_tag = None

  # determine author and committer
  git_author_name = check_call_out((git_path, 'log', '-1', '--pretty=format:%an'))
  git_author_email = check_call_out((git_path, 'log', '-1', '--pretty=format:%ae'))
  git_author = '%s <%s>' % (git_author_name, git_author_email)
  git_committer_name = check_call_out((git_path, 'log', '-1', '--pretty=format:%cn'))
  git_committer_email = check_call_out((git_path, 'log', '-1', '--pretty=format:%ce'))
  git_committer = '%s <%s>' % (git_committer_name, git_committer_email)

  # refresh index
  check_call_out((git_path, 'update-index', '-q', '--refresh'))

  # check working copy for changes
  status_output = subprocess.call((git_path, 'diff-files', '--quiet'))
  if status_output != 0:
    git_status = 'UNCLEAN: Modified working tree'
  else:
    # check index for changes
    status_output = subprocess.call((git_path, 'diff-index', '--cached',
      '--quiet', 'HEAD'))
    if status_output != 0:
      git_status = 'UNCLEAN: Modified index'
    else:
      git_status = 'CLEAN: All modifications committed'

  # print details in a directly importable form
  print >>fileobj, 'id = "%s"' % git_id
  print >>fileobj, 'date = "%s"' % git_date
  print >>fileobj, 'branch = "%s"' % git_branch
  if git_tag is None:
      print >>fileobj, 'tag = None'
  else:
      print >>fileobj, 'tag = "%s"' % git_tag
  print >>fileobj, 'author = "%s"' % git_author
  print >>fileobj, 'author_name = "%s"' % git_author_name
  print >>fileobj, 'author_email = "%s"' % git_author_email
  print >>fileobj, 'builder = "%s"' % git_builder
  print >>fileobj, 'builder_name = "%s"' % git_builder_name
  print >>fileobj, 'builder_email = "%s"' % git_builder_email
  print >>fileobj, 'committer = "%s"' % git_committer
  print >>fileobj, 'committer_name = "%s"' % git_committer_name
  print >>fileobj, 'committer_email = "%s"' % git_committer_email
  print >>fileobj, 'status = "%s"' % git_status
  print >>fileobj, 'version = id'

  # add a verbose report for convenience
  print >>fileobj, 'verbose_msg = """%s"""' % \
  """Branch: %s
Tag: %s
Id: %s

Builder: %s
Build date: %s
Repository status: %s""" \
  % (git_branch, git_tag, git_id, git_builder, build_date, git_status)

def write_empty_git_version(fileobj):
  """
  A fallback function that can populate a git_version.py with null values.
  See write_git_version() above.
  """
  print >>fileobj, "id = None"
  print >>fileobj, "date = None"
  print >>fileobj, "branch = None"
  print >>fileobj, "tag = None"
  print >>fileobj, "author = None"
  print >>fileobj, "author_name = None"
  print >>fileobj, "author_email = None"
  print >>fileobj, "builder = None"
  print >>fileobj, "builder_name = None"
  print >>fileobj, "builder_email = None"
  print >>fileobj, "committer = None"
  print >>fileobj, "committer_name = None"
  print >>fileobj, "committer_email = None"
  print >>fileobj, "status = None"
  print >>fileobj, "version = None"
  print >>fileobj, "verbose_msg = \"No version information available; not built from a git repository\""

if __name__=="__main__":
  import sys
  write_git_version(sys.stdout)
