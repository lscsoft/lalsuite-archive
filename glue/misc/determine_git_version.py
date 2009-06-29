#!/usr/bin/env python
#
# git_version.py - determine git version info
#
# Copyright (C) 2009, Adam Mercer <adam.mercer@ligo.org>,
#                     Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>
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
import re
import subprocess
import time

#
# process management functions
#

class GitInvocationError(exceptions.LookupError):
  pass

# return output from running given command
def run_external_command(command, honour_ret_code=True):
  # start external command process
  p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, close_fds=True)

  # get outputs
  try:
    output = p.stdout.read()
    err = p.stderr.read()
    ret_code = p.wait()
  finally:
    p.stdout.close()
    p.stderr.close()

  # throw exception if process failed
  if honour_ret_code and (ret_code != 0):
    raise GitInvocationError, 'failed to run "%s"' % command

  return ret_code, output

#
# primary functions
#
def in_git_repository():
  """
  Return True if we are in a git repository and False if not.

  NB: Unfortunately there is a magic number without any documentation to back
  it up. It turns out that git status returns non-zero exit codes for all sorts
  of success conditions, but I cannot find any documentation of them. 128 was
  determined empirically.  I sure hope that it's portable.
  The return code of 127 indicates that the git command cannot be found.
  """
  ret_code = run_external_command('git status', honour_ret_code=False)[0]
  return ret_code not in (127, 128)

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
  # determine current time and treat it as the build time
  build_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())

  # determine builder
  builder_name_cmd = 'git config user.name'
  builder_email_cmd = 'git config user.email'
  git_builder_name = run_external_command(builder_name_cmd,
    honour_ret_code=False)[1].strip()
  git_builder_email = run_external_command(builder_email_cmd,
    honour_ret_code=False)[1].strip()
  git_builder = "%s <%s>" % (git_builder_name, git_builder_email)

  # determine git id
  id_cmd = 'git log -1 --pretty="format:%H"'
  git_id = run_external_command(id_cmd)[1].strip()

  # determine commit date, iso utc
  date_cmd = 'git log -1 --pretty="format:%ct"'
  git_udate = float(run_external_command(date_cmd)[1].strip())
  git_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime(git_udate))

  # determine branch
  branch_cmd = 'git branch --no-color'
  branch_regexp = re.compile(r"\* ((?!\(no branch\)).*)", re.MULTILINE)
  branch_match = branch_regexp.search(run_external_command(branch_cmd)[1].strip())
  if branch_match is None:
    git_branch = None
  else:
    git_branch = branch_match.group(1)

  # determine tag
  tag_cmd = 'git describe --exact-match --tags %s' % git_id
  status, git_tag = run_external_command(tag_cmd, honour_ret_code=False)
  git_tag = git_tag.strip()
  if status != 0:
    git_tag = None

  # determine author and committer
  author_name_cmd = 'git log -1 --pretty="format:%an"'
  author_email_cmd = 'git log -1 --pretty="format:%ae"'
  committer_name_cmd = 'git log -1 --pretty="format:%cn"'
  committer_email_cmd = 'git log -1 --pretty="format:%ce"'
  git_author_name = run_external_command(author_name_cmd)[1].strip()
  git_author_email = run_external_command(author_email_cmd)[1].strip()
  git_author = '%s <%s>' % (git_author_name, git_author_email)
  git_committer_name = run_external_command(committer_name_cmd)[1].strip()
  git_committer_email = run_external_command(committer_email_cmd)[1].strip()
  git_committer = '%s <%s>' % (git_committer_name, git_committer_email)

  # determine tree status
  status_cmd = 'git status'
  status_output = run_external_command(status_cmd, honour_ret_code=False)[1]\
                  .strip()
  if ("# Changed but not updated:\n" in status_output) or \
     ("# Changes to be committed:\n" in status_output):
    git_status = 'UNCLEAN: Some modifications not committed'
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
