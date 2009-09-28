#!/usr/bin/env python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
#
"""
This checklist script is responsible for creating a MoinMoin Wiki
page.  This checklist page is meant to replace the original HTML based
page created by the script makeCheckList.py.  The output of this
script is a text file that can be cut and paste into a MoinMoin Wiki
editor, resulting in a Wikified checklist under MoinMoin version
control.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__date__    = '$Date$'
__version__ = '$Revision$'
__prog__    = 'makeCheckListWiki.py'


####################################################################
# Read in the standard python modules and the LSC python modules
# needed to create the MoinMoin file
####################################################################
import copy
import numpy
import optparse
import os
import random
import socket
import sys
import time
import urllib
sys.path.append('@PYTHONLIBDIR@')


####################################################################
# Custom wiki class to make writing MoinMoin text simpler
####################################################################
class wiki(object):
  def __init__(self,open_box=False,fname="wiki.txt"):
    if open_box: fname = "open_box_" + fname
    self.fname = fname
    self.file = open(fname,"w")

  def image_link(self,path,webserver):
    thumb = "thumb_" + path
    command = 'convert ' + path + ' -resize 300x300 -antialias ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    s = '[[ImageLink('+webserver+'/'+thumb+','+webserver+'/'+path+',width=300][,alt=none])]]'
    self.file.write(s)

  def image_table(self,image_list, webserver):
    if not image_list: return
    for j, i in enumerate(image_list):
      if not (j) % 3: self.file.write("\n\n||")
      self.image_link(i, webserver)
      self.file.write("||")
    self.file.write("\n\n")
  
  def image_glob(self, pat):
    image_list = []
    for image in glob.glob(pat):
      if 'thumb' in image: continue
      else: image_list.append(image)
    image_list.sort()
    return image_list

  def section(self,title):
    s = "=== "+title.strip()+" ===\n"
    self.file.write(s)

  def write(self,val):
    self.file.write(val)

  def finish(self):
    self.file.close()  

def parse_command_line():
  parser = OptionParser(version = "%prog CVS $Id$", usage = "%prog [options] [file ...]", description = "%prog computes mass/mass upperlimit")
  parser.add_option("--webserver", help = "Set the webserver path.  Required.  Example https://ldas-jobs.ligo.caltech.edu/~channa/highmass_months_23-24_summary_page")
  parser.add_option("--open-box", action = "store_true", help = "Produce open box page")
  parser.add_option("--output-name-tag", default = "", metavar = "name", help = "Set the basename for image search")
  opts, filenames = parser.parse_args()

  if not opts.webserver:
    print >>sys.stderr, "must specify a webserver"
    sys.exit(1)
  return opts, filenames
####################################################################
# End Custom wiki class
####################################################################



####################################################################
# Main part of script
####################################################################

####################################################################
# The command line arguements used to construct this page
# the options drive the creation of the page
####################################################################
usage = """uage: %prog [options]"""

parser = OptionParser(usage)
#Add all options to setup the checklist items
