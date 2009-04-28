# $Id$
#
# Copyright (C) 2006  Duncan A. Brown
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

import sys

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from pylal import llwapp
from pylal.date import LIGOTimeGPS
from pylal import SnglInspiralUtils

__author__ = "Duncan Brown <dbrown@ligo.caltech.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Preparation
#
# =============================================================================
#

def get_tables(doc):
  snglinspiraltable = table.get_table(
    doc, lsctables.SnglInspiralTable.tableName)

  input_times = None
  output_times = None
  try:
    searchsummtable = table.get_table(
      doc, lsctables.SearchSummaryTable.tableName)
    input_times = searchsummtable.get_inlist().extent()
    output_times = searchsummtable.get_outlist().extent()
  except ValueError:
    pass
    
  return input_times, output_times, snglinspiraltable


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

def append_process(doc, **kwargs):
  process = llwapp.append_process(
    doc, program = "ligolw_sicluster", version = __version__, 
    cvs_repository = "lscsoft", cvs_entry_time = __date__, 
    comment = kwargs["comment"])

  ligolw_process.append_process_params(doc, process, 
    [("--cluster-window", "lstring", kwargs["cluster_window"])])
  if kwargs["snr_threshold"] > 0:
    ligolw_process.append_process_params(doc, process, 
      [("--snr-threshold", "lstring", kwargs["snr_threshold"])])
  if kwargs["sort_descending_snr"]:
    ligolw_process.append_process_params(doc, process, 
      [("--sort-descending-snr", "lstring", " ")])
  if kwargs["sort_ascending_snr"]:
    ligolw_process.append_process_params(doc, process, 
      [("--sort-ascending-snr", "lstring", " ")])

  return process


#
# =============================================================================
#
#                             Clustering Algorithm
#
# =============================================================================
#

def SnglInspiralCluster(a, b):
  """
  Replace a with a cluster constructed from a and b. 
  """
  if b.snr >= a.snr:
    return b
  else:
    return a


def ClusterSnglInspiralTable(triggers, testfunc, clusterfunc, 
  twindow, bailoutfunc = None):
  """
  Cluster the triggers in the list.  testfunc should accept a pair of
  triggers, and return 0 if they should be clustered.  clusterfunc
  should accept a pair of triggers, and replace the contents of the
  first with a cluster constructed from the two.  If bailoutfunc is
  provided, the triggers will be sorted using testfunc as a
  comparison operator, and then only pairs of triggers for which
  bailoutfunc returns 0 will be considered for clustering.
  """
  while True:
    did_cluster = False

    if bailoutfunc:
      triggers.sort(testfunc)

    i = 0
    while i < len(triggers):
      j = i + 1
      while j < len(triggers):
        if not testfunc(triggers[i], triggers[j],twindow):
          triggers[i] = clusterfunc(triggers[i], triggers[j])
          del triggers[j]
          did_cluster = True
        else:
          if bailoutfunc:
            if bailoutfunc(triggers[i], triggers[j]):
              break
          j += 1
      i += 1

    if not did_cluster:
      return


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_sicluster(doc, **kwargs):
  # Extract segments and tables
  inseg, outseg, snglinspiraltable = get_tables(doc)

  # Add process information
  try:
    process = append_process(doc, **kwargs)
  except ValueError:
    process = None

  # Delete all triggers below threshold
  if kwargs["snr_threshold"] > 0:
    thresh = float(kwargs["snr_threshold"])
    if kwargs["verbose"]:
      print >>sys.stderr, "discarding triggers with snr < %f..." % \
        kwargs["snr_threshold"]
    for i in range(len(snglinspiraltable) - 1, -1, -1):
      if snglinspiraltable[i].snr <= thresh:
        del snglinspiraltable[i]

  # Cluster
  if kwargs["verbose"]:
    print >>sys.stderr, "clustering..."
  ClusterSnglInspiralTable(snglinspiraltable, 
    kwargs["testfunc"], kwargs["clusterfunc"],
    LIGOTimeGPS(kwargs["cluster_window"]), kwargs["bailoutfunc"])

  # Sort by signal-to-noise ratio
  if kwargs["sort_ascending_snr"] or kwargs["sort_descending_snr"]:
    if kwargs["verbose"]:
      print >>sys.stderr, "sorting by snr..."
    snglinspiraltable.sort(SnglInspiralUtils.CompareSnglInspiralBySnr)
    if kwargs["sort_descending_snr"]:
      snglinspiraltable.reverse()

  # Add search summary information
  if process and inseg and outseg:
    llwapp.append_search_summary(doc, process, inseg = inseg, outseg = outseg, 
      nevents = len(snglinspiraltable))
  if process:
    llwapp.set_process_end_time(process)

  return doc
