# $Id$
#
# Copyright (C) 2008  Nickolas Fotopoulos
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

import sqlite3
import sys

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"[11:-2]

#
# DAGManLog class, to represent the state of a DAGMan run
#

class DAGManLog(object):  
  # specify state as (state_id, name)
  states = ((0, "not yet run"),
            (1, "idle"),
            (2, "running"),
            (3, "held"),
            (4, "succeeded"),
            (5, "failed"))

  # state transition map
  # (grep_key, state_id, node_id_col)
  # more verbosely: (text that, if found, triggers the state transition,
  #                  state_id of destination state,
  #                  index to space-delimited column that contains node name)
  state_map = (("ULOG_SUBMIT", 1, 7),
               ("ULOG_EXECUTE", 2, 7),
               ("ULOG_JOB_HELD", 3, 7),
               ("ULOG_JOB_RELEASED", 1, 7),
               ("completed successfully", 4, 3),
               ("failed with status", 5, 3))

  def __init__(self, db_filename=":memory:"):
    # create database
    self.conn = sqlite3.connect(db_filename)
    c = self.conn.cursor()

    # create tables
    c.execute("""CREATE TABLE IF NOT EXISTS states
    (state_id INTEGER PRIMARY KEY, name TEXT UNIQUE)""")
    # c.execute("""CREATE TABLE IF NOT EXISTS state_map
    # (grep_key TEXT, state_id INTEGER, node_id_col INTEGER)""")
    c.execute("""CREATE TABLE IF NOT EXISTS jobs
    (job_id INTEGER PRIMARY KEY, sub_file TEXT UNIQUE)""")
    c.execute("""CREATE TABLE IF NOT EXISTS nodes
    (node_id TEXT UNIQUE, job_id INTEGER, state_id INTEGER)""")

    # populate state and state_map
    c.executemany("INSERT OR IGNORE INTO states (state_id, name) VALUES (?, ?)", self.states)
    # c.executemany("INSERT OR IGNORE INTO state_map (grep_key, state_id, node_id_col) VALUES (?, ?, ?)", self.state_map)
    self.conn.commit()

  def read_dag(self, dag_fileobj):
    """
    Populate jobs and initialize nodes to never run state.
    """
    c = self.conn.cursor()
    for line in dag_fileobj:
      if line.startswith("JOB"):
        tup = line.split()
        node_id = tup[1]
        sub_file = tup[2]
        c.execute("INSERT OR IGNORE INTO jobs (sub_file) VALUES (?)", (sub_file,))
        job_id = c.execute("SELECT job_id FROM jobs WHERE sub_file = ?", (sub_file,)).next()[0]
        c.execute("INSERT OR REPLACE INTO nodes (node_id, job_id, state_id) VALUES (?, ?, ?)", (node_id, job_id, 0))
    self.conn.commit()

  def get_job_dict(self):
    """
    Return a dictionary of counts for each sub file.
    """
    c = self.conn.cursor()
    return dict(c.execute("SELECT sub_file, count(*) FROM jobs JOIN nodes ON nodes.job_id = jobs.job_id GROUP BY sub_file"))

  def read_dagman_out(self, dagman_out_fileobj):
    """
    Parse .dag.dagman.out file and update the status of each node.
    """
    c = self.conn.cursor()
    for line in dagman_out_fileobj:
      # determine what state, if any, this line denotes, and record it
      for grep_key, state_id, node_id_col in self.state_map:
        if grep_key in line:
          node_id = line.split()[node_id_col]
          c.execute("UPDATE nodes SET state_id = ? WHERE node_id = ?",
                    (state_id, node_id))
          break
    self.conn.commit()

  def get_state_subfile_dict(self):
    """
    Return a set of nested dictionaries of the form:
    {state: {sub_file: number in state}}.
    """
    c = self.conn.cursor()
    d = {}
    for state_id, state_name in self.states:
      d[state_name] = dict(c.execute("SELECT sub_file, count(*) FROM jobs JOIN nodes ON nodes.job_id = jobs.job_id WHERE nodes.state_id = ? GROUP BY sub_file", (state_id,)))
    return d
      
  def get_state_dict(self):
    """
    Return a dictionary of {node state: number in state}.
    """
    c = self.conn.cursor()
    return dict(c.execute("SELECT name, count(*) FROM nodes JOIN states WHERE nodes.state_id = states.state_id GROUP BY states.state_id"))


#
# Text summary functions
#

def print_dict_of_counts(d):
  """
  Print the entries of a dict of counts in some standard way:

      value[0] key[0]
      value[1] key[1]
      ...
  """
  for key, value in d.iteritems():
    print "%10d %s" % (value, key)

def print_nested_dict_of_counts(d):
  """
  Print the entries of a nested dict of counts in some standard way:
  d = {key1: {key2: value}}

  key1[0]:
      value[0] key2[0]
      value[1] key2[1]
  ...
  """
  for key1, sub_dict in d.iteritems():
    print key1 + ":"
    print_dict_of_counts(sub_dict)