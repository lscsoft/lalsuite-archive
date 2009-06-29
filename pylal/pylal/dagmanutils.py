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

try:
  import sqlite3
except ImportError:
  from pysqlite2 import dbapi2 as sqlite3

import datetime
import hashlib
import os
import sys
import time

itertools = __import__("itertools")

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"[11:-2]

#
# Enums
#
STATE_NOT_YET_RUN = 0
STATE_IDLE = 1
STATE_RUNNING = 2
STATE_HELD = 3
STATE_SUCCEEDED = 4
STATE_FAILED = 5
STATE_TERMINATED = 6

#
# DAGManLog class, to represent the state of a DAGMan run
#

class DAGManLog(object):
  # specify state as (state_id, name)
  states = ((STATE_NOT_YET_RUN, "not yet run"),
            (STATE_IDLE, "idle"),
            (STATE_RUNNING, "running"),
            (STATE_HELD, "held"),
            (STATE_SUCCEEDED, "succeeded"),
            (STATE_FAILED, "failed"),
            (STATE_TERMINATED, "terminated (unknown status)"))

  # state transition map
  # (grep_key, state_id, node_id_col)
  # more verbosely: (text that, if found, triggers the state transition,
  #                  state_id of destination state,
  #                  index to space-delimited column that contains node name)
  state_map = (("ULOG_SUBMIT", STATE_IDLE, 7),
               ("ULOG_EXECUTE", STATE_RUNNING, 7),
               ("ULOG_JOB_HELD", STATE_HELD, 7),
               ("ULOG_JOB_RELEASED", STATE_IDLE, 7),
               ("ULOG_JOB_EVICTED", STATE_IDLE, 7),
               ("ULOG_JOB_TERMINATED", STATE_TERMINATED, 7),
               ("ULOG_POST_TERMINATED", STATE_TERMINATED, 7))

  def __init__(self, db_filename=":memory:"):
    # create database
    self.conn = sqlite3.connect(db_filename,
        detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = self.conn.cursor()

    # create tables
    c.execute("""CREATE TABLE IF NOT EXISTS states
    (state_id INTEGER PRIMARY KEY, name TEXT UNIQUE)""")
    c.execute("""CREATE TABLE IF NOT EXISTS jobs
    (job_id INTEGER PRIMARY KEY, sub_file TEXT UNIQUE)""")
    c.execute("""CREATE TABLE IF NOT EXISTS nodes
    (node_id TEXT UNIQUE, job_id INTEGER)""")
    c.execute("""CREATE TABLE IF NOT EXISTS transitions
    (ts TIMESTAMP, node_id TEXT, state_id INTEGER)""")
    c.execute("""CREATE TABLE IF NOT EXISTS hash (hash TEXT)""")

    # populate states
    c.executemany("INSERT OR IGNORE INTO states (state_id, name) VALUES (?, ?)", self.states)
    self.conn.commit()

  def read_dag(self, dag_fileobj):
    """
    Populate jobs and create an initial transition to never run state.
    """
    c = self.conn.cursor()
    for line in dag_fileobj:
      ts = datetime.datetime(1900, 1, 1, 0, 0, 0, 0)
      if line.startswith("JOB"):
        tup = line.split()
        node_id = tup[1]
        sub_file = tup[2]
        c.execute("INSERT OR IGNORE INTO jobs (sub_file) VALUES (?)", (sub_file,))
        job_id = c.execute("SELECT job_id FROM jobs WHERE sub_file = ?", (sub_file,)).next()[0]
        c.execute("INSERT OR REPLACE INTO nodes (node_id, job_id) VALUES (?, ?)", (node_id, job_id))
        c.execute("INSERT INTO transitions (ts, node_id, state_id) VALUES (?, ?, ?)", (ts, node_id, 0))
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

    # try to infer year from dagman.out, else use current year.
    if hasattr(dagman_out_fileobj, "name"):
        start = datetime.datetime.fromtimestamp(\
            os.path.getctime(dagman_out_fileobj.name))
    else:
        start = datetime.datetime.now()
    year = str(start.year)

    dagman_out_hash = hashlib.md5()
    dagman_out_iter = iter(dagman_out_fileobj)
    for line in dagman_out_iter:
      dagman_out_hash.update(line)
      # determine what state, if any, this line denotes, and record it
      for grep_key, state_id, node_id_col in self.state_map:
        if grep_key in line: # we have a match
          tup = line.split()
          node_id = tup[node_id_col]
          ts = datetime.datetime(*time.strptime(\
            year + "/" + tup[0] + " " + tup[1],
            "%Y/%m/%d %H:%M:%S")[:6])
          # if ts is before log file creation, must increment year
          if ts < start:
            year = str(start.year + 1)
            ts = datetime.datetime(*time.strptime(year + "/" + tup[0] + " "\
                + tup[1], "%Y/%m/%d %H:%M:%S")[:6])
          if state_id == STATE_TERMINATED:
              line = dagman_out_iter.next() # job status on next line
              dagman_out_hash.update(line)
              if "completed successfully" in line:
                state_id = STATE_SUCCEEDED
              elif "failed" in line:
                state_id = STATE_FAILED
              else:
                sys.stderr.write("Warning: Cannot determine success or "\
                                 "failure: " + line + "\n")

          c.execute("INSERT INTO transitions (ts, node_id, state_id) VALUES (?, ?, ?)", (ts, node_id, state_id))
          break

    # write hash of the file we just parsed
    c.execute("INSERT INTO hash (hash) VALUES (?)", (dagman_out_hash.hexdigest(),))
    self.conn.commit()

  def get_state_subfile_dict(self):
    """
    Return a set of nested dictionaries of the form:
    {state: {sub_file: number in state}}.
    """
    c = self.conn.cursor()
    d = {}
    for state_id, state_name in self.states:
      d[state_name] = dict(c.execute("SELECT sub_file, count(*) FROM (SELECT * FROM transitions GROUP BY transitions.node_id) LEFT NATURAL JOIN nodes NATURAL JOIN jobs WHERE state_id=? GROUP BY sub_file", (state_id,)))
    return d

  def get_state_dict(self):
    """
    Return a dictionary of {node state: number in state}.
    """
    c = self.conn.cursor()
    return dict(c.execute("SELECT name, count(*) FROM (SELECT * FROM transitions GROUP BY transitions.node_id) NATURAL JOIN states GROUP BY states.state_id"))

  def get_runtime_dict(self):
    """
    Return a dictionary of {node_id: total run time} for nodes that completed
    successfully.  The times included are those in STATE_RUNNING.
    """
    c = self.conn.cursor()

    # only allow durations from nodes that have completed
    durations = dict((key, 0) for (key,) in \
      c.execute("SELECT node_id FROM transitions WHERE state_id=?",
                (STATE_SUCCEEDED,)))

    for node_id in durations:
      # find starts and stops
      starts = c.execute("SELECT ts FROM transitions WHERE node_id=? AND state_id=? ORDER BY ts", (node_id, STATE_RUNNING))
      stops = c.execute("SELECT ts FROM transitions WHERE node_id=? AND state_id NOT IN (?, ?) ORDER BY ts LIMIT -1 OFFSET 1", (node_id, STATE_RUNNING, STATE_NOT_YET_RUN))

      # do the work; NB: Condor only has second precision at the moment.
      for (start,), (stop,) in itertools.izip(starts, stops):
        if stop < start:
          raise ValueError, "stop before start"
        dur = stop - start
        durations[node_id] += 86400 * dur.days + dur.seconds

    # now add together all nodes that belong to a given job
    # XXX: This can perhaps be put in a separate function.
    # job_durations = dict((job_name, [0, datetime.timedelta(0)]) \
    #   for job_name in c.execute("SELECT name FROM jobs"))
    # for node_id, dur in durations.iteritems():
    #   sub_file = c.execute("SELECT sub_file FROM (SELECT * FROM nodes WHERE node_id=?) LEFT NATURAL JOIN jobs", node_id)
    #   job_durations[sub_file][0] += 1
    #   job_durations[sub_file][1] += dur

    return durations

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

