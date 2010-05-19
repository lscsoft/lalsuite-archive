"""
A simple module for acquiring pidfile locks (e.g., for use by daemons).

Copyright (C) 2010 by Peter F. Couvares, Syracuse University
mailto: pfcouvar@syr.edu
"""

import os
import sys
import errno
import fcntl
import time

# from glue import git_version

# __author__ = 'Peter Couvares <pfcouvar@syr.edu>'
# __version__ = "git id %s" % git_version.id
#__date__ = git_version.date


# inspired by http://stackoverflow.com/questions/1005972
def pid_exists(pid):
    """ Returns true if the given pid exists, false otherwise. """
    try:
        # signal 0 is harmless and can be safely used to probe pid existence
        os.kill(pid, 0)
    except OSError, e:
        # "permission denied" proves existence; otherwise, no such pid
        return e.errno == errno.EPERM
    else:
        return True


def get_lock(lockfile):
    """
    Tries to write a lockfile containing the current pid.  Excepts if
    the lockfile already contains the pid of a running process.

    Although this should prevent a lock from being granted twice, it
    can theoretically deny a lock unjustly in the unlikely event that
    the original process is gone but another unrelated process has
    been assigned the same pid by the OS.
    """

    pidfile = open(lockfile, "a+")

    # here we do some meta-locking by getting an exclusive lock on the
    # pidfile before reading it, to prevent two daemons from seeing a
    # stale lock at the same time, and both trying to run
    try:
        fcntl.flock(pidfile.fileno(), fcntl.LOCK_EX|fcntl.LOCK_NB)
    except IOError,e:
        raise RuntimeError, "failed to lock %s: %s" % (pidfile_path, e)

    # we got the file lock, so check the pid therein
    pidfile.seek(0)
    pidfile_pid = pidfile.readline().strip()
    
    if pidfile_pid.isdigit() and pid_exists(int(pidfile_pid)):
        raise RuntimeError, ("pidfile %s contains pid (%s) of a running "
                             "process" % (lockfile, pidfile_pid))
    else:
        print ("pidfile %s contains stale pid %s; writing new lock" %
               (lockfile, pidfile_pid))

    # the pidfile didn't exist or was stale, so grab a new lock
    pidfile.truncate(0)
    pidfile.write("%d\n" % os.getpid())
    pidfile.close()

    # should be entirely unecessary, but paranoia always served me well
    confirm_lock(lockfile)
    return True


def confirm_lock(lockfile):
    """
    Confirm that the given lockfile contains our pid.
    Should be entirely unecessary, but paranoia always served me well.
    """
    pidfile = open(lockfile, "r")
    pidfile_pid = pidfile.readline().strip()
    pidfile.close()
    if int(pidfile_pid) != os.getpid():
        raise RuntimeError, ("pidfile %s contains pid %s; expected pid %s!" %
                             (lockfile, os.getpid(), pidfile_pid))
    return True
