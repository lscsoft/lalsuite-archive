"""
This module defines the some useful functions for GSI socket servers.
"""

__author__ = "Scott Koranda <skoranda@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import os
import sys

class Gridmap(object):
  """
  Class to represent the grid-mapfile used for the server. 
  Methods __getitem__ and __setitem__ are added for convenience
  so that an instance can be treated like a dictionary.

  Entries in the grid-mapfile should be on a single line and
  surrounded by double quotes. No mapping to a local user ID 
  is necessary.
  """
  def __init__(self, filePath=None, logger=None):
    """
    Initialize a new instance. If filePath is None then check
    the environment for the variable GRIDMAP and uses the 
    path found as the path to the grid-mapfile to open and parse.
    
    @return: an instance that can be treated like a dictionary
    """
    if not filePath:
      try: 
        path = os.environ["GRIDMAP"]
      except:
        if logger:
          logger.error("Error parsing GRIDMAP from environment: %s" % e)
          logger.error("Authentication to server will fail for all subjects")
    else:   
      path = filePath
      self.path = path
      self.logger = logger


  def __getitem__(self, key):
    f = None
    try:
      f = open(self.path, 'r')
    except Exception, e:
      if self.logger:
        self.logger.error("Error opening grid-mapfile %s: %s" % (self.path, e))
        self.logger.error("Authentication to server will fail for all subjects")
                
    try:
      subjects = []
      lines = [ s.strip() for s in f.readlines()]
      for s in lines:
        a = s.split('"')
        for b in a:
          if b == '': continue
          else:
            subjects.append(b)
            break
    except Exception, e:
      if self.logger:
        self.logger.error("Error parsing grid-mapfile %s: %s" % (self.path, e))
        self.logger.error("Authentication to server may fail for some subjects")

    self.d = {} 
    for s in subjects:
      self.d[s] = 1

    if f:
      f.close()

    if self.d.has_key(key): return 1
    else: return 0


class AuthCallback(object):
  def __init__(self, gridmap, logger, callback=None):
    self.gridmap = gridmap
    self.logger = logger
    self.callback = callback

  def __call__(self, arg, handle, identity, context):
    """
    Callback function for GSI authentication. Each time this 
    function is called a new instance of the GridMap class is created
    that reads the grid-mapfile. If the subject passed in as the identity
    is in the grid-mapfile then 1 is returned and the connection is
    authenticated and allowed to continue. If the subject passed is
    not in grid-mapfile then a 0 is returned and the connection is
    broken.

    This function is supplied to help create an instance of io.AuthData
    which is then supplied to help create an instance of io.TCPIOAttr,
    which is then supplied to help create the instance of 
    io.ThreadingGSITCPSocketServer. See the
    U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>} documentation.

    Everything is within a try/except clause since if an exception is raised
    within this function experience shows that the connection will be
    allowed to authenticate!

    @param logger:  an instance of a logger to which the identity can be
      logged, or None
    @param handle: handle to the io.AuthData instance, not used here
    @param identity: the identity of the party trying to authenticate
    @type identity: string
    @param context: the GSI context, not used here

    @return: 1 to accept and 0 to reject
    """
    try:
      # load gridmap file
      if self.gridmap[identity]:
        if self.logger:
          self.logger.info("Accepted connection from %s" % identity)
        if self.callback:
          self.callback(arg, 1)
        return 1
      else:
        if self.logger:
          self.logger.info("Rejected connection from %s" % identity)
        if self.callback:
          self.callback(arg, 0)
        return 0
    except Exception, e:
      try:
        self.logger.error("Error within authentication callback: %s" % e)
      except Exception, e:
        pass
      return 0


def daemon():
  """
  Common code for making a process into a daemon. 

  Fork once, change directory to /, get a new session, set the umask to 0,
  then fork again.

  @return: Parent never returns. The child gets None.
  """
  # do the first fork
  try:
    pid = os.fork()
    if pid > 0:
      # exit first parent
      sys.exit(0)
  except OSError, e:
    print >>sys.stderr, "fork #1 failed: %d (%s)" % (e.errno, e.strerror)
    sys.exit(1)
        
  # decouple from parent environment
  os.chdir("/")
  os.setsid()
  os.umask(0)
        
  # do the second fork
  try:
    pid = os.fork()
    if pid > 0:
      # exit from second parent
      sys.exit(0)
  except OSError, e:
    print >>sys.stderr, "fork #2 failed: %d (%s)" % (e.errno, e.strerror)
    sys.exit(1)
        
  # send stdout and stderr to /dev/null
  try:
    devnull = open("/dev/null", "w")
    sys.stdout = devnull
    sys.stderr = devnull
  except Exception, e:
    print >>sys.__stderr__, "Unable to direct to /dev/null: %s" % e
    sys.exit(1)

