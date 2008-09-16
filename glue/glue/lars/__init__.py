
import glue.lal
from subprocess import Popen, PIPE
from urlparse import urlsplit
import urllib2
import os
import sys
from tempfile import mkdtemp
from socket import gethostbyaddr, gethostname
import logging
import pickle

__all__ = ["Cache", "CacheEntry", "initializeCli", "initializeLogging"]

DEFAULT_DB = "cvs::pserver:gravity.phys.uwm.edu:2401/usr/local/cvs/larsdb?file=inspiral/advertising.cache"

config = {}

def initializeLogging(defaultLogLevel=logging.WARNING):
    # Initialize the logging facility
    #  XXX This is very crude.  Output looks ugly.
    logging.basicConfig(level=defaultLogLevel)
    return logging.getLogger('lars')


def initializeCli():
    # Default database -- in some ini or something...  ??

    global config

    if os.environ.has_key('LARS_LIB'):
        sys.path.insert(0, os.environ['LARS_LIB'])

    larsDbUrl = os.environ.get("LARS_DB") or DEFAULT_DB

    if os.environ.has_key("LARS_DIR"):
        basedir = os.environ["LARS_DIR"]
    else:
        basedir = os.path.join(os.environ["HOME"], ".lars")

    mountdir = os.path.join(basedir, "mnt")
    tmpdir = os.path.join(basedir, "tmp")

    # Make sure all lars-y directories exist and are accessible.
    #

    if not os.access(basedir,os.F_OK):
      if verbose: print "Creating .lars directory hierarchy at " + basedir
      try:
        os.makedirs(basedir)
      except:
        print "Directory " + basedir + " does not exist and could not be created"
        sys.exit(1)
      try:    
        os.mkdir(mountdir)
      except: 
        print "Directory " + mountdir + " does not exist and could not be created"
        sys.exit(1)
      try:    
        os.mkdir(basedir+"/tmp")
      except: 
        print "Directory " + mountdir + " does not exist and could not be created"
        sys.exit(1)
    # make the above things exist now

    config['larsDbUrl'] = larsDbUrl
    config['basedir']   = basedir
    config['mountdir']  = mountdir
    config['tmpdir']    = tmpdir

    return config

CacheEntry = glue.lal.CacheEntry

class LarsBadPasswordException(Exception): pass

class Cache(glue.lal.Cache):

    def getSearch(cls, searchDirUrl):
        # retrieve all cache files for a given search
        (scheme, netloc, path, query, frag) = urlsplit(searchDirUrl)
        if scheme == "file":
            filePattern = os.path.join(path, "*.cache")
            if netloc == "localhost" or islocalhost(netloc):
                #print "netloc (%s) is local" % netloc
                #print "looking for cache(s): %s" % filePattern
                p = Popen(["cat %s"% filePattern],
                           shell=True, stdin=open("/dev/null"),
                           stdout=PIPE, close_fds=True)
            else:
                # assuming no port or username in netloc and query and frag are nil
                p = Popen(["ssh", netloc, "cat", filePattern],
                           stdin=open("/dev/null"),
                           stdout=PIPE, close_fds=True)
            retval = cls.fromfile(p.stdout)
            p.stdout.close()
            del p
            return retval
        else:
            raise Exception("unknown scheme '%s' in '%s'" % (scheme, searchDirUrl))

    getSearch = classmethod(getSearch)

    def get(cls, url):
        (scheme, netloc, path, query, frag) = urlsplit(url)
        if scheme == "file":
            if islocalhost(netloc):
                rv = cls.fromfile(open(path,"r"))
            else:
                # ugh
                p = Popen(["ssh", netloc, "cat", path], stdout=PIPE, close_fds=True)
                rv = cls.fromfile(p)
        elif scheme == "cvs":
            # Assuming format is ":pserver:/loc/of/repo?file=DBFILE"
            f = CvsFile(url[4:])
            rv = cls.fromfile(f)
        elif scheme == "http":
            f = urllib2.urlopen(url)
            rv = cls.fromfile(f)
        else:
            raise Exception("Unknown scheme: %s" % scheme)
        return rv
    get = classmethod(get)

    def writeToUrl(self, url):
        (scheme, netloc, path, query, frag) = urlsplit(url)
        if scheme == "file":
            if netloc == "localhost":
                f = open(path, "w")
            else:
                p = Popen(["ssh", netloc, "cat", ">", path],
                           stdin=PIPE, close_fds=True)
                f = p.stdin
        if scheme == "cvs":
            f = CvsFile(url[4:], "w")
        else:
            raise Exception("Unknown scheme: %s" % scheme)
        for entry in self:
            print >>f, str(entry)
        f.close()

class CvsFile(file):
    def __init__(self, cvspath, perms="r"):
        self._cvstmpdir = mkdtemp()
	self._loginattempts = 4
        if not cvspath.startswith(":pserver:"):
            raise Exception("cannot parse cvsroot '%s'" % cvspath)
        try: 
            breakpt = cvspath.find("?")
            self._cvsroot = cvspath[:breakpt]
            (key, val) = cvspath[breakpt+1:].split("=")
            if key != "file": pass #this is really an error
            self._cvsfname = val
        except Exception, e:
            raise Exception(e, "cannot parse cvsroot '%s'" % cvspath)
        self._cvsop("checkout")
        self._cvsfpath = os.path.join(self._cvstmpdir, self._cvsfname)
        file.__init__(self, self._cvsfpath, perms)

    def _cvsop(self, op, arg=None, trylogin=True):
        if not isinstance(op, list):
            op = [op]
	if not arg:
	    arg = [self._cvsfname]
	if not isinstance(arg, list):
	    arg = [arg]
        p = Popen(["cvs", "-d", self._cvsroot] + op + arg,
                  cwd=self._cvstmpdir,
                  stdout=open("/dev/null","rw"),
		  stderr=PIPE)
        p.wait()

	if ( p.returncode != 0 ):
	    if "cvs login" in p.stderr.read() and trylogin and sys.stdin.isatty():
	        for i in range(0, self._loginattempts):
		    try:
		        self._cvsop("login", trylogin=False)
			try:
			    self._cvsop(op, arg, trylogin=False)
			except: raise
		        return
	    	    except:
		        pass
	        raise LarsBadPasswordException()
	    else:
	        raise Exception("Problem with cvs")

    def __del__(self):
        # XXX "change" should be something better
	try:
	    self._cvsop(["commit", "-m", "change"], trylogin=False)
        except: pass
        for root, dirs, files in os.walk(self._cvstmpdir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(self._cvstmpdir)


def islocalhost(hostname):
    localhostname, aliases, ips = gethostbyaddr(gethostname())
    localhostnames = ["localhost", localhostname, "127.0.0.0"] + aliases + ips
    return hostname in localhostnames

# From the man page.
rsyncReturnCodes = {
      0  :   "Success",
      1  :   "Syntax or usage error",
      2  :   "Protocol incompatibility",
      3  :   "Errors selecting input/output files, dirs",
      4  :   "Requested  action  not supported: an attempt was made to manipulate 64-bit files on a platform that cannot support them; or  an  option was  specified  that  is  supported  by  the  client and not by the server.",
      5  :   "Error starting client-server protocol",
      6  :   "Daemon unable to append to log-file",
      10 :   "Error in socket I/O",
      11 :   "Error in file I/O",
      12 :   "Error in rsync protocol data stream",
      13 :   "Errors with program diagnostics",
      14 :   "Error in IPC code",
      20 :   "Received SIGUSR1 or SIGINT",
      21 :   "Some error returned by waitpid()",
      22 :   "Error allocating core memory buffers",
      23 :   "Partial transfer due to error",
      24 :   "Partial transfer due to vanished source files",
      25 :   "The --max-delete limit stopped deletions",
      30 :   "Timeout in data send/receive",}

def copyDir(src, dest): # XXX should go into lars.util

    log = logging.getLogger("lars.util")
    log.debug("copying directory %s to %s" % (src, dest))

    (src_scheme,  src_netloc,  src_path, _, _)  = urlsplit(src)
    (dest_scheme, dest_netloc, dest_path, _, _) = urlsplit(dest)

    if not src_netloc:
        src_netloc = "localhost"

    if src_scheme and src_scheme != "file":
        raise Exception("Unknown scheme in directory transfer source URL: %s" % src)
    if src_scheme and src_scheme != "file":
        raise Exception("Unknown scheme in directory transfer destination URL: %s" % dest)
    if dest_netloc and not islocalhost(dest_netloc):
        raise Exception("Cannot copyDir to non local directory: %s" % dest)

    rsyncCommand = ["rsync", "-r", "%s:%s" % (src_netloc, src_path), os.path.abspath(dest_path)]
    log.debug("Rsync command: %s" % rsyncCommand)

    p = Popen(rsyncCommand,
              stdin=open("/dev/null", "rw"),
              stdout=open("/dev/null", "rw"),
              #stderr=open("/dev/null", "rw"),
              close_fds=True
             )
    p.wait()
    log.debug("Rsync completed with return code: %s" % str(p.returncode))

    if p.returncode:
        if p.returncode in [23,24]:
            log.warn("Rsync problem: %s" % rsyncReturnCodes[p.returncode])
        else:
            log.error("Rsync problem: %s" % rsyncReturnCodes[p.returncode])
            raise Exception(rsyncReturnCodes[p.returncode])


def copyFile(src, dest):
    # lame -- not robust or flexible.  way too many assumptions
    (scheme, netloc, path, _, _) = urlsplit(src)
    # should at least check the scheme. and if we could open the out file
    destFile = open(dest, "w")
    p = Popen(["ssh", netloc, "cat", path],
	      stdout=destFile,
              stderr=open("/dev/null", "rw"),
	      close_fds=True)
    p.wait()
    destFile.close()
    if p.returncode != 0:
        # SHOULD BE LOGGED NOT PRINTED
        print "Could not copy: ", src
        os.unlink(dest)

class SshFsManager:
    def __init__(self, mountdir=None):
        global config
        self.log = logging.getLogger("lars.SshFsManager")
        mountdir = mountdir or config.get("mountdir")
        if not mountdir:
            raise Exception("SshFsManager: no mount directory")
        self.mountdir = mountdir
        # XXX make sure mountdir exists, is a directory, is read/writable

        self.mountinfo = {}
        try:
            self.infofile = open(os.path.join(mountdir, "mount.info"), "rb")
            self._readMountInfo()
        except IOError:
            pass

    def _readMountInfo(self):
        self.mountinfo = pickle.load(self.infofile)

    def _writeMountInfo(self):
        # XXX We have some locking/concurrency issues here...
        # For now we shall ignore them.
        self.mountinfo = pickle.load(self.infofile)
        outfile = open(self.infofile, "wb")
        pickle.dump(self.mountinfo, outfile)
        outfile.close()

    def mount(self, srcUrl, destDir):
        (scheme, netloc, path, _, _) = urlsplit(srcUrl)
        fs = "%s:%s" % (netloc, path)
        try:    os.mkdir(destDir)
        except: pass
        self.log.debug("Mounting '%s' on '%s'" % (fs, destDir))
        p = Popen(["sshfs", fs, os.path.abspath(destDir)], cwd=self.mountdir)
        p.wait()
        del p

