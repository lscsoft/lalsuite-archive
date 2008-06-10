
import glue.lal
from subprocess import Popen, PIPE
from urlparse import urlsplit
import urllib2
import os
import sys
from tempfile import mkdtemp
from socket import gethostbyaddr, gethostname

__all__ = ["Cache", "CacheEntry"]

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
