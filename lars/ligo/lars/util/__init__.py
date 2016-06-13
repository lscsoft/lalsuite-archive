
import os
from urlparse import urlsplit
import logging
from subprocess import Popen

log = logging.getLogger("lars.util")

def mountSearch(searchUrl, targetDir):
    (scheme, netloc, path, query, frag) = urlsplit(searchUrl)
    # netloc is just host -- we didn't put user/port in there right?
    fs = "%s:%s" % (netloc, path)
    targetDir = os.path.abspath(os.path.expanduser(targetDir))
    try:    os.mkdir(targetDir)
    except: pass
    log.info("Mounting '%s' on '%s'" % (fs, targetDir))
    #p = Popen(["sshfs", fs, os.path.abspath(targetDir)], cwd=basedir)
    p = Popen(["sshfs", fs, os.path.abspath(targetDir)])
    p.wait()
    del p
