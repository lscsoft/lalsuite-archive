
import os
from urlparse import urlsplit
import logging

log = logging.getLogger("lars.util")

def mountSearch(searchUrl, targetDir):
    (scheme, netloc, path, query, frag) = urlsplit(searchUrl)
    # netloc is just host -- we didn't put user/port in there right?
    fs = "%s:%s" % (netloc, path)
    try:    os.mkdir(targetDir)
    except: pass
    log.info("Mounting '%s' on '%s'" % (fs, targetDir))
    p = Popen(["sshfs", fs, os.path.abspath(targetDir)], cwd=basedir)
    p.wait()
    del p
