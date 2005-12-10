#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()
import shutil
import sys

form = cgi.FieldStorage()

filename = form.getfirst("filename")

if filename[:5] != "/tmp/" or filename[-4:] != ".png":
	sys.exit(0)

print >>sys.stdout, "Content-Type: image/png\n"
shutil.copyfileobj(file(form.getfirst("filename")), sys.stdout)
