#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()
import shutil
import sys

from glue import lal
from glue import segments
from glue import segmentsUtils

from eventdisplay import *

#
# Init
#

now = runtconvert(TconvertCommand("now"))

default_start = "now"
default_duration = str(-8 * 3600)


#
# Parse display request
#

class Query(object):
	segment = None
	ratewidth = None
	band = None
	instrument = None
	cache = None

form = cgi.FieldStorage()

query = Query()

start = form.getfirst("start", default_start).lower()
duration = lal.LIGOTimeGPS(form.getfirst("dur", default_duration))
if start == "now":
	query.segment = segments.segment(now, now + duration)
	refreshtag = """<meta http-equiv="refresh" content="600"></meta>"""
else:
	query.segment = segments.segment(lal.LIGOTimeGPS(start), lal.LIGOTimeGPS(start) + duration)
	refreshtag = ""

query.ratewidth = float(form.getfirst("ratewidth", "60"))
query.band = segments.segment(float(form.getfirst("lofreq", "0")), float(form.getfirst("hifreq", "2500")))
query.instrument = form.getfirst("inst", "H1")
query.cache = "/home/kipp/cgi-bin/" + query.instrument + "/filelist.cache"


#
# An error message
#

def errormsg(msg):
	return """<center>Error: %s</center>""" % msg


#
# Plot markup
#

def plotmarkup(filename):
	src = "image.cgi?filename=" + filename
	return """<a href="%s"><img src="%s" width="500"></a>""" % (src, src)


#
# Form markup
#

def formmarkup(query):
	def instrumentlist(default):
		s = """<select name="inst"><option>""" + default + """</option>"""
		for inst in [inst for inst in ["G1", "H1", "H2", "L1"] if inst != default]:
			s += "<option>" + inst + "</option>"
		return s + "</select>"

	s = """<form action="eventdisplay.cgi" method="get">\n"""
	s += """<label for="inst">Instrument: </label>""" + instrumentlist(query.instrument) + "<br>\n"
	s += """<label for="start">GPS Time Range: start=</label><input type="text" name="start" value=\"""" + form.getfirst("start", start) + """\"> s, duration=<input type="text" name="dur" value=\"""" + form.getfirst("dur", "") + """\"> s<br>\n<label for="ratewidth">Trigger Rate Window Width: </label><input type="text" name="ratewidth" value=\"""" + str(query.ratewidth) + """\"> s<br>\n<label for="lofreq">Frequency Band: </label><input type="text" name="lofreq" value=\"""" + str(query.band[0]) + """\"> Hz to <input type="text" name="hifreq" value=\"""" + str(query.band[1]) + """\"> Hz<br>\n<center><input type="submit" value="Submit"></center>
</form>"""
	return s


#
# Testing...
#

seglist = segmentsUtils.fromlalcache(file(query.cache), coltype=lal.LIGOTimeGPS).coalesce()

print "Content-Type: text/html\n"

print """<html>"""

print refreshtag

print "<h1>Excess Power Event Interface</h1>"
print formmarkup(query)

if query.segment.duration() > 24 * 3600:
	# Time interval too long error
	print errormsg("Requested segment is too long (24 hour max)")
else:
	print "<h2>%s s Starting At %s</h2>" % (duration, start.title())

	# Make plots
	tfplotdesc = TFPlotDescription(query.instrument, query.segment, query.band, seglist)
	rateplotdesc = RatePlotDescription(query.instrument, query.ratewidth, query.segment, seglist)

	table = gettriggers(query.cache, tfplotdesc.trig_segment() | rateplotdesc.trig_segment())

	maketfplot(tfplotdesc, table)
	makerateplot(rateplotdesc, table)

	print plotmarkup(rateplotdesc.filename)
	print plotmarkup(tfplotdesc.filename)

print """</html>"""
