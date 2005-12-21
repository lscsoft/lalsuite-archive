#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()

from glue import lal
from glue import segments

import eventdisplay

#
# Init
#

now = eventdisplay.runtconvert(eventdisplay.TconvertCommand("now"))

default_start = "now"
default_duration = str(-1 * 3600)


#
# Parse display request
#

class Query(object):
	segment = None
	ratewidth = None
	band = None
	instrument = None

form = cgi.FieldStorage()

query = Query()

start = form.getfirst("start", default_start).lower()
duration = lal.LIGOTimeGPS(form.getfirst("dur", default_duration))
if start == "now":
	query.segment = segments.segment(now, now + duration)
	refresh = """<meta http-equiv="refresh" content="%d"></meta>""" % (abs(duration) / 3600 * 60 + 60)
else:
	query.segment = segments.segment(lal.LIGOTimeGPS(start), lal.LIGOTimeGPS(start) + duration)
	refresh = ""

query.ratewidth = float(form.getfirst("ratewidth", "60"))
query.band = segments.segment(float(form.getfirst("lofreq", "0")), float(form.getfirst("hifreq", "2500")))
query.instrument = form.getfirst("inst", "H1")


#
# An error message
#

def errormsg(msg):
	return """<center>Error: %s</center>""" % msg


#
# Plot markup
#

def plotmarkup(name, query):
	src = "%s?inst=%s&start=%s&dur=%s&ratewidth=%s&lofreq=%s&hifreq=%s" % (name, query.instrument, query.segment[0], query.segment.duration(), query.ratewidth, query.band[0], query.band[1])
	return """<a href="%s"><img src="%s" width="800"></a>""" % (src, src)


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
# Excess Power Event Interface
#

print "Content-Type: text/html\n"

print """<html>"""

print refresh

print "<h1>Excess Power Event Interface</h1>"
print formmarkup(query)

if query.segment.duration() > 24 * 3600:
	# Time interval too long error
	print errormsg("Requested segment is too long (24 hour max)")
else:
	print "<h2>%s s Starting At %s</h2>" % (duration, start.title())

	# Include plots
	print plotmarkup("rateplot.cgi", query)
	print plotmarkup("conf_vs_time.cgi", query)
	print plotmarkup("tfplot.cgi", query)
	print plotmarkup("conf_vs_freq.cgi", query)
	print plotmarkup("rate_vs_freq.cgi", query)
	print plotmarkup("rate_vs_conf.cgi", query)

print """</html>"""
