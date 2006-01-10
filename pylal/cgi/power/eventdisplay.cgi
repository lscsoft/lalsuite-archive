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
	freqwidth = None
	band = None
	instrument = None

form = cgi.FieldStorage()

query = Query()

start = form.getfirst("start", default_start).lower()
duration = lal.LIGOTimeGPS(form.getfirst("dur", default_duration))
if start == "now":
	query.segment = segments.segment(now, now + duration)
	refresh = """<meta http-equiv="refresh" content="%d"></meta>""" % (abs(duration) / 3600 * 180 + 180)
else:
	query.segment = segments.segment(lal.LIGOTimeGPS(start), lal.LIGOTimeGPS(start) + duration)
	refresh = ""

query.ratewidth = float(form.getfirst("ratewidth", "60"))
query.freqwidth = float(form.getfirst("freqwidth", "16"))
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

def _imgsrc(name, query):
	return "%s?inst=%s&start=%s&dur=%s&ratewidth=%s&freqwidth=%s&lofreq=%s&hifreq=%s" % (name, query.instrument, query.segment[0], query.segment.duration(), query.ratewidth, query.freqwidth, query.band[0], query.band[1])

def plot_pngthumbnail(name, query):
	src = _imgsrc(name, query) + "&format=png"
	return """<a href="%s"><img src="%s" width="800"></a>""" % (src, src)

def plot_epslink(name, query):
	src = _imgsrc(name, query) + "&format=eps"
	return """<a href="%s">EPS</a>""" % src

def plot_table_row(name, query):
	markup = "\t<td>" + plot_pngthumbnail(name, query) + "</td>"
	markup += "\t<td>" + plot_epslink(name, query) + "</td>"
	return markup


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
	s += """<label for="start">GPS Time Range: start=</label><input type="text" name="start" value=\"""" + form.getfirst("start", start) + """\"> s, duration=<input type="text" name="dur" value=\"""" + form.getfirst("dur", "") + """\"> s<br>\n"""
	s += """<label for="ratewidth">Triggers per Second Window: </label><input type="text" name="ratewidth" value=\"""" + str(query.ratewidth) + """\"> s<br>\n"""
	s += """<label for="freqwidth">Triggers per Hz Window: </label><input type="text" name="freqwidth" value=\"""" + str(query.freqwidth) + """\"> Hz<br>\n"""
	s += """<label for="lofreq">Frequency Band: </label><input type="text" name="lofreq" value=\"""" + str(query.band[0]) + """\"> Hz to <input type="text" name="hifreq" value=\"""" + str(query.band[1]) + """\"> Hz<br>\n"""
	s += """<center><input type="submit" value="Submit"></center>
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
	# Table of plots
	print "<center>"
	print "<h2>%s s Starting At %s</h2>" % (duration, start.title())
	print "<table>"
	print "<tr>" + plot_table_row("rateplot.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("conf_vs_time.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("tfplot.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("rate_vs_freq.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("conf_vs_freq.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("rate_vs_conf.cgi", query) + "</tr>"
	print "<tr>" + plot_table_row("rate_vs_snr.cgi", query) + "</tr>"
	print "</table>"
	print "</center>"

print """</html>"""
