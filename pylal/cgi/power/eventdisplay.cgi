#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()

from glue import lal
from glue import segments
from glue import segmentsUtils

import eventdisplay

#
# Init
#

now = lal.LIGOTimeGPS(eventdisplay.runtconvert(eventdisplay.TconvertCommand("now")))

default_start = "now"
default_duration = str(-1 * 3600)

G1seglist = segmentsUtils.fromlalcache(file("G1/filelist.cache"), coltype = lal.LIGOTimeGPS).coalesce()
H1seglist = segmentsUtils.fromlalcache(file("H1/filelist.cache"), coltype = lal.LIGOTimeGPS).coalesce()
H2seglist = segmentsUtils.fromlalcache(file("H2/filelist.cache"), coltype = lal.LIGOTimeGPS).coalesce()
L1seglist = segmentsUtils.fromlalcache(file("L1/filelist.cache"), coltype = lal.LIGOTimeGPS).coalesce()


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
	s = "<tr>\n"
	s += "\t<td>" + plot_pngthumbnail(name, query) + "</td>\n"
	s += "\t<td>" + plot_epslink(name, query) + "</td>\n"
	s += "</tr>"
	return s


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
	s += """<table>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="inst">Instrument:</label></td>\n"""
	s += """	<td>""" + instrumentlist(query.instrument) + """</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>GPS Time Range:</td>\n"""
	s += """	<td><label for="start">start=</label><input type="text" name="start" value=\"""" + form.getfirst("start", start) + """\"> s, <label for="dur">duration=</label><input type="text" name="dur" value=\"""" + form.getfirst("dur", "") + """\"> s</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="ratewidth">Triggers per Second Window:</label></td>\n"""
	s += """	<td><input type="text" name="ratewidth" value=\"""" + str(query.ratewidth) + """\"> s</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="freqwidth">Triggers per Hz Window:</label></td>\n"""
	s += """	<td><input type="text" name="freqwidth" value=\"""" + str(query.freqwidth) + """\"> Hz</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="lofreq">Frequency Band:</label></td>\n"""
	s += """	<td><input type="text" name="lofreq" value=\"""" + str(query.band[0]) + """\"> Hz to <input type="text" name="hifreq" value=\"""" + str(query.band[1]) + """\"> Hz</td>\n"""
	s += """</tr>\n"""
	s += """</table>\n"""
	s += """<center><input type="submit" value="Submit"></center>
</form>"""
	return s


#
# Live time info
#

def live_time_row_markup(inst, seglist):
	livetime = float(seglist.duration())
	if livetime > 0.0:
		rate = livetime / float(seglist.extent().duration())
	else:
		rate = 0.0
	s = """<tr>\n"""
	s += """	<td align="left">%s</td>\n""" % inst
	s += """	<td align="right">%.2f h</td>\n""" % (livetime / 60.0 / 60.0,)
	s += """	<td align="center">%.3f</td>\n""" % (rate,)
	s += """</tr>\n"""
	return s

def live_time_markup():
	s = """<table frame="box" rules="all">\n"""
	s += """<thead><tr>\n"""
	s += """	<th>Instruments</th>\n"""
	s += """	<th>Live Time</th>\n"""
	s += """	<th>Rate</th>\n"""
	s += """</tr></thead>\n"""
	s += """<tbody>\n"""
	s += live_time_row_markup("G1", G1seglist)
	s += live_time_row_markup("H1", H1seglist)
	s += live_time_row_markup("H2", H2seglist)
	s += live_time_row_markup("L1", L1seglist)
	s += live_time_row_markup("G1 & H1", G1seglist & H1seglist)
	s += live_time_row_markup("G1 & H2", G1seglist & H2seglist)
	s += live_time_row_markup("G1 & L1", G1seglist & L1seglist)
	s += live_time_row_markup("H1 & H2", H1seglist & H2seglist)
	s += live_time_row_markup("H1 & L1", H1seglist & L1seglist)
	s += live_time_row_markup("H2 & L1", H2seglist & L1seglist)
	s += live_time_row_markup("G1 & H1 & H2", G1seglist & H1seglist & H2seglist)
	s += live_time_row_markup("G1 & H1 & L1", G1seglist & H1seglist & L1seglist)
	s += live_time_row_markup("G1 & H2 & L1", G1seglist & H2seglist & L1seglist)
	s += live_time_row_markup("H1 & H2 & L1", H1seglist & H2seglist & L1seglist)
	s += live_time_row_markup("G1 & H1 & H2 & L1", G1seglist & H1seglist & H2seglist & L1seglist)
	s += """</tbody>\n"""
	s += """</table>"""
	return s

def s5_live_time_summary(seglist):
	s5length = 1.0 * 365.25 * 24.0 * 60.0 * 60.0	# 1 year
	livetime = float(seglist.duration())
	if livetime > 0.0:
		rate = livetime / float(seglist.extent().duration())
	else:
		rate = 0.0
	s5end = now + (s5length - livetime) / rate
	s = """<table>\n"""
	s += """<tr>\n"""
	s += """	<td>H1 & H2 & L1 hours required for S5</td>\n"""
	s += """	<td>%.2f h</td>\n""" % (s5length / 60.0 / 60.0,)
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>Fraction of S5 completed</td>\n"""
	s += """	<td>%.1f%%</td>\n""" % (100.0 * livetime / s5length,)
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>Estimated S5 termination</td>\n"""
	s += """	<td>%s (GPS %d s)</td>\n""" % (eventdisplay.runtconvert(eventdisplay.TconvertCommand(str(int(s5end)))), int(s5end))
	s += """</tr>\n"""
	s += """</table>\n"""
	return s


#
# Excess Power Event Interface
#

print "Content-Type: text/html\n"

print """<html>"""

print refresh

print "<center><h1>Excess Power Event Interface</h1></center>"
print "<center>(Patience Required)</center>"
print """<p>You can <a href="http://www.lsc-group.phys.uwm.edu/cgi-bin/cvs/viewcvs.cgi/pylal/cgi/power/?cvsroot=lscsoft">browse the source code for these web pages</a>.</p>"""
print "<p>"
print formmarkup(query)
print "</p>"

print """<hr width="90%">"""
print "<center>"
print "<h2>S5 Live Times To Date</h2>"
print live_time_markup()
print s5_live_time_summary(H1seglist & H2seglist & L1seglist)
print "</center>"
print """<hr width="90%">"""

if query.segment.duration() > 24 * 3600:
	# Time interval too long error
	print errormsg("Requested segment is too long (24 hour max)")
else:
	# FIXME: doesn't work because user apache doesn't have LSC grid
	# credentials.
	#print str(eventdisplay.getsegments(eventdisplay.SegFindConfigH2, "Science", segments.segment(eventdisplay.s5start, now)))

	# Table of plots
	print "<center>"
	print "<h2>%s s Starting At %s</h2>" % (duration, start.title())
	print "<p>"
	print "<table>"
	print plot_table_row("rateplot.cgi", query)
	print plot_table_row("conf_vs_time.cgi", query)
	print plot_table_row("tfplot.cgi", query)
	print plot_table_row("rate_vs_freq.cgi", query)
	print plot_table_row("conf_vs_freq.cgi", query)
	print plot_table_row("rate_vs_conf.cgi", query)
	print plot_table_row("rate_vs_snr.cgi", query)
	print "</table>"
	print "</p>"
	print "</center>"

print """</html>"""
