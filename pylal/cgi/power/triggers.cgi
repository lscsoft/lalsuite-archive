#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()

from glue import lal

import eventdisplay
import webplot


#
# Hack PlotDescription class from webplot.py to describe the output document
#

class Description(webplot.PlotDescription):
	def __init__(self):
		webplot.PlotDescription.__init__(self)
		self.set_format("xml")
		return self

	def parse_form(self):
		webplot.PlotDescription.parse_form(self)
		self.set_format("xml")
		return self


#
# Generate and send trigger file.
#

description = Description().parse_form()

command = eventdisplay.LLAddCommand(webplot.CacheURLs(description.cache, description.segment), output = description.filename)

eventdisplay.runlladd(command)

webplot.SendImage(description)
