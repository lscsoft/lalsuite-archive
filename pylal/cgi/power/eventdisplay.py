#!/usr/bin/python

import os
import popen2

from glue import lal

#
# How to run tconvert
#

class TconvertCommand(object):
	def __init__(self, tspec = None):
		self._exec = "/home/kipp/local/bin/lalapps_tconvert"
		self.tspec = tspec

	def __str__(self):
		s = self._exec
		if self.tspec:
			s += " " + self.tspec
		return s

def runtconvert(command):
	if type(command) != TconvertCommand:
		raise ValueError, "invalid argument to runtconvert(command): command must type TconvertCommand"
	child = popen2.Popen3(str(command), True)
	for line in child.childerr:
		pass
	for line in child.fromchild:
		result = lal.LIGOTimeGPS(line)
	status = child.wait()
	if not os.WIFEXITED(status) or os.WEXITSTATUS(status):
		raise Exception, "failure running \"" + str(command) + "\""
	return result
