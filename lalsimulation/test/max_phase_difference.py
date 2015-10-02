#!/usr/bin/python

import numpy
import scipy.interpolate
import sys

orig_fn = sys.argv[1]
current_fn = sys.argv[2]

orig = numpy.loadtxt(orig_fn)
current = numpy.loadtxt(current_fn)

tck = scipy.interpolate.splrep(current[:,0], current[:,1], s=0)
phase_interpolated = scipy.interpolate.splev(orig[:,0], tck, der=0)

maxdiff = numpy.amax(numpy.abs(phase_interpolated-orig[:,1]))

print "# max phase difference: %e" % maxdiff
