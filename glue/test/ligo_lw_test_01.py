import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy
import sys

from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import utils

xmldoc = utils.load_filename("ligo_lw_test_01.xml", verbose = True)

for n, a in enumerate(xmldoc.getElementsByTagName(ligolw.Array.tagName)):
	print >>sys.stderr, "found %s array '%s'" % ("x".join(map(str, a.array.shape)), a.getAttribute("Name"))
	fig = figure.Figure()
	FigureCanvas(fig)
	axes = fig.gca()
	axes.loglog()
	axes.grid(True)
	for i in range(1, a.array.shape[0]):
		axes.plot(numpy.fabs(a.array[0]), numpy.fabs(a.array[i]))
	axes.set_title(a.getAttribute("Name"))
	print >>sys.stderr, "saving as 'ligo_lw_test_01_%d.png' ..." % n
	fig.savefig("ligo_lw_test_01_%d.png" % n)
	print >>sys.stderr, "done."
