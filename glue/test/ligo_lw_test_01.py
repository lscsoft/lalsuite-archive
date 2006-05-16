import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

from glue.ligolw import ligolw
from glue.ligolw import array

print "Loading ligo_lw_test_01.xml..."
doc = ligolw.Document()
ligolw.make_parser(array.LIGOLWContentHandler(doc)).parse(file("ligo_lw_test_01.xml"))

for n, a in enumerate(doc.getElementsByTagName(ligolw.Array.tagName)):
	print "Found %s array \"%s\"..." % ("x".join(map(str, a.array.shape)), a.getAttribute("Name")),
	fig = figure.Figure()
	FigureCanvasAgg(fig)
	fig.set_figsize_inches(16, 8)
	axes = fig.gca()
	axes.loglog()
	axes.grid(True)
	for i in range(1, a.array.shape[0]):
		axes.plot(numarray.fabs(a.array[0]), numarray.fabs(a.array[i]))
	axes.set_title(a.getAttribute("Name"))
	fig.savefig("ligo_lw_test_01_%d.png" % n)
	print "saved as ligo_lw_test_01_%d.png" % n
