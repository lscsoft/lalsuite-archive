import doctest
import sys
from glue.ligolw import _ilwd
from glue.ligolw import ilwd

failures = doctest.testmod(_ilwd)[0]
failures += doctest.testmod(ilwd)[0]

sys.exit(bool(failures))
