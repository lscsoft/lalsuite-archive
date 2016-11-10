#!/usr/bin/env python

import doctest
import sys
from glue.ligolw import ligolw

if __name__ == '__main__':
	failures = doctest.testmod(ligolw)[0]
	sys.exit(bool(failures))
