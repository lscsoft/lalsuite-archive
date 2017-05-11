#!/usr/bin/env python

import doctest
import sys
from glue.ligolw import array

if __name__ == '__main__':
	failures = doctest.testmod(array)[0]
	sys.exit(bool(failures))
