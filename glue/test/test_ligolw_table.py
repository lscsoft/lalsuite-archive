#!/usr/bin/env python

import doctest
import sys
from glue.ligolw import table

if __name__ == '__main__':
	failures = doctest.testmod(table)[0]
	sys.exit(bool(failures))
