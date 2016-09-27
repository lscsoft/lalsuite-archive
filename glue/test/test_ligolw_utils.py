#!/usr/bin/env python

import doctest
import sys
from glue.ligolw import utils

if __name__ == '__main__':
	failures = doctest.testmod(utils)[0]
	sys.exit(bool(failures))
