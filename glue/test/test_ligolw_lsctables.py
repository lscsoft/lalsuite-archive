#!/usr/bin/env python

import doctest
import sys
from glue.ligolw import lsctables

if __name__ == '__main__':
	failures = doctest.testmod(lsctables)[0]
	sys.exit(bool(failures))
