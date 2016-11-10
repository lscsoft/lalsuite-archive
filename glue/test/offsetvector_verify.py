import doctest
import sys
from glue import offsetvector

failures = doctest.testmod(offsetvector)[0]

sys.exit(bool(failures))
