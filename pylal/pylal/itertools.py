"""
This code has been moved to the iterutils module in Glue.
"""


import warnings


warnings.warn("The pylal.itertools module is deprecated, the code has been migrated to the glue.iterutils module.  This is only a warning, your program will continue to work, but please update your code to import the new Glue module instead because the pylal.itertools module will be deleted in the near future.", DeprecationWarning)


from glue.iterutils import *
