
import os

from distutils.core import setup

version = "1.1.0"

setup(
  name = "ligo-lars",
  version = version,
  maintainer = "Brian Moe",
  maintainer_email = "brian.moe@ligo.org",
  description = "LIGO Archiving Service",
  long_description = "Lars is a protype for an archival service for LIGO",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/lars.html",
  license = 'GPL',
  provides = ['ligo.lars'],
  packages = [ 'ligo.lars', 'ligo.lars.cli'],

  requires = ['ligo', 'glue.segments', 'M2Crypto'],

  scripts = [
    os.path.join('bin','lars'),
    os.path.join('bin','lars_add'),
    os.path.join('bin','lars_search'),
  ],

)
