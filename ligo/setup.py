
from distutils.core import setup

version = "1.0.1"

setup(
  name = "ligo-common",
  version = version,
  maintainer = "Brian Moe",
  maintainer_email = "brian.moe@ligo.org",
  description = "Empty LIGO modules",
  long_description = "Empty module placeholder for other LIGO modules",
  license = 'GPL',
  provides = ['ligo'],
  packages = [ 'ligo'],
)

