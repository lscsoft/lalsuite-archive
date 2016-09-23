from setuptools import setup

version = "1.0.2"

setup(
  name = "ligo-common",
  version = version,
  maintainer = "Brian Moe",
  maintainer_email = "brian.moe@ligo.org",
  description = "Empty LIGO modules",
  long_description = "Empty module placeholder for other LIGO modules",
  license = 'GPL',
  namespace_packages = ['ligo'],
  provides = ['ligo'],
  packages = [ 'ligo'],
)

