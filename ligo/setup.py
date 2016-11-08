from setuptools import setup

version = "1.0.2"

setup(
  name = "ligo-common",
  version = version,
  maintainer = "Tanner Prestegard",
  maintainer_email = "tanner.prestegard@ligo.org",
  description = "Empty LIGO modules",
  long_description = "Empty module placeholder for other LIGO modules",
  license = 'GPL',
  provides = ['ligo'],
  packages = [ 'ligo'],
)

