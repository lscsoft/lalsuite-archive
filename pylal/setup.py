# $Id$
# 
# setup for lgen

import os
from distutils.core import setup, Extension

def stripfirsttwo(string):
	return string[2:]

lallibs = map(stripfirsttwo, os.popen("pkg-config lal lalmetaio --libs-only-l").read().split())

lallibdirs = map(stripfirsttwo, os.popen("pkg-config lal --libs-only-L").read().split())
lallibdirs = lallibdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --libs-only-L").read().split())

lalincdirs = map(stripfirsttwo, os.popen("pkg-config lal --cflags-only-I").read().split())
lalincdirs = lalincdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --cflags-only-I").read().split())


setup(
	name = "lgen",
	version = "0.1",
	author = "Patrick Brady",
	author_email = "patrick@gravity.phys.uwm.edu",
	description = "LSC Graphics Toolkit",
	url = "http://www.lsc-group.phys.uwm.edu/daswg/",
	license = "See file LICENSE",
	packages = ["lgen"],
	ext_modules = [
		Extension("lgen.support", ["src/support.c"],
			include_dirs = lalincdirs,
			libraries = lallibs,
			library_dirs = lallibdirs)
	],
	scripts = [
		os.path.join("bin", "plotburst"),
		os.path.join("bin", "plotsiminspiral"),
		os.path.join("bin", "plotgbb"),
		os.path.join("bin", "plotinspiral"),
		os.path.join("bin", "plotinspinj")
	],
	data_files = [
		("etc", [
			os.path.join("etc", "lgen-user-env.sh"),
			os.path.join("etc", "lgen-user-env.csh")
		])
	]
)
