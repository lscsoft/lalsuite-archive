# $Id$
# 
# setup for pylal

import os
from distutils.core import setup, Extension

def stripfirsttwo(string):
	return string[2:]

lallibs = map(stripfirsttwo, os.popen("pkg-config lal lalmetaio lalsupport --libs-only-l").read().split())

lallibdirs = map(stripfirsttwo, os.popen("pkg-config lal --libs-only-L").read().split())
lallibdirs = lallibdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --libs-only-L").read().split())

lalincdirs = map(stripfirsttwo, os.popen("pkg-config lal --cflags-only-I").read().split())
lalincdirs = lalincdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --cflags-only-I").read().split())


setup(
	name = "pylal",
	version = "0.1",
	author = "Patrick Brady",
	author_email = "patrick@gravity.phys.uwm.edu",
	description = "LSC Graphics Toolkit",
	url = "http://www.lsc-group.phys.uwm.edu/daswg/",
	license = "See file LICENSE",
	packages = ["pylal"],
	ext_modules = [
		Extension("pylal.support", ["src/support.c"],
			include_dirs = lalincdirs,
			libraries = lallibs,
			library_dirs = lallibdirs,
			runtime_library_dirs = lallibdirs)
	],
	scripts = [
		os.path.join("bin", "plotburst"),
		os.path.join("bin", "plotsiminspiral"),
		os.path.join("bin", "plotgbb"),
		os.path.join("bin", "plotinspiral"),
		os.path.join("bin", "plotinspinj"),
    os.path.join("bin", "plotinspdiff"),
    os.path.join("bin", "plotinspmissed"),
    os.path.join("bin", "plotthinca"),
    os.path.join("bin", "plotwindow"),
    os.path.join("bin", "plotcoincwindow"),
    os.path.join("bin", "lalapps_pire")

	],
	data_files = [
		("etc", [
			os.path.join("etc", "pylal-user-env.sh"),
			os.path.join("etc", "pylal-user-env.csh")
		])
	]
)
