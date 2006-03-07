# $Id$
# 
# setup for pylal

import os
from distutils.core import setup, Extension

def stripfirsttwo(string):
	return string[2:]

lallibs = map(stripfirsttwo, os.popen("pkg-config lal lalframe lalmetaio lalsupport --libs-only-l").read().split())

lallibdirs = map(stripfirsttwo, os.popen("pkg-config lal --libs-only-L").read().split())
lallibdirs = lallibdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --libs-only-L").read().split())
lallibdirs = lallibdirs + map(stripfirsttwo, os.popen("pkg-config libframe --libs-only-L").read().split())

lalincdirs = map(stripfirsttwo, os.popen("pkg-config lal --cflags-only-I").read().split())
lalincdirs = lalincdirs + map(stripfirsttwo, os.popen("pkg-config libmetaio --cflags-only-I").read().split())
lalincdirs = lalincdirs + map(stripfirsttwo, os.popen("pkg-config libframe --cflags-only-I").read().split())


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
			runtime_library_dirs = lallibdirs),
		Extension("pylal.frgetvect", ["src/frgetvect.c"],
			include_dirs = lalincdirs,
			libraries = lallibs,
			library_dirs = lallibdirs,
			runtime_library_dirs = lallibdirs)
	],
	scripts = [
		os.path.join("bin", "plotbinj"),
		os.path.join("bin", "plotburst"),
		os.path.join("bin", "plotburstrate"),
		os.path.join("bin", "plotchannel"),
		os.path.join("bin", "plotsiminspiral"),
		os.path.join("bin", "plotgbb"),
		os.path.join("bin", "plotinspiral"),
		os.path.join("bin", "plotinspinj"),
		os.path.join("bin", "plotinspdiff"),
		os.path.join("bin", "plotinspmissed"),
		os.path.join("bin", "plotnumtemplates"),
		os.path.join("bin", "plotinspiralrange"),
		os.path.join("bin", "plotcoincseglength"),
		os.path.join("bin", "plotthinca"),
		os.path.join("bin", "plotwindow"),
		os.path.join("bin", "plotcoincwindow"),
		os.path.join("bin", "ploteffdistcut"),
		os.path.join("bin", "plotefficiency"),
		os.path.join("bin", "plotsnrchi"),
		os.path.join("bin", "plotdistance"),
		os.path.join("bin", "s3_statistic"),
		os.path.join("bin", "lalapps_ll2cache"),
		os.path.join("bin", "lalapps_path2cache"),
		os.path.join("bin", "lalapps_pire"),
		os.path.join("bin", "ligolw_binjfind"),
		os.path.join("bin", "ligolw_bucluster"),
		os.path.join("bin", "ligolw_bucut"),
		os.path.join("bin", "ligolw_tisi")
	],
	data_files = [
		("etc", [
			os.path.join("etc", "pylal-user-env.sh"),
			os.path.join("etc", "pylal-user-env.csh")
		])
	]
)
