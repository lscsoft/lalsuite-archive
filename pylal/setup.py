# $Id$
# 
# setup for pylal

import os
from distutils.core import setup, Extension

class PkgConfig(object):
	def __init__(self, names):
		def stripfirsttwo(string):
			return string[2:]
		self.libs = map(stripfirsttwo, os.popen("pkg-config --libs-only-l %s" % names).read().split())
		self.libdirs = map(stripfirsttwo, os.popen("pkg-config --libs-only-L %s" % names).read().split())
		self.incdirs = map(stripfirsttwo, os.popen("pkg-config --cflags-only-I %s" % names).read().split())

full_lal_pkg_config = PkgConfig("lal lalframe lalmetaio lalsupport")
lal_pkg_config = PkgConfig("lal")


setup(
	name = "pylal",
	version = "0.1",
	author = "Patrick Brady",
	author_email = "patrick@gravity.phys.uwm.edu",
	description = "LSC Graphics Toolkit",
	url = "http://www.lsc-group.phys.uwm.edu/daswg/",
	license = "See file LICENSE",
	packages = ["pylal", "pylal.xlal"],
	ext_modules = [
		Extension("pylal.support", ["src/support.c"],
			include_dirs = full_lal_pkg_config.incdirs,
			libraries = full_lal_pkg_config.libs,
			library_dirs = full_lal_pkg_config.libdirs,
			runtime_library_dirs = full_lal_pkg_config.libdirs),
		Extension("pylal.frgetvect", ["src/frgetvect.c"],
			include_dirs = full_lal_pkg_config.incdirs,
			libraries = full_lal_pkg_config.libs,
			library_dirs = full_lal_pkg_config.libdirs,
			runtime_library_dirs = full_lal_pkg_config.libdirs),
		Extension("pylal.xlal.date", ["src/xlal/date.c"],
			include_dirs = lal_pkg_config.incdirs,
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs)
	],
	scripts = [
		os.path.join("bin", "plotbinj"),
		os.path.join("bin", "plotburca"),
		os.path.join("bin", "plotburst"),
		os.path.join("bin", "plotburstrate"),
		os.path.join("bin", "plotchannel"),
		os.path.join("bin", "plotsiminspiral"),
		os.path.join("bin", "plotnumgalaxies"),
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
		os.path.join("bin", "plotinsppop"),
		os.path.join("bin", "plottisi"),
		os.path.join("bin", "s3_statistic"),
		os.path.join("bin", "hipecoire"),
		os.path.join("bin", "lalapps_ll2cache"),
		os.path.join("bin", "lalapps_path2cache"),
		os.path.join("bin", "lalapps_pire"),
		os.path.join("bin", "ligolw_binjfind"),
		os.path.join("bin", "ligolw_bucluster"),
		os.path.join("bin", "ligolw_bucut"),
		os.path.join("bin", "ligolw_burca"),
		os.path.join("bin", "ligolw_cafe"),
		os.path.join("bin", "ligolw_segments"),
		os.path.join("bin", "ligolw_sschunk"),
		os.path.join("bin", "ligolw_tisi")
	],
	data_files = [
		("etc", [
			os.path.join("etc", "pylal-user-env.sh"),
			os.path.join("etc", "pylal-user-env.csh")
		])
	]
)
