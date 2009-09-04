# $Id$
# 
# setup for pylal


import os
from misc import determine_git_version
from distutils.core import setup, Extension
from distutils.command import install
from distutils.command import build_py
from distutils.command import sdist
from distutils import log
from sys import version_info
from numpy.lib.utils import get_include as numpy_get_include


class PkgConfig(object):
	def __init__(self, names):
		def stripfirsttwo(string):
			return string[2:]
		self.libs = map(stripfirsttwo, os.popen("pkg-config --libs-only-l %s" % names).read().split())
		self.libdirs = map(stripfirsttwo, os.popen("pkg-config --libs-only-L %s" % names).read().split())
		self.incdirs = map(stripfirsttwo, os.popen("pkg-config --cflags-only-I %s" % names).read().split())
		self.extra_cflags = os.popen("pkg-config --cflags-only-other %s" % names).read().split()

lal_pkg_config = PkgConfig("lal")
lalframe_pkg_config = PkgConfig("lalframe")

def remove_root(path, root):
	if root:
		return os.path.normpath(path).replace(os.path.normpath(root), "")
	return os.path.normpath(path)

class pylal_build_py(build_py.build_py):
	def run(self):
		# create the git_version module
		if determine_git_version.in_git_repository():
			try:
				log.info("generating pylal/git_version.py")
				git_version_fileobj = open("pylal/git_version.py", "w")
				determine_git_version.write_git_version(git_version_fileobj)
			finally:
				git_version_fileobj.close()
		elif os.path.exists("pylal/git_version.py"):
			# We're probably being built from a release tarball; don't overwrite
			log.info("not in git checkout; using existing pylal/git_version.py")
		else:
			log.info("not in git checkout; writing empty pylal/git_version.py")
			try:
				git_version_fileobj = open("pylal/git_version.py", "w")
				determine_git_version.write_empty_git_version(git_version_fileobj)
			finally:
				git_version_fileobj.close()

		# resume normal build procedure
		build_py.build_py.run(self)

class pylal_install(install.install):
	def run(self):
		# create the user env scripts
		if self.install_purelib == self.install_platlib:
			pylal_pythonpath = self.install_purelib
		else:
			pylal_pythonpath = self.install_platlib + ":" + self.install_purelib

		pylal_prefix = remove_root(self.prefix, self.root)
		pylal_install_scripts = remove_root(self.install_scripts, self.root)
		pylal_pythonpath = remove_root(pylal_pythonpath, self.root)
		pylal_install_platlib = remove_root(self.install_platlib, self.root)

		if not os.path.isdir("etc"):
			os.mkdir("etc")
		log.info("creating pylal-user-env.sh script")
		env_file = open(os.path.join("etc", "pylal-user-env.sh"), "w")
		print >> env_file, "# Source this file to access PYLAL"
		print >> env_file, "PYLAL_PREFIX=" + pylal_prefix
		print >> env_file, "export PYLAL_PREFIX"
		print >> env_file, "PATH=" + pylal_install_scripts + ":${PATH}"
		print >> env_file, "PYTHONPATH=" + pylal_pythonpath + ":${PYTHONPATH}"
		print >> env_file, "LD_LIBRARY_PATH=" + pylal_install_platlib + ":${LD_LIBRARY_PATH}"
		print >> env_file, "DYLD_LIBRARY_PATH=" + pylal_install_platlib + ":${DYLD_LIBRARY_PATH}"
		print >> env_file, "export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH"
		env_file.close()

		log.info("creating pylal-user-env.csh script")
		env_file = open(os.path.join("etc", "pylal-user-env.csh"), "w")
		print >> env_file, "# Source this file to access PYLAL"
		print >> env_file, "setenv PYLAL_PREFIX " + pylal_prefix
		print >> env_file, "setenv PATH " + pylal_install_scripts + ":${PATH}"
		print >> env_file, "if ( $?PYTHONPATH ) then"
		print >> env_file, "  setenv PYTHONPATH " + pylal_pythonpath + ":${PYTHONPATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv PYTHONPATH " + pylal_pythonpath
		print >> env_file, "endif"
		print >> env_file, "if ( $?LD_LIBRARY_PATH ) then"
		print >> env_file, "  setenv LD_LIBRARY_PATH " + pylal_install_platlib + ":${LD_LIBRARY_PATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv LD_LIBRARY_PATH " + pylal_install_platlib
		print >> env_file, "endif"
		print >> env_file, "if ( $?DYLD_LIBRARY_PATH ) then"
		print >> env_file, "  setenv DYLD_LIBRARY_PATH " + pylal_install_platlib + ":${DYLD_LIBRARY_PATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv DYLD_LIBRARY_PATH " + pylal_install_platlib
		print >> env_file, "endif"
		env_file.close()

		# now run the installer
		install.install.run(self)


class pylal_sdist(sdist.sdist):
	def run(self):
		# remove the automatically generated user env scripts
		for script in ["pylal-user-env.sh", "pylal-user-env.csh"]:
			log.info("removing " + script )
			try:
				os.unlink(os.path.join("etc", script))
			except:
				pass

		# create the git_version module
		if determine_git_version.in_git_repository():
			log.info("generating pylal/git_version.py")
			try:
				git_version_fileobj = open("pylal/git_version.py", "w")
				determine_git_version.write_git_version(git_version_fileobj)
			finally:
				git_version_fileobj.close()
		else:
			log.info("not in git checkout; writing empty pylal/git_version.py")
			try:
				git_version_fileobj = open("pylal/git_version.py", "w")
				determine_git_version.write_empty_git_version(git_version_fileobj)
			finally:
				git_version_fileobj.close()

		# now run sdist
		sdist.sdist.run(self)


setup(
	name = "pylal",
	version = "0.1",
	author = "Patrick Brady",
	author_email = "patrick@gravity.phys.uwm.edu",
	description = "LSC Graphics Toolkit",
	url = "http://www.lsc-group.phys.uwm.edu/daswg/",
	license = "See file LICENSE",
	packages = [
		"pylal",
		"pylal.xlal"
	],
	cmdclass = {
		"build_py": pylal_build_py,
		"install": pylal_install,
		"sdist": pylal_sdist
	},
	ext_modules = [
		Extension(
			"pylal.Fr",
			["src/Fr.c"],
			include_dirs = lalframe_pkg_config.incdirs + [numpy_get_include()],
			libraries = lalframe_pkg_config.libs,
			library_dirs = lalframe_pkg_config.libdirs,
			runtime_library_dirs = lalframe_pkg_config.libdirs,
			extra_compile_args = lalframe_pkg_config.extra_cflags
		),
		Extension(
			"pylal.tools",
			["src/tools.c"],
			include_dirs = lal_pkg_config.incdirs,
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		),
		Extension(
			"pylal.xlal.date",
			["src/xlal/date.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include()],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		),
		Extension(
			"pylal.xlal.inject",
			["src/xlal/inject.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		),
		Extension(
			"pylal.xlal.tools",
			["src/xlal/tools.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		),
		Extension(
			"pylal.xlal.window",
			["src/xlal/window.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include()],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		),
		Extension(
			"pylal.xlal.burstsearch",
			["src/xlal/burstsearch.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include()],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs
		)
	],
	scripts = [
		os.path.join("bin", "analyseQscan.py"),
		os.path.join("bin", "distrib_fu_qscan_results.py"),
		os.path.join("bin", "submit_remote_scan.py"),
		os.path.join("bin", "exttrig_likelihood_pipe"),
		os.path.join("bin", "fup_triggers.py"),
		os.path.join("bin", "grbSelect"),
		os.path.join("bin", "galaxies_in_polygon"),
		os.path.join("bin", "lal_query_cache"),
		os.path.join("bin", "makeCheckList.py"),
		os.path.join("bin", "paste_insp_triggers"),
		os.path.join("bin", "plotbank"),
		os.path.join("bin", "plotbinj"),
		os.path.join("bin", "plotburca"),
		os.path.join("bin", "plotburca2"),
		os.path.join("bin", "plotburst"),
		os.path.join("bin", "plotburstrate"),
		os.path.join("bin", "plotchannel"),
		os.path.join("bin", "plotcohsnr"),
		os.path.join("bin", "plotcoincmissed"),
		os.path.join("bin", "plotchiatimeseries"),
		os.path.join("bin", "plotdetresponse"),
		os.path.join("bin", "plotgrbl"),
		os.path.join("bin", "plotlalseries"),
		os.path.join("bin", "plotsiminspiral"),
		os.path.join("bin", "plotnumgalaxies"),
		os.path.join("bin", "calcMassCut"),
		os.path.join("bin", "upperlimit.py"),
		os.path.join("bin", "write_iul_page"),
		os.path.join("bin", "lalapps_compute_posterior"),
		os.path.join("bin", "plotulvsmass"), 
		os.path.join("bin", "plotifar"),
		os.path.join("bin", "plotinjnum"),
		os.path.join("bin", "plotinspfound"),
		os.path.join("bin", "plotinspiral"),
		os.path.join("bin", "plotinspinj"),
		os.path.join("bin", "plotinspdiff"),
		os.path.join("bin", "plotinspmissed"),
		os.path.join("bin", "plotnumtemplates"),
		os.path.join("bin", "plotinspiralrange"),
		os.path.join("bin", "plotcoincseglength"),
		os.path.join("bin", "plotsegments"),
		os.path.join("bin", "plotthinca"),
		os.path.join("bin", "pylal_cache_to_mvsc.py"),
		os.path.join("bin", "pylal_mvsc_player.py"),
		os.path.join("bin", "mvsc_plots.py"),
		os.path.join("bin", "mvsc_plot_cuts.py"),
		os.path.join("bin", "mvsc_htmlwriter.py"),
		os.path.join("bin", "pylal_combine_posteriors"),
		os.path.join("bin", "pylal_followup_missed"),
		os.path.join("bin", "followupRescueHtml"),
		os.path.join("bin", "pylal_exttrig_summary"),
		os.path.join("bin", "pylal_grblikelihood"),
		os.path.join("bin", "pylal_grbUL"),
		os.path.join("bin", "pylal_grbtimeslide_stats"),
		os.path.join("bin", "pylal_exttrig_llmonitor"),
		os.path.join("bin", "pylal_exttrig_llsummary"),
		os.path.join("bin", "pylal_query_dq"),
		os.path.join("bin", "pylal_relic"),
		os.path.join("bin", "plotethinca"),
		os.path.join("bin", "plotwindow"),
		os.path.join("bin", "plotcoincwindow"),
		os.path.join("bin", "ploteffdistcut"),
		os.path.join("bin", "plotefficiency"),
		os.path.join("bin", "plotsnrchi"),
		os.path.join("bin", "frame_check"),
		os.path.join("bin", "IFOstatus_check"),
		os.path.join("bin", "plotsnrchisq_pipe"),
		os.path.join("bin", "plotmcmc.py"),
		os.path.join("bin", "plotinsppop"),
		os.path.join("bin", "plottisi"),
		os.path.join("bin", "query_dagman_log"),
		os.path.join("bin", "s3_statistic"),
		os.path.join("bin", "antime"),
		os.path.join("bin", "septime"),
		os.path.join("bin", "hipecoire"),
		os.path.join("bin", "inspiral_likelihood"),
		os.path.join("bin", "lalapps_cbc_plotsummary"),
		os.path.join("bin", "lalapps_excesspowerfinal"),
		os.path.join("bin", "lalapps_ll2cache"),
		os.path.join("bin", "lalapps_newcorse"),
		os.path.join("bin", "lalapps_path2cache"),
		os.path.join("bin", "lalapps_stringfinal"),
		os.path.join("bin", "ligolw_binjfind"),
		os.path.join("bin", "ligolw_bucluster"),
		os.path.join("bin", "ligolw_bucut"),
		os.path.join("bin", "ligolw_burca"),
		os.path.join("bin", "ligolw_burca_tailor"),
		os.path.join("bin", "ligolw_cafe"),
		os.path.join("bin", "ligolw_conv_inspid"),
		os.path.join("bin", "ligolw_inspinjfind"),
		os.path.join("bin", "ligolw_rinca"),
		os.path.join("bin", "ligolw_segments"),
		os.path.join("bin", "ligolw_sschunk"),
		os.path.join("bin", "ligolw_sicluster"),
		os.path.join("bin", "ligolw_tisi"),
		os.path.join("bin", "ligolw_thinca"),
		os.path.join("bin", "ligolw_thinca_to_coinc"),
		os.path.join("bin", "ligolw_veto"),
		os.path.join("bin", "ligolw_cbc_hardware_inj_page"),
		os.path.join("bin", "ligolw_fr_to_science"),
		os.path.join("bin", "inspiral_likelihood"),
		os.path.join("bin", "inspiral_likelihood_hipe"),
		os.path.join("bin", "KW_veto_setup"),
		os.path.join("bin", "KW_veto_calc"),
		os.path.join("bin", "KW_veto_plots"),
		os.path.join("bin", "KW_veto_channelPage"),
		os.path.join("bin", "KW_veto_reportPage"),
		os.path.join("bin", "KW_veto_qscanSetup"),
		os.path.join("bin", "pylal_plot_inspiral_skymap"),
		os.path.join("bin", "plotskypoints"),
		os.path.join("bin", "upper_limit_results"),
		os.path.join("bin", "pylal_expose"),
		os.path.join("bin", "ligolw_cbc_dbsimplify"),
		os.path.join("bin", "ligolw_cbc_dbaddinj"),
		os.path.join("bin", "ligolw_cbc_printlc"),
		os.path.join("bin", "ligolw_cbc_cluster_coincs"),
		os.path.join("bin", "ligolw_cbc_cfar"),
		os.path.join("bin", "ligolw_cbc_plotslides"),
		os.path.join("bin", "ligolw_cbc_plotifar"),
		os.path.join("bin", "ligolw_cbc_compute_durations"),
		os.path.join("bin", "extractCommand"),
		os.path.join("bin", "OddsPostProc.py"),
		os.path.join("bin", "prepare_sendback.py"),
		os.path.join("bin", "qsub_wscan.sh"),
		os.path.join("bin", "qsub_wscanlite.sh"),
		os.path.join("bin", "virgo_qscan_in2p3.py"),
		os.path.join("bin", "wscan_in2p3.sh"),
		os.path.join("bin", "wscanlite_in2p3.sh")
	],
	data_files = [ ("etc", [
		os.path.join("etc", "pylal-user-env.sh"),
		os.path.join("etc", "pylal-user-env.csh")
		] ) ]
)
