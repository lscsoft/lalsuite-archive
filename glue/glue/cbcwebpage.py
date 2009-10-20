# Copyright (C) 2009 Chad Hanna
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

from glue import markup
from glue.markup import oneliner as e
from glue import git_version

import subprocess
import os, sys, time, socket
import shutil

__author__ = "Chad Hanna <channa@caltech.edu>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


###############################################################################
##### UTILITY FUNCTIONS #######################################################
###############################################################################
def web_path_to_url(path):
	host = socket.getfqdn()
	pl = path.rstrip('/').split('/')

	#FIXME add more hosts as you need them
	if 'ligo.caltech.edu' in host: return "https://ldas-jobs.ligo.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'ligo-la.caltech.edu' in host: return "https://ldas-jobs.ligo-la.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'ligo-wa.caltech.edu' in host: return "https://ldas-jobs.ligo-wa.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'phys.uwm.edu' in host: return "https://ldas-jobs.phys.uwm.edu/~" + pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'phy.syr.edu' in host: return "https://sugar-jobs.phy.syr.edu/~" + pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'aei.uni-hannover.de' in host: return "https://atlas.atlas.aei.uni-hannover.de/~" + pl[pl.index('WWW')-1] + '/' + '/'.join(pl[pl.index('WWW')+1:])
	print sys.stderr, "WARNING: could not find web server, returning empty string"
	return ''


def create_toggle(filename="toggle.js"):
	"""
This function is just an alias to create a javascript for the toggle on/off. 
  
@return: nothing
	"""
	fname = open(filename, "w")
	fname.write("""
function toggle2(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
    		ele.style.display = "none";
  	}
	else {
		ele.style.display = "block";
	}
}
function afterLoadFrame() {
	$('#iframecontent a[rel="external"]').attr('target','_blank');
	$('#iframecontent input').hide();
	$('#iframecontent p:first').hide(); 
	}
function loadFrame(sourceURL) {
	$("#iframecontent").load(sourceURL,{},afterLoadFrame);
	/* Remove the last two arguments to disable toggling from the title. */
	}
	""")
	fname.close()
	return filename

def script_dict():
	script = {}
	tog = create_toggle()
	script[tog] = 'javascript'
	script['http://ajax.googleapis.com/ajax/libs/jquery/1.2.6/jquery.min.js'] = 'javascript'
	return script


def copy_ihope_style(stylefile="cbcwebpage.css", base_dir="."):

	# FIXME this is a stupid way to find the path... changes to build scripts, set env var?
	path = which('ligo_data_find')
	if path: path = os.path.split(path)[0]
	else: 
		print >>sys.stderr, "COULD NOT FIND STYLE FILES %s IN %s, ABORTING" % (stylefile, path)
		raise ValueError
		sys.exit(1)
	out = path.replace('bin','etc') + '/' + stylefile
	if not os.path.isfile(out):
		print >>sys.stderr, "COULD NOT FIND STYLE FILES %s IN %s, ABORTING" % (stylefile, path)
		raise ValueError
		sys.exit(1)
	shutil.copy(out, base_dir)
	
	return base_dir + '/' + os.path.split(out.rstrip('/'))[1]

def which(prog):
	which = subprocess.Popen(['which',prog], stdout=subprocess.PIPE)
	out = which.stdout.read().strip()
	if not out:
		print >>sys.stderr, "ERROR: could not find %s in your path, have you built the proper software and source the proper env. scripts?" % (prog,prog)
		raise ValueError
		sys.exit(1)
	return out

def user_and_date():
	tmstr = "/".join([str(i) for i in time.gmtime()[0:3]])
        tmstr += " " + ":".join([str(i) for i in time.gmtime()[3:5]])
	return "%s - %s" % (os.environ['USER'], tmstr)

###############################################################################
##### CBC WEB PAGE CLASSES ####################################################
###############################################################################

# PROBABLY DOES NOT EVER NEED TO BE USED DIRECTLY, BUT IS USED IN cbcpage
class _subpage_id(object):
	def __init__(self, id, link_text, closed_flag=0):
		self.id = id
		self.link_text = link_text
		self.closed_flag = closed_flag

class _imagelink(markup.page):
	def __init__(self, imageurl, thumburl, tag="img", width=240):
		markup.page.__init__(self, mode="strict_html")
		self.add("<a href=%s><img src=%s width=%d></a>" % (imageurl, thumburl, width))

	def get_content(self):
		return self.content
		
class _imagelinkcpy(markup.page):
	def __init__(self, imagepath, thumbpath=None, tag="img", width=240):
		markup.page.__init__(self, mode="strict_html")
		try: os.mkdir('Images')
		except: pass
		#So that you can give it a url
		imagepath.replace('file://localhost','').strip()
		imgname = os.path.split(imagepath.rstrip('/'))[1]
		shutil.copy(imagepath, 'Images/')
		if not thumbpath:
			thumbname = 'Images/' + "thumb_" + imgname
			command = 'convert ' + 'Images/' + imgname + ' -resize 300x300 -antialias ' + thumbname
			popen = subprocess.Popen(command.split())
			popen.communicate()
			status = popen.returncode
			imgname = 'Images/' + imgname
		else:
			thumbpath.replace('file://localhost','').strip()
			thumbname = os.path.split(thumbpath.rstrip('/'))[1]
			shutil.copy(thumbpath, 'Images/')
			thumbname = 'Images/' + thumbname
			imgname = 'Images/' + imgname
		self.add("<a href=%s><img src=%s width=%d></a>" % (imgname, thumbname, width))

	def get_content(self):
		return self.content
	

class _table(markup.page):
	def __init__(self, two_d_data, title="", caption="", tag="table", num="1"):
		markup.page.__init__(self, mode="strict_html")
		self.add("<br>")
		if title: 
			self.b("%s. %s" %(num, title.upper()) )
	
		self.table()
		for row in two_d_data:
			self.tr
			tdstr = ""
			for col in row:
				tdstr += "<td>%s</td>" % (str(col),)
			self.add(tdstr)
			self.tr.close()
		self.table.close()
		if self.caption: self.i("%s. %s" %(num, caption))
		self.add("<br>")

	def get_content(self):
		return self.content			

# PROBABLY DOES NOT EVER NEED TO BE USED DIRECTLY, BUT IS USED IN cbcpage
class _section(markup.page):
	def __init__(self, tag, title="", secnum="1", level=2):
		markup.page.__init__(self, mode="strict_html")
		self.secnum = secnum
		self._title = title
		self.sections = {}
		self.section_ids = []
		self.level = level 
		self.tag = tag
		self.id = tag + self.secnum
		self.tables = 0
		self.add('<div class="contenu"><h%d id="toggle_%s" onclick="javascript:toggle2(\'div_%s\', \'toggle_%s\');"> %s. %s </h%d>' % (level, self.id, secnum, self.id, secnum, title, level) )
		self.div(id="div_"+secnum , style='display:none;')

	def add_section(self, tag, title=""):
		secnum = "%s.%d" % (self.secnum, len(self.sections.values())+1)
		self.sections[tag] = _section(tag, title=title, secnum=secnum, level=self.level+1)
		self.section_ids.append([len(self.sections.values()), tag])

	def get_content(self):
		self.section_ids.sort()
		out = self.content
		self.div.close()
		self.div.close()
		for num, key in self.section_ids:
			out.extend(self.sections[key].get_content())
		return out
	
	def add_table(self, two_d_data, title="", caption="", tag="table", num=0):
		self.tables += 1
		tabnum = self.secnum+"T:"+str(self.tables)
		table = _table(two_d_data, title=title, caption=caption, tag="table", num=tabnum)
		self.content.extend(table.get_content())



# MAIN CBC WEB PAGE CLASS 
class cbcpage(markup.page):

	def __init__(self, title="cbc web page", path='./', css=None, script=None, verbose=False):
		"""
		"""
		if not css: css = copy_ihope_style()
		if not script: script = script_dict()
		self.verbose = verbose
		self._style = css
		self._title = title
		self._script = script
		self.path = path

		markup.page.__init__(self, mode="strict_html")	
		self._escape = False
		doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
		doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""
		self.init(title=self._title, css=self._style, script=self._script , doctype=doctype)
		
		self.subpages = {}
		self.subpage_ids = [] 
		self.external_frames = []	
		
		self.sections = {}
		self.section_ids = []

		self.tables = 0
		
	def add_subpage(self, tag, title, link_text=None):
		""" 
		"""

		subpage_num = len(self.subpages.values())
		if not link_text: link_text=str(subpage_num)

		# tuple including number so it can be sorted later
		self.subpages[tag] = cbcpage(title=title,css=self._style,script=self._script)
		self.subpage_ids.append( [subpage_num, _subpage_id(tag, link_text)] )

	def close_subpage(self,id=None):
		#SECTIONS WILL AUTOMATICALLY BE CLOSED IN WRITE METHOD IF NOT DONE EXPLICITELY
		self.subpage_ids.sort()
		if not id: id = subpage_ids[-1][1].id

		self.subpages[id].div.close()
		self.subpages[id].add("<!-- close div contenu-->")
		self.subpages[id].div.close()
		self.subpages[id].add_footer()

	def add_footer(self):
		#FIXME DO SOMETHING
		pass

	def add_external_frame(self, linkurl, linktext):
		self.external_frames.append([linkurl, linktext])

	def write(self, file_name="index"):

		if self.subpage_ids:

			self.div(id_="wrapper")
			self.div(id_="menubar")
			self.div(id_="menu")
			self.subpage_ids.sort()
			# do the sub pages
			for num,secid in self.subpage_ids:
				id = secid.id
				if secid.closed_flag == 0: self.close_subpage(id)
				secfname = file_name + "_" + id
				self.subpages[id].write(secfname)
				self.div(class_="menuitem")
				self.add('\t<a class="menulink" href="javascript:loadFrame(\'%s.html\');"> %d: %s </a>\n' % (secfname, num, secid.link_text) )
				self.div.close()
			for i, ext_frame in enumerate(self.external_frames):
				self.div(class_="menuitem")
				self.add('\t<a class="menulink" href="javascript:loadFrame(\'%s\');"> %d: %s </a>\n' % (ext_frame[0], num+i, ext_frame[1]) )
				self.div.close()
			self.div.close()
			self.div.close()
			self.div(id_="ihope")
			self.add('<h2> CBC </h2>')
			self.add('<img width=90 src="http://upload.wikimedia.org/wikipedia/commons/thumb/5/5e/BH_LMC.png/180px-BH_LMC.png">')
			self.div.close()
			self.div(id_='header')
			self.add('<h1>' + self._title  +' </h1>')
			self.add('<h3> ' + user_and_date() + ' </h3>')
			self.div.close()
			self.div(id_='iframecontent')
			self.add('<p id="placeholder">Please select a report section on the left.</p>')
			self.div.close()
			self.div.close()

		# do the sections
		self.section_ids.sort()
		for num, key in self.section_ids:
			self.content.extend(self.sections[key].get_content())
		pagefile = file('%s/%s.html' % (self.path, file_name), 'w')
		pagefile.write(str(self))
		pagefile.close()
			
	def add_section(self, tag, title="", level=2):
		"""
		"""
		secnum = len(self.sections.values()) + 1
		self.section_ids.append([secnum, tag])
		self.sections[tag] = _section(title=title, tag=tag, secnum=str(secnum), level=level)

	def add_table(self, two_d_data, title="", caption="", tag="table"):
		self.tables += 1
		table = _table(two_d_data, title=title, caption=caption, tag="table", num="T:"+str(self.tables))
		self.content.extend(table.get_content())
