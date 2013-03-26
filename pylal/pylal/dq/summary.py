#!/usr/bin/env python

# Copyright (C) 2012 Duncan M. Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

"""
This module defines the classes used to generate interferometer summary screens
from data.
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division

import re
import calendar
import datetime
import os
import sys
import numpy
import copy
import StringIO
import lal

from dateutil.relativedelta import relativedelta

# set display
from matplotlib import use
use("Agg")

from pylal import git_version
from pylal import htmlutils
from pylal import plotdata
from pylal import plotsegments
from pylal import plottriggers
from pylal import seriesutils
from pylal.dq import dqTriggerUtils
from pylal.plotutils import display_name as latex

from glue import markup
from glue import segments
from glue import segmentsUtils
from glue.ligolw import table

# set metadata
__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# set regex
_r_cchar  = re.compile("[\W_]+")
rindex = re.compile("index.html\Z")

# set mode enum
SUMMARY_MODE_GPS   = 0
SUMMARY_MODE_DAY   = 1
SUMMARY_MODE_WEEK  = 2
SUMMARY_MODE_MONTH = 3
SUMMARY_MODE_YEAR  = 5
MODE_NAME = {SUMMARY_MODE_GPS: "GPS",\
             SUMMARY_MODE_DAY: "DAY",\
             SUMMARY_MODE_WEEK: "WEEK",\
             SUMMARY_MODE_MONTH: "MONTH",\
             SUMMARY_MODE_YEAR: "YEAR"}

# =============================================================================
# Classes for tabs
# =============================================================================

class SummaryTab(object):
    """
    The default meta-class for the summary tabs. It provides the basics so
    developers should subclass SummaryTab to generate something useful.
    """
    mode = SUMMARY_MODE_GPS
    ifo  = None
    def __init__(self, name, **kwargs):
        """
        Initialise a SummaryTab object. At least a name must be given.

        @param name: a name string for the tab
        @type name: C{str}
        """ 
        # set name
        self.name = name

        # set defaults
        kwargs.setdefault("plots", list())
        kwargs.setdefault("subplots", list())
        kwargs.setdefault("children", list())
        kwargs.setdefault("parent", None)
        kwargs.setdefault("state", None)
        kwargs.setdefault("information", None)

        # set all other values
        for key,val in kwargs.iteritems():
            setattr(self, key, val) 
 
    # ================
    # basic properties

    @property
    def name(self):
        """Descriptive string for this Tab."""
        if not hasattr(self, "_name"):
            self.name = ""
        return self._name
    @name.setter
    def name(self, value):
        self._name = str(value)
    @name.deleter
    def name(self):
        del self._name

    @property
    def directory(self):
        """path to the directory for this Tab."""
        if not hasattr(self, "_directory"):
            n = _r_cchar.sub("_", self.name.lower())
            if hasattr(self, "parent") and self.parent is not None:
                self.directory = os.path.join(self.parent.directory, n)
            else:
                self.directory = os.path.join(os.path.curdir, n)
        return self._directory
    @directory.setter
    def directory(self, d):
        self._directory = os.path.normpath(d)
    @directory.deleter
    def directory(self):
        del self._directory

    @property
    def start_time(self):
        """GPS start time for this Tab."""
        return self._start_time
    @start_time.setter
    def start_time(self, gps):
        self._start_time = lal.LIGOTimeGPS(gps)
    @start_time.deleter
    def start_time(self):
        del self._start_time

    @property
    def end_time(self):
        """GPS end time for this Tab."""
        return self._end_time
    @end_time.setter
    def end_time(self, gps):
        self._end_time = lal.LIGOTimeGPS(gps)
    @end_time.deleter
    def end_time(self):
        del self._end_time

    @property
    def span(self):
        """GPS [start_time, stop_time) segment for this Tab."""
        return segments.segment(self.start_time, self.end_time)
    @span.setter
    def span(self, seg):
        self.start_time = lal.LIGOTimeGPS(seg[0])
        self.end_time = lal.LIGOTimeGPS(seg[1])
    @span.deleter
    def span(self):
        del self._span

    @property
    def segments(self):
        """glue.segments.segmentlist describing the analysis segments for
        this Tab.
        """
        return self._segments
    @segments.setter
    def segments(self, seglist):
        self._segments =\
            segments.segmentlist([segments.segment(map(float, s))\
                                  for s in seglist])
    @segments.deleter
    def segments(self):
        del self._segments

    # ===============
    # plot properties

    @property
    def plots(self):
        """list of plots for this Tab."""
        return self._plotlist
    @plots.setter
    def plots(self, plotlist):
        self._plotlist = list(map(os.path.normpath, plotlist))
    @plots.deleter
    def plots(self):
        del self._plotlist

    @property
    def subplots(self):
        """list of subplots for this Tab."""
        return self._subplotlist
    @subplots.setter
    def subplots(self, subplotlist):
        self._subplotlist = list(map(os.path.normpath, subplotlist))
        return self._subplotlist
    @subplots.deleter
    def subplots(self):
        del self._subplotlist

    # ===============
    # HTML properties

    @property
    def href(self):
        """path str to the frame for this Tab."""
        if not hasattr(self, "_href"):
            basename = "%s.html" % _r_cchar.sub("_", self.state.name).lower()
            self.href = os.path.join(self.directory, basename)
        return self._href
    @href.setter
    def href(self, url):
        self._href = os.path.normpath(url)
    @href.deleter
    def href(self):
        del self._href

    @property
    def index(self):
        """URL to the index for this Tab."""
        return os.path.normpath("%s%s" % (self.directory, os.path.sep))

    #
    # class methods
    #

    def add_child(self, tab):
        self.children.append(tab)
        self.plots.append(tab.plots[0])

    def get_child(self, name):
        names = [c.name for c in self.children]
        try:
            idx = names.index(name)
        except ValueError,e:
            raise RunTimeError("Parent tab has no child named \"%s\"." % name)
        else:
            return self.children[idx]
        
    def write_menu(self, jobdir, startdate=None, weekday=calendar.MONDAY):
        """
        Write the HTML menu for this tab.

        @param jobdir: output directory of this job
        @type jobdir: C{str}
        @param startdate: date at which to start the menubar calendar
        @type startdate: C{datetime.datetime}
        @param weekday: day on which to start week in calendar (0=monday)
        @type weekday: C{int}
        """
        self.menu = markup.page()
        self.menu.div(id_="menubar")

        # work out which sections we want
        if self.mode == SUMMARY_MODE_DAY:
            sections = ["Today", "Previous day", "Next day", "Calendar"]
            marker   = "today"
            strf     = "%Y%m%d"
            delta    = datetime.timedelta(days=1)
        elif self.mode == SUMMARY_MODE_WEEK:
            sections = ["This week", "Previous week", "Next week",\
                        "Calendar"]
            marker   = "thisweek"
            strf     = "%Y%m%d"
            delta    = datetime.timedelta(days=7)
        elif self.mode == SUMMARY_MODE_MONTH:
            sections = ["This month", "Previous month", "Next month",\
                        "Calendar"]
            marker   = "thismonth"
            strf     = "%Y%m"
            delta    = relativedelta(months=1)
        elif self.mode == SUMMARY_MODE_YEAR:
            sections = ["This year", "Previous year", "Next year",\
                        "Calendar"]
            marker   = "thisyear"
            strf     = "%Y"
            delta    = relativedelta(years=1)
        else:
            sections = ["Full data"]
        sections.extend(["About", "Glossary"])

        # link different run states
        if self.name == "Online" and len(self.states):
            self.menu.div(class_="statetoggle")
            self.menu.ul(id_="id_statetoggle")
            for i,s in enumerate(self.states):
                s2 = _r_cchar.sub("-", s.name.lower())
                setcls = """$("#li_%s").attr("class","%s");$("#div_%s").%s()"""
                onclick = []
                for j,tag in enumerate(self.states):
                    tag    = _r_cchar.sub("-", tag.name.lower())
                    class_ = j==i and "open" or "closed"
                    show   = j==i and "show" or "hide"
                    onclick.append(setcls % (tag, class_, tag, show))
                class_ = i==0 and "open" or "closed"
                id_    = "li_%s" % s2
                href   = "#%s"   % s2
                self.menu.li(s.name, class_=class_, id_=id_, href=href,\
                             onclick=";".join(onclick))
            self.menu.ul.close()
            self.menu.div.close()
        else:
            run_states = [child.state for child in self.children\
                          if child.state is not None]
            if len(run_states):
                self.menu.div(class_="statetoggle")
                self.menu.ul(id_="statetoggle")
                for i,state in enumerate(run_states):
                    title   = "%s time" % state.name.title()
                    class_  = i==0 and "open" or "closed"
                    id_     = "li_%s" % _r_cchar.sub("-", state.name.lower())
                    onclick = "$(this).loadrunstate(\"%s\");"\
                              % self.children[i].href
                    self.menu.li(title, class_=class_, id_=id_,\
                                 onclick=onclick)
                self.menu.ul.close()
                self.menu.div.close()
        
        # write section links
        if self.mode != SUMMARY_MODE_GPS:
            cday = datetime.datetime(*lal.GPSToUTC(int(self.start_time))[:6])
            cstrf = cday.strftime(strf)
        for i,sec in enumerate(sections):
            if re.match("(Today|This )", sec):
                if re.search(cstrf, self.index):
                    link =\
                        re.sub("%s%s%s" % (os.path.sep, cstrf, os.path.sep),\
                               "%s%s%s" % (os.path.sep, marker, os.path.sep),\
                               self.index)
                else:
                    link = os.path.join(jobdir, os.pardir, marker)
            elif sec.startswith("Previous "):
                previous = (cday-delta).strftime(strf)
                if re.search(cstrf, self.index):
                    link = self.index
                else:
                    link = jobdir
                link = link.replace(cstrf, previous)
            elif sec.startswith("Next "):
                next_ = (cday+delta).strftime(strf)
                if re.search(cstrf, self.index):
                    link = self.index
                else:
                    link = jobdir
                link = link.replace(cstrf, next_)
            elif sec == "About" and self.mode != SUMMARY_MODE_GPS:
                link = os.path.join(jobdir, sec.lower())
            else:
                link = sec.lower()
            link = "%s%s"\
                   % (os.path.normpath(rindex.sub("", link)), os.path.sep)
            self.menu.a(sec, id_="link_%d" % i, class_="menulink", href=link)

        # write menucalendar
        if self.mode == SUMMARY_MODE_GPS or self.name == "Calendar":
            self.menu.div.close()
        else:
            now = int(lal.GPSTimeNow())
            enddate = datetime.datetime(*lal.GPSToUTC(now)[:6])
            if self.name in ["Summary", "Glossary", "Calendar"]:
                menucal = calendar_page(startdate, enddate, weekday=weekday,\
                                        path=None, ncol=1, reverse=True)
                menupath = os.path.join("html", "menucal.html")
            else:
                menucal = calendar_page(startdate, enddate, weekday=weekday,\
                                        path=self.directory, ncol=1,\
                                        reverse=True)
                menupath = os.path.join("html", "%s_menucal.html"\
                            % _r_cchar.sub("_", re.split("\d\d\d\d%s"\
                                                         % os.path.sep,\
                                                         self.directory)[-1]))
            with open(menupath, "w") as menuf:
                menuf.write(re.sub(" class=\"\"", "", str(menucal)))
            self.menu.div("",id_="menucalendar")
            self.menu.script("$(\"div#menucalendar\").load(\"%s\");"\
                             % menupath, type="text/javascript")
            self.menu.div.close()
 
    def write_ifo_links(self, baselist, thisifo, tablink):
        """
        Write the links to summary pages for other interferometers.

        @param baselist: list of (ifo, url) pairs defining the base for each
            IFO.
        @type baselist: list of tuples
        @param thisifo: IFO prefix for this Tab.
        @type thisifo: C{str}
        @param tablink: relative URL of this Tab to append to all IFO base
            URLs.
        @type tablink: C{str}
        """
        if re.search("None", str(tablink)):
            tablink = ""
        else:
            tablink = rindex.sub("", tablink)
        self.ifobar = markup.page()
        for ifo,base in baselist:
            base = os.path.normpath(base)
            if ifo.upper() == thisifo.upper():
                self.ifobar.h1(markup.oneliner.a(ifo, href=tablink))
            else:
                self.ifobar.h3(markup.oneliner.a(ifo.upper(),\
                                          href=os.path.join(base, tablink)))

    def frametohtml(self):
        """
        Write this Tab's frame to file as HTML. This allows it to be jquery
        loaded into a parent page.
        """
        with open(self.href, "w") as framef:
            framef.write(re.sub(" class=\"\"", "", str(self.frame)))

    def build(self, **initargs):
        """
        Write this Tab's full HTML to the index file.

        @keyword **initargs: keyword arguments to pass to markup.page.init.
        """
        # default banner
        if isinstance(self, OnlineSummaryTab):
            banner = markup.oneliner.h1("Online data")
        else:
            date = datetime.date(*lal.GPSToUTC(int(self.start_time))[:3])
            if self.mode == SUMMARY_MODE_DAY:
                banner = "%s: %d-%d" % (date.strftime("%B %d %Y"),\
                                        self.start_time, self.end_time)
            elif self.mode == SUMMARY_MODE_WEEK:
                banner = "Week of %s: %d-%d" % (date.strftime("%B %d %Y"),\
                                                self.start_time, self.end_time)
            elif self.mode == SUMMARY_MODE_MONTH:
                banner = "%s: %d-%d" % (date.strftime("%B %Y"),\
                                        self.start_time, self.end_time)
            elif self.mode == SUMMARY_MODE_YEAR:
                banner = "%s: %d-%d" % (date.strftime("%Y"),\
                                        self.start_time, self.end_time)
            else:
                banner = "%d-%d" % (self.start_time, self.end_time)
            banner = markup.oneliner.h1(banner)

        # default online button
        if self.mode != SUMMARY_MODE_GPS:
            homebutton = markup.page()
            homebutton.ul(class_="buttons")
            homebutton.li(id_="online")
            homebutton.a("Online", href="")
            homebutton.li.close()
            homebutton.ul.close()
        else:
            homebutton = False

        # set init defaults
        initargs.setdefault("title", "Summary: %s" % self.name)

        page = htmlutils.build_page(icon=self.ifobar, banner=banner,\
                                    homebutton=homebutton, tabs=self.tabs,\
                                    menu=self.menu, frame=self.frame,\
                                    **initargs)
        with open(os.path.join(self.index, "index.html"), "w") as indexf:
            indexf.write(re.sub(" class=\"\"", "", str(page)))

    def add_information(cp, cpsec):
        """
        Add an "Information" section to this SummaryTab.

        @param cp: INI configuration for this job.
        @type cp: C{configparser.configparser}
        @param cpsec: INI-file section name from which to read the
            'html-information' option.
        @type cpsec: C{str}
        """
        if cp.has_option(cpsec, "html-information"):
            self.information = cp.get(cpsec, "html-information")
        else:
            sys.stderr.write("Warning. No HTML information found for %s. "\
                             "Please add this.\n" % self.name)
            sys.stderr.flush()


class SectionSummaryTab(SummaryTab):
    """
    SummaryTab representing an overview Tab for a set of child Tabs.
    """
    def write_tabs(self, sectiontabs):
        """
        Write the tabbar for this Tab.

        @param sectiontabs: list of SummaryTab objects to include in tabbar
        @type sectiontabs: C{list}
        """
        self.tabs = markup.page()

        # sort with "Summary" at the top, and "Misc." at the bottom
        sectiontabs.sort(key=lambda t: t.name == "Summary" and 1 or\
                                                 "Misc."   and 3 or 2)

        # build ordered list of all tabs
        alltabs = []
        for tab in sectiontabs:
            alltabs.append(tab)
            if tab.name != "Summary":
                alltabs.extend(tab.children)
        allnames = [t.name for t in alltabs]
        if self.name in allnames:
            current = allnames.index(self.name)
        else:
            current = None

        # draw left, right arrows
        self.tabs.div(class_="right")
        self.tabs.ul(class_="tabs", id_="arrows")
        if current is not None:
            self.tabs.li(title="Previous", class_="arrow")
            href = rindex.sub("", alltabs[(current-1)%len(alltabs)].index)
            self.tabs.a("&#x2190",href=href)
            self.tabs.li.close()
            self.tabs.li(title="Next", class_="arrow")
            href = rindex.sub("", alltabs[(current+1)%len(alltabs)].index)
            self.tabs.a("&#x2192",href=href)
            self.tabs.li.close()
        self.tabs.ul.close()
        self.tabs.div.close()

        # open tablist
        self.tabs.ul(class_="tabs")

        # write each SectionSummaryTab
        subtabs = None
        for i,tab in enumerate(sectiontabs):
            if self.name == tab.name\
            or (self.parent and self.parent.name == tab.name):
                class_  = "open"
                subtabs = tab.children
                parent  = tab
            else:
                class_  = "closed"
            self.tabs.li(id_=_r_cchar.sub("_", tab.name.lower()),\
                         class_=class_, title=self.name)
            self.tabs.a(tab.name, href=rindex.sub("", tab.index))
            if class_=="closed" and tab.name == "Summary":
                self.tabs.ul("",class_="dropdown")
            elif class_=="closed":
                self.tabs.ul(class_="dropdown")
                for subtab in tab.children:
                    self.tabs.li(id_=_r_cchar.sub("_", subtab.name.lower()),\
                            class_="closed", title="Open %s" % subtab.name)
                    self.tabs.a(subtab.name, href=rindex.sub("", subtab.index))
                    self.tabs.li.close()
                self.tabs.ul.close()
            self.tabs.li.close()
        self.tabs.ul.close()

        self.tabs.div(class_="clear")
        self.tabs.div.close()

        if self.name != "Summary" and subtabs:
            self.tabs.ul(class_="subtabs",\
                         id_=_r_cchar.sub("_", parent.name.lower())) 
            for tab in subtabs:
                if self.name == tab.name:
                    class_="open"
                else:
                    class_="closed"
                self.tabs.li(id_=_r_cchar.sub("_", tab.name.lower()),\
                             class_=class_, title="Open %s" % tab.name)
                self.tabs.a(tab.name, href=rindex.sub("", tab.index))
                self.tabs.li.close()
        else:
            self.tabs.ul(class_="subtabs")
        self.tabs.ul.close()

        # clear divs
        self.tabs.div("", class_="clear")

    def finalize(self):
        """
        Write the HTML frame for this SectionSummary.
        """
        self.frame = markup.page()

        self.frame.h1(self.name, id_="h1_0")
        self.frame.div(id_="div_0")
 
        if self.name == "Summary":
            order = ["Sensitivity", "Triggers", "Segments"]
            self.children.sort(key=lambda x: x.name in order\
                                             and order.index(x.name)+1 or 1000)
            children = []
            for tab in self.children:
                children.extend(tab.children)
        else:
            children = self.children
        n = len(children) > 1 and 2 or 1
        
        self.frame.table(style="table-layout: fixed; width: 100%;")
        for i,tab in enumerate(children):
            if i % n == 0:
                self.frame.tr()
            self.frame.td()
            self.frame.h3(tab.name, class_="summary")
            self.frame.a(href=tab.index, title=tab.name)
            self.frame.img(src=tab.plots[0][0], alt=tab.name, class_="full")
            self.frame.a.close()
            self.frame.td.close()
            if i % n == (n-1):
                self.frame.tr.close()
        self.frame.table.close()

        self.frame.div.close()

    def frametohtml(self):
        raise AttributeError("The frametohtml function should not be used"+\
                             " for SegmentSummaryTabs.")


class MetaStateTab(SectionSummaryTab):
    """
    This meta-class should stand between a SectionSummaryTab and each of it's
    children, allowing to hold the information for running the same tab in
    a variety of different run states.
    """
    def finalize(self):
        """
        Write the HTML page for this MetaStateTab.
        """
        self.frame = markup.page()
        state = self.children[0].state
        html = self.children[0].href
        self.frame.div.close()
        self.frame.script("$(\"#li_%s\").loadrunstate(\"%s\");"\
                          % (_r_cchar.sub("-", state.name).lower(), html),\
                          type_="text/javascript")
        self.frame.div()


class SegmentSummaryTab(SummaryTab):
    """
    Object representing a summary of segments.
    """
    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.flags = list()
        self.segdict = segments.segmentlistdict()
        self.segment_files = dict()

    # ==========
    # properties

    @property
    def segdict(self):
        """glue.segments.segmentlistdict containing glue.segments.segmentlists
        for each flag used in this SegmentSummaryTab."""
        return self._segdict
    @segdict.setter
    def segdict(self, segdict):
        # set type
        for key in segdict.keys():
            segdict[key] = type(segdict[key])([map(float) for s in\
                                               segdict[key]])
        self._segdict = segments.segmentlistdict(segdict)
    @segdict.deleter
    def segdict(self):
        del self._segdict

    #
    # class methods 
    #

    def add_segments(self, flag, seglist):
        """
        Add a glue.segments.segmentlist for a given flag to this object

        @param flag: name of flag to record
        @type flag: C{str}
        @param seglist: list of segments to record
        @type seglist: C{glue.segments.segmentlist}
       
        """
        if self.segdict.has_key(flag):
            self.segdict[flag] += segments.segmentlist(seglist)
        else:
            self.flags.append(flag)
            self.segdict[flag] = segments.segmentlist(seglist)

    def tosegwizard(self, flag, filename=None):
        """
        Write the segments know for the given flag to a segwizard format file.
        If an output filename is not given, one is generated.

        @param flag: name of flag whose segments are to be written
        @type flag: C{str}
        @param filename: path string of file to be written
        @type filename: C{str}
        """
        if not filename:
            name = _r_cchar.sub("-", _r_cchar.sub("_", flag.upper()), 1)
            filename = os.path.join(self.directory, "%s-%d-%d.txt"\
                       % (name, self.start_time, abs(self.span)))
        with open(filename, "w") as segf:
            segmentsUtils.tosegwizard(segf, self.segdict[flag])
        self.segment_files[flag] = filename

    def plotsegments(self, outfile, subplot=False, **kwargs):
        """
        Plot the segments associated with this SegmentSummaryTab.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotsegments.plotsegmentlistdict
        """
        if not subplot:
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        kwargs.setdefault("keys", list(self.flags))
        if re.search("segments\Z", self.name, re.I):
            kwargs.setdefault("title", self.name)
        else:
            kwargs.setdefault("title", "%s segments" % self.name)
        kwargs.setdefault("subtitle", "%d-%d" % (self.start_time,self.end_time))
        desc = kwargs.pop("description", None)
        plotsegments.plotsegmentlistdict(self.segdict, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotduration(self, outfile, flag, subplot=False, **kwargs):
        """
        Plot the segment durations for the given flag.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param flag: name of flag whose segments are to be plotted
        @type flag: C{str}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotsegments.plotsegmentlistdict
        """
        desc = kwargs.pop("description", None)
        if not subplot:
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        kwargs.setdefault("flag", flag)
        if re.search("segments\Z", self.name, re.I):
            kwargs.setdefault("title", self.name)
        else:
            kwargs.setdefault("title", "%s segments" % self.name)
        plotsegments.plotduration(self.segdict[flag], outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plothistogram(self, outfile, subplot=False, **kwargs):
        """
        Plot a histogram of the segments associated with this SegmentSummaryTab.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotsegments.plotsegmentlistdict
        """
        desc = kwargs.pop("description", None)
        kwargs.setdefault("keys", list(self.flags))
        if re.search("segments\Z", self.name, re.I):
            kwargs.setdefault("title", self.name)
        else:
            kwargs.setdefault("title", "%s segments" % self.name)
        plotsegments.plothistogram(self.segdict, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotdutycycle(self, outfile, subplot=False, **kwargs):
        """
        Plot the duty cycle of segments associated with this SegmentSummaryTab.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotsegments.plotsegmentlistdict
        """
        desc = kwargs.pop("description", None)
        kwargs.setdefault("keys", list(self.flags))
        if not subplot:
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        if re.search("segments\Z", self.name, re.I):
            kwargs.setdefault("title", self.name)
        else:
            kwargs.setdefault("title", "%s segments" % self.name)
        if self.mode == SUMMARY_MODE_DAY:
            kwargs.setdefault("binlength", 3600)
        elif self.mode == SUMMARY_MODE_WEEK or self.mode == SUMMARY_MODE_MONTH:
            kwargs.setdefault("binlength", 86400)
        elif self.mode == SUMMARY_MODE_YEAR:
            kwargs.setdefault("binlength", 365/12*86400)
        plotsegments.plotdutycycle(self.segdict, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def finalize(self):
        """
        Generate a markup.page object representing the HTML summary of the 
        segments associated with this SegmentSummaryTab.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises %s segments." % self.name,\
                     class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                     rel="full", class_="fancybox-button")
        self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],class_="full")
        self.frame.a.close()
 
        # construct livetime table
        uptime = abs(self.span)
        headers = ["Flags", "Livetime (s)", "Duty cycle (\%)"]
        data    = [[flag, float(abs(self.segdict[flag])),\
                    100*float(abs(self.segdict[flag]))/float(uptime)]\
                   for flag in self.flags]
        self.frame.add(htmlutils.write_table(headers, data, {"table":"full"})())
        self.frame.div.close()

        # other plots
        self.frame.h2("Plots", id_="h2_0-2", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-2")
        pclass = len(self.plots[1:]) == 1 and "full" or "half"
        for plot,desc in self.plots[1:]:
            self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                         rel="full")
            self.frame.img(src=plot, alt=desc, class_=pclass)
            self.frame.a.close()
        self.frame.div.close()

        # segment lists
        self.frame.h2("Segment lists", id_="h2_0-3", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-3")
        for i,flag in enumerate(self.flags):
            self.frame.h3(flag, id_="h3_0-3-%d" % i, class_="closed",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-3-%d" % i, style="display: none;")
            segfile = self.segment_files.get(flag, None)
            if segfile is not None:
                self.frame.p("The full segment list can be downloaded from %s."\
                             % markup.oneliner.a("this file", href=segfile),\
                             class_="line")
            segwizard = StringIO.StringIO()
            segmentsUtils.tosegwizard(segwizard, self.segdict[flag])
            self.frame.pre(segwizard.getvalue())
            segwizard.close()
            self.frame.div.close()
        self.frame.div.close()

        # subplots
        if len(self.subplots):
            self.frame.h2("Subplots", id_="h2_0-4", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            for plot,desc in self.subplots:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="subplots")
                self.frame.img(src=plot, alt=desc, class_="quarter")
            self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-5", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # dereference
        for attr in ["segments", "segdict"]:
            if hasattr(self, attr):
                try:
                    delattr(self, attr)
                except AttributeError:
                    raise


class DataSummaryTab(SummaryTab):
    """
    SummaryTab representing the summary of data extracted from frames.
    """
    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.timeseries  = dict()
        self.sampling    = dict()
        self.spectrum    = dict()
        self.minspectrum = dict()
        self.maxspectrum = dict()
        self.designspectrum = dict()
        self.referencespectrum = dict()
        self.spectrogram = dict()
        self.channels    = []
        self.framefiles  = dict()

    #
    # class methods
    #

    def add_timeseries(self, channel, series):
        """
        Add a glue.ligolw.table for a given channel to this object

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: TimeSeries data object to store
        @type series: C{lal.XXXXTimeSeries}
        """
        self.timeseries[channel] = series
        self.sampling[channel] = 1/series.deltaT
        if channel not in self.channels:
            self.channels.append(channel)

    def add_spectrogram(self, channel, *sequencelist):
        """
        Add a list of REAL8VectorSequence objects as the spectrogram for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param sequencelist: spectrogram data object to store, or list of
            spectrograms
        """
        if not self.spectrogram.has_key(channel):
            self.spectrogram[channel] = []
        if hasattr(sequencelist, "__contains__"):
            self.spectrogram[channel].extend(sequencelist)
        else:
            self.spectrogram[channel].append(sequencelist)
        if channel not in self.channels:
            self.channels.append(channel)

    def add_median_spectrum(self, channel, series):
        """
        Add a FrequencySeries object as the median spectrum for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: FrequencySeries data object to store
        @type series: C{lal.XXXXFrequencySeries}
        """
        self.spectrum[channel] = series
        if channel not in self.channels:
            self.channels.append(channel)

    def add_min_spectrum(self, channel, series):
        """
        Add a FrequencySeries object as the minimum spectrum for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: FrequencySeries data object to store
        @type series: C{lal.XXXXFrequencySeries}
        """
        if not series.name.endswith("min"):
            series.name = "%s_min" % series.name
        self.minspectrum[channel] = series
        if channel not in self.channels:
            self.channels.append(channel)

    def add_max_spectrum(self, channel, series):
        """
        Add a REAL8FrequencySeries object as the maximum spectrum for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: FrequencySeries data object to store
        @type series: C{lal.XXXXFrequencySeries}
        """
        if not series.name.endswith("max"):
            series.name = "%s_max" % series.name
        self.maxspectrum[channel] = series
        if channel not in self.channels:
            self.channels.append(channel)
        
    def add_design_spectrum(self, channel, series):
        """
        Add a FrequencySeries object as the maximum spectrum for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: FrequencySeries data object to store
        @type series: C{lal.XXXXFrequencySeries}
        """
        self.designspectrum[channel] = series
        
    def add_reference_spectrum(self, channel, series):
        """
        Add a FrequencySeries object as the maximum spectrum for this
        channel.

        @param channel: name of channel to store
        @type channel: C{str}
        @param series: FrequencySeries data object to store
        @type series: C{lal.XXXXFrequencySeries}
        """
        self.referencespectrum[channel] = series
        
    def plottimeseries(self, outfile, channels=None, subplot=False, **kwargs):
        """
        Plot the timeseries for this TriggerSummary, one column against another.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot, default: all of them
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotdata.plottimeseries
        """
        if not channels:
            channels = self.channels
        desc = kwargs.pop("description", None)
        if not subplot:
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        data = [seriesutils.duplicate(self.timeseries[channel])\
                for channel in channels]
        for i,channel in enumerate(channels):
            data[i].name = channel
        plotdata.plottimeseries(data, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plothistogram(self, outfile, channels=None, subplot=False, **kwargs):
        """
        Plot the histogram of data contained with the timeseries' for this
        DataSummary.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot, default: all of them
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotdata.plothistogram
        """
        if not channels:
           channels = self.channels
        desc = kwargs.pop("description", None)
        data = [self.timeseries[channel] for channel in channels]
        for i,channel in enumerate(channels):
            data[i].name = channel
        if len(channels) == 1:
            data[0].name = "_"
        plotdata.plothistogram(data, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotspectrogram(self, outfile, channel=None, ratio=None, subplot=False,\
                        **kwargs):
        """
        Plot the spectrogram of the given channel

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot, default: all of them
        @type channels: C{list}
        @param ratio: description of ratio to calculate on-the-fly,
            e.g. "median", or "design"
        @type ratio: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotdata.plotspectrogram
        """
        if not channel:
            channel = self.channels[0]
        desc = kwargs.pop("description", None)
        if not subplot:
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        kwargs.setdefault("title", latex(channel))
 
        # get data
        if len(self.spectrogram[channel]):
            data,epoch,deltaT,f0,deltaF = map(list,\
                                              zip(*self.spectrogram[channel]))
        else:
            data = []
            epoch = lal.LIGOTimeGPS(self.start_time)
            deltaT = 1
            f0 = 0
            deltaF = 1

        # get ratio
        if ratio:
            if ratio == "median":
                spectrum = self.spectrum[channel].data.data
            elif ratio == "max":
                spectrum = self.maxspectrum[channel].data.data
            elif ratio == "min":
                spectrum = self.minspectrum[channel].data.data
            elif ratio == "design":
                spectrum = self.designspectrum[channel].data.data
            elif ratio == "reference":
                spectrum = self.referencespectrum[channel].data.data
            else:
                raise ValueError("Unknown ratio mode \"%s\"" % ratio)
            for i,sequence in enumerate(data):
                data[i] =\
                    lal.CreateREAL8VectorSequence(sequence.length,\
                                                      sequence.vectorLength)
                data[i].data = sequence.data/spectrum
                numpy.putmask(data[i].data, numpy.isnan(data[i].data), 1e100)

        plotdata.plotspectrogram(data, outfile, epoch=epoch, deltaT=deltaT,\
                                 f0=f0, deltaF=deltaF, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotspectrum(self, outfile, channel=None, psd=False, subplot=False,\
                     **kwargs):
        """
        Plot the spectrum of the given channel

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot, default: all of them
        @type channels: C{list}
        @param psd: plot spectrum as a PSD, default False
        @type psd: C{bool}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plotdata.plotfrequencyseries
        """
        if not channel:
            channel = self.channels[0]
        desc = kwargs.pop("description", None)
        serieslist = [self.spectrum[channel], self.minspectrum[channel],\
                      self.maxspectrum[channel]]
        if self.designspectrum.has_key(channel):
            serieslist.append(self.designspectrum[channel])
        if self.referencespectrum.has_key(channel):
            serieslist.append(self.referencespectrum[channel])
        if psd:
            for i,series in serieslist:
                serieslist[i] =\
                    seriesutils.fromarray(series.data.data**2,\
                                          name=series.name, epoch=series.epoch,\
                                          deltaT=series.deltaF, f0=series.f0,\
                                          sampleUnits=series.sampleUnits,\
                                          frequencyseries=True)
        plotdata.plotfrequencyseries(serieslist, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def finalize(self):
        """
        Generate a markup.page object representing the HTML summary of the 
        data associated with this DataSummaryTab.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises %s data." % self.name,\
                     class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        if len(self.plots):
            self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                         rel="full", class_="fancybox-button")
            self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],\
                           class_="full")
            self.frame.a.close()
 
        # construct info table
        headers = ["Channel", "Sampling rate"]
        data    = [[channel, self.sampling[channel]]\
                   for channel in self.channels\
                   if not re.search("[-._](min|max)\Z", channel)]
        self.frame.add(htmlutils.write_table(headers, data, {"table":"full"})())
        self.frame.div.close()

        # other plots
        self.frame.h2("Plots", id_="h2_0-2", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-2")
        pclass = len(self.plots[1:]) == 1 and "full" or "half"
        for plot,desc in self.plots[1:]:
            self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                         rel="full")
            self.frame.img(src=plot, alt=desc, class_=pclass)
            self.frame.a.close()
        self.frame.div.close()

        # subplots
        if len(self.subplots):
            self.frame.h2("Subplots", id_="h2_0-4", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            for plot,desc in self.subplots:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="subplots")
                self.frame.img(src=plot, alt=desc, class_="quarter")
            self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-5", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # dereference
        for attr in ["timeseries", "spectrum", "minspectrum", "maxspectrum",\
                     "spectrogram"]:
            if hasattr(self, attr):
                delattr(self, attr)


class RangeSummaryTab(DataSummaryTab):
    """
    SummaryTab representing a summary of detection range extracted from
    frames directly or calculated from h(t).
    """
    def __init__(self, *args, **kwargs):
        DataSummaryTab.__init__(self, *args, **kwargs)
        self.sources = list()

    def add_source(self, sourcedict):
        self.sources.append(sourcedict)
        self.sources.sort(key=lambda s: (s["type"],s["name"]))

    def finalize(self):
        """
        Generate a markup.page object representing the HTML summary of the 
        data associated with this RangeSummaryTAb.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises %s range data." % self.name,\
                     class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        if len(self.plots):
            self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                         rel="full", class_="fancybox-button")
            self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],\
                           class_="full")
            self.frame.a.close()
 
        # construct info table
        headers = ["Source", "Type", "Mean range (Mpc)", "Min range (Mpc)",\
                   "Max range (Mpc)"]
        td = []
        for s in self.sources:
            data = self.timeseries[s["name"]].data.data
            data = data[data!=0]
            if data.size:
                td.append([s["name"], s["type"], data.mean(), data.min(),\
                             data.max()])
            else:
                td.append([s["name"], s["type"], "-", "-", "-"])
        self.frame.add(htmlutils.write_table(headers, td, {"table":"full"})())
        self.frame.div.close()

        # other plots
        self.frame.h2("Plots", id_="h2_0-2", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-2")
        pclass = len(self.plots[1:]) == 1 and "full" or "half"
        for plot,desc in self.plots[1:]:
            self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                         rel="full")
            self.frame.img(src=plot, alt=desc, class_=pclass)
            self.frame.a.close()
        self.frame.div.close()

        # subplots
        if len(self.subplots):
            self.frame.h2("Subplots", id_="h2_0-4", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            for plot,desc in self.subplots:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="subplots")
                self.frame.img(src=plot, alt=desc, class_="quarter")
            self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-5", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # dereference
        for attr in ["timeseries"]:
            if hasattr(self, attr):
                delattr(self, attr)


class TriggerSummaryTab(SummaryTab):
    """
    Object representing a summary of triggers.
    """
    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.triggers = dict()
        self.channels = []
        self.trigger_files = dict()

    def add_triggers(self, channel, trigtable):
        """
        Add a glue.ligolw.table for a given channel to this object

        @param channel: name of channel to record
        @type channel: C{str}
        @param trigtable: LIGOLw table of triggers to record
        @type trigtable: C{glue.ligolw.table.Table}
        """
        if self.triggers.has_key(channel):
            self.triggers[channel].extend(trigtable)
        else:
            self.triggers[channel] = trigtable
        if channel not in self.channels:
            self.channels.append(channel)

    def toxml(self, channel, filename=None):
        """
        Write the trigtable for the given channel to an xml file.

        @param channel: name of channel to record
        @type channel: C{str}
        @param filename: path of file to write, preferrably xml.gz extension
        @type filename: C{str}
        """
        if not filename:
            name = _r_cchar.sub("-", _r_cchar.sub("_", channel.upper()), 1)
            filename = os.path.join(self.directory, "%s-%d-%d.txt"\
                       % (name, self.start_time, abs(self.span)))
        with open(filename, "w") as trigf:
            dqTriggerUtils.totrigxml(trigf, self.triggers[channel],\
                                     program=sys.argv[0])

    def plottable(self, outfile, channels=None, subplot=False,\
                  **kwargs):
        """
        Plot the triggers for this TriggerSummary, one column against another.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plottriggers.plottable
        """
        desc = kwargs.pop("description", None)
        if kwargs.get("xcolumn", None) == "time":
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        if channels and len(channels) == 1:
            kwargs.setdefault("title", "%s (%s)" % (latex(channels[0]),\
                                                    latex(self.etg)))
        if channels is not None:
            trigs = dict((key,val) for (key,val) in self.triggers.iteritems()\
                         if key in channels)
            plottriggers.plottable(trigs, outfile, **kwargs)
        if len(self.triggers.keys()) == 1:
            kwargs.setdefault("title","%s (%s)" % (latex(self.channels[0]),\
                                                   latex(self.etg)))
            plottriggers.plottable({"_":self.triggers.values()[0]}, outfile,\
                                   **kwargs)
        else:
            plottriggers.plottable(self.triggers, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plothistogram(self, outfile, channels=None, subplot=False, **kwargs):
        """
        Plot a histogram of the triggers for this TriggerSummary.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plottriggers.plothistogram
        """
        desc = kwargs.pop("description", None)
        if kwargs.get("cumulative", True) and not kwargs.has_key("normalize"):
            kwargs["normalize"] = float(abs(self.span))
        if channels is not None:
            trigs = dict((key,val) for key in self.triggers.keys() if key\
                         in channels)
            plottriggers.plothistogram(trigs, outfile, **kwargs)
        elif len(self.triggers.keys()) == 1:
            plottriggers.plothistogram({"_":self.triggers.values()[0]},\
                                       outfile, **kwargs)
        else:
            plottriggers.plothistogram(self.triggers, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotrate(self, outfile, channels=None, subplot=False, **kwargs):
        """
        Plot the rate of triggers for this TriggerSummary.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plottriggers.plotrate
        """
        desc = kwargs.pop("description", None)
        kwargs.setdefault("xlim", [self.start_time, self.end_time])
        if channels:
            trigs = dict((key,val) for key in self.triggers.keys() if key\
                         in channels)
            plottriggers.plotrate(trigs, outfile, **kwargs)
        elif len(self.triggers.keys()) == 1:
            plottriggers.plotrate({"_":self.triggers.values()[0]}, outfile,\
                                   **kwargs)
        else:
            plottriggers.plotrate(self.triggers, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotautocorrelation(self, outfile, channels=None, subplot=False,\
                            **kwargs):
        """
        Plot the auto correlation function of the trigers for this
        TriggerSummary.

        @param outfile: filename for plot
        @type outfile: C{str}
        @param channels: list of channels to plot
        @type channels: C{list}
        @param subplot: record plot as a subplot, default: False
        @type subplot: C{bool}
        @keyword **kwargs: other arguments to pass to
            plottriggers.plotautocorrelation
        """
        desc = kwargs.pop("description", None)
        if kwargs.get("xcolumn", None) == "time":
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        if channels:
            trigs = dict((key,val) for key in self.triggers.keys() if key\
                         in channels)
            plottriggers.plotautocorrelation(trigs, outfile, **kwargs)
        elif len(self.triggers.keys()) == 1:
            plottriggers.plotautocorrleation({"_":self.triggers.values()[0]},\
                                             outfile, **kwargs)
        else:
            plottriggers.plotautocorrelation(self.triggers, outfile, **kwargs)
        if subplot:
            self.subplots.append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def finalize(self, snrs=[5,8,10,20,50,100,500,1000]):
        """
        Generate a markup.page object summarising the triggers in this
        TriggerSummaryTab.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises %s triggers." % self.name,\
                     class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                     rel="full", class_="fancybox-button")
        self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],class_="full")
        self.frame.a.close()

        # trig stats
        th = ["Channel"]+["SNR >= %d" % s for s in snrs]
        td = []
        for chan in self.channels:
            snr = numpy.asarray(self.triggers[chan].getColumnByName("snr"))
            td.append([chan]+[(snr>s).sum() for s in snrs]) 
        self.frame.add(htmlutils.write_table(th, td, {"table":"full"})())
        self.frame.div.close()

        # other plots
        self.frame.h2("Plots", id_="h2_0-2", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-2")
        pclass = len(self.plots[1:]) == 1 and "full" or "half"
        for plot,desc in self.plots[1:]:
            self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                         rel="full")
            self.frame.img(src=plot, alt=desc, class_=pclass)
            self.frame.a.close()
        self.frame.div.close()
 
        # loudest triggers
        self.frame.h2("Loudest triggers", id_="h2_0-3", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-3")
        for i,chan in enumerate(self.channels):
            self.frame.h3(chan, id_="h3_0-3-%d" % i, class_="closed",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-3-%d" % i, style="display: none;")
            trigfile = self.trigger_files.get(chan, None)
            if trigfile is not None:
                self.frame.p("The full segment list can be downloaded from %s."\
                             % markup.oneliner.a("this file", href=trigfile),\
                             class_="line")
            # FIXME write list of loudest triggers
            self.triggers[chan].sort(key=lambda t: t.snr, reverse=True)
            trigs = table.new_from_template(self.triggers[chan])
            trigs.extend(self.triggers[chan][:20])
            data = []
            data.append(plottriggers.get_column(trigs, "time"))
            head = ["Time"]
            if "central_freq" in trigs.validcolumns.keys():
                data.append(plottriggers.get_column(trigs, "duration"))
                data.append(plottriggers.get_column(trigs, "central_freq"))
                data.append(plottriggers.get_column(trigs, "bandwidth"))
                head.extend(["Duration", "Frequency", "Bandwidth"])
            else:
                data.append(plottriggers.get_column(trigs, "template_duration"))
                data.append(plottriggers.get_column(trigs, "mchirp"))
                head.extend(["Template duration", "Chirp mass"])
            data.append(plottriggers.get_column(trigs, "snr"))
            head.append("SNR")
            data = map(list, zip(*data))
            self.frame.add(htmlutils.write_table(head,data,{"table":"full"})())
            self.frame.div.close()
        self.frame.div.close()

        # subplots
        if len(self.subplots):
            self.frame.h2("Subplots", id_="h2_0-4", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            for plot,desc in self.subplots:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="subplots")
                self.frame.img(src=plot, alt=desc, class_="quarter")
            self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-5", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # dereference
        for attr in ["triggers"]:
            if hasattr(self, attr):
                delattr(self, attr)


class AuxTriggerSummaryTab(TriggerSummaryTab):
    """
    Object representing a summary of auxiliary channel triggers.
    """
    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.triggers      = dict()
        self.channels      = []
        self.mainchannel   = []
        self.trigger_files = dict()
        self.coincs        = dict()
        self.numcoincs     = dict()
        self.sigma         = dict()
        self.auxplots      = dict()

    def get_coincs(self, channel1, channel2, dt=0.03):
        """
        Find triggers for channel1 within dt of a trigger for channel2.
        """
        self.coincs[(channel1, channel2)] =\
            dqTriggerUtils.get_coincs(self.triggers[channel1],\
                                      self.triggers[channel2], dt=dt)
        self.numcoincs[(channel1, channel2)] =\
            len(self.coincs[(channel1, channel2)])

    def compute_coinc_significance(self, channel, dt=0.03, shifts=[0]):
        """
        Calculate the statistical significance of the number of triggers in
        channel coincident with a trigger in the mainchannel.
        Give a list of numerical time shifts to calculate for multiple events.
        """
        trigs1 = self.triggers[channel]
        trigs2 = self.triggers[self.mainchannel]
        t = plottriggers.get_column(trigs1, "time")

        _reburst  = re.compile("burst", re.I)
        _rering   = re.compile("ring", re.I)
        _recbc    = re.compile("inspiral", re.I)
        shifts    = numpy.asarray(shifts)
        ncoinc    = numpy.zeros(len(shifts))
        for i,shift in enumerate(shifts):
            if shift==0 and self.numcoincs.has_key((self.mainchannel, channel)):
                ncoinc[i] = self.numcoincs[(self.mainchannel, channel)]
            else:
                ncoinc[i] = dqTriggerUtils.get_number_coincs(trigs2, trigs1,\
                                                             dt=dt,\
                                                             timeshift=shift)
        
        # calculate numsigma for each slide
        mean  = ncoinc[shifts!=0].mean()
        std   = ncoinc[shifts!=0].std()
        self.sigma[channel] = dict(zip(shifts, numpy.fabs(ncoinc-mean)/std))

    def plottable(self, outfile, channel, **kwargs):
        """
        Plot the triggers for the given channel.
        """
        desc = kwargs.pop("description", None)
        if kwargs.get("xcolumn", None) == "time":
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        kwargs.setdefault("title", "%s (%s)" % (latex(channel),\
                                                latex(self.etg)))
        plottriggers.plottable({"_":self.triggers[channel]}, outfile,\
                               **kwargs)
        self.auxplots[channel].append((outfile, desc))

    def plothistogram(self, outfile, channel, **kwargs):
        desc = kwargs.pop("description", None)
        plottriggers.plothistogram({"_":self.triggers[channel]}, outfile,\
                                   **kwargs)
        self.auxplots[channel].append((outfile, desc))

    def plotrate(self, outfile, channel, **kwargs):
        desc = kwargs.pop("description", None)
        kwargs.setdefault("xlim", [self.start_time, self.end_time])
        plottriggers.plotrate({"_":self.triggers[channel]}, outfile,\
                               **kwargs)
        self.auxplots[chan].append((outfile, desc))

    def plotautocorrelation(self, outfile, channel, **kwargs):
        desc = kwargs.pop("description", None)
        if kwargs.get("xcolumn", None) == "time":
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        plottriggers.plotautocorrleation({"_":self.triggers[channel]},\
                                         outfile, **kwargs)
        self.auxplots[chan].append((outfile, desc))

    def plotcoincs(self, outfile, channel, **kwargs):
        """
        Plot the coincident triggers between the given channel and the
        mainchannel.
        """
        desc = kwargs.pop("description", None)
        if kwargs.get("xcolumn", None) == "time":
            kwargs.setdefault("xlim", [self.start_time, self.end_time])
        kwargs.setdefault("title", "Coincident %s and %s (%s)"\
                                   % (latex(channel), latex(self.mainchannel),\
                                      latex(self.etg)))
        trigs = {channel:self.coincs[(channel, self.mainchannel)],\
                 self.mainchannel:self.coincs[(self.mainchannel, channel)]}
        plottriggers.plottable(trigs, outfile, **kwargs)
        self.auxplots[channel].append((outfile, desc))
        
    def plotsigma(self, outfile, channel=None, **kwargs):
        """
        Plot the statistical significance of the number of coincidences between
        the mainchannel and all auxiliary channels.
        """
        desc = kwargs.pop("description", None)
        xdata = []
        ydata = []

        # get data for one channel
        if channel:
            xdata,ydata = zip(*sorted(self.sigma[channel].items(),\
                                      key=lambda (x,y): x))
        # or data for all channels
        else:
            for chan in self.channels:
                if chan == self.mainchannel:
                    continue
                if self.sigma[chan].has_key(0.0):
                    ydata.append(self.sigma[chan][0.0])
                    xdata.append(chan)

        # get axis params
        kwargs.pop("xlim", None)
        ylim = kwargs.pop("ylim", None)
        kwargs.pop("logx", False)
        logy = kwargs.pop("logy", False)

        # get labels
        if channel:
            xlabel = kwargs.pop("xlabel", "Time shift (s)")
        else:
            xlabel = kwargs.pop("xlabel", "")
        ylabel = kwargs.pop("ylabel", r"Coincidence significance "\
                                       "($\mathrm{\sigma}$)")
        if channel:
            title  = kwargs.pop("title", "Coincident %s and %s (%s)"\
                                % (latex(channel), latex(self.mainchannel),\
                                   latex(self.etg)))
        else:
            title  = kwargs.pop("title",\
                                "Significance of coincidences with %s (%s)"\
                                % (latex(self.mainchannel), latex(self.etg)))
        subtitle = kwargs.pop("subtitle", "")

        # get misc params
        bbox = kwargs.pop("bbox_inches", None)
        cbar = kwargs.pop("hidden_colorbar", None)

        # make plot
        if channel:
            plot = plotdata.plotutils.BarPlot(xlabel, ylabel, title, subtitle)
        else:
            plot = plotdata.plotutils.BarPlot(xlabel, ylabel, title, subtitle,\
                                              figsize=[24,6])
        plot.add_content(numpy.arange(len(xdata)), ydata, **kwargs)
        plot.finalize()
        if cbar:
            plotdata.plotutils.add_colorbar(plot.ax, visible=False)      
 
        # set xticks to channel names and rotate
        plot.ax.set_xlim(-1, len(xdata))
        plot.ax.set_xticks(numpy.arange(0,len(xdata)))
        if channel:
            plot.ax.set_xticklabels(map(lambda x: "$%s$"\
                                        % plotdata.plotutils.float_to_latex(x),\
                                        xdata))
        else:
            plot.ax.set_xticklabels(map(latex, xdata))
            for i,t in enumerate(plot.ax.get_xticklabels()):
                t.set_rotation(315)
                t.set_verticalalignment('top')
                t.set_horizontalalignment('left')
                t.set_fontsize("smaller")
            plot.fig.subplots_adjust(bottom=0.3)

        if ylim:
            plot.ax.set_ylim(ylim)

        plot.savefig(outfile, bbox_inches=bbox)
        plot.close() 
        if channel:
            self.auxplots[channel].append((outfile, desc))
        else:
            self.plots.append((outfile, desc))

    def plotnumslides(self, outfile, channel, **kwargs):
        """
        Plot the number of coincidences between this channel and the
        mainchannel for all slides and zerolag.
        """
        desc = kwargs.pop("description", None)
        data = sorted(self.sigma[channel].items(), key=lambda x: x[0])
        shifts,sigma = map(numpy.asarray, zip(*data))

        # get axis params        
        kwargs.pop("xlim")
        xlim = [shifts.min() - abs(shifts.min())*0.01,\
                shifts.max() + abs(shifts.max())+0.01]
        ylim = kwargs.pop("ylim", None)
        kwargs.pop("logx")
        logy = kwargs.pop("logy")

        # get labels
        xlabel = kwargs.pop("xlabel", "Time shift (s)")
        ylabel = kwargs.pop("ylabel", "Number of coincidences")
        title  = kwargs.pop("title", "")
        subtitle = kwargs.pop("subtitle", "")

        # get misc params
        bbox = kwargs.pop("bbox_inches", None)
        cbar = kwargs.pop("hidden_colorbar", None)

        # make plot
        plot = plotdata.plotutils.SimplePlot(xlabel, ylabel, title, subtitle)
        plot.add_content(shift, sigma, **kwargs)
        plot.finalize()
        if logy:
            plot.ax.sey_yscale("log")
        plot.ax.set_xlim(xlim)
        if logy:
            plot.ax.set_ylim(ylim)
        if cbar:
            plotdata.plotutils.add_colorbar(plot.ax, visible=False)      
 
        plotdata.plotutils.set_ticks(plot.ax, x=False, y=True)
        plot.savefig(outfile, bbox_inches=bbox)
        plot.close() 
        self.auxplots[chan].append((outfile, desc))

    def finalize(self):
        """
        Generate a glue.markup.page summarising the auxiliary channel triggers
        for this AuxTriggerSummaryTab.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises auxiliary channel %s triggers."\
                     % self.name, class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        if len(self.plots):
            self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                         class_="fancybox-button", rel="full")
            self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],\
                           class_="full")
            self.frame.a.close()
  
        # trig stats
        if len(self.numcoincs.keys()) != 0:
            th = ["Channel", "Num. coinc. with<br>%s" % self.mainchannel,\
                  "Num. %s<br>coinc. with aux." % self.mainchannel,\
                  "Zero shift coinc. &sigma;"]
            td = []
            for chan in self.channels:
                if chan == self.mainchannel:
                    continue
                td.append([chan])
                if (chan, self.mainchannel) in self.numcoincs.keys():
                    td[-1].extend([self.numcoincs[(chan, self.mainchannel)],\
                                   self.numcoincs[(self.mainchannel, chan)]])
                else:
                    td[-1].extend(["-", "-"])
                if self.sigma.has_key(chan) and self.sigma[chan].has_key(0.0):
                    td[-1].append("%.2f" % self.sigma[chan][0.0])
                else:
                    td[-1].append("-")
            self.frame.add(htmlutils.write_table(th, td, {"table":"full"})())
         
        self.frame.div.close()

        # channel information
        self.frame.h2("Auxiliary channels", id_="h2_0-2", class_="closed",\
                      onclick="toggleVisible();")
        self.frame.div(id_="div_0-2", style="display: none;")

        for i,chan in enumerate(self.channels):
            if chan == self.mainchannel:
                continue
            self.frame.h3(chan, id_="h3_0-2-%d" % i, class_="closed",\
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-2-%d" % i, style="display: none;")
            if (chan, self.mainchannel) in self.numcoincs:
                th = ["Num. coinc with<br>%s" % (self.mainchannel),\
                      "Num %s<br>coinc with aux." % (self.mainchannel),\
                      "Zero time-shift coincidence &sigma;"]
                td = list()
                if (chan, self.mainchannel) in self.numcoincs.keys():
                    td.extend([self.numcoincs[(chan, self.mainchannel)],\
                               self.numcoincs[(self.mainchannel, chan)]])
                else:
                    td.extend(["-", "-"])
                if self.sigma.has_key(chan) and self.sigma[chan].has_key(0.0):
                    td.append("%.2f" % self.sigma[chan][0.0])
                else:
                    td.append("-")
                self.frame.add(htmlutils.write_table(th,td,{"table":"full"})())
            class_ = plotclass(len(self.auxplots[chan]))
            for p,d in self.auxplots[chan]:
                self.frame.a(href=p, title=d, class_="fancybox-button", rel=class_)
                self.frame.img(src=p, alt=d,  class_=class_)
                self.frame.a.close()
            self.frame.div.close()
        self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-3", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-3")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # deference
        for attr in ["triggers", "coincs"]:
            if hasattr(self, attr):
                delattr(self, attr)


class StateVectorSummaryTab(SegmentSummaryTab):
    """
    Object representing the summary of a bitmasked channel.
    """
    def __init__(self, *args, **kwargs):
        SegmentSummaryTab.__init__(self, *args, **kwargs)
        self.bitmask = dict()
        self.derived_bitmask = dict()

    def add_bitmask(self, bitmask):
        self.bitmask = dict(bitmask)
        self.bits = sorted(self.bitmask.keys())

    def add_derived_bits(self, derived_bitmask):
        self.derived_bitmask = dict(derived_bitmask)
        self.derived_bits = sorted(self.derived_bitmask.keys())

    def finalize(self):
        """
        Generate a markup.page object representing the HTML summary of the 
        segments associated with this StateVectorSummaryTab.
        """
        # opening
        self.frame = markup.page()
        self.frame.h1(self.name, id_="h1_0", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0")
        self.frame.p("This page summarises %s state vector." % self.name,\
                     class_="line")

        # summary
        self.frame.h2("Summary", id_="h2_0-1", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-1")
        self.frame.a(href=self.plots[0][0], title=self.plots[0][1],\
                     class_="fancybox-button", rel="full")
        self.frame.img(src=self.plots[0][0], alt=self.plots[0][1],\
                       class_="full")
        self.frame.a.close()
 
        # construct livetime table
        uptime = float(abs(self.span))
        headers = ["Bit", "Flag", "Livetime (s)", "Duty cycle (%%)"]
        data    = [[bit, self.bitmask[bit],\
                    abs(self.segdict[self.bitmask[bit]]),\
                    "%.2f" %\
                        (100*(abs(self.segdict[self.bitmask[bit]])/uptime))]\
                    for bit in self.bits]
        self.frame.add(htmlutils.write_table(headers, data, {"table":"full"})())
        self.frame.div.close()

        # other plots
        if len(self.plots) > 1:
            self.frame.h2("Plots", id_="h2_0-2", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-2")
            for plot,desc in self.plots[1:]:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="half")
                self.frame.img(src=plot, alt=desc, class_="half")
                self.frame.a.close()
            self.frame.div.close()

        # segment lists
        self.frame.h2("Segment lists", id_="h2_0-3", class_="open",\
                     onclick="toggleVisible();")
        self.frame.div(id_="div_0-3")
        for i,bit in enumerate(self.bits):
            flag = self.bitmask[bit]
            self.frame.h3(flag, id_="h3_0-3-%d" % i,\
                          class_="closed", onclick="toggleVisible();")
            self.frame.div(id_="div_0-3-%d" % i, style="display: none;")
            segfile = self.segment_files.get(flag, None)
            if segfile is not None:
                self.frame.p("The full segment list can be downloaded from %s."\
                             % markup.oneliner.a("this file", href=segfile),\
                             class_="line")
            segwizard = StringIO.StringIO()
            segmentsUtils.tosegwizard(segwizard, self.segdict[flag])
            self.frame.pre(segwizard.getvalue())
            segwizard.close()
            self.frame.div.close()
        self.frame.div.close()

        # subplots
        if len(self.subplots):
            self.frame.h2("Subplots", id_="h2_0-4", class_="open",\
                         onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            for plot,desc in self.subplots:
                self.frame.a(href=plot, title=desc, class_="fancybox-button",\
                             rel="subplots")
                self.frame.img(src=plot, alt=desc, class_="quarter")
            self.frame.div.close()

        # information
        if self.information:
            self.frame.h2("Information", id_="h2_0-5", class_="open",
                          onclick="toggleVisible();")
            self.frame.div(id_="div_0-4")
            self.frame.add(self.information)
            self.frame.div.close()

        self.frame.div.close()

        # deference
        for attr in ["segdict"]:
            if hasattr(self, attr):
                delattr(self, attr)


class GlossaryTab(SummaryTab):

    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.tabs = ""
        self.definitions = dict()

    def add_entry(key, val):
        self.definitions[key] = val

    def add_entries(*args):
        for key,val in args:
            self.add_entry(key, val)

    def finalize(self):
        """
        Write the GlossaryTab HTML frame.
        """
        self.frame = htmlutils.write_glossary(self.definitions)


class AboutTab(SummaryTab):

    version = None
    files   = list()

    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.tabs = ""

    def add_file(self, name, path):
        self.files.append((name, path))

    def finalize(self):
        """
        Write the AboutTab HTML frame.
        """
        self.frame = htmlutils.about_page(sys.argv[0], sys.argv[1:],\
                                          self.version, self.files)


class CalendarTab(SummaryTab):

    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.weekday = 0
        self.tabs = ""

    @property
    def start_date(self):
        """datetime.datetime start date for this CalendarTab."""
        return self._start_date
    @start_date.setter
    def start_date(self, date):
        if isinstance(date, str):
            if re.search("/", startdate):
               self._start_date = datetime.datetime.strptime(date, "%Y/%m/%d")
            else:
               self._start_date = datetime.datetime.strptime(date, "%Y%m%d")
        elif isinstance(date, int):
            self._start_date = datetime.date(*lal.GPSToUTC(int(date))[:3])
        elif isinstance(date, datetime.datetime)\
        or isinstance(date, datetime.date):
            self._start_date = date

    @property
    def end_date(self):
        """datetime.datetime end date for this CalendarTab."""
        return self._end_date
    @end_date.setter
    def end_date(self, date):
        if isinstance(date, str):
            if re.search("/", enddate):
               self._end_date = datetime.datetime.strptime(date, "%Y/%m/%d")
            else:
               self._end_date = datetime.datetime.strptime(date, "%Y%m%d")
        elif isinstance(date, int):
            self._end_date = datetime.date(*lal.GPSToUTC(int(date))[:3])
        elif isinstance(date, datetime.datetime)\
        or isinstance(date, datetime.date):
            self._end_date = date

    #
    # class methods
    #

    def finalize(self):
        """
        Returns a glue.markup.page containing a nested HTML table of links
        for the given SummaryTab.
        """
        # get startday from configuration
        years = numpy.arange(self.start_date.year, self.end_date.year+1)

        # make page
        self.frame = markup.page()
        self.frame.h1('Calendar')
        for i,y in enumerate(years[::-1]):
            startdate = datetime.date(y, 1, 1)
            enddate   = min(self.end_date, datetime.date(y, 12, 31))
            self.frame.h2(y, id="h2_%d" % i, onclick="toggleVisible();",\
                    class_="calendar %s" % (i==0 and 'open' or 'closed'))
            self.frame.div(id="div_%d" % i, class_="calendar",\
                     style="display: %s;" % (i==0 and 'block' or 'none'))
            self.frame.add(calendar_page(startdate, enddate,\
                                         weekday=self.weekday)())
            self.frame.div.close()


class OnlineSummaryTab(SummaryTab):
    """
    SummaryTab representing the online data page.
    """
    def __init__(self, *args, **kwargs):
        SummaryTab.__init__(self, *args, **kwargs)
        self.tabs     = ""
        self.states   = list()
        self.sections = list()
        self.plots    = dict()
        self.refresh  = None

    @property
    def plots(self):
        """dict of plots for this OnlineSummaryTab."""
        return self._plotdict
    @plots.setter
    def plots(self, plotdict):
        self._plotdict = plotdict
        if isinstance(plotdict, dict):
            for key,val in plotdict.iteritems():
                for key2,val2 in val.iteritems():
                    self._plotdict[key][key2] =\
                        map(os.path.normpath(plotdict[key][key2]))
    @plots.deleter
    def plots(self):
        del self._plotdict

    def add_states(self, statelist):
        self.states.append(statelist)
        for state in statelist:
            if not state in self.plots.keys():
                self.plots[state] = dict()

    def add_plots(self, state, plotlist):
        state = SummaryState(str(state))
        if not state in self.states:
            self.states.append(state)
            self.plots[state] = dict()
        for key,plot in plotlist:
            if key not in self.sections:
                self.sections.append(key)
            if key not in self.plots[state].keys():
                self.plots[state][key] = []
            self.plots[state][key].extend(plot.split(','))

    def write_tabs(self):
        """
        Write the tabbar used for this OnlineSummaryTab.
        """
        self.tabs = markup.page()
        self.tabs.ul(class_="buttons")
        for i,section in enumerate(sorted(self.sections)):
            self.tabs.li(section, class_="open", onclick="toggleVisible();",\
                         id_="button_%s" % _r_cchar.sub("-", section), \
                         title="Hide online %s plots" % section)
        self.tabs.ul.close()

    def finalize(self):
        """
        Write the HTML frame for this OnlineSummaryTab.
        """
        if not self.tabs:
            self.write_tabs()
        self.frame = markup.page()
        if len(self.states) == 0:
            self.frame.p("No online monitors configured.", class_="line")
        else:
            self.frame.p("This frame shows the current status of this "+\
                        "instrument. Select which sections to view using the "+\
                        "buttons above.", class_="line", id_="online")
            if self.refresh:
                self.frame.p("Plots will auto refresh every %d seconds"\
                       % self.refresh, class_="line")
            for i,state in enumerate(self.states):
                style = i==0 and "block" or "none"
                self.frame.div(id_="div_%s"\
                               % _r_cchar.sub("-", state.name.lower()),\
                               style="display: %s;" % style)
                for j,section in enumerate(sorted(self.plots[state].keys())):
                    self.frame.div(class_="toggle_%s"% _r_cchar.sub("-", section),\
                                   style="display: block;")
                    self.frame.h2(section, id_="h2_%d" % i, class_="open",\
                            onclick="toggleVisible();")
                    self.frame.div(id_="div_%d" % i, style="display: block;")
                    c = plotclass(len(self.plots[state][section]))
                    for k,plot in enumerate(self.plots[state][section]):
                        self.frame.a(id_="img_%d-%d" % (i,k), href=plot, rel=c)
                        self.frame.img(src=plot, id_="imd_%d-%d" % (i,k),\
                                       class_="%s online" % c)
                        self.frame.a.close()
                    self.frame.div.close()
                    self.frame.div.close()
                self.frame.div.close()
        if self.refresh:
            self.frame.script("var t=%s; refreshImages(t);",\
                             type="text/javascript")

# =============================================================================
# Define run state object
# =============================================================================

class SummaryState(object):
    """
    Object representing a choice of IFO state.
    """
    def __init__(self, name):
        """
        Define a new SummaryState with the given name

        @param name: descriptive name for this SummaryState
        @type name: C{str}
        @return: a new SummaryState object
        @rtype: C{summary.SummaryState}
        """
        # set default variables
        self.name = name
        self.definition = None
        self.segments = segments.segmentlist()
        self.set = False

        # define the tag and match
        self.tag = _r_cchar.sub("_", name.lower()).upper()
        self.match = re.compile("%s\Z" % self.tag, re.I).match

    @property
    def start_time(self):
        """GPS start time for this Tab."""
        return self._start_time
    @start_time.setter
    def start_time(self, gps):
        self._start_time = lal.LIGOTimeGPS(gps)
    @start_time.deleter
    def start_time(self):
        del self._start_time

    @property
    def end_time(self):
        """GPS end time for this Tab."""
        return self._end_time
    @end_time.setter
    def end_time(self, gps):
        self._end_time = lal.LIGOTimeGPS(gps)
    @end_time.deleter
    def end_time(self):
        del self._end_time

    @property
    def span(self):
        """GPS [start_time, stop_time) segment for this Tab."""
        return segments.segment(self.start_time, self.end_time)
    @span.setter
    def span(self, seg):
        self.start_time = lal.LIGOTimeGPS(seg[0])
        self.end_time = lal.LIGOTimeGPS(seg[1])
    @span.deleter
    def span(self):
        del self._span

    @property
    def segments(self):
        """glue.segments.segmentlist describing the valid segments for
        this SummaryState.
        """
        return self._segments
    @segments.setter
    def segments(self, seglist):
        self._segments =\
            segments.segmentlist([segments.segment(map(float, s))\
                                  for s in seglist])
    @segments.deleter
    def segments(self):
        del self._segments

# =============================================================================
# Help functions
# =============================================================================

def calendar_page(startdate, enddate, path=None, jobdir=".",\
                  weekday=calendar.MONDAY, ncol=4, reverse=False):
    """
    Write an HTML calendar for the given [gpsstart, gpsend) interval.
    One table per month, one link per day.
    """
    d = datetime.date(startdate.year, 1, 1)
    calendar.setfirstweekday(weekday)

    if reverse:
        m = enddate
    else:
        m = datetime.date(startdate.year, 1, 1)
 
    page = markup.page()

    if not path:
        dir1 = os.path.sep
        dir2 = os.path.sep
    else:
        path = os.path.normpath(path)
        dir1 = os.path.join("", *re.split(os.path.sep, path)[-2:])
        dir2 = os.path.join("", *re.split(os.path.sep, path)[-1:])
        if re.match("%s\d\d\d\d\d\d" % os.path.sep, dir1):
            dir1 = dir2

    # loop over months
    i = 0
    if ncol==1:
        page.table(class_="calendar")
    else:
        page.table(class_="calendar year")
    while startdate.year <= m.year < enddate.year+1:
        # link year
        if  i==0 or (i % ncol == 0 and ((reverse and m.month==12)\
                                        or (not reverse and m.month==1))):
            page.tr()
            Y = m.year
            href = None
            for e,t in zip(["html", "php"], [dir1, dir2]):
                idx = os.path.join(jobdir, "archive_yearly",\
                                   "%d%s" % (m.year, t), "index.%s" % e)
                if os.path.isfile(os.path.expanduser(idx)):
                    href = "%s%s" % (os.path.split(idx)[0], os.path.sep)
                    break
            if href:
                page.th(markup.oneliner.a(Y, href=href),
                        class_="year", colspan="100%")
            else:
                page.th(Y, class_="year", colspan="100%")
            page.tr.close()
 
            
        # open new table
        if i % ncol == 0:
            page.tr()
        # get delta
        if reverse and m.month == 1:
            delta = datetime.timedelta(days=-calendar.monthrange(m.year-1,\
                                                                 12)[1])
        elif reverse:
            delta = datetime.timedelta(days=-calendar.monthrange(m.year,\
                                                                m.month-1)[1])
        else:
            delta = datetime.timedelta(days=calendar.monthrange(m.year,\
                                                                m.month)[1])

        # if month is outside range, make a blank entry for it
        if (m.year == startdate.year and m.month < startdate.month)\
        or (m.year == enddate.year and m.month > enddate.month):
            if ncol==1:
                page.td(str())
            else:
                page.td(str(), class_="calendar month")

        # other wise make a month table for it
        else:
            if ncol==1:
                page.td()
                page.table(class_="calendar")
            else:
                page.td(class_="calendar month")
                page.table(class_="calendar month")
            # get month as calendar
            month = calendar.monthcalendar(m.year, m.month)

            # get month link as table header
            page.tr()
            H = "%s %s" % (calendar.month_name[m.month], m.year)
            href = None
            for e,t in zip(["html", "php"], [dir1, dir2]):
                idx = os.path.join(jobdir, "archive_monthly",\
                                   "%s%s" % (m.strftime("%Y%m"), t),\
                                   "index.%s" % e)
                if os.path.isfile(os.path.expanduser(idx)):
                    href = "%s%s" % (os.path.split(idx)[0], os.path.sep)
                    break
            if href:
                page.th(markup.oneliner.a(H, class_="day", href=href),
                        colspan="100%")
            else:
                page.th(H, colspan="100%") 
            page.tr.close()

            # loop over weeks
            for week in month:
                page.tr()
                for idx,day in enumerate(week):
                    if day != 0:
                        break
                w = (datetime.date(m.year, m.month, day)\
                     - datetime.timedelta(days=idx)).strftime("%Y%m%d")
                href = None
                for e,t in [(e,t) for e in ["html", "php"]\
                                  for t in [dir1, dir2]]:
                    idx = os.path.join(jobdir, "archive_weekly",\
                                       "%s%s" % (w, t), "index.%s" % e)
                    if os.path.isfile(os.path.expanduser(idx)):
                        href = "%s%s" % (os.path.split(idx)[0], os.path.sep)
                        break
                if href:
                    page.td(markup.oneliner.a("w", class_="day", href=href))
                else:
                    page.td("w")

                for day in week:
                    # format gaps in tabls
                    if day == 0:
                        page.td("")
                        continue
                    # find day page and link
                    d = datetime.date(m.year, m.month, day).strftime("%Y%m%d")
                    href = None
                    for e,t in [(e,t) for e in ["html", "php"]\
                                      for t in [dir1, dir2]]:
                        idx = os.path.join(jobdir, "archive_daily",\
                                           "%s%s" % (d, t), "index.%s" % e)
                        if os.path.isfile(os.path.expanduser(idx)):
                            href = "%s%s" % (os.path.split(idx)[0], os.path.sep)
                            break
                    if href:
                        page.td(markup.oneliner.a(str(day), class_="day",\
                                                  href=href))
                    else:
                        page.td(str(day))
                page.tr.close()
            page.td.close()
            page.table.close()

        # close row
        if i % ncol == ncol -1:
             page.tr.close()

        # increment
        m += delta
        i += 1

    page.table.close()
    return page

def plotclass(n):
    """
    Guess the plot class to use from the number of plots
    """
    if n % 5 == 0:
        return "fifth"
    elif n % 4 == 0:
        return "quarter"
    elif n % 3 == 0:
        return "third"
    elif n % 2 == 0:
        return "half"
    else:
        return "full"
