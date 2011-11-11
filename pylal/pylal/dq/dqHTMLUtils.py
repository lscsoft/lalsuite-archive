#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re
from glue import markup

from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a few extensions to glue.markup to streamline GEO/LIGO detector characterisation tools that use very similar web interfaces
"""

# =============================================================================
# Write table
# =============================================================================

def write_table(page, headers, data, cl=''):

  """
    Write table into glue.markup.page object. headers are written with <th>,
    multiple columns of data are written with <td>. Classes include:
      * "", default class writes standard table with single column of headers
            and multiple columns of data
      * "list", writes table of 'header[i]: data[i]' definition-style entries

    Arguments:

      page : glue.markup.page
        page object into which to write table
      headers : list
        list of table header elements
      data : list
        list (or nested list) of table data elements, list of lists used for
        multiple rows
 
    Keyword arguments:

      cl : string
        name for HTML table class, cl='list' treats special case above

  """

  # open table
  page.table(class_=cl)

  # list: print two vertical columns of header:data pairs
  if cl=='list':
    for i in range(len(headers)):

      page.tr()
      page.th()
      page.add(str(headers[i]))
      page.th.close()
      page.td()
      page.add(str(data[i]))
      page.td.close()
      page.tr.close()

  # otherwise print 'standard' table with single header row and multiple data
  # rows
  else:
    page.tr()
    if len(headers)==1:
      page.th(colspan="100%")
      page.add(str(headers[0]))
      page.th.close()
    else:
      for n in headers:
        page.th()
        page.add(str(n))
        page.th.close()
    page.tr.close()

    if data and not re.search('list',str(type(data[0]))):
      data = [data]

    for row in data:
      page.tr()
      for item in row:
        page.td()
        page.add(str(item))
        page.td.close()
      page.tr.close()

  page.table.close()

  return page

# =============================================================================
# Write <head>
# =============================================================================

def write_head(title, css, js, base=None, refresh=None, jquery=True):

  """
    Returns glue.markup.page object with <head> tag filled.

    Arguments:

      title : string
        text for <title> tag
      css : string
        relative path to style sheet
      js : string
        relative path to javascript

    Keyword arguments:

      base : string
        absolute http(s) path of url base
      refresh : int
        number of seconds after which to refresh page automatically
      jquery : [ True | False ]
        import jquery AJAX script in header, default: True
  """

  page = markup.page(mode="strict_html")
  page._escape = False

  doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
  doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""
  page.add(doctype)
  page.head()
  if base:
    page.base(href=base)
  if refresh:
    page.add('<meta http-equiv="refresh" content="%s">' % refresh)
  page.link(media="all", href=css, type="text/css", rel="stylesheet")
  page.title()
  page.add(title)
  page.title.close()

  if jquery:
    page.script(src="http://ajax.googleapis.com/ajax/libs/jquery/1.2.6"\
                    "/jquery.min.js", type="text/javascript")
    page.script.close()
  page.script(src=js, type="text/javascript")
  page.script.close()
  page.head.close()

  return page

# =============================================================================
# Write <div id="header">
# =============================================================================

def write_banner(title, text='&nbsp;'):

  """
    Returns glue.markup.page object for <div id="header">
  """

  page = markup.page(mode="strict_html")
  page._escape = False

  page.div(id="header")
  page.h1()
  page.add(title)
  page.h1.close()
  page.h3()
  page.add(text)
  page.h3.close()

  page.hr(class_="short")
  page.hr(class_="long")

  page.div.close()

  return page

# =============================================================================
# Write <div id="menubar">
# =============================================================================

def write_menu(sections, pages, current=None):

  """
    Returns glue.markup.page object for <div id="menubar">, constructing menu
    in HTML.

    Arguments:

      sections : list
        ordered list of menu entry names
      pages : dict
        dict of section:href pairs holding link paths for each element of
        sections list
  """

  page = markup.page(mode="strict_html")
  page._escape = False

  page.div(id="menubar")

  for i,sec in enumerate(sections):
    if sec==current:
      cl = "menulink selected"
    else:
      cl = "menulink"
    page.a(id="link_%s" % i, class_=cl, href="%s" % pages[sec])
    page.add(sec)
    page.a.close()

  page.script(type="text/javascript")
  page.script.close()

  page.div.close()

  return page

# =============================================================================
# Initialise page
# =============================================================================

def init_page(head, banner, menu):

  """
    Initialise html into markup page, including <head> tag, banner and menu.
  """

  # write html
  page = markup.page()
  page._escape = False

  # initialise page
  page.html(lang="en")
  page.add(head())

  # open body
  page.body()

  # open container for page (needed to position footer)
  page.div(id="container")
  # add banner
  page.add(banner())
  # open content (tab below banner and above footer)
  page.div(id="content")
  # print menu
  page.add(menu())
  # initialise maintab
  page.div(id="maintab")

  return page

# =============================================================================
# Close page
# =============================================================================

def close_page(page, footer=False):

  """
    Close content, maintab, and container divs, write footer and close
    <body> and <html> tags.
  """

  # close maintab
  page.div.close()
  # close content tab
  page.div.close()
  # close container
  page.div.close()
  # add footer
  if footer:
    page.add(foot())

  # close page
  page.body.close()
  page.html.close()

  return page

# =============================================================================
# Write glossary
# =============================================================================

def write_glossary(page, terms):

  """
    Write a glossary of DQ terms into the glue.markup.page object page using the
    list of (term,definition) tuples terms.
  """

  page.h2()
  page.add('Glossary')
  page.h2.close()
  page.add('This section gives a glossary of terms used in this analysis')

  i=1

  # write section
  page = write_h(page, 'Glossary', [i], cl=3)
  page.div(id="div_%s" % (i), style="display: none;")

  page.add('The LIGO-Virgo acronym wiki can be found on')
  page.a(href="https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/Acronyms")
  page.add('this page')
  page.a.close()

  # write glossary table
  tmp = {}
  for (term, text) in terms:
    tmp[term.title()] = text
  th = sorted(tmp.keys())
  td = [tmp[t] for t in th]
  page = write_table(page, th, td, 'list')

  page.div.close()

  return page

# =============================================================================
# Write <hX>
# =============================================================================

def write_h(page, title, idx, cl=3):

  """
    Write hX header into glue.markup.page object page with toggling. Text
    contained within title is printed while link is constructed and toggled
    using javascipt toggleVisible function
  """

  if not isinstance(idx, list):
    idx = [idx]

  idx = map(str, idx)

  page.input(id="input_%s" % '.'.join(idx), type="button", class_="h%s" % cl,\
             value=title, onclick="toggleVisible('%s')" % '.'.join(idx))

  return page

# =============================================================================
# Link image
# =============================================================================

def link_image(page, href, src, alt, title, class_=None):

  """
    Link image into glue.markup.page object page with standard options
  """

  page.a(href=href, title=title, rel="external")
  page.img(class_=class_, src=src, alt=alt)
  page.a.close()

  return page

# =============================================================================
# Link file
# =============================================================================

def link_file(page, href, text):

  """
    Link file into glue.markup.page object page, with associated text and
    options.
  """

  page.a(href="%s" % href, rel="external")
  page.add(text)
  page.a.close()

  return page

