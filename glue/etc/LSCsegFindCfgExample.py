"""
The LSCsegFindCfgExample.py module contains example configuration settings for LSCsegFind
"""
__author__ = "Greg Mendell: Example configuration file for LSCsegFind"
__date__ = '$Date$'
__version__ = '$Revision$'[0:0]

# Note that this line is important; it defines cfgDefaults as a python dictionary:
cfgDefaults = {}

# Set defaults for each segment NAME to use with the LSCsegFind --type=NAME option here.

# Syntax is cfgDefaults['NAME'] = {'desc':'STRING','url':'URL','coalesce':BOOLEAN}
# cfgDefaults['NAME'] : replace NAME a name to use with the LSCsegFind --type=NAME option
# "desc":'STRING'     : replace STRING with a description of the segments
# "url":'URL'         : replace URL with the url to get segments from of this type.
# "coalesce":BOOLEAN  : replace BOOLEAN with True or False; True if segments need to be coalesced after retrieval from the URL

cfgDefaults['S3H1'] = { 'desc':'Conlog Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/scirun/S3/LockStatistics/S3segments/S3H1v00_segs.txt', 'coalesce':False }

cfgDefaults['S3H2'] = { 'desc':'Conlog Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/scirun/S3/LockStatistics/S3segments/S3H2v00_segs.txt', 'coalesce':False }

cfgDefaults['S3L1'] = { 'desc':'Conlog Science-mode Segments', 'url':'http://london.ligo-la.caltech.edu/scirun/S3/LockStatistics/S3segments/S3L1v00_segs.txt', 'coalesce':False }

cfgDefaults['M5H1'] = { 'desc':'Conlog Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/scirun/M5/LockStatistics/M5segments/M5H1v00_segs.txt', 'coalesce':False }

cfgDefaults['M5H2'] = { 'desc':'Conlog Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/scirun/M5/LockStatistics/M5segments/M5H2v00_segs.txt', 'coalesce':False }

cfgDefaults['H1DMT'] = { 'desc':'DMT Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/gds/monitor_reports/SegGener/H1-ScienceMode.txt', 'coalesce':True }

cfgDefaults['H2DMT'] = { 'desc':'DMT Science-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/gds/monitor_reports/SegGener/H2-ScienceMode.txt', 'coalesce':True }

cfgDefaults['H1DMT_CM'] = { 'desc':'DMT Common-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/gds/monitor_reports/SegGener/H1-CommonMode.txt', 'coalesce':True }

cfgDefaults['H2DMT_CM'] = { 'desc':'DMT Common-mode Segments', 'url':'http://blue.ligo-wa.caltech.edu/gds/monitor_reports/SegGener/H2-CommonMode.txt', 'coalesce':True }