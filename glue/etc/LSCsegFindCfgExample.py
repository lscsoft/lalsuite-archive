"""
The LSCsegFindCfgExample.py module contains configuration settings for LSCsegFind
"""
__author__ = "Greg Mendell: Example configuration file for LSCsegFind"
__date__ = '$Date$'
__version__ = '$Revision$'[0:0]

# Note that this line is important; it defines urlDefaults as a dictionary:
urlDefaults = {}
# Set default URLs for any IFO(s) here:
urlDefaults["H1"] = "http://blue.ligo-wa.caltech.edu/scirun/S3/LockStatistics/S3segments/S3H1v00_segs.txt"

urlDefaults["H2"] = "http://blue.ligo-wa.caltech.edu/scirun/S3/LockStatistics/S3segments/S3H2v00_segs.txt"

urlDefaults["L1"] = "http://london.ligo-la.caltech.edu/scirun/S3/LockStatistics/S3segments/S3L1v00_segs.txt"

