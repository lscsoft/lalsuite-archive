#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

'''
A collection of utilities to assist in printing out information from an xmldoc.
'''

import sys
import time
import datetime

from glue.ligolw.utils import print_tables
from glue.ligolw import ligolw
from glue.ligolw import table

from pylal.xlal.date import XLALGPSToUTC
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
try:
    XLALGPSToUTC(LIGOTimeGPS(987654321,0))
except TypeError:
    # s6 code
    from pylal.xlal.date import LIGOTimeGPS


__author__ = "Collin Capano <cdcapano@physics.syr.edu>"
__date__ = "$Date$" 
__version__ = "$Revision$"



# =============================================================================
#
#                           Utilities
#
# =============================================================================

def generic_get_pyvalue(obj):
    if obj.value is None:
        return None
    return ligolwtypes.ToPyType[obj.type or "lstring"](obj.value)


def get_columns_to_print(xmldoc, tableName):
    """
    Retrieves canonical columns to print for the given tableName.
    Returns a columnList, row_span and rspan_break lists.
    """
    tableName = tableName.endswith(":table") and tableName or tableName+":table"
    summTable = table.get_table(xmldoc, tableName )
    # get rankname
    rankname = [col.getAttribute("Name").split(":")[-1]
        for col in summTable.getElementsByTagName(u'Column') if "rank" in col.getAttribute("Name")][0]

    if tableName == "loudest_events:table":
        durname = [col.getAttribute("Name").split(":")[-1]
            for col in summTable.getElementsByTagName(u'Column') if "duration" in col.getAttribute("Name")][0]
        columnList = [
            rankname,
            'combined_far',
            'fap',
            'fap_1yr',
            'snr',
            'end_time',
            'end_time_ns',
            'end_time_utc__Px_click_for_daily_ihope_xP_',
            'mass',
            'mchirp',
            'ifos__Px_click_for_elog_xP_',
            'instruments_on',
            'mini_followup',
            durname]
        row_span_columns = rspan_break_columns = [durname]
    elif tableName == "selected_found_injections:table":
        durname = [col.getAttribute("Name").split(":")[-1]
            for col in summTable.getElementsByTagName(u'Column') if "duration" in col.getAttribute("Name")][0]
        columnList = [
            rankname,
            'mini_followup',
            'injected_end_time',
            'injected_end_time_ns',
            'injected_end_time_utc__Px_click_for_daily_ihope_xP_',
            'injected_eff_dist_h',
            'injected_eff_dist_l',
            'injected_eff_dist_v',
            'injected_mchirp',
            'injected_mass1',
            'injected_mass2',
            'recovered_match_rank',
            'recovered_ifos',
            'instruments_on__Px_click_for_elog_xP_',
            'recovered_combined_far',
            'recovered_fap',
            'recovered_fap_1yr',
            'recovered_snr',
            'recovered_end_time',
            'recovered_end_time_ns',
            'recovered_mchirp',
            'recovered_mass',
            durname]
        row_span_columns = rspan_break_columns = [
            rankname,
            'mini_followup',
            'injected_end_time',
            'injected_end_time_ns', 
            'injected_end_time_utc__Px_click_for_daily_ihope_xP_',
            'injected_eff_dist_h',
            'injected_eff_dist_l',
            'injected_eff_dist_v',
            'injected_mchirp',
            'injected_mass1',
            'injected_mass2']
    elif tableName == "close_missed_injections:table":
        columnList = [
            'rank',
            'decisive_distance',
            'end_time',
            'end_time_ns',
            'end_time_utc__Px_click_for_daily_ihope_xP_',
            'mchirp',
            'mass1',
            'mass2',
            'eff_dist_h',
            'eff_dist_l',
            'eff_dist_v',
            'instruments_on__Px_click_for_elog_xP_',
            'mini_followup'
            ]
        row_span_columns = rspan_break_columns = []
    else:
        # unrecognized table, just return all the columns in the table
        columnList = [col.getAttribute("Name").split(":")[-1] for col in summTable.getElementsByTagName(u'Column')]
        row_span_columns = rspan_break_columns = []
        
    return columnList, row_span_columns, rspan_break_columns


#
#   Some helper functions for manipulating times
#

def get_dst_start_end(ifo, year):
    """
    Figures out what dates daylight savings time starts and ends at a given site on a given year.
    """
    # in the United Stats, prior to 2007, DST began on the first Sunday in April
    # and ended on the last Sunday in October (http://aa.usno.navy.mil/faq/docs/daylight_time.php)
    if ("H" in ifo  or "L" in ifo) and year < 2007:
        for ii in range(1,28):
            dst_start = datetime.datetime(year, 4, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday':
                break
        for ii in range(31,0,-1):
            dst_end = datetime.datetime(year, 10, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    # in the US, starting in 2007, DST begins on the second Sunday in March and ends on the first
    # Sunday in November
    elif ("H" in ifo  or "L" in ifo) and year >= 2007:
        nn = 1
        for ii in range(1,31):
            dst_start = datetime.datetime(year, 3, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday' and nn == 2:
                break
            elif dst_start.strftime('%A') == 'Sunday':
                nn += 1
        for ii in range(1,28):
            dst_end = datetime.datetime(year, 11, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    # in Europe, DST begins on the last Sunday of March and ends on the last Sunday of October
    # source: http://www.timeanddate.com/news/time/europe-dst-starts-march-29-2009.html
    elif ("V" in ifo or "G" in ifo):
        for ii in range(31,0,-1):
            dst_start = datetime.datetime(year, 3, ii, 2, 0, 0)
            if dst_start.strftime('%A') == 'Sunday':
                break
        for ii in range(31,0,-1):
            dst_end = datetime.datetime(year, 10, ii, 2, 0, 0)
            if dst_end.strftime('%A') == 'Sunday':
                break
    else:
        raise ValueError, "unrecognized ifo %s" % ifo
    
    return dst_start, dst_end
        

def get_sitelocaltime_from_gps(ifo, gpstime):
    # get the utc time in datetime.datetime format
    utctime = XLALGPSToUTC(LIGOTimeGPS(gpstime, 0))
    utctime = datetime.datetime(utctime[0],utctime[1],utctime[2],utctime[3],utctime[4],utctime[5],utctime[6])
    # figure out if daylight savings time was on or not
    dst_start, dst_end = get_dst_start_end(ifo, utctime.year)
    # figure out the appropriate time offset
    if "H" in ifo:
        toffset = datetime.timedelta(hours=-7)
    elif "L" in ifo:
        toffset = datetime.timedelta(hours=-5)
    elif ("V" in ifo or "G" in ifo):
        toffset = datetime.timedelta(hours=+2)
    # apply the dst time offset to see if daylight savings was on; if not, adjust the toffset
    if not (utctime + toffset >= dst_start and utctime + toffset < dst_end):
        toffset = toffset + datetime.timedelta(hours=-1)

    return utctime + toffset


def format_end_time_in_utc(gps_sec):
    return time.strftime("%a %d %b %Y %H:%M:%S", XLALGPSToUTC(LIGOTimeGPS(gps_sec, 0)))


def get_elog_page(ifo, gpstime):
    # set site_address
    if "H" in ifo:
        site_address = "http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?group=detector"
    elif "L" in ifo:
        site_address = "http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?group=detector"
    elif "V" in ifo:
        #FIXME: What's the site address and format for Virgo log book?
        site_address = "https://pub3.ego-gw.it/logbook/"
    # get local time at the site
    site_localtime = get_sitelocaltime_from_gps(ifo, gpstime)
    # set the address
    if "H" in ifo or "L" in ifo:
        site_address = "%s&date_to_view=%s" % ( site_address, site_localtime.strftime("%m/%d/%Y") )

    return site_address

def get_daily_ihope_page(gpstime, pages_location = "https://ldas-jobs.ligo.caltech.edu/~cbc/ihope_daily"):
    utctime = XLALGPSToUTC(LIGOTimeGPS(gpstime, 0))
    return "%s/%s/%s/" %(pages_location, time.strftime("%Y%m", utctime), time.strftime("%Y%m%d", utctime))


def create_hyperlink(address, link):
    return '<a href="%s">%s</a>' % (address, link)
