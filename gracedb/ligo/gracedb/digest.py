#!/usr/bin/env python
#
# Copyright (C) 2012  Leo Singer <leo.singer@ligo.org>
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
"""
Prepare a summary or digest of a gravitational wave candidate in GraCEDb,
and provide functions to write it in different formats, including FITS HEALPix
image, VOEvent, JSON, and a 'mad libs' format suitable for the body of a
GCN circular.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Configure logging.
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger("digest")


def _download_to_strio(client, graceid, filename):
    """Download a file from GraCEDb into an in-memory string, wrapped in a
    file-like object. Raise an IOError if anything goes wrong with the data
    transfer."""
    from cStringIO import StringIO
    strio = StringIO()
    response = client.download(graceid, filename, strio)
    if response:
        raise IOError(response)
    return strio


def get_digest(graceid, *args, **kwargs):
    """Look up an event in GraCEDb by its event id and return a digest as dictionary."""
    from ligo import gracedb
    from cStringIO import StringIO
    import glue.ligolw.utils as ligolw_utils
    import glue.ligolw.table as ligolw_table
    import glue.ligolw.lsctables as ligolw_lsctables
    import healpy

    client = gracedb.Client(*args, **kwargs)
    log.info("started GraCEDb client")

    # Read the coinc.xml file.
    filename = "coinc.xml"
    try:
        strio = _download_to_strio(client, graceid, filename)
    except IOError:
        log.exception("%s:could not retrieve %s", graceid, filename)
        raise
    log.info("%s:retrieved %s", graceid, filename)

    # Now rewind to the beginning of the file.
    strio.seek(0)

    # Parse the file as LIGO-LW XML.
    xmldoc, digest = ligolw_utils.load_fileobj(strio)

    # Retrieve the tables that we need.
    process_table = ligolw_table.get_table(xmldoc, ligolw_lsctables.ProcessTable.tableName)
    coinc_inspiral_table = ligolw_table.get_table(xmldoc, ligolw_lsctables.CoincInspiralTable.tableName)
    sngl_inspiral_table = ligolw_table.get_table(xmldoc, ligolw_lsctables.SnglInspiralTable.tableName)

    # Retrieve values from the table that we need.
    snr = coinc_inspiral_table[0].snr
    false_alarm_rate = coinc_inspiral_table[0].false_alarm_rate
    geocent_end_time = coinc_inspiral_table[0].end_time + 1e-9 * coinc_inspiral_table[0].end_time_ns
    mchirp = coinc_inspiral_table[0].mchirp
    mtotal = coinc_inspiral_table[0].mass

    # Construct digest.
    digest = {
        "graceid": graceid,
        "detectors": [],
        "snr": snr,
        "false_alarm_rate": false_alarm_rate,
        "geocent_end_time": geocent_end_time,
        "mchirp": mchirp,
        "mtotal": mtotal
    }

    # Now try to retrieve skymap.
    filename = "general/bayestar/skymap.fits"
    try:
        strio = _download_to_strio(client, graceid, filename)

        # Now rewind to the beginning of the file.
        strio.seek(0)
    except IOError:
        log.exception("%s:could not retrieve %s", graceid, filename)
    else:
        log.info("%s:retrieved %s", graceid, filename)
        digest["skymap"] = healpy.read_map(strio)
        log.info("extracted HEALPix map")

    return digest


def put_digest_fits(fileobj, digest):
    """Format the digest as a FITS file and write to fileobj."""
    # FIXME: This just writes a FITS header, not a whole FITS file.
    # For the sky map, do we want to use the 'HEALPix' format or
    # the WCS 'HPX' format? Find out which is more commonly used.

    from astropy.io import fits
    import lal
    import datetime

    date_obs = datetime.datetime(*lal.GPSToUTC(long(digest["geocent_end_time"]))[:-2]).isoformat()
    date = datetime.datetime.now().isoformat()

    header = fits.Header()
    header.append(("OBS_ID", digest["graceid"], "unique observation ID"))
    header.append(("SNR", digest["snr"], "network signal to noise ratio (dimensionless)"))
    header.append(("FAR", digest["false_alarm_rate"], "false alarm rate (Hz)")) # FIXME: what units?
    header.append(("DATE-OBS", date_obs, "date of the observation"))
    header.append(("DATE", date, "date of file creation"))
    header.append(("MCHIRP", digest["mchirp"], "chirp mass (solar masses)"))
    header.append(("MTOTAL", digest["mtotal"], "total mass (solar masses)"))

    log.info("writing fits")
    print >>fileobj, str(header).strip()


def put_digest_voevent(fileobj, digest):
    """Format the digest as a VOEvent and write to fileobj."""
    log.info("writing voevent")
    raise NotImplementedError


def put_digest_madlibs(fileobj, digest):
    """Format the digest in a human-readable form with commentary, suitable for
    the content of a GCN circular, and write to fileobj."""
    log.info("writing madlib")
    raise NotImplementedError


def put_digest_json(fileobj, digest):
    """Format the digest as JSON and write to fileobj."""
    import json
    log.info("writing json")
    json.dump(dict(digest), fileobj)


if __name__ == "__main__":
    import optparse
    import sys

    # Dictionary that maps command line arguments -> functions.
    putters = {
        "fits": put_digest_fits,
        "voevent": put_digest_voevent,
        "madlibs": put_digest_madlibs,
        "json": put_digest_json
    }

    default_putter = "fits"

    # Command line interface.
    import optparse
    parser = optparse.OptionParser(description=__doc__,
        usage="%prog [options] [-o OUTPUT] GRACEID")
    parser.add_option("-o", "--output", metavar="FILE", help="Output file (default=stdout)")
    parser.add_option("-f", "--format", metavar="|".join(putters.keys()),
        choices=putters.keys(), default=default_putter,
        help="Digest format (default=%s)" % default_putter)

    opts, args = parser.parse_args()

    if len(args) > 1:
        parser.error("too many command line arguments")
    elif len(args) == 0:
        parser.error("not enough command line arguments")
    else:
        graceid = args[0]

    digest = get_digest(graceid)
    putter = putters[opts.format]

    if opts.output is None:
        fileobj = sys.stdout
    else:
        fileobj = open(opts.output, "w")

    putter(fileobj, digest)
