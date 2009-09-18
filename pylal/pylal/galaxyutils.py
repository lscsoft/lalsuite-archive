#
# Copyright (C) 2007  Nickolas Fotopoulos
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

from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import math
itertools = __import__("itertools")  # system-wide itertools

import numpy

##############################################################################
# math
##############################################################################

def _make_within_pi(angles):
    """
    Make an array of angles completely between -pi and pi.
    """
    return angles - numpy.sign(angles) * \
        numpy.round(abs(angles) / (2 * numpy.pi)) * 2 * numpy.pi

def is_inside_polygon(point, vertices):
    """
    Return True if the (2-D) point is inside the (2-D) polygon defined by the
    vertices.
  
    Adapted from:
    http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/ (solution 2)
    """
    point = numpy.array(point)
    centered_vertices = numpy.empty(shape=(len(vertices) + 1, 2), dtype=float)
    centered_vertices[:-1] = vertices - point
    centered_vertices[-1] = centered_vertices[0]
    angles = numpy.arctan2(centered_vertices[:, 1], centered_vertices[:, 0])
    diff_angles = _make_within_pi(angles[1:] - angles[:-1])
    angle_sum = diff_angles.sum()
    return abs(angle_sum) >= numpy.pi

def is_within_distances(gal_dist, dmin, dmax):
    """
    Return True if dmax > gal_dist > dmin
    """
    if (dmin > dmax):
        raise ValueError, "minimum distance is greater than maximum distance " + str(dmin) + " > "+ str(dmax)
    if (0 > dmin) or (0 > dmax):
        raise ValueError, "negative distance " + str(dmin) + " " + str(dmax)
    return (dmax > gal_dist) and (dmin < gal_dist)

def plot_points(points, polygon):
    """
    Return the handle to a figure containing a scatter plot of a set of
    points, each colored by is_inside_polygon's determination of whether it is
    inside or outside the given polygon.
    
    This is an excellent way to test is_inside_polygon.
    """
    found_list = [is_inside_polygon(point, polygon) for point in points]
    found_points = numpy.array([point for point,found in \
        itertools.izip(points, found_list) if found])
    missed_points = numpy.array([point for point,found in \
        itertools.izip(points, found_list) if not found])
    
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.plot(found_points[:, 0], found_points[:, 1], "bx", label="found")
    ax.plot(missed_points[:, 0], missed_points[:, 1], "ro", label="missed")
    
    closed_poly = numpy.vstack((polygon, polygon[0]))
    ax.plot(closed_poly[:, 0], closed_poly[:, 1], "k-", label="_nolabel_")
    return fig

# visual test code for is_inside_polygon

# import pylab
# points = numpy.random.uniform(low=0.0, high=5.0, size=(5000,2))
# test_polygon = numpy.array([[2, 2], [2, 3], [3, 2]], dtype=float)
# plot_points(points, test_polygon)
# pylab.show()
# sys.exit()

##############################################################################
# unit conversion
##############################################################################

def hms2rad(ra_sex):
    """
    Convert right ascension and from a sexagesimal string to floating point
    radians.  Handles h:m:s or h:m.
    """
    tup = ra_sex.split(":")
    h = int(tup[0])
    m = int(tup[1])
    s = float(tup[2])
    
    if (h < 0 or h > 23) or (m < 0 or m > 60) or (s < 0 or s >= 60):
        raise ValueError, "hour, minute, or second out of bounds " + ra_sex
    return 2 * numpy.pi * (h + (m + s / 60) / 60) / 24

def dms2rad(dec_sex):
    """
    Convert declination from a colon-delimited sexagesimal string to floating
    point radians.  Handles d:m:s.
    """
    tup = dec_sex.split(":")
    if tup[0].startswith("-"):
        d = int(tup[0][1:])
        sign = -1
    else:
        d = int(tup[0])
        sign = +1
    m = int(tup[1])
    s = float(tup[2])
  
    if (d > 89) or (m < 0 or m > 60) or (s < 0 or s >= 60):
        raise ValueError, "degree, minute, or second out of bounds: " + dec_sex
    return numpy.pi * sign * (d + (m + s / 60) / 60) / 180

def hm2rad(ra_sex):
    """
    Convert right ascension and from a sexagesimal string to floating point
    radians.  Handles h:m.
    """
    tup = ra_sex.split(":")
    h = int(tup[0])
    m = float(tup[1])
  
    if (h < 0 or h > 23) or (m < 0 or m > 60):
        raise ValueError, "hour or minute out of bounds " + ra_sex
    return 2 * numpy.pi * (h + m / 60) / 24

def dm2rad(dec_sex):
  """
  Convert declination from a colon-delimited sexagesimal string to floating
  point radians.  Handles d:m.
  """
  tup = dec_sex.split(":")
  if tup[0].startswith("-"):
      d = int(tup[0][1:])
      sign = -1
  else:
      d = int(tup[0])
      sign = 1
  m = float(tup[1])

  
  if (d > 89) or (m < 0 or m > 60):
    raise ValueError, "degree or minute out of bounds: " + dec_sex
  return numpy.pi * sign * (d + m / 60) / 180

##############################################################################
# galaxy and galaxy catalog representations
##############################################################################

def column_cmp(k1, k2):
    """
    Comparison function to sort column names by their column index.
    """
    return cmp(GalaxyCatalog.valid_columns[k1][0], GalaxyCatalog.valid_columns[k2][0])

class GalaxyCatalog(list):
    valid_columns =  {
        "name": (0, str),
        "ra": (1, hm2rad),
        "dec": (2, dm2rad),
        # distance measured in kpc
        "distance_kpc": (3, float),
        # luminosity measured in milky way equivalent galaxies
        "luminosity_mwe": (4, float),
        "metal_correction": (5, float),
        "magnitude_error": (6, float),
        "distance_error": (7, float),
        }
    
    def entry_from_line(cls, line, load_columns):
        # create blank entry
        row = cls.entry_class()
        
        # parse line
        tup = line.split()
        
        # fill the entry
        for col_name in load_columns:
            col_index, col_type = cls.valid_columns[col_name]
            setattr(row, col_name, col_type(tup[col_index]))
        
        return row
    entry_from_line = classmethod(entry_from_line)
    
    def from_file(cls, fileobj, load_columns=None):
        # set/validate columns to load
        if load_columns is None:
            load_columns = cls.valid_columns.keys()
            load_columns.sort(column_cmp)
        else:
            for col in load_columns:
                if col not in cls.valid_columns:
                    raise ValueError, "no such column exists"
        cls.entry_class.__slots__ = load_columns
        
        # load them
        return cls([cls.entry_from_line(line, load_columns) for line \
                    in fileobj if not line.startswith("#")])
    from_file = classmethod(from_file)

    def within_polygon(self, vertices):
        return self.__class__([gal for gal in self if is_inside_polygon(gal.coords, vertices)])

    def within_distances(self, dmin, dmax):
        return self.__class__([gal for gal in self if is_within_distances(gal.distance_kpc, dmin, dmax)])

    def __repr__(self):
        return "\n".join(itertools.imap(str, self))

class Galaxy(object):
    """
    A galaxy object that knows how to initialize itself from a line in a text
    file and consumes a minimum of memory.
    """
    __slots__ = GalaxyCatalog.valid_columns.keys()
    __slots__.sort(column_cmp)

    def __str__(self):
        return "\t".join([str(getattr(self, slot)) for slot in self.__slots__])

    def __repr__(self):
        return "Galaxy(\"" + str(self) + "\")"

    def _coords_getter(self):
        return (self.ra, self.dec)
    coords = property(fget=_coords_getter)

GalaxyCatalog.entry_class = Galaxy
