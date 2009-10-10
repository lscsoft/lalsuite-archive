# Copyright (C) 2009  Kipp Cannon
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
This is a convenience module providing all the LAL datatype wrappings from
the pylal.xlal.datatypes subpackage in a single import.  It is recommended
that you import this module unless you require only exactly one or a few
specific types in which case you can import the individual modules if that
is easier.

Example:

>>> from pylal import datatypes as laltypes
>>> x = laltypes.REAL8TimeSeries()
"""


import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# remember to keep this list up to date
#


from xlal.datatypes.complex16fftplan import *
from xlal.datatypes.complex16frequencyseries import *
from xlal.datatypes.lalunit import *
from xlal.datatypes.ligotimegps import *
from xlal.datatypes.real8fftplan import *
from xlal.datatypes.real8frequencyseries import *
from xlal.datatypes.real8timeseries import *
from xlal.datatypes.real8window import *
