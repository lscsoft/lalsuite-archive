# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
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
High-level document manipulation utilities.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import ligolw


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#

import warnings
class PartialLIGOLWContentHandler(ligolw.PartialLIGOLWContentHandler):
	def __init__(self, document, filter):
		warnings.warn("docutils.PartialLIGOLWContentHandler() class is deprecated:  use ligolw.PartialLIGOLWContentHandler() instead", DeprecationWarning)
		ligolw.PartialLIGOLWContentHandler.__init__(self, document, filter)
