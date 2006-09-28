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
LSC Table definitions.  These must be kept synchronized with the official
definitions in the LDAS CVS repository at
http://www.ldas-sw.ligo.caltech.edu/cgi-bin/cvsweb.cgi/ldas/dbms/db2/sql.
Maintainership of the table definitions is left as an excercise to
interested users.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

from xml import sax

from glue import lal
from glue import segments
import ligolw
import table
import types
import ilwd


#
# =============================================================================
#
#                            Convenience Functions
#
# =============================================================================
#

def New(Type, columns = None):
	"""
	Convenience function for constructing pre-defined LSC tables.  The
	optional columns argument is a list of the names of the columns the
	table should be constructed with.  If columns = None, then the
	table is constructed with all valid columns included (pass columns
	= [] to create a table with no columns).

	Example:

	>>> import lsctables
	>>> new = lsctables.New(lsctables.ProcessTable)
	"""
	new = Type(sax.xmlreader.AttributesImpl({u"Name": Type.tableName}))
	colnamefmt = ":".join(Type.tableName.split(":")[:-1]) + ":%s"
	if columns is not None:
		for key in columns:
			if key not in new.validcolumns:
				raise ligolw.ElementError, "New(): invalid Column '%s' for Table '%s'" % (key, new.tableName)
			new.appendChild(table.Column(sax.xmlreader.AttributesImpl({u"Name": colnamefmt % key, u"Type": new.validcolumns[key]})))
	else:
		for key, value in new.validcolumns.items():
			new.appendChild(table.Column(sax.xmlreader.AttributesImpl({u"Name": colnamefmt % key, u"Type": value})))
	new._end_of_columns()
	new.appendChild(table.TableStream(sax.xmlreader.AttributesImpl({u"Name": Type.tableName})))
	return new


def IsTableElement(Type, elem):
	"""
	Convenience function to check that an element is a Table of type
	Type.
	"""
	if elem.tagName != ligolw.Table.tagName:
		return False
	return table.CompareTableNames(elem.getAttribute("Name"), Type.tableName) == 0


def IsTableProperties(Type, tagname, attrs):
	"""
	Convenience function to check that the given tag name and
	attributes match those of a Table of type Type.
	"""
	if tagname != ligolw.Table.tagName:
		return False
	return table.CompareTableNames(attrs["Name"], Type.tableName) == 0


def getTablesByType(elem, Type):
	"""
	Return a list of tables of type Type under elem.
	"""
	return table.getTablesByName(elem, Type.tableName)


def HasNonLSCTables(elem):
	"""
	Return True if the document tree below elem contains non-LSC
	tables, otherwise return False.
	"""
	for t in elem.getElementsByTagName(ligolw.Table.tagName):
		if table.StripTableName(t.getAttribute("Name")) not in TableByName:
			return True
	return False


#
# =============================================================================
#
#                     Table Class With a Mapping Protocol
#
# =============================================================================
#

class LSCTableRowDict(object):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows when multiple rows share the same key.
	"""
	def __init__(self, table_elem):
		"""
		Initialize the mapping on the list of rows.
		"""
		self._table = table_elem
		self._keys = table_elem.getColumnByName(table_elem.ids.column_name)

	def __len__(self):
		"""
		Return the number of unique IDs.
		"""
		return len(self.keys())

	def __getitem__(self, key):
		"""
		Return a list of the rows whose ID equals key.
		"""
		l = [self._table[i] for i, k in enumerate(self._keys) if k == key]
		if not len(l):
			raise KeyError, key
		return l

	def __setitem__(self, key, values):
		"""
		Replace the rows whose ID equals key with the rows in the
		list of values, appending to the table if there are no rows
		with that ID.
		"""
		# FIXME: should we assign the key to rows?
		del self[key]
		map(self._table.append, values)

	def __delitem__(self, key):
		"""
		Delete all the rows whose ID equals key.
		"""
		for i in xrange(len(self._keys), -1, -1):
			if self._keys[i] == key:
				del self._table[i]

	def __iter__(self):
		"""
		Iterate over the unique IDs.
		"""
		return {}.fromkeys(self._keys).iterkeys()

	iterkeys = __iter__

	def __contains__(self, key):
		"""
		Return True if the table contains a row whose ID equals
		key, otherwise return False.
		"""
		return key in self._keys

	has_key = __contains__

	def keys(self):
		"""
		Return a list of the unique IDs.
		"""
		return {}.fromkeys(self._keys).keys()

	def iteritems(self):
		"""
		Iterate over (key, value) pairs.
		"""
		for key in self:
			yield key, self[key]

	def itervalues(self):
		"""
		Return an iterator over rows.
		"""
		for key in self:
			yield self[key]


class LSCTable(table.Table):
	"""
	A table with a mapping protocol for rows.
	"""
	def _end_of_columns(self):
		table.Table._end_of_columns(self)
		self.dict = LSCTableRowDict(self)


#
# =============================================================================
#
#                                process:table
#
# =============================================================================
#

class ProcessIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "process", "process_id", n)


class ProcessTable(LSCTable):
	tableName = "process:table"
	validcolumns = {
		"program": "lstring",
		"version": "lstring",
		"cvs_repository": "lstring",
		"cvs_entry_time": "int_4s",
		"comment": "lstring",
		"is_online": "int_4s",
		"node": "lstring",
		"username": "lstring",
		"unix_procid": "int_4s",
		"start_time": "int_4s",
		"end_time": "int_4s",
		"jobid": "int_4s",
		"domain": "lstring",
		"ifos": "lstring",
		"process_id": "ilwd:char"
	}
	ids = ProcessIDs()

	def get_ids_by_program(self, program):
		"""
		Return a sorted list of the process IDs for rows whose
		program string equals the given program.
		"""
		ids = [row.process_id for row in self if row.program == program]
		ids.sort()
		return ids


class Process(object):
	__slots__ = ProcessTable.validcolumns.keys()


ProcessTable.RowType = Process


#
# =============================================================================
#
#                                lfn:table
#
# =============================================================================
#

class LfnIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "lfn", "lfn_id", n)


class LfnTable(LSCTable):
	tableName = "lfn:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"lfn_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"start_time": "int_4s",
		"end_time": "int_4s",
	}
	ids = LfnIDs()


class Lfn(object):
	__slots__ = LfnTable.validcolumns.keys()

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in LfnTable.validcolumns.keys():
			if key == "lfn_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0


LfnTable.RowType = Lfn


#
# =============================================================================
#
#                             process_params:table
#
# =============================================================================
#

class ProcessParamsTable(table.Table):
	tableName = "process_params:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"param": "lstring",
		"type": "lstring",
		"value": "lstring"
	}

	def append(self, row):
		if row.type not in types.Types:
			raise ligolw.ElementError, "ProcessParamsTable.append(): unrecognized type '%s'" % row.type
		table.Table.append(self, row)


class ProcessParams(object):
	__slots__ = ProcessParamsTable.validcolumns.keys()


ProcessParamsTable.RowType = ProcessParams


#
# =============================================================================
#
#                             search_summary:table
#
# =============================================================================
#

class SearchSummaryTable(table.Table):
	tableName = "search_summary:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"shared_object": "lstring",
		"lalwrapper_cvs_tag": "lstring",
		"lal_cvs_tag": "lstring",
		"comment": "lstring",
		"ifos": "lstring",
		"in_start_time": "int_4s",
		"in_start_time_ns": "int_4s",
		"in_end_time": "int_4s",
		"in_end_time_ns": "int_4s",
		"out_start_time": "int_4s",
		"out_start_time_ns": "int_4s",
		"out_end_time": "int_4s",
		"out_end_time_ns": "int_4s",
		"nevents": "int_4s",
		"nnodes": "int_4s"
	}

	def get_inlist(self):
		"""
		Return a segmentlist object describing the times spanned by
		the input segments of all rows in the table.
		"""
		return segments.segmentlist([row.get_in() for row in self])

	def get_outlist(self):
		"""
		Return a segmentlist object describing the times spanned by
		the output segments of all rows in the table.
		"""
		return segments.segmentlist([row.get_out() for row in self])

	def get_inprocs(self, seglist):
		"""
		Return a list of the process IDs for the processes whose
		input segments intersect some part of the segmentlist
		seglist.
		"""
		return [row.process_id for row in self if segments.segmentlist([row.get_in()]) & seglist]

	def get_outprocs(self, seglist):
		"""
		Return a list of the process IDs for the processes whose
		output segments intersect some part of the segmentlist
		seglist.
		"""
		return [row.process_id for row in self if segments.segmentlist([row.get_out()]) & seglist]

	def get_out_segmentlistdict(self, process_ids = None):
		"""
		Return a segmentlistdict mapping instrument to out segment
		list.  If process_ids is a list of process IDs, then only
		rows with matching IDs are included otherwise all rows are
		included.
		"""
		seglistdict = segments.segmentlistdict()
		for row in self:
			if process_ids is None or row.process_id in process_ids:
				for ifo in row.ifos.split(","):
					if ifo in seglistdict:
						seglistdict[ifo].append(row.get_out())
					else:
						seglistdict[ifo] = segments.segmentlist([row.get_out()])
		return seglistdict.coalesce()


class SearchSummary(object):
	__slots__ = SearchSummaryTable.validcolumns.keys()

	def get_in(self):
		"""
		Return the input segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.in_start_time, self.in_start_time_ns), lal.LIGOTimeGPS(self.in_end_time, self.in_end_time_ns))

	def set_in(self, seg):
		"""
		Set the input segment.
		"""
		self.in_start_time, self.in_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.in_end_time, self.in_end_time_ns = seg[1].seconds, seg[1].nanoseconds

	def get_out(self):
		"""
		Get the output segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.out_start_time, self.out_start_time_ns), lal.LIGOTimeGPS(self.out_end_time, self.out_end_time_ns))

	def set_out(self, seg):
		"""
		Set the output segment.
		"""
		self.out_start_time, self.out_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.out_end_time, self.out_end_time_ns = seg[1].seconds, seg[1].nanoseconds


SearchSummaryTable.RowType = SearchSummary


#
# =============================================================================
#
#                            search_summvars:table
#
# =============================================================================
#

class SearchSummVarsTable(table.Table):
	tableName = "search_summvars:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"name": "lstring",
		"string": "lstring",
		"value": "real_8"
	}


class SearchSummVars(object):
	__slots__ = SearchSummVarsTable.validcolumns.keys()


SearchSummVarsTable.RowType = SearchSummVars


#
# =============================================================================
#
#                               sngl_burst:table
#
# =============================================================================
#

class SnglBurstIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sngl_burst", "event_id", n)


class SnglBurstTable(LSCTable):
	tableName = "sngl_burst:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"filter_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"stop_time": "int_4s",
		"stop_time_ns": "int_4s",
		"duration": "real_4",
		"flow": "real_4",
		"fhigh": "real_4",
		"central_freq": "real_4",
		"bandwidth": "real_4",
		"amplitude": "real_4",
		"snr": "real_4",
		"confidence": "real_4",
		"tfvolume": "real_4",
		"hrss": "real_4",
		"time_lag": "real_4",
		"peak_time": "int_4s",
		"peak_time_ns": "int_4s",
		"peak_frequency": "real_4",
		"peak_strain": "real_4",
		"peak_time_error": "real_4",
		"peak_frequency_error": "real_4",
		"peak_strain_error": "real_4",
		"ms_start_time": "int_4s",
		"ms_start_time_ns": "int_4s",
		"ms_stop_time": "int_4s",
		"ms_stop_time_ns": "int_4s",
		"ms_duration": "real_4",
		"ms_flow": "real_4",
		"ms_fhigh": "real_4",
		"ms_bandwidth": "real_4",
		"ms_hrss": "real_4",
		"ms_snr": "real_4",
		"ms_confidence": "real_4",
		"param_one_name": "lstring",
		"param_one_value": "real_8",
		"param_two_name": "lstring",
		"param_two_value": "real_8",
		"param_three_name": "lstring",
		"param_three_value": "real_8",
		"event_id": "ilwd:char"
	}
	ids = SnglBurstIDs()


class SnglBurst(object):
	__slots__ = SnglBurstTable.validcolumns.keys()
	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def get_stop(self):
		return lal.LIGOTimeGPS(self.stop_time, self.stop_time_ns)

	def set_stop(self, gps):
		self.stop_time, self.stop_time_ns = gps.seconds, gps.nanoseconds

	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

	def set_peak(self, gps):
		self.peak_time, self.peak_time_ns = gps.seconds, gps.nanoseconds

	def get_period(self):
		start = lal.LIGOTimeGPS(self.start_time, self.start_time_ns)
		return segments.segment(start, start + self.duration)

	def set_period(self, period):
		self.start_time, self.start_time_ns = period[0].seconds, period[0].nanoseconds
		self.duration = float(period.duration())

	def get_band(self):
		low = self.central_freq - self.bandwidth / 2
		return segments.segment(low, low + self.bandwidth)

	def set_band(self, band):
		self.central_freq = (band[0] + band[1])/2.0
		self.bandwidth = band.duration()


SnglBurstTable.RowType = SnglBurst


#
# =============================================================================
#
#                             sngl_inspiral:table
#
# =============================================================================
#

class SnglInspiralIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sngl_inspiral", "event_id", n)


class SnglInspiralTable(table.Table):	# FIXME: should be LSCTable
	tableName = "sngl_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"template_duration": "real_8",
		"event_duration": "real_8",
		"amplitude": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"mtotal": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"f_final": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "int_4s",
		"sigmasq": "real_8",
		"rsqveto_duration": "real_4",
		"Gamma0": "real_4",
		"Gamma1": "real_4",
		"Gamma2": "real_4",
		"Gamma3": "real_4",
		"Gamma4": "real_4",
		"Gamma5": "real_4",
		"Gamma6": "real_4",
		"Gamma7": "real_4",
		"Gamma8": "real_4",
		"Gamma9": "real_4",
		"event_id": "int_8s"	# FIXME: column should be ilwd
	}
	# FIXME:  event_id column needs to be changed to ilwd
	#ids = SnglInspiralIDs()

	def get_column(self,column):
		if column == 'effective_snr':
			return self.get_effective_snr()
		if column == 'snr_over_chi':
			return self.get_snr_over_chi()
		elif column == 'chirp_distance':
			return self.get_chirp_dist()
		else:
			return self.getColumnByName(column).asarray()

	def get_effective_snr(self):    
		snr = self.get_column('snr')
		chisq = self.get_column('chisq')
		chisq_dof = self.get_column('chisq_dof')
		return snr/ (1 + snr**2/250)**(0.25) / (chisq/(2*chisq_dof - 2) )**(0.25)
    
	def get_chirp_distance(self,ref_mass = 1.40):
		mchirp = self.get_column('mchirp')
		eff_dist = self.get_column('eff_distance')
		return eff_dist * (2.**(-1./5) * ref_mass / mchirp)**(5./6)

	def get_snr_over_chi(self):
		return self.get_column('snr')/self.get_column('chisq')**(1./2)
		
	def ifocut(self,ifo):
		ifoTrigs = table.new_from_template(self)
		for row in self:
			if row.ifo == ifo:
				ifoTrigs.append(row)
		return ifoTrigs

	def veto(self,seglist):
		vetoed = table.new_from_template(self)
		keep = table.new_from_template(self)
		for row in self:
			time = row.get_end()
			if time in seglist:
				vetoed.append(event)
			else:
				keep.append(event)

		return keep
	
	def getslide(self,slide_num):
		"""
		Return the triggers with a specific slide number.
		@param slide_num: the slide number to recover (contained in the event_id)
		"""
		slideTrigs = table.new_from_template(self)
		for row in self:
			if ( (row.event_id % 1000000000) / 100000 ) == slide_num:
				slideTrigs.append(row)
     
		return slideTrigs


class SnglInspiral(object):
	__slots__ = SnglInspiralTable.validcolumns.keys()

	def get_end(self):
		return lal.LIGOTimeGPS(self.end_time, self.end_time_ns)

	def set_end(self, gps):
		self.end_time, self.end_time_ns = gps.seconds, gps.nanoseconds

	def get_effective_snr(self):
		return self.snr/ (1 + self.snr**2/250)**(0.25)/(self.chisq/(2*self.chisq_dof - 2) )**(0.25) 


SnglInspiralTable.RowType = SnglInspiral


#
# =============================================================================
#
#                             sngl_ringdown:table
#
# =============================================================================
#

class SnglRingDownIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sngl_ringdown", "event_id", n)


class SnglRingDownTable(table.Table):	# FIXME: should be LSCTable
	tableName = "sngl_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"frequency": "real_4",
		"quality": "real_4",
		"mass": "real_4",
		"spin": "real_4",
		"snr": "real_4",
		"eff_distance": "real_4",
		"sigma_sq": "real_8",
		"event_id": "int_8s"	# FIXME: column should be ilwd
	}
	# FIXME:  event_id column needs to be changed to ilwd
	#ids = SnglRingDownIDs()


class SnglRingDown(object):
	__slots__ = SnglRingDownTable.validcolumns.keys()


SnglRingDownTable.RowType = SnglRingDown


#
# =============================================================================
#
#                             multi_inspiral:table
#
# =============================================================================
#

class MultiInspiralIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "multi_inspiral", "event_id", n)


class MultiInspiralTable(LSCTable):
	tableName = "multi_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifos": "lstring",
		"search": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"amplitude": "real_4",
		"ifo1_eff_distance": "real_4",
		"ifo2_eff_distance": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"ifo1_snr": "real_4",
		"ifo2_snr": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "real_4",
		"sigmasq": "real_4",
		"ligo_axis_ra": "real_4",
		"ligo_axis_dec": "real_4",
		"ligo_angle": "real_4",
		"ligo_angle_sig": "real_4",
		"inclination": "real_4",
		"polarization": "real_4",
		"event_id": "ilwd:char"
	}
	ids = MultiInspiralIDs()


class MultiInspiral(object):
	__slots__ = MultiInspiralTable.validcolumns.keys()


MultiInspiralTable.RowType = MultiInspiral


#
# =============================================================================
#
#                              sim_inspiral:table
#
# =============================================================================
#

class SimInspiralIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sim_inspiral", "simulation_id", n)


class SimInspiralTable(LSCTable):
	tableName = "sim_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_end_time": "int_4s",
		"geocent_end_time_ns": "int_4s",
		"h_end_time": "int_4s",
		"h_end_time_ns": "int_4s",
		"l_end_time": "int_4s",
		"l_end_time_ns": "int_4s",
		"g_end_time": "int_4s",
		"g_end_time_ns": "int_4s",
		"t_end_time": "int_4s",
		"t_end_time_ns": "int_4s",
		"v_end_time": "int_4s",
		"v_end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"source": "lstring",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"distance": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"inclination": "real_4",
		"coa_phase": "real_4",
		"polarization": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"spin1x": "real_4",
		"spin1y": "real_4",
		"spin1z": "real_4",
		"spin2x": "real_4",
		"spin2y": "real_4",
		"spin2z": "real_4",
		"theta0": "real_4",
		"phi0": "real_4",
		"f_lower": "real_4",
		"f_final": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"eff_dist_g": "real_4",
		"eff_dist_t": "real_4",
		"eff_dist_v": "real_4",
		"simulation_id": "ilwd:char"
	}
	ids = SimInspiralIDs()

	def get_column(self,column):
		if 'chirp_dist' in column:
			site = column[-1]
			return self.get_chirp_dist(site)
		elif column == 'spin1':
			return self.get_spin_mag(1)
		elif column == 'spin2':
			return self.get_spin_mag(2)
		else:
			return self.getColumnByName(column).asarray()

	def get_chirp_dist(self,site,ref_mass = 1.40):
		mchirp = self.get_column('mchirp')
		eff_dist = self.get_column('eff_dist_' + site)
		return eff_dist * (2.**(-1./5) * ref_mass / mchirp)**(5./6)

	def get_spin_mag(self,objectnumber):
		sx = self.get_column('spin' + str(objectnumber) + 'x')
		sy = self.get_column('spin' + str(objectnumber) + 'y')
		sz = self.get_column('spin' + str(objectnumber) + 'z')
		return (sx**2 + sy**2 + sz**2)**(0.5)

	def veto(self,seglist,site=None):
		keep = table.new_from_template(self)
		for row in self:
			time = row.get_end(site)
			if time not in seglist:
				keep.append(row)

		return keep


class SimInspiral(object):
	__slots__ = SimInspiralTable.validcolumns.keys()

	def get_end(self,site = None):
		if not site:
			return lal.LIGOTimeGPS(self.geocent_end_time, self.geocent_end_time_ns)
		else:
			return lal.LIGOTimeGPS(getattr(self,site + 'end_time'), getattr(self,site + 'end_time_ns'))


SimInspiralTable.RowType = SimInspiral


#
# =============================================================================
#
#                               sim_burst:table
#
# =============================================================================
#

class SimBurstIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sim_burst", "simulation_id", n)


class SimBurstTable(LSCTable):
	tableName = "sim_burst:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_peak_time": "int_4s",
		"geocent_peak_time_ns": "int_4s",
		"h_peak_time": "int_4s",
		"h_peak_time_ns": "int_4s",
		"l_peak_time": "int_4s",
		"l_peak_time_ns": "int_4s",
		"peak_time_gmst": "real_8",
		"dtminus": "real_4",
		"dtplus": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"coordinates": "lstring",
		"polarization": "real_4",
		"hrss": "real_4",
		"hpeak": "real_4",
		"distance": "real_4",
		"freq": "real_4",
		"tau": "real_4",
		"zm_number": "int_4s",
		"simulation_id": "ilwd:char"
	}
	ids = SimBurstIDs()


class SimBurst(object):
	__slots__ = SimBurstTable.validcolumns.keys()

	def cmp(self, other):
		"""
		Return 0 if self and other describe the same injection,
		non-0 otherwise.
		"""
		a = (
			self.waveform,
			self.geocent_peak_time,
			self.geocent_peak_time_ns,
			self.h_peak_time,
			self.h_peak_time_ns,
			self.l_peak_time,
			self.l_peak_time_ns,
			self.peak_time_gmst,
			self.dtminus,
			self.dtplus,
			self.longitude,
			self.latitude,
			self.coordinates,
			self.polarization,
			self.hrss,
			self.hpeak,
			self.distance,
			self.freq,
			self.tau,
			self.zm_number
		)
		b = (
			other.waveform,
			other.geocent_peak_time,
			other.geocent_peak_time_ns,
			other.h_peak_time,
			other.h_peak_time_ns,
			other.l_peak_time,
			other.l_peak_time_ns,
			other.peak_time_gmst,
			other.dtminus,
			other.dtplus,
			other.longitude,
			other.latitude,
			other.coordinates,
			other.polarization,
			other.hrss,
			other.hpeak,
			other.distance,
			other.freq,
			other.tau,
			other.zm_number
		)
		return cmp(a, b)

	def get_geocent_peak(self):
		return lal.LIGOTimeGPS(self.geocent_peak_time, self.geocent_peak_time_ns)

	def set_geocent_peak(self, gps):
		self.geocent_peak_time, self.geocent_peak_time_ns = gps.seconds, gps.nanoseconds

	def get_peak(self, instrument):
		observatory = instrument[0]
		if observatory == "H":
			return lal.LIGOTimeGPS(self.h_peak_time, self.h_peak_time_ns)
		if observatory == "L":
			return lal.LIGOTimeGPS(self.l_peak_time, self.l_peak_time_ns)
		raise ValueError, instrument

	def set_peak(self, instrument, gps):
		observatory = instrument[0]
		if observatory == "H":
			self.h_peak_time, self.h_peak_time_ns = gps.seconds, gps.nanoseconds
		if observatory == "L":
			self.l_peak_time, self.l_peak_time_ns = gps.seconds, gps.nanoseconds
		raise ValueError, instrument


SimBurstTable.RowType = SimBurst


#
# =============================================================================
#
#                              sim_ringdown:table
#
# =============================================================================
#

class SimRingDownIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sim_ringdown", "simulation_id", n)


class SimRingDownTable(LSCTable):
	tableName = "sim_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"coordinates": "lstring",
		"geocent_start_time": "int_4s",
		"geocent_start_time_ns": "int_4s",
		"h_start_time": "int_4s",
		"h_start_time_ns": "int_4s",
		"l_start_time": "int_4s",
		"l_start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"mass": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"distance": "real_4",
		"inclination": "real_4",
		"polarization": "real_4",
		"epsilon": "real_4",
		"spin": "real_4",
		"frequency": "real_4",
		"quality": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"h0": "real_4",
		"hrss": "real_4",
		"hrss_h": "real_4",
		"hrss_l": "real_4",
		"simulation_id": "ilwd:char"
	}
	ids = SimRingDownIDs()


class SimRingDown(object):
	__slots__ = SimRingDownTable.validcolumns.keys()


SimRingDownTable.RowType = SimRingDown


#
# =============================================================================
#
#                               summ_value:table
#
# =============================================================================
#

class SummValueTable(table.Table):
	tableName = "summ_value:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"ifo": "lstring",
		"name": "lstring",
		"value": "real_4",
		"comment": "lstring"
	}


class SummValue(object):
	__slots__ = SummValueTable.validcolumns.keys()


SummValueTable.RowType = SummValue


#
# =============================================================================
#
#                            sim_inst_params:table
#
# =============================================================================
#

class SimInstParamsIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "sim_inst_params", "simulation_id", n)


class SimInstParamsTable(LSCTable):
	tableName = "sim_inst_params:table"
	validcolumns = {
		"simulation_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"value": "real_8"
	}
	ids = SimInstParamsIDs()


class SimInstParams(object):
	__slots__ = SimInstParamsTable.validcolumns.keys()


SimInstParamsTable.RowType = SimInstParams


#
# =============================================================================
#
#                               stochastic:table
#
# =============================================================================
#

class StochasticTable(table.Table):
	tableName = "stochastic:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"duration": "int_4s",
		"duration_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"cc_stat": "real_8",
		"cc_sigma": "real_8"
	}


class Stochastic(object):
	__slots__ = StochasticTable.validcolumns.keys()


StochasticTable.RowType = Stochastic


#
# =============================================================================
#
#                               stochsumm:table
#
# =============================================================================
#

class StochSummTable(table.Table):
	tableName = "stochsumm:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"y_opt": "real_8",
		"error": "real_8"
	}


class StochSumm(object):
	__slots__ = StochSummTable.validcolumns.keys()


StochSummTable.RowType = StochSumm


#
# =============================================================================
#
#                            external_trigger:table
#
# =============================================================================
#

# FIXME: this table looks broken to me.  There is no unique ID column, thus
# it is not possible to refer to entries in this table from other tables.
# There is a "notice_id" column, but that is for recording the native
# identifier as used by the source of the trigger.  It cannot be relied
# upon to be unique within this table (two different sources might *happen*
# to use the same identifier format, like "event001").

class ExtTriggersTable(table.Table):
	tableName = "external_trigger:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"det_alts": "lstring",
		"det_band": "lstring",
		"det_fluence": "lstring",
		"det_fluence_int": "lstring",
		"det_name": "lstring",
		"det_peak": "lstring",
		"det_peak_int": "lstring",
		"det_snr": "lstring",
		"email_time": "int_4s",
		"event_dec": "real_4",
		"event_dec_err": "real_4",
		"event_epoch": "lstring",
		"event_err_type": "lstring",
		"event_ra": "real_4",
		"event_ra_err": "real_4",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"event_type": "lstring",
		"event_z": "real_4",
		"event_z_err": "real_4",
		"notice_comments": "lstring",
		"notice_id": "lstring",
		"notice_sequence": "lstring",
		"notice_time": "int_4s",
		"notice_type": "lstring",
		"notice_url": "lstring",
		"obs_fov_dec": "real_4",
		"obs_fov_dec_width": "real_4",
		"obs_fov_ra": "real_4",
		"obs_fov_ra_width": "real_4",
		"obs_loc_ele": "real_4",
		"obs_loc_lat": "real_4",
		"obs_loc_long": "real_4",
		"ligo_fave_lho": "real_4",
		"ligo_fave_llo": "real_4",
		"ligo_delay": "real_4",
		"event_number_gcn": "int_4s",
		"event_number_grb": "lstring",
		"event_status": "int_4s"
	}


class ExtTriggers(object):
	__slots__ = ExtTriggersTable.validcolumns.keys()


ExtTriggersTable.RowType = ExtTriggers


#
# =============================================================================
#
#                                 filter:table
#
# =============================================================================
#

class FilterTable(table.Table):
	tableName = "filter:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"program": "lstring",
		"start_time": "int_4s",
		"filter_name": "lstring",
		"comment": "lstring"
	}


class Filter(object):
	__slots__ = FilterTable.validcolumns.keys()


FilterTable.RowType = Filter


#
# =============================================================================
#
#                                segment:table
#
# =============================================================================
#

class SegmentIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "segment", "segment_id", n)


class SegmentTable(LSCTable):
	tableName = "segment:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"active": "int_4s",
		"segnum": "int_4s",
		"insertion_time": "int_4s"
	}
	ids = SegmentIDs()


class Segment(object):
	__slots__ = SegmentTable.validcolumns.keys()

	def get(self):
		"""
		Return the segment described by this row.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.start_time, self.start_time_ns), lal.LIGOTimeGPS(self.end_time, self.end_time_ns))

	def set(self, segment):
		"""
		Set the segment described by this row.
		"""
		self.start_time, self.start_time_ns = segment[0].seconds, segment[0].nanoseconds
		self.end_time, self.end_time_ns = segment[1].seconds, segment[1].nanoseconds

	def get_active(self):
		"""
		Return True if the segment is active, False if the segment
		is inactive and None if neither is the case.
		"""
		if self.active > 0:
			return True
		if self.active < 0:
			return False
		return None

	def set_active(self, active):
		"""
		Sets the segment to active if active is True, to inactive
		if active if False, and undefined if active is None.
		"""
		if active is None:
			self.active = 0
		elif active:
			self.active = 1
		else:
			self.active = -1
		return self


SegmentTable.RowType = Segment


#
# =============================================================================
#
#                            segment_def_map:table
#
# =============================================================================
#

class SegmentDefMapIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "segment_def_map", "seg_def_map_id", n)


class SegmentDefMapTable(LSCTable):
	tableName = "segment_def_map:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"seg_def_map_id": "ilwd:char",
		"segment_cdb": "int_4s",
		"segment_id": "ilwd:char",
		"segment_def_cdb": "int_4s",
		"segment_def_id": "ilwd:char",
		"state_vec_map": "int_4s",
		"insertion_time": "int_4s"
	}
	ids = SegmentDefMapIDs()


class SegmentDefMap(object):
	__slots__ = SegmentDefMapTable.validcolumns.keys()


SegmentDefMapTable.RowType = SegmentDefMap


#
# =============================================================================
#
#                            segment_definer:table
#
# =============================================================================
#

class SegmentDefIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "segment_definer", "segment_def_id", n)


class SegmentDefTable(LSCTable):
	tableName = "segment_definer:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_def_id": "ilwd:char",
		"run": "lstring",
		"ifos": "lstring",
		"name": "lstring",
		"version": "int_4s",
		"comment": "lstring",
		"state_vec_major": "int_4s",
		"state_vec_minor": "int_4s",
		"insertion_time": "int_4s"
	}
	ids = SegmentDefIDs()


class SegmentDef(object):
	__slots__ = SegmentDefTable.validcolumns.keys()


SegmentDefTable.RowType = SegmentDef


#
# =============================================================================
#
#                               time_slide:table
#
# =============================================================================
#

class TimeSlideIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "time_slide", "time_slide_id", n)


class TimeSlideTable(LSCTable):
	tableName = "time_slide:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"instrument": "lstring",
		"offset": "real_8"
	}
	ids = TimeSlideIDs()

	def get_offset_dict(self, id):
		"""
		Return a dictionary of instrument/offset pairs as described
		by the rows having the given ID.
		"""
		d = {}
		for row in self.dict[id]:
			if row.instrument in d:
				raise KeyError, "%s: duplicate instrument %s" % (id, row.instrument)
			d[row.instrument] = row.offset
		return d

	def is_null(self, id):
		"""
		Test that a time slide ID identifies an all-zero time
		slide.
		"""
		for offset in self.get_offset_dict(id).itervalues():
			if offset:
				return False
		return True


class TimeSlide(object):
	__slots__ = TimeSlideTable.validcolumns.keys()


TimeSlideTable.RowType = TimeSlide


#
# =============================================================================
#
#                             coinc_definer:table
#
# =============================================================================
#

class CoincDefIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "coinc_definer", "coinc_def_id", n)


class CoincDefTable(LSCTable):
	tableName = "coinc_definer:table"
	validcolumns = {
		"coinc_def_id": "ilwd:char",
		"table_name": "lstring"
	}
	ids = CoincDefIDs()

	def get_contributors(self, id):
		"""
		Return a list of contributing table names for the given ID.
		"""
		l = [row.table_name for row in self.dict[id]]
		l.sort()
		return l


class CoincDef(object):
	__slots__ = CoincDefTable.validcolumns.keys()


CoincDefTable.RowType = CoincDef


#
# =============================================================================
#
#                              coinc_event:table
#
# =============================================================================
#

class CoincIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "coinc_event", "coinc_event_id", n)


class CoincTable(LSCTable):
	tableName = "coinc_event:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"coinc_def_id": "ilwd:char",
		"coinc_event_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"nevents": "int_4u"
	}
	ids = CoincIDs()


class Coinc(object):
	__slots__ = CoincTable.validcolumns.keys()


CoincTable.RowType = Coinc


#
# =============================================================================
#
#                            coinc_event_map:table
#
# =============================================================================
#

class CoincMapTable(table.Table):
	tableName = "coinc_event_map:table"
	validcolumns = {
		"coinc_event_id": "ilwd:char",
		"event_id": "ilwd:char"
	}


class CoincMap(object):
	__slots__ = CoincMapTable.validcolumns.keys()


CoincMapTable.RowType = CoincMap


#
# =============================================================================
#
#                               ligolw_mon:table
#
# =============================================================================
#

class LIGOLWMonIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "ligolw_mon", "event_id", n)


class LIGOLWMonTable(LSCTable):
	tableName = "ligolw_mon:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"time": "int_4s",
		"time_ns": "int_4s",
		"amplitude": "real_8",
		"confidence": "real_8",
		"frequency": "real_8",
		"event_id": "ilwd:char",
		"insertion_time": "int_4s"
	}
	ids = LIGOLWMonIDs()


class LIGOLWMon(object):
	__slots__ = LIGOLWMonTable.validcolumns.keys()

	def get_time(self):
		return lal.LIGOTimeGPS(self.time, self.time_ns)

	def set_time(self, gps):
		self.time, self.time_ns = gps.seconds, gps.nanoseconds


LIGOLWMonTable.RowType = LIGOLWMon


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#

#
# Table name ---> table type mapping.
#

TableByName = {
	table.StripTableName(ProcessTable.tableName): ProcessTable,
	table.StripTableName(LfnTable.tableName): LfnTable,
	table.StripTableName(ProcessParamsTable.tableName): ProcessParamsTable,
	table.StripTableName(SearchSummaryTable.tableName): SearchSummaryTable,
	table.StripTableName(SearchSummVarsTable.tableName): SearchSummVarsTable,
	table.StripTableName(SnglBurstTable.tableName): SnglBurstTable,
	table.StripTableName(SnglInspiralTable.tableName): SnglInspiralTable,
	table.StripTableName(SnglRingDownTable.tableName): SnglRingDownTable,
	table.StripTableName(MultiInspiralTable.tableName): MultiInspiralTable,
	table.StripTableName(SimInspiralTable.tableName): SimInspiralTable,
	table.StripTableName(SimBurstTable.tableName): SimBurstTable,
	table.StripTableName(SimRingDownTable.tableName): SimRingDownTable,
	table.StripTableName(SummValueTable.tableName): SummValueTable,
	table.StripTableName(SimInstParamsTable.tableName): SimInstParamsTable,
	table.StripTableName(StochasticTable.tableName): StochasticTable,
	table.StripTableName(StochSummTable.tableName): StochSummTable,
	table.StripTableName(ExtTriggersTable.tableName): ExtTriggersTable,
	table.StripTableName(FilterTable.tableName): FilterTable,
	table.StripTableName(SegmentTable.tableName): SegmentTable,
	table.StripTableName(SegmentDefMapTable.tableName): SegmentDefMapTable,
	table.StripTableName(SegmentDefTable.tableName): SegmentDefTable,
	table.StripTableName(TimeSlideTable.tableName): TimeSlideTable,
	table.StripTableName(CoincDefTable.tableName): CoincDefTable,
	table.StripTableName(CoincTable.tableName): CoincTable,
	table.StripTableName(CoincMapTable.tableName): CoincMapTable,
	table.StripTableName(LIGOLWMonTable.tableName): LIGOLWMonTable
}


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#


#
# Override portions of the ligolw.LIGOLWContentHandler class
#

__parent_startTable = ligolw.LIGOLWContentHandler.startTable

def startTable(self, attrs):
	name = table.StripTableName(attrs["Name"])
	if name in TableByName:
		return TableByName[name](attrs)
	return __parent_startTable(self, attrs)

ligolw.LIGOLWContentHandler.startTable = startTable
