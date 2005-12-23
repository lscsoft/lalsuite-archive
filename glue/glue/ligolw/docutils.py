"""
High-level Table element manipulation utilities.
"""

from ligolw import ElementError
import lsctables


def MergeElements(elem1, elem2):
	"""
	Move the children of elem2 to elem1, and unlink elem2 from its
	parent.  The return value is elem1.  If the two elements are
	tables, then the contents of the second table into the first table,
	and unlink the second table from the document tree.  This does
	*not* check that this operation is OK, i.e. if the two tables are
	compatible.
	"""
	# FIXME: what about attributes?
	if elem1.tagName != elem2.tagName:
		raise ElementError, "MergeElements(elem1, elem2): elements must have same names"
	if elem1.tagName == "LIGO_LW":
		map(elem1.appendChild, elem2.childNodes)
	elif elem1.tagName == "Table":
		elem1.rows.extend(elem2.rows)
	if elem2.parentNode:
		elem2.parentNode.removeChild(elem2)
	return elem1


def TablesAreCompatible(table1, table2):
	"""
	Return True if the two tables have compatible columns.  This
	implies that a row object from one can be placed in the other.
	"""
	names1 = [column.getAttribute("Name") for column in table1.getElementsByTagName("Column")]
	names1.sort()
	names2 = [column.getAttribute("Name") for column in table2.getElementsByTagName("Column")]
	names2.sort()
	return names1 == names2


def MergeCompatibleTables(elem):
	"""
	Below the given element, find all Tables whose structure is
	described in lsctables, and merge compatible ones of like type.
	That is, merge all SnglBurstTables that have the same columns into
	a single table, etc..
	"""
	for tname in lsctables.TableByName.keys():
		tables = [t for t in elem.getElementsByTagName("Table") if t.getAttribute("Name") == tname]
		for i in range(1, len(tables)):
			if TablesAreCompatible(tables[0], tables[i]):
				MergeElements(tables[0], tables[i])


def RemapProcessIDs(elem, mapping):
	"""
	Recurse over all Table elements whose structure is described in the
	module lsctables, and remap the process IDs of all rows according
	to mapping.  Rows with process IDs not named in the mapping are
	ignored.
	"""
	for table in [t for t in elem.getElementsByTagName("Table") if t.getAttribute("Name") in lsctables.TableByName.keys()]:
		for row in table.rows:
			try:
				row.process_id = mapping[row.process_id]
			except KeyError:
				pass


def GetProcessParams(table, id):
	"""
	Turn the rows matching the process ID process_id in the process
	params table "table" into a dictionary of parameter/value pairs.
	"""
	if table.getAttribute("Name") != lsctables.ProcessParamsTable.tableName:
		raise ValueError, "GetProcessParams(table, id): table is not a process params table"
	params = {}
	for row in table.rows:
		if row.process_id != id:
			pass
		elif params.has_key(row.param):
			raise ElementError, "duplicate process parameter %s for process ID %s" % (row.param, id)
		elif row.type in metaio.StringTypes:
			params[row.param] = str(row.value)
		elif row.type in metaio.IntTypes:
			params[row.param] = int(row.value)
		elif row.type in metaio.FloatTypes:
			params[row.param] = float(row.value)
		else:
			raise ElementError, "unrecognized Type attribute %s" % row.type
	return params
