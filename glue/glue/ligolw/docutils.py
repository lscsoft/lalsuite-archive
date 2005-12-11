"""
High-level Table element manipulation utilities.
"""

import lsctables


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


def MergeTables(table1, table2):
	"""
	Transfer the contents of the second table into the first table, and
	unlink the second table from the document tree.  This does *not*
	check that this operation is OK, i.e. if the two tables are
	compatible.
	"""
	table1.rows.extend(table2.rows)
	if table2.parentNode:
		table2.parentNode.removeChild(table2)


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
				MergeTables(tables[0], tables[i])
