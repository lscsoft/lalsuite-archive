"""
High-level Table element manipulation utilities.
"""


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
