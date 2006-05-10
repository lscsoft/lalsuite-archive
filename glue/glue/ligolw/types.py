#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
Definitions of type strings found in LIGO Light Weight XML files.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                               Type Information
#
# =============================================================================
#

StringTypes = ["char_s", "ilwd:char", "ilwd:char_u", "lstring", "string"]
IntTypes = ["int_2s", "int_2u", "int_4s", "int_4u", "int_8s", "int_8u", "int"]
FloatTypes = ["real_4", "real_8", "float", "double"]

Types = StringTypes + IntTypes + FloatTypes

ToNumArrayType = {
	"int_2s": "Int16",
	"int_2u": "UInt16",
	"int_4s": "Int32",
	"int_4u": "UInt32",
	"int_8s": "Int64",
	"int_8u": "UInt64",
	"int": "Int32",
	"real_4": "Float32",
	"real_8": "Float64",
	"float": "Float64",
	"double": "Float64"
}
