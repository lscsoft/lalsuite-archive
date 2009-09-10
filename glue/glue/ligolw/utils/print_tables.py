# xml_convert.py

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
A collection of utilities to convert xml-tables to other formats, such as
wiki or html.
"""
import sys

from glue.ligolw import ligolw
from glue.ligolw import table

__author__ = "Collin Capano <cdcapano@ligo.caltech.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Utilities
#
# =============================================================================
#

def set_output_format( output_format ):
    """
    Sets output format; returns standard bits of table. These are:
        ttx: how to start a title for a set of tables
        xtt: how to end a title for a set of tables
        tx: how to start a table
        xt: how to close a table
        capx: how to start a caption for the table
        xcap: how to close a caption for the table
        rx: how to start a row and the first cell in the row
        xr: how to close a row and the last cell in the row
        xccx: how to close a cell and open a new one
    """
    if output_format == 'wiki':
        ttx = '== '
        xtt = ' =='
        tx = ''
        xt = ''
        capx = "'''"
        xcap = "'''"
        rx = '||'
        xr = '||'
        xccx = '||'

    elif output_format == "html":
        ttx = '<b>'
        xtt = '</b><hr>'
        tx = '<table border = "1">'
        xt = '</table><br><br>'
        capx = '<caption>'
        xcap = '</caption>'
        rx = '<tr><td>'
        xr = '</td></tr>'
        xccx = '</td><td>'

    else:
        raise ValueError, "unrecognized output_format %s" % output_format

    return ttx, xtt, tx, xt, capx, xcap, rx, xr, xccx


def smart_round( val, decimal_places = 2):
    """
    For floats >= 10.**-(decimal_places - 1), rounds off to the valber of decimal places specified.
    For floats < 10.**-(decimal_places - 1), puts in exponential form then rounds off to the decimal
    places specified.

    @val: value to round. If val is not a float, just returns val
    @decimal_places: number of decimal places to round to
    """
    if isinstance(val, float) and val != 0.0:
        if val >= 10.**-(decimal_places - 1):
            conv_str = ''.join([ '%.', str(decimal_places), 'f' ])
        else:
            conv_str = ''.join([ '%.', str(decimal_places), 'e' ])
        val = float( conv_str % val )

    return val


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def print_tables(xmldoc, output, output_format, tableList = [], round_floats = True, decimal_places = 2, title = None):
    """
    Method to print tables in an xml file in other formats.
    Input is an xmldoc, output is a file object containing the
    tables.

    @xmldoc: document to convert
    @output: file object to write output to; if None, will write to stdout
    @output_format: format to convert to
    @tableList: only convert the listed tables. Default is
     to convert all the tables found in the xmldoc. Tables
     not converted will not be included in the returned file
     object.
    @round_floats: If turned on, will smart_round floats to specifed
     number of places.
    @decimal_places: If round_floats turned on, will smart_round to this
     number of decimal places.
    """
    # get the tables to convert
    if tableList == []:
        tableList = [tb.getAttribute("Name") for tb in xmldoc.childNodes[0].getElementsByTagName(u'Table')]

    # set the output
    if output is None:
        output = sys.stdout

    # get table bits
    ttx, xtt, tx, xt, capx, xcap, rx, xr, xccx = set_output_format( output_format )

    # set the title if desired
    if title is not None:
        print >> output, "%s%s%s" %(ttx,str(title),xtt)
    # cycle over the tables in the xmldoc
    for table_name in tableList:
        this_table = table.get_table(xmldoc, table_name)
        col_names = [col.getAttribute("Name").split(":")[-1] for col in this_table.getElementsByTagName(u'Column')]
        # start the table and print table name
        print >> output, tx
        print >> output, "%s%s%s" %(capx, table_name, xcap)
        print >> output, "%s%s%s" %(rx, xccx.join(col_names), xr)

        # print the data in the table
        for row in this_table:
            if round_floats:
                out_row = [ str(smart_round( getattr(row, col_name), decimal_places = decimal_places )) for col_name in col_names ]
            else:
                out_row = [ str(getattr(row, col_name)) for col_name in col_names ]
            print >> output, "%s%s%s" %(rx, xccx.join(out_row), xr)

        # close the table and go on to the next
        print >> output, xt
