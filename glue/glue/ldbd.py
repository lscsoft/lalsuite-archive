"""
lightweight database dumper

Copyright (C) 2003 Duncan Brown
 
This file is part of the lightweight datapase dumper (ldbd)

The ldbd module provides classes for manipulating LIGO metadata database
tables.

References:
http://www.ligo.caltech.edu/docs/T/T990101-02.pdf
http://www.ligo.caltech.edu/docs/T/T990023-01.pdf
http://ldas-sw.ligo.caltech.edu/doc/db2/doc/html/ilwdformat.html
"""

__author__ = 'Duncan Brown <dbrown@ligo.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]
# $Source$

import os
import sys
import string
import re
import csv
try:
  import mx.ODBC.DB2 as mxdb
  from mx.ODBC.DB2 import SQL
except:
  pass

try:
  import thread
except:
  """
  We're running on a single-threaded OS (or the Python interpreter has
  not been compiled to support threads) so return a standard dictionary.
  """
  _tss = {}
  def get_thread_storage():
    return _tss
  def destroy_thread_storage():
    del _tss
else:
  _tss = {}
  _tss_lock = thread.allocate_lock()
  def get_thread_storage():
    """Return a thread-specific storage dictionary."""
    thread_id = thread.get_ident()
    tss = _tss.get(thread_id)
    if tss is None:
      try:
        _tss_lock.acquire()
        _tss[thread_id] = tss = {}
      finally:
        _tss_lock.release()
    return tss
  def destroy_thread_storage():
    """Destroy a thread-specific storage dictionary."""
    thread_id = thread.get_ident()
    tss = _tss.get(thread_id)
    try:
      _tss_lock.acquire()
      del tss
      del _tss[thread_id]
    finally:
      _tss_lock.release()
    return

class LIGOLwParseError(Exception):
  """Error parsing LIGO lightweight XML file"""
  def __init__(self,args=None):
    self.args = args

class LIGOLwDBError(Exception):
  """Error interacting with database"""
  def __init__(self,args=None):
    self.args = args

class LIGOMetadataDatabase:
  """
  Contains a tuple of tables in the order that insertions should
  ocour and a dictionary of mandatory uniquw id fields for each
  table that must be populated.
  """
  def __init__(self,database):
    """
    database = the name of the LIGO database to initalize
    """
    self.database = database
    self.uniqueids = {}
    conn = mxdb.Connect(database)
    curs = conn.cursor()
    curs.execute("SELECT tabname FROM syscat.tables WHERE definer<>'SYSIBM'")
    self.tables = curs.fetchall()
    curs.execute("SELECT tabname, colname FROM syscat.columns " +
      "WHERE typename = 'CHARACTER' AND length = 13 AND nulls = 'N'")
    for tab, col in curs.fetchall():
      tab = tab.lower()
      col = col.lower()
      try:
        self.uniqueids[tab][col] = 'ilwd:char'
      except KeyError:
        self.uniqueids[tab] = {}
        self.uniqueids[tab][col] = 'ilwd:char'
    curs.close()
    conn.close()

class UniqueIds:
  """
  Contains a dictionary of unique ids which can be queried based
  on name. If a unique id does not exist in the dictionaty, one
  is fetched from the database. The unique id dictionary is specifc
  to the thread of execution so that one ligolw parser can be used
  in muitiple threads.
  """
  def __init__(self,curs):
    """
    curs = database cursor to the currently open database
    """
    get_thread_storage()['uqids'] = {}
    get_thread_storage()['curs'] = curs

  def __del__(self):
    destroy_thread_storage()

  def lookup(self,istring):
    """
    istring = the ilwd:char string corresponding to a unique id
    """
    try:
      return get_thread_storage()['uqids'][istring]
    except KeyError:
      curs = get_thread_storage()['curs']
      curs.execute('VALUES GENERATE_UNIQUE()')
      get_thread_storage()['uqids'][istring] = curs.fetchone()[0]
      return get_thread_storage()['uqids'][istring]

class LIGOLwParser:
  """
  Provides methods for parsing the data from a LIGO lightweight XML
  file parsed with pyRXP into a dictionary
  """

  def __init__(self):
    """
    Initializes a LIGO lightweight XML parser with the necessary 
    regular expressions and function for tuple translation

    Before calling parsteuple() the user must tell the class
    which unique id dictionary to use for the file by:
    
      p = LIGOLwParser()
      p.unique = UniqueIds(mycursor)
    """
    self.tabrx = re.compile(r'(\A[a-z0-9_]+:|\A)([a-z0-9_]+):table\Z')
    self.colrx = re.compile(r'(\A[a-z0-9_]+:|\A)([a-z0-9_]+:|\A)([a-z0-9_]+)\Z')
    self.llsrx = re.compile(r'\A\s*"')
    self.rlsrx = re.compile(r'"\s*\Z')
    self.licrx = re.compile(r'\A\s+"')
    self.ricrx = re.compile(r'"*\s*\Z')
    self.octrx = re.compile(r'\A\\[0-9][0-9][0-9]')
    self.dlmrx = re.compile(r'\\,')
    self.cp = csv.parser(0,',',1)
    self.unique = None
    self.types = {
      'int_2s' : int,
      'int_4s' : int,
      'real_4' : float,
      'real_8' : float,
      'lstring' : self.__lstring,
      'ilwd:char' : self.__ilwdchar,
      'ilwd:char_u' : self.__ilwdchar
    }

  def __lstring(self,lstr):
    """
    Returns a parsed lstring by stripping out and instances of
    the escaped delimiter. Sometimes the raw lstring has whitespace
    and a double quote at the beginning or end. If present, these
    are removed.
    """
    lstr = self.llsrx.sub('',lstr.encode('ascii'))
    lstr = self.rlsrx.sub('',lstr)
    lstr = self.dlmrx.sub(',',lstr)
    return lstr
    
  def __ilwdchar(self,istr):
    """
    If the ilwd:char field contains octal data, it is translated
    to a binary string and returned. Otherwise a lookup is done
    in the unique id dictionary and a binary string containing the
    correct unique id is returned.
    """
    istr = self.licrx.sub('',istr.encode('ascii'))
    istr = self.ricrx.sub('',istr)
    if self.octrx.match(istr):
      exec "istr = '"+istr+"'"
    else:
      try:
        istr = self.unique.lookup(istr)
      except AttributeError:
        raise LIGOLwParseError, 'unique id table has not been initialized'
    return istr

  def parsetuple(self,xmltuple):
    """
    Parse an XML tuple returned by pyRXP into a dictionary
    of LIGO metadata elements. The dictionary contains one
    entry for each table found in the XML tuple.
    """
    # first extract all the table and columns from the tuple from the
    # children of the ligo lightweight parent tuple
    table = {}
    tupleidx = 0
    for tag in xmltuple[2]:
      if tag[0] == 'Table' or tag[0] == 'table':
        tab = tag[1]['Name'].encode('ascii').lower()
        try:
          tab = self.tabrx.match(tab).group(2)
        except AttributeError:
          raise LIGOLwParseError, 'unable to parse a valid table name '+tab
        # initalize the table dictionary for this table
        table[tab] = { 
          'pos' : tupleidx,
          'column' : {},
          'stream' : (), 
          'query' : ''
          }
        # parse for columns in the tables children
        # look for the column name and type in the attributes
        # store the index in which the columns were found as
        # we need this to decode the stream later
        for subtag in tag[2]:
          if subtag[0] == 'Column' or subtag[0] == 'column':
            col = subtag[1]['Name'].encode('ascii').lower()
            try:
              col = self.colrx.match(col).group(3)
            except AttributeError:
              raise LIGOLwParseError, 'unable to parse a valid column name '+col
            try:
              typ = subtag[1]['Type'].encode('ascii').lower()
            except KeyError:
              raise LIGOLwParseError, 'type is missing for column '+col
            table[tab]['column'][col] = typ
            table[tab].setdefault('orderedcol',[]).append(col)
      tupleidx += 1

    # now iterate the dictionary of tables we have created looking for streams
    for tab in table.keys():
      for tag in xmltuple[2][table[tab]['pos']][2]:
        if tag[0] == 'Stream' or tag[0] == 'stream':
          # store the stream delimiter and create the esacpe regex
          try:
            delim = tag[1]['Delimiter'].encode('ascii')
          except KeyError:
            raise LIGOLwParseError, 'stream is missing delimiter'
          if delim != ',':
            raise LIGOLwParseError, 'unable to handle stream delimiter: '+delim

          # strip newlines from the stream and parse it
          stream = self.cp.parse(re.sub(r'\n','',tag[2][0]))

          # turn the csv stream into a list of lists
          slen = len(stream)
          ntyp = len(table[tab]['column'])
          mlen, lft = divmod(slen,ntyp)
          if lft != 0:
            raise LIGOLwParseError, 'invalid stream length for given columns'
          lst = [[None] * ntyp for i in range(mlen)]

          # translate the stream data to the correct data types
          for i in range(slen):
            j, k = divmod(i,ntyp)
            try:
              thiscol = table[tab]['orderedcol'][k]
              lst[j][k] = self.types[table[tab]['column'][thiscol]](stream[i])
            except KeyError, ValueError:
              raise LIGOLwParseError, 'error translating data in stream'
          table[tab]['stream'] = map(tuple,lst)

    # return the created table to the caller
    return table
          

class LIGOMetadata:
  """
  LIGO Metadata object class. Contains methods for parsing a LIGO
  lightweight XML file and inserting it into a database, executing
  and SQL query to retrive data from the database and writing it
  to a LIGO lightweight XML file
  """
  def __init__(self,ldb,xmlparser,lwtparser):
    """
    Connects to the database and creates a cursor. Initializes the unique
    id table for this LIGO lw document.

    ldb = LIGOMetadataDatabase object
    xmlparser = pyRXP XML to tuple parser object
    lwtparser = LIGOLwParser object (tuple parser)
    """
    self.ldb = ldb
    self.dbcon = mxdb.Connect(self.ldb.database)
    self.dbcon.setconnectoption(SQL.AUTOCOMMIT, SQL.AUTOCOMMIT_OFF)
    self.curs = self.dbcon.cursor()
    self.xmlparser = xmlparser
    self.lwtparser = lwtparser
    self.lwtparser.unique = None
    self.table = {}

  def __del__(self):
    if self.lwtparser.unique:
      del self.lwtparser.unique
    self.curs.close()
    self.dbcon.close()

  def reset(self):
    """Clear any existing table"""
    if self.table:
      del self.table
    self.table = {}

  def parse(self,xml):
    """
    Parses an XML document into a form read for insertion into the database

    xml = the xml document to be parsed
    """
    ligolwtup = self.xmlparser(xml)
    self.lwtparser.unique = UniqueIds(self.curs)
    self.table = self.lwtparser.parsetuple(ligolwtup)

  def add_lfn(self,lfn):
    """
    Add an LFN table to a parsed LIGO_LW XML document.

    lfn = lfn to be added
    """
    # get the process_id from the process table
    pid_col = self.table['process']['orderedcol'].index('process_id')
    pid = self.table['process']['stream'][0][pid_col]
    try:
      self.table['lfn']['stream'].append((pid,lfn))
    except KeyError:
      self.table['lfn'] = {
        'pos' : 0,
        'column' : {'process_id' : 'ilwd:char', 'lfn' : 'lstring'},
        'stream' : [(pid, lfn)],
        'query' : '',
        'orderedcol' : ['process_id', 'lfn' ]
        }

  def set_dn(self,dn):
    """
    Add an gridcert table to a parsed LIGO_LW XML document.

    dn = dn to be added
    """
    # get the process_id from the process table
    pid_col = self.table['process']['orderedcol'].index('process_id')
    pid = self.table['process']['stream'][0][pid_col]
    self.table['gridcert'] = {
      'pos' : 0,
      'column' : {'process_id' : 'ilwd:char', 'dn' : 'lstring'},
      'stream' : [(pid, dn)],
      'query' : '',
      'orderedcol' : ['process_id', 'dn' ]
      }
    
  def insert(self):
    """Insert the object into the database"""
    if len(self.table) == 0:
      raise LIGOLwDBError, 'attempt to insert empty table'
    for tab in self.table.keys():
      # find and add any missing unique ids
      generate = []
      missingcols = [k for k in self.ldb.uniqueids[tab] 
        if k not in self.table[tab]['column']]
      for m in missingcols:
        generate.append(',GENERATE_UNIQUE()')
        self.table[tab]['orderedcol'].append(m)
      # and construct the sql query
      self.table[tab]['query'] = ' '.join( 
        ['INSERT INTO', tab, '(', ','.join(self.table[tab]['orderedcol']), 
        ') VALUES (', ','.join(['?' for x in self.table[tab]['column']]) , 
        ''.join(generate), ')'])
    for tabtup in self.ldb.tables:
      tab = tabtup[0].lower()
      try:
        try: 
          self.curs.execute(self.table[tab]['query'],
            self.table[tab]['stream'])
          rowcount = self.curs.rowcount
        except mxdb.Error, e:
          self.dbcon.rollback()
          raise LIGOLwDBError, e[2]
        except mxdb.Warning, e:
          self.dbcon.rollback()
          raise LIGOLwDBError, e[2]
      except KeyError:
        pass
    self.dbcon.commit()
    return rowcount


  def select(self,sql):
    """
    Execute an SQL select statement and stuff the results into a
    dictionary.

    sql = the (case sensitve) SQL statment to execute
    """
    if len(self.table) != 0:
      raise LIGOLwDBError, 'attempt to fill non-empty table from database'
    ligolw = ''
    self.table = {}
    sqltypes = {
      -2 : 'ilwd:char_u',
      1 : 'lstring',
      4  : 'int_4s',
      7 : 'real_4',
      8 : 'real_8',
      12 : 'lstring',
      93 : 'lstring', 
      }
    try:
      tab = re.compile(r'[Ff][Rr][Oo][Mm]\s+([A-Za-z0-0_]+)([,\s]+|$)').search(sql).group(1)
    except AttributeError:
      raise LIGOLwDBError, 'could not find table name in query ' + str(sql)
    self.table[tab] = {
      'pos' : 0,
      'column' : {},
      'stream' : (), 
      'query' : sql
      }
    try:
      self.curs.execute(sql)
    except mxdb.Error, e:
      raise LIGOLwDBError, e[2]
    desc = self.curs.description
    for col,typ,disp,intsz,prec,sca,nul in desc:
      try:
        self.table[tab]['column'][col] = sqltypes[typ]
      except KeyError:
        raise LIGOLwDBError, 'unknown type returned by database ' + str(typ)
      self.table[tab].setdefault('orderedcol',[]).append(col)

    try:
      self.table[tab]['stream'] = self.curs.fetchall()
    except mxdb.Error, e:
      raise LIGOLwDBError, e[2]

    return len(self.table[tab]['stream'])

  def xml(self):
    """Convert a table dictionary to LIGO lightweight XML"""
    if len(self.table) == 0:
      raise LIGOLwDBError, 'attempt to convert empty table to xml'
    ligolw = """\
<?xml version='1.0' encoding='utf-8' ?>
<!DOCTYPE LIGO_LW [
<!ELEMENT LIGO_LW ((LIGO_LW|Comment|Param|Table|Array|Stream|IGWDFrame|AdcData|AdcInterval|Time|Detector)*)>
<!ATTLIST LIGO_LW
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED>

<!ELEMENT Comment (#PCDATA)>

<!ELEMENT Param (#PCDATA|Comment)*>
<!ATTLIST Param 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Start CDATA #IMPLIED
          Scale CDATA #IMPLIED
          Unit CDATA #IMPLIED
          DataUnit CDATA #IMPLIED>

<!ELEMENT Table (Comment?,Column*,Stream?)>
<!ATTLIST Table 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED>

<!ELEMENT Column EMPTY>
<!ATTLIST Column
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Unit CDATA #IMPLIED>

<!ELEMENT Array (Dim*,Stream?)>
<!ATTLIST Array 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Unit CDATA #IMPLIED>

<!ELEMENT Dim (#PCDATA)>
<!ATTLIST Dim 
          Name  CDATA #IMPLIED
          Unit CDATA #IMPLIED
          Start CDATA #IMPLIED
          Scale CDATA #IMPLIED>

<!ELEMENT Stream (#PCDATA)>
<!ATTLIST Stream 
          Name      CDATA #IMPLIED
          Type      (Remote|Local) "Local"
          Delimiter CDATA ","
          Encoding  CDATA #IMPLIED
          Content   CDATA #IMPLIED>

<!ELEMENT IGWDFrame ((Comment|Param|Time|Detector|AdcData|LIGO_LW|Stream?|Array|IGWDFrame)*)>
<!ATTLIST IGWDFrame 
          Name CDATA #IMPLIED>

<!ELEMENT Detector ((Comment|Param|LIGO_LW)*)>
<!ATTLIST Detector 
          Name CDATA #IMPLIED>

<!ELEMENT AdcData ((AdcData|Comment|Param|Time|LIGO_LW|Array)*)>
<!ATTLIST AdcData 
          Name CDATA #IMPLIED>

<!ELEMENT AdcInterval ((AdcData|Comment|Time)*)>
<!ATTLIST AdcInterval 
          Name CDATA #IMPLIED
          StartTime CDATA #IMPLIED
          DeltaT CDATA #IMPLIED>

<!ELEMENT Time (#PCDATA)>
<!ATTLIST Time 
          Name CDATA #IMPLIED
          Type (GPS|Unix|ISO-8601) "ISO-8601">
]>

<LIGO_LW>
"""

    for tab in self.table.keys():
      ligolw += '   <Comment>'+self.table[tab]['query']+'</Comment>\n'
      ligolw += '   <Table Name="'+tab+':table">\n'
      for col in self.table[tab]['orderedcol']:
        ligolw +='      <Column Name="'+tab.lower()+':'+col.lower()+'" Type="'+self.table[tab]['column'][col].lower()+'"/>\n'
      ligolw += '      <Stream Name="'+tab.lower()+':table" Type="Local" Delimiter=",">\n'
      stridx = 0
      ligolw += '      '
      for tup in self.table[tab]['stream']:
        if stridx != 0:
          ligolw += ',\n      '
        colidx = 0
        for tupi in tup:
          if tupi is not None:
            coltype = self.table[tab]['column'][self.table[tab]['orderedcol'][colidx]]
            if re.match(r'\Ailwd:char_u\Z',coltype):
              ligolw += '"'
              for ch in str(tupi):
                ligolw += '\\%.3o' % (ord(ch))
              ligolw += '"'
            elif re.match(r'\Alstring\Z',coltype):
              ligolw += '"'+str(tupi)+'"'
            elif re.match(r'\Areal_4\Z',coltype):
              ligolw += '%13.7e' % tupi
            elif re.match(r'\Areal_8\Z',coltype):
              ligolw += '%22.16e' % tupi
            else:
              ligolw += str(tupi)
          else:
            ligolw += ''
          if colidx < (len(self.table[tab]['column']) - 1):
            ligolw += ','
          colidx += 1
        stridx += 1
      ligolw += '\n      </Stream>\n'
      ligolw += '   </Table>\n'
    ligolw += '</LIGO_LW>'

    return ligolw
