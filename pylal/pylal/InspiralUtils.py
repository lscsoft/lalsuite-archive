from glue import lal
from glue import segments
import socket, os

def WriteProcessParams(page, name, version, command_line): 
  """
  function to write out the process params that the code was called with
  """
  try:
    import markup
    from markup import oneliner as e
  except: 
    raise ImportError("Require markup.py to generate the html page")
  text = "Figures produced with "+name+", " \
      + version[1:len(version)-1] + ", invoked with arguments:\n\n" \
      + name
  for arg in command_line:
    text += " " +  arg
  page.pre( text )
  return page

def AddFileToCache(fname, cache):
  """
  Add the given file to the lal.Cache
  """
  file_name = fname.split('.')[0].split('-')
  cache.append(lal.CacheEntry( file_name[0], file_name[1],
    segments.segment(int(file_name[2]), 
      int(file_name[2]) + int(file_name[3])),
    'file://' + socket.gethostbyaddr(socket.gethostname())[0] + \
     os.getcwd() + '/' + fname))

def GenerateCache(fileList):
  """
  Generate a lal.Cache for the list of files
  """
  cache = lal.Cache()
  for file in fileList:
    AddFileToCache(file, cache)
  return(cache)
