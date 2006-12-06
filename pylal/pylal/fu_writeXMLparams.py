from pylal.fu_utils import *

##############################################################################
# Function to write the xml of the triggers to an HTML table 
##############################################################################
#class HTMLTable:
#  def __init__(self):
#    self.columns = []
#    self.headers = []
#    
#  def add_column(self,rows,header):
#    self.columns.append(rows)
#    self.headers.append(header)
#
#  def write(self,file):
#    file.write('\n<br><table><tr>')
#    cnt = 0
#    for i in self.columns:
#      file.write('<td><b>'+ self.headers[cnt]+'</b><table>')
#      cnt +=1
#      for j in i:
#        file.write('<tr><td>'+str(j)+'</td></tr>\n')
#      file.write('</table></td>')
#    file.write('\n</tr></table><br>\n')

def getSlots(xml):
  values = []
  for i in xml.__slots__:
    temp = None
    try: temp = getattr(xml,i)
    except: pass
    if temp:
      values.append('<i><font color=blue>'+str(i)+'</i></font>=' + str(temp))
  return values

def writeXMLparams(trig):
  container = HTMLcontainer(trig,__name__)
  print os.getcwd()
  tableFile = open(container.locallink,'w')
  writeIULHeader(tableFile)
  #tableFile = open(container.link,'w')
  tableFile.write('<h3>Follow up of trigger [' +str(trig.eventID) +']</h3>\n')
  container.text = "click here for a link to the xml parameters"
  #if trig.is_trigs():
  #  pass
  table = HTMLTable()
  for ifo in trig.gpsTime:
    if trig.gpsTime[ifo]: 
      xml = getattr(trig.coincs, str(ifo))
      table.add_column(getSlots(xml),str(ifo)) 
  table.write(tableFile)
  tableFile.close()
  return container


  
