from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ilwd
from glue.ligolw import types

################################################################################
# Definitions of some of the classes and functions that are used in the main
# program.
################################################################################
class InspiralLikelihoodTable(table.Table):
  tableName = "Inspiral_Likelihood:table"
  validcolumns = {
    "ifo": "lstring",
    "false_alarm_prob": "real_4",
    "detect_prob": "real_4",
    "likelihood": "real_4",
    "bg_frac": "real_4",
    "detect_bg_frac_ratio": "real_4",
    "event_id": "ilwd:char"
  }
  def get_column(self,column):
    return self.getColumnByName(column).asarray()

  def getslide(self,slide_num):
    """
    Return the triggers with a specific slide number.
    @param slide_num: the slide number to recover (contained in the event_id)
    """
    slideTrigs = lsctables.New(InspiralLikelihoodTable)
    slideTrigs.extend([r for r in self if r.get_slide_number() == slide_num])
    return slideTrigs

  def coinctype(self, ifo):
    """
    Return the events with the given IFO.
    @param ifo: parameter specifying the IFO ( e.g. "H1H2L1")  
    """
    ifo_triggers=lsctables.New(InspiralLikelihoodTable)

    for row in self:
      if row.ifo == ifo:
        ifo_triggers.append(row)

    return ifo_triggers

class InspiralLikelihood(object):
  __slots__ = InspiralLikelihoodTable.validcolumns.keys()

  def get_id_parts(self):
    """
    Return the three pieces of the int_8s-style sngl_inspiral
    event_id.
    """
    int_event_id = int(self.event_id)
    x = int_event_id // 1000000000
    slidenum = (int_event_id % 1000000000) // 100000
    y = int_event_id % 100000
    return x, slidenum, y

  def get_slide_number(self):
    """
    Return the slide-number for this trigger
    """
    x, slidenum, y = self.get_id_parts()
    slide_number = slidenum
    if slide_number > 5000:
      slide_number = 5000 - slide_number
    return slide_number

InspiralLikelihoodTable.RowType = InspiralLikelihood


def ReadInspiralLikelihoodFromFiles(fileList):
  """
  Read the InspiralLikelihoodTables from a list of files
  @param fileList: list of input files
  """
  FullTable = lsctables.New(InspiralLikelihoodTable)
  for file in fileList:
    doc = utils.load_url(file, verbose = False, gz = file[-3:] == ".gz",
        xmldoc = ligolw.Document())
    try:
      SubTable = ligolw.table.get_table(doc, InspiralLikelihoodTable.tableName)
    except:
      SubTable = None
    if SubTable:
      for row in SubTable:
        FullTable.append(row)
  return FullTable

def generate_prefix_and_suffix(ifo_times="", program_name="", ifo_tag="",
  user_tag="", gps_start_time=0, gps_end_time=0):
  """
  Generates the string according to inspiral pipeline conventions, that can
  be used as a name of a file.
  """
  prefix = program_name
  if ifo_times:
    prefix = ifo_times + "-" + prefix
  if ifo_tag:
    prefix = prefix + "_" + ifo_tag
  if user_tag:
    prefix = prefix + "_" + user_tag

  if gps_start_time > 0 and gps_end_time > 0:
    suffix = "_" +str(gps_start_time) + "-" + \
      str(gps_end_time - gps_start_time)
  else:
    suffix = ""
  return prefix, suffix

