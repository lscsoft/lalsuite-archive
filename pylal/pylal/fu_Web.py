#!/usr/bin/env @PYTHONPROG@
"""
followup Web page classes

$Id$

This
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]


# Think of a web/wiki page as an object that has organized content
# The organization of this content is simple in this case, the content can
# be organized by sections, subsections, or it can be a free for all.
# All content however will be held inside of a table.  There can be of course
# many tables for each web page.  This makes sense since
# any thing one would want to write can be made in a one cell table.  This
# extra layer helps make presentation uniform and coding easier.  
class Content:
  
  def __init__(self):
    self.contentList = []
    self.table = []
    #self.root = ''

  def link(self,link,text):
    thisLink = Link(link,text)
    self.contentList.append(thisLink)

  def verbatim(self,text):
    thisVerbatim = Verbatim(text)
    self.contentList.append(thisVerbatim)

  def text(self,text):
    thisText = Text(text)
    self.contentList.append(thisText)

  def list(self,list):
    thisList = List(list)
    self.contentList.append(thisList)

  def image(self, image, link=None, width = 400):
    thisImage = Image(image, self.root, link, width)
    self.contentList.append(thisImage)

  def linebreak(self, times = 1):
    thisBreak = Break(times)
    self.contentList.append(thisBreak)
 
  def appendTable(self,rows,columns,border=0,width=800):
    number = len(self.table)
    thisTable = Table(rows,columns,border,width,number,self.root)
    self.contentList.append(thisTable)
    self.table.append(thisTable)



#Add sub pages and provide slick mechanism for writing the whole tree!?!
#this will help with the followup dag web generation...
# ADD CONCEPT OF RELATIVE PATH TO GIVE MORE FLEXIBILITY TO DIRECTORY STRUCTURE
class WebPage(Content):
  """ 
  Class to store the web page output of a followup program
  """
  def __init__(self,title,filename,root=''):
    #Content will be written before sections, which themselves may have content
    self.contentList = []
    self.table = []
    self.section = []
    self.title = title
    self.root = root
    self.subPage = []
    self.filename = filename

  def appendSection(self,heading):
    number = len(self.section)
    self.section.append( Section( heading,number,self.root ) )  

  
  def appendSubPage(self,title, file, root=''):
    self.subPage.append(WebPage(title, file, root))

  def linkNewSubPage(self,title,file, text='', root=''):
    if text:
      self.link(root+file,text)
    else:
      self.link(root+file, title)
    self.appendSubPage(title,file, root)


  def write(self,type):
    self.writeHeader(self.file,type)
    self.writeTitle(self.file,type)
    self.writeTableOfContents(self.file,type)
    # first do the content in this page
    for content in self.contentList:
      content.write(self.file,type)
    for section in self.section:
      section.write(self.file,type)
    # now do the sub pages recursively
    for page in self.subPage:
      page.cleanWrite(type)
    

  def writeTitle(self, file, type):
    if type == 'IUL':
      pass # cause this is done in the header for IUL type

  def writeHeader(self, file, type):
    if type == 'IUL':
      file.write('<%method title>' + self.title + '</%method>\n')
      file.write('<%method headline>' + self.title + '</%method>\n')
      file.write('<%method cvsid>$Id$</%method>\n')
 
  def writeTableOfContents(self,file,type):
    if type == 'IUL':
      file.write('<h3 id="toc">Table of contents</h3>\n') 
      sectionTOC  = [] 

      for section in self.section:
        link = section.toclink
        subSectionTOC = []
        for subsection in section.subSection:
          subSectionTOC.append( [Link(subsection.toclink, subsection.heading)])
        if subSectionTOC:
          sectionTOC.append( [Link(section.toclink, section.heading), List(subSectionTOC)] )
        else: 
          sectionTOC.append( [Link(section.toclink, section.heading)] )
      TOC = List(sectionTOC)
      TOC.write(file,type)
          
#  def cleanWrite(self, filename, type):
  def cleanWrite(self,type):
    self.file = open(self.filename,'w')
    self.write(type)          
    self.file.close()

# This class shouldn't really be used without a webpage as it's parent
class Section(Content):
  """
  Class to store a section of a webpage
  """
  def __init__(self,heading,secNumber,root=''):
    self.contentList = []
    self.table = []
    self.subSection = []
    self.heading = heading
    self.secNumber = secNumber 
    self.toclink = root + '#section' + str(self.secNumber)
    self.root = root

  def appendSubSection(self,heading):
    number = len(self.subSection)
    self.subSection.append( SubSection( heading,self.secNumber,number, self.root ) )

  def write(self,file,type):
    self.writeSectionHeader(file,type)
    for content in self.contentList:
      content.write(file,type)
    for subSection in self.subSection:
      subSection.write(file,type)
  
  def writeSectionHeader(self,file,type):
    if type == 'IUL':
      file.write('<h3 id="section'+str(self.secNumber)+'">'+str(self.secNumber)+'.  ' + self.heading+'\n')
      file.write('<a href="'+self.root+'#toc">[Back to TOC]</a></h3>\n')
      
# This class shouldn't really be used without a section as its parent, which
# itself has a webpage as its parent
class SubSection(Content):
  """
  Class to store a subsection of a webpage
  """
  def __init__(self,heading,secNumber,subNumber, root=''):
    self.contentList = []
    self.table = []
    self.heading = heading
    self.secNumber = secNumber
    self.subNumber = subNumber
    self.root = root
    self.toclink = root + '#subsection' + str(self.secNumber) + '.' + str(self.subNumber)

  def write(self,file,type):
    self.writeSubSectionHeader(file,type)
    for content in self.contentList:
      content.write(file,type)

  def writeSubSectionHeader(self,file,type):
    if type == 'IUL':
      file.write('<h4 id="subsection'+str(self.secNumber)+'.'+str(self.subNumber)+'">'+str(self.secNumber)+'.'+str(self.subNumber)+'.  '+self.heading+'\n')
      file.write('<a href="'+self.root+'#toc">[Back to TOC]</a></h4>\n')




# here is where the real 'content' is.  Remember that all content is
# contained inside of a table. Table allows you to create rows, columns
# and cells, there shouldn't be any normal reason to use those classes
# by themselves  - currently this just does tables of form (m x n) 
class Table:
  """
  Class to store the web page output of a followup program
  """
  def __init__(self,rows,columns,border,width,number,root):
    self.row = []
    self.number = number
    self.border = border
    self.width = width
    self.root = root
    for i in range(rows):
      self.row.append(Row(columns,root))

  def write(self,file,type):
    if type == 'IUL':
      file.write('<table border=' + str(self.border)+' width='+str(self.width)+'>\n')
      for row in self.row:
        row.write(file,type)
      file.write('</table><br>\n')

# You shouldn't need to use this class - best to use Table
class Row:
  """
  Class to store a table row
  """
  def __init__(self,columns,root):
    self.cell = []
    self.root = root
    for j in range(columns):
      self.cell.append(Cell(root))

  def write(self,file,type):
    if type == 'IUL':
      file.write('<tr>')
      for cell in self.cell:
        cell.write(file,type)
      file.write('</tr>\n')

# You shouldn't need to use this class - best to use Table
class Cell(Content):
  """
  Class to store a table cell
  """
  def __init__(self,root):
    self.contentList = []
    self.table = []
    self.root = root

  def write(self, file, type):
    if type == 'IUL':
      file.write('<td>')
      for content in self.contentList:
        content.write(file,type)
      file.write('</td>')


# finally something that does something. (but not very much)
class Link:

  def __init__(self, link, text):
    self.link = link
    self.text = text

  def write(self, file, type):
    if type == 'IUL':
      file.write('<a href="' + self.link + '">' + self.text + '</a>')
      
class Verbatim:

  def __init__(self, text):
    self.text = text

  def write(self, file, type):
    if type == 'IUL':
      file.write('\n<br><pre>' + self.text + '</pre><br>\n')


class Text:
  
  def __init__(self,text):
    self.text = text
 
  def write(self, file, type):
    if type == 'IUL':
      file.write('<p>')
      file.write(self.text)
      file.write('</p>\n')

class List:
 
  def __init__(self,list):
    # valid objects in this list are lists of Link, Verbatim, Text, List, Image
    # or any custom object with a compliant write method.
    # remember this expects a list of lists!
    self.list = list
  
  def write(self, file, type):
    if type == 'IUL':
      file.write('<ol>\n')
      for item in self.list:
        file.write('<li>')
        for element in item:
          element.write(file,type)
        file.write('\n')
      file.write('</ol>\n')

class Image:

  def __init__(self,image,root,link=None,width=400):
    self.link = link
    self.image = root+image
    self.width = width
  def write(self, file, type):
    if type == 'IUL':
      if self.link:
        file.write('<a href="'+self.link+'"><img src="'+self.image+'" width = '+str(self.width)+'></a>')
      else: 
        file.write('<img src="'+self.image+'" width = '+str(self.width)+'>')   

class Break:

  def __init(self, times = 1):
    self.times = range(times)

  def write(self, file, type):
    if type == 'IUL':
      for time in self.times:
        file.write('<br>')
    
## TEST PROGRAM
###########################

testPage = WebPage("This web page made itself",'index.html',"http://www.lsc-group.phys.uwm.edu/iulgroup/investigations/s5/people/channa/fu_Web/")
testPage.linkNewSubPage("This web page made itself",'index2.html',"click here!","http://www.lsc-group.phys.uwm.edu/iulgroup/investigations/s5/people/channa/fu_Web/")
testPage.appendSection("What is up with fu_Web.py?")
testPage.section[0].appendSubSection("motivation for fu_Web.py")
testPage.section[0].subSection[0].text("fu_Web.py is a module in pylal that allows for a simple web output scheme.  It provides a standard way for a developer to produce simple web pages within their codes - for example this page")
testPage.section[0].appendSubSection("theory of fu_Web.py")
testPage.section[0].subSection[1].text("We need a simple formalism to abstract simple web design to make presenting information easier.  fu_Web.py provides organization with the hope that it will facilitate the conveyance of important results")
testPage.section[0].subSection[1].text("Web.py is made to allow flexible write methods that can turn the same code into HTML, Wiki, (Latex - if your a masochist) [of course right now it is only HTML]")
testPage.section[0].appendSubSection("WHERE DO I GET IT ?!!")
testPage.section[0].subSection[2].link("http://www.lsc-group.phys.uwm.edu/cgi-bin/cvs/viewcvs.cgi/pylal/pylal/fu_Web.py?cvsroot=lscsoft", "TRY HERE!")

testPage.appendSection("Using fu_Web.py")

testPage.section[1].appendSubSection("It's good for more than just follow-ups!")

testPage.section[1].subSection[0].text("Everyone who writes python based code to produce end data products like plots and tables that belong on the web can use this code")
testPage.section[1].appendSubSection("Sections, subSections")
testPage.section[1].subSection[1].text("Hopefully you will all choose to organize your work based on sections and subsections.  If you do the table of contents and links will be automatically generated when the write method is called")
testPage.section[1].appendSubSection("What can I do with this")
testPage.section[1].subSection[2].appendTable(3,2,1)
testPage.section[1].subSection[2].table[0].row[0].cell[0].text("Well you could try a table")
testPage.section[1].subSection[2].table[0].row[0].cell[1].link(testPage.root+"#toc","or a link")
testPage.section[1].subSection[2].table[0].row[1].cell[0].list([[Text("or a list")],[Link(testPage.root+"#toc","of links")]])
testPage.section[1].subSection[2].table[0].row[1].cell[1].image("L1_793808911_704_snr.png")

verbatim = '''
testPage = WebPage("This web page made itself",
  "http://www.lsc-group.phys.uwm.edu/iulgroup/investigations/s5/people/channa/fu_Web/")
testPage.appendSection("What is up with fu_Web.py?")
testPage.section[0].appendSubSection("motivation for fu_Web.py")
......
testPage.section[1].subSection[2].appendTable(3,2,1)
testPage.section[1].subSection[2].table[0].row[0].cell[0].text("Well you could try a table")
......
testPage.cleanWrite('index.html', 'IUL')
'''
testPage.section[1].subSection[2].table[0].row[2].cell[0].text("Some verbatim?")
testPage.section[1].subSection[2].table[0].row[2].cell[1].verbatim(verbatim)


testPage.appendSection("Next Steps")
testPage.section[2].text("I hope that people will begin to code using this set of classes, and furthermore will start to improve it.")
testPage.section[2].text("The next major step is to integrate this into the followup pipeline.  I'll need to write a parser/metadata output generator for the pipeline so that each step can talk back to the code that writes the final web page")

testPage.cleanWrite('IUL')
