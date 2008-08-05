#!/usr/bin/env @PYTHONPROG@
#
"""
followup web page scraping utilities

$Id$

"""

__author__ = 'Cristina Valeria Torres <cristina.torres@ligo.org>'
__date__ = '$Date: 2008/08/15 14:00'
__version__ = '$Revision$'

######################################################################

import os
import string
import sys
import time

class scrapePage:
    """
    This class is responisble for taking in out expected html
    formatted file and allowing us to manipulate the central table of
    interest while keeping the rest of the html available for later
    writing to a disk.
    """
    def __init__(self):
        self.filename=''
        self.fileRead=False
        self.originalText=list()
        self.startKey="<h3>Follow-up tests</h3>\n"
        self.endKey="<h3>Parameter estimation</h3>\n"
        self.saveLines=list()
        self.topOfPage=list()
        self.middleOfPage=list()
        self.endOfPage=list()
        self.tableObject=list()
        self.rowNames=list()
        self.colNames=list()
        self.removeKey=["<table","<tbody"]
        ignoreKeysMatch=list()
        for match in self.removeKey:
            ignoreKeysMatch.append("/"+match.strip("<"))
        self.removeKey.extend(ignoreKeysMatch)
    #End Init

    def setContextKeys(self,newStartKey="",newEndKey=""):
        """
        Calling self.setContectKeys will allow you to specify two new
        context keys to select a single table from a parsed html
        file.  The two arguments for this function require that you
        specify a key which is one entire line long from the source
        html file that you want to extract the table from.  This will
        allow the code to save the surrounding html and allow you to
        manipulate the table in a more natural manner of
        table[row][col] to edit the entries.  Do not set the keys to
        partial line matches or non-printing characters, this will
        almost ensure the failure of the html parser.
        """
        if newStartKey != "":
            self.startKey=newStartKey
        if newEndKey != "":
            self.endKey=newEndKey
    #End setContextKeys()

    def readfile(self,filename):
        """
        Reads in a text html file given a filename.
        """
        fp=open(filename)
        self.originalText=fp.readlines()
        fp.close()
        startCut=0
        endCut=0
        foundStartKey=False
        foundEndKey=False
        currentLine=0
        while (currentLine < self.originalText.__len__()):
            if not(foundStartKey):
                self.topOfPage.append(self.originalText[currentLine])
            if foundStartKey and not(foundEndKey):
                self.middleOfPage.append(self.originalText[currentLine])
            if foundStartKey and foundEndKey:
                self.endOfPage.append(self.originalText[currentLine])
            if self.originalText[currentLine].__contains__(self.startKey):
                foundStartKey=True
            if self.originalText[currentLine].__contains__(self.endKey):
                foundEndKey=True
                self.endOfPage.append(self.originalText[currentLine])
            currentLine=currentLine+1
        cleantext=list()
        for line in self.middleOfPage:
            foundKey=0
            for key in self.removeKey:
                if line.lower().__contains__(key):
                    foundKey=foundKey+1
            if foundKey == 0:
                cleantext.append(line)
            else:
                self.saveLines.append(line)
        text=str().join(cleantext)
        tableObject=list()
        for row in text.replace("<tr","<MARK><tr").split("<MARK>"):
            self.tableObject.append(row.replace("<td","<MARK><td").split("<MARK>"))
        rowNames=list()
        rowNumber=0
        for row in self.tableObject:
            if row.__len__() > 2:
                self.rowNames.append([rowNumber,row[1]])
                rowNumber=rowNumber+1
            else:
                self.rowNames.append([rowNumber,"\n"])
                rowNumber=rowNumber+1
        colNumber=0
        for col in range(0,self.tableObject[1].__len__()):
            self.colNames.append([colNumber,self.tableObject[1][col]])
            colNumber=colNumber+1
#End read method

    def getColumnByText(self,textString='',colNum=1):
        """
        Given a text string expected in Column #1 we select the
        specified column given as an argument here. If there was
        nothing found return empty string.
        """
        currentRow=0
        rowCount=self.rowNames.__len__()-1
        foundRow=-1
        while currentRow <= rowCount:
            if self.__compareKeyWords__(textString.lower(),self.rowNames[currentRow][1].lower()):
                foundRow=self.rowNames[currentRow][0]
                currentRow=rowCount+1
            currentRow=currentRow+1
        if (foundRow > -1):
            outputData=self.tableObject[foundRow][colNum]
            return outputData            
        else:
            return""
    #End getColumnByText()

    def showRows(self):
        """
        Call this method after reading the html to literally see 
        the row labels inside the HTML table we are manipulating.
        """
        for row in self.rowNames:
            sys.stdout.write("Row %i, %s"%(int(row[0]),str(row[1])))
            sys.stdout.flush()
     #End showRows()

    def getRowList(self):
        """
        This method gets the list of rows in the table for that 
        htmlPage() instance.  The data returned in a list of two element
        lists. Like [[a,b],[c,d],...,[y,z]]
        """
        return self.rowNames

    def showCols(self):
        """
        Call this method after reading the html to literally see
        the column labels inside of the html table we are manipulating.
        """
        colNum=1
        for col in self.colNames:
            sys.stdout.write("Col %i, %s"%(int(col[0]),str(col[1])))
            sys.stdout.flush()

    def getColumnByCoord(self,RowNum,ColNum):
        """
        Given a row number and column number return that element in
        the table. If the coords do not exist return empty string.
        """
        return self.tableObject[RowNum][ColNum]
    #End getColumnByCoord()

    def insertTextAtCoord(self,RowNum,ColNum,Text):
        """
        Given a row number and column number insert the argument text
        over what currently exists.  If the RowNum and ColNum is out
        of bounds do nothing.
        """
        self.tableObject[RowNum][ColNum]=Text
    #End insertTextAtCoord

    def insertTextGivenText(self,matchText,colNum,Text):
        """
        Looks for given row matching column 1 to given text.  It then
        inserts the Text into the column specified by ColNum.  If
        there is no match or ColNum is out of bound nothing is done.
        """
        currentRow=0
        rowCount=self.rowNames.__len__()-1
        foundRow=-1
        while currentRow <= rowCount:
            if self.__compareKeyWords__(matchText.lower(),self.rowNames[currentRow][1].lower()):
                foundRow=self.rowNames[currentRow][0]
                currentRow=rowCount+1
            currentRow=currentRow+1
        if foundRow > -1:
            self.tableObject[foundRow][colNum]=Text
    #End insertTextGivenText()

    def __buildMiddleOfPage__(self):
        """
        This method should not be called explicity.  It will rebuild
        the table object variable into a chunk of html for writing to
        the disk.
        """
        tmpMiddle=list()
        saveCount=self.saveLines.__len__()/2
        topSave=list()
        bottomSave=list()
        for i in range(0,saveCount-1):
            topSave.append(self.saveLines[i])
        for i in range(saveCount,self.saveLines.__len__()):
            bottomSave.append(self.saveLines[i])
        middleTop=str().join(topSave)
        middleMid=str().join([str().join(x) for x in self.tableObject])
        middleBot=str().join(bottomSave)
        self.middleOfPage=list()
        self.middleOfPage.append(str().join([middleTop,middleMid,middleBot]))
    #End __buildMiddleOfPage__()

    def __stripHTMLTags__(self,stringIN):
        """
        Take input string and remove all tags inside of < >
        delimiters.
        """
        leftD="<"
        rightD=">"
        input=stringIN
        result=''
        maxloop=0
        ignoreKeys=["http","em","br","h1","h2","h3","hr"]
        ignoreKeysMatch=list()
        for match in ignoreKeys:
            ignoreKeysMatch.append("/"+match.strip("<"))
        ignoreKeys.extend(ignoreKeysMatch)
        foundKeys=0
        output=list()
        while ((input.__contains__("<") and input.__contains__(">")) and (maxloop < 100)):
            maxloop=maxloop+1
            tag=input.__getslice__(input.find("<"),input.find(">")+1)
            foundKeys=0
            for key in ignoreKeys:
                if tag.lower().__contains__(key):
                    foundKeys=foundKeys+1
            if (not(tag.__contains__(" ")) and foundKeys == 0):
                output.append(input.split(tag,1)[0])
                input=input.split(tag,1)[1]
            if (not(tag.__contains__(" ")) and foundKeys > 0):
                output.append(input.split(tag,1)[0])
                input=input.split(tag,1)[1]
        output.append(input)
        result=str().join(output)
        return result
    #End __stripHTMLTags__()

    def __stripRowNumber__(self,stringA):
        """
        Takes the string representing the table row number.  It strips
        the number strip from the front. The input string is assumed
        to have the form #?? Word Words More Words
        where the only number is #??  Ideally this method should only
        be called by self.__compareKeyWords__()
        """
        delimiter="#"
        if (stringA.find(delimiter) == -1):
            return stringA
        [startTXT,middleTXT]=stringA.split("#",1)
        middleTXT=middleTXT.split(" ",1)[1]
        return startTXT+middleTXT
        
    def __compareKeyWords__(self,stringA="",stringB="",exact=False):
        """
        Break stringA into keywords minus html tags.  Then take these
        words and make sure they exist inside of stringB.
        If the exact key is True then strings like
        Big blue bird     will not match Big blue pretty bird
        if the string is left as default (False) then we allow the
        above string to be matched since all the words in the first
        string are contained in the second string.
        """
        if (
            (self.__stripHTMLTags__(self.__stripRowNumber__(stringA)).isspace())
            or
            (self.__stripHTMLTags__(self.__stripRowNumber__(stringB)).isspace())
            ):
            return False
        keyWordList=self.__stripHTMLTags__(self.__stripRowNumber__(stringA)).lower().split()
        matchCount=0
        match=False
        stringB=self.__stripHTMLTags__(self.__stripRowNumber__(stringB)).lower()
        for word in stringB.split():
            for key in keyWordList:
                if (word.__contains__(key)):
                    match=True
            if match:
                matchCount=matchCount+1
            match=False
        if exact:
            if (
                (matchCount==keyWordList.__len__()
                 and
                 (keyWordList.__len__() == list(stringB.split()).__len__())
                 )):
                return True
            else:
                return False
        if matchCount>=keyWordList.__len__():
            return True
        else:
            return False
    #End __compareKeyWords__()

    def writeHTML(self,filename):
        """
        Writes out the html that was manipulated to the file filaname.
        """
        fp=open(filename,'w')
        outputData=list()
        self.__buildMiddleOfPage__()
        outputData.extend(self.topOfPage)
        outputData.extend(self.middleOfPage)
        outputData.extend(self.endOfPage)
        fp.writelines(outputData)
        fp.close()
    #End writeHTML()

#End CLASS scrapPage

                  
