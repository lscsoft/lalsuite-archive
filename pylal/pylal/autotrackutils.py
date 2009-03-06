#!/usr/bin/env python
"""
  * Copyright (C) 2004, 2005 Cristina V. Torres
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with with program; see the file COPYING. If not, write to the
  *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  *  MA  02111-1307  USA
"""
__author__ = 'Cristina Torres <cristina.torres@ligo.org>'
__date__ = '$Date$'
__version__ = ''

import os
import numpy
import copy
import sys
import sqlite3 as sqlite

if os.getenv("DISPLAY") == None:
    #Non-interactive
    try:
        import matplotlib
        matplotlib.use("Agg")
        import pylab
    except Exception, errorInfo: #RuntimeError,ImportError:
        disableGraphics=True
        sys.stderr.write("Error trying to import NON-INTERACTIVE pylab!\n")
        sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
        sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
        sys.stderr.write("Pylab functionality unavailable!\n")
else:
    #Interactive
    try:
        import pylab
    except Exception, errorInfo: #RuntimeError,ImportError:
        disableGraphics=True
        sys.stderr.write("Error trying to import INTERACTIVE pylab!\n")
        sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
        sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
        sys.stderr.write("Pylab functionality unavailable!\n")

"""
This file is intended to provide the autotrack utilities required to
process the results of the tracksearch hybrid MDGC codes.
"""


def generateNhoodListFromFile(filename="",timestamp=str("0")):
    fp=open(filename,'rt')
    rawData=fp.readlines()
    fp.close()
    linesPerNhood=2
    if rawData.__len__().__mod__(linesPerNhood) != 0:
        raise autotrackError("File appears inconsistent!")
    eCount=rawData.__len__().__div__(linesPerNhood)
    nhoodList=list()
    for index in range(0,eCount):
        thisNhood=nhood()
        thisNhood.createFromString(rawData[index*2],rawData[index*2+1])
        thisNhood.__setBirthDate__(timestamp)
        nhoodList.append(thisNhood)
    return nhoodList
    #End

class autotrackError(Exception):
    """
    Generic class raise this flag if we do not have a better
    flag name.
    """
    #End class

class autotrackDefineMismatchError(Exception):
    """
    Custom error exceptions for these class objects.
    """
    #End class autotrackDefineError

class autotrackDefineIdMismatchError(Exception):
    """
    Custom error handling of mismatched ID
    """
    #End class autotrackDefineIdMismatchError

class nhood:
    """
    This class provides the definition of a single defined autotrack
    defined neighborhood.  We use these neighborhoods to track
    how the instrument behavior groups change, appear or disappear.
    """
    def __init__(self):
        """ This method initializes an empty neighborhood.
        To populate the neightborhood invoke the proper method
        depending on the source of the data, ie ascii mysql etc
        """
        self.idNum=float(-1)
        self.density=float(0)
        self.memberCount=0
        self.volume=0
        self.colCount=0
        self.birthdate=str("-1")
        self.discoverer=str()
        self.lastSeen=float(0)
        self.center=None
        self.bound=None
    #End __init__ method

    def __setID__(self,inputArg):
        """
        Should set the ID numeric field
        """
        if type(float(0)) != type(inputArg):
            inputArg=float(inputArg)
        self.idNum=inputArg
    #End

    def getID(self):
        """
        Fetches the numeric ID assigned to this nhood instance
        """
        return self.idNum
    #End

    def __setBirthDate__(self,bDate=str("")):
        """
        Set a text string which is the GPS birthdate,
        as closely as possible for this neighborhood.
        """
        if type(str()) != type(bDate):
            bDate=str(bDate)
        self.birthdate=bDate
    #End

    def getBirthDateText(self):
        """
        This retrieves the text string birthdate stored in nhood
        instance.
        """
        return self.birthdate
    #End

    def createFromString(self,mString,bString,delim=None):
        """
        This input method assumes you've opened a text file container
        and will input the data from text strings.  The assumed
        delimiter is a space but can be set in the method call.
        Old file spec had 17 cols new has 19 cols
        """
        mData=str(mString).lower().split(delim)
        bData=str(bString).lower().split(delim)
        mCol=mData.__len__()
        sCol=bData.__len__()
        if (mData.__len__()) != (bData.__len__()):
            raise autotrackDefineMismatch("Array lengths %i:%i"%(mData.__len__(),bData.__len__()))
        self.colCount=mCol
        #Break off first few elements before creating arrays
        mID=mData.pop(0)
        mDensity=mData.pop(0)
        if mCol >= 17:
            mCount=mData.pop(0)
            mVolume=mData.pop(0)
        #
        sID=bData.pop(0)
        sDensity=bData.pop(0)
        if sCol >=17:
            mCount=mData.pop(0)
            mVolume=mData.pop(0)
        if mID != sID:
            raise autotrackDefineIdMismatchError("Group labels do not match!")
        if mDensity != sDensity:
            raise autotrackDefineIdMismatchError("Group density values do not match!")
        if mCount != sCount:
            raise autotrackDefineIdMismatchError("Group count values do not match!")
        if mVolume==sVolume:
            raise autotrackDefineIdMismatchError("Group volume measures do not match!")
        self.__setID__(float(mID))
        self.__setDensity__(float(mDensity))
        self.__setMemberCount__(float(mCount))
        self.__setVolume__(float(mVolume))
        self.center=numpy.array(mData,'float64')
        self.bound=numpy.array(bData,'float64')
    #End createFromString

    def __setMemberCount__(self,mCount=0):
        """
        Sets the tally of members in the group.
        """
        self.memberCount=mCount
    #End

    def getMemberCount(self):
        """
        get the registered members listed for this nhood
        """
        return self.memberCount
    #End

    def __setVolume__(self.volume=0):
        """
        Sets the volume value of this grouping.
        """
        self.volume=volume
    #End

    def getVolume(self):
        """
        Gets the registered volume for a given nhood.
        """
        return self.volume
    #End

    def exportToString(self):
        """
        Create a text string that can be directly inserted into the
        text field of the sqlite database table TGN.
        """
        delimiter=":"
        outputString=""
        outputString=outputString+":%s"%(self.getID())
        outputString=outputString+":%s"%(self.getDensity())
        if self.colCount >= 17:
            outputString=outputString+":%s"%(self.getMemberCount)
            outputString=outputString+":%s"%(self.getVolume)
        for elem in self.center:
            outputString=outputString+":%s"%(elem)
        outputString=outputString+"\n"
        outputString=outputString+":%s"%(self.getID())
        outputString=outputString+":%s"%(self.getDensity())
        if self.colCount >= 17:
            outputString=outputString+":%s"%(self.getMemberCount)
            outputString=outputString+":%s"%(self.getVolume)
        for elem in self.bound:
            outputString=outputString+":%s"%(elem)
        return outputString
    #end exportToString()
    def getDensity(self):
        """
        Returns the value of nhood density set.
        """
        return self.density
    #End

    def __setDensity__(self,inputArg):
        """
        Sets the input density value to the nhood instance.
        """
        if type(float(0)) != type(inputArg):
            inputArg=float(inputArg)
        self.density=inputArg
    #End

    def isNULL(self):
        """
        Check the defined properties, returns true if they are all
        zeroed out which is NULL according the the matlab generator.
        """
        isNull=bool(False)
        cV=self.getCenterVector()
        bV=self.getBoundVector()
        if numpy.any(cV==0) and numpy.any(bV==0):
            isNull=bool(True)
        return isNull
    #End isNull

    def getBoundVector(self):
        """
        Gets the variance components of this nhood instance.
        """
        return self.bound
    #End getBoundVector

    def getCenterVector(self):
        """
        Gets the center of the neighborhood
        """
        return self.center
    #End getCenterVector

    def isSame(self,nhood):
        """
        Checks to see if self instance is IDENTICAL to 
        nhood instance given as argument!
        """
        samePoint=bool(False)
        samePoint=self.checkOverlap(nhood,0)
        if samePoint and self.idNum==nhood.idNum and self.density==nhood.density:
            return bool(True)
        return bool(False)
            
    def checkOverlap(self,nhood,boundSize=0.5):
        """
        Check to see if SELF neighborhood overlaps with other input
        NHOOD class.  We define them as over lapping if they are
        withing boundSize stddevs of the center of SELF compared
        to the center of argument NHOOD.
        """
        stepVector=boundSize*numpy.array(self.getBoundVector()).__abs__()
        diffVector=numpy.array(
            self.getCenterVector()
            -nhood.getCenterVector()).__abs__()
        if numpy.array(diffVector<=stepVector).all():
            return bool(True)
        else:
            return bool(False)
        #End checkOverlap

    def getSeparation(self,nhood):
        """
        Gets the resultant seperation vectors and normalizes this 
        value by the boundVector, then use this to compute a 
        normalized magnitude of the vector.
        """ 
        diffVector=numpy.array(
            self.getCenterVector()
            -nhood.getCenterVector()).__abs__()
        bv=self.getBoundVector()
        sepVector=diffVector.__div__(bv)
#         if numpy.isnan(sepVector).any():
#             print diffVector
#             print bv
#             print sepVector
        mag=numpy.sqrt(numpy.inner(sepVector,sepVector))
        return mag
    #End getSeperation
    
    def nearestNhoodID(self,nhoodList,boundSize=0.5):
        """
        wrapper
        """
        myN=self.nearestNhood(nhoodList,boundSize)
        myID=myN.getID()
        return myID

    def nearestNhood(self,nhoodList,boundSize=0.5):
        """
        Takes a list of nhoods and compares it to self
        to determine which is the closest one, this ideally
        is the same group and we can associate the group IDs.
        """
        if type(list())!=type(nhoodList):
            raise autotrackError("Type of input to method nearestNhoodID is wrong!")
        distanceKey=list()
        for index in range(0,nhoodList.__len__()):
            #Ignore entries in the list that all NULL
            if not nhoodList[index].isNULL():
                if self.checkOverlap(nhoodList[index],boundSize):
                    dist=self.getSeparation(nhoodList[index])
                    idVal=nhoodList[index].getID()
                    distanceKey.append([dist,idVal])
        distanceKey.sort()
        try:
            findID=int(distanceKey[0][1])
        except IndexError:
            findID=-1
        if findID > -1:
            for nhd in nhoodList:
                if nhd.getID() == findID:
                    foundNhood=nhd
        else:
            foundNhood=nhood()
        return foundNhood
    #End nearestNhood
#End class nhood
        
class autotrackSQL:
    """
    This class provides sqlite table creation query deletion etc,
    related funtions for working the the autotrack databases.  This
    """
    def __init__(self,dbName="autotrack_default.sqlite"):
        """
        Initializes the variables associated with the autotrackSQL
        database manipulations. Setup is {TableName,Bool,{Definition}}
        """
        self.dbFile=""
        self.dbSocket=None
        self.defineTables=(
            ('tgn',bool(True),
             ('group_serial_number',
              'group_birthdate',
              'discoverer',
              'group_label',
              'group_density',
              'group_member_count',
              'group_last_seen',
              'statistics')
             ),
            ('scientist_entry',bool(True),
             ('group_serial_number',
              'entry_date',
              'scientist',
              'channels_of_interest',
              'URLs_of_interest',
              'description_text',
              'solution_text',
              'misc_information',
              'extra_field')
             ),
            ('plot_locations',bool(False),
             ('group_serial_number',
              'line_plot',
              'aux_plots')
             ),
            ('sightings',bool(True),
             ('group_serial_number',
              'parent_group_serial_number',
              'earliest_sighting',
              'sighting_list')
             ),
            ('historical_tgn',bool(False),
             ('group_serial_number',
              'group_birthdate',
              'discoverer',
              'group_label',
              'group_density',
              'group_member_count',
              'group_last_seen',
              'statistics')
             )
            )
        #Look for the table if it does not exist create it
        #otherwise load it into this object!
        self.__selectDB__(dbName)
        #End __init__()

    def __selectDB__(self,dbName='autotrack_default.sqlite'):
        """
        Selects the specified if it exists is reads the contents
        or if not is issues a warning and tells you the db 
        needs to be created.
        """
        #Code up how to check for a db and either create the tables
        #or read them
        #End self.__selectDB__()
        self.dbSocket=sqlite.connect(dbName)
        self.dbFile=dbName
        #Try seeing if there are tables that we want in there already
        try:
            self.dbSocket.execute("select * from %s"%(self.defineTables[0][0]))
        except:
            sys.stdout.write("It appears that we need to create the tables.\n")
            self.createTables()
        #End __selectDB__()

    def getTableIndex(self,name=""):
        """
        Given a string searches the tuple self.defineTables to
        determine the index of that table to various function calls.
        """
        answer=int(-1)
        for index in range(self.defineTables.__len__()):
            if self.defineTables[index][0].lower() == name.lower():
                answer=index
        if not(type(int()) == type(answer)):
            raise autotrackError("getTableIndex type invalid!\n %s"(answer))
        return answer
    #End getTableIndex

    def __createSingleTable__(self,name):
        """
        This method will generate a table from the list of tables
        and definitions specified by self.autotrackTableDef or
        self.auxTableDef
        """
        tableIndex=self.getTableIndex(name)
        thisTable=self.defineTables[tableIndex]
        tableName=thisTable[0]
        required=thisTable[1]
        colDefs=thisTable[2]
        rowString=""
        for col in colDefs:
            rowString=rowString+"%s "%(col)
        commandString="create table %s (%s)"%(tableName,rowString)
        try:
            self.dbSocket.execute(commandString)
        except:
            sys.stderr.write("Error trying to create table.\n")
            self.dbSocket.commit()
            return
        self.dbSocket.commit()
        #End __createSingleTable__()
        
    def getSocket(self):
        """
        Returns an object which is a raw handle to the dqlite DB.
        This is equivalent to being returned X from
        X=sqlite.connect('dbfilename')
        """
        return self.dbSocket
        #End getSocket

    def createTables(self,overwrite=bool(False)):
        """
        This function call will create the specified sqlite db unless
        there exists one already.  In that case we throw an error
        unless we explicitly want to overwrite that table.
        """
        for index in range(self.defineTables.__len__()):
            thisTableName=self.defineTables[index][0]
            sys.stdout.write("Creating table %s\n"%(thisTableName))
            self.__createSingleTable__(thisTableName)
        self.dbSocket.commit()
        #End createTables()
