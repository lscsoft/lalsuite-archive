#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       lalinference.py
#
#       Copyright 2011
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#===============================================================================
# Preamble
#===============================================================================

"""
This module exposes the data structures used in LALInference in convenience
classes.
"""

from pylal import git_version

from pylal._lalinference import BaseLALVariables,BaseLALIFOData,BarycenterInput,PosVelAcc,EphemerisData

__author__="Ben Aylott <benjamin.aylott@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

#===============================================================================
# Class definitions
#===============================================================================

class LALVariables(BaseLALVariables):
    """
    This class dresses what is found BaseLALVariables with special python
    methods. It is the preferred method of interacting with LALVariables
    data structures. 
    """
    def __init__(self,initDict=None):
        BaseLALVariables.__init__(self)

        if initDict is not None:
            self.addVariables(initDict)

    def __getitem__(self,key):
        return self.getVariable(key)

    def __len__(self):
        return self.getVariableDimension()

    def addVariables(self,inputDict):

        for name,data in inputDict.items():
            val,varytype=data
            self.addVariable(name,val,varytype)
        
class LALIFOData(BaseLALIFOData):
    """
    This class dresses what is found BaseLALIFOData with special python
    methods. It is the preferred method of interacting with LALIFOData
    data structures. 
    """
    def __init__(self):
        BaseLALIFOData.__init__(self)

    

if __name__ == '__main__':
    #Run unit tests
    pass
