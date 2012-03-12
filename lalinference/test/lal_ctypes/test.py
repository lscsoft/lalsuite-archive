#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  test.py
#  
#  Copyright 2012 Ben Aylott <beaylott@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

from ctypes import *

from lalinference import *

def main():
    varss=LALInferenceVariables()

    value=REAL8(90.)
    valuep=pointer(value)
    
    liblalinference.LALInferenceAddVariable(byref(varss),c_char_p("ASD"),valuep,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR)
    ret=liblalinference.LALInferenceGetVariable(byref(varss),c_char_p("ASD"))
    ret2=cast(ret,POINTER(REAL8))
    print ret2.contents.value
    return 0
if __name__ == '__main__':
    main()

