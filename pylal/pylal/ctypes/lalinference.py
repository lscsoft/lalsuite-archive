#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  lalinference.py
#  
#  Part of the ctypes wrapper library for LAL.
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

import pylal.ctypes

from pylal.ctypes.datatypes.primitives import *
from pylal.ctypes.datatypes.complex import COMPLEX8,COMPLEX16
from pylal.ctypes.datatypes.vector import REAL8Vector,UINT4Vector,LIGOTimeGPSVector
from pylal.ctypes.datatypes.real8timeseries import REAL8TimeSeries
from pylal.ctypes.datatypes.real8frequencyseries import REAL8FrequencySeries
from pylal.ctypes.datatypes.complex16frequencyseries import COMPLEX16FrequencySeries
from pylal.ctypes.datatypes.complex16timeseries import COMPLEX16TimeSeries
from pylal.ctypes.datatypes.window import REAL8Window
from pylal.ctypes.datatypes.fftplan import REAL8FFTPlan
from pylal.ctypes.datatypes.laldetector import LALDetector
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS

from pylal.ctypes.utils import make_enum_typedef,_set_types,_set_xlal_types

VARNAME_MAX=1024
   
#LALInferenceParamVaryType

globals().update(make_enum_typedef(
    "LALInferenceParamVaryType",
    [
        "LALINFERENCE_PARAM_LINEAR",
        "LALINFERENCE_PARAM_CIRCULAR",
        "LALINFERENCE_PARAM_FIXED",
        "LALINFERENCE_PARAM_OUTPUT"
    ]
))

#LALInferenceVariableType

class gslMatrix(Structure):
    pass
    
gsl_matrix=gslMatrix
    
class string(c_char_p):
    pass

class void_ptr(c_void_p):
    pass

LALInferenceVariableType_list=[
    "INT4",
    "INT8",
    "UINT4",
    "REAL4", 
    "REAL8", 
    "COMPLEX8", 
    "COMPLEX16", 
    "gslMatrix",
    "REAL8Vector",
    "UINT4Vector",
    "string",
    "void_ptr"
]

class mapped_enum_typedef(c_uint):
    def __init__(self,enum_value,litype):
        c_uint.__init__(self,enum_value)
        self.litype=litype

def __make_enum_litype_mapped_typedef(enum_names):
    vdict={}
    new_enum_typedef=type("LALInferenceVariableType",(mapped_enum_typedef,),{})
    vdict["LALInferenceVariableType"]=new_enum_typedef
    i=0
    for enum_name in enum_names:
        vdict["LALINFERENCE_"+enum_name+"_t"]=new_enum_typedef(i,eval(enum_name))
        i+=1
    globals().update(vdict)
    
__make_enum_litype_mapped_typedef(LALInferenceVariableType_list)

#LALInferenceVariableType

globals().update(make_enum_typedef(
    "LALInferenceDomain",
    [
        "LALINFERENCE_DOMAIN_TIME", 
        "LALINFERENCE_DOMAIN_FREQUENCY"
    ]
))

#LALInferenceDomain

globals().update(make_enum_typedef(
    "LALInferenceDomain",
    [
        "LALINFERENCE_DOMAIN_TIME", 
        "LALINFERENCE_DOMAIN_FREQUENCY"
    ]
))

#LALInferenceApplyTaper

globals().update(make_enum_typedef(
    "LALInferenceApplyTaper",
    [
        "LALINFERENCE_TAPER_NONE",
        "LALINFERENCE_TAPER_START",
        "LALINFERENCE_TAPER_END",
        "LALINFERENCE_TAPER_STARTEND",
        "LALINFERENCE_TAPER_NUM_OPTS",
        "LALINFERENCE_RING",
        "LALINFERENCE_SMOOTH"
    ]
))

#LALInferenceRunState
class LALInferenceRunState(Structure): pass

#LALInferenceVariableItem
class LALInferenceVariableItem(Structure): pass

#LALInferenceVariables
class LALInferenceVariables(Structure):
    _fields_ = [
        ("dimension",INT4),
        ("head",POINTER(LALInferenceVariableItem))
    ]
    
    def __init__(self):
        Structure.__init__(self)
        self.ptr=pointer(self)
        
    def addVariable(self,name,value,li_variable_type,li_param_vary_type):
        value=li_variable_type.litype(value)
        valuep=pointer(value)
        value_voidp=cast(valuep,c_void_p)
        LALInferenceAddVariable(byref(self),c_char_p(name),value_voidp,li_variable_type,li_param_vary_type)
    
    def removeVariable(self,name):
        LALInferenceRemoveVariable(self.ptr,c_char_p(name))
    
    def getVariable(self,name):
        
        get=LALInferenceGetVariable(self.ptr,c_char_p(name))
         
        if get:
            vtype=self.getVariableType(name)
            return cast(get,POINTER(vtype.litype)).contents.value
            
        else:
            raise KeyError
            
    def getVariableType(self,name):
        typei=LALInferenceGetVariableType(self.ptr,c_char_p(name))
        enum_name=LALInferenceVariableType_list[typei.value]
        return eval("LALINFERENCE_"+enum_name+"_t")
    
    def getVariableTypeByIndex(self,index):
        typei=LALInferenceGetVariableTypeByIndex(self.ptr,c_int(index))
        enum_name=LALInferenceVariableType_list[typei.value]
        return eval("LALINFERENCE_"+enum_name+"_t")
        
    def getVariableDimension(self):
        return int(LALInferenceGetVariableDimension(self.ptr).value)
    
    def getVariableDimensionNonFixed(self):
        return int(LALInferenceGetVariableDimensionNonFixed(self.ptr).value)
    
    def setVariable(self,name,value):
        valuep=pointer(value)
        value_voidp=cast(valuep,c_void_p)
        return LALInferenceSetVariable(self.ptr,c_char_p(name),value_voidp)
    
    def printVariables(self):
        return LALInferencePrintVariables(self.ptr)
    
    def checkVariable(self,name):
        return bool(LALInferenceCheckVariable(pointer(self),c_char_p(name)))
    
    def checkVariableNonFixed(self,name):
        return bool(LALInferenceCheckVariableNonFixed(pointer(self),c_char_p(name)))
    
    def __copy__(self):
        return self.__deepcopy__()
        
    def __deepcopy__(self):
        return self.copyVariables()
    
    def copyVariables(self):
        new=LALInferenceVariables()
        LALInferenceCopyVariables(self.ptr,byref(new))
        return new
    
    def setVariable(self,name,value):
        typeo=self.getVariableType(name)
        a=cast(pointer(typeo.litype(value)),c_void_p)
        LALInferenceSetVariable(self.ptr,c_char_p(name),a)
        
    def removeVariable(self,name):
        LALInferenceRemoveVariable(self.ptr,c_char_p(name))
        
    def __eq__(self,other):
        return bool(LALInferenceCompareVariables(self.ptr,byref(other)))
        
    def __neq__(self,other):
        return not self.__eq__(other)
        
    def __str__(self):
        return str(self.printVariables())
        

LALInferenceVariableItem._fields_ = [
    ("name",c_char_p),
    ("value",c_void_p),
    ("type",LALInferenceVariableType),
    ("vary",LALInferenceParamVaryType),
    ("next",POINTER(LALInferenceVariableItem))
]

#LALInferenceVariables
LALInferenceVariables_func_table=[
    ["LALInferenceCheckVariable",c_int,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCheckVariableNonFixed",c_int,[POINTER(LALInferenceVariables),c_char_p]],
]

LALInferenceVariables_xlal_func_table=[
    ["LALInferenceAddVariable",None,[POINTER(LALInferenceVariables),c_char_p,c_void_p,LALInferenceVariableType,LALInferenceParamVaryType]],
    ["LALInferenceGetVariable",c_void_p,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetVariableDimension",INT4,[POINTER(LALInferenceVariables)]],
    ["LALInferenceGetVariableDimensionNonFixed",INT4,[POINTER(LALInferenceVariables)]],
    ["LALInferenceGetVariableTypeByIndex",LALInferenceVariableType,[POINTER(LALInferenceVariables),c_int]],
    ["LALInferenceGetVariableType",LALInferenceVariableType,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetVariableVaryType",LALInferenceParamVaryType,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetVariableName",c_char_p,[POINTER(LALInferenceVariables),c_int]],
    ["LALInferenceSetVariable",None,[POINTER(LALInferenceVariables),c_char_p,c_void_p]],
    ["LALInferenceRemoveVariable",None,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCheckVariable",c_int,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCheckVariableNonFixed",c_int,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCopyVariables",None,[POINTER(LALInferenceVariables),POINTER(LALInferenceVariables)]],
    ["LALInferencePrintVariables",None,[POINTER(LALInferenceVariables)]],
    ["LALInferenceCompareVariables",c_int,[POINTER(LALInferenceVariables),POINTER(LALInferenceVariables)]]
]

#LALInferencePrior
class LALInferencePrior(object):
    def __init__(self):
        self._as_parameter_=LALInferenceVariables()
        self._ptr=pointer(self._as_parameter_)
        
    def addMinMaxPrior(self,name,min_value,max_value,li_variable_type):
        
        min_value_p=pointer(REAL8(min_value))
        max_value_p=pointer(REAL8(max_value))
        
        LALInferenceAddMinMaxPrior(self._ptr,c_char_p(name),min_value_p,max_value_p,li_variable_type)
        
    def removeMinMaxPrior(self,name):
        
        LALInferenceRemoveMinMaxPrior(self._ptr,c_char_p(name))
    
    def checkMinMaxPrior(self,name):
        
        return bool(LALInferenceCheckMinMaxPrior(self._ptr,c_char_p(name)))
    
    def getMinMaxPrior(self,name):
        min_p=pointer(REAL8())
        max_p=pointer(REAL8())
    
        LALInferenceGetMinMaxPrior(self._ptr,c_char_p(name),min_p,max_p)
        
        return min_p.contents.value,max_p.contents.value
    
    def addGaussianPrior(self,name,mu_value,sigma_value,li_variable_type):
        
        mu_value_p=pointer(REAL8(mu_value))
        sigma_value_p=pointer(REAL8(sigma_value))
        
        LALInferenceAddGaussianPrior(self._ptr,c_char_p(name),mu_value_p,sigma_value_p,li_variable_type)
        
    def removeGaussianPrior(self,name):
        
        LALInferenceRemoveGaussianPrior(self._ptr,c_char_p(name))
    
    def checkGaussianPrior(self,name):
        
        return bool(LALInferenceCheckGaussianPrior(self._ptr,c_char_p(name)))
    
    def getGaussianPrior(self,name):
        mu_p=pointer(REAL8())
        sigma_p=pointer(REAL8())
    
        LALInferenceGetGaussianPrior(self._ptr,c_char_p(name),mu_p,sigma_p)
        
        return mu_p.contents.value,sigma_p.contents.value
    
    def addCorrelatedPrior(self,name,cor_matrix,idx):
    
        cor_matrix_pp=pointer(pointer(cor_matrix))
        idx_p=pointer(UINT(idx))
        
        LALInferenceAddCorrelatedPrior(self._ptr,c_char_p(name),cor_matrix_pp,idx_p)
        
    def removeCorrelatedPrior(self,name):

        LALInferenceRemoveCorrelatedPrior(self._ptr,c_char_p(name))
        
    def checkCorrelatedPrior(self,name):
        
        return bool(LALInferenceCheckCorrelatedPrior(self._ptr,c_char_p(name)))
        
    def getCorrelatedPrior(self,name):
        pass
        

LALInferencePrior_func_table=[
    ["LALInferenceCheckMinMaxPrior",c_int,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCheckGaussianPrior",c_int,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceCheckCorrelatedPrior",c_int,[POINTER(LALInferenceVariables),c_char_p]],
]

LALInferencePrior_xlal_func_table=[
    ["LALInferenceAddMinMaxPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(REAL8),POINTER(REAL8),LALInferenceVariableType]],
    ["LALInferenceRemoveMinMaxPrior",None,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetMinMaxPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(REAL8),POINTER(REAL8)]],
    ["LALInferenceAddGaussianPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(REAL8),POINTER(REAL8),LALInferenceVariableType]],
    ["LALInferenceRemoveGaussianPrior",None,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetGaussianPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(REAL8),POINTER(REAL8)]],
    ["LALInferenceAddCorrelatedPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(POINTER(gsl_matrix)),POINTER(UINT4)]],
    ["LALInferenceRemoveCorrelatedPrior",None,[POINTER(LALInferenceVariables),c_char_p]],
    ["LALInferenceGetCorrelatedPrior",None,[POINTER(LALInferenceVariables),c_char_p,POINTER(POINTER(gsl_matrix)),POINTER(UINT4)]],
]



#LALInferenceIFOData
class LALInferenceIFOData(Structure): pass

#Function types
LALInferenceAlgorithm=CFUNCTYPE(None,POINTER(LALInferenceRunState))
LALInferenceEvolveOneStepFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState))
LALInferencePriorFunction=CFUNCTYPE(REAL8,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))
LALInferenceTemplateFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState))
LALInferenceLikelihoodFunction=CFUNCTYPE(REAL8,POINTER(LALInferenceVariables),POINTER(LALInferenceIFOData),POINTER(LALInferenceTemplateFunction))
LALInferenceProposalFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))
LALInferenceLogFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))

class BarycenterInput(Structure): pass
class EphemerisData(Structure): pass

#LALInferenceIFOData
LALInferenceIFOData._fields_ = [
    ('name',c_char_p),
    ('timeData',POINTER(REAL8TimeSeries)),
    ('timeModelhPlus',POINTER(REAL8TimeSeries)),
    ('timeModelhCross',POINTER(REAL8TimeSeries)),
    ('whiteTimeData',POINTER(REAL8TimeSeries)), 
    ('windowedTimeData',POINTER(REAL8TimeSeries)),
    ('nullloglikelihood',REAL8),
    ('loglikelihood',REAL8),
    ('acceptedloglikelihood',REAL8),
    ('fPlus',REAL8),
    ('fCross',REAL8),
    ('timeshift',REAL8),
    ('freqData',POINTER(COMPLEX16FrequencySeries)),
    ('freqModelhCross',POINTER(COMPLEX16FrequencySeries)),
    ('freqModelhPlus',POINTER(COMPLEX16FrequencySeries)),
    ('whiteFreqData',POINTER(COMPLEX16FrequencySeries)),
    ('compTimeData',POINTER(COMPLEX16TimeSeries)),
    ('compModelData',POINTER(COMPLEX16TimeSeries)),
    ('dataTimes',POINTER(LIGOTimeGPSVector)),
    ('modelParams',POINTER(LALInferenceVariables)),
    ('dataParams',POINTER(LALInferenceVariables)),
    ('modelDomain',LALInferenceDomain),
    ('oneSidedNoisePowerSpectrum',POINTER(REAL8FrequencySeries)),
    ('timeDomainNoiseWeights',POINTER(REAL8TimeSeries)),
    ('window',POINTER(REAL8Window)),
    ('timeToFreqFFTPlan',POINTER(REAL8FFTPlan)),
    ('fLow',REAL8),
    ('fHigh',REAL8),
    ('detector',POINTER(LALDetector)),
    ('bary',POINTER(BarycenterInput)),
    ('ephem',POINTER(EphemerisData)),
    ('epoch',LIGOTimeGPS),
    ('SNR',REAL8),
    ('STDOF',REAL8),
    ('next',POINTER(LALInferenceIFOData))
]

#ProcessParamsTable
class ProcessParamsTable(Structure): pass

#gsl_rng
class gsl_rng(Structure): pass

#LALInferenceRunState
LALInferenceRunState._fields_ = [
  ('commandLine',POINTER(ProcessParamsTable)),
  ('algorithm',POINTER(LALInferenceAlgorithm)),
  ('evolve',POINTER(LALInferenceEvolveOneStepFunction)),
  ('prior',POINTER(LALInferencePriorFunction)),
  ('likelihood',POINTER(LALInferenceLikelihoodFunction)),
  ('proposal',POINTER(LALInferenceProposalFunction)),
  ('template',POINTER(LALInferenceTemplateFunction)),
  ('logsample',POINTER(LALInferenceLogFunction)),
  ('data',POINTER(LALInferenceIFOData)),
  ('currentParams',POINTER(LALInferenceVariables)),
  ('priorArgs',POINTER(LALInferenceVariables)),
  ('proposalArgs',POINTER(LALInferenceVariables)),
  ('algorithmParams',POINTER(LALInferenceVariables)),
  ('livePoints',POINTER(POINTER(LALInferenceVariables))),
  ('differentialPoints',POINTER(POINTER(LALInferenceVariables))),
  ('differentialPointsLength',c_size_t),
  ('differentialPointsSize',c_size_t),
  ('currentLikelihood',REAL8),
  ('currentPrior',REAL8),
  ('GSLrandom',gsl_rng)
]

def __create_lalinference_functions(tables):
    for table in tables:
        for funcname,restype,argtypes in table:
            
            globals()[funcname]=_set_types(pylal.ctypes.liblalinference,funcname,restype,argtypes) 


def __create_lalinference_xlal_functions(tables):
    for table in tables:
        for funcname,restype,argtypes in table:
            globals()[funcname]=_set_xlal_types(pylal.ctypes.liblalinference,funcname,restype,argtypes) 

__create_lalinference_functions([LALInferenceVariables_func_table])
__create_lalinference_xlal_functions([LALInferenceVariables_xlal_func_table])

__create_lalinference_functions([LALInferencePrior_func_table])
__create_lalinference_xlal_functions([LALInferencePrior_xlal_func_table])
