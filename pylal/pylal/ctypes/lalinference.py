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

from pylal.datatypes import *

VARNAME_MAX=1024
   
#LALInferenceParamVaryType

globals().update(make_li_enum_typedef(
    "LALInferenceParamVaryType",
    [
        "LALINFERENCE_PARAM_LINEAR",
        "LALINFERENCE_PARAM_CIRCULAR",
        "LALINFERENCE_PARAM_FIXED",
        "LALINFERENCE_PARAM_OUTPUT"
    ]
))

#LALInferenceVariableType

globals().update(make_li_enum_typedef(
    "LALInferenceVariableType",
    [
        "LALINFERENCE_INT4_t",
        "LALINFERENCE_INT8_t",
        "LALINFERENCE_UINT4_t",
        "LALINFERENCE_REAL4_t", 
        "LALINFERENCE_REAL8_t", 
        "LALINFERENCE_COMPLEX8_t", 
        "LALINFERENCE_COMPLEX16_t", 
        "LALINFERENCE_gslMatrix_t",
        "LALINFERENCE_REAL8Vector_t",
        "LALINFERENCE_UINT4Vector_t",
        "LALINFERENCE_string_t",
        "LALINFERENCE_void_ptr_t"
    ]
))

#LALInferenceVariableType

globals().update(make_li_enum_typedef(
    "LALInferenceDomain",
    [
        "LALINFERENCE_DOMAIN_TIME", 
        "LALINFERENCE_DOMAIN_FREQUENCY"
    ]
))

#LALInferenceDomain

globals().update(make_li_enum_typedef(
    "LALInferenceDomain",
    [
        "LALINFERENCE_DOMAIN_TIME", 
        "LALINFERENCE_DOMAIN_FREQUENCY"
    ]
))

#LALInferenceApplyTaper

globals().update(make_li_enum_typedef(
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
        valuep=pointer(value)
        value_voidp=cast(valuep,c_void_p)
        LALInferenceAddVariable(byref(self),c_char_p(name),value_voidp,li_variable_type,li_param_vary_type)
    
    def removeVariable(self,name):
        LALInferenceRemoveVariable(self.ptr,c_char_p(name))
    
    def getVariable(self,name):
        print LALInferenceGetVariable(self.ptr,c_char_p(name))
        
    def getVariableDimension(self):
        return LALInferenceGetVariableDimension(self.ptr)
        
    def setVariable(self,name,value):
        valuep=pointer(value)
        value_voidp=cast(valuep,c_void_p)
        return LALInferenceSetVariable(self.ptr,c_char_p(name),value_voidp)
    
    def printVariables(self):
        return LALInferencePrintVariables(self.ptr)
        
    def printSample(self):
        pass        
    

LALInferenceVariableItem._fields_ = [
    ("name",c_char_p),
    ("value",c_void_p),
    ("type",LALInferenceVariableType),
    ("vary",LALInferenceParamVaryType),
    ("next",POINTER(LALInferenceVariableItem))
]

#LALInferenceVariables
LALInferenceVariables_table=[
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
    ["LALInferenceCopyVariables",None,[POINTER(LALInferenceVariables),POINTER(LALInferenceVariables)]],
    ["LALInferencePrintVariables",None,[POINTER(LALInferenceVariables)]]
]

def __create_lalinference_functions(tables):
    for table in tables:
        for funcname,restype,argtypes:
            globals[funcname]=_set_types(pylal.ctypes.liblalinference,funcname,restype,argtypes) 

##LALInferenceIFOData
#class LALInferenceIFOData(Structure): pass

##Function types
#LALInferenceAlgorithm=CFUNCTYPE(None,POINTER(LALInferenceRunState))
#LALInferenceEvolveOneStepFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState))
#LALInferencePriorFunction=CFUNCTYPE(REAL8,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))
#LALInferenceTemplateFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState))
#LALInferenceLikelihoodFunction=CFUNCTYPE(REAL8,POINTER(LALInferenceVariables),POINTER(LALInferenceIFOData),POINTER(LALInferenceTemplateFunction))
#LALInferenceProposalFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))
#LALInferenceLogFunction=CFUNCTYPE(None,POINTER(LALInferenceRunState),POINTER(LALInferenceVariables))

#class BarycenterInput(Structure): pass
#class EphemerisData(Structure): pass

##LALInferenceIFOData
#LALInferenceIFOData._fields_ = [
    #('name',c_char_p),
    #('timeData',POINTER(REAL8TimeSeries)),
    #('timeModelhPlus',POINTER(REAL8TimeSeries)),
    #('timeModelhCross',POINTER(REAL8TimeSeries)),
    #('whiteTimeData',POINTER(REAL8TimeSeries)), 
    #('windowedTimeData',POINTER(REAL8TimeSeries)),
    #('nullloglikelihood',REAL8),
    #('loglikelihood',REAL8),
    #('acceptedloglikelihood',REAL8),
    #('fPlus',REAL8),
    #('fCross',REAL8),
    #('timeshift',REAL8),
    #('freqData',POINTER(COMPLEX16FrequencySeries)),
    #('freqModelhCross',POINTER(COMPLEX16FrequencySeries)),
    #('freqModelhPlus',POINTER(COMPLEX16FrequencySeries)),
    #('whiteFreqData',POINTER(COMPLEX16FrequencySeries)),
    #('compTimeData',POINTER(COMPLEX16TimeSeries)),
    #('compModelData',POINTER(COMPLEX16TimeSeries)),
    #('dataTimes',POINTER(LIGOTimeGPSVector)),
    #('modelParams',POINTER(LALInferenceVariables)),
    #('dataParams',POINTER(LALInferenceVariables)),
    #('modelDomain',LALInferenceDomain),
    #('oneSidedNoisePowerSpectrum',POINTER(REAL8FrequencySeries)),
    #('timeDomainNoiseWeights',POINTER(REAL8TimeSeries)),
    #('window',POINTER(REAL8Window)),
    #('timeToFreqFFTPlan',POINTER(REAL8FFTPlan)),
    #('fLow',REAL8),
    #('fHigh',REAL8),
    #('detector',POINTER(LALDetector)),
    #('bary',POINTER(BarycenterInput)),
    #('ephem',POINTER(EphemerisData)),
    #('epoch',LIGOTimeGPS),
    #('SNR',REAL8),
    #('STDOF',REAL8),
    #('next',POINTER(LALInferenceIFOData))
#]

##ProcessParamsTable
#class ProcessParamsTable(Structure): pass


##LALInferenceRunState
#LALInferenceRunState._fields_ = [
  #('commandLine',POINTER(ProcessParamsTable)),
  #('algorithm',POINTER(LALInferenceAlgorithm)),
  #('evolve',POINTER(LALInferenceEvolveOneStepFunction)),
  #('prior',POINTER(LALInferencePriorFunction)),
  #('likelihood',POINTER(LALInferenceLikelihoodFunction)),
  #('proposal',POINTER(LALInferenceProposalFunction)),
  #('template',POINTER(LALInferenceTemplateFunction)),
  #('logsample',POINTER(LALInferenceLogFunction)),
  #('data',POINTER(LALInferenceIFOData)),
  #('currentParams',POINTER(LALInferenceVariables)),
  #('priorArgs',POINTER(LALInferenceVariables)),
  #('proposalArgs',POINTER(LALInferenceVariables)),
  #('algorithmParams',POINTER(LALInferenceVariables)),
  #('livePoints',POINTER(POINTER(LALInferenceVariables))),
  #('differentialPoints',POINTER(POINTER(LALInferenceVariables))),
  #('differentialPointsLength',c_size_t),
  #('differentialPointsSize',c_size_t),
  #('currentLikelihood',REAL8),
  #('currentPrior',REAL8),
  #('GSLrandom',gsl_rng)
#]
