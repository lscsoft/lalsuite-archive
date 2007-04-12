#!/usr/bin/env python
"""
MetaArray - a subclass of ndarray that holds metadata and preserves it across
            array operations.
Metadata - a class for metadata stored in MetaArray
MetaArrayList - a subclass of list that ollows element-wise MetaArray
                operations

Spectrum - a subclass of MetaArray as demonstration
SpectrumMetadata - a subclass of MetaData that strictly checks for metadata
                   compatibility between two Spectra
SpectrumList - subclass of MetaArrayList that has some nice features specific
               to SpectrumMetadata

$Id$
"""
from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"[11:-2]

import os
import sys
import copy

import numpy
from glue import segments

##############################################################################
# Method wrappers
##############################################################################

class _arraymethod(object):
    """
    Used to attach methods to MetaArrays.  Wrap ndarray methods that return
    a single array.  Merge metadata of all input Spectra.
    
    __init__ gets called when we define the MetaArray (sub-)class and attach
        methods.
    __get__ gets called on the object at method call time.
    __call__ is called immediately after __get__.
    """
    def __init__(self, methname):
        self._name = methname
        self.obj = None
        self.__doc__ = getattr(numpy.ndarray, self._name, None).__doc__
    
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    
    def __call__(self, *args, **params):
        data = self.obj.A
        cls = type(self.obj)
        result = getattr(data, self._name)(*args, **params).view(cls)
        result.metadata = self.obj.metadata.copy()
        for arg in args:
            result.metadata |= getattr(arg, 'metadata', None)
        return result

class _elementwise_method(object):
    """
    Used to attach methods to MetaArrayLists.  Wrap ndarray methods that
    operate upon arrays and apply them to each element of the list.
    
    __init__ gets called when we define the MetaArrayList class and attach
        methods.
    __get__ gets called on the list object at method call time.
    __call__ is called immediately after __get__.
    """
    def __init__(self, methname):
        self._name = methname
        self.list = None
        self.listtype = None
        self.itemtype = None
        # I don't think I can get the docstring, unfortunately
        # self.__doc__ = getattr(itemtype, self._name, None).__doc__

    def __get__(self, list, listtype=None):
        self.list = list
        self.listtype = type(list)
        self.itemtype = list._itemtype
        return self

    def __call__(self, *others, **params):
        result = self.listtype([])
        cls = self.itemtype
        # if we have two alike lists, just take pairwise action
        if len(others) == 1 and isinstance(others[0], self.listtype):
            result.extend([getattr(a, self._name)(b, **params).view(cls) \
                for a, b in zip(self.list, others[0])])
        # otherwise, just try passing the args along
        else:
            result.extend([getattr(a, self._name)(*others, **params).view(cls) for a in self.list])
        return result

##############################################################################
# Generic classes
##############################################################################

class Metadata(object):
    """
    Abstract class to hold metadata
    """
    # specify all allowed metadata here; type is required
    typemap = {}
    __slots__ = typemap.keys() + ['typemap']
    
    def __new__(cls, metadata=None):
        if isinstance(metadata, Metadata):
            return metadata
        return super(Metadata, cls).__new__(cls)
    
    def __init__(self, metadata):
        # Recycle if possible
        if isinstance(metadata, Metadata):
            return
        # user forgot to provide metadata
        elif metadata is None:
            return None
        # create Metadata from dictionary
        elif isinstance(metadata, dict):
            for key, val in metadata.items():
                try:
                    setattr(self, key, self.typemap[key](val))
                except KeyError, e:
                    raise KeyError, \
                        "This key is not in the typemap of %s: %s" \
                        % (type(self), str(e))
            
            # all fields must be filled
            for slot in self.__slots__:
                if not hasattr(self, slot):
                    raise ValueError, \
                        "Not enough metadata supplied; missing %s" % slot
        else:
            raise NotImplemented
    
    def __str__(self):
        return str(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))
    
    def __repr__(self):
        return repr(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))
    
    def __ior__(self, other):
        """
        Merge metadata; this must be subclassed.
        """
        if self is None:
            return other
        if other is None:
            return self
        raise NotImplemented
    
    def __or__(self, other):
        return self.copy().__ior__(other)
    
    def copy(self):
        return type(self)(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))

class MetaArray(numpy.ndarray):
    """
    An array containing a data and metadata.  Intended to be
    subclassed.
    
    On b = MetaArray(a) where a is a MetaArray, metadata is copied, but data
    is not.
    """
    __array_priority__ = 10.1  # ufuncs mixing array types return MetaArray
    _metadata_type = Metadata
    
    def __new__(subtype, data=None, metadata=None, dtype=None, copy=False,
        subok=True):
        # Case 1: data is an MetaArray, it has metadata already.
        if isinstance(data, MetaArray):
            if dtype is None:
                dtype = data.dtype
            else:
                dtype = numpy.dtype(dtype)
            if not copy and dtype==data.dtype and metadata is None:
                return data
            elif metadata is not None:
                # share memory, but not metadata
                new = numpy.array(data, dtype=dtype, copy=copy, subok=True)
                new.metadata = subtype._metadata_type(metadata)
                return new
            else:
                # copy, update metadata, then return
                new = data.astype(dtype)
                new.metadata = metadata
                new._baseclass = _baseclass
                return new
        
        # All other cases, create a new array and attach metadata
        # Unless you specified otherwise, we'll reuse memory from existing
        # arrays.
        new = numpy.array(data, dtype=dtype, copy=copy, subok=True)
        _baseclass = type(new)
        new = new.view(subtype)
        new.metadata = subtype._metadata_type(metadata)
        new._baseclass = _baseclass
        return new
    
    def __array_finalize__(self, obj):
        """
        Called anytime a MetaArray is returned; make sure that metadata
        is set to something.
        """
        self.metadata = getattr(obj, "metadata", None)
        self._baseclass = getattr(obj, '_baseclass', type(obj))
    
    def __array_wrap__(self, obj):
        """
        Called anytime a ufunc operates on a MetaArray and another object.
        The result of the ufunc is obj.  The MetaArray operand is self.
        """
        result = obj.view(type(self))
        result.metadata = self.metadata.copy()
        return result
    
    def __repr__(self):
        return "%s, %s)" % (numpy.ndarray.__repr__(self)[:-1], repr(self.metadata))
    
    def __str__(self):
        return "%s %s" % (numpy.ndarray.__str__(self), str(self.metadata))
    
    # methods that return an array, wrapped to return a MetaArray
    __abs__ = _arraymethod('__abs__')
    __add__ = _arraymethod('__add__')
    __and__ = _arraymethod('__and__')
    __copy__ = _arraymethod('__copy__')
    __deepcopy__ = _arraymethod('__deepcopy__')
    __div__ = _arraymethod('__div__')
    __divmod__ = _arraymethod('__divmod__')
    __floordiv__ = _arraymethod('__floordiv__')
    __hex__ = _arraymethod('__hex__')
    __iadd__ = _arraymethod('__iadd__')
    __iand__ = _arraymethod('__iand__')
    __idiv__ = _arraymethod('__idiv__')
    __ifloordiv__ = _arraymethod('__ifloordiv__')
    __ilshift__ = _arraymethod('__ilshift__')
    __imod__ = _arraymethod('__imod__')
    __imul__ = _arraymethod('__imul__')
    __invert__ = _arraymethod('__invert__')
    __ior__ = _arraymethod('__ior__')
    __ipow__ = _arraymethod('__ipow__')
    __irshift__ = _arraymethod('__irshift__')
    __isub__ = _arraymethod('__isub__')
    __itruediv__ = _arraymethod('__itruediv__')
    __ixor__ = _arraymethod('__ixor__')
    __lshift__ = _arraymethod('__lshift__')
    __mul__ = _arraymethod('__mul__')
    __rmod__ = _arraymethod('__rmod__')
    __rmul__ = _arraymethod('__rmul__')
    __ror__ = _arraymethod('__ror__')
    __rpow__ = _arraymethod('__rpow__')
    __rrshift__ = _arraymethod('__rrshift__')
    __rshift__ = _arraymethod('__rshift__')
    __rsub__ = _arraymethod('__rsub__')
    __rtruediv__ = _arraymethod('__rtruediv__')
    __rxor__ = _arraymethod('__rxor__')
    __sub__ = _arraymethod('__sub__')
    __truediv__ = _arraymethod('__truediv__')
    __xor__ = _arraymethod('__xor__')
    astype = _arraymethod('astype')
    byteswap = _arraymethod('byteswap')
    choose = _arraymethod('choose')
    clip = _arraymethod('clip')
    compress = _arraymethod('compress')
    conj = _arraymethod('conj')
    conjugate = _arraymethod('conjugate')
    copy = _arraymethod('copy')
    cumprod = _arraymethod('cumprod')
    cumsum = _arraymethod('cumsum')
    diagonal = _arraymethod('diagonal')
    fill = _arraymethod('fill')
    flat = _arraymethod('flat')
    flatten = _arraymethod('flatten')
    repeat = _arraymethod('repeat')
    squeeze = _arraymethod('squeeze')
    transpose = _arraymethod('transpose')
    
    T = property(fget = lambda self: self.transpose())
    H = property(fget = lambda self: self.T.conj()) # Hermitian transpose
    
    def _get_data(self):
        return self.view(self._baseclass)
    A = property(fget=_get_data) # get at the underlying Array data

class MetaArrayList(list):
    _itemtype = MetaArray
    
    # these methods will act upon each element of this list
    __abs__ = _elementwise_method('__abs__')
    __add__ = _elementwise_method('__add__')
    __and__ = _elementwise_method('__and__')
    __copy__ = _elementwise_method('__copy__')
    __deepcopy__ = _elementwise_method('__deepcopy__')
    __div__ = _elementwise_method('__div__')
    __divmod__ = _elementwise_method('__divmod__')
    __floordiv__ = _elementwise_method('__floordiv__')
    __hex__ = _elementwise_method('__hex__')
    __iadd__ = _elementwise_method('__iadd__')
    __iand__ = _elementwise_method('__iand__')
    __idiv__ = _elementwise_method('__idiv__')
    __ifloordiv__ = _elementwise_method('__ifloordiv__')
    __ilshift__ = _elementwise_method('__ilshift__')
    __imod__ = _elementwise_method('__imod__')
    __imul__ = _elementwise_method('__imul__')
    __invert__ = _elementwise_method('__invert__')
    __ior__ = _elementwise_method('__ior__')
    __ipow__ = _elementwise_method('__ipow__')
    __irshift__ = _elementwise_method('__irshift__')
    __isub__ = _elementwise_method('__isub__')
    __itruediv__ = _elementwise_method('__itruediv__')
    __ixor__ = _elementwise_method('__ixor__')
    __lshift__ = _elementwise_method('__lshift__')
    __mul__ = _elementwise_method('__mul__')
    __rmod__ = _elementwise_method('__rmod__')
    __rmul__ = _elementwise_method('__rmul__')
    __ror__ = _elementwise_method('__ror__')
    __rpow__ = _elementwise_method('__rpow__')
    __rrshift__ = _elementwise_method('__rrshift__')
    __rshift__ = _elementwise_method('__rshift__')
    __rsub__ = _elementwise_method('__rsub__')
    __rtruediv__ = _elementwise_method('__rtruediv__')
    __rxor__ = _elementwise_method('__rxor__')
    __sub__ = _elementwise_method('__sub__')
    __truediv__ = _elementwise_method('__truediv__')
    __xor__ = _elementwise_method('__xor__')
    astype = _elementwise_method('astype')
    byteswap = _elementwise_method('byteswap')
    choose = _elementwise_method('choose')
    clip = _elementwise_method('clip')
    compress = _elementwise_method('compress')
    conj = _elementwise_method('conj')
    conjugate = _elementwise_method('conjugate')
    copy = _elementwise_method('copy')
    cumprod = _elementwise_method('cumprod')
    cumsum = _elementwise_method('cumsum')
    diagonal = _elementwise_method('diagonal')
    fill = _elementwise_method('fill')
    flat = _elementwise_method('flat')
    flatten = _elementwise_method('flatten')
    repeat = _elementwise_method('repeat')
    squeeze = _elementwise_method('squeeze')
    sum = _elementwise_method('sum')
    transpose = _elementwise_method('transpose')
    
    T = property(fget = lambda self: self.transpose()) # Transpose
    H = property(fget = lambda self: self.T.conj()) # Hermitian transpose
    
    # special case a few attribute accessors
    def _get_real(self):
        return type(self)([x.real for x in self])
    real = property(fget=_get_real)
    def _get_imag(self):
        return type(self)([x.imag for x in self])
    imag = property(fget=_get_imag)
    
    def _get_data(self):
        """
        Return a regular list of arrays.
        """
        return [x.view(x._baseclass) for x in self]
    A = property(fget=_get_data)

##############################################################################
# Define Spectrum and associated structures
##############################################################################

class SpectrumMetadata(Metadata):
    """
    Hold the metadata associated with a spectrum, including frequency
    resolution, a segmentlist indicating what times were involved in taking
    the spectrum, a channel name, etc.
    """
    # specify all allowed metadata here; type is required
    typemap = {"df": float,
               "f_low": float,
               "f_high": float,
               "num_bins": int,
               "segments": segments.segmentlist,
               "comments": list}
    __slots__ = typemap.keys() + ['typemap']
    
    def __ior__(self, other):
        """
        Merge metadata.  No repeats.  Throw error on incompatible spectra.
        Let None act as identity.
        """
        if other is None:
            return self
        if self is None:
            return other
        
        # check that metadata are compatible for merging
        assert numpy.alltrue([getattr(self, attr) == getattr(other, attr) \
            for attr in self.__slots__ \
            if attr not in ('segments', 'comments', 'typemap')])
        
        # union segments
        self.segments |= other.segments
        
        # add only new comments
        self.comments.extend([comment for comment in other.comments \
            if comment not in self.comments])
        return self

class Spectrum(MetaArray):
    """ This is a MetaArray, but with the metadata typemap specified. """
    _metadata_type = SpectrumMetadata
    pass

class SpectrumList(MetaArrayList):
    """
    A list of Spectra with a few convenience functions defined for all lists
    of spectra.  Intended to be subclassed.
    """
    _itemtype = Spectrum
    
    def segments(self):
        """
        Return the (uncoalesced) list of segments represented by the Spectra
        in this SpectrumList.
        """
        segs = segments.segmentlist()
        for spectrum in self:
            segs.extend(spectrum)
        return segs
    
    def extent(self):
        return self.segments().extent()

class SpectrumDict(dict):
    """
    A dictionary allowing access to FFTs of different data streams and
    providing convenient batch operations.
    """
    pass