import numpy as np
from ctypes import c_double, POINTER, cast, _Pointer, Array
from typing import Any

def dp_array_to_pointer(array: np.ndarray, copy = False) -> _Pointer:
    """
    Convert a numpy array to a pointer to a double array.
    """
    # Check if the array is contiguous
    if not array.flags['F_CONTIGUOUS'] or array.dtype != np.float64:
        array = np.asfortranarray(array, dtype= np.float64)
    # If copy is True, make a copy of the array
    if copy:
        array = np.array(array, copy= True, order= 'F').reshape(array.shape)
    return array.ctypes.data_as(POINTER(c_double))

def pointer_to_dp_array(
        pointer: _Pointer, 
        shape: tuple[int, ...],
        copy: bool = False
    ) -> np.ndarray:
    """
    Convert a pointer to a double array to a numpy array.
    """
    data_dtype: type[c_double] | type[Array[Any]] = (c_double)
    for dim in shape:
        data_dtype = data_dtype * dim
    data_ptr_type = data_dtype

    data_ptr = cast(pointer, POINTER(data_ptr_type)).contents 
    data_ = np.frombuffer(
        data_ptr, 
        dtype =np.float64,
        count = np.prod(shape),
    ).reshape(shape, order='F')
    return data_

def _ctype_ndarray(element_type, shape):
    """ Create an ndarray of the given element type and shape """
    for dim in shape[::-1]:
        element_type = dim * element_type
        # prevent the type name include np.ctypeslib
        element_type.__module__ = None
    return element_type