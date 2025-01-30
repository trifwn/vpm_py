import numpy as np
from ctypes import c_double, POINTER, cast, _Pointer

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
    p_arr_type = POINTER(_ctype_ndarray(c_double, shape))

    obj = cast(pointer, p_arr_type).contents
    data = np.asarray(obj, dtype= np.float64, order= 'F').reshape(shape)
    if copy:
        data = np.array(data, copy= True, order= 'F')
    return data

def _ctype_ndarray(element_type, shape):
    """ Create an ndarray of the given element type and shape """
    for dim in shape[::-1]:
        element_type = dim * element_type
        # prevent the type name include np.ctypeslib
        element_type.__module__ = None
    return element_type