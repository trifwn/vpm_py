import numpy as np
from ctypes import c_double, POINTER, cast, _Pointer, Array


def dp_array_to_pointer(array: np.ndarray, copy = False) -> _Pointer:
    """
    Convert a numpy array to a pointer to a double array.
    """
    # Check if the array is contiguous
    if not array.flags['C_CONTIGUOUS']:
        array = np.ascontiguousarray(array, dtype= np.float64)   
    
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
    size_array = (c_double)
    for dim in shape:
        size_array *= dim
    
    pointer = cast(pointer, POINTER(size_array))
    data = np.ctypeslib.as_array(pointer, shape=shape)
    # data = np.array(data, copy= copy, order= 'F').reshape(shape)
    return data
