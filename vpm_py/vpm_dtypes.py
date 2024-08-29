import numpy as np
from ctypes import c_double, POINTER, cast, _Pointer, Array
np.set_printoptions(precision= 3, suppress= True)

def dp_array_to_pointer(array: np.ndarray, copy = False) -> _Pointer:
    """
    Convert a numpy array to a pointer to a double array.
    """
    # Check if the array is contiguous
    if not array.flags['F_CONTIGUOUS']:
        array = np.asfortranarray(array)
    # Check if the array is double precision
    if array.dtype != np.float64:
        array = array.astype(np.float64)
    # If copy is True, make a copy of the array
    if copy:
        array = np.array(array, copy= True, order= 'F').reshape(array.shape)
    # if array.ndim == 2:
    #     print(f"Got array with shape {array.shape}")
    #     print(f"array[0,:]")
    #     print(array[0,:])
    #     print(f"array[1,:]")
    #     print(array[1,:])
    #     print(f"array[2,:]")
    #     print(array[2,:])
    return array.ctypes.data_as(POINTER(c_double))

def pointer_to_dp_array(
        pointer: _Pointer, 
        shape: tuple[int, ...],
        copy: bool = False
    ) -> np.ndarray:
    """
    Convert a pointer to a double array to a numpy array.
    """
    pointer = cast(pointer, POINTER(c_double))
    data = np.ctypeslib.as_array(pointer, shape= shape)
    data = data.reshape(shape, order= 'F')
    if copy:
        data = np.array(data, copy= True, order= 'F')
    return data
