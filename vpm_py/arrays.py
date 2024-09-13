import numpy as np
import ctypes
from ctypes import c_int,  byref, POINTER, c_double, cast, cdll, Array
from ctypes.util import find_library
from ctypes import CDLL
from tempfile import NamedTemporaryFile
import os
import glob
from shutil import copy2
from typing import Any

from vpm_py.vpm_io import print_red, print_green, print_IMPORTANT
from vpm_py.vpm_dtypes import dp_array_to_pointer, pointer_to_dp_array

here = os.path.abspath(os.path.dirname(__file__))
lib_locations = os.path.join(here, 'shared_libs')
lib_arrays_path = glob.glob(os.path.join(lib_locations, 'libvpm_*.so'))[0]
lib_arrays_ext = lib_arrays_path[lib_arrays_path.rfind('.'):]

class F_Array_Struct(ctypes.Structure):
    _fields_ = [
        ("ndims", c_int),              # Integer for the number of dimensions
        ("total_size", c_int),         # Integer for the total size of the data array
        ("shape_ptr", POINTER(c_int)),  # Pointer to the shape array (array of integers)
        ("data_ptr", POINTER(c_double)) # Pointer to the data array (array of doubles)
    ]

    def __init__(self, ndims=0, total_size=0, shape_ptr=None, data_ptr=None):
        super().__init__()
        self.ndims = c_int(ndims)
        self.total_size = c_int(total_size)
        self.shape_ptr = shape_ptr
        self.data_ptr = data_ptr

    def __str__(self) -> str:
        _str = f"ND_Array with {self.ndims} dimensions and total size {self.total_size}\n"
        # Print the memory address of the shape and data pointers
        # Check for null pointers
        if not self.shape_ptr:
            _str += f"Shape pointer is null\n"
        else:
            _str += f"Shape pointer address: {ctypes.addressof(self.shape_ptr.contents)}\n"
        if not self.data_ptr:
            _str += f"Data pointer is null\n"
        else:
            _str += f"Data pointer address: {ctypes.addressof(self.data_ptr.contents)}\n"
        _str += f"Shape pointer: {self.shape_ptr}\n"
        _str += f"Data pointer: {self.data_ptr}\n"
        return _str
    
    def is_null(self):
        """
        Should return True if the shape or data pointers are null
        """
        if not self.shape_ptr:
            return True
        if not self.data_ptr:
            return True
        return False

    @classmethod
    def null(cls, ndims, total_size):
        return cls(
            ndims= ndims,
            total_size= total_size,
            data_ptr= cast(0, POINTER(c_double)),
            shape_ptr= cast(0, POINTER(c_int))
        )

class F_Array(object):
    _lib_array: CDLL = None  
    _lib_array_path = None
    libc = None

    def __init__(
            self, 
            shape: tuple[int,...] | np.ndarray[Any, np.dtype[np.int32]], 
            data_container : np.ndarray | None =None, 
            name: str | None = None
        ):
        """
        An F_Array is a data container that can easily interoperate with fortran.

        Args:
            shape (tuple[int,...] | np.ndarray[Any, np.dtype[np.int32]]): The array shape
            data_container (np.ndarray | None, optional): The actual data. If provided will be used otherwise will be initialized to 0.
            name (str | None, optional): Debugging parameter to track arrays. Will be removed.
        """
        self.name = name
        if not F_Array._lib_array:
            F_Array._load_library()
        self.ndims = len(shape)
        self.total_size = int(np.prod(shape)) # Total size of the data array
        self.shape_container = np.array(shape, dtype=np.int32, order='F')

        if data_container is not None:
            self.data_container = None
            self.owns_data = False
        else:
            data_container = np.zeros(shape, dtype=np.float64, order='F')
            self.data_container = data_container
            self.owns_data = True
        data_dtype: type[c_double] | type[Array[Any]] = (c_double)
        for dim in shape:
            data_dtype = data_dtype * dim
        self.data_ptr_type = data_dtype
        self.shape_ptr = self.shape_container.ctypes.data_as(POINTER(c_int * self.ndims))
        self.data_ptr =  data_container.ctypes.data_as(POINTER(data_dtype))

        data_ptr_zero = cast(self.data_ptr, POINTER(c_double))
        shape_ptr_zero = cast(self.shape_ptr, POINTER(c_int))
        _ = self._lib_array.create_dtype(
            c_int(self.ndims), c_int(self.total_size), shape_ptr_zero, data_ptr_zero
        )

        # print_red(f"Created array {self.name}")
        # print_green(f"\tThe array owns the data: {self.owns_data}")
        # print_green(f"\tThe array has shape {self.shape}")
        # print_green(f"\tThe array data are in {data_container.__array_interface__["data"]} ")

    @classmethod
    def _load_library(cls):
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_arrays_ext)
        tmp.close()
        F_Array._lib_array = cdll.LoadLibrary(copy2(lib_arrays_path, tmp.name))
        F_Array._lib_array_path = tmp.name
        F_Array._setup_function_signatures()
        F_Array.libc = CDLL(find_library('c'))

    
    @classmethod
    def _setup_function_signatures(cls):
        cls._lib_array.create_dtype.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_double)]
        cls._lib_array.create_dtype.restype = F_Array_Struct
        cls._lib_array.get_data_ptr.argtypes = []
        cls._lib_array.get_data_ptr.restype = POINTER(c_double)
        cls._lib_array.get_shape_ptr.argtypes = []
        cls._lib_array.get_shape_ptr.restype = POINTER(c_int)
        cls._lib_array.add.argtypes = [POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(F_Array_Struct)]
        cls._lib_array.sub.argtypes = [POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(F_Array_Struct)]
        cls._lib_array.mul.argtypes = [POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(F_Array_Struct)]
        cls._lib_array.div.argtypes = [POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(F_Array_Struct)]
        cls._lib_array.print_array.argtypes = [POINTER(F_Array_Struct)]
        cls._lib_array.free_array.argtypes = [POINTER(F_Array_Struct)]

    @classmethod
    def from_ctype(
        self, 
        array_struct: F_Array_Struct, 
        ownership= False, 
        name=None
    ):
        # Get the shape from the shape pointer
        ndims = array_struct.ndims
        shape_ptr = array_struct.shape_ptr
        shape_ptr = cast(shape_ptr, POINTER(c_int * ndims))
        shape = np.frombuffer(
            shape_ptr.contents, 
            dtype=np.int32,
            count=ndims
        )
        data_ptr = array_struct.data_ptr
        data_ptr = cast(data_ptr, POINTER(c_double*array_struct.total_size))
        data = np.asarray(data_ptr.contents, dtype=np.float64).reshape(shape, order='F')
        
        arr = F_Array(shape, data_container=data, name=name)
        return arr
    
    @classmethod
    def from_ndarray(self, np_array):
        # If the array is not Fortran ordered, make it Fortran ordered
        if not np_array.flags['F_CONTIGUOUS']:
            np_array = np.asfortranarray(np_array, dtype=np.float64)
        arr = F_Array(np_array.shape, data_container=np_array)
        return arr
    
    def to_numpy(self):
        return self.data

    def copy(self):
        data = np.array(self.data, copy=True)
        return F_Array.from_ndarray(data)
    
    def transfer_data_ownership(self):
        """
        This function transfers the ownership of the data to the caller so that 
        the data are not deleted when the object is deleted. If the data are not
        owned by the object, the data view is returned to the caller.

        Returns:
            ndarray: The data array
        """
        if self.owns_data:
            # Ensure the data are not deleted and just returned to the caller
            self.owns_data = False
            ret = self.data_container
            self.data_container = None
        else:
            ret = self.data
        return ret

    def to_ctype(self):
        return F_Array_Struct(self.ndims, self.total_size, self.shape_ptr[0], self.data_ptr[0][0])

    def print_in_fortran(self):
        array_struct = self.to_ctype()
        self._lib_array.print_array(byref(array_struct))

    def _call_fortran_op(self, other: "F_Array", op: str):
        # Call Fortran operation
        result_array = F_Array(np.array(self.shape, copy=True, dtype=np.float64), name=f"{self.name}_{op}_{other.name}")

        if op == 'add':
            self._lib_array.add(byref(self.to_ctype()), byref(other.to_ctype()), byref(result_array.to_ctype()))
        elif op == 'sub':
            self._lib_array.sub(byref(self.to_ctype()), byref(other.to_ctype()), byref(result_array.to_ctype()))
        elif op == 'mul':
            self._lib_array.mul(byref(self.to_ctype()), byref(other.to_ctype()), byref(result_array.to_ctype()))
        elif op == 'div':
            self._lib_array.div(byref(self.to_ctype()), byref(other.to_ctype()), byref(result_array.to_ctype()))
        else:
            raise ValueError(f"Operation {op} not supported")
        return F_Array.from_ctype(result_array)

    def __check_elementwise_op(self, other):
        is_array = isinstance(other, F_Array)
        is_np_array = isinstance(other, np.ndarray)
        if not is_array and not is_np_array:
            raise ValueError("Can only add Array or numpy array")
        if len(other.shape) != self.ndims:
            raise ValueError(f"Arrays must have the same number of dimensions. Got arrays with shapes: {self.shape} and {other.shape}")
        if np.any(self.shape != other.shape):
            raise ValueError(f"Arrays must have the same shape. Got arrays with shapes: {self.shape} and {other.shape}")
        return is_array, is_np_array

    def __add__(self, other):
        other_is_array, other_is_np_array = self.__check_elementwise_op(other)
        if other_is_np_array:
            result = self.data + other
        if other_is_array:
            result = self._call_fortran_op(other, 'add')
        return result
        
    def __sub__(self, other):
        is_array, is_np_array = self.__check_elementwise_op(other)
        if is_np_array:
            result = self.data - other
        if is_array:
            result = self._call_fortran_op(other, 'sub')
        return result
    
    def __mul__(self, other):
        is_array, is_np_array = self.__check_elementwise_op(other)
        if is_np_array:
            result = self.data * other
        if is_array:
            result = self._call_fortran_op(other, 'mul')
        return result
    
    def __truediv__(self, other):
        is_array, is_np_array = self.__check_elementwise_op(other)
        if is_np_array:
            result = self.data / other
        if is_array:
            result = self._call_fortran_op(other, 'div')
        return result
    
    def __matmul__(self, other):
        raise NotImplementedError("Matrix multiplication not implemented yet")

    def __del__(self):
        try:
            # print_red(f"Deleting array {self.name}")
            # print_green(f"\tThe array owns the data: {self.owns_data}")
            # print_green(f"\tThe array has shape {self.shape}")
            # print_green(f"\tThe array data are in {self.data.__array_interface__["data"]} ")
            # Ensure to call free_array correctly
            # Free the memory allocated by the Fortran code
            # self.libc.free(self.data_ptr)
            del self.data_container
            del self.shape_container
            del self.shape_ptr
            del self.data_ptr
            del self.ndims
            del self.total_size
            # if self._lib_array_path:
                # os.remove(self._lib_array_path)
            del self
        except Exception as e:
            print(f"Error during deletion: {e}")

    def __getitem__(self, index):
        "Return the data array with the given index"
        return np.asfortranarray(self.data)[index]

    def __setitem__(self, index, value):
        is_array = isinstance(value, F_Array)
        is_np_array = isinstance(value, np.ndarray)
        is_scalar = np.isscalar(value)

        if is_scalar:
            self.data[index] = value
        elif is_array:
            self.data[index] = value.data
        elif is_np_array:
            self.data[index] = value
    
    @property
    def shape(self) -> np.ndarray[Any, np.dtype[np.int32]]:
        shape = np.frombuffer(
            cast(self.shape_ptr, POINTER(c_int * self.ndims)).contents, 
            dtype=np.int32,
            count=self.ndims
        )
        return shape
    
    @property
    def data(self)-> np.ndarray[Any, np.dtype[np.floating]]:
        data_ = np.frombuffer(
            cast(self.data_ptr, POINTER(self.data_ptr_type)).contents,
            dtype =np.float64,
            count = self.total_size,
        ).reshape(self.shape, order='F')
        return data_

    def __repr__(self):
        return f"Array with shape {self.shape}"
    
    def __str__(self):
        string = f"Array with shape {self.shape}\n"
        string += f"Data:\n"
        string += np.array2string(
            self.data, separator=',', prefix='\t',
            formatter={'float': lambda x: "%.3f" % x}
        )
        return string

def benchmark(name, func, iter = 1000, *args, **kwargs):
    import time
    start = time.time()
    for _ in range(iter):
        func(*args, **kwargs)
    end = time.time()
    print(f"-"*100)
    print_green(f"Benchmarking {name}")
    print(f"\tfinished in {(end-start) * 1000} ms for {iter} iterations")
    print(f"\tAverage time: {(end-start) * 1000 / iter} ms")
    print(f"-"*100)


def main():
    # Call the Fortran function
    arr = F_Array((2,2))

    arr2 = F_Array((2,2))
    arr2[:,:] = np.random.randint(0,100, size=(2,2))
    arr[:,:] = np.random.randint(0,100, size=(2,2))

    print("Testing the Fortran API")
    
    print(f"Testing operations between arrays")
    print(f"Array 1:")
    arr.print_in_fortran()
    print(f"Array 2:")
    arr2.print_in_fortran()
    print(f"Adding the arrays")
    res = arr + arr2
    print(f"Result:")
    res.print_in_fortran()

    print(f"SO COOL! Now let's change the data in the array")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            seed = np.random.randint(0,100)
            print(f"Setting data[{i},{j}] = {seed} was {arr.data[i,j]}")
            arr.data[i,j] = seed 
            arr.print_in_fortran()

    # BENCMARKING
    print_IMPORTANT("Starting benchmarking")
    array_shape = (100,10)
    iterations = 1000


    # ARRAY CREATION
    def create_array(shape):
        arr = F_Array(shape)
        del(arr)
    
    def create_np_array(shape):
        arr = np.zeros(shape)
        del(arr)
    
    benchmark("Array creation", create_array, iter=iterations, shape=array_shape)
    benchmark("Numpy array creation", create_np_array, iter=iterations, shape=array_shape)
    
    # ARRAY WRITING
    def assign_values(arr, seed):
        arr.data[:] = seed
        assert np.all(arr.data == seed)
        
    def assign_values_np(arr, seed):
        arr = seed
        assert np.all(arr == seed)
    
    seed = np.random.randint(1,100, size=array_shape)
    benchmark("Fortran array writing", assign_values, iter=iterations, arr=F_Array(array_shape), seed=seed)
    benchmark("Numpy array writing", assign_values_np, iter=iterations, arr=np.zeros(array_shape), seed=seed)
    
    # ARRAY OPERATIONS
    arr1 = F_Array(array_shape)
    arr1[:] = seed
    np_arr1 = np.zeros(array_shape)
    np_arr1[:] = seed

    arr2 = F_Array(array_shape)
    arr2[:] = seed
    np_arr2 = np.zeros(array_shape)
    np_arr2[:] = seed


    def add_arrays(arr1, arr2):
        res = arr1 + arr2
        assert np.all(res.data == 2*seed)
    
    def add_operation():
        arr1._lib_array.add(byref(arr1.to_ctype()), byref(arr2.to_ctype()), byref(arr1.to_ctype()))

    def add_np_arrays(arr1, arr2):
        res = arr1 + arr2
        assert np.all(res == 2*seed)

    benchmark("Fortran array addition", add_arrays, iter=iterations, arr1=arr1, arr2=arr2)
    benchmark("Numpy array addition", add_np_arrays, iter=iterations, arr1=np_arr1, arr2=np_arr2)
    
    def sub_arrays(arr1, arr2):
        res = arr1 - arr2
        assert np.all(res.data == 0)
    
    def sub_np_arrays(arr1, arr2):
        res = arr1 - arr2
        assert np.all(res == 0)

    benchmark("Fortran array subtraction", sub_arrays, iter=iterations, arr1 = arr1 , arr2 = arr2) 
    benchmark("Numpy array subtraction", sub_np_arrays, iter=iterations, arr1=np_arr1, arr2=np_arr2)

    def mul_arrays(arr1, arr2):
        res = arr1 * arr2
        assert np.all(res.data == seed*seed)
    
    def mul_np_arrays(arr1, arr2):
        res = arr1 * arr2
        assert np.all(res == seed*seed)

    benchmark("Fortran array multiplication", mul_arrays, iter=iterations, arr1 = arr1, arr2 = arr2)
    benchmark("Numpy array multiplication", mul_np_arrays, iter=iterations, arr1=np_arr1, arr2=np_arr2)

    def div_arrays(arr1, arr2):
        res = arr1 / arr2
        assert np.all(res.data == 1)
    
    def div_np_arrays(arr1, arr2):
        res = arr1 / arr2
        assert np.all(res == 1)

    benchmark("Fortran array division", div_arrays, iter=iterations, arr1 = arr1, arr2 = arr2)
    benchmark("Numpy array division", div_np_arrays, iter=iterations, arr1=np_arr1, arr2=np_arr2)


if __name__ == "__main__":
    main()