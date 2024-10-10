import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref
from .vpm_dtypes import dp_array_to_pointer
from .arrays import F_Array_Struct, F_Array
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))
lib_locations = os.path.join(here, 'shared_libs')

lib_vpm_path = glob.glob(os.path.join(lib_locations, 'liboperators*.so'))[0]
lib_vpm_ext = lib_vpm_path[lib_vpm_path.rfind('.'):]

class OperatorsLib:
    "Singleton class to load the shared library"
    _lib: CDLL = None 
    _lib_path = None

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.load_lib()
            cls.instance = super(OperatorsLib, cls).__new__(cls)
        return cls.instance
    
    @classmethod
    def load_lib(cls):
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_vpm_ext)
        tmp.close()
        OperatorsLib._lib = cdll.LoadLibrary(copy2(lib_vpm_path, tmp.name))
        OperatorsLib._lib_path = tmp.name
        OperatorsLib._setup_function_signatures()

    @classmethod
    def _setup_function_signatures(cls):
        # Calculate Derivatives
        cls._lib.derivative.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(c_int),
            POINTER(F_Array_Struct)
        ]
        cls._lib.derivative.restype = None

        # Calculate Divergence
        cls._lib.divergence.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.divergence.restype = None

        # Calculate Curl
        cls._lib.curl.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.curl.restype = None

        # Calculate Laplacian
        cls._lib.laplacian.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.laplacian.restype = None

        # Calculate Gradient
        cls._lib.gradient.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.gradient.restype = None

        # Calculate Vector Laplacian
        cls._lib.vector_laplacian.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.vector_laplacian.restype = None

        # Calculate Hessian
        cls._lib.hessian.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.hessian.restype = None

        # Calculate Jacobian
        cls._lib.jacobian.argtypes = [
            POINTER(F_Array_Struct), POINTER(c_int), POINTER(c_double),
            POINTER(F_Array_Struct)
        ]
        cls._lib.jacobian.restype = None

    @classmethod
    def calc_derivative(cls, field, dx, dy, dz, order, direction):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)

        result = F_Array_Struct.null(ndims=4, total_size=1)

        cls._lib.derivative(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, 
            byref(c_int(order)), byref(c_int(direction)), byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_divergence(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)

        result = F_Array_Struct.null(ndims=4, total_size=1)

        cls._lib.divergence(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_curl(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)

        result = F_Array_Struct.null(ndims=4, total_size=1)

        cls._lib.curl(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_laplacian(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        
        dpm_pointer = dp_array_to_pointer(np.array([dx, dy, dz]))
        result = F_Array_Struct.null(ndims=field_f.ndims , total_size=1)
        
        cls._lib.laplacian(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_gradient(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)

        result = F_Array_Struct.null(ndims=3, total_size=1)
        cls._lib.gradient(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_vector_laplacian(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)
        original_dims = field_f.ndims
        result = F_Array_Struct.null(ndims= original_dims, total_size=1)
        cls._lib.vector_laplacian(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_hessian(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)
        original_dims = field_f.ndims
        result = F_Array_Struct.null(ndims= original_dims + 2, total_size=1) 
        cls._lib.hessian(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np

    @classmethod
    def calc_jacobian(cls, field, dx, dy, dz):
        field = np.asfortranarray(field)
        field_f = F_Array(shape=field.shape, data_container=field)
        dpm = np.array([dx, dy, dz])
        dpm_pointer = dp_array_to_pointer(dpm)

        result = F_Array_Struct.null(ndims= field_f.ndims + 1, total_size=1) 
        cls._lib.jacobian(
            byref(field_f.to_ctype()), byref(c_int(3)), dpm_pointer, byref(result)
        )
        result_arr = F_Array.from_ctype(result, ownership=True)
        result_np = result_arr.transfer_data_ownership()
        return result_np