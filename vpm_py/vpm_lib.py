import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL
from .arrays import F_Array_Struct
here = os.path.abspath(os.path.dirname(__file__))
lib_locations = os.path.join(here, 'shared_libs')

lib_vpm_path = glob.glob(os.path.join(lib_locations, 'libvpm_*.so'))[0]
lib_vpm_ext = lib_vpm_path[lib_vpm_path.rfind('.'):]

class VPM_Lib:
    "Singleton class to load the shared library"
    _lib_vpm: CDLL = None 
    _lib_vpm_path = None

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.load_lib()
            cls.instance = super(VPM_Lib, cls).__new__(cls)
        return cls.instance
    
    @classmethod
    def load_lib(cls):
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_vpm_ext)
        tmp.close()
        VPM_Lib._lib_vpm = cdll.LoadLibrary(copy2(lib_vpm_path, tmp.name))
        VPM_Lib._lib_vpm_path = tmp.name
        VPM_Lib._setup_function_signatures()

    @classmethod
    def _setup_function_signatures(cls):
        # API.init
        cls._lib_vpm.init.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int),
            POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_int),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
            POINTER(c_int)
        ]
        cls._lib_vpm.init.restype = None

        # API.vpm
        cls._lib_vpm.vpm.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(F_Array_Struct),
            POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(c_int),
            POINTER(c_double), POINTER(c_int)
        ]
        cls._lib_vpm.vpm.restype = None

        # API.finalize
        cls._lib_vpm.finalize.argtypes = []
        cls._lib_vpm.finalize.restype = None
        
        # API.remesh_particles_3d
        cls._lib_vpm.remesh_particles_3d.argtypes = [
            POINTER(c_int), POINTER(F_Array_Struct), POINTER(F_Array_Struct),
            POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(c_int)
        ]
        cls._lib_vpm.remesh_particles_3d.restype = None
