import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL
from .arrays import F_Array_Struct
here = os.path.abspath(os.path.dirname(__file__))
lib_locations = os.path.join(here, 'shared_libs')

lib_vpm_path = glob.glob(os.path.join(lib_locations, 'vpm_py_api.so'))[0]
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
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), 
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), 
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), 
            POINTER(c_int), POINTER(c_int)
        ]
        cls._lib_vpm.init.restype = None

        # API.vpm
        cls._lib_vpm.vpm.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), 
            POINTER(F_Array_Struct), POINTER(F_Array_Struct), POINTER(c_int),
            POINTER(c_double), POINTER(c_int)
        ]
        cls._lib_vpm.vpm.restype = None

        # Add bindings for vpm_project_solve
        cls._lib_vpm.vpm_project_solve.argtypes = [
            POINTER(c_int),          # NTIME_in
            POINTER(c_double),       # XP_in
            POINTER(c_double),       # QP_in
            POINTER(c_int),          # NVR_in
            POINTER(c_int),          # NVR_size_in
            POINTER(c_int),          # neqpm_in
            POINTER(F_Array_Struct)  # RHS_pm_out
        ]
        cls._lib_vpm.vpm_project_solve.restype = None

        # Add bindings for vpm_define
        cls._lib_vpm.vpm_define.argtypes = [
            POINTER(c_int),    # NTIME_in
            POINTER(c_double), # XP_in
            POINTER(c_double), # QP_in
            POINTER(c_int),    # NVR_in
            POINTER(c_int),    # NVR_size_in
            POINTER(c_int)     # neqpm_in
        ]
        cls._lib_vpm.vpm_define.restype = None

        # Add bindings for vpm_solve_velocity
        cls._lib_vpm.vpm_solve_velocity.argtypes = [
            POINTER(c_int),          # NTIME_in
            POINTER(c_double),       # XP_in
            POINTER(c_double),       # QP_in
            POINTER(c_double),       # UP_in
            POINTER(c_double),       # GP_in
            POINTER(c_int),          # NVR_in
            POINTER(c_int),          # NVR_size_in
            POINTER(c_int),          # neqpm_in
            POINTER(F_Array_Struct), # RHS_pm_out
            POINTER(F_Array_Struct)  # Vel_out
        ]
        cls._lib_vpm.vpm_solve_velocity.restype = None

        # Add bindings for vpm_solve_velocity_deformation
        cls._lib_vpm.vpm_solve_velocity_deformation.argtypes = [
            POINTER(c_int),          # NTIME_in
            POINTER(c_double),       # XP_in
            POINTER(c_double),       # QP_in
            POINTER(c_double),       # UP_in
            POINTER(c_double),       # GP_in
            POINTER(c_int),          # NVR_in
            POINTER(c_int),          # NVR_size_in
            POINTER(c_int),          # neqpm_in
            POINTER(F_Array_Struct), # RHS_pm_out
            POINTER(F_Array_Struct), # Vel_out
            POINTER(F_Array_Struct)  # Deform_out
        ]
        cls._lib_vpm.vpm_solve_velocity_deformation.restype = None

        # Add bindings for vpm_interpolate
        cls._lib_vpm.vpm_interpolate.argtypes = [
            POINTER(c_int),          # NTIME_in
            POINTER(c_double),       # XP_in
            POINTER(c_double),       # QP_in
            POINTER(c_double),       # UP_in
            POINTER(c_double),       # GP_in
            POINTER(c_int),          # NVR_in
            POINTER(c_int),          # NVR_size_in
            POINTER(c_int),          # neqpm_in
            POINTER(F_Array_Struct)  # RHS_pm_out
        ]
        cls._lib_vpm.vpm_interpolate.restype = None

        # Add bindings for vpm_diffuse
        cls._lib_vpm.vpm_diffuse.argtypes = [
            POINTER(c_double),       # NI_in
            POINTER(c_double),       # XP_in
            POINTER(c_double),       # QP_in
            POINTER(c_double),       # UP_in
            POINTER(c_double),       # GP_in
            POINTER(c_int),          # NVR_in
            POINTER(c_int),          # NVR_size_in
            POINTER(c_int),          # neqpm_in
            POINTER(F_Array_Struct)  # RHS_pm_out
        ]
        cls._lib_vpm.vpm_diffuse.restype = None

        # Add bindings for vpm_correct_vorticity
        cls._lib_vpm.vpm_correct_vorticity.argtypes = [
            POINTER(c_double), # XP_in
            POINTER(c_double), # QP_in
            POINTER(c_int),    # NVR_in
            POINTER(c_int),    # neqpm_in
            POINTER(c_int),    # NVR_size_in
        ]
        cls._lib_vpm.vpm_correct_vorticity.restype = None

        # Add bindings for vpm_solve_pressure
        cls._lib_vpm.vpm_solve_pressure.argtypes = [
            POINTER(F_Array_Struct),
            POINTER(F_Array_Struct),
            POINTER(F_Array_Struct),
        ]

        # API.finalize
        cls._lib_vpm.finalize.argtypes = []
        cls._lib_vpm.finalize.restype = None
        
        # API.remesh_particles_3d
        cls._lib_vpm.remesh_particles_3d.argtypes = [
            POINTER(c_int), POINTER(c_int), 
            POINTER(F_Array_Struct), POINTER(F_Array_Struct),
            POINTER(F_Array_Struct), POINTER(F_Array_Struct), 
            POINTER(c_int), POINTER(c_double)
        ]
        cls._lib_vpm.remesh_particles_3d.restype = None

        cls._lib_vpm.write_particle_mesh_solution.argtypes = []
        cls._lib_vpm.write_particle_mesh_solution.restype = None
        cls._lib_vpm.write_particles.argtypes = []
        cls._lib_vpm.write_particles.restype = None

        # # Calculate Derivatives
        # cls._lib_vpm.calc_derivative.argtypes = [
        #     POINTER(F_Array_Struct), POINTER(c_double), POINTER(c_int), POINTER(c_int),
        #     POINTER(F_Array_Struct)
        # ]
        # cls._lib_vpm.calc_derivative.restype = None
