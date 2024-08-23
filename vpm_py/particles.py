import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref, cast
from .arrays import F_Array, F_Array_Struct
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))
lib_locations = os.path.join(here, 'shared_libs')

lib_vpm_path = glob.glob(os.path.join(lib_locations, 'libvpm_*.so'))[0]
lib_vpm_ext = lib_vpm_path[lib_vpm_path.rfind('.'):]


class Particles:
    def __init__(self, number_equations: int) -> None:
        self.number_equations = number_equations
        self.load_lib()
    
    def load_lib(self):
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_vpm_ext)
        tmp.close()
        self._lib_particles: CDLL = cdll.LoadLibrary(copy2(lib_vpm_path, tmp.name))
        self._lib_particles_path = tmp.name
        self.setup_function_signatures()

    def setup_function_signatures(self):
        self._lib_particles.get_num_particles.argtypes = [POINTER(c_int)]
        self._lib_particles.get_num_particles.restype = None
        self._lib_particles.get_particle_positions.argtypes = [POINTER(F_Array_Struct)]
        self._lib_particles.get_particle_positions.restype = None
        self._lib_particles.get_particle_strengths.argtypes = [POINTER(F_Array_Struct)]
        self._lib_particles.get_particle_strengths.restype = None
        self._lib_particles.get_particle_velocities.argtypes = [POINTER(F_Array_Struct)]
        self._lib_particles.get_particle_velocities.restype = None
        self._lib_particles.get_particle_deformation.argtypes = [POINTER(F_Array_Struct)]
        self._lib_particles.get_particle_deformation.restype = None
        self._lib_particles.set_particle_positions.argtypes = [POINTER(c_double)]
        self._lib_particles.set_particle_positions.restype = None
        self._lib_particles.set_particle_strengths.argtypes = [POINTER(c_double)]
        self._lib_particles.set_particle_strengths.restype = None
        self._lib_particles.set_particle_velocities.argtypes = [POINTER(c_double)]
        self._lib_particles.set_particle_velocities.restype = None
        self._lib_particles.set_particle_deformation.argtypes = [POINTER(c_double)]
        self._lib_particles.set_particle_deformation.restype = None
        self._lib_particles.print_particles.argtypes = []
        self._lib_particles.print_particles.restype = None
    
    def print_particles(self):
        self._lib_particles.print_particles()
    
    @property
    def XP(self):
        """
        Particle positions
        """
        NVR = self.NVR
        XP_arr = F_Array((3, NVR))
        XP_struct = XP_arr.to_ctype()
        self._lib_particles.get_particle_positions(byref(XP_struct))
        XP_arr = F_Array.from_ctype(XP_struct)
        return XP_arr.data

    @XP.setter
    def XP(self, XP):
        XP = np.ascontiguousarray(XP, dtype=np.float64)
        XP_ptr = XP.ctypes.data_as(POINTER(c_double))
        self._lib_particles.set_particle_positions(XP_ptr)

    @property
    def QP(self):
        """
        Particle strengths
        """
        NVR = self.NVR
        QP_arr = F_Array((self.number_equations + 1, NVR))
        QP_struct = QP_arr.to_ctype()
        self._lib_particles.get_particle_strengths(byref(QP_struct))
        QP_arr = F_Array.from_ctype(QP_struct)
        return QP_arr.data
        

    @QP.setter
    def QP(self, QP):
        QP = np.ascontiguousarray(QP, dtype=np.float64)
        QP_ptr = QP.ctypes.data_as(POINTER(c_double))
        self._lib_particles.set_particle_strengths(QP_ptr)

    @property
    def UP(self):
        """
            Particle velocities
        """
        NVR = self.NVR
        UP_arr = F_Array((3, NVR))
        UP_struct = UP_arr.to_ctype()
        self._lib_particles.get_particle_velocities(byref(UP_struct))
        UP_arr = F_Array.from_ctype(UP_struct)
        return UP_arr.data
        
    @UP.setter
    def UP(self, UP):
        UP = np.ascontiguousarray(UP, dtype=np.float64)
        UP_ptr = UP.ctypes.data_as(POINTER(c_double))
        self._lib_particles.set_particle_velocities(UP_ptr)
    
    @property
    def GP(self):
        """
            Particle deformations
        """
        NVR = self.NVR
        GP_arr = F_Array((3, NVR))
        GP_struct = GP_arr.to_ctype()
        self._lib_particles.get_particle_deformation(byref(GP_struct))
        GP_arr = F_Array.from_ctype(GP_struct)
        return GP_arr.data
    
    @GP.setter
    def GP(self, GP):
        GP = np.ascontiguousarray(GP, dtype=np.float64)
        GP_ptr = GP.ctypes.data_as(POINTER(c_double))
        self._lib_particles.set_particle_deformation(GP_ptr)

    @property
    def NVR(self):
        NVR = c_int()
        self._lib_particles.get_num_particles(byref(NVR))
        return NVR.value
    
    def __del__(self):
        self._lib_particles.finalize()
        os.remove(self._lib_particles_path)
