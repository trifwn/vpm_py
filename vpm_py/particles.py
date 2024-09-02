import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref
import numpy as np

from vpm_py.arrays import F_Array, F_Array_Struct
from vpm_py.vpm_io import print_IMPORTANT
from vpm_py.vpm_lib import VPM_Lib

class Particles:
    def __init__(self, number_equations: int) -> None:
        self.number_equations = number_equations
        self.load_lib()

        self.particle_positions = np.array([])
        self.particle_strengths = np.array([])
        self.particle_velocities = np.array([])
        self.particle_deformations = np.array([])
    
    def load_lib(self)-> None:
        lib = VPM_Lib()
        self._lib_particles = lib._lib_vpm
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
    
    def print_info_fortran(self):
        self._lib_particles.print_particles()
    
    def print_positions_fortran(self):
        self._lib_particles.print_particle_positions()
    
    # def update_positions(self, XP):
    
    @property
    def XP(self):
        """
        Particle positions
        """
        NVR = self.NVR
        XP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        self._lib_particles.get_particle_positions(byref(XP_struct))
        XP_arr = F_Array.from_ctype(XP_struct, name = "from fortran XP")
        # self.particle_positions = XP_arr.data
        return XP_arr

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
        neq = self.number_equations
        QP_struct = F_Array_Struct.null(ndims=2, total_size=(neq + 1)*NVR)
        self._lib_particles.get_particle_strengths(byref(QP_struct))
        QP_arr = F_Array.from_ctype(QP_struct, name ="from fortran QP")
        return QP_arr

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
        UP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        self._lib_particles.get_particle_velocities(byref(UP_struct))
        UP_arr = F_Array.from_ctype(UP_struct,name = "from fortran UP")
        return UP_arr
        
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
        GP_struct = F_Array_Struct.null(ndims=2, total_size = 3*NVR)
        self._lib_particles.get_particle_deformation(byref(GP_struct))
        GP_arr = F_Array.from_ctype(GP_struct,name = "from fortran GP")
        return GP_arr
    
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
    
    def print_data_ownership(self):
        print_IMPORTANT(f"owners of the data:", color_divider="yellow", color_text="yellow")
        print(f"\tXP data owner: {self.particle_positions.__array_interface__['data']}")
        print(f"\tQP data owner: {self.particle_strengths.__array_interface__['data']}")
        print(f"\tUP data owner: {self.particle_velocities.__array_interface__['data']}")
        print(f"\tGP data owner: {self.particle_deformations.__array_interface__['data']}")

    # def __del__(self):
    #     self._lib_particles.finalize()
    #     os.remove(self._lib_particles_path)
