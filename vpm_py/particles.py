import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref, c_char_p, create_string_buffer
import numpy as np

from vpm_py.arrays import F_Array, F_Array_Struct
from vpm_py.console_io import print_IMPORTANT
from vpm_py.vpm_lib import VPM_Lib

class Particles:
    def __init__(self, number_equations: int) -> None:
        self.number_equations = number_equations
        self.load_lib()

        self.particle_positions = np.array([])
        self.particle_charges = np.array([])
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
        Particle Charges
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
    
    def save_to_file(self, filename: str = "particles", folder: str = "results", filetype: str = "h5"):
        """
        Save the particles to a file named: folder/____filename where the ____ indicates the timestep

        Args:
            folder (str, optional): Folder to save the file. Defaults to "results".
            filename (str, optional): Filename. Defaults to "particles".
            filetype (str, optional): Filetype. Defaults to hdf5 format "h5" .
        """
        if not folder.endswith("/"):
            folder += "/"

        filename_b = filename.encode('utf-8')
        folder_b = folder.encode('utf-8')

        # If folder does not exist, create it
        os.makedirs(folder, exist_ok=True)

        # Call the Fortran routine using ctypes
        if filetype == "h5":
            print(f"\tSaving particles to HDF5 file")
            self._lib_particles.write_particles_hdf5(folder_b, filename_b)
        elif filetype == "dat":
            print(f"\tSaving particles to formatted file")
            self._lib_particles.write_particles(folder_b, filename_b)
        else:
            raise ValueError(f"Invalid filetype: {filetype}")