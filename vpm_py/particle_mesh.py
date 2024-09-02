import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref
import numpy as np

from vpm_py.arrays import F_Array, F_Array_Struct
from vpm_py.vpm_io import print_IMPORTANT
from vpm_py.vpm_lib import VPM_Lib


class ParticleMesh: 
    def __init__(
        self, 
        Nx: int= 10,
        Ny: int= 10,
        Nz: int= 10,
        num_equations: int = 3
    ) -> None:

        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.number_equations = num_equations
        
        self.load_lib()
        self.Ux = np.zeros((Nx, Ny, Nz))
        self.Uy = np.zeros((Nx, Ny, Nz))
        self.Uz = np.zeros((Nx, Ny, Nz))
        self.RHS = np.zeros((num_equations + 1, Nx, Ny, Nz))
    
    def update_U(self, Ux, Uy, Uz):
        self.Ux = Ux
        self.Uy = Uy
        self.Uz = Uz
    
    def load_lib(self):
        lib = VPM_Lib()
        self._lib_pmgrid = lib._lib_vpm
        self.setup_function_signatures()

    def setup_function_signatures(self):
        pass
   