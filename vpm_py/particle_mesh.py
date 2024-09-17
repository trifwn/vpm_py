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
        self.load_lib()
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
        self._lib = lib._lib_vpm
        self.setup_function_signatures()

    def setup_function_signatures(self):
        self._lib.get_NX_pm_coarse.argtypes = [POINTER(c_int)]
        self._lib.get_NX_pm_coarse.restype = None
        self._lib.get_NY_pm_coarse.argtypes = [POINTER(c_int)]
        self._lib.get_NY_pm_coarse.restype = None
        self._lib.get_NZ_pm_coarse.argtypes = [POINTER(c_int)]
        self._lib.get_NZ_pm_coarse.restype = None
        self._lib.get_NN.argtypes = [POINTER(c_int * 3)]
        self._lib.get_NN.restype = None
        self._lib.get_NN_bl.argtypes = [POINTER(c_int * 6)]
        self._lib.get_NN_bl.restype = None
        self._lib.get_Xbound.argtypes = [POINTER(c_double * 6)]
        self._lib.get_Xbound.restype = None
        self._lib.get_neqpm.argtypes = [POINTER(c_int)]
        self._lib.get_neqpm.restype = None
        self._lib.set_RHS_pm.argtypes = [POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
        self._lib.set_RHS_pm.restype = None
        # self._lib.get_dpm.argtypes = [POINTER(F_Array_Struct)]
        # self._lib.get_dpm.restype = None

    @property
    def NX_pm(self):
        NX_pm = c_int()
        self._lib.get_NX_pm_coarse(byref(NX_pm))
        return NX_pm.value
    
    @property
    def NY_pm(self):
        NY_pm = c_int()
        self._lib.get_NY_pm_coarse(byref(NY_pm))
        return NY_pm.value
    
    @property
    def NZ_pm(self):
        NZ_pm = c_int()
        self._lib.get_NZ_pm_coarse(byref(NZ_pm))
        return NZ_pm.value

    @property
    def nn(self):       
        """
            NN(6) is the number of cells in each direction
        """
        NN = (c_int * 3)() 
        self._lib.get_NN(byref(NN))
        return np.array(NN)
    
    @property
    def nn_bl(self):
        """
        NN_bl(6) is the start and end indices of the cells assigned to the processor.
        """
        NN_bl = (c_int * 6)()
        self._lib.get_NN_bl(byref(NN_bl))
        return np.array(NN_bl) 
    
    @property
    def xbound(self):
        """
        Xbound(6) is the boundary of the domain.
        """
        xbound = (c_double * 6)()
        self._lib.get_Xbound(byref(xbound))
        return np.array(xbound)

    @property
    def Xmin(self):
        return self.xbound[0]
    
    @property
    def Xmax(self):
        return self.xbound[3]
    
    @property
    def Ymin(self):
        return self.xbound[1]
    
    @property
    def Ymax(self):
        return self.xbound[4]
    
    @property
    def Zmin(self):
        return self.xbound[2]
    
    @property
    def Zmax(self):
        return self.xbound[5]

    def get_num_equations(self):
        """
            Get the number of equations
        """
        neqpm = c_int()
        self._lib.get_neqpm(byref(neqpm))
        return neqpm.value
   
    def set_rhs_pm(self, RHS_pm: np.ndarray):
        """
        Set the right-hand side of the particle mesh.
        """
        RHS_pm = np.asfortranarray(RHS_pm, dtype=np.float64)
        # Pass as array of shape (num_equations, NXB, NYB, NZB) not pointer
        RHS_pm_ = RHS_pm.ctypes.data_as(POINTER(c_double))
        sizes = np.array(RHS_pm.shape, dtype=np.int32)
        size1 = sizes[0]
        size2 = sizes[1]
        size3 = sizes[2]
        size4 = sizes[3]
        self._lib.set_RHS_pm(
            RHS_pm_, 
            byref(c_int(size1)), 
            byref(c_int(size2)), 
            byref(c_int(size3)), 
            byref(c_int(size4))
        )
    
    def save_to_file(self, folder: str = "results", filename: str = "particle_mesh.dat"):
        """
            Write the particle mesh solution
        """
        self._lib.write_particle_mesh_solution()

    # def get_dpm(self):
    #     """
    #         Get the particle positions
    #     """
    #     dpm = F_Array_Struct.null(ndims=1, total_size = 3)
    #     self._lib.get_dpm(byref(dpm))
    #     return dpm