import os
import glob
from shutil import copy2
from tempfile import NamedTemporaryFile
from ctypes import c_double, c_int, POINTER, cdll, CDLL, byref
import numpy as np

from vpm_py.arrays import F_Array, F_Array_Struct
from vpm_py.console_io import print_IMPORTANT
from vpm_py.vpm_lib import VPM_Lib


class ParticleMesh: 
    def __init__(
        self, 
        num_equations: int = 3
    ) -> None:
        self.load_lib()
        self.number_equations = num_equations
        
        self.load_lib()
        self.U = np.zeros((3, 1, 1, 1))
        self.deformation = np.zeros((3, 1, 1, 1))
        self.pressure = np.zeros((1, 1, 1))
        self.q_pressure = np.zeros((1, 1, 1))
        self.u_pressure = np.zeros((1, 1, 1))
        self.SOL = np.zeros((num_equations, 1, 1, 1))
        self.RHS = np.zeros((num_equations , 1, 1, 1))
    
    def load_lib(self):
        lib = VPM_Lib()
        self._lib = lib._lib_vpm
        self.setup_function_signatures()

    def setup_function_signatures(self):
        self._lib.get_NX_pm.argtypes = [POINTER(c_int)]
        self._lib.get_NX_pm.restype = None
        self._lib.get_NY_pm.argtypes = [POINTER(c_int)]
        self._lib.get_NY_pm.restype = None
        self._lib.get_NZ_pm.argtypes = [POINTER(c_int)]
        self._lib.get_NZ_pm.restype = None
        self._lib.get_NN.argtypes = [POINTER(c_int * 3)]
        self._lib.get_NN.restype = None
        self._lib.get_NN_bl.argtypes = [POINTER(c_int * 6)]
        self._lib.get_NN_bl.restype = None
        self._lib.get_fine_bounds.argtypes = [POINTER(c_double * 6)]
        self._lib.get_fine_bounds.restype = None
        self._lib.get_Xbounds.argtypes = [POINTER(c_double), POINTER(c_double)]
        self._lib.get_Xbounds.restype = None
        self._lib.get_Ybounds.argtypes = [POINTER(c_double), POINTER(c_double)]
        self._lib.get_Ybounds.restype = None
        self._lib.get_Zbounds.argtypes = [POINTER(c_double), POINTER(c_double)]
        self._lib.get_Zbounds.restype = None
        self._lib.get_neqpm.argtypes = [POINTER(c_int)]
        self._lib.get_neqpm.restype = None
        self._lib.set_RHS_pm.argtypes = [POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
        self._lib.set_RHS_pm.restype = None
        # self._lib.get_dpm.argtypes = [POINTER(F_Array_Struct)]
        # self._lib.get_dpm.restype = None
    
    @property
    def grid_positions(self):

        # First, we need to get the number of grid points in each direction
        NX_pm = self.NX_pm
        NY_pm = self.NY_pm
        NZ_pm = self.NZ_pm

        Xmin = c_double()
        Xmax = c_double()
        self._lib.get_Xbounds(byref(Xmin), byref(Xmax))
        Xmin = Xmin.value
        Xmax = Xmax.value

        Ymin = c_double()
        Ymax = c_double()
        self._lib.get_Ybounds(byref(Ymin), byref(Ymax))
        Ymin = Ymin.value
        Ymax = Ymax.value

        Zmin = c_double()
        Zmax = c_double()
        self._lib.get_Zbounds(byref(Zmin), byref(Zmax))
        Zmin = Zmin.value
        Zmax = Zmax.value

        # Get dx, dy, dz
        dx = (Xmax - Xmin) / (NX_pm - 1)
        dy = (Ymax - Ymin) / (NY_pm - 1)
        dz = (Zmax - Zmin) / (NZ_pm - 1)
        
        X = np.linspace(Xmin, Xmax, NX_pm)
        Y = np.linspace(Ymin, Ymax, NY_pm)
        Z = np.linspace(Zmin, Zmax, NZ_pm)
        X,Y,Z = np.meshgrid(X, Y, Z, indexing='ij')
        return np.array([X, Y, Z])

    @property
    def NX_pm(self):
        NX_pm = c_int()
        self._lib.get_NX_pm(byref(NX_pm))
        return NX_pm.value
    
    @property
    def NY_pm(self):
        NY_pm = c_int()
        self._lib.get_NY_pm(byref(NY_pm))
        return NY_pm.value
    
    @property
    def NZ_pm(self):
        NZ_pm = c_int()
        self._lib.get_NZ_pm(byref(NZ_pm))
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
        self._lib.get_fine_bounds(byref(xbound))
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
    
    def save_to_file(self,  filename: str = "particle_mesh", folder: str = "results", filetype: str = "h5"):
        """
            Write the particle mesh solution

            Parameters
            ----------
            filename : str
                The filename to write to
            folder : str
                The folder to write to
            filetype : str
                The filetype to write to
        """
        if not folder.endswith("/"):
            folder += "/"

        filename_b = filename.encode('utf-8')
        folder_b = folder.encode('utf-8')

        os.makedirs(folder, exist_ok=True)
        os.makedirs(os.path.join(folder,'results'), exist_ok=True)

        if filetype == "h5":
            self._lib.write_particle_mesh_solution_hdf5(folder_b, filename_b)
        elif filetype == "dat":
            self._lib.write_particle_mesh_solution(folder_b, filename_b)
        else:
            raise ValueError(f"Filetype {filetype} not recognized")

    def save_pressure_to_file(self, filename: str = "pressure", folder: str = "results", filetype: str = "h5"):
        """
            Write the pressure solution

            Parameters
            ----------
            filename : str
                The filename to write to
            folder : str
                The folder to write to
            filetype : str
                The filetype to write to
        """
        # if not folder.endswith("/"):
        #     folder += "/"

        # filename_b = filename.encode('utf-8')
        # folder_b = folder.encode('utf-8')

        # os.makedirs(folder, exist_ok=True)
        # os.makedirs(os.path.join(folder,'results'), exist_ok=True)

        # if filetype == "h5":
        #     # Convert pressure to F_Array
        #     pressure = F_Array.from_ndarray(self.pressure)
        #     self._lib.write_pressure_hdf5(folder_b, filename_b, byref(pressure.to_ctype()))
        # else:
        #     raise ValueError(f"Filetype {filetype} not recognized")

        import h5py as h5
        if not folder.endswith("/"):
            folder += "/"
        filename = os.path.join(folder, filename)
        with h5.File(filename, 'w') as f:
            f.create_dataset(filename, data=self.pressure)

        

    # def get_dpm(self):
    #     """
    #         Get the particle positions
    #     """
    #     dpm = F_Array_Struct.null(ndims=1, total_size = 3)
    #     self._lib.get_dpm(byref(dpm))
    #     return dpm