import os
from ctypes import c_double, c_int, POINTER, byref
import numpy as np
from vpm_py.vpm_lib import VPM_Lib
from vpm_py.arrays import F_Array


class ParticleMesh: 
    def __init__(
        self, 
        num_equations: int = 3
    ) -> None:
        self.load_lib()
        self.number_equations = num_equations
        
        self.load_lib()
        self._velocity = np.zeros((3,0,0,0), dtype=np.float64, order='F')
        self._deformation = np.zeros((3,0,0,0), dtype=np.float64, order='F')
        self._pressure = np.zeros((0,0,0), dtype=np.float64, order='F')
        self._q_pressure = np.zeros((0,0,0), dtype=np.float64, order='F')
        self._u_pressure = np.zeros((0,0,0,0), dtype=np.float64, order='F')
        self._SOL = np.array([])
        self._RHS = np.array([])
    
    @property
    def velocity(self):
        if self._velocity is None or self._velocity.size == 0:
            return None
        return self._velocity
    
    @velocity.setter
    def velocity(self, value: np.ndarray):
        self._velocity = self._convert_to_numpy(value) 
    
    @property
    def deformation(self):
        if self._deformation is None or self._deformation.size == 0:
            return None
        return self._deformation
    
    @deformation.setter
    def deformation(self, value: np.ndarray):
        self._deformation = self._convert_to_numpy(value)
    
    @property
    def pressure(self):
        if self._pressure is None or self._pressure.size == 0:
            return None
        return self._pressure
    
    @pressure.setter
    def pressure(self, value: np.ndarray):
        self._pressure = self._convert_to_numpy(value)
    
    @property
    def q_pressure(self):
        if self._q_pressure is None or self._q_pressure.size == 0:
            return None
        return self._q_pressure
    
    @q_pressure.setter
    def q_pressure(self, value: np.ndarray):
        self._q_pressure = self._convert_to_numpy(value)
    
    @property
    def u_pressure(self):
        if self._u_pressure is None or self._u_pressure.size == 0: 
            return None
        return self._u_pressure
    
    @u_pressure.setter
    def u_pressure(self, value: np.ndarray):
        self._u_pressure = self._convert_to_numpy(value)

    @property
    def SOL(self):
        if self._SOL is None or self._SOL.size == 0:
            return None
        return self._SOL
    
    @SOL.setter
    def SOL(self, value: np.ndarray):
        self._SOL = self._convert_to_numpy(value)
    
    @property
    def RHS(self):
        if self._RHS is None or self._RHS.size == 0: 
            return None
        return self._RHS
    
    @RHS.setter
    def RHS(self, value: np.ndarray):
        self._RHS = self._convert_to_numpy(value)
    
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
        self._lib.get_dpm.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self._lib.get_dpm.restype = None
    
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
        
        X = np.linspace(Xmin, Xmax, NX_pm)
        Y = np.linspace(Ymin, Ymax, NY_pm)
        Z = np.linspace(Zmin, Zmax, NZ_pm)
        X,Y,Z = np.meshgrid(X, Y, Z, indexing='ij')
        return np.array([X, Y, Z])

    @property
    def grid_spacing(self):
        dx = c_double()
        dy = c_double()
        dz = c_double()
        self._lib.get_dpm(byref(dx), byref(dy), byref(dz))
        return np.array([dx.value, dy.value, dz.value])

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
    def grid_size(self):       
        """
            NN(3) is the number of cells in each direction
        """
        NN = (c_int * 3)() 
        self._lib.get_NN(byref(NN))
        return np.array(NN)
    
    @property
    def grid_limits(self):
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
    
    def store_mesh_results(
        self,
        rhs: np.ndarray | F_Array | None = None,
        sol: np.ndarray | F_Array |  None = None,
        velocity: np.ndarray | F_Array  | None = None,
        deformation: np.ndarray | F_Array | None = None,
        pressures: np.ndarray | F_Array |  None = None,
    ):
        if rhs is not None:
            if isinstance(rhs, F_Array):
                rhs = rhs.to_numpy()
            if not np.shares_memory(rhs, self._RHS):
                # del self._RHS
                self.RHS = rhs

        if sol is not None:
            if isinstance(sol, F_Array):
                sol = sol.to_numpy()
            if not np.shares_memory(sol, self._SOL):
                # del self._SOL
                self.SOL = sol

        if velocity is not None:
            if not np.shares_memory(velocity, self._velocity):
                # del self._velocity
                self.velocity = velocity

        if deformation is not None:
            if isinstance(deformation, F_Array):
                deformation = deformation.to_numpy()
            if not np.shares_memory(deformation, self._deformation):
                # del self._deformation
                self.deformation = deformation

        if pressures is not None:
            if isinstance(pressures, F_Array):
                pressures = pressures.to_numpy()
                u_pressure = pressures[0, :, :, :]
                q_pressure = pressures[1, :, :, :]
                pressure   = pressures[2, :, :, :]

            if not np.shares_memory(pressure, self._pressure):
                del self._pressure
                self.pressure = pressure
            
            if not np.shares_memory(u_pressure, self._u_pressure):
                del self._u_pressure
                self.u_pressure = u_pressure
            
            if not np.shares_memory(q_pressure, self._q_pressure):
                del self._q_pressure
                self.q_pressure = q_pressure
 
    def _convert_to_numpy(self, array: np.ndarray | F_Array) -> np.ndarray:
        """Convert F_Array to numpy array if needed."""
        if isinstance(array, F_Array):
            array = array.data
        
        if not isinstance(array, np.ndarray):
            raise ValueError("Invalid array type")
        
        # Make sure the array has Fortran order
        if not array.flags['F_CONTIGUOUS']:
            array = np.asfortranarray(array, dtype=np.float64)
        return array 
    
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
        import h5py as h5
        if not folder.endswith("/"):
            folder += "/"
        filename = os.path.join(folder, filename)
        with h5.File(filename, 'w') as f:
            f.create_dataset(filename, data=self.pressure)
