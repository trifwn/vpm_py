import os
from ctypes import c_int, POINTER, byref, _Pointer
import numpy as np

from vpm_py.arrays import F_Array, F_Array_Struct
from vpm_py.vpm_lib import VPM_Lib
from .vpm_dtypes import dp_array_to_pointer

from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()

class Particles:
    def __init__(self, number_equations: int) -> None:
        self.number_equations = number_equations
        self.load_lib()

        # Assing 2D numpy arrays to store particle data
        self._particle_positions = np.zeros((3, 0), dtype=np.float64, order='F')
        self._particle_charges = np.zeros((self.number_equations + 1, 0), dtype=np.float64, order='F')
        self._particle_velocities = np.zeros((3, 0), dtype=np.float64, order='F') 
        self._particle_deformations = np.zeros((3, 0), dtype=np.float64, order='F')

    
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
        self._lib_particles.print_particles.argtypes = []
        self._lib_particles.print_particles.restype = None
    
    def print_info_fortran(self):
        self._lib_particles.print_particles()
    
    def print_positions_fortran(self):
        self._lib_particles.print_particle_positions()
    
    @property
    def particle_positions(self):
        return self._particle_positions
    
    @particle_positions.setter
    def particle_positions(self, value: np.ndarray):
        self._particle_positions = self._convert_to_numpy(value)

    @property
    def particle_charges(self):
        return self._particle_charges
    
    @particle_charges.setter
    def particle_charges(self, value: np.ndarray):
        self._particle_charges = self._convert_to_numpy(value)

    @property
    def particle_velocities(self):
        return self._particle_velocities
    
    @particle_velocities.setter
    def particle_velocities(self, value: np.ndarray):
        self._particle_velocities = self._convert_to_numpy(value)

    @property
    def particle_deformations(self):
        return self._particle_deformations
    
    @particle_deformations.setter
    def particle_deformations(self, value: np.ndarray):
        self._particle_deformations = self._convert_to_numpy(value)

    @property
    def XP(self):
        """
        Particle positions
        """
        NVR = self.NVR
        XP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        self._lib_particles.get_particle_positions(byref(XP_struct))
        XP_arr = F_Array.from_ctype(XP_struct)
        return XP_arr

    @property
    def QP(self):
        """
        Particle Charges
        """
        NVR = self.NVR
        neq = self.number_equations
        QP_struct = F_Array_Struct.null(ndims=2, total_size=(neq + 1)*NVR)
        self._lib_particles.get_particle_strengths(byref(QP_struct))
        QP_arr = F_Array.from_ctype(QP_struct)
        return QP_arr

    @property
    def UP(self):
        """
            Particle velocities
        """
        NVR = self.NVR
        UP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        self._lib_particles.get_particle_velocities(byref(UP_struct))
        UP_arr = F_Array.from_ctype(UP_struct)
        return UP_arr
        
    @property
    def GP(self):
        """
            Particle deformations
        """
        NVR = self.NVR
        GP_struct = F_Array_Struct.null(ndims=2, total_size = 3*NVR)
        self._lib_particles.get_particle_deformation(byref(GP_struct))
        GP_arr = F_Array.from_ctype(GP_struct)
        return GP_arr
    
    @property
    def NVR(self):
        NVR = c_int()
        self._lib_particles.get_num_particles(byref(NVR))
        return NVR.value
    
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
    
    def get_array_pointers(self, return_all = False) -> list[_Pointer]:
        """Create pointers for arrays with copy."""
        arrays = [
            self.particle_positions,
            self.particle_charges,
        ]
        if return_all:
            arrays.extend([
                self.particle_velocities,
                self.particle_deformations
            ])
        return [dp_array_to_pointer(arr, copy=False) for arr in arrays]

    def validate_stored_particle_counts(self, NVR: int | None = None) -> None:
        """Validate particle counts match across arrays and return NVR_size."""
        if NVR is None:
            NVR = self.NVR

        positions = self.particle_positions
        charges = self.particle_charges
        velocities = self.particle_velocities
        deformations = self.particle_deformations

        if not (NVR == positions.shape[1]):
            print(f"Warning: NVR size mismatch: {NVR} != {positions.shape[1]}, for positions")
            self.particle_positions = np.zeros((3, NVR), dtype=np.float64, order='F')
        
        if not (NVR == charges.shape[1]):
            print(f"Warning: NVR size mismatch: {NVR} != {charges.shape[1]}, for charges")
            self.particle_charges = np.zeros((self.number_equations + 1, NVR), dtype=np.float64, order='F')

        if velocities is not None and not (NVR == velocities.shape[1]):
            print(f"Warning: NVR size mismatch: {NVR} != {velocities.shape[1]}, for velocities")
            self.particle_velocities = np.zeros((3, NVR), dtype=np.float64, order='F')

        if deformations is not None and not (NVR == deformations.shape[1]):
            print(f"Warning: NVR size mismatch: {NVR} != {deformations.shape[1]} for deformations")
            self.particle_deformations = np.zeros((3, NVR), dtype=np.float64, order='F')

    def store_particles(
            self, 
            positions: np.ndarray | F_Array | None = None, 
            charges: np.ndarray | F_Array | None = None, 
            velocities: np.ndarray | F_Array | None = None, 
            deformations: np.ndarray | F_Array | None = None
        ) -> None:
        """Store particle results back to class attributes."""
        if positions is not None:
            del self._particle_positions 
            self.particle_positions = positions

        if charges is not None:
            del self._particle_charges
            self.particle_charges = charges

        if velocities is not None:
            del self._particle_velocities
            self.particle_velocities = velocities

        if deformations is not None:
            del self._particle_deformations
            self.particle_deformations = deformations

        if rank == 0:
            positions = self.particle_positions
            print(f"Stored particles: {self.NVR} == {positions.shape[1]}")
            mystring = f"Particle positions: type: {type(positions)} "
            mystring += f"shape: {positions.shape}\n"
            if positions.size > 1:
                mystring += f"\tmax: {np.max(positions[:]):.6f}, min: {np.min(positions[:]):.6f} size in MB: {positions.nbytes/1024/1024:.6f}"
            print(mystring)

            charges = self.particle_charges
            mystring = f"Particle charges: type: {type(charges)} "
            mystring += f"shape: {charges.shape}\n"
            if charges.size > 1:
                mystring += f"\tmax: {np.max(charges[:]):.6f}, min: {np.min(charges[:]):.6f} size in MB: {charges.nbytes/1024/1024:.6f}"
            print(mystring)

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
        os.makedirs(os.path.join(folder,'results'), exist_ok=True)

        # Call the Fortran routine using ctypes
        if filetype == "h5":
            self._lib_particles.write_particles_hdf5(folder_b, filename_b)
        elif filetype == "dat":
            self._lib_particles.write_particles(folder_b, filename_b)
        else:
            raise ValueError(f"Invalid filetype: {filetype}")