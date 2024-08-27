import numpy as np
import os

from ctypes import c_int,  byref, POINTER, cdll, c_double, cast
from shutil import copy2
from tempfile import NamedTemporaryFile
from mpi4py import MPI

# Local imports
from . import ParticleMesh
from . import Particles
from .vpm_io import  print_IMPORTANT, print_green, print_blue, print_red
from .vpm_lib import VPM_Lib
from .utils import divide_processors
from .arrays import F_Array
from .vpm_dtypes import dp_array_to_pointer, pointer_to_dp_array

class VPM(object):
    """
    Interface to the VPM Fortran routines.
    """

    def __init__(
        self,
        number_of_equations: int = 3,
        number_of_processors: int = 1,
        max_particle_num: int = 1000,
        rank: int = 0
    ):
        lib = VPM_Lib()
        self._lib = lib._lib_vpm
        self.vpm_lib = lib

        self.rank = rank
        self.num_processors = number_of_processors
        self.max_particle_num = max_particle_num
        
        # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
        NBI, NBJ, NBK = divide_processors(number_of_processors)
        if self.rank == 0:
            print(f"Number of processors: {self.num_processors}")
            print(f"NBI: {NBI}, NBJ: {NBJ}, NBK: {NBK}")
        
        self.dpm = np.array([4.0, 4.0, 4.0])
        self.NBI = NBI
        self.NBJ = NBJ
        self.NBK = NBK
        self.num_equations = number_of_equations

        # Initialize the VPM
        self.initialize(self.dpm[0], self.dpm[1], self.dpm[2], NBI, NBJ, NBK)
        self.particle_mesh = ParticleMesh()
        self.particles = Particles(number_equations= self.num_equations)

    def initialize(
        self, 
        DXpm: float, 
        DYpm: float, 
        DZpm: float, 
        NBI: int, 
        NBJ: int, 
        NBK: int, 
        ncoarse = 8, 
        projection_type: int = 4, 
        boundary_condition: int =2, 
        variable_volume: bool = True, 
        EPSVOL: float = 1e-6, 
        remesh: bool = True, 
        ncell_rem: int = 1, 
        use_tree: bool = True, 
        ilevmax: int = 1,
        OMPTHREADS: int = 1,
        is_box_fixed: bool = False, 
        slice: bool = True, 
        IPMWRITE: int = 1, 
        IPMWSTART: int = 1, 
        IPMWSTEPS: int = 150
    ):
        """Wrapper function for calling Fortran subroutine `init`.

        Args:
            DXpm (float): Grid spacing in x-direction.
            DYpm (float): Grid spacing in y-direction.
            DZpm (float): Grid spacing in z-direction.
            NBI (int): Number of cells in x-direction.
            NBJ (int): Number of cells in y-direction.
            NBK (int): Number of cells in z-direction.
            ncoarse (int, optional): Number of coarse cells. Defaults to 8.
            projection_type (int, optional): Projection type. Defaults to 4.
            boundary_condition (int, optional): Boundary condition. Defaults to 2.
            variable_volume (bool, optional): Variable volume. Defaults to True.
            EPSVOL (float, optional): Volume tolerance. Defaults to 1e-6.
            remesh (bool, optional): Remesh. Defaults to True.
            ncell_rem (int, optional): Number of cells to remesh. Defaults to 1.
            use_tree (bool, optional): Use tree. Defaults to True.
            ilevmax (int, optional): Maximum level. Defaults to 1.
            OMPTHREADS (int, optional): Number of threads. Defaults to 1.
            is_box_fixed (bool, optional): Is box fixed. Defaults to True.
            slice (bool, optional): Slice. Defaults to True.
            IPMWRITE (int, optional): Write IPM. Defaults to 1.
            IPMWSTART (int, optional): Write IPM start. Defaults to 1.
            IPMWSTEPS (int, optional): Write IPM steps. Defaults to 1. 
        """
        self.dpm = np.array([DXpm, DYpm, DZpm])
        self._lib.init(
            byref(c_double(DXpm)), byref(c_double(DYpm)), byref(c_double(DZpm)), byref(c_int(projection_type)),
            byref(c_int(boundary_condition)), byref(c_int(variable_volume)), byref(c_double(EPSVOL)), byref(c_int(ncoarse)),
            byref(c_int(NBI)), byref(c_int(NBJ)), byref(c_int(NBK)), byref(c_int(remesh)), byref(c_int(ncell_rem)),
            byref(c_int(use_tree)), byref(c_int(ilevmax)), byref(c_int(OMPTHREADS)), byref(c_int(is_box_fixed)),
            byref(c_int(slice)), byref(c_int(IPMWRITE)), byref(c_int(IPMWSTART)), byref(c_int(IPMWSTEPS))
        )
        if(self.rank == 0):
            print_green(f"Finished initializing VPM {self.rank}:")
            # Print the arguments passed
            print(f"\tDXpm= {DXpm}")
            print(f"\tDYpm= {DYpm}")
            print(f"\tDZpm= {DZpm}")
            print(f"\tinterf_iproj= {projection_type}")
            print(f"\tibctyp= {boundary_condition}")
            print(f"\tIDVPM= {variable_volume}")
            print(f"\tEPSVOL= {EPSVOL}")
            print(f"\tncoarse= {ncoarse}")
            print(f"\tNBI= {NBI}")
            print(f"\tNBJ= {NBJ}")
            print(f"\tNBK= {NBK}")
            print(f"\tNREMESH= {remesh}")
            print(f"\tncell_rem= {ncell_rem}")
            print(f"\tiyntree= {use_tree}")
            print(f"\tilevmax= {ilevmax}")
            print(f"\tOMPTHREADS= {OMPTHREADS}")
            print(f"\tidefine= {is_box_fixed}")
            print(f"\tiynslice= {slice}")
            print(f"\tIPMWRITE= {IPMWRITE}")
            print(f"\tIPMWSTART= {IPMWSTART}")
            print(f"\tIPMWSTEPS= {IPMWSTEPS}")

    def vpm(
        self, 
        num_particles: int, 
        num_equations: int,
        mode: int,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
        RHS_PM: np.ndarray,
        timestep: int,
        viscosity: float,
    ):
        """_summary_

        Args:
            num_particles (int): Number of particles
            num_equations (int): Number of equations to model
            mode (int): 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse
            particle_positions (np.ndarray): Particle positions array of shape (3, NVR_in)
            particle_strengths (np.ndarray): Particle strenghts array of shape (num_equations + 1, NVR_in)
            particle_velocities (np.ndarray): Particle velocities array of shape (3, NVR_in)
            particle_deformations (np.ndarray): Particle deformation array of shape (3, NVR_in)
            RHS_PM (np.ndarray): Forcing term of $\\nabla^2 u_i = RHS_i$ of shape (num_equations, NXB, NYB, NZB)
            timestep (int): Timestep
            viscosity (float): Viscosity term for the diffusion equation
            max_particle_num (float): Maximum number of particles allowed
        """
        # Ensure numpy arrays are contiguous
        XP_ptr = dp_array_to_pointer(particle_positions   , copy = True)
        QP_ptr = dp_array_to_pointer(particle_strengths   , copy = True)
        UP_ptr = dp_array_to_pointer(particle_velocities  , copy = True)
        GP_ptr = dp_array_to_pointer(particle_deformations, copy = True)
        RHS_pm_ptr = dp_array_to_pointer(RHS_PM)

        # Get the pointers to the numpy arrays
        Velx_ptr = dp_array_to_pointer(self.particle_mesh.Ux)
        Vely_ptr = dp_array_to_pointer(self.particle_mesh.Uy)
        Velz_ptr = dp_array_to_pointer(self.particle_mesh.Uz)
        
        self._lib.vpm(
            XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(num_equations)),
            byref(c_int(mode)), RHS_pm_ptr, Velx_ptr, Vely_ptr, Velz_ptr,
            byref(c_int(timestep)), byref(c_double(viscosity)), 
            byref(c_int(self.max_particle_num))
        )

        # store the results
        self.particles.particle_positions = pointer_to_dp_array(XP_ptr, particle_positions.shape)
        self.particles.particle_strengths = pointer_to_dp_array(QP_ptr, particle_strengths.shape)
        self.particles.particle_velocities = pointer_to_dp_array(UP_ptr, particle_velocities.shape)
        self.particles.particle_deformations = pointer_to_dp_array(GP_ptr, particle_deformations.shape)
        # self._store5 = RHS_pm_ptr
        
    def remesh_particles_3d(self, iflag: int):
        """Remesh particles in 3D

        Args:
            iflag (int): Flag to remesh particles
        """
        NVR = self.particles.NVR
        print_green(f"\tNumber of particles before remeshing: {NVR}", self.rank)
        XP_arr = F_Array((3, NVR))
        QP_arr = F_Array((self.num_equations + 1, NVR))
        UP_arr = F_Array((3, NVR))
        GP_arr = F_Array((3, NVR))

        NVR = c_int(NVR)
        XP_struct = XP_arr.to_ctype()
        QP_struct = QP_arr.to_ctype()
        UP_struct = UP_arr.to_ctype()
        GP_struct = GP_arr.to_ctype()
        self._lib.remesh_particles_3d(
            byref(c_int(iflag)), byref(XP_struct), byref(QP_struct),
            byref(UP_struct), byref(GP_struct), byref(NVR)
        )
        XP_arr = F_Array.from_ctype(XP_struct)
        QP_arr = F_Array.from_ctype(QP_struct)
        UP_arr = F_Array.from_ctype(UP_struct)
        GP_arr = F_Array.from_ctype(GP_struct)
        return XP_arr, QP_arr, UP_arr, GP_arr

    def print_pmesh_parameters(self):
        """
            Print the parameters of the particle mesh
        """
        print_IMPORTANT(f"Pmesh parameters {self.rank}/{self.num_processors - 1}")
        self._lib.print_pmeshpar()

    def print_projlib_parameters(self):
        """
            Print the parameters of the projection library
        """
        print_IMPORTANT(f"Proj lib parameters {self.rank}/{self.num_processors - 1}")
        self._lib.print_projlib()

    def print_pmgrid_parameters(self):
        """
            Print the parameters of the particle mesh grid
        """
        print_IMPORTANT(f"PM grid parameters {self.rank}/{self.num_processors - 1}")
        self._lib.print_pmgrid()

    def print_vpm_vars(self):
        """
            Print the variables of the VPM
        """
        print_IMPORTANT(f"VPM variables {self.rank}/{self.num_processors - 1}")
        self._lib.print_vpm_vars()
    
    def print_vpm_size(self):
        """
            Print the size of the VPM
        """
        print_IMPORTANT(f"VPM size {self.rank}/{self.num_processors - 1}")
        self._lib.print_vpm_size()

    def print_parvar(self):
        """
            Print the particle variables
        """
        print_IMPORTANT(f"Particle variables {self.rank}/{self.num_processors - 1}")
        self._lib.print_parvar()
    
    def set_rhs_pm(self, RHS_PM: np.ndarray):
        """
        Set the right-hand side of the particle mesh.
        """
        RHS_PM = np.ascontiguousarray(RHS_PM, dtype=np.float64)
        # Pass as array of shape (num_equations, NXB, NYB, NZB) not pointer
        RHS_PM_ = RHS_PM.ctypes.data_as(POINTER(c_double))
        sizes = np.array(RHS_PM.shape, dtype=np.int32)
        size1 = sizes[0]
        size2 = sizes[1]
        size3 = sizes[2]
        size4 = sizes[3]
        self._lib.set_RHS_pm(
            RHS_PM_, 
            byref(c_int(size1)), 
            byref(c_int(size2)), 
            byref(c_int(size3)), 
            byref(c_int(size4))
        )
    
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
        self._lib.get_Xbound(byref(xbound))
        return np.array(xbound)

    def get_size_XP(self):
        """
            Get the size of the particle positions
        """
        size_XP = (c_int*2)()
        self._lib.get_size_XP(byref(size_XP))
        return size_XP.value
    
    def get_num_equations(self):
        """
            Get the number of equations
        """
        neqpm = c_int()
        self._lib.get_neqpm(byref(neqpm))
        return neqpm.value
    
    def finalize(self):
        self._lib.finalize()
        if self.rank == 0:
            print_green(f"Finalized VPM. {self.rank}")

    def __del__(self):
        self.finalize()
        os.remove(self.vpm_lib._lib_vpm_path)
