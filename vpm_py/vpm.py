import numpy as np
from mpi4py import MPI
import os

from ctypes import c_int,  byref, c_double, POINTER, c_char, c_char_p, create_string_buffer, cast
# Local imports
from . import ParticleMesh
from . import Particles
from .console_io import  print_IMPORTANT, print_green, print_blue, print_red
from .vpm_lib import VPM_Lib
from .utils import divide_processors
from .arrays import F_Array, F_Array_Struct
from .vpm_dtypes import dp_array_to_pointer, pointer_to_dp_array

class VPM(object):
    """
    Interface to the VPM Fortran routines.
    """

    def __init__(
        self,
        number_of_equations: int = 3,
        number_of_processors: int = 1,
        rank: int = 0,
        verbocity: int = 1,
        dx_particle_mesh: float = 0.2,
        dy_particle_mesh: float = 0.2,
        dz_particle_mesh: float = 0.2,
    ):
        lib = VPM_Lib()
        self._lib = lib._lib_vpm
        self.vpm_lib = lib

        self.rank = rank
        self.num_processors = number_of_processors
        
        # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
        NBI, NBJ, NBK = divide_processors(number_of_processors)
        if self.rank == 0:
            print(f"Number of processors: {self.num_processors}")
            print(f"NBI: {NBI}, NBJ: {NBJ}, NBK: {NBK}")
        
        self.dpm = np.array([dx_particle_mesh, dy_particle_mesh, dz_particle_mesh])
        self.NBI = NBI
        self.NBJ = NBJ
        self.NBK = NBK
        self.num_equations = number_of_equations

        self.set_verbosity(verbocity)
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
        remesh: bool = True, 
        use_tree: bool = True, 
        ilevmax: int = 4,
        OMPTHREADS: int = 1,
        is_domain_fixed: bool = False, 
        IPMWRITE: int = 0, 
        IPMWSTART: int = -1, 
        IPMWSTEPS: int = -1,
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
            remesh (bool, optional): Remesh. Defaults to True.
            use_tree (bool, optional): Use tree. Defaults to True.
            ilevmax (int, optional): Maximum level. Defaults to 4.
            OMPTHREADS (int, optional): Number of threads. Defaults to 1.
            is_domain_fixed (bool, optional): Is box fixed. Defaults to True.
            IPMWRITE (int, optional): Write IPM. Defaults to 1.
            IPMWSTART (int, optional): Write IPM start. Defaults to 1.
            IPMWSTEPS (int, optional): Write IPM steps. Defaults to 1.
        """
        self.dpm = np.array([DXpm, DYpm, DZpm])
        self._lib.init(
            byref(c_double(DXpm)), byref(c_double(DYpm)), byref(c_double(DZpm)), byref(c_int(projection_type)),
            byref(c_int(boundary_condition)), byref(c_int(variable_volume)), byref(c_int(ncoarse)),
            byref(c_int(NBI)), byref(c_int(NBJ)), byref(c_int(NBK)), byref(c_int(remesh)), 
            byref(c_int(use_tree)), byref(c_int(ilevmax)), byref(c_int(OMPTHREADS)), byref(c_int(is_domain_fixed)),
            byref(c_int(IPMWRITE)), byref(c_int(IPMWSTART)), byref(c_int(IPMWSTEPS)),
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
            print(f"\tncoarse= {ncoarse}")
            print(f"\tNBI= {NBI}")
            print(f"\tNBJ= {NBJ}")
            print(f"\tNBK= {NBK}")
            print(f"\tNREMESH= {remesh}")
            print(f"\tiyntree= {use_tree}")
            print(f"\tilevmax= {ilevmax}")
            print(f"\tOMPTHREADS= {OMPTHREADS}")
            print(f"\tidefine= {is_domain_fixed}")
            print(f"\tIPMWRITE= {IPMWRITE}")
            print(f"\tIPMWSTART= {IPMWSTART}")
            print(f"\tIPMWSTEPS= {IPMWSTEPS}")
    
    def vpm(
        self, 
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        mode: int,
        timestep: int,
        viscosity: float,
        num_particles: int | None = None,
    ):
        """_summary_

        Args:
            num_equations (int): Number of equations to model
            particle_positions (np.ndarray): Particle positions array of shape (3, NVR_in)
            particle_charges (np.ndarray): Particle charges array of shape (num_equations + 1, NVR_in)
            mode (int): 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse
            timestep (int): Timestep
            viscosity (float): Viscosity term for the diffusion equation
            num_particles (int, optional): Number of points to treat as particles. When None defaults
                        to all the points in the arrays. If different from None, the rest of the points
                        will be treated as source terms. Defaults to None.
        """
        # Check if the arrays are F_Arrays
        if isinstance(particle_positions, F_Array):
            particle_positions = particle_positions.data
        if isinstance(particle_charges, F_Array):
            particle_charges = particle_charges.data

        # Check that the arrays have the same number of particles
        NVR_size = particle_positions.shape[1]
        if not (particle_charges.shape[1] == NVR_size):
            raise ValueError("Number of particles in particle_charges does not match particle_positions")
        
        # Create F_Arrays to store the results
        particle_velocities = np.zeros_like(particle_positions, dtype=np.float64, order='K')  
        particle_deformations = np.zeros_like(particle_positions, dtype=np.float64, order='K')

        # Get the pointers to arrays for the particles
        XP_ptr = dp_array_to_pointer(particle_positions, copy = True)
        UP_ptr = dp_array_to_pointer(particle_velocities, copy = True)
        QP_ptr = dp_array_to_pointer(particle_charges, copy = True)
        GP_ptr = dp_array_to_pointer(particle_deformations, copy = True)

        # Get the pointers to arrays for the grid values
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size= 1)
        Velx_ptr = F_Array_Struct.null(ndims=3, total_size= 1)
        Vely_ptr = F_Array_Struct.null(ndims=3, total_size= 1)
        Velz_ptr = F_Array_Struct.null(ndims=3, total_size= 1)

        if num_particles is None:
            num_particles = NVR_size

        self._lib.vpm(
            XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(num_equations)),byref(c_int(mode)), 
            byref(RHS_pm_ptr), byref(Velx_ptr), byref(Vely_ptr), byref(Velz_ptr),
            byref(c_int(timestep)), byref(c_double(viscosity)),byref(c_int(NVR_size))
        )
        # store the results of the particles
        self.particles.particle_positions = pointer_to_dp_array(XP_ptr, particle_positions.shape)
        self.particles.particle_charges = pointer_to_dp_array(QP_ptr, particle_charges.shape)
        self.particles.particle_velocities = pointer_to_dp_array(UP_ptr, particle_velocities.shape)
        self.particles.particle_deformations = pointer_to_dp_array(GP_ptr, particle_deformations.shape)
                
        # store the results of the particle mesh
        neq = self.num_equations

        self.particle_mesh.number_equations = neq
        if not Velx_ptr.is_null():
            Ux_arr = F_Array.from_ctype(Velx_ptr, ownership=True, name = "Ux")
            self.particle_mesh.Ux = Ux_arr.transfer_data_ownership()
        
        if not Vely_ptr.is_null():
            Uy_arr = F_Array.from_ctype(Vely_ptr, ownership=True, name = "Uy")
            self.particle_mesh.Uy = Uy_arr.transfer_data_ownership()
        
        if not Velz_ptr.is_null():
            Uz_arr = F_Array.from_ctype(Velz_ptr, ownership=True, name = "Uz")
            self.particle_mesh.Uz = Uz_arr.transfer_data_ownership()

        if not RHS_pm_ptr.is_null():
            RHS_arr = F_Array.from_ctype(RHS_pm_ptr, ownership=True, name = "RHS")
            self.particle_mesh.RHS = RHS_arr.transfer_data_ownership()

    def remesh_particles(self, project_particles: bool, particles_per_cell: int= 1, cut_off: float = 1e-9):
        """Remesh particles in 3D

        Args:
            project_particles (bool): Whether to project the particles or use RHS_pm 
            particles_per_cell (int): Number of particles per cell
            cut_off (float): Cut off value for the remeshing
        """
        NVR = self.particles.NVR
        neq = self.num_equations
        XP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        UP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        QP_struct = F_Array_Struct.null(ndims=2, total_size=(neq + 1)*NVR)
        GP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        NVR = c_int(NVR)
        self._lib.remesh_particles_3d(
            byref(c_int(project_particles)), byref(c_int(particles_per_cell)),
            byref(XP_struct), byref(QP_struct),
            byref(UP_struct), byref(GP_struct),
            byref(NVR), byref(c_double(cut_off))
        )
        XP_arr = F_Array.from_ctype(XP_struct, ownership=True, name = "XP_remesh")
        QP_arr = F_Array.from_ctype(QP_struct, ownership=True, name = "QP_remesh")
        UP_arr = F_Array.from_ctype(UP_struct, ownership=True, name = "UP_remesh")
        GP_arr = F_Array.from_ctype(GP_struct, ownership=True, name = "GP_remesh")

        # store the results
        self.particles.particle_positions = XP_arr.transfer_data_ownership()
        self.particles.particle_velocities = UP_arr.transfer_data_ownership()
        self.particles.particle_charges = QP_arr.transfer_data_ownership()
        self.particles.particle_deformations = GP_arr.transfer_data_ownership()
        return XP_arr, QP_arr

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
    
    def set_verbosity(self, verbosity: int):
        """
            Set the verbosity of the VPM
        """
        self._lib.set_verbosity(byref(c_int(verbosity)))

    def finalize(self):
        self._lib.finalize()
        if self.rank == 0:
            print_green(f"Finalized VPM. {self.rank}")

    def __del__(self):
        self.finalize()
        os.remove(self.vpm_lib._lib_vpm_path)
