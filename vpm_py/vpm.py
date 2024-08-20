import numpy as np
import os
import glob

from ctypes import c_int,  byref, POINTER, cdll, c_double
from shutil import copy2
from tempfile import NamedTemporaryFile
from mpi4py import MPI

# Local imports
from . import ParticleMesh
from . import Particles

here = os.path.abspath(os.path.dirname(__file__))
lib_path = glob.glob(os.path.join(here, 'libvpm_*.so'))[0]
lib_ext = lib_path[lib_path.rfind('.'):]

def print_blue(text):
    print(f"\033[94m{text}\033[00m")

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_IMPORTANT(text):
    print(f"\033[93m{'-'*100}\033[00m")
    print(f"\033[91m{text}\033[00m")
    print(f"\033[93m{'-'*100}\033[00m")

def print0(rank, text):
    """Print only if rank is 0"""
    if rank == 0:
        print(text) 

class VPM(object):
    """Interface to the VPM Fortran routines.

    Attributes
   
    """

    def __init__(
        self,
        num_equations: int = 3,
    ):
        super().__init__()
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_ext)
        tmp.close()
        self._lib_path = tmp.name
        self.original_lib = lib_path
        copy2(lib_path, self._lib_path)
        self._lib = cdll.LoadLibrary(self._lib_path)
        
        # Print the name of the so used
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()
        if self.rank == 0:
            print_IMPORTANT(f"Using {lib_path}")

        # Define argument types for the Fortran subroutines
        # API.init
        self._lib.init.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int),
                                    POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_int),
                                    POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                    POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                    POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                    POINTER(c_int)]
        # API.vpm
        self._lib.vpm.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_double),
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int),
            POINTER(c_double), POINTER(c_int)
        ]
        # API.finalize
        self._lib.finalize.argtypes = []
        
        # API.remesh_particles_3d
        self._lib.remesh_particles_3d.argtypes = [POINTER(c_int)]


        # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
        # Find the factors of the number of processors
        factors = np.unique(np.array([i for i in range(1, self.num_procs + 1) if self.num_procs % i == 0]))
        # Find the subsets of factors that multiply to the number of processors
        subsets = []
        for i in range(len(factors)):
            for j in range(i, len(factors)):
                for k in range(j, len(factors)):
                    if factors[i] * factors[j] * factors[k] == self.num_procs:
                        subsets.append((factors[i], factors[j], factors[k]))
        # Find the subset that has the smallest sum
        min_sum = np.inf
        for subset in subsets:
            if sum(subset) < min_sum:
                min_sum = sum(subset)
                NBI, NBJ, NBK = subset

        # Print the number of processors in each direction
        if self.rank == 0:
            print(f"Number of processors: {self.num_procs}")
            print(f"NBI: {NBI}, NBJ: {NBJ}, NBK: {NBK}")
        
        self.dpm = np.zeros(3)
        # Initialize the VPM
        self.initialize(
            DXpm=4.0,
            DYpm=4.0,
            DZpm=4.0,
            NBI=NBI,
            NBJ=NBJ,
            NBK=NBK
        )
        self.particle_mesh = ParticleMesh()
        self.particles = Particles()

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
        is_box_fixed: bool = True, 
        slice: bool = True, 
        IPMWRITE: int = 1, 
        IPMWSTART: int = 1, 
        IPMWSTEPS: int = 1
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
        max_particle_num: int = 1000
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
        particle_positions = np.ascontiguousarray(particle_positions, dtype=np.float64)
        particle_strengths = np.ascontiguousarray(particle_strengths, dtype=np.float64)
        particle_velocities = np.ascontiguousarray(particle_velocities, dtype=np.float64)
        particle_deformations = np.ascontiguousarray(particle_deformations, dtype=np.float64)
        
        RHS_PM = np.ascontiguousarray(RHS_PM, dtype=np.float64)

        # Get the pointers to the numpy arrays
        u_x = self.particle_mesh.Ux
        u_y = self.particle_mesh.Uy
        u_z = self.particle_mesh.Uz
        u_x = np.ascontiguousarray(u_x, dtype=np.float64)
        u_y = np.ascontiguousarray(u_y, dtype=np.float64)
        u_z = np.ascontiguousarray(u_z, dtype=np.float64)
        Velx_ptr = u_x.ctypes.data_as(POINTER(c_double))
        Vely_ptr = u_y.ctypes.data_as(POINTER(c_double))
        Velz_ptr = u_z.ctypes.data_as(POINTER(c_double))

        XP_ptr = particle_positions.ctypes.data_as(POINTER(c_double))
        QP_ptr = particle_strengths.ctypes.data_as(POINTER(c_double))
        UP_ptr = particle_velocities.ctypes.data_as(POINTER(c_double))
        GP_ptr = particle_deformations.ctypes.data_as(POINTER(c_double))
        RHS_pm_ptr = RHS_PM.ctypes.data_as(POINTER(c_double))

        self._lib.vpm(XP_ptr, QP_ptr, UP_ptr, GP_ptr, byref(c_int(num_particles)),
                        byref(c_int(num_equations)), byref(c_int(mode)), RHS_pm_ptr,
                        Velx_ptr, Vely_ptr, Velz_ptr, byref(c_int(timestep)),
                        byref(c_double(viscosity)), byref(c_int(max_particle_num))
                    )
        
        print(f"\t\tVelx_ptr: {Velx_ptr}")
        print(f"\t\tVely_ptr: {Vely_ptr}")
        print(f"\t\tVelz_ptr: {Velz_ptr}")
        if mode in [1,2]:
            # Get the values from the ptrs to the numpy arrays
            shape_pm = (self.NX_pm, self.NY_pm, self.NZ_pm)
            Ux = np.ctypeslib.as_array(Velx_ptr, shape=shape_pm)
            Uy = np.ctypeslib.as_array(Vely_ptr, shape=shape_pm)
            Uz = np.ctypeslib.as_array(Velz_ptr, shape=shape_pm)
            # Convert to numpy arrays as they are fortran arrays
            Ux = np.array(Ux, dtype=np.float64).T
            Uy = np.array(Uy, dtype=np.float64).T
            Uz = np.array(Uz, dtype=np.float64).T

            self.particle_mesh.update_U(Ux, Uy, Uz)

            if self.rank == 0:
                print_IMPORTANT(f"Getting the values of Ux, Uy, Uz\nGot shape: {shape_pm}")
                        # Check arrays for nan and print the number of nans and their positions
                if np.isnan(Ux).any():
                    print(f"Ux has {np.sum(np.isnan(Ux))} nans")
                    print(np.argwhere(np.isnan(Ux)))
                
                if np.isnan(Uy).any():
                    print(f"Uy has {np.sum(np.isnan(Uy))} nans")
                    print(np.argwhere(np.isnan(Uy)))

                if np.isnan(Uz).any():
                    print(f"Uz has {np.sum(np.isnan(Uz))} nans")
                    print(np.argwhere(np.isnan(Uz)))

    def remesh_particles_3d(self, iflag: int):
        """Remesh particles in 3D

        Args:
            iflag (int): Flag to remesh particles
        """
        self._lib.remesh_particles_3d(byref(c_int(iflag)))
        # print_green(f"Remeshed particles {self.rank}/{self.num_procs - 1} completed successfully.")

    def print_pmesh_parameters(self):
        """
            Print the parameters of the particle mesh
        """
        print_IMPORTANT(f"Pmesh parameters {self.rank}/{self.num_procs - 1}")
        self._lib.print_pmeshpar()

    def print_projlib_parameters(self):
        """
            Print the parameters of the projection library
        """
        print_IMPORTANT(f"Proj lib parameters {self.rank}/{self.num_procs - 1}")
        self._lib.print_projlib()

    def print_pmgrid_parameters(self):
        """
            Print the parameters of the particle mesh grid
        """
        print_IMPORTANT(f"PM grid parameters {self.rank}/{self.num_procs - 1}")
        self._lib.print_pmgrid()

    def print_vpm_vars(self):
        """
            Print the variables of the VPM
        """
        print_IMPORTANT(f"VPM variables {self.rank}/{self.num_procs - 1}")
        self._lib.print_vpm_vars()
    
    def print_vpm_size(self):
        """
            Print the size of the VPM
        """
        print_IMPORTANT(f"VPM size {self.rank}/{self.num_procs - 1}")
        self._lib.print_vpm_size()

    def print_parvar(self):
        """
            Print the particle variables
        """
        print_IMPORTANT(f"Particle variables {self.rank}/{self.num_procs - 1}")
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

    def get_particle_positions(self):
        """
            Get the particle positions
        """
        XP = (c_double * 3 * self.NVR)()
        self._lib.get_particle_positions(byref(XP))
        return np.array(XP).T
    
    @property
    def XP(self):
        return self.get_particle_positions()
    
    @XP.setter
    def XP(self, XP):
        XP = np.ascontiguousarray(XP, dtype=np.float64)
        XP_ptr = XP.ctypes.data_as(POINTER(c_double))
        self._lib.set_particle_positions(XP_ptr)
    
    ##  PARTICLE STRENGTHS
    def get_particle_strengths(self):
        """
            Get the particle strengths
        """
        neq = c_int()
        self._lib.get_neqpm(byref(neq))

        QP = (c_double * (neq.value + 1)  * self.NVR)()
        self._lib.get_particle_strengths(byref(QP))
        return np.array(QP).T

    @property
    def QP(self):
        return self.get_particle_strengths()
    
    @QP.setter
    def QP(self, QP):
        QP = np.ascontiguousarray(QP, dtype=np.float64)
        QP_ptr = QP.ctypes.data_as(POINTER(c_double))
        self._lib.set_particle_strengths(QP_ptr)
    
    ## PARTICLE VELOCITIES
    def get_particle_velocities(self):
        """
            Get the particle velocities
        """
        UP = (c_double * 3 * self.NVR)()
        self._lib.get_particle_velocities(byref(UP))
        return np.array(UP).T
    
    @property
    def UP(self):
        return self.get_particle_velocities()
    
    @UP.setter
    def UP(self, UP):
        UP = np.ascontiguousarray(UP, dtype=np.float64)
        UP_ptr = UP.ctypes.data_as(POINTER(c_double))
        self._lib.set_particle_velocities(UP_ptr)
    
    ## PARTICLE DEFORMATIONS
    def get_particle_deformations(self):
        """
            Get the particle deformations
        """
        GP = (c_double * 3 * self.NVR)()
        self._lib.get_particle_deformation(byref(GP))
        return np.array(GP).T
    
    @property
    def GP(self):
        return self.get_particle_deformations()
    
    @GP.setter
    def GP(self, GP):
        GP = np.ascontiguousarray(GP, dtype=np.float64)
        GP_ptr = GP.ctypes.data_as(POINTER(c_double))
        self._lib.set_particle_deformation(GP_ptr)

    def get_NVR(self):
        """
            Get the number of particles
        """
        NVR = c_int()
        self._lib.get_num_particles(byref(NVR))
        return NVR.value

    @property
    def NVR(self):
        return self.get_NVR()

    def get_size_XP(self):
        """
            Get the size of the particle positions
        """
        size_XP = (c_int*2)()
        self._lib.get_size_XP(byref(size_XP))
        return size_XP.value
    
    @property
    def num_particles(self):
        return self.get_NVR()

    def get_num_equations(self):
        """
            Get the number of equations
        """
        neqpm = c_int()
        self._lib.get_neqpm(byref(neqpm))
        return neqpm.value
    
    @property
    def num_equations(self):
        return self.get_num_equations()

    def get_NN(self):
        """
            Get the number of cells in each direction
        """
        NN = (c_int * 3)() 
        self._lib.get_NN(byref(NN))
        return np.array(NN)
    
    def get_NN_bl(self):
        """
            Get the number of blocks in each direction
        """
        # NN is an array of 6 integers
        NN_bl = (c_int * 6)()
        self._lib.get_NN_bl(byref(NN_bl))
        return np.array(NN_bl) 

    def get_xbound(self):
        """
            Get the xbound
        """
        xbound = (c_double * 6)()
        self._lib.get_Xbound(byref(xbound))
        return np.array(xbound)
    
    @property
    def nn(self):
        return self.get_NN()
    
    @property
    def nn_bl(self):
        return self.get_NN_bl()
    
    @property
    def xbound(self):
        return self.get_xbound()
 
    def finalize(self):
        self._lib.finalize()
        if self.rank == 0:
            print_green(f"Finalized VPM. {self.rank}")

    def __del__(self):
        self.finalize()
        os.remove(self._lib_path)
