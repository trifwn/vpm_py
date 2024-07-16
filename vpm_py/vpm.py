import numpy as np
import os
import glob

from ctypes import c_bool, c_int, c_float, byref, POINTER, cdll, c_double
from shutil import copy2
from tempfile import NamedTemporaryFile
from mpi4py import MPI

here = os.path.abspath(os.path.dirname(__file__))
lib_path = glob.glob(os.path.join(here, 'libvpm*.so'))[0]
lib_ext = lib_path[lib_path.rfind('.'):]

def print_blue(text):
    print(f"\033[94m{text}\033[00m")

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_IMPORTANT(text):
    print(f"\033[93m{"-"*100}\033[00m")
    print(f"\033[91m{text}\033[00m")
    print(f"\033[93m{"-"*100}\033[00m")

class VPM(object):
    """Interface to the VPM Fortran routines.

    Attributes
   
    """

    def __init__(self):
        super().__init__()
        tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_ext)
        tmp.close()
        self._lib_path = tmp.name
        copy2(lib_path, self._lib_path)
        self._lib = cdll.LoadLibrary(self._lib_path)

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
        self.finalize = self._lib.finalize
        
        # API.remesh_particles_3d
        self._lib.remesh_particles_3d.argtypes = [POINTER(c_int)]

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
        # For now, assume NBI = NBJ = NBK
        NBI = NBJ = NBK = 1
        while NBI * NBJ * NBK < self.num_procs:
            if NBI == NBJ == NBK:
                NBI *= 2
            elif NBI == NBJ:
                NBK *= 2
            else:
                NBJ *= 2
        
        # Print the number of processors in each direction
        if self.rank == 0:
            print(f"Number of processors: {self.num_procs}")
            print(f"NBI: {NBI}, NBJ: {NBJ}, NBK: {NBK}")

        # Initialize the VPM
        self.initialize(
            DXpm=4.0,
            DYpm=4.0,
            DZpm=4.0,
            NBI=NBI,
            NBJ=NBJ,
            NBK=NBK
        )

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
        U_PM: np.ndarray,
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
            U_PM (np.ndarray): Velocity field in the particle mesh of shape (3, NXB, NYB, NZB)
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

        u_x = U_PM[0, :, :, :]
        u_y = U_PM[1, :, :, :]
        u_z = U_PM[2, :, :, :]
        u_x = np.ascontiguousarray(u_x, dtype=np.float64)
        u_y = np.ascontiguousarray(u_y, dtype=np.float64)
        u_z = np.ascontiguousarray(u_z, dtype=np.float64)

        XP_ptr = particle_positions.ctypes.data_as(POINTER(c_double))
        QP_ptr = particle_strengths.ctypes.data_as(POINTER(c_double))
        UP_ptr = particle_velocities.ctypes.data_as(POINTER(c_double))
        GP_ptr = particle_deformations.ctypes.data_as(POINTER(c_double))
        RHS_pm_ptr = RHS_PM.ctypes.data_as(POINTER(c_double))
        Velx_ptr = u_x.ctypes.data_as(POINTER(c_double))
        Vely_ptr = u_y.ctypes.data_as(POINTER(c_double))
        Velz_ptr = u_z.ctypes.data_as(POINTER(c_double))

        self._lib.vpm(XP_ptr, QP_ptr, UP_ptr, GP_ptr, byref(c_int(num_particles)),
                                  byref(c_int(num_equations)), byref(c_int(mode)), RHS_pm_ptr,
                                  Velx_ptr, Vely_ptr, Velz_ptr, byref(c_int(timestep)),
                                  byref(c_double(viscosity)), byref(c_int(max_particle_num)))
        
        print_green(f"WhatToDo: {mode} from {self.rank + 1}/{self.num_procs} completed successfully.")

    def remesh_particles_3d(self, ncell_rem: int):
        """Remesh particles in 3D

        Args:
            ncell_rem (int): Number of cells to remesh
        """
        self._lib.remesh_particles_3d(byref(c_int(ncell_rem)))

        print_green(f"Remeshed particles {self.rank + 1}/{self.num_procs} completed successfully.")
 
    def finalize(self):
        self._lib.finalize()
        if self.rank == 0:
            print_green(f"Finalized VPM. {self.rank}")

    def __del__(self):
        self.finalize()
        os.remove(self._lib_path)

# Example usage:
if __name__ == "__main__":
    # try:
    #     import vpm_f2py
    # except ImportError as e:
    #     print(f"Error: {e}")
    #     # exit(1)
    vpm = VPM()


    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    if rank == 0:
        print_blue(f"Number of processes: {comm.Get_size()}")
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()
    
    
    # Define the input parameters
    NVR = 200
    neqpm = 4
    NXB = 10
    NYB = 10
    NZB = 10

    # Distribute the particles uniformly in the domain
    divisor = int((np.power(NVR, 1/3)))
    x = np.linspace(0, 1, divisor)
    y = np.linspace(0, 1, divisor)
    z = np.linspace(0, 1, divisor)
    NVR = len(x) * len(y) * len(z)

    XP_in = np.array(np.meshgrid(x, y, z)).reshape(3, -1)
    
    # Define the particle quantities as random numbers
    par_strenghts = np.ones((neqpm, NVR))

    par_velocities = np.zeros((3, NVR), dtype=np.float64)
    par_velocities[0,:] = 1.0

    # Deformation of the particles
    par_deformations = np.zeros((3, NVR), dtype=np.float64)
    
    RHS_pm_in = np.zeros((neqpm, NXB, NYB, NZB), dtype=np.float64)
    Velx = np.zeros((NXB, NYB, NZB), dtype=np.float64)
    Vely = np.zeros((NXB, NYB, NZB), dtype=np.float64)
    Velz = np.zeros((NXB, NYB, NZB), dtype=np.float64)
    NTIME = 0
    NI = 0.1
    NVRM = NVR + 10

    ##################################### MODE TESTING #####################################
    for WhatToDo in [0, 1, 2, 3, 5, 4]:
    # 0: Initialize
    # 1: Solve
    # 2: Convect
    # 3: Project
    # 4: Back
    # 5: Diffuse
        if rank == 0:
            print_IMPORTANT(f"Calling the vpm subroutine. WhatToDo = {WhatToDo}")
            print_green(f"Number of particles: {NVR}")
            print_green(f"Size of Velx: {Velx.shape}")
            print_green(f"Size of Vely: {Vely.shape}")
            print_green(f"Size of Velz: {Velz.shape}")
            print_green(f"Size of RHS_pm: {RHS_pm_in.shape}")
        vpm.vpm(
            num_particles= NVR,
            num_equations= neqpm,
            mode= WhatToDo,
            particle_positions= XP_in,
            particle_strengths= par_strenghts,
            particle_velocities= par_velocities,
            particle_deformations= par_deformations,
            RHS_PM= RHS_pm_in,
            U_PM= np.array([Velx, Vely, Velz]),
            timestep= NTIME,
            viscosity= NI,
            max_particle_num= NVRM
        )
        comm.Barrier()
        if rank == 0:
            print_IMPORTANT(f"Completed successfully. WhatToDo = {WhatToDo}")
            print_green(f"Number of particles: {NVR}")
            print_green(f"Size of Velx: {Velx.shape}")
            print_green(f"Size of Vely: {Vely.shape}")
            print_green(f"Size of Velz: {Velz.shape}")
            print_green(f"Size of RHS_pm: {RHS_pm_in.shape}")
            print('\n\n\n')

        # Remesh the particles
        comm.Barrier()
        vpm.remesh_particles_3d(1)
        comm.Barrier()

    if rank == 0:
        # Plot the velocity fields
        import matplotlib.pyplot as plt
        # Make 3 3D plots for the velocity fields
        fig , axs = plt.subplots(1, 3, figsize=(15, 5))
        axs[0].imshow(Velx[:, :, 0])
        axs[0].set_title("Velx")
        axs[1].imshow(Vely[:, :, 0])
        axs[1].set_title("Vely")
        axs[2].imshow(Velz[:, :, 0])
        axs[2].set_title("Velz")
        plt.show()

        # Plot the particle positions, and color them by the particle quantities and quiver plot for the particle velocities
        fig, axs = plt.subplots(1, 4, figsize=(10, 5), subplot_kw={'projection': '3d'})
        for i in range(3):
            axs[i].scatter(XP_in[0], XP_in[1], XP_in[2], c=par_strenghts[i])
            axs[i].set_title(f"Strength in eq_{["x", "y", "z"][i]}")
        axs[3].quiver(XP_in[0], XP_in[1], XP_in[2], par_velocities[0], par_velocities[1], par_velocities[2])
        axs[3].set_title("Particle velocities")
        plt.show()
    comm.Barrier()
    print_green(f"Rank: {rank} completed successfully.")
    del(vpm)