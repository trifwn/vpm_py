from vpm_py import VPM
from time import sleep
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def print_blue(text):
    print(f"\033[94m{text}\033[00m")

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_red(text):
    print(f"\033[91m{text}\033[00m")

def print_IMPORTANT(text):
    print(f"\033[93m{'-'*100}\033[00m")
    print(f"\033[91m{text}\033[00m")
    print(f"\033[93m{'-'*100}\033[00m")

def main():
    # try:
    #     import vpm_f2py
    # except ImportError as e:
    #     print(f"Error: {e}")
    #     # exit(1)
    vpm = VPM()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    ################TESTING################
    if rank == 0:
        lib = vpm._lib
        import ctypes
        import numpy as np
        from ctypes import c_int,  byref, POINTER, cdll, c_double, c_void_p

        # Declare argument and return types (updated)
        lib.print_array.argtypes = []
        lib.init_array.argtypes = [c_int]
        lib.get_module_array.restype = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')

        # Declare argument and return types for NDArray functions
        lib.free_ndarray.argtypes = []
        lib.get_ndarray_data_ptr.argtypes = []
        lib.get_ndarray_data_ptr.restype = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        lib.print_ndarray.argtypes = []
        lib.init_ndarray.argtypes = [POINTER(c_int), POINTER(c_int)] 

        # ARRAY TESTING
        print_IMPORTANT("Testing the Fortran API - ARRAY")

        # Initialize the array in Fortran
        n = 5
        lib.init_array(n)

        # Print the initial array in Fortran
        lib.print_array()

        # Get a NumPy array view of the Fortran module array
        print(f"Getting the array from Fortran")
        fortran_array = lib.get_module_array()
        arr = (ctypes.c_double * n)
        nd_ptr = ctypes.cast(fortran_array, ctypes.POINTER(arr))
        np_array = np.ctypeslib.as_array(nd_ptr.contents, shape=(n,))
        print(f"NumPy array: {np_array}")

        for i in range(n):
            seed = np.random.randint(0, 100)
            np_array[i] = seed
            print(f"Seed: {seed}")
            # Print the modified array in Fortran
            lib.print_array()

        lib.init_array(n)
        lib.change_array()
        print("After reinitialization")
        print(f"Arry in Py: {np_array}")
        lib.print_array()

        # MATRIX TESTING
        print_IMPORTANT("Testing the Fortran API - MATRIX")
        # Create a 2D array in Fortran
        shape = (3, 4, 5)  # Python tuple representing the shape
        shape = np.ascontiguousarray(shape, dtype=np.int32)
        shape_ptr = shape.ctypes.data_as(POINTER(c_int))
        ndims = len(shape)
        lib.init_ndarray(
            shape_ptr, byref(c_int(ndims)) 
        )

        # Print the initial array in Fortran
        lib.print_ndarray()
        # Get the data pointer and convert to NumPy array
        data_ptr = lib.get_ndarray_data_ptr()
        arr = (ctypes.c_double * np.prod(shape))
        nd_ptr = ctypes.cast(data_ptr, ctypes.POINTER(arr))
        np_array = np.ctypeslib.as_array(nd_ptr.contents, shape=shape)
        print(f"NumPy array: {np_array}")

        # # Modify the array in Python
        # for i in range(shape[0]):
        #     for j in range(shape[1]):
        #         seed = np.random.randint(0, 100)
        #         np_array[i, j] = seed
        #         print(f"Seed at ({i}, {j}): {seed}")
        #         lib.print_ndarray()  # Print after each modification

        # # Free the Fortran array
        # lib.free_ndarray()

        #######################################
    exit(1)

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
    NTIME = 0
    NI = 0.1
    NVRM = NVR + 10

    ##################################### MODE TESTING #####################################
    for WhatToDo in [0, 1, 2, 3, 4, 2 ,5]:
    # 0: Initialize
    # 1: Solve
    # 2: Convect
    # 3: Project
    # 4: Back
    # 5: Diffuse
        if rank == 0:
            print_IMPORTANT(f"Calling the vpm subroutine. WhatToDo = {WhatToDo}")
            print_green(f"Number of particles: {NVR}")
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
            timestep= NTIME,
            viscosity= NI,
            max_particle_num= NVRM
        )
        comm.Barrier()
        if rank == 0:
            print_IMPORTANT(f"Completed successfully. WhatToDo = {WhatToDo}")
            print_green(f"Number of particles: {NVR}")
            print('\n\n\n')

        # DEBUGGING
        # if rank == 0:
        #     vpm.print_projlib_parameters()
        #     vpm.print_pmgrid_parameters()
        #     vpm.print_pmesh_parameters()
        #     vpm.print_vpm_vars()
        #     vpm.print_vpm_size()
        #     vpm.print_parvar()
        #     print('\n\n')
        # sleep(1)

        NVR = vpm.NVR
        XP = vpm.XP
        QP = vpm.QP
        UP = vpm.UP
        GP = vpm.GP
        # Convection
        if WhatToDo == 2:
            DT = 0.1
            XP = XP + UP * DT
            QP = QP - 1. *GP * DT
            vpm.XP = XP
            vpm.QP = QP


        # Remesh the particles
        comm.Barrier()
        vpm.remesh_particles_3d(1)
        comm.Barrier()

        NX_pm = vpm.NX_pm
        NY_pm = vpm.NY_pm
        NZ_pm = vpm.NZ_pm
        neqpm = vpm.num_equations
        NN = vpm.nn
        NN_bl = vpm.nn_bl
        Xbound = vpm.xbound
        Dpm = vpm.dpm
        if rank == 0:
            print_IMPORTANT(f"Remeshing completed successfully. WhatToDo = {WhatToDo}")
            print_green(f"Dpm: {Dpm}")
            print_green(f"Size of PM grid: {NX_pm} x {NY_pm} x {NZ_pm}")
            print_green(f"Number of equations: {neqpm}")
            print_green(f"Number of particles: {NVR}")
            print_green(f"Size of RHS_pm: {RHS_pm_in.shape}")
            print_green(f"Size of XP: {XP.shape}")
            print_green(f"Size of QP: {QP.shape}")
            print_green(f"Size of UP: {UP.shape}")
            print_green(f"Size of GP: {GP.shape}")
            print_green(f"Size of NN: {NN.shape}")
            print_green(f"Size of NN_bl: {NN_bl.shape}")
            print_green(f"Size of Xbound: {Xbound.shape}")
        
    comm.Barrier()
    print_green(f"Rank: {rank} completed successfully.")
    del(vpm)

if __name__ == "__main__":
    main()