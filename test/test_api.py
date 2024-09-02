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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    number_of_processors = comm.Get_size()

    vpm = VPM(
        number_of_equations= 3,
        max_particle_num= 1000,
        rank= rank,
        number_of_processors= number_of_processors,
    )
    
    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    if rank == 0:
        print_blue(f"Number of processes: {number_of_processors}")
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
    divisor = 5# int((np.power(NVR, 1/3)))
    x = np.linspace(0, 1, divisor)
    y = np.linspace(0, 1, divisor)
    z = np.linspace(0, 1, divisor)
    NVR = len(x) * len(y) * len(z)

    XP_in = np.array(np.meshgrid(x, y, z)).reshape(3, -1)
    # Define the particle quantities as random numbers
    par_strenghts = np.ones((neqpm, NVR))

    par_velocities = np.zeros((3, NVR), dtype=np.float64)
    par_velocities[0,:] = -1.0
    par_velocities[1,:] = -1.0/ 2.
    par_velocities[2,:] = -1.0/ 3.

    # Deformation of the particles
    par_deformations = np.zeros((3, NVR), dtype=np.float64)
    
    RHS_pm_in = np.zeros((neqpm, NXB, NYB, NZB), dtype=np.float64)
    NTIME = 0
    NI = 0.1
    DT = 0.1

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
        vpm.vpm(
            num_particles= NVR,
            num_equations= neqpm,
            mode= WhatToDo,
            particle_positions= XP_in,
            particle_strengths= par_strenghts,
            particle_velocities= par_velocities,
            particle_deformations= par_deformations,
            timestep= NTIME,
            viscosity= NI,
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

        XP = vpm.particles.XP
        QP = vpm.particles.QP
        UP = vpm.particles.UP
        GP = vpm.particles.GP
        NVR = vpm.particles.NVR

        # Convection
        if WhatToDo == 2:
            for i in range(NVR):
                XP[:, i] = XP[:, i] + UP[:, i] * DT
                QP[:3, i] = QP[:3, i] - 1. * GP[:3, i] * DT
            
        # Remesh the particles
        comm.Barrier()
        (
            XP_new, QP_new, UP_new, GP_new
        ) = vpm.remesh_particles_3d(1)
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
    MPI.Finalize()