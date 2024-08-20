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
    vpm = VPM()
    print_IMPORTANT(
        vpm.original_lib
    )
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    np_procs = comm.Get_size()

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    if rank == 0:
        print_blue(f"Number of processes: {np_procs}")
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()
     
    DT = 0.1
    NI = 0.1
    neq = 3 

    Vref = np.float64(1.) # Convert Vref to 64-bit floating point

    # Create particles
    NVR = np.int32(100)
    XPR = np.zeros((3, NVR), dtype=np.float64)
    QPR = np.ones((neq + 1, NVR), dtype=np.float64)
    MAX_PARTICLES = 10000
    
    # Initialization VPM
    from hill_spherical_vortex import hill_assign
    _, RHS_pm_in = hill_assign(
        Dpm= vpm.dpm,
        NN= vpm.nn,
        NN_bl= vpm.nn_bl,
        Xbound= vpm.xbound,
        neqpm= vpm.num_equations,
    )
    
    comm.Barrier()
    vpm.vpm(
        num_particles=NVR,
        num_equations=neq,
        mode = 0,
        particle_positions= XPR, 
        particle_strengths=QPR, 
        particle_velocities=np.zeros((3, NVR)), 
        particle_deformations=np.zeros((3, NVR)), 
        RHS_PM=RHS_pm_in,
        timestep=0,
        viscosity=NI,
        max_particle_num=MAX_PARTICLES
    )
    comm.Barrier()

    # Remeshing
    if rank == 0:
        st = MPI.Wtime()

    _, RHS_pm_hill = hill_assign(
        Dpm= vpm.dpm,
        NN= vpm.nn,
        NN_bl= vpm.nn_bl,
        Xbound= vpm.xbound,
        neqpm= vpm.num_equations,
    )
    vpm.set_rhs_pm(RHS_pm_hill)
    if rank == 0:
        st = MPI.Wtime()
        print_red(f"Remeshing")
    vpm.remesh_particles_3d(-1) 
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing took {int((et - st) / 60)}m {int(et - st) % 60}s")

    if rank == 0:
        XPR = vpm.XP
        QPR = vpm.QP
        UPR = np.zeros_like(XPR)
        GPR = np.zeros_like(XPR)
        size_XPR = XPR.shape
    else:
        size_XPR = None
    
    size_XPR = comm.bcast(size_XPR, root=0)
    if rank != 0:
        XPR = np.zeros(size_XPR, dtype=np.float64)
        QPR = np.zeros((neq + 1, size_XPR[1]), dtype=np.float64)
        UPR = np.zeros_like(XPR)
        GPR = np.zeros_like(XPR)
    
    # Create the plot to live update the particles
    if rank == 0:
        fig = plt.figure()
        ax: Axes3D = fig.add_subplot(111, projection='3d')
        sc = ax.scatter(XPR[0,:], XPR[1,:], XPR[2,:], c=QPR[0,:], marker='o')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.ion()
        plt.show(block = False)

    comm.Barrier()
    # call vpm with mode = 0
    vpm.vpm(
        num_particles= vpm.num_particles,
        num_equations= vpm.num_equations,
        mode = 0,
        particle_positions= XPR,
        particle_strengths=QPR,
        particle_velocities=UPR,
        particle_deformations=GPR,
        RHS_PM=RHS_pm_in,
        timestep=0,
        viscosity=NI,
        max_particle_num=MAX_PARTICLES
    )

    # Main loop
    NVR = vpm.NVR
    T = 0
    max_iter = 20
    for i in range(1, max_iter):
        comm.Barrier()
        T += DT
        if rank == 0:
            print(f"---------------------------------")
            print_green(f"ITERATION= {i} of {max_iter}")
            print_green(f"T={DT*i}")
            print_green(f"DT={DT}")
            print(f"---------------------------------")

        UPR = np.zeros((3, NVR), dtype=np.float64)
        GPR = np.zeros((3, NVR), dtype=np.float64)

        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 2,
            particle_positions= XPR, 
            particle_strengths=QPR, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            timestep=i,
            viscosity=NI,
            max_particle_num=MAX_PARTICLES
        )

        if rank == 0:
            # Print the size of the particles
            print("\n\n")
            print(f"{'-'*100}")
            print(f"Number of particles: {NVR}")
            print(f"Number of equations: {neq}")
            print('\n')

            print(f"XPR:")
            print(f"Mean: {np.mean(XPR, axis=1)}")
            print(f"Max: {np.max(XPR, axis=1)}")
            print(f"Min: {np.min(XPR, axis=1)}")
            print('\n')

            print(f"UPR:")
            print(f"Mean: {np.mean(UPR, axis=1)}")
            print(f"Max: {np.max(UPR, axis=1)}")
            print(f"Min: {np.min(UPR, axis=1)}")
            print('\n')

            print(f"GPR:")
            print(f"Mean: {np.mean(GPR, axis=1)}")
            print(f"Max: {np.max(GPR, axis=1)}")
            print(f"Min: {np.min(GPR, axis=1)}")
            print('\n')

            # print(f"RHS:")
            # RHS = vpm.particle_mesh.RHS.flatten()
            # print(f"Mean: {np.mean(RHS)}")
            # print(f"Max: {np.max(RHS)}")
            # print(f"Min: {np.min(RHS)}")

            # print(f"U_PM:")
            # ux = vpm.particle_mesh.Ux
            # uy = vpm.particle_mesh.Uy
            # uz = vpm.particle_mesh.Uz
            # V_mag = np.sqrt(ux**2 + uy**2 + uz**2).flatten()
            # print(f"Mean: {np.mean(V_mag)}")
            # print(f"Max: {np.max(V_mag)}")
            # print(f"Min: {np.min(V_mag)}")
            
            print(f"{'-'*100}")
            print('\n\n')
            
            # Update the plot
            sc.set_offsets(np.c_[XPR[0,:], XPR[1,:]])
            sc.set_3d_properties(XPR[2,:], 'z')
            Q_MAG = np.sqrt(QPR[0,:]**2 + QPR[1,:]**2 + QPR[2,:]**2)
            sc.set_array(Q_MAG)
            fig.suptitle(f"{NVR} Particles at T={T}")
            
            # Update the particles
            fig.canvas.flush_events()
            fig.canvas.draw()
            plt.pause(0.0001)
        sleep(0.5)
        comm.Barrier()

        # Move the particles
        if rank == 0:
            XPR = vpm.XP
            UPR = vpm.UP
            QPR = vpm.QP
            GPR = vpm.GP
            for j in range(vpm.NVR):
                XPR[:, j] += (UPR[:, j]) * DT
                FACDEF = 1.0
                QPR[:3, j] -= FACDEF * GPR[:3, j] * DT
            vpm.XP = XPR
            vpm.UP = UPR
            vpm.QP = QPR
            vpm.GP = GPR
            print_IMPORTANT(f"Convected Particles")

        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 0,
            particle_positions= XPR, 
            particle_strengths=QPR, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            timestep=i,
            viscosity=NI,
            max_particle_num=MAX_PARTICLES
        )
        comm.Barrier()

        # Remeshing
        # if rank == 0:
        #     st = MPI.Wtime()
        #     print_red(f"Remeshing")
        # vpm.remesh_particles_3d(1)
        # if rank == 0:
        #     et = MPI.Wtime()
        #     print(f"\tRemeshing took {int((et - st) / 60)}m {int(et - st) % 60}s")

    MPI.Finalize()

if __name__ == "__main__":
    main()
