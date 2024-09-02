from vpm_py import VPM
from time import sleep
import numpy as np
from mpi4py import MPI
import numpy as np

from vpm_py.vpm_io import print_IMPORTANT, print_red, print_green, print_blue
from vpm_py.visualization import Particle3DPlot
from hill_spherical_vortex import hill_assign

def main():
    # Initialize MPI
    comm = MPI.COMM_WORLD
    start_time = MPI.Wtime()
    rank = comm.Get_rank()
    np_procs = comm.Get_size()
    
    # Initialize VPM
    vpm = VPM(
        max_particle_num= 10000,
        number_of_equations= 3,
        number_of_processors= np_procs,
        rank= rank,
        dx_particle_mesh= 0.1,
        dy_particle_mesh= 0.1,
        dz_particle_mesh= 0.1
    )

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    print_blue(f"Number of processes: {np_procs}", rank)
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()
     
    DT = 1 * 0.1
    NI = -0.1
    neq = 3 
    UINF = np.array([0., 0., 0.])

    # Create particles
    NVR = 100
    XPR = np.zeros((3, NVR), dtype=np.float64)
    QPR = np.ones((neq + 1, NVR), dtype=np.float64)
    UPR = np.zeros((3, NVR), dtype=np.float64)
    GPR = np.zeros((3, NVR), dtype=np.float64)
    
    # Initialization VPM
    comm.Barrier()
    vpm.vpm(
        num_particles=NVR,
        num_equations=neq,
        mode = 0,
        particle_positions= XPR, 
        particle_strengths= QPR, 
        particle_velocities= UPR, 
        particle_deformations= GPR,
        timestep=0,
        viscosity=NI,
    )
    comm.Barrier()

    # Remeshing
    if rank == 0:
        st = MPI.Wtime()

    print_IMPORTANT(f"Hill vortex initialization", rank)
    _, RHS_pm_hill = hill_assign(
        Dpm= vpm.dpm,
        NN= vpm.nn,
        NN_bl= vpm.nn_bl,
        Xbound= vpm.xbound,
        neqpm= vpm.num_equations,
        a = 2.0,
        us = 1.0,
        z0 = 0.0,
    )
    vpm.set_rhs_pm(RHS_pm_hill)
    print_red(f"Setting RHS_PM as computed from the hill vortex", rank)
    
    if rank == 0:
        st = MPI.Wtime()
        print_red(f"Remeshing")
    XPR, QPR, GPR, UPR = vpm.remesh_particles_3d(-1) 
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing took {int((et - st) / 60)}m {int(et - st) % 60}s\n")

    print_IMPORTANT(f"Particles initialized", rank)
    # Get the particles
    XPR = vpm.particles.XP
    UPR = vpm.particles.UP
    QPR = vpm.particles.QP
    GPR = vpm.particles.GP

    UPR[:,:] = 0
    GPR[:,:] = 0
    if rank != 0:
        XPR[:,:] = 0
        QPR[:,:] = 0
    
    # Create the plot to live update the particles
    if rank == 0:
        plotter = Particle3DPlot()
        plotter.update(
            x = XPR[0,:],
            y = XPR[1,:],
            z = XPR[2,:],
            c = np.sqrt(QPR[0,:]**2 + QPR[1,:]**2 + QPR[2,:]**2)
        )

    comm.Barrier()
    # call vpm with mode = 0
    vpm.vpm(
        num_particles= vpm.particles.NVR,
        num_equations= vpm.num_equations,
        mode = 0,
        particle_positions    =  XPR,
        particle_strengths    =  QPR,
        particle_velocities   =  UPR,
        particle_deformations =  GPR,
        timestep=0,
        viscosity=NI,
    )
    # Main loop
    T = 0
    max_iter = 100
    for i in range(1, max_iter+1):
        comm.Barrier()
        NVR = vpm.particles.NVR
        T += DT
        print_IMPORTANT(
            f"Iteration= {i} of {max_iter}\nT={T}\nDT={DT}",
            rank = rank,
            color_divider="green",
            color_text="green"
        )
        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 2,
            particle_positions    =  XPR,
            particle_strengths    =  QPR,
            particle_velocities   =  UPR,
            particle_deformations =  GPR,
            timestep=i,
            viscosity=NI,
        )

        print_IMPORTANT(f"INFO", rank)
        if rank == 0:
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP
            UPR = vpm.particles.UP
            GPR = vpm.particles.GP
            # Print the size of the particles
            # print(f"Number of particles: {NVR}")
            # print(f"Number of equations: {neq}")
            # print('\n')

            # print_green(f"XPR:")
            # print(f"Mean: {np.mean(XPR.data, axis=1)}")
            # print(f"Max: {np.max(XPR.data, axis=1)}")
            # print(f"Min: {np.min(XPR.data, axis=1)}")
            # print('\n')

            # print_green(f"UPR:")
            # print(f"Mean: {np.mean(UPR.data, axis=1)}")
            # print(f"Max: {np.max(UPR.data, axis=1)}")
            # print(f"Min: {np.min(UPR.data, axis=1)}")
            # print('\n')
            
            # print_green(f"QPR:")
            # print(f"Mean: {np.mean(QPR.data, axis=1)}")
            # print(f"Max: {np.max(QPR.data, axis=1)}")
            # print(f"Min: {np.min(QPR.data, axis=1)}")
            # print('\n')

            # print_green(f"GPR:")
            # print(f"Mean: {np.mean(GPR.data, axis=1)}")
            # print(f"Max: {np.max(GPR.data, axis=1)}")
            # print(f"Min: {np.min(GPR.data, axis=1)}")
            # print('\n')

            # print(f"Particle Mesh Values:")
            # print_green(f"RHS:")
            # RHS = vpm.particle_mesh.RHS
            # print(f"Mean: {np.mean(RHS)}")
            # print(f"Max: {np.max(RHS)}")
            # print(f"Min: {np.min(RHS)}")
            # print('\n')

            # ux = vpm.particle_mesh.Ux
            # uy = vpm.particle_mesh.Uy
            # uz = vpm.particle_mesh.Uz
            # for u, name in zip([ux, uy, uz], ["Ux", "Uy", "Uz"]):
            #     print_green(f"{name}:")
            #     print(f"Mean: {np.mean(u)}")
            #     print(f"Max: {np.max(u)}")
            #     print(f"Min: {np.min(u)}")
            #     print('\n')
            
            # V_mag = np.sqrt(ux**2 + uy**2 + uz**2)
            # print_green(f"Velocity Magnitude:")
            # print(f"Mean: {np.mean(V_mag)}")
            # print(f"Max: {np.max(V_mag)}")
            # print(f"Min: {np.min(V_mag)}")
            # print('\n')

            print_IMPORTANT(f"Convecting Particles", rank)
            # # Move the particles
            for j in range(vpm.particles.NVR):
                XPR[:3, j] += (UPR[:3, j] + UINF) * DT
                FACDEF = 1.0
                QPR[:3, j] -= FACDEF * GPR[:3, j] * DT

            # Update the plot
            plotter.update(
                x = XPR[0,:],
                y = XPR[1,:],
                z = XPR[2,:],
                c = np.sqrt(QPR[0,:]**2 + QPR[1,:]**2 + QPR[2,:]**2)
            ) 

        sleep(0.5)
        print_IMPORTANT(f"Redefine Bounds", rank)

        comm.Barrier()
        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 0,
            particle_positions    =  XPR,
            particle_strengths    =  QPR,
            particle_velocities   =  UPR,
            particle_deformations =  GPR,
            timestep=i,
            viscosity=NI,
        )
        comm.Barrier()

        if i<5:
            print_IMPORTANT(f"Remeshing", rank)
            # Remeshing
            if rank == 0:
                st = MPI.Wtime()
                print_red(f"Remeshing")
            XPR, QPR, UPR, GPR = vpm.remesh_particles_3d(1)
            if rank == 0:
                et = MPI.Wtime()
                print(f"\tRemeshing took {int((et - st) / 60)}m {int(et - st) % 60}s")

    MPI.Finalize()
    end_time = MPI.Wtime()
    print(f"Time taken: {end_time - start_time}")

if __name__ == "__main__":
    main()
