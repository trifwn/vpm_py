from vpm_py import VPM
from time import sleep
import numpy as np
from mpi4py import MPI
import numpy as np

from vpm_py.console_io import print_IMPORTANT, print_red, print_green, print_blue
from vpm_py.visualization import StandardVisualizer
from test_hill_spherical_vortex import hill_assign_parallel

def main():
    # Initialize MPI
    comm = MPI.COMM_WORLD
    start_time = MPI.Wtime()
    rank = comm.Get_rank()
    np_procs = comm.Get_size()
    
    # Initialize VPM
    vpm = VPM(
        number_of_equations= 3,
        number_of_processors= np_procs,
        rank= rank,
        verbocity= 0,
        dx_particle_mesh= 0.1,
        dy_particle_mesh= 0.1,
        dz_particle_mesh= 0.1
    )
    if rank == 0:
        plotter = StandardVisualizer(plot_particles=("charge","magnitude"))

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
    XPR_zero = np.zeros((3, NVR), dtype=np.float64)
    XPR_zero[:, 0] = np.array([-1.5, -1.5, -1.5])
    XPR_zero[:, 1] = np.array([ 1.5,  1.5,  1.5])
    XPR_zero[:, 0] = np.array([-1.956, -1.9832, -2.04])
    XPR_zero[:, 1] = np.array([2.01, 2.1, 1.89])
    QPR_zero = np.ones((neq + 1, NVR), dtype=np.float64)

    # Initialization VPM
    comm.Barrier()
    vpm.vpm(
        num_equations=neq,
        mode = 0,
        particle_positions= XPR_zero, 
        particle_charges= QPR_zero, 
        timestep=0,
        viscosity=NI,
    )
    comm.Barrier()

    if rank == 0:
        st = MPI.Wtime()

    print_IMPORTANT(f"Hill vortex initialization", rank)
    _, RHS_pm_hill = hill_assign_parallel(
        Dpm= vpm.dpm,
        NN= vpm.particle_mesh.nn,
        NN_bl= vpm.particle_mesh.nn_bl,
        Xbound= vpm.particle_mesh.xbound,
        neqpm= vpm.num_equations,
        sphere_radius = 1.0,
        u_freestream = 1.0,
        sphere_z_center = 0.0,
    )
    vpm.particle_mesh.set_rhs_pm(RHS_pm_hill)
    print_red(f"Setting RHS_pm as computed from the hill vortex", rank)
    
    if rank == 0:
        st = MPI.Wtime()
        print_red(f"Remeshing")
    XPR, QPR = vpm.remesh_particles(project_particles=False) 
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing finished in {int((et - st) / 60)}m {int(et - st) % 60}s\n")

    print_IMPORTANT(f"Particles initialized", rank)

    
    # Create the plot to live update the particles
    if rank == 0:
        plotter.update_particle_plots(
            iteration=0,
            particle_positions= XPR[:,:],
            particle_charges= QPR[:,:],
            particle_velocities= vpm.particles.UP[:,:],
            particle_deformations= vpm.particles.GP[:,:]
        )

    comm.Barrier()
    # call vpm with mode = 0
    vpm.vpm(
        num_equations= vpm.num_equations,
        mode = 0,
        particle_positions    =  XPR,
        particle_charges    =  QPR,
        timestep=0,
        viscosity=NI,
    )
    # Main loop
    T = 0
    max_iter = 500
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
            num_equations=neq,
            mode = 2,
            particle_positions    =  XPR,
            particle_charges    =  QPR,
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

            print_green(f"UPR:")
            print(f"Mean: {np.mean(UPR.data, axis=1)}")
            print(f"Max: {np.max(UPR.data, axis=1)}")
            print(f"Min: {np.min(UPR.data, axis=1)}")
            print('\n')
            
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
                # Translate the particles
                XPR[:3, j] += (UPR[:3, j] + UINF) * DT
                # Diffusion of vorticity
                # FACDEF = 1.0
                # QPR[:3, j] -= FACDEF * GPR[:3, j] * DT

            # Update the plot
            plotter.update_particle_plots(
                iteration=i,
                particle_positions= XPR[:,:],
                particle_charges= QPR[:,:],
                particle_velocities= UPR[:,:],
                particle_deformations= GPR[:,:]
            )
            vpm.particles.save_to_file(filename= f"particles_test", folder="results_test")
            vpm.particle_mesh.save_to_file(filename= f"pm_test", folder="results_test")

        sleep(0.5)
        print_IMPORTANT(f"Redefine Bounds", rank)

        comm.Barrier()
        vpm.vpm(
            num_equations=neq,
            mode = 0,
            particle_positions    =  XPR,
            particle_charges    =  QPR,
            timestep=i,
            viscosity=NI,
        )
        comm.Barrier()

        # if i<5:
        #     print_IMPORTANT(f"Remeshing", rank)
        #     # Remeshing
        #     if rank == 0:
        #         st = MPI.Wtime()
        #         print_red(f"Remeshing")
        #     XPR, QPR, UPR, GPR = vpm.remesh_particles_3d(1)
        #     if rank == 0:
        #         et = MPI.Wtime()
        #         print(f"\tRemeshing finished in {int((et - st) / 60)}m {int(et - st) % 60}s")

    MPI.Finalize()
    end_time = MPI.Wtime()
    print(f"Time taken: {end_time - start_time}")

if __name__ == "__main__":
    main()
