from vpm_py import VPM
from time import sleep
import numpy as np
from mpi4py import MPI
import numpy as np

from vpm_py.console_io import print_IMPORTANT, print_red, print_green, print_blue
from vpm_py.visualization import StandardVisualizer
from test_hill_problem import hill_assign_parallel

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
        verbocity= 2,
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
     
    DT =  0.1
    NI = -0.1
    neq = 3 
    UINF = np.array([0., 0., 0.])

    # Create particles
    NVR = 100
    XPR_zero = np.zeros((3, NVR), dtype=np.float64)
    XPR_zero[:, 0] = np.array([-2., -2., -2.])
    XPR_zero[:, 1] = np.array([2., 2., 2.])
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
        sphere_radius = 1.5,
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
        print(f"\tRemeshing finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

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
            particle_charges      =  QPR,
            timestep=i,
            viscosity=NI,
        )

        if rank == 0:
            print_IMPORTANT(f"INFO", rank)
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP
            UPR = vpm.particles.UP
            GPR = vpm.particles.GP
            # Print the size of the particles
            print(f"Number of particles: {NVR}")
            print(f"Number of equations: {neq}")
            print('\n')

            print_green(f"UPR:")
            print(f"Mean: {np.mean(UPR.data, axis=1)}")
            print(f"Max: {np.max(UPR.data, axis=1)}")
            print(f"Min: {np.min(UPR.data, axis=1)}")
            print('\n')
            
            U_PM = vpm.particle_mesh.U
            for name, u in zip(["Ux", "Uy", "Uz"], U_PM): 
                print_green(f"{name}:")
                print(f"Mean: {np.mean(u)}")
                print(f"Max: {np.max(u)}")
                print(f"Min: {np.min(u)}")
                print('\n')

            print_IMPORTANT(f"Convecting Particles", rank)
            
            st = MPI.Wtime()
            # # Move the particles
            XPR[:3, :] += (UPR[:3, :]) * DT
            # for j in range(vpm.particles.NVR):
            #     # Translate the particles
            #     XPR[:3, j] += (UPR[:3, j] + UINF) * DT
            #     # Diffusion of vorticity
            #     # FACDEF = 1.0
            #     # QPR[:3, j] -= FACDEF * GPR[:3, j] * DT
            et = MPI.Wtime()

            print(f"\tConvection finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
            print_IMPORTANT(f"Updating the plot", rank)

            st = MPI.Wtime()
            # Update the plot
            plotter.update_particle_plots(
                iteration=i,
                particle_positions= XPR[:,:],
                particle_charges= QPR[:,:],
                particle_velocities= UPR[:,:],
                particle_deformations= GPR[:,:]
            )
            plotter.update_mesh_plots(
                iteration=i,
                pm_positions= vpm.particle_mesh.grid_positions,
                pm_velocities= vpm.particle_mesh.U,
                pm_charges= vpm.particle_mesh.RHS,
                pm_deformations= vpm.particle_mesh.deformation
            )
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
            print_IMPORTANT(f"Saving the particles and particle mesh", rank)

            st = MPI.Wtime()
            vpm.particles.save_to_file(filename= f"particles_test", folder="results_test")
            vpm.particle_mesh.save_to_file(filename= f"particle_mesh_test", folder="results_test")
            et = MPI.Wtime()

            print(f"\tSaving the particles and particle mesh finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

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

    MPI.Finalize()
    end_time = MPI.Wtime()
    print_IMPORTANT(f"Time taken: {int((end_time - start_time) / 60)}m {int(end_time - start_time) % 60}s", rank=rank) 

if __name__ == "__main__":
    main()
