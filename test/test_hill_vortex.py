from time import sleep

import numpy as np
from mpi4py import MPI
from test_hill_problem import hill_assign_parallel

from vpm_py import VPM
from vpm_py.console_io import (print_blue, print_green, print_IMPORTANT,
                               print_red)
from vpm_py.visualization import StandardVisualizer


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
        verbocity= 1,
        dx_particle_mesh= 0.1,
        dy_particle_mesh= 0.1,
        dz_particle_mesh= 0.1
    )
    if rank == 0:
        plotter = StandardVisualizer(plot_particles=("charge","magnitude"))
        vpm.attach_plotter(plotter)

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    print_blue(f"Number of processes: {np_procs}", rank)
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()
     
    DT =  0.5
    neq = 3 

    # Create particles
    NVR = 100
    XPR_zero = np.zeros((3, NVR), dtype=np.float64)
    XPR_zero[:, 0] = np.array([-2., -2., -2.])
    XPR_zero[:, 1] = np.array([2., 2., 2.])
    QPR_zero = np.ones((neq + 1, NVR), dtype=np.float64)

    # Initialization VPM
    comm.Barrier()
    vpm.vpm_define(
        num_equations=neq,
        particle_positions= XPR_zero, 
        particle_charges= QPR_zero, 
    )
    comm.Barrier()

    # Initialize Hill Vortex
    if rank == 0:
        st = MPI.Wtime()
        print_IMPORTANT("Hill vortex initialization", rank)
    _, RHS_pm_hill = hill_assign_parallel(
        Dpm= vpm.dpm,
        NN= vpm.particle_mesh.nn,
        NN_bl= vpm.particle_mesh.nn_bl,
        Xbound= vpm.particle_mesh.xbound,
        neqpm= vpm.num_equations,
        sphere_radius = 2.0,
        u_freestream = 1.0,
        sphere_z_center = 0.0,
    )
    vpm.particle_mesh.set_rhs_pm(RHS_pm_hill)
    print_red("Setting RHS_pm as computed from the hill vortex", rank)
    
    if rank == 0:
        st = MPI.Wtime()
        print_red("Remeshing")
    XPR_hill, QPR_hill = vpm.remesh_particles(project_particles=False) 
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

    print_IMPORTANT("Particles initialized", rank)
    vpm.update_plot(0)

    comm.Barrier()
    vpm.vpm_define(
        num_equations= vpm.num_equations,
        particle_positions  =  XPR_hill[:,:],
        particle_charges    =  QPR_hill[:,:],
    )
    # Main loop
    T = 0
    max_iter = 500
    XPR = XPR_hill.copy()
    QPR = QPR_hill.copy()
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
        vpm.vpm_solve_velocity_deformation(
            timestep=i,
            num_equations=neq,
            particle_positions    =  XPR,
            particle_charges      =  QPR,
        )

        if rank == 0:
            print_IMPORTANT("INFO", rank)
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP
            UPR = vpm.particles.UP
            GPR = vpm.particles.GP
            U_PM = vpm.particle_mesh.U

            for name, u in zip(["Ux", "Uy", "Uz"], U_PM): 
                print_green(f"{name}:")
                print(f"Mean: {np.mean(u)}")
                print(f"Max: {np.max(u)}")
                print(f"Min: {np.min(u)}")
                print('\n')

            print_IMPORTANT("Convecting Particles", rank)
            
            st = MPI.Wtime()
            # Convect the particles and apply vortex stretching
            for j in range(vpm.particles.NVR):
                XPR[:3, j] = XPR[:3, j] + UPR[:3, j] * DT
                QPR[:3, j] = QPR[:3, j] - GPR[:3, j] * DT
            et = MPI.Wtime()

            print(f"\tConvection finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
            print_IMPORTANT("Updating the plot", rank)

            st = MPI.Wtime()
            vpm.update_plot(i)
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Saving the particles and particle mesh", rank)
            st = MPI.Wtime()
            vpm.particles.save_to_file(filename= "particles_test", folder="vpm_case")
            vpm.particle_mesh.save_to_file(filename= "particle_mesh_test", folder="vpm_case")
            et = MPI.Wtime()
            print(f"\tSaving the particles and particle mesh finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

        sleep(0.5)
        print_IMPORTANT("Redefine Bounds", rank)

        comm.Barrier()
        vpm.vpm_define(
            num_equations=neq,
            particle_positions  =  XPR,
            particle_charges    =  QPR,
            timestep=i,
        )

        # vpm.vpm_diffuse(
        #     num_equations=neq,
        #     particle_positions= XPR,
        #     particle_charges= QPR,
        #     particle_velocities= UPR,
        #     particle_deformations= GPR,
        #     viscosity= -0.1
        # )

        XPR, QPR = vpm.remesh_particles(project_particles=True)
        comm.Barrier()

        if i == 1:
            repeate_correction = 20
        else:
            repeate_correction = 1
        for _ in range(repeate_correction):
            if rank == 0:
                XPR = vpm.particles.XP
                QPR = vpm.particles.QP
                NVR = vpm.particles.NVR
            else:
                NVR = None
            vpm.vpm_correct_vorticity(
                num_equations=neq,
                particle_positions=XPR,
                particle_charges=QPR,
                num_particles=NVR,
            )
            comm.Barrier()

    MPI.Finalize()
    end_time = MPI.Wtime()
    print_IMPORTANT(f"Time taken: {int((end_time - start_time) / 60)}m {int(end_time - start_time) % 60}s", rank=rank) 

if __name__ == "__main__":
    main()
