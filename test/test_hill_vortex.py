from time import sleep

import numpy as np
from mpi4py import MPI
from utils.hill_problem import hill_assign_parallel

from vpm_py import VPM
from vpm_py.console_io import (print_blue, print_green, print_IMPORTANT,
                               print_red)
from vpm_py.visualization import StandardVisualizer


def main():
    # PROBLEM STATEMENT
    DT =  0.01 
    UINF = np.array([0.0, 0.0, 1.0])
    REYNOLDS_NUMBER = 100
    SPHERE_RADIUS = 2.0
    # Reynolds number = U * L / nu , where U is the velocity, L is the radius of the sphere and nu is the kinematic viscosity
    # nu = U * L / REYNOLDS_NUMBER
    VISCOSITY = np.linalg.norm(UINF) * 2.0 / REYNOLDS_NUMBER
    TIMESTEPS = 1500

    # CASE FOLDER
    CASE_FOLDER = "/mnt/c/Users/tryfonas/Data/hill_vortex/"

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
        dz_particle_mesh= 0.1,
        case_folder= CASE_FOLDER,
    )
    if rank == 0:
        plotter = StandardVisualizer(
            plot_particles=("charge","magnitude"), 
            plot_slices=("velocity","magnitude") 
        )
        vpm.attach_visualizer(plotter)
        vpm.setup_animation_writer(
            filename=f'{CASE_FOLDER}hill_vortex.mp4',
            fps=10,
        )

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    print_blue(f"Number of processes: {np_procs}", rank)
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()

    # Print Problem parameters
    print_red(f"Reynolds number: {REYNOLDS_NUMBER}", rank)
    print_red(f"Viscosity: {VISCOSITY}", rank)
    print_red(f"Sphere radius: {SPHERE_RADIUS}", rank)
    print_red(f"UINF: {UINF}", rank)
    print_red(f"DT: {DT}", rank)
    print_red(f"Approximate CFL: {np.sum(np.linalg.norm(UINF) * DT / vpm.dpm)}", rank)
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
        sphere_radius = SPHERE_RADIUS,
        u_freestream = UINF[2],
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
    # vpm.update_plot(0)

    comm.Barrier()
    vpm.vpm_define(
        num_equations= vpm.num_equations,
        particle_positions  =  XPR_hill[:,:],
        particle_charges    =  QPR_hill[:,:],
    )
    # Main loop
    T = 0
    XPR = XPR_hill.copy()
    QPR = QPR_hill.copy()
    for i in range(1, TIMESTEPS+1):
        comm.Barrier()
        NVR = vpm.particles.NVR
        T += DT
        print_IMPORTANT(
            f"Iteration= {i} of {TIMESTEPS}\nT={T}\nDT={DT}",
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
            XPR[:3,:] += UPR * DT
            QPR[:3,:] -= GPR * DT
            et = MPI.Wtime()
            print(f"\tConvection finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Updating the plot", rank)
            st = MPI.Wtime()
            vpm.update_plot(f"Time: {T:.2f}s | Iteration: {i}/{TIMESTEPS}")
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Saving the particles and particle mesh", rank)
            st = MPI.Wtime()
            vpm.particles.save_to_file(filename= "particles_test", folder=CASE_FOLDER)
            vpm.particle_mesh.save_to_file(filename= "particle_mesh_test", folder=CASE_FOLDER)
            et = MPI.Wtime()
            print(f"\tSaving the particles and particle mesh finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

        print_IMPORTANT("Redefine Bounds", rank)
        comm.Barrier()
        vpm.vpm_define(
            num_equations=neq,
            particle_positions  =  XPR,
            particle_charges    =  QPR,
            timestep=i,
        )

        vpm.vpm_diffuse(
            num_equations=neq,
            particle_positions= XPR,
            particle_charges= QPR,
            viscosity= -VISCOSITY
        )

        XPR, QPR = vpm.remesh_particles(project_particles=True)
        # if i == 1:
        #     repeate_correction = 20
        # else:
        #     repeate_correction = 1
        # for _ in range(repeate_correction):
        #     if rank == 0:
        #         XPR = vpm.particles.XP
        #         QPR = vpm.particles.QP
        #         NVR = vpm.particles.NVR
        #     else:
        #         NVR = None
        #     vpm.vpm_correct_vorticity(
        #         num_equations=neq,
        #         particle_positions=XPR,
        #         particle_charges=QPR,
        #         num_particles=NVR,
        #     )
        comm.Barrier()

    # Finalize
    end_time = MPI.Wtime()
    print_IMPORTANT(f"Time taken: {int((end_time - start_time) / 60)}m {int(end_time - start_time) % 60}s", rank=rank) 
    MPI.Finalize()

if __name__ == "__main__":
    main()
