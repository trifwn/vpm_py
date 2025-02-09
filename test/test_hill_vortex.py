import numpy as np
from mpi4py import MPI
from test_problems.hill_vortex import hill_assign_parallel

from vpm_py import VPM
from vpm_py.console_io import (print_blue, print_green, print_IMPORTANT,
                               print_red)
from vpm_py.visualization import StandardVisualizer

import psutil
import os

def get_memory_usage():
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return memory_info.rss  # Resident Set Size (RSS) in bytes

def main():
    # PROBLEM STATEMENT
    UINF = np.array([0.0, 0.0, 1.0])
    SPHERE_RADIUS = 1.0
    REYNOLDS_NUMBER = 10.0 
    # Reynolds number = U * L / nu , where U is the velocity, L is the radius of the sphere and nu is the kinematic viscosity
    # nu = U * L / REYNOLDS_NUMBER
    VISCOSITY = np.linalg.norm(UINF) * SPHERE_RADIUS / REYNOLDS_NUMBER
    # DT should be set according to the CFL condition: CFL = U * DT / dx < 1
    dpm = np.array([0.05, 0.05, 0.05])
    DT = 0.5 * (dpm[0]) / np.linalg.norm(UINF)
    TIMESTEPS = 1000
    CFL_LIMITS = [0.3 , 0.9]
    CFL_TARGET = 0.5
    REMESH_FREQUENCY = 40

    # OPTIONS
    remesh = True
    apply_vorticity_correction = True

    # CASE FOLDER
    CASE_FOLDER = "/mnt/c/Users/tryfonas/Data/hill_vortex"
    if REYNOLDS_NUMBER == np.inf:
        CASE_FOLDER += "_Re=inf"
    else:
        CASE_FOLDER += f"_Re={REYNOLDS_NUMBER}"

    if apply_vorticity_correction:
        CASE_FOLDER += "_correct"
    else:
        CASE_FOLDER += "_nocorrect"
    
    if not remesh:
        CASE_FOLDER += "_no_remesh"
    else:
        CASE_FOLDER += f"_remesh_{REMESH_FREQUENCY}"
    
    CASE_FOLDER += "/"

    INITIALIZE_CASE = True 
    # Initialize MPI
    comm = MPI.COMM_WORLD
    start_time = MPI.Wtime()
    rank = comm.Get_rank()
    np_procs = comm.Get_size()
    
    if rank == 0:
        print(f"Case folder: {CASE_FOLDER}")

    # Initialize VPM
    vpm = VPM(
        number_of_equations= 3,
        number_of_processors= np_procs,
        rank= rank,
        verbocity= 1,
        dx_particle_mesh= dpm[0],
        dy_particle_mesh= dpm[1],
        dz_particle_mesh= dpm[2],
        case_folder= CASE_FOLDER,
    )
    if rank == 0:
        plotter = StandardVisualizer(
            plot_particles= ("charge","magnitude"), 
            plot_slices   = [
                 ('velocity', 'magnitude'),
                 ('charge', 'magnitude'),
                #  ("q_pressure","Q"), 
                #  ("u_pressure","U"), 
                 ("pressure","P"), 
            ],
            # Figure size should be 1920x1080
            figure_size= (18.2, 10.0),
        )
        vpm.attach_visualizer(plotter)
        vpm.setup_animation_writer(
            filename=f'{CASE_FOLDER}hill_vortex_pressure.mp4',
            fps=10,
        )
        pass

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    print_blue(f"Number of processes: {np_procs}", rank)
    print_blue(f"Rank: {rank}")

    # Print Problem parameters
    print_red(f"Reynolds number: {REYNOLDS_NUMBER}", rank)
    print_red(f"Viscosity: {VISCOSITY}", rank)
    print_red(f"Sphere radius: {SPHERE_RADIUS}", rank)
    print_red(f"UINF: {UINF}", rank)
    print_red(f"DT: {DT}", rank)
    print_red(f"Approximate CFL: {np.sum(np.linalg.norm(UINF) * DT / vpm.dpm)}", rank)
    neq = 3 

    if INITIALIZE_CASE:
        # Initialize the particles
        XPR_zero, QPR_zero, NVR = initialize_hill_vortex(
            vpm= vpm,
            neq= neq,
            SPHERE_RADIUS= SPHERE_RADIUS,
            UINF= UINF,
            comm= comm,
            rank= rank,
        )
        print_IMPORTANT(f"Initalized {NVR} particles", rank)
        start_iter = 0
    else:
        from vpm_py.file_io import get_latest_particle_file
        XPR_zero, _, QPR_zero, _, start_iter = get_latest_particle_file(
            folder= CASE_FOLDER
        )
        # Load the particles
        NVR = XPR_zero.shape[1]
        comm.Barrier()
        print_IMPORTANT(f"Loaded particles from file: {NVR} particles", rank)

    XPR = XPR_zero.copy()
    QPR = QPR_zero.copy()
    comm.Barrier()
    vpm.vpm_define(
        num_equations= vpm.num_equations,
        particle_positions  =  XPR[:,:],
        particle_charges    =  QPR[:,:],
    )

    # for i in range(1, 2):
    #     vpm.vpm_correct_vorticity(
    #         particle_positions=XPR,
    #         particle_charges=QPR,
    #         num_particles=NVR,
    #     )

    # Main loop
    T = 0
    # Open the file to write the memory usage
    if rank == 0:
        with open(f"{CASE_FOLDER}memory_usage.txt", "w") as f:
            f.write("".join([f"Process {i}, " for i in range(1, np_procs)]) + "\n")  

    for i in range(start_iter, start_iter +  TIMESTEPS+1):
        memory_usage = get_memory_usage()
        if rank == 0:
            memory_usages = np.zeros(np_procs)
            memory_usages[0] = memory_usage / 1024**2
            for j in range(1, np_procs):
                memory_usages[j] = comm.recv(source=j) / 1024**2
            # Append the memory usage to the file
            with open(f"{CASE_FOLDER}memory_usage.txt", "a") as f:
                # Write the values separated by commas
                f.write(",".join([f"{memory:.2f}" for memory in memory_usages]) + "\n")
        else:
            comm.send(memory_usage, dest=0)


        NVR = vpm.particles.NVR
        comm.Barrier()
        grid_dimensions = vpm.particle_mesh.grid_size
        print_IMPORTANT(
            f"Iteration= {i} of {TIMESTEPS}\nT={T}\nDT={DT}\nNumber of particles: {NVR}\n" +
            f"PM Cells: {np.prod(grid_dimensions)}\n" +
            f"CASE_FOLDER: {CASE_FOLDER}\n" ,
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

        vpm.vpm_solve_pressure()

        if rank == 0:
            print_IMPORTANT("INFO", rank)
            XPR = vpm.particles.particle_positions
            QPR = vpm.particles.particle_charges
            UPR = vpm.particles.particle_velocities
            GPR = vpm.particles.particle_deformations
            U_PM = vpm.particle_mesh.velocity
            PRESSURE_PM = vpm.particle_mesh.pressure

            # Position statistics
            for name, x in zip(["X", "Y", "Z"], XPR):
                print_green(f"{name} Statistics:")
                print(f"  Mean: {np.mean(x):.6e}  Max:  {np.max(x):.6e} Min:  {np.min(x):.6e}")
                print()

            # Velocity statistics
            for name, u in zip(["Ux", "Uy", "Uz"], U_PM):
                print_green(f"{name} Statistics:")
                print(f"  Mean: {np.mean(u):.6e}  Max:  {np.max(u):.6e} Min:  {np.min(u):.6e}")
                print()  

            # Pressure statistics
            if PRESSURE_PM is not None:
                for name, p in zip(["P", "Q"], PRESSURE_PM):
                    print_green(f"{name} Statistics:")
                    print(f"  Mean: {np.mean(p):.6e} Max:  {np.max(p):.6e} Min:  {np.min(p):.6e}")
                    print()

            DT = stability_check(DT, VISCOSITY, U_PM, vpm.dpm, CFL_LIMITS, CFL_TARGET)

            print_IMPORTANT("Convecting Particles", rank)
            st = MPI.Wtime()
            # Convect the particles and apply vortex stretching
            XPR[:3,:] += UPR[:3, :] * DT
            QPR[:3,:] -= GPR[:3, :] * DT

            et = MPI.Wtime()
            print(f"\tConvection finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Updating the plot", rank)
            st = MPI.Wtime()
            remesh_str = 'remesh = True' if remesh else 'remesh = False'
            correct_str = 'correction = True' if apply_vorticity_correction else 'correction = False'
            vpm.update_plot(
                f"Reynolds {REYNOLDS_NUMBER} |  Time: {T + DT:.2f}s | Iteration: {i}/{TIMESTEPS} | {remesh_str} | {correct_str}",
                dt = DT
            )
            
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Saving the particles and particle mesh", rank)
            st = MPI.Wtime()
            vpm.particles.save_to_file(filename= "particles", folder=CASE_FOLDER, metadata={"timestep": i, "time": T, "dt": DT})
            vpm.particle_mesh.save_to_file(filename= "particle_mesh", folder=CASE_FOLDER)
            vpm.particle_mesh.save_pressure_to_file(filename= f"{i:06d}_pressure", folder=f"{CASE_FOLDER}/results")
            
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

        if remesh and i % REMESH_FREQUENCY == 0:
            print_IMPORTANT("Remeshing", rank)
            XPR, QPR = vpm.remesh_particles(project_particles=True, cut_off=1e-9)
    
        if apply_vorticity_correction:
            NVR = vpm.particles.NVR
            print_IMPORTANT("Applying Vorticity Correction", rank)
            vpm.vpm_correct_vorticity(
                    particle_positions=XPR,
                    particle_charges=QPR,
                    num_particles=NVR,
                )
        T += DT
        if REYNOLDS_NUMBER != np.inf:
            print_IMPORTANT("Applying Diffusion", rank)
            vpm.vpm_diffuse(viscosity= VISCOSITY)
            if rank == 0:
                QPR = vpm.particles.particle_charges
                GPR = vpm.particles.particle_deformations
                # Diffuse the particles
                QPR[:3,:] = QPR[:3,:] - GPR[:, :] *  DT
                print('Diffusion finished: with viscosity = ', VISCOSITY)
        
        XPR = vpm.particles.particle_positions
        QPR = vpm.particles.particle_charges
        UPR = vpm.particles.particle_velocities
        GPR = vpm.particles.particle_deformations
        
    # Finalize
    end_time = MPI.Wtime()
    print_IMPORTANT(f"Time taken: {int((end_time - start_time) / 60)}m {int(end_time - start_time) % 60}s", rank=rank) 
    MPI.Finalize()

def initialize_hill_vortex(
    vpm: VPM, 
    neq: int, 
    SPHERE_RADIUS: float, 
    UINF, 
    comm, 
    rank
):
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
        NN= vpm.particle_mesh.grid_size,
        NN_bl= vpm.particle_mesh.grid_limits,
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
    NVR = vpm.particles.NVR
    return XPR_hill, QPR_hill, NVR

def stability_check(DT, VISCOSITY, U_PM, dpm, CFL_LIMITS, CFL_TARGET):
    # Calculate CFL components
    CFL_x = np.max(np.abs(U_PM[0, :, :, :])) * DT / dpm[0]
    CFL_y = np.max(np.abs(U_PM[1, :, :, :])) * DT / dpm[1]
    CFL_z = np.max(np.abs(U_PM[2, :, :, :])) * DT / dpm[2]
    CFL = CFL_x + CFL_y + CFL_z

    # Check CFL stability
    print(f"CFL_x: {CFL_x:.3f}, CFL_y: {CFL_y:.3f}, CFL_z: {CFL_z:.3f}")
    if CFL > 1:
        print_red(f"CFL: {CFL:.3f}") 
        print_red("CFL stability criterion violated\n")
    else:
        print_green(f"CFL: {CFL:.3f}")
        print_green("CFL stability criterion satisfied\n")

    # Calculate diffusion stability criterion
    stability_criterion = 1 / (6 * VISCOSITY) / (1 / dpm[0]**2 + 1 / dpm[1]**2 + 1 / dpm[2]**2)
    print('Stability criterion for diffusion:')
    print(f'\tDT < 1/(2ν) * [1/(1/dx² + 1/dy² + 1/dz²)] = {stability_criterion:.3e}')

    # Determine diffusion safety limit
    if DT > stability_criterion:
        print_red('\tNot Satisfied')
        print(f"\tCurrent DT: {DT:.3e} > Stability limit: {stability_criterion:.3e}")
        DT_diffusion = 0.7 * stability_criterion  # Apply safety factor
    else:
        print_green('\tSatisfied')
        print(f"\tCurrent DT: {DT:.3e} < Stability limit: {stability_criterion:.3e}")
        DT_diffusion = stability_criterion  # Use full stability limit

    # Calculate CFL-based timestep boundaries
    DT_CFL_LOW = (CFL_LIMITS[0] / CFL) * DT
    DT_CFL_HIGH = (CFL_LIMITS[1] / CFL) * DT
    DT_CFL_TARGET = (CFL_TARGET / CFL) * DT
    DT_initial = DT  # Store initial DT for logging

    # Determine dominant constraint
    if DT_CFL_HIGH < DT_diffusion:
        print("\nCFL condition is stricter than diffusion")
    elif DT_CFL_LOW > DT_diffusion:
        print("\nDiffusion is stricter than CFL condition")

    # Adjust DT for CFL constraints
    if DT >= DT_CFL_HIGH:
        print_red("\nAdjusting timestep to satisfy CFL upper limit")
        print_red(f"\tDT: {DT:.3e} → {DT_CFL_TARGET:.3e} (CFL target = {CFL_TARGET})")
        DT = DT_CFL_TARGET
    elif DT <= DT_CFL_LOW:
        print_red("\nAdjusting timestep to satisfy CFL lower limit")
        print_red(f"\tDT: {DT:.3e} → {DT_CFL_TARGET:.3e} (CFL target = {CFL_TARGET})")
        DT = DT_CFL_TARGET

    # Adjust DT for diffusion constraints after CFL adjustment
    if DT >= DT_diffusion:
        new_DT = 0.8 * DT_diffusion  # Additional safety factor
        print_red("\nAdjusting timestep for diffusion stability")
        print_red(f"\tDT after CFL: {DT:.3e} → {new_DT:.3e} (0.8 x diffusion limit)")
        DT = new_DT
    else:
        print_green("\nDiffusion stability confirmed")
        print_green(f"\tCurrent DT: {DT:.3e} < Diffusion limit: {DT_diffusion:.3e}")

    # Final validation
    assert DT <= DT_CFL_HIGH, "CFL upper limit violated"
    assert DT <= DT_diffusion, "Diffusion limit violated"
    print(f"\nFinal DT: {DT:.3e} (Initial: {DT_initial:.3e})")
    return DT


if __name__ == "__main__":
    main()
