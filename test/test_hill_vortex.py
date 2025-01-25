import numpy as np
from mpi4py import MPI
from hill_vortex.hill_problem import hill_assign_parallel

from vpm_py import VPM
from vpm_py.console_io import (print_blue, print_green, print_IMPORTANT,
                               print_red)
from vpm_py.visualization import StandardVisualizer

# import matplotlib
# matplotlib.use('Agg')


def main():
    # PROBLEM STATEMENT
    UINF = np.array([0.0, 0.0, 1.0])
    REYNOLDS_NUMBER = 10. #np.inf 
    SPHERE_RADIUS = 2.0
    # Reynolds number = U * L / nu , where U is the velocity, L is the radius of the sphere and nu is the kinematic viscosity
    # nu = U * L / REYNOLDS_NUMBER
    VISCOSITY = np.linalg.norm(UINF) * SPHERE_RADIUS / REYNOLDS_NUMBER
    # DT should be set according to the CFL condition: CFL = U * DT / dx < 1
    DT = 0.5 * 0.1 / np.linalg.norm(UINF)
    TIMESTEPS = 1000
    CFL_LIMITS = [0.3 , 0.9]
    CFL_TARGET = 0.5

    # OPTIONS
    remesh = True
    apply_vorticity_correction = False

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
    
    CASE_FOLDER += "/"

    INITIALIZE_CASE = False 
    
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
        dx_particle_mesh= 0.1,
        dy_particle_mesh= 0.1,
        dz_particle_mesh= 0.1,
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
    else:
        from vpm_py.file_io import get_latest_particle_file
        XPR_zero, _, QPR_zero, _ = get_latest_particle_file(
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

    for i in range(1, 2):
        vpm.vpm_correct_vorticity(
            particle_positions=XPR,
            particle_charges=QPR,
            num_particles=NVR,
        )

    # Main loop
    T = 0
    for i in range(1, TIMESTEPS+1):
        comm.Barrier()
        NVR = vpm.particles.NVR
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

        vpm.vpm_solve_pressure()

        if rank == 0:
            print_IMPORTANT("INFO", rank)
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP
            UPR = vpm.particles.UP
            GPR = vpm.particles.GP
            U_PM = vpm.particle_mesh.U
            PRESSURE_PM = vpm.particle_mesh.pressure

            for name, u in zip(["Ux", "Uy", "Uz"], U_PM): 
                print_green(f"{name}:")
                print(f"Mean: {np.mean(u)}")
                print(f"Max: {np.max(u)}")
                print(f"Min: {np.min(u)}")
                print('\n')

            for name, p in zip(["P", "Q"], PRESSURE_PM):
                print_green(f"{name}:")
                print(f"Mean: {np.mean(p)}")
                print(f"Max: {np.max(p)}")
                print(f"Min: {np.min(p)}")
                print('\n')
            
            # Calculate CFL
            CFL_x = np.max(np.abs(U_PM[0, :, :, :])) * DT / vpm.dpm[0]
            CFL_y = np.max(np.abs(U_PM[1, :, :, :])) * DT / vpm.dpm[1]
            CFL_z = np.max(np.abs(U_PM[2, :, :, :])) * DT / vpm.dpm[2]
            CFL   = CFL_x + CFL_y + CFL_z
            print_green(f"CFL: {CFL}")
            print(f"CFL_x: {CFL_x}")
            print(f"CFL_y: {CFL_y}")
            print(f"CFL_z: {CFL_z}")
            print('\n')

            # THE STABILITY CRITERION FOR THE DIFFUSION IS: DT < dx^2 / (2 * nu)
            print('Stability criterion for diffusion:')
            print('\tDT < 1 / (2 * nu) * 1 / (1/dx^2 + 1/dy^2 + 1/dz^2)')
            if DT > 1 / (2 * VISCOSITY) / (1 / vpm.dpm[0]**2 + 1 / vpm.dpm[1]**2 + 1 / vpm.dpm[2]**2):
                print_red('\tNot Satisfied')
                print(f"\tDT: {DT} > {1 / (2 * VISCOSITY) / (1 / vpm.dpm[0]**2 + 1 / vpm.dpm[1]**2 + 1 / vpm.dpm[2]**2)}")
                DT_diffusion = 1 / (2 * VISCOSITY) / (1 / vpm.dpm[0]**2 + 1 / vpm.dpm[1]**2 + 1 / vpm.dpm[2]**2)
            else:
                print_green('\tSatisfied')
                print(f"\tDT: {DT} < {1 / (2 * VISCOSITY) / (1 / vpm.dpm[0]**2 + 1 / vpm.dpm[1]**2 + 1 / vpm.dpm[2]**2)}")
                DT_diffusion = DT
            print('\n')
            
            # Adjust the timestep
            DT = adjust_CFL(CFL, CFL_LIMITS, CFL_TARGET, DT, DT_diffusion)

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
            )
            
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

            print_IMPORTANT("Saving the particles and particle mesh", rank)
            st = MPI.Wtime()
            vpm.particles.save_to_file(filename= "particles", folder=CASE_FOLDER)
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

        if remesh and i % 20 == 0 and i != 0:
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

        if rank == 0:
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP

        if REYNOLDS_NUMBER != np.inf:
            print_IMPORTANT("Applying Diffusion", rank)
            vpm.vpm_diffuse(
                viscosity= VISCOSITY,
                particle_positions= XPR,
                particle_charges= QPR,
            )
            if rank == 0:
                QPR = vpm.particles.QP
                GPR = vpm.particles.GP
                print('Old max(QPR) = ', np.max(np.abs(QPR[:,:])))
                print('Old min(QPR) = ', np.min(np.abs(QPR[:,:])))
                # Diffuse the particles
                QPR[:3,:] = QPR[:3,:] + GPR[:, :] *  DT
                print('New max(QPR) = ', np.max(np.abs(QPR[:,:])))
                print('New min(QPR) = ', np.min(np.abs(QPR[:,:])))
                print('Diffusion finished: with viscosity = ', VISCOSITY)
                print('min (GPR) = ', np.min(GPR[:,:]))
                print('mean(GPR) = ', np.mean(GPR[:,:]))
                print('max (GPR) = ', np.max(GPR[:,:]))

        if rank == 0:
            print('\n')
            print('-----------------------------------------')
            XPR = vpm.particles.XP
            QPR = vpm.particles.QP
        
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
    return XPR_hill, QPR_hill, NVR

def adjust_CFL(CFL, CFL_LIMITS, CFL_TARGET, DT, DT_Limit):
    if CFL > max(CFL_LIMITS) or CFL < min(CFL_LIMITS) :
        DT_old = DT
        DT = CFL_TARGET / CFL * DT
        print_red(f"Adjusting the timestep so that the CFL condition is satisfied and the new CFL is {CFL_TARGET}")
        print_red(f"DT: was {DT_old} -> Adjusting to {DT}")
    elif DT > DT_Limit:
        print_red("Adjusting the timestep so difussion is stable")
        print_red(f"DT: was {DT} -> Adjusting to {DT_Limit}")
        DT = DT_Limit
    else:
        print_green("CFL condition is satisfied")
    print('\n')
    return DT


if __name__ == "__main__":
    main()
