import numpy as np
from mpi4py import MPI
from test_problems.ricker_wavelet import ricker_wavelet_3d_with_gradient 
from vpm_py import VPM
from vpm_py.console_io import (print_blue, print_IMPORTANT, print_red)
from vpm_py.visualization import StandardVisualizer

# import matplotlib
# matplotlib.use('Agg')


def main():
    # PROBLEM STATEMENT
    UINF = np.array([0.0, 0.0, 1.0])
    REYNOLDS_NUMBER = 10. #np.inf 
    SPHERE_RADIUS = 1.0
    # Reynolds number = U * L / nu , where U is the velocity, L is the radius of the sphere and nu is the kinematic viscosity
    # nu = U * L / REYNOLDS_NUMBER
    VISCOSITY = 100 #np.linalg.norm(UINF) * SPHERE_RADIUS / REYNOLDS_NUMBER
    dpm = np.array([0.1, 0.1, 0.1])
    # DT should be set according to the CFL condition: CFL = U * DT / dx < 1
    DT =  1 / (6 * VISCOSITY) / (1 / dpm[0]**2 + 1 / dpm[1]**2 + 1 / dpm[2]**2)

    TIMESTEPS = 500

    # OPTIONS
    remesh = True

    # CASE FOLDER
    CASE_FOLDER = f"/mnt/c/Users/tryfonas/Data/ricker_wavelet_nu={VISCOSITY}"

    if not remesh:
        CASE_FOLDER += "_no_remesh"
    
    CASE_FOLDER += "/"

    INITIALIZE_CASE = True #False 
    
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
                #  ('velocity', 'magnitude'),
                 ('charge', 'magnitude'),
                 ('charge', 'x'),
                 ('charge', 'y'),
                 ('charge', 'z'),
                #  ("q_pressure","Q"), 
                #  ("u_pressure","U"), 
                #  ("pressure","P"), 
            ],
            # Figure size should be 1920x1080
            figure_size= (18.2, 10.0),
        )
        vpm.attach_visualizer(plotter)
        vpm.setup_animation_writer(
            filename=f'{CASE_FOLDER}rick_wavelet_diffusion.mp4',
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
        XPR_zero, QPR_zero, NVR = initialize_ricker_wavelet(
            vpm= vpm,
            neq= neq,
            radius= SPHERE_RADIUS,
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
            st = MPI.Wtime()
            vpm.update_plot(f"Nu = {VISCOSITY} |  Time: {T + DT:.2f}s | Iteration: {i}/{TIMESTEPS}")
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

        if remesh and i % 20 == 0 and i != 0:
            print_IMPORTANT("Remeshing", rank)
            XPR, QPR = vpm.remesh_particles(project_particles=True, cut_off=1e-9)
    
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
                QPR[:3,:] = QPR[:3,:] - GPR[:, :] *  DT
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

def initialize_ricker_wavelet(
    vpm: VPM, 
    neq: int, 
    radius: float, 
    comm, 
    rank
):
    # Create particles
    NVR = 100
    XPR_zero = np.zeros((3, NVR), dtype=np.float64)
    XPR_zero[:, 0] = -5.0 * radius 
    XPR_zero[:, 1] = 5.0 * radius
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
    st = MPI.Wtime()
    print_IMPORTANT("Ricker Wavelet initialization", rank)
    _, RHS_pm_ricker = ricker_wavelet_3d_with_gradient(
        NN= vpm.particle_mesh.grid_size,
        Xbound= vpm.particle_mesh.xbound,
        Dpm= vpm.dpm,
        radius = radius,
    )
    et = MPI.Wtime()
    print(f"\tRicker Wavelet initialized in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
    vpm.particle_mesh.set_rhs_pm(RHS_pm_ricker)
    print_red("Setting RHS_pm as computed from the hill vortex", rank)
    
    if rank == 0:
        st = MPI.Wtime()
        print_red("Remeshing")
    NVR, XPR_hill, QPR_hill = vpm.remesh_particles(project_particles=False) 
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

    print_IMPORTANT("Particles initialized", rank)
    return XPR_hill, QPR_hill, NVR

if __name__ == "__main__":
    main()
