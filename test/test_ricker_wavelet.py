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
    REYNOLDS_NUMBER = 1. #np.inf 
    SPHERE_RADIUS = 1.0
    # Reynolds number = U * L / nu , where U is the velocity, L is the radius of the sphere and nu is the kinematic viscosity
    # nu = U * L / REYNOLDS_NUMBER
    VISCOSITY = 100 #np.linalg.norm(UINF) * SPHERE_RADIUS / REYNOLDS_NUMBER
    dpm = np.array([0.2, 0.2, 0.2])
    # DT should be set according to the CFL condition: CFL = U * DT / dx < 1
    DT =  1 / (12 * VISCOSITY) / (1 / dpm[0]**2 + 1 / dpm[1]**2 + 1 / dpm[2]**2)

    TIMESTEPS = 1000

    # OPTIONS
    remesh = True
    REMESH_FREQUENCY = 40
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
                 ('charge', 'magnitude'),
                 ('charge', 'x'),
                 ('charge', 'y'),
                 ('charge', 'z'),
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

        if rank == 0:
            st = MPI.Wtime()
            vpm.update_plot(f"Nu = {VISCOSITY} |  Time: {T + DT:.2f}s | Iteration: {i}/{TIMESTEPS}")
            et = MPI.Wtime()
            print(f"\tUpdating the plot finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
            vpm.particles.save_to_file(filename= "particles", folder=CASE_FOLDER, metadata={"timestep": i, "time": T, "dt": DT})
            vpm.particle_mesh.save_to_file(filename= "particle_mesh", folder=CASE_FOLDER)


        if remesh and i % REMESH_FREQUENCY == 0 and i != 0:
            print_IMPORTANT("Remeshing", rank)
            XPR, QPR = vpm.remesh_particles(project_particles=True, cut_off=1e-9)
    
        T += DT
        if rank == 0:
            XPR = vpm.particles.particle_positions
            QPR = vpm.particles.particle_charges

        if REYNOLDS_NUMBER != np.inf:
            print_IMPORTANT("Applying Diffusion", rank)
            vpm.vpm_diffuse(
                viscosity= VISCOSITY,
                particle_positions= XPR,
                particle_charges= QPR,
            )
            if rank == 0:
                QPR = vpm.particles.particle_charges
                GPR = vpm.particles.particle_deformations
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
            XPR = vpm.particles.particle_positions
            QPR = vpm.particles.particle_charges
        
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
    wavelet_field, gradient_field = ricker_wavelet_3d_with_gradient(
        NN= vpm.particle_mesh.grid_size,
        Xbound= vpm.particle_mesh.xbound,
        Dpm= vpm.dpm,
        t = 0.0,
        nu = 0.0,
        sigma= radius 
    )

    et = MPI.Wtime()
    print(f"\tRicker Wavelet initialized in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")
    vpm.particle_mesh.set_rhs_pm(gradient_field)
    print_red("Setting RHS_pm as computed from the hill vortex", rank)
    
    if rank == 0:
        st = MPI.Wtime()
        print_red("Remeshing")
    XPR_hill, QPR_hill = vpm.remesh_particles(project_particles=False) 
    NVR = vpm.particles.NVR
    if rank == 0:
        et = MPI.Wtime()
        print(f"\tRemeshing finished in {int((et - st) / 60)}m {(et - st) % 60:.2f}s\n")

    print_IMPORTANT("Particles initialized", rank)
    return XPR_hill, QPR_hill, NVR

if __name__ == "__main__":
    main()
