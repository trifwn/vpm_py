from vpm_py import VPM
import numpy as np
from mpi4py import MPI
import numpy as np

from vpm_py.console_io import print_IMPORTANT, print_red, print_green, print_blue
from vpm_py.visualization import StandardVisualizer
from vpm_py.arrays import F_Array
from test_hill_spherical_vortex import hill_assign_parallel, visualize_vorticity


# Initialize MPI
comm = MPI.COMM_WORLD
start_time = MPI.Wtime()
rank = comm.Get_rank()
np_procs = comm.Get_size()

# Initialize VPM
vpm = VPM(
    number_of_equations=3,
    number_of_processors=np_procs,
    rank=rank,
    verbocity=2,
    dx_particle_mesh=0.1,
    dy_particle_mesh=0.1,
    dz_particle_mesh=0.1,
)

if rank == 0:
    plotter = StandardVisualizer(plot_particles=("charge", 'magnitude'))#, plot_mesh=("velocity", "magnitude"))

# PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
print_blue(f"Number of processes: {np_procs}", rank)
comm.Barrier()
print_blue(f"Rank: {rank}")
comm.Barrier()


DT = 0.1
NI = -0.1
neq = 3
UINF = np.array([0.0, 0.0, 0.0])

# Create particles
NVR = 100
XPR_zero = np.zeros((3, NVR), dtype=np.float64)
XPR_zero[:, 0] = np.array([-2, -2, -2])
XPR_zero[:, 1] = np.array([2, 2, 2])
QPR_zero = np.ones((neq + 1, NVR), dtype=np.float64)

# Initialization VPM
comm.Barrier()
vpm.vpm(
    num_equations=neq,
    mode=0,
    particle_positions=XPR_zero,
    particle_charges=QPR_zero,
    timestep=0,
    viscosity=NI,
)
comm.Barrier()


# # Initialize Hill Vortex
if rank == 0:
    st = MPI.Wtime()
print_IMPORTANT(f"Hill vortex initialization", rank)
_, RHS_pm_hill = hill_assign_parallel(
    Dpm=vpm.dpm,
    NN=vpm.particle_mesh.nn,
    NN_bl=vpm.particle_mesh.nn_bl,
    Xbound=vpm.particle_mesh.xbound,
    neqpm=vpm.num_equations,
    sphere_radius=1.5,
    u_freestream=1.0,
    sphere_z_center=0.0,
)
# visualize_vorticity(RHS_pm_hill, vpm.nn_bl)

vpm.particle_mesh.set_rhs_pm(RHS_pm_hill)
print_red(f"Setting RHS_PM as computed from the hill vortex", rank)

if rank == 0:
    st = MPI.Wtime()
    print_red(f"Remeshing")
XPR_hill, QPR_hill = vpm.remesh_particles(project_particles=False)
if rank == 0:
    et = MPI.Wtime()
    print(f"\tRemeshing finished in {int((et - st) / 60)}m {int(et - st) % 60}s\n")

print_IMPORTANT(f"Particles initialized", rank)

# Create the plot to live update the particles
if rank == 0:
    plotter.update_particle_plots(
        iteration= 0,
        particle_positions= XPR_hill[:,:],
        particle_charges= QPR_hill[:,:],
        particle_velocities= vpm.particles.UP.to_numpy(copy=True),
        particle_deformations= vpm.particles.GP.to_numpy(copy=True),
    )

# # Define Timestep 
def check_redefine(i: int, XPR: F_Array, QPR: F_Array):
    min_particle_poistion_x = np.min(XPR.data[0,:])
    max_particle_poistion_x = np.max(XPR.data[0,:])
    min_particle_poistion_y = np.min(XPR.data[1,:])
    max_particle_poistion_y = np.max(XPR.data[1,:])
    min_particle_poistion_z = np.min(XPR.data[2,:])
    max_particle_poistion_z = np.max(XPR.data[2,:])
    redefine = 0
    if ((min_particle_poistion_x < vpm.particle_mesh.Xmin + 2 * vpm.dpm[0]) or 
        (min_particle_poistion_y < vpm.particle_mesh.Ymin + 2 * vpm.dpm[1]) or 
        (min_particle_poistion_z < vpm.particle_mesh.Zmin + 2 * vpm.dpm[2]) or
        (max_particle_poistion_x > vpm.particle_mesh.Xmax - 2 * vpm.dpm[0]) or
        (max_particle_poistion_y > vpm.particle_mesh.Ymax - 2 * vpm.dpm[1]) or
        (max_particle_poistion_z > vpm.particle_mesh.Zmax - 2 * vpm.dpm[2])):
        redefine = 1
        print_IMPORTANT("Redefining the particle mesh because the particles are out of bounds", rank)
    # Broadcast the redefine flag
    comm = MPI.COMM_WORLD
    redefine = comm.bcast(redefine, root=0)

    if redefine == 1:
        vpm.vpm(
            num_equations=neq,
            mode = 0,
            particle_positions = XPR,
            particle_charges = QPR,
            timestep=i,
            viscosity=NI,
        )
        comm.Barrier()


max_iter = 500
def solve(i: int, T: float, XPR: F_Array, QPR: F_Array):
    # Check if the particles are out of bounds
    check_redefine(i, XPR, QPR)
    
    NVR = vpm.particles.NVR
    print_IMPORTANT(
        f"Iteration= {i} of {max_iter}\nT={T}\nDT={DT}\nNumber of particles={NVR}",
        rank=rank,
        color_divider="green",
        color_text="green",
    )
    vpm.vpm(
        num_equations=neq,
        mode=2,
        particle_positions=XPR,
        particle_charges=QPR,
        timestep=i,
        viscosity=NI,
    )

    if rank == 0:
        UPR = vpm.particles.UP.to_numpy(copy=True)
        GPR = vpm.particles.GP.to_numpy(copy=True)
        # Update the plot
        plotter.update_particle_plots(
            iteration=i,
            particle_positions=XPR[:,:],
            particle_charges=QPR[:,:],
            particle_velocities=UPR,
            particle_deformations=GPR,
        )
        # plotter.update_mesh_plots(
        #     iteration=i,
        #     mesh_velocity=vpm.particle_mesh.
        #     mesh_vorticity=vpm.particle_mesh.
        # )
    else:
        UPR = vpm.particles.particle_velocities
        GPR = vpm.particles.particle_deformations

    comm.Barrier()
    return UPR 

def timestep(
    i:int, 
    t:float, 
    XPR: F_Array, 
    QPR: F_Array, 
):
    XPR_TMP = XPR.copy()
    QPR_TMP = QPR.copy()
    U1 = solve(i, t, XPR_TMP, QPR_TMP)
    if rank == 0:
        if np.any(np.isnan(U1)):
            print_red(f"U1 has NaN values", rank)
            raise ValueError("U1 has NaN values")

    XPR_TMP.data[:, :] = XPR.data[:, :] + (DT / 2) * U1[:, :]
    U2 = solve(i, t + DT / 2, XPR_TMP, QPR_TMP)
    if rank == 0:
        if np.any(np.isnan(U2)):
            print_red(f"U2 has NaN values", rank)
            raise ValueError("U2 has NaN values")

    XPR_TMP.data[:, :] = XPR.data[:, :] + (DT / 2) * U2[:, :]
    U3 = solve(i, t + DT / 2, XPR_TMP, QPR_TMP)
    if rank == 0:
        if np.any(np.isnan(U3)):
            print_red(f"U3 has NaN values", rank)
            raise ValueError("U3 has NaN values")

    XPR_TMP.data[:, :] = XPR.data[:, :] + DT * U3[:, :]
    U4 = solve(i, t + DT, XPR_TMP, QPR_TMP)
    if rank == 0:
        if np.any(np.isnan(U4)):
            print_red(f"U4 has NaN values", rank)
            raise ValueError("U4 has NaN values")

    # U_mean = U1
    if rank == 0:
        U_mean = 1/6 * (U1[:, :] + 2 * U2[:, :] + 2 * U3[:, :] + U4[:, :])
        if np.any(np.isnan(U_mean)):
            print_red(f"U_mean has NaN values", rank)
            print_red(f"U1: {U1}", rank)
            print_red(f"U2: {U2}", rank)
            print_red(f"U3: {U3}", rank)
            print_red(f"U4: {U4}", rank)
            raise ValueError("U_mean has NaN values")
        print_IMPORTANT(f"Convecting Particles", rank)
        # We do not convect the vorticity
        for j in range(vpm.particles.NVR):
            XPR[:, j] = XPR[:, j] + U_mean[:,j] * DT

        print_IMPORTANT(f"Saving to file", rank)
        vpm.particles.save_to_file(filename=f"particles")
        vpm.particle_mesh.save_to_file(filename=f"particle_mesh")

    print_IMPORTANT(f"Redefine Bounds", rank)
    vpm.vpm(
        num_equations=neq,
        mode=0,
        particle_positions=XPR,
        particle_charges=QPR,
        timestep=i,
        viscosity=NI,
    )
    comm.Barrier()
    return XPR, QPR

# # Simulate
T = 0.
max_iter = 500
XPR = XPR_hill 
QPR = QPR_hill
for i in range(max_iter):
    XPR, QPR = timestep(i, T, XPR, QPR)
    T += DT
