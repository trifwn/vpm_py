

from vpm_py import VPM
import numpy as np
from mpi4py import MPI
import numpy as np

from vpm_py.vpm_io import print_IMPORTANT, print_red, print_green, print_blue
from vpm_py.visualization import Particle3DPlot
from hill_spherical_vortex import hill_assign, visualize_vorticity


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
    verbocity= 2,
    dx_particle_mesh= 0.1,
    dy_particle_mesh= 0.1,
    dz_particle_mesh= 0.1,
)
if rank == 0:
    plotter = Particle3DPlot()

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
XPR_zero[:, 0] = np.array([-2, -2, -2])
XPR_zero[:, 1] = np.array([2, 2, 2])
QPR_zero = np.ones((neq + 1, NVR), dtype=np.float64)
UPR_zero = np.zeros((3, NVR), dtype=np.float64)
GPR_zero = np.zeros((3, NVR), dtype=np.float64)


# Initialization VPM
comm.Barrier()
vpm.vpm(
    num_equations=neq,
    mode = 0,
    particle_positions= XPR_zero, 
    particle_strengths= QPR_zero, 
    particle_velocities= UPR_zero, 
    particle_deformations= GPR_zero,
    timestep=0,
    viscosity=NI,
)
comm.Barrier()

if rank == 0:
    st = MPI.Wtime()

print_IMPORTANT(f"Hill vortex initialization", rank)
_, RHS_pm_hill = hill_assign(
    Dpm= vpm.dpm,
    NN= vpm.nn,
    NN_bl= vpm.nn_bl,
    Xbound= vpm.xbound,
    neqpm= vpm.num_equations,
    sphere_radius = 1.0,
    u_freestream = 1.0,
    sphere_z_center = 0.0,
)


# visualize_vorticity(RHS_pm_hill, vpm.nn_bl)
vpm.set_rhs_pm(RHS_pm_hill)
print_red(f"Setting RHS_PM as computed from the hill vortex", rank)

if rank == 0:
    st = MPI.Wtime()
    print_red(f"Remeshing")
XPR, QPR, GPR, UPR = vpm.remesh_particles_3d(-1) 
if rank == 0:
    et = MPI.Wtime()
    print(f"\tRemeshing finished in {int((et - st) / 60)}m {int(et - st) % 60}s\n")

print_IMPORTANT(f"Particles initialized", rank)

# Create the plot to live update the particles
if rank == 0:
    plotter.update(
        x = XPR[0,:],
        y = XPR[1,:],
        z = XPR[2,:],
        c = np.sqrt(QPR[0,:]**2 + QPR[1,:]**2 + QPR[2,:]**2)
    )

T = 0.
max_iter = 100

from vpm_py.arrays import F_Array

def solve(i,T,XPR):
    comm.Barrier()
    NVR = vpm.particles.NVR
    XPR = vpm.particles.XP
    QPR = vpm.particles.QP
    UPR = vpm.particles.UP
    GPR = vpm.particles.GP
    vpm.vpm(
        num_equations=neq,
        mode = 2,
        particle_positions    =  XPR,
        particle_strengths    =  QPR,
        particle_velocities   =  UPR,
        particle_deformations =  GPR,
        timestep=i,
        viscosity=NI,
    )

    return UPR

def timestep(i, t,
    XPR: F_Array, 
    QPR: F_Array, 
    UPR: F_Array, 
    GPR: F_Array
):
    NVR = vpm.particles.NVR
    print_IMPORTANT(
        f"Iteration= {i} of {max_iter}\nT={T}\nDT={DT}\nNumber of particles={NVR}",
        rank = rank,
        color_divider="green",
        color_text="green"
    )
    print_green(f"RK4 step number 0", rank)
    XPR_TMP = XPR.copy()
    U1 = solve(i, t, XPR)
    
    print_green(f"RK4 step number 1", rank)
    XPR_TMP[:,:] = XPR[:,:] + (DT / 2) * U1[:,:]
    U2 = solve(i, t + DT / 2, XPR_TMP)
    
    print_green(f"RK4 step number 2", rank)
    XPR_TMP[:,:] = XPR[:,:] + (DT / 2) * U2[:,:]
    U3 = solve(i, t + DT / 2, XPR_TMP)
    
    print_green(f"RK4 step number 3", rank)
    XPR_TMP[:,:] = XPR[:,:] + DT * U3[:,:]
    U4 = solve(i, t + DT, XPR_TMP)

    DX = (DT / 6) * (U1[:,:] + 2 * U2[:,:] + 2 * U3[:,:] + U4[:,:])
    

    NVR = vpm.particles.NVR
    if rank == 0:
        XPR = vpm.particles.XP
        QPR = vpm.particles.QP
        UPR = vpm.particles.UP
        GPR = vpm.particles.GP
        # Print the size of the particles

        print_IMPORTANT(f"Convecting Particles", rank)
        # # Move the particles
        print_green(f"U:")
        print(f"Mean: {np.mean(DX/DT, axis=1)}")
        print(f"Max: {np.max(DX/DT, axis=1)}")
        print(f"Min: {np.min(DX/DT, axis=1)}")
        print('\n')
        XPR[:3,:] = XPR[:3,:] + DX[:3,:] 
        
        # Update the plot
        plotter.update(
            x = XPR[0,:],
            y = XPR[1,:],
            z = XPR[2,:],
            c = np.sqrt(QPR[0,:]**2 + QPR[1,:]**2 + QPR[2,:]**2)
        ) 

    print_IMPORTANT(f"Redefine Bounds", rank)
    comm.Barrier()
    vpm.vpm(
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
    # if i%10 == 0:
    #     print_IMPORTANT(f"Remeshing", rank)
    #     # Remeshing
    #     if rank == 0:
    #         st = MPI.Wtime()
    #         print_red(f"Remeshing")
    #     XPR, QPR, UPR, GPR = vpm.remesh_particles_3d(1)
    #     if rank == 0:
    #         et = MPI.Wtime()
    #         print(f"\tRemeshing finished in {int((et - st) / 60)}m {int(et - st) % 60}s")
    return XPR, UPR, QPR, GPR


for i in range(max_iter):
    XPR, UPR, QPR, GPR = timestep(i, T, XPR, UPR, QPR, GPR)
    T += DT
