from vpm_py import VPM
from time import sleep
from scipy.io import FortranFile
import numpy as np
from mpi4py import MPI
import numpy as np

def print_blue(text):
    print(f"\033[94m{text}\033[00m")

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_red(text):
    print(f"\033[91m{text}\033[00m")

def print_IMPORTANT(text):
    print(f"\033[93m{'-'*100}\033[00m")
    print(f"\033[91m{text}\033[00m")
    print(f"\033[93m{'-'*100}\033[00m")

def main():
    vpm = VPM()
    print_IMPORTANT(
        vpm.original_lib
    )
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    np_procs = comm.Get_size()

    # PRINT THE RANK OF THE PROCESS AND DETERMINE HOW MANY PROCESSES ARE RUNNING
    if rank == 0:
        print_blue(f"Number of processes: {np_procs}")
    comm.Barrier()
    print_blue(f"Rank: {rank}")
    comm.Barrier()
     
    DT = 0.1
    NI = 0.1
    neq = 3 

    # Create particles
    NVR_ext = np.int32(100)
    Vref = np.float64(1.)
    Vref = np.float64(Vref)      # Convert Vref to 64-bit floating point

    XPR = np.zeros((3, NVR_ext), dtype=np.float64)
    QPR = np.ones((neq + 1, NVR_ext), dtype=np.float64)
    
    # Initialize particle
    x_pos_min = -10
    x_pos_max = 10 
    for i in range(NVR_ext):
        x = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        y = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        z = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        XPR[:, i] = np.array([x,y,z], dtype=np.float64)
        QPR[:3, i] = np.array([1.,2.,1.], dtype=np.float64) * 10
    

    if rank == 0: 
        # Open the file for writing in binary mode
        with FortranFile("particles.bin", 'w') as f:
            # Write NVR_ext and Vref
            f.write_record(NVR_ext)
            f.write_record(Vref)
            
            # Write XPR and QPR arrays
            for i in range(NVR_ext):
                f.write_record(XPR[:, i])
                f.write_record(QPR[:, i])
        print_green(f"Successfully wrote {NVR_ext} particles to 'particles.bin'.")
 

    MAX_PARTICLES = NVR_ext + 100

    # Initialization VPM
    UINF = np.zeros(3, dtype=np.float64)
    UINF[0] = Vref/10
    RHS_pm_in = np.zeros((neq, 10, 10, 10), dtype=np.float64) # Placeholder, adjust as needed
    U_PM = np.zeros((3,10,10,10), dtype=np.float64) # Placeholder, adjust as needed
    
    comm.Barrier()
    vpm.vpm(
        num_particles=NVR_ext,
        num_equations=neq,
        mode = 0,
        particle_positions= XPR, 
        particle_strengths=QPR, 
        particle_velocities=np.zeros((3, NVR_ext)), 
        particle_deformations=np.zeros((3, NVR_ext)), 
        RHS_PM=RHS_pm_in,
        U_PM=U_PM,
        timestep=0,
        viscosity=NI,
        max_particle_num=MAX_PARTICLES
    )
    comm.Barrier()

    # Remeshing
    if rank == 0:
        st = MPI.Wtime()

    vpm.remesh_particles_3d(1) 

    if rank == 0:
        et = MPI.Wtime()
        print(f"Remeshing {int((et - st) / 60)}m {int(et - st) % 60}s")

    # Main loop
    NVR = NVR_ext
    T = 0
    TMAX = 20
    for i in range(1, 100):
        comm.Barrier()
        T += DT
        if rank == 0:
            print(f"---------------------------------")
            print_green(f"ITERATION= {i} of {TMAX}")
            print_green(f"T={DT*i}")
            print_green(f"DT={DT}")
            print(f"---------------------------------")

        UPR = np.zeros((3, NVR), dtype=np.float64)
        GPR = np.zeros((3, NVR), dtype=np.float64)

        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 2,
            particle_positions= XPR, 
            particle_strengths=QPR, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            U_PM=U_PM,
            timestep=i,
            viscosity=NI,
            max_particle_num=NVR_ext
        )

        if rank == 0:

            print(f"XPR:")
            # Print mean/max/min of XPR
            print(f"Mean: {np.mean(XPR, axis=1)}")
            print(f"Max: {np.max(XPR, axis=1)}")
            print(f"Min: {np.min(XPR, axis=1)}")
            print('\n\n')
            print(f"UPR:")
            # Print mean/max/min of UPR
            print(f"Mean: {np.mean(UPR, axis=1)}")
            print(f"Max: {np.max(UPR, axis=1)}")
            print(f"Min: {np.min(UPR, axis=1)}")
            print('\n\n')
            print(f"GPR:")
            # Print mean/max/min of GPR
            print(f"Mean: {np.mean(GPR, axis=1)}")
            print(f"Max: {np.max(GPR, axis=1)}")
            print(f"Min: {np.min(GPR, axis=1)}")
            print('\n\n')
            print(f"U_PM:")
            # Print mean/max/min of U_PM
            U = np.sqrt(U_PM[0,:,:,:]**2 + U_PM[1,:,:,:]**2 + U_PM[2,:,:,:]**2).flatten()
            print(f"Mean: {np.mean(U)}")
            print(f"Max: {np.max(U)}")
            print(f"Min: {np.min(U)}")
            print('\n\n')
        sleep(5)
        comm.Barrier()

        if rank == 0:
            for j in range(NVR_ext):
                XPR[:, j] += (UPR[:, j] + UINF) * DT
                FACDEF = 1.0
                QPR[:3, j] -= FACDEF * GPR[:3, j] * DT

        vpm.vpm(
            num_particles=NVR,
            num_equations=neq,
            mode = 0,
            particle_positions= XPR, 
            particle_strengths=QPR, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            U_PM=U_PM,
            timestep=i,
            viscosity=NI,
            max_particle_num=MAX_PARTICLES
        )
        comm.Barrier()

        if rank == 0:
            st = MPI.Wtime()
            print_red(f"Remeshing")
        vpm.remesh_particles_3d(1)
        if rank == 0:
            et = MPI.Wtime()
            print(f"\tRemeshing took {int((et - st) / 60)}m {int(et - st) % 60}s")

    MPI.Finalize()

if __name__ == "__main__":
    main()
