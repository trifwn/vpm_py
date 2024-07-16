from vpm_py import VPM
import numpy as np

from mpi4py import MPI
import numpy as np

def print_blue(text):
    print(f"\033[94m{text}\033[00m")

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_IMPORTANT(text):
    print(f"\033[93m{'-'*100}\033[00m")
    print(f"\033[91m{text}\033[00m")
    print(f"\033[93m{'-'*100}\033[00m")

def main():
    vpm = VPM()
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
    
    # Initialize variables
    NI_in = 0.1
    DT_in = 0.05
    neq = 3

    # Initialize variables for remeshing
    Vref = None
    NVR_ext = None
    XPR = None
    QPR = None
     
    if rank < 1e15:
        # Read particles
        # with open('particles.bin', 'rb') as f:
        #     NVR_ext = np.fromfile(f, dtype=np.int32, count=1)[0]
        #     Vref = np.fromfile(f, dtype=np.float64, count=1)[0]
        #     XPR = np.zeros((3, NVR_ext), dtype=np.float64)
        #     QPR = np.zeros((neq + 1, NVR_ext), dtype=np.float64)
        #     for i in range(NVR_ext):
        #         XPR[:, i] = np.fromfile(f, dtype=np.float64, count=3)
        #         QPR[:3, i] = np.fromfile(f, dtype=np.float64, count=3)
        
        # Create 2 particles
        NVR_ext = np.int32(100)
        Vref = np.float64(0.1)

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
            QPR[:3, i] = np.array([10,10,10], dtype=np.float64)
        
        Vref = np.float64(Vref)      # Convert Vref to 64-bit floating point

        if rank == 0: 
            # Open the file for writing in binary mode
            from scipy.io import FortranFile
            with FortranFile("particles.bin", 'w') as f:
                # Write NVR_ext and Vref
                f.write_record(NVR_ext)
                f.write_record(Vref)
                
                # Write XPR and QPR arrays
                for i in range(NVR_ext):
                    f.write_record(XPR[:, i])
                    f.write_record(QPR[:, i])
                    print(f"Successfully wrote {NVR_ext} particles to 'particles.bin'.")
 
        RMETM = 0.001
        NVR_sources = 0
        XSOUR = np.zeros((3, NVR_sources), dtype=np.float64)
        QSOUR = np.zeros((neq + 1, NVR_sources), dtype=np.float64)

    # Broadcast NVR_ext and Vref
    NVR_ext = comm.bcast(NVR_ext, root=0)
    MAX_PARTICLES = NVR_ext + 100
    Vref = comm.bcast(Vref, root=0)

    # Initialization VPM
    UINF = np.zeros(3, dtype=np.float64)
    UINF[0] = 7
    RHS_pm_in = np.zeros((neq, 10, 10, 10), dtype=np.float64) # Placeholder, adjust as needed
    velx = np.zeros((10, 10, 10), dtype=np.float64) # Placeholder, adjust as needed
    vely = np.zeros((10, 10, 10), dtype=np.float64) # Placeholder, adjust as needed
    velz = np.zeros((10, 10, 10), dtype=np.float64) # Placeholder, adjust as needed
    U_PM = np.zeros((3,10,10,10), dtype=np.float64) # Placeholder, adjust as needed
    # Assuming vpm function is wrapped from Fortran and available in Python
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
        viscosity=NI_in,
        max_particle_num=NVR_ext
    )
    comm.Barrier()

    # Remeshing
    if rank == 0:
        st = MPI.Wtime()

    vpm.remesh_particles_3d(1)  # Assuming this function is available in Python

    if rank == 0:
        et = MPI.Wtime()
        print(f"Remeshing {int((et - st) / 60)}m {int(et - st) % 60}s")

    # Main loop
    NVR_all = NVR_ext
    T = 0
    for i in range(1, 6001):
        T += DT_in

        if rank < 1e15:
            NVR_all = NVR_ext #+ NVR_sources
            XP_all = np.zeros((3, NVR_all), dtype=np.float64)
            QP_all = np.zeros((neq + 1, NVR_all), dtype=np.float64)
            UPR = np.zeros((3, NVR_all), dtype=np.float64)
            GPR = np.zeros((3, NVR_all), dtype=np.float64)
            XP_all[:, :NVR_ext] = XPR
            QP_all[:, :NVR_ext] = QPR
            XP_all[:, NVR_ext:NVR_all] = XSOUR
            QP_all[:, NVR_ext:NVR_all] = QSOUR
            print_green(
                f"NVR_all: {NVR_all}, XP_all: {XP_all.shape}, QP_all: {QP_all.shape}, UPR: {UPR.shape}, GPR: {GPR.shape}"
            )
        NVR_all = comm.bcast(NVR_all, root=0)

        vpm.vpm(
            num_particles=NVR_all,
            num_equations=neq,
            mode = 1,
            particle_positions= XP_all, 
            particle_strengths=QP_all, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            U_PM=U_PM,
            timestep=i,
            viscosity=NI_in,
            max_particle_num=NVR_ext
        )

        if rank == 0:
            max_val = np.max(UPR)
            print(max_val)

            for j in range(NVR_ext):
                XPR[:, j] += (UPR[:, j] + UINF) * DT_in
                FACDEF = 1.0
                QPR[:3, j] -= FACDEF * GPR[:3, j] * DT_in

            for j in range(NVR_sources):
                XSOUR[:, j] += UINF * DT_in

        vpm.vpm(
            num_particles=NVR_all,
            num_equations=neq,
            mode = 0,
            particle_positions= XP_all, 
            particle_strengths=QP_all, 
            particle_velocities=UPR, 
            particle_deformations=GPR, 
            RHS_PM=RHS_pm_in,
            U_PM=U_PM,
            timestep=i,
            viscosity=NI_in,
            max_particle_num=NVR_ext
        )

        if rank == 0:
            st = MPI.Wtime()

        if i % 1 == 0:
            vpm.remesh_particles_3d(1)

        if rank == 0:
            et = MPI.Wtime()
            print(f"Remeshing {int((et - st) / 60)}m {int(et - st) % 60}s")

    MPI.Finalize()

def definevort(RHS_pm, Xbound, Dpm, NN, NN_bl):
    # Constants
    xc1, xc2 = 0.0, 0.0
    yc1, yc2 = -1.5, 1.5
    PI = 4.0 * np.arctan(1.0)
    
    # Initialize analytic_sol array
    analytic_sol = np.zeros((1, NN[0], NN[1], NN[2]))
    
    # Iterate over indices
    for j in range(NN_bl[1], NN_bl[4]+1):
        for i in range(NN_bl[0], NN_bl[3]+1):
            xi = Xbound[0] + (i-1) * Dpm[0]
            yi = Xbound[1] + (j-1) * Dpm[1]
            ksi1 = np.sqrt((xi - xc1)**2 + (yi - yc1)**2)
            ksi2 = np.sqrt((xi - xc2)**2 + (yi - yc2)**2)
            th1 = np.arctan2((yi - yc1), (xi - xc1))
            th2 = np.arctan2((yi - yc2), (xi - xc2))
            if th1 < 0.0:
                th1 += 2.0 * PI
            if th2 < 0.0:
                th2 += 2.0 * PI
            w1 = (2.0 - ksi1**2) * np.exp(0.5 * (1.0 - ksi1**2))
            w2 = -(2.0 - ksi2**2) * np.exp(0.5 * (1.0 - ksi2**2))
            RHS_pm[0, i-1, j-1, 0] = -(w1 + w2)
            analytic_sol[0, i-1, j-1, 0] = np.exp(0.5 * (1.0 - ksi1**2)) - np.exp(0.5 * (1.0 - ksi2**2))
    return analytic_sol, 

if __name__ == "__main__":
    main()
