import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI

def hills_vortex(control_point, sphere_radius, u_freestream, z_0):    
    """

    Args:
        CP (_type_): Control Point
        a: The radius of the sphere 
        us: The constant free stream velocity far away from the sphere
        z0: The z-coordinate of the sphere center

    Returns:
        _type_: _description_
    """
    # Unpacking control point
    x, y, z = control_point
    
    # Adjust z by z0
    z -= z_0
    
    # Compute distances
    r_3d = np.sqrt(x**2 + y**2 + z**2)
    r = np.sqrt(x**2 + y**2)
    
    # Initialize vectors
    e_phi = np.zeros(3, dtype=np.float64)
    e_rho = np.zeros(3, dtype=np.float64)
    
    # Tangential and radial coordinates
    if r/sphere_radius < 1e-12:
        e_phi[:] = 0.0
        e_rho[:] = 0.0
    else:
        e_phi[0] = -y / r
        e_phi[1] = x / r
        e_rho[0] = x / r
        e_rho[1] = y / r
    
    if r_3d < sphere_radius:
        # Inside the sphere
        u_z = 3.0/5.0 * u_freestream * (1.0 - (2.0 * r**2 + z**2) / sphere_radius**2) + 2.0/5.0 * u_freestream
        u_rho = 3.0/5.0 * u_freestream * (r * z) / sphere_radius**2
        omega_phi = 3.0 * u_freestream * r / sphere_radius**2
        deformation_phi = 9.0/5.0 * u_freestream**2 / sphere_radius**2 * r * z
    else:
        # Outside the sphere
        u_z = 2.0/5.0 * u_freestream * (sphere_radius**2 / (z**2 + r**2))**(5.0/2.0) * (2.0 * z**2 - r**2) / (2.0 * sphere_radius**2)
        u_rho = 3.0/5.0 * u_freestream * r * z / sphere_radius**2 * (sphere_radius**2 / (z**2 + r**2))**(5.0/2.0)
        omega_phi = 0.0
        deformation_phi = 0.0
    
    # Induced velocity
    u_induced = np.zeros(3, dtype=np.float64)
    u_induced[0] = u_rho * e_rho[0]
    u_induced[1] = u_rho * e_rho[1]
    u_induced[2] = u_z
    
    # Deformation
    deformation = np.zeros(3, dtype=np.float64)
    deformation[0] = deformation_phi * e_phi[0]
    deformation[1] = deformation_phi * e_phi[1]
    deformation[2] = 0.0
    
    # Vorticity
    vorticity = np.zeros(3, dtype=np.float64)
    vorticity[0] = omega_phi * e_phi[0]
    vorticity[1] = omega_phi * e_phi[1]
    vorticity[2] = 0.0

    return u_induced, deformation, vorticity

def hill_assign(NN, NN_bl, Xbound, Dpm, neqpm, sphere_radius=1.0, u_freestream=-1.0, sphere_z_center=0.0):
    RHS_pm_bl = np.zeros((
        neqpm,
        NN_bl[3] - NN_bl[0] + 1, 
        NN_bl[4] - NN_bl[1] + 1, 
        NN_bl[5] - NN_bl[2] + 1
    ), order='F', dtype=float)
    # The RHS_pm array contains the right-hand side of the Poisson equation
    # Meaning it contains the negative vorticity values 
    analytic_sol = np.zeros((6, NN[0], NN[1], NN[2]))
    # Analytic solution for the Hill's spherical vortex contains:
    # 3 velocity components (u, v, w) and 3 deformation components (du/dx, dv/dy, dw/dz)

    # Unpack NN_bl for clarity
    i_start, j_start, k_start = NN_bl[0], NN_bl[1], NN_bl[2]
    i_end, j_end, k_end = NN_bl[3], NN_bl[4], NN_bl[5]

    # Create range arrays for i, j, k
    i_range = np.arange(i_start-1, i_end - 1)
    j_range = np.arange(j_start-1, j_end - 1)
    k_range = np.arange(k_start-1, k_end - 1)

    # Create a meshgrid to replace the loops
    i_grid, j_grid, k_grid = np.meshgrid(i_range, j_range, k_range, indexing='ij')

    # Vectorized computation of CP (Coordinates of Points)
    CP = np.stack([
        Xbound[0] + (i_grid) * Dpm[0],
        Xbound[1] + (j_grid) * Dpm[1],
        Xbound[2] + (k_grid) * Dpm[2]
    ], axis=-1)


    with open('results/hill_spherical_vortex.dat', 'w') as f:
        f.write('Variables = "x", "y", "z", "u", "v", "w", "omega_x", "omega_y", "omega_z"\n')
        f.write(f'Zone T="Hill\'s Spherical Vortex", I={NN[0]-1}, J={NN[1]-1}, K={NN[2]-1}, F=POINT\n')
        for k in range(k_start -1, k_end -1):
            for j in range(j_start -1, j_end -1):
                for i in range(i_start -1, i_end -1):
                        # Call fUi_HillsVortex_1 function
                        Uind, Defm, Vort = hills_vortex(
                            CP[i,j,k,:], sphere_radius, u_freestream, sphere_z_center
                        )
                        # Update RHS_pm_bl with negative vorticity values and analytic_sol with the computed values
                        RHS_pm_bl[:3, i, j, k] = -Vort[:3]
                        analytic_sol[:3, i, j, k] = Uind
                        analytic_sol[3:6, i, j, k] = Defm
                        # Write the data to the file
                        f.write(
                            f'{CP[i,j,k,0]:.8e} {CP[i,j,k,1]:.8e} {CP[i,j,k,2]:.8e} {Uind[0]:.8e} {Uind[1]:.8e} {Uind[2]:.8e} {-Vort[0]:.8e} {-Vort[1]:.8e} {-Vort[2]:.8e}\n'
                        )
    return analytic_sol, RHS_pm_bl

def hill_assign_parallel(NN, NN_bl, Xbound, Dpm, neqpm, sphere_radius=1.0, u_freestream=-1.0, sphere_z_center=0.0):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    RHS_pm_bl = np.zeros((
        neqpm,
        NN_bl[3] - NN_bl[0] + 1, 
        NN_bl[4] - NN_bl[1] + 1, 
        NN_bl[5] - NN_bl[2] + 1
    ), order='F', dtype=float)

    analytic_sol = np.zeros((6, NN[0], NN[1], NN[2]))

    i_start, j_start, k_start = NN_bl[0], NN_bl[1], NN_bl[2]
    i_end, j_end, k_end = NN_bl[3], NN_bl[4], NN_bl[5]

    # Distribute k-range among processors
    k_range = np.array_split(range(k_start - 1, k_end - 1), size)[rank]

    local_RHS_pm_bl = np.zeros_like(RHS_pm_bl)
    local_analytic_sol = np.zeros_like(analytic_sol)

    CP = np.zeros((i_end - i_start + 1, j_end - j_start + 1, len(k_range), 3))
    for idx, k in enumerate(k_range):
        for j in range(j_start - 1, j_end - 1):
            for i in range(i_start - 1, i_end - 1):
                CP[i-(i_start-1), j-(j_start-1), idx] = [
                    Xbound[0] + i * Dpm[0],
                    Xbound[1] + j * Dpm[1],
                    Xbound[2] + k * Dpm[2]
                ]

    # if rank == 0:
    #     file = open('results/hill_spherical_vortex.dat', 'w')
    #     file.write(f'Variables = "x", "y", "z", "u", "v", "w", "omega_x", "omega_y", "omega_z"\n')
    #     file.write(f'Zone T="Hill\'s Spherical Vortex", I={NN[0]-1}, J={NN[1]-1}, K={NN[2]-1}, F=POINT\n')
    # else:
    #     file = None

    for idx, k in enumerate(k_range):
        for j in range(j_start - 1, j_end - 1):
            for i in range(i_start - 1, i_end - 1):
                Uind, Defm, Vort = hills_vortex(
                    CP[i-(i_start-1), j-(j_start-1), idx], sphere_radius, u_freestream, sphere_z_center
                )
                local_RHS_pm_bl[:3, i-(i_start-1), j-(j_start-1), idx] = -Vort[:3]
                local_analytic_sol[:3, i, j, k] = Uind
                local_analytic_sol[3:6, i, j, k] = Defm

                # if file:
                #     file.write(
                #         f'{CP[i-(i_start-1),j-(j_start-1),idx,0]:.8e} {CP[i-(i_start-1),j-(j_start-1),idx,1]:.8e} '
                #         f'{CP[i-(i_start-1),j-(j_start-1),idx,2]:.8e} {Uind[0]:.8e} {Uind[1]:.8e} {Uind[2]:.8e} '
                #         f'{-Vort[0]:.8e} {-Vort[1]:.8e} {-Vort[2]:.8e}\n'
                #     )

    # if file:
    #     file.close()

    # Gather results on process 0
    gathered_RHS_pm_bl = comm.gather(local_RHS_pm_bl, root=0)
    gathered_analytic_sol = comm.gather(local_analytic_sol, root=0)

    if rank == 0:
        for i, part in enumerate(gathered_RHS_pm_bl):
            k_slice = np.array_split(range(k_start - 1, k_end - 1), size)[i]
            RHS_pm_bl[:, :, :, k_slice - (k_start - 1)] = part[:, :, :, :len(k_slice)]

        for i, part in enumerate(gathered_analytic_sol):
            k_slice = np.array_split(range(k_start - 1, k_end - 1), size)[i]
            analytic_sol[:, :, :, k_slice] = part[:, :, :, k_slice]

    # Broadcast final results to all processes
    RHS_pm_bl = comm.bcast(RHS_pm_bl, root=0)
    analytic_sol = comm.bcast(analytic_sol, root=0)

    return analytic_sol, RHS_pm_bl

def hill_error(NN, NN_bl, Xbound, Dpm, SOL_pm, vel_pm, analytic_sol):
    # Initialize the error array
    error = np.zeros((7, NN[0], NN[1], NN[2]), dtype=float)
    velvrx_pm = vel_pm[0]
    velvry_pm = vel_pm[1]
    velvrz_pm = vel_pm[2]

    max_err = np.zeros(7, dtype=float)
    mean_err = np.zeros(7, dtype=float)

    # Main computation loop
    for k in range(NN_bl[2] + 1, NN_bl[5]):
        for j in range(NN_bl[1] + 1, NN_bl[4]):
            for i in range(NN_bl[0] + 1, NN_bl[3]):
                CP = np.array([Xbound[0] + (i - 1) * Dpm[0], 
                               Xbound[1] + (j - 1) * Dpm[1], 
                               Xbound[2] + (k - 1) * Dpm[2]], dtype=float)

                error[0, i, j, k] = abs(velvrx_pm[i, j, k] - analytic_sol[0, i, j, k])
                error[1, i, j, k] = abs(velvry_pm[i, j, k] - analytic_sol[1, i, j, k])
                error[2, i, j, k] = abs(velvrz_pm[i, j, k] - analytic_sol[2, i, j, k])

                error[3, i, j, k] = abs(SOL_pm[0, i, j, k] - analytic_sol[3, i, j, k])
                error[4, i, j, k] = abs(SOL_pm[1, i, j, k] - analytic_sol[4, i, j, k])
                error[5, i, j, k] = abs(SOL_pm[2, i, j, k] - analytic_sol[5, i, j, k])

                vmag_num = np.sqrt(velvrx_pm[i, j, k]**2 + velvry_pm[i, j, k]**2 + velvrz_pm[i, j, k]**2)
                vmag_anal = np.sqrt(analytic_sol[0, i, j, k]**2 + analytic_sol[1, i, j, k]**2 + analytic_sol[2, i, j, k]**2)
                error[6, i, j, k] = abs(vmag_num - vmag_anal)

                max_err = np.maximum(max_err, error[:, i, j, k])
                mean_err += error[:, i, j, k]
    mean_err /= (NN[0] * NN[1] * NN[2])
    print(f'----Maximum Velocity Error-----: {max_err[6] * 100:.2f}%')
    print(f'----Mean Velocity Error-----: {mean_err[6] * 100:.2f}%')

def visualize_vorticity(RHS_pm_bl, NN_bl):
    """
    Plots the vorticity components (vorticity_x, vorticity_y, vorticity_z) from RHS_pm_bl in 3D.

    Parameters:
    RHS_pm_bl (np.ndarray): Array containing the vorticity data with shape (3, nx, ny, nz).
                            The first dimension corresponds to the 3 vorticity components (x, y, z).
    NN_bl (list): List containing the start and end indices for the 3D domain (ix, iy, iz).
    """
    # Get the indices for the 3D grid
    i_start, j_start, k_start = NN_bl[0], NN_bl[1], NN_bl[2]
    i_end, j_end, k_end = NN_bl[3], NN_bl[4], NN_bl[5]

    # Create meshgrid for i, j, k
    i_range = np.arange(i_start, i_end )
    j_range = np.arange(j_start, j_end )
    k_range = np.arange(k_start, k_end )

    i_grid, j_grid, k_grid = np.meshgrid(i_range, j_range, k_range, indexing='ij')

    # Flatten the grids for plotting
    i_flat = i_grid.flatten()
    j_flat = j_grid.flatten()
    k_flat = k_grid.flatten()

    # Flatten the vorticity components for each x, y, z
    vorticity_x = RHS_pm_bl[0].flatten()
    vorticity_y = RHS_pm_bl[1].flatten()
    vorticity_z = RHS_pm_bl[2].flatten()
    vorticity_mag = vorticity_x**2 + vorticity_y**2 + vorticity_z**2

    # Create 3D scatter plots for each vorticity component
    fig = plt.figure(figsize=(15, 5))


    # # Plot vorticity_x
    # ax1 = fig.add_subplot(131, projection='3d')
    # ax1.scatter(i_flat, j_flat, k_flat, c=vorticity_x, cmap='coolwarm', marker='o')
    # ax1.set_title("Vorticity X")
    # ax1.set_xlabel("i")
    # ax1.set_ylabel("j")
    # ax1.set_zlabel("k")

    # # Plot vorticity_y
    # ax2 = fig.add_subplot(132, projection='3d')
    # ax2.scatter(i_flat, j_flat, k_flat, c=vorticity_y, cmap='coolwarm', marker='o')
    # ax2.set_title("Vorticity Y")
    # ax2.set_xlabel("i")
    # ax2.set_ylabel("j")
    # ax2.set_zlabel("k")

    # Plot vorticity_z

    # Create a mask of all points where vorticity_mag > 1e-6
    mask = vorticity_mag > 1e-12
    ax3 = fig.add_subplot(111, projection='3d')
    sc = ax3.scatter(i_flat[mask], j_flat[mask], k_flat[mask], c=vorticity_mag[mask], cmap='viridis', marker='o')
    ax3.set_title("Vorticity Magnitude")
    ax3.set_xlabel("i")
    ax3.set_ylabel("j")
    ax3.set_zlabel("k")
    # Add colorbar
    fig.colorbar(sc, ax=ax3, orientation='vertical', label='Vorticity Magnitude') 

    # Adjust layout
    plt.tight_layout()
    plt.show()
