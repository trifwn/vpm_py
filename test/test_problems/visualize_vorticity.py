import numpy as np
import matplotlib.pyplot as plt
from vpm_py import VPM

def visualize_vorticity(
    vpm: VPM, 
    vorticity_pm: np.ndarray,
):
    """
    Plots the vorticity components (vorticity_x, vorticity_y, vorticity_z) from RHS_pm_bl in 3D.

    Parameters:
    RHS_pm_bl (np.ndarray): Array containing the vorticity data with shape (3, nx, ny, nz).
                            The first dimension corresponds to the 3 vorticity components (x, y, z).
    NN_bl (list): List containing the start and end indices for the 3D domain (ix, iy, iz).
    """

    Dpm= vpm.dpm,
    NN= vpm.particle_mesh.grid_size,
    Xbound= vpm.particle_mesh.xbound,
    
    nx = NN[0]
    ny = NN[1]
    nz = NN[2]

    xs = Xbound[0] + np.arange(nx) * Dpm[0]
    ys = Xbound[1] + np.arange(ny) * Dpm[1]
    zs = Xbound[2] + np.arange(nz) * Dpm[2]


    xs = Xbound[0] + np.arange(nx) * Dpm[0]
    ys = Xbound[1] + np.arange(ny) * Dpm[1]
    zs = Xbound[2] + np.arange(nz) * Dpm[2]
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

    # Get the indices for the 3D grid
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()

    # Flatten the vorticity components for each x, y, z
    vorticity_x = vorticity_pm[0, :, :, :].flatten()
    vorticity_y = vorticity_pm[1, :, :, :].flatten()
    vorticity_z = vorticity_pm[2, :, :, :].flatten()
    vorticity_mag = vorticity_x**2 + vorticity_y**2 + vorticity_z**2

    # Create 3D scatter plots for each vorticity component
    fig = plt.figure(figsize=(15, 5))
    
    # Create a mask of all points where vorticity_mag > 1e-6
    mask = vorticity_mag > 1e-12
    ax3 = fig.add_subplot(111, projection='3d')
    sc = ax3.scatter(X[mask], Y[mask], Z[mask], c=vorticity_mag[mask], cmap='viridis', marker='o')
    ax3.set_title("Vorticity Magnitude")
    ax3.set_xlabel("i")
    ax3.set_ylabel("j")
    ax3.set_zlabel("k")
    # Add colorbar
    fig.colorbar(sc, ax=ax3, orientation='vertical', label='Vorticity Magnitude') 

    # Adjust layout
    plt.tight_layout()
    plt.show()
