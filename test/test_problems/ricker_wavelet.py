import numpy as np

def ricker_wavelet_3d(control_point, sphere_radius):
    """
    Compute the 3D Ricker wavelet scalar field at a given control point.

    Args:
        control_point (tuple): (x, y, z) coordinates of the control point.
        sphere_radius (float): The characteristic radius of the wavelet.

    Returns:
        float: The Ricker wavelet value at the control point.
    """
    # Unpack control point
    x, y, z = control_point
    # Compute 3D radial distance squared
    r_squared = x**2 + y**2 + z**2
    sigma_squared = sphere_radius**2
    sigma_4 = sphere_radius**4

    # Compute the Ricker wavelet value
    if r_squared / sigma_squared <= 1:
        wavelet_value = (1/(np.pi * sigma_4)) * (1 - 0.5 * r_squared / sigma_squared) * np.exp(-r_squared / (2 * sigma_squared))
    else:
        wavelet_value = 0.0

    return wavelet_value


def ricker_wavelet_3d_field(NN, Xbound, Dpm, sphere_radius=1.0):
    """
    Compute the Ricker wavelet scalar field over a 3D grid.

    Args:
        NN (tuple): Number of points in each dimension (nx, ny, nz).
        Xbound (tuple): Starting coordinates (x0, y0, z0).
        Dpm (tuple): Grid spacing in each dimension (dx, dy, dz).
        sphere_radius (float): The characteristic radius of the wavelet.

    Returns:
        ndarray: 3D array containing the wavelet values.
    """
    nx, ny, nz = NN
    x0, y0, z0, _, _ , _ = Xbound
    dx, dy, dz = Dpm

    # Initialize 3D array for the wavelet values
    wavelet_field = np.zeros((nx, ny, nz))

    # Iterate over the grid and compute the wavelet value for each point
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Compute coordinates of the control point
                x = x0 + i * dx
                y = y0 + j * dy
                z = z0 + k * dz

                # Compute the Ricker wavelet value
                wavelet_field[i, j, k] = ricker_wavelet_3d((x, y, z), sphere_radius)

    return wavelet_field


def ricker_wavelet_3d_with_gradient(NN, Xbound, Dpm, radius=1.0 ):
    """
    Compute the Ricker wavelet scalar field and gradient over a 3D grid.

    Args:
        NN (tuple): Number of points in each dimension (nx, ny, nz).
        Xbound (tuple): Starting coordinates (x0, y0, z0).
        Dpm (tuple): Grid spacing in each dimension (dx, dy, dz).
        sphere_radius (float): The characteristic radius of the wavelet.
        z_0 (float): The z-coordinate of the wavelet center.

    Returns:
        tuple: 
            - wavelet_field (ndarray): 3D array containing the wavelet values.
            - gradient_field (ndarray): 4D array containing the gradient values (x, y, z) at each point.
    """
    nx, ny, nz = NN
    x0, y0, z0, _, _, _ = Xbound
    dx, dy, dz = Dpm

    # Initialize 3D array for the wavelet values and a 4D array for gradients
    wavelet_field = np.zeros((nx, ny, nz))
    gradient_field = np.zeros((3, nx, ny, nz))

    # Iterate over the grid and compute the wavelet value and gradient for each point
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Compute coordinates of the control point
                x = x0 + i * dx
                y = y0 + j * dy
                z = z0 + k * dz

                r = np.array([x, y, z])
                # Compute the wavelet value
                wavelet_value = ricker_wavelet_3d((x, y, z), radius)
                wavelet_field[i, j, k] = wavelet_value

                # Compute gradient
                r_norm = np.linalg.norm(r)
                r_3d = r_norm if r_norm > 1e-12 else 1e-12

                if r_norm / radius <= 1:
                    grad = -wavelet_value * r / r_3d
                else:
                    grad = np.array([0.0, 0.0, 0.0])

                gradient_field[:, i, j, k] = grad

    return wavelet_field, gradient_field
