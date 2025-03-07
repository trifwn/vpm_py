import numpy as np

def ricker_wavelet_3d(control_point, sigma):
    """
    Compute the 3D Ricker wavelet scalar field at a given control point.

    Args:
        control_point (tuple): (x, y, z) coordinates of the control point.
        sigma (float): Standard deviation (characteristic radius) of the wavelet.

    Returns:
        float: The Ricker wavelet value at the control point.
    """
    x, y, z = control_point
    r_squared = x**2 + y**2 + z**2
    sigma_squared = sigma**2

    # Compute the Ricker wavelet value using the analytical expression:
    # ψ(r) = (3 - r^2/σ^2) (1/(2πσ^2)^(3/2)) exp(-r^2/(2σ^2))
    wavelet_value = (3 - r_squared/sigma_squared) * (1 / (2 * np.pi * sigma_squared)**(3/2)) \
                    * np.exp(-r_squared / (2 * sigma_squared))
    
    return wavelet_value

def ricker_wavelet_3d_with_gradient(NN, Xbound, Dpm, sigma=1.0, t=0.0, nu=1.0):
    """
    Compute the Ricker wavelet scalar field and gradient over a 3D grid
    with time evolution under diffusion.

    Args:
        NN (tuple): Number of points in each dimension (nx, ny, nz).
        Xbound (tuple): Starting coordinates (x0, y0, z0).
        Dpm (tuple): Grid spacing in each dimension (dx, dy, dz).
        sigma (float): The characteristic radius (standard deviation) of the wavelet.
        t (float): Time at which to evaluate the field.
        nu (float): Diffusivity coefficient.

    Returns:
        tuple: 
            - wavelet_field (ndarray): 3D array containing the wavelet values.
            - gradient_field (ndarray): 4D array containing the gradient values (x, y, z) at each point.
    """
    nx, ny, nz = NN
    x0, y0, z0 = Xbound[:3]
    dx, dy, dz = Dpm

    # Initialize 3D array for the wavelet values and a 4D array for gradients
    wavelet_field = np.zeros((nx, ny, nz))
    gradient_field = np.zeros((3, nx, ny, nz))

    # Time evolution parameter: tau = sigma^2 + 2*nu*t
    tau = sigma**2 + 2 * nu * t

    # Precompute the diffused Gaussian coefficient
    gaussian_coeff = 1 / (2 * np.pi * tau)**(3/2)

    # Iterate over the grid and compute the wavelet value and its gradient at each point
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Compute coordinates of the control point
                x = x0 + i * dx
                y = y0 + j * dy
                z = z0 + k * dz

                # Compute r^2 and r
                r_squared = x**2 + y**2 + z**2
                r = np.sqrt(r_squared)

                # Compute the scalar field (for t=0, this recovers the original Ricker wavelet)
                wavelet_field[i, j, k] = sigma**2 * (3 /tau - r_squared / tau**2) * gaussian_coeff \
                                         * np.exp(-r_squared / (2 * tau))
                
                # Compute the time-evolved Gaussian: G(r,t) = 1/[(2πτ)^(3/2)] exp(-r^2/(2τ))
                G_rt = gaussian_coeff * np.exp(-r_squared / (2 * tau))
                
                # Compute the gradient using the analytical expression:
                # ∇ψ(r,t) = -σ^2 * (r/τ^2) * (5 - r^2/τ) * G(r,t)  · (r̂)
                # Avoid division by zero by ensuring r is not exactly zero.
                r_unit = np.array([x, y, z]) / (r if r > 1e-12 else 1e-12)
                grad_magnitude = - sigma**2 * (r / tau**2) * (5 - r_squared / tau) * G_rt

                # Store the gradient (vector field is purely radial)
                gradient_field[:, i, j, k] = grad_magnitude * r_unit

    return wavelet_field, gradient_field
