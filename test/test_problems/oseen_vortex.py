import numpy as np
from vpm_py import VPM
from scipy.special import gammainc

def oseen_vortex(
    control_points: np.ndarray, 
    G: float, 
    viscosity: float,
    density: float,
    t: float
):  
    x = control_points[:, :, :, 0]
    y = control_points[:, :, :, 1]
    
    # Initialize velocity and vorticity arrays
    velocity = np.zeros_like(control_points)
    vorticity = np.zeros_like(control_points)
    pressure = np.zeros_like(x)
    
    # Convert to spherical coordinates
    r = np.sqrt(x**2 + y**2 )
    theta = np.arctan2(y, x)

    exponent = r**2 / (4 * viscosity * t)
    D = G / (2 * np.pi)

    u_theta = (D/r) *(1 - np.exp(-exponent))

    r_safe = np.where(r == 0, 1e-10, r)
    
    # Compute the pressure field
    pressure = density * D**2 * (
        1/(2 * r_safe**2) * (1 - np.exp(-exponent))**2 
        - 1/(4 * viscosity * t) * (
            gammainc(0, exponent) - gammainc(0, 2 * exponent)
        )
    )

    # Calculate velocity and vorticity
    velocity[:, :, :, 0] = u_theta * np.sin(theta) 
    velocity[:, :, :, 1] = -u_theta * np.cos(theta)
    velocity[:, :, :, 2] = 0

    vorticity[:, :, :, 0] = 0
    vorticity[:, :, :, 1] = 0
    vorticity[:, :, :, 2] = G / (4 * np.pi * viscosity * t) * np.exp(-exponent) 
    return velocity, vorticity, pressure

def oseen_assign(
    vpm: VPM, 
    viscosity: float,
    density: float,
    t: float,
    gamma: float, 
):

    XYZ = vpm.particle_mesh.grid_positions
    X = XYZ[0, :, :, :]
    Y = XYZ[1, :, :, :]
    Z = XYZ[2, :, :, :]
 
    CP = np.stack([X, Y, Z], axis=-1)
    shape_CP = CP.shape

    analytical_velocity, analytical_vorticity, analytical_pressure = oseen_vortex(
        control_points= CP, 
        G= gamma,
        viscosity= viscosity, 
        density= density,
        t= t
    )
    # Expand the arrays to the original shape
    analytical_velocity = analytical_velocity.reshape(shape_CP)
    analytical_vorticity = analytical_vorticity.reshape(shape_CP)
    analytical_pressure = analytical_pressure.reshape(shape_CP[:3])

    # Move axis
    analytical_vorticity = np.moveaxis(analytical_vorticity, [0, 1, 2, 3], [1, 2, 3, 0])
    analytical_velocity = np.moveaxis(analytical_velocity, [0, 1, 2, 3], [1, 2, 3, 0])

    analytical_vorticity = np.asfortranarray(analytical_vorticity)
    analytical_velocity = np.asfortranarray(analytical_velocity)
    analytical_pressure = np.asfortranarray(analytical_pressure)
    return analytical_velocity, analytical_vorticity, analytical_pressure
