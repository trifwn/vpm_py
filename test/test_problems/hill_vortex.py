import numpy as np
from vpm_py import VPM

def hill_vortex(
    control_points, 
    sphere_radius, 
    U, 
    density = 1.225,
    z_0 = 0
):    
    """
    Args:
        CP (_type_): Control Point
        a: The radius of the sphere 
        U: The constant free stream velocity far away from the sphere
        z0: The z-coordinate of the sphere center

    Returns:
        _type_: _description_
    """
    density = 1.225
    x = control_points[:, :, :, 0]
    y = control_points[:, :, :, 1]
    z = control_points[:, :, :, 2] - z_0
    
    # Initialize velocity and vorticity arrays
    velocity = np.zeros_like(control_points)
    vorticity = np.zeros_like(control_points)
    pressure = np.zeros_like(x)
    
    # Convert to spherical coordinates
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)
    eta = np.arctan2(np.sqrt(x**2 + y**2), z)

    inside_mask = r < sphere_radius
    outside_mask = r >= sphere_radius
    a = sphere_radius

    # Allocate the arrays
    u_eta = np.zeros_like(r)
    u_phi = np.zeros_like(r)
    u_z = np.zeros_like(r)
    omega_phi = np.zeros_like(r)
    pressure = np.zeros_like(r)

    eta_inside = eta[inside_mask]
    z_inside = z[inside_mask]
    eta_outside = eta[outside_mask]
    z_outside = z[outside_mask]

    # Velocity inside and outside the sphere
    u_eta[inside_mask] = 3 * U / (2 * a**2) * eta_inside * z_inside
    u_phi[inside_mask] = 0.0
    u_z[inside_mask] = 3 * U / (2 * a**2) *  (a**2  - 2 * eta_inside**2 - z_inside**2)

    u_eta[outside_mask] = (3 * a **3 * U ) / 2 * (eta_outside * z_outside) / (eta_outside**2 + z_outside**2)**(5/2)
    u_phi[outside_mask] = 0.0
    u_z[outside_mask] = -U - (a **3 * U) / 2 * (eta_outside **2 - 2 * z_outside ** 2) / (eta_outside ** 2 + z_outside ** 2) ** (5/2)

    # Vorticity inside and outside the sphere
    omega_phi[inside_mask] = 3 * U * np.sin(eta_inside) / (2 * a) * (1/2 + 2 * (r[inside_mask]/a)**2)
    omega_phi[outside_mask] = 0.

    # Pressure inside and outside the sphere
    pressure[inside_mask] = - (9 * U**2 * density) / (8 * a**4) * (
                        - eta_inside **4 + z_inside**4 + (a**2) *(eta_inside**2 - 2 * z_inside**2)
                    )
    pressure[outside_mask] = U**2 * density / 8 * (
        5 - (4 * a**3 * (eta_outside**2 - 2* z_outside**2) ) / (eta_outside**2 + z_outside**2)**(5/2)
        - a**6 * (eta_outside**2 + 4 * z_outside**2) / (eta_outside**2 + z_outside**2)**(4)
    )

    # Unit vectors (Converting to Cartesian coordinates)
    velocity[:,:,:,0] = u_eta * np.cos(phi) * np.cos(eta) - u_phi * np.sin(phi) - u_z * np.sin(eta) * np.cos(phi)
    velocity[:,:,:,1] = u_eta * np.sin(phi) * np.cos(eta) + u_phi * np.cos(phi) - u_z * np.sin(eta) * np.sin(phi)
    velocity[:,:,:,2] = u_eta * np.sin(eta) + u_z * np.cos(eta)
    
    
    vorticity[:,:,:,0] = - omega_phi * np.sin(phi)
    vorticity[:,:,:,1] = omega_phi * np.cos(phi)
    vorticity[:,:,:,2] = 0.0
    # NaN to zero
    velocity = np.nan_to_num(velocity)
    vorticity = np.nan_to_num(vorticity)

    return velocity, vorticity, pressure

def hill_assign(
    vpm: VPM, 
    sphere_radius: float =1.0, 
    u_freestream: float =-1.0, 
    sphere_z_center: float =0.0,
    density: float = 1.225
):

    print("Getting analytical solution")
    XYZ = vpm.particle_mesh.grid_positions
    X = XYZ[0, :, :, :]
    Y = XYZ[1, :, :, :]
    Z = XYZ[2, :, :, :]

    CP = np.stack([X, Y, Z], axis=-1)
    shape_CP = CP.shape
    print(f"CP shape: {shape_CP}")
    
    analytical_velocity, analytical_vorticity, analytical_pressure = hill_vortex(
        control_points= CP, 
        sphere_radius= sphere_radius,
        U= u_freestream,
        z_0= sphere_z_center,
        density= density
    )
    # PRINT VOLOCITY
    print(f"Analytical vorticity: {analytical_vorticity.shape}")

    # Expand the arrays to the original shape
    analytical_velocity = analytical_velocity.reshape(shape_CP)
    analytical_vorticity = analytical_vorticity.reshape(shape_CP)
    analytical_pressure = analytical_pressure.reshape(shape_CP[:3])

    # Move axis
    analytical_vorticity = np.moveaxis(analytical_vorticity, [0, 1, 2, 3], [1, 2, 3, 0])
    # analytical_velocity = np.moveaxis(analytical_velocity, [0, 1, 2, 3], [3, 2, 1, 0])
    # analytical_pressure = np.moveaxis(analytical_pressure, [0, 1, 2], [2, 1, 0])

    analytical_vorticity = np.asfortranarray(analytical_vorticity)
    # analytical_pressure = np.asfortranarray(analytical_pressure)
    # analytical_velocity = np.asfortranarray(analytical_velocity)
    return analytical_velocity, analytical_vorticity, analytical_pressure
