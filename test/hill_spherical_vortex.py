import numpy as np

def fUi_HillsVortex_1(CP, a, us, z0):    
    # Unpacking control point
    x, y, z = CP
    
    # Adjust z by z0
    z -= z0
    
    # Compute distances
    r3d = np.sqrt(x**2 + y**2 + z**2)
    rho = np.sqrt(x**2 + y**2)
    
    # Initialize vectors
    e_phi = np.zeros(3, dtype=np.float64)
    e_rho = np.zeros(3, dtype=np.float64)
    
    # Tangential and radial coordinates
    if rho/a < 1e-12:
        e_phi[:] = 0.0
        e_rho[:] = 0.0
    else:
        e_phi[0] = -y / rho
        e_phi[1] = x / rho
        e_rho[0] = x / rho
        e_rho[1] = y / rho
    
    if r3d < a:
        # Inside the sphere
        uz = 3.0/5.0 * us * (1.0 - (2.0 * rho**2 + z**2) / a**2) + 2.0/5.0 * us
        urho = 3.0/5.0 * us * (rho * z) / a**2
        om_phi = 3.0 * us * rho / a**2
        defm_phi = 9.0/5.0 * us**2 / a**2 * rho * z
        Grad = np.zeros(9, dtype=np.float64)
    else:
        # Outside the sphere
        uz = 2.0/5.0 * us * (a**2 / (z**2 + rho**2))**(5.0/2.0) * (2.0 * z**2 - rho**2) / (2.0 * a**2)
        urho = 3.0/5.0 * us * rho * z / a**2 * (a**2 / (z**2 + rho**2))**(5.0/2.0)
        om_phi = 0.0
        defm_phi = 0.0
        Grad = np.zeros(9, dtype=np.float64)
    
    # Induced velocity
    Uind = np.zeros(3, dtype=np.float64)
    Uind[0] = urho * e_rho[0]
    Uind[1] = urho * e_rho[1]
    Uind[2] = uz
    
    # Deformation
    Defm = np.zeros(3, dtype=np.float64)
    Defm[0] = defm_phi * e_phi[0]
    Defm[1] = defm_phi * e_phi[1]
    
    # Vorticity
    Vort = np.zeros(3, dtype=np.float64)
    Vort[0] = om_phi * e_phi[0]
    Vort[1] = om_phi * e_phi[1]
    
    return Uind, Grad, Defm, Vort

def hill_assign(NN, NN_bl, Xbound, Dpm, RHS_pm_bl, neqpm):
    # Initialize RHS array
    RHS_pm_bl[:3, :, :, :] = 0.0
    # Allocate analytic_sol
    analytic_sol = np.zeros((6, NN[0], NN[1], NN[2]))

    # Main computation loop
    for k in range(NN_bl[2], NN_bl[5]):
        for j in range(NN_bl[1], NN_bl[4]):
            for i in range(NN_bl[0], NN_bl[3]):
                CP = np.array([Xbound[0] + (i - 1) * Dpm[0], 
                               Xbound[1] + (j - 1) * Dpm[1], 
                               Xbound[2] + (k - 1) * Dpm[2]], dtype=float)
                
                Uind = np.zeros(3, dtype=float)
                Grad = np.zeros(9, dtype=float)
                Defm = np.zeros(3, dtype=float)
                Vort = np.zeros(3, dtype=float)

                # Call fUi_HillsVortex_1 function (assumed to be defined similarly in Python)
                Uind, Grad, Defm, Vort = fUi_HillsVortex_1(CP, 1.0, -1.0, 0.0)

                # print(i,j,k)
                # Update RHS_pm_bl with negative vorticity values
                RHS_pm_bl[:3, i, j, k] = -Vort[:3]

                # Store analytical solution
                analytic_sol[:3, i, j, k] = Uind
                analytic_sol[3:6, i, j, k] = Defm

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

# Example Usage
def main():
    import matplotlib.pyplot as plt
    sizing = 60
    NN = [sizing, sizing, sizing]  # Grid size
    NN_bl = [0,0,0, sizing, sizing, sizing]  # Block limits
    Xbound = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]  # Bounds in each direction
    Dpm = [(Xbound[3] - Xbound[0]) / (NN[0] - 1),
        (Xbound[4] - Xbound[1]) / (NN[1] - 1),
        (Xbound[5] - Xbound[2]) / (NN[2] - 1)]  # Grid spacing
    neqpm = 3  # Number of equations per grid point

    # Initialize the RHS_pm_bl array
    RHS_pm_bl = np.zeros((neqpm, NN[0], NN[1], NN[2]), dtype=float)

    # Call hill_assign to compute the fields
    hill_assign(NN, NN_bl, Xbound, Dpm, RHS_pm_bl, neqpm)
    print("Hill vortex computation complete")
    # Extract the computed vorticity field from RHS_pm_bl
    vorticity_field = -RHS_pm_bl[:3]

    # hill_error(NN, NN_bl, Xbound, Dpm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm, analytic_sol)


    # Plotting the vorticity field
    # We'll plot a slice of the vorticity field at the mid-plane (z=0.5)
    for mid_index in [int(NN[2] / 2), int(NN[2] / 2) + 10]: 
        X, Y = np.meshgrid(np.linspace(Xbound[0], Xbound[3], NN[0]), np.linspace(Xbound[1], Xbound[4], NN[1]))
        fig, ax = plt.subplots(1, 3, figsize=(18, 6))

        # Vorticity components
        titles = ['Vorticity X', 'Vorticity Y', 'Vorticity Z']
        for i in range(3):
            # ax[i].contourf(X, Y, vorticity_field[i, :, :, mid_index], cmap='jet', levels=100)
            ax[i].matshow(vorticity_field[i, :, :, mid_index], cmap='turbo')
            #  Scatter the grid points
            ax[i].set_title(titles[i])
            ax[i].set_xlabel('X (m)')
            ax[i].set_ylabel('Y (m)')

        plt.show()

if __name__ == "__main__":
    main()