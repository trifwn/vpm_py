
import numpy as np

def definevort(RHS_pm, Xbound, Dpm, NN, NN_bl):
    # Parameters for the vorticity calculation
    xc1, xc2 = 0.0, 0.0
    yc1, yc2 = -1.5, 1.5
    PI = 4 * np.arctan(1.0)
    
    # Allocate and initialize analytic_sol
    analytic_sol = np.zeros((1, NN[0], NN[1], NN[2]), dtype=float)
    
    # Loop through the grid
    for j in range(NN_bl[1], NN_bl[4]):
        for i in range(NN_bl[0], NN_bl[3]):
            xi = Xbound[0] + (i - 1) * Dpm[0]
            yi = Xbound[1] + (j - 1) * Dpm[1]
            ksi1 = np.sqrt((xi - xc1)**2 + (yi - yc1)**2)
            ksi2 = np.sqrt((xi - xc2)**2 + (yi - yc2)**2)
            th1 = np.arctan2(yi - yc1, xi - xc1)
            th2 = np.arctan2(yi - yc2, xi - xc2)
            if th1 < 0:
                th1 += 2 * PI
            if th2 < 0:
                th2 += 2 * PI

            w1 = (2.0 - ksi1**2) * np.exp(0.5 * (1.0 - ksi1**2))
            w2 = -(2.0 - ksi2**2) * np.exp(0.5 * (1.0 - ksi2**2))
            RHS_pm[0, i, j, 0] = -(w1 + w2)
            analytic_sol[0, i, j, 0] = np.exp(0.5 * (1.0 - ksi1**2)) - np.exp(0.5 * (1.0 - ksi2**2))
            
def vort_error(NN, NN_bl, Xbound, Dpm, SOL_pm, analytic_sol):
    # Parameters for error calculation
    analytic_sol = np.zeros((1, NN[0], NN[1], NN[2]), dtype=float)  # Assuming analytic_sol is available globally
    error = np.zeros((1, NN[0], NN[1], NN[2]), dtype=float)
    
    # Calculate the maximum value of the analytical solution for normalization
    analytic_max = np.max(np.abs(analytic_sol[0, :, :, :]))
    
    max_err = 0.0
    for j in range(NN_bl[1], NN_bl[4] + 1):
        for i in range(NN_bl[0], NN_bl[3] + 1):
            xi = Xbound[0] + (i - 1) * Dpm[0]
            yi = Xbound[1] + (j - 1) * Dpm[1]
            # Assume CP(3) is zero for 2D case
            CP = [xi, yi, 0]
            
            error_val = np.abs(SOL_pm[0, i, j, 0] - analytic_sol[0, i, j, 0])
            error[0, i, j, 0] = error_val
            
            max_err = max(max_err, error_val)
            mean_err += error_val

    mean_err /= (NN_bl[3] - NN_bl[0] + 1) * (NN_bl[4] - NN_bl[1] + 1)
    
    # Output the maximum error percentage
    print(f'----Maximum Phi Error-----: {max_err / analytic_max * 100:.2f}%')
    print(f"----Mean Phi Error-----: {mean_err * 100:.2f}%")

# Example Usage
def main():
    import matplotlib.pyplot as plt
    sizing = 60
    NN = [sizing, sizing, sizing]  # Grid size
    NN_bl = [0,0,0, sizing, sizing, sizing]  # Block limits
    Xbound = [-1, -1, -1, 1.0, 1.0, 1.0]  # Bounds in each direction
    Dpm = [(Xbound[3] - Xbound[0]) / (NN[0] - 1),
        (Xbound[4] - Xbound[1]) / (NN[1] - 1),
        (Xbound[5] - Xbound[2]) / (NN[2] - 1)]  # Grid spacing
    neqpm = 3  # Number of equations per grid point

    # Initialize the RHS_pm_bl array
    RHS_pm_bl = np.zeros((neqpm, NN[0], NN[1], NN[2]), dtype=float)

    # Call hill_assign to compute the fields
    definevort( RHS_pm_bl, Xbound, Dpm, NN, NN_bl)
    print("Hill vortex computation complete")
    # Extract the computed vorticity field from RHS_pm_bl
    vorticity_field = -RHS_pm_bl[:3]

    # vort_error(NN, NN_bl, Xbound, Dpm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm, analytic_sol)


    # Plotting the vorticity field
    # We'll plot a slice of the vorticity field at the mid-plane (z=0.5)
    for mid_index in [0]: 
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