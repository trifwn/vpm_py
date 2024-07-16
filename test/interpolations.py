import numpy as np
import matplotlib.pyplot as plt

def projection_fun(itype, x):
    xabs = np.abs(x)
    
    if itype == 2:
        return np.where(xabs <= 1, 1 - xabs, 0)
    
    elif itype == 3:
        return np.where(xabs > 1.5, 0,
               np.where(xabs >= 0.5, 0.5 * (1.5 - xabs)**2,
               0.5 * (xabs + 1.5)**2 - 1.5 * (xabs + 0.5)**2))
    
    elif itype == 4:
        return np.where(xabs > 2, 0,
               np.where(xabs >= 1, 0.5 * (2 - xabs)**2 * (1 - xabs),
               1 - 2.5*xabs**2 + 1.5*xabs**3))
    
    elif itype == 5:
        return np.where(xabs > 1.5, 0,
               np.where(xabs > 0.5, 0.5 * (1 - xabs) * (2 - xabs),
               1 - xabs**2))
    
    elif itype == 6:
        return np.where(xabs > 2, 0,
               np.where(xabs > 1, 0.5 * (1 - xabs) * (2 - xabs),
               1 - xabs**2))
    
    else:
        raise ValueError(f"No such projection function: {itype}")

# Create x values
x = np.linspace(-3, 3, 1000)

# Set up the plot
plt.figure(figsize=(12, 8))

# Plot each function
for itype in range(2, 7):
    y = projection_fun(itype, x)
    plt.plot(x, y, label=f'Type {itype}')

plt.title('Projection Functions')
plt.xlabel('x')
plt.ylabel('W(x)')
plt.legend()
plt.grid(True)
plt.axhline(y=0, color='k', linestyle='--')
plt.axvline(x=0, color='k', linestyle='--')

plt.show()