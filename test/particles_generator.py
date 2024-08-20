from scipy.io import FortranFile
import numpy as np

NVR = np.int32(100)
Vref = np.float64(1.) # Convert Vref to 64-bit floating point
neq = 3

# Create particles
XPR = np.zeros((3, NVR), dtype=np.float64)
QPR = np.ones((neq + 1, NVR), dtype=np.float64)


with FortranFile("particles.bin", 'w') as f:
   # Initialize particle
    x_pos_min = -10
    x_pos_max = 10 
    for i in range(NVR):
        x = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        y = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        z = (x_pos_max - x_pos_min) / 100 * i + x_pos_min
        XPR[:, i] = np.array([x,y,z], dtype=np.float64)
        QPR[:3, i] = np.array([1.,2.,1.], dtype=np.float64) * 10
    

    # Open the file for writing in binary mode
    # Write NVR_ext and Vref
    f.write_record(NVR)
    f.write_record(Vref)
    
    # Write XPR and QPR arrays
    for i in range(NVR):
        f.write_record(XPR[:, i])
        f.write_record(QPR[:, i])
    print(f"Successfully wrote {NVR} particles to 'particles.bin'.")