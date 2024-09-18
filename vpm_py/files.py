import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.io import FortranFile

def process_particle_ouput_file(filename, folder):
    """Process a single particle file and return the data arrays."""
    data = {
        "XS": [], "YS": [], "ZS": [], "QXS": [], "QYS": [], "QZS": [], "QMAG": [],
                                      "UXS": [], "UYS": [], "UZS": [], "UMAG": []
    }
    with open(folder + filename) as file:
        lines = file.readlines()
        for l in lines[2:]:
            l = l.split()
            try:
                data["XS"].append(float(l[1]))
                data["YS"].append(float(l[2]))
                data["ZS"].append(float(l[3]))
                data["UXS"].append(float(l[4]))
                data["UYS"].append(float(l[5]))
                data["UZS"].append(float(l[6]))
                data["QXS"].append(float(l[7]))
                data["QYS"].append(float(l[8]))
                data["QZS"].append(float(l[9]))
                data["UMAG"].append(np.sqrt(data["UXS"][-1]**2 + data["UYS"][-1]**2 + data["UZS"][-1]**2))
                data["QMAG"].append(np.sqrt(data["QXS"][-1]**2 + data["QYS"][-1]**2 + data["QZS"][-1]**2))
            except ValueError:
                print('Error in line: ', l)
                continue
    # Concatenate particle positions into a single array of shape (3, N)
    XS = np.array(data["XS"])
    YS = np.array(data["YS"])
    ZS = np.array(data["ZS"])
    particle_positions = np.vstack((XS, YS, ZS))
    # Concatenate particle velocities into a single array of shape (3, N)
    UXS = np.array(data["UXS"])
    UYS = np.array(data["UYS"])
    UZS = np.array(data["UZS"])
    particle_velocities = np.vstack((UXS, UYS, UZS))
    # Concatenate particle vorticity into a single array of shape (3, N)
    QXS = np.array(data["QXS"])
    QYS = np.array(data["QYS"])
    QZS = np.array(data["QZS"])
    particle_strengths = np.vstack((QXS, QYS, QZS))
    particle_deformations = np.zeros_like(particle_positions)
    return particle_positions, particle_velocities, particle_strengths, particle_deformations

def process_pm_output_file(f, folder):
    with open(folder + f) as file:
        lines = file.readlines()
        vars = [v.replace('"','') for v in lines[0].split()[2:]]
        sizes = [i for i in lines[1].split()]
        IS, JS, KS = [int(sizes[i][2:]) for i in range(1, 4)]
    df = pd.read_csv(folder + f, sep=r'\s+', header=None, skiprows=3) 
    df.columns = vars
    mesh_data = {
        "XS": df['X'].values.reshape(KS, JS, IS, order='C'),
        "YS": df['Y'].values.reshape(KS, JS, IS, order='C'),
        "ZS": df['Z'].values.reshape(KS, JS, IS, order='C'),
        "UXS": df['U'].values.reshape(KS, JS, IS, order='C'),
        "UYS": df['V'].values.reshape(KS, JS, IS, order='C'),
        "UZS": df['W'].values.reshape(KS, JS, IS, order='C'),
        "QXS": df['VORTX'].values.reshape(KS, JS, IS, order='C'),
        "QYS": df['VORTY'].values.reshape(KS, JS, IS, order='C'),
        "QZS": df['VORTZ'].values.reshape(KS, JS, IS, order='C')
    }
    mesh_data["UMAG"] = np.sqrt(mesh_data["UXS"]**2 + mesh_data["UYS"]**2 + mesh_data["UZS"]**2)
    mesh_data["QMAG"] = np.sqrt(mesh_data["QXS"]**2 + mesh_data["QYS"]**2 + mesh_data["QZS"]**2)

    return mesh_data

def process_multiple_files(files, folder, process_func):
    with mp.Pool() as pool:
        results = pool.starmap(process_func, [(f, folder) for f in files])
    return results

def process_multiple_pm_files(files, folder):
    return process_multiple_files(files, folder, process_pm_output_file)

def process_multiple_particle_files(files, folder):
    return process_multiple_files(files, folder, process_particle_ouput_file)

def write_particle_to_bin(
    XPR, 
    QPR, 
    Vref=1.0,
    filename="particles.bin"
):
    NVR = XPR.shape[1]
    neq = QPR.shape[0] - 1

    with FortranFile(filename, 'w') as f:
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