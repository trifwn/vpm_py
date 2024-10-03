import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.io import FortranFile
import h5py
from typing import Any  

def process_particle_output_file_h5(filename: str, folder: str | None= None):
    """Process a single particle file stored in HDF5 format and return the data arrays."""
    if folder:
        filename = folder + filename

    with h5py.File(filename, 'r') as f:
        # Read the dataset groups

        try:
            # Concatenate particle positions into a single array of shape (3, N)
            particle_positions = f['XPR'][:].T
            # Concatenate particle velocities into a single array of shape (3, N)
            particle_velocities = f['UPR'][:].T
            # Concatenate particle vorticity into a single array of shape (3, N)
            particle_charges = f['QPR'][:].T
            # Create a placeholder for particle deformations, as the original code assumes this array
            particle_deformations = f['GPR'][:].T
        except KeyError:
            print('Error reading the dataset groups for file: ', filename)
            return None

    return particle_positions, particle_velocities, particle_charges, particle_deformations

def process_pm_output_file_h5(filename: str, folder: str | None = None):
    """Process a single PM solution file stored in HDF5 format and return the data arrays."""
    if folder:
        filename = folder + filename

    with h5py.File(filename, 'r') as f:
        # Read the dataset groups
        XS = f['X'][:]
        YS = f['Y'][:]
        ZS = f['Z'][:]
        UXS = f['U'][:]
        UYS = f['V'][:]
        UZS = f['W'][:]
        RHS = f['RHS'][:]
        SOL = f['SOL'][:]

        try:
            DEFORM = f['DEFORMX'][:]
            DEFORMY = f['DEFORMY'][:]
            DEFORMZ = f['DEFORMZ'][:]
        except KeyError:
            pass 

        XS = np.moveaxis(XS, [0, 1, 2], [2, 1, 0])
        YS = np.moveaxis(YS, [0, 1, 2], [2, 1, 0])
        ZS = np.moveaxis(ZS, [0, 1, 2], [2, 1, 0])
        UXS = np.moveaxis(UXS, [0, 1, 2], [2, 1, 0])
        UYS = np.moveaxis(UYS, [0, 1, 2], [2, 1, 0])
        UZS = np.moveaxis(UZS, [0, 1, 2], [2, 1, 0])
        RHS = np.moveaxis(RHS, [0, 1, 2, 3], [3, 2, 1, 0])
        SOL = np.moveaxis(SOL, [0, 1, 2, 3], [3, 2, 1, 0])

    # Mesh grid needs to have shape adjusted (as done in the original function)
    neq = RHS.shape[0]
    mesh_positions = np.array([XS, YS, ZS])
    mesh_velocities = np.array([UXS, UYS, UZS])
    mesh_charges = np.array(RHS)
    mesh_solution = np.array(SOL)
    # Create a placeholder for mesh deformations, as the original code assumes this array
    mesh_deformations = np.zeros_like(mesh_positions)

    return neq, mesh_positions, mesh_velocities, mesh_charges, mesh_deformations, mesh_solution

def process_particle_ouput_file_dat(filename: str , folder: str | None = None):
    """Process a single particle file and return the data arrays."""

    data: dict[str, Any] = {
        "XS": [], "YS": [], "ZS": [], "QXS": [], "QYS": [], "QZS": [], "QMAG": [],
                                      "UXS": [], "UYS": [], "UZS": [], "UMAG": []
    }
    if folder:
        filename = folder + filename
    with open(filename) as file:
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
    particle_charges = np.vstack((QXS, QYS, QZS))
    particle_deformations = np.zeros_like(particle_positions)
    return particle_positions, particle_velocities, particle_charges, particle_deformations

def process_pm_output_file_dat(filename: str, folder: str | None = None):
    if folder:
        filename = folder + filename
    with open(filename) as f:
        lines = f.readlines()
        vars = [v.replace('"','') for v in lines[0].split()[2:]]
        sizes = [i for i in lines[1].split()]
        IS, JS, KS = [int(sizes[i][2:]) for i in range(1, 4)]
    df = pd.read_csv(filename, sep=r'\s+', header=None, skiprows=3) 
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
    for k, v in mesh_data.items():
        mesh_data[k] = np.moveaxis(v, [0, 1, 2], [2, 1, 0])

    mesh_positions = np.array([mesh_data["XS"], mesh_data["YS"], mesh_data["ZS"]])
    mesh_velocities = np.array([mesh_data["UXS"], mesh_data["UYS"], mesh_data["UZS"]])
    mesh_charges = np.array([mesh_data["QXS"], mesh_data["QYS"], mesh_data["QZS"]])
    mesh_deformations = np.zeros_like(mesh_positions)

    return mesh_positions, mesh_velocities, mesh_charges, mesh_deformations

def process_pm_output_file(filename: str, folder: str | None = None):
    if filename.endswith('.h5'):
        return process_pm_output_file_h5(filename, folder)
    else:
        return process_pm_output_file_dat(filename, folder)

def process_particle_ouput_file(filename: str, folder: str | None = None):
    if filename.endswith('.h5'):
        return process_particle_output_file_h5(filename, folder)
    else:
        return process_particle_ouput_file_dat(filename, folder)

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