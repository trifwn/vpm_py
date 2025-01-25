import numpy as np
import pandas as pd
import h5py

def process_pm_file_h5(filename: str, folder: str | None = None):
    """Process a single PM solution file stored in HDF5 format and return the data arrays."""
    if folder:
        filename = folder + filename

    with h5py.File(filename, 'r') as f:
        # Read the dataset groups
        XS = f['X'][:]
        YS = f['Y'][:]
        ZS = f['Z'][:]
        RHS = f['RHS'][:]
        SOL = f['SOL'][:]

        XS = np.moveaxis(XS, [0, 1, 2], [2, 1, 0])
        YS = np.moveaxis(YS, [0, 1, 2], [2, 1, 0])
        ZS = np.moveaxis(ZS, [0, 1, 2], [2, 1, 0])
        RHS = np.moveaxis(RHS, [0, 1, 2, 3], [3, 2, 1, 0])
        SOL = np.moveaxis(SOL, [0, 1, 2, 3], [3, 2, 1, 0])

        # Mesh grid needs to have shape adjusted (as done in the original function)
        neq = RHS.shape[0]
        mesh_positions = np.array([XS, YS, ZS])
        mesh_solution = np.array(SOL)
        mesh_charges = np.array(RHS)
        try:
            UPM = f['VEL'][:]
            UPM = np.moveaxis(UPM, [0, 1, 2, 3], [3, 2, 1, 0])
            mesh_velocities = np.array(UPM) 
        except KeyError:
            mesh_velocities = np.zeros_like(mesh_positions)

        try:
            STRETCHING = f['VORTEXSTRETCH'][:]
            STRETCHING = np.moveaxis(STRETCHING, [0, 1, 2, 3], [3, 2, 1, 0])
            mesh_vortex_stretching = np.array(STRETCHING)
        except KeyError:
            mesh_vortex_stretching = np.zeros_like(mesh_positions)

        try:
            P = f['P'][:]
            Q = f['Q'][:]
            P = np.moveaxis(P, [0, 1, 2], [2, 1, 0])
            Q = np.moveaxis(Q, [0, 1, 2], [2, 1, 0])
            mesh_p = np.array(P)
            mesh_q = np.array(Q)
        except KeyError:
            mesh_p = np.zeros_like(mesh_positions)
            mesh_q = np.zeros_like(mesh_positions)
        mesh_u = np.zeros_like(mesh_positions)
    return neq, mesh_positions, mesh_velocities, mesh_charges, mesh_vortex_stretching, mesh_solution, mesh_p, mesh_q, mesh_u

def process_pm_file_dat(filename: str, folder: str | None = None):
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
    mesh_vortex_stretching = np.zeros_like(mesh_positions)

    return mesh_positions, mesh_velocities, mesh_charges, mesh_vortex_stretching
