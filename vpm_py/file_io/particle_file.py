import h5py
import numpy as np
from typing import Any

def process_particle_file_h5(filename: str, folder: str | None= None):
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


def process_particle_file_dat(filename: str , folder: str | None = None):
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
