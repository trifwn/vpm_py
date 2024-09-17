import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes

folder = 'results_euler/'

# Get the files in the folder
import os
files = os.listdir(folder)

cap = 400

files_vr = [f for f in files if f.endswith('particles.dat')][:cap]
files_vr.sort()
files_sol = [f for f in files if f.endswith('pm_output.dat')][:cap] 
files_sol.sort()

ALL_XS = []
ALL_YS = []
ALL_ZS = []
ALL_UXS = []
ALL_UYS = []
ALL_UZS = []
ALL_QXS = []
ALL_QYS = []
ALL_QZS = []

def process_file(filename):
    """Process a single file and return the data arrays."""
    data = {
        "XS": [],
        "YS": [],
        "ZS": [],
        "UXS": [],
        "UYS": [],
        "UZS": [],
        "QXS": [],
        "QYS": [],
        "QZS": []
    }

    with open(folder + filename) as file:
        lines = file.readlines()
        vars = lines[0].split()

        # From the file name get the NTIME
        ntime = int(filename[:5])

        # Initialize lists for each variable
        xs, ys, zs, uxs, uys, uzs, qxs, qys, qzs = [], [], [], [], [], [], [], [], []

        for l in lines[2:]:
            l = l.split()
            try:
                id = int(l[0])
                x = float(l[1])
                y = float(l[2])
                z = float(l[3])
                u = float(l[4])
                v = float(l[5])
                w = float(l[6])
                vort_x = float(l[7])
                vort_y = float(l[8])
                vort_z = float(l[9])
            except ValueError:
                print('Error in line: ', l)
                continue

            xs.append(x)
            ys.append(y)
            zs.append(z)
            uxs.append(u)
            uys.append(v)
            uzs.append(w)
            qxs.append(vort_x)
            qys.append(vort_y)
            qzs.append(vort_z)

        # Convert lists to NumPy arrays
        data["XS"] = np.array(xs)
        data["YS"] = np.array(ys)
        data["ZS"] = np.array(zs)
        data["UXS"] = np.array(uxs)
        data["UYS"] = np.array(uys)
        data["UZS"] = np.array(uzs)
        data["QXS"] = np.array(qxs)
        data["QYS"] = np.array(qys)
        data["QZS"] = np.array(qzs)

    return data

# Process files in parallel
with mp.Pool() as pool:
    results = pool.map(process_file, files_vr)

# Collect results
for result in results:
    ALL_XS.append(result["XS"])
    ALL_YS.append(result["YS"])
    ALL_ZS.append(result["ZS"])
    ALL_UXS.append(result["UXS"])
    ALL_UYS.append(result["UYS"])
    ALL_UZS.append(result["UZS"])
    ALL_QXS.append(result["QXS"])
    ALL_QYS.append(result["QYS"])
    ALL_QZS.append(result["QZS"])

ALL_PM_XS = []
ALL_PM_YS = []
ALL_PM_ZS = []
ALL_PM_UXS = []
ALL_PM_UYS = []
ALL_PM_UZS = []
ALL_PM_QXS = []
ALL_PM_QYS = []
ALL_PM_QZS = []
ALL_PM_PSIX = []
ALL_PM_PSIY = []
ALL_PM_PSIZ = []
ALL_PM_DEFORMX = []
ALL_PM_DEFORMY = []
ALL_PM_DEFORMZ = []

def process_file_pm(f):
    with open(folder + f) as file:
        lines = file.readlines()
        vars = lines[0].split()
        vars = vars[2:]
        vars = [v.replace('"','') for v in vars]
        
        # From the file name get the NTIME
        ntime = int(f[:5])
        sizes = [i for i in lines[1].split()]
        IS = int(sizes[1][2:])
        JS = int(sizes[2][2:])
        KS = int(sizes[3][2:])
        
    df = pd.read_csv(folder + f,sep=r'\s+',header=None, skiprows=3) 
    # Make the vars the column names
    df.columns = vars
    return df, IS, JS, KS

# Process files in parallel
with mp.Pool() as pool:
    results = pool.map(process_file_pm, files_sol)

# Collect results
for result in results:
    df, IS, JS, KS = result

    XS = df['X'].values.reshape(KS, JS, IS, order='C')
    YS = df['Y'].values.reshape(KS, JS, IS, order='C')
    ZS = df['Z'].values.reshape(KS, JS, IS, order='C')
    UXS = df['U'].values.reshape(KS, JS, IS, order='C')
    UYS = df['V'].values.reshape(KS, JS, IS, order='C')
    UZS = df['W'].values.reshape(KS, JS, IS, order='C')
    QXS = df['VORTX'].values.reshape(KS, JS, IS, order='C')
    QYS = df['VORTY'].values.reshape(KS, JS, IS, order='C')
    QZS = df['VORTZ'].values.reshape(KS, JS, IS, order='C')
    PSIX = df['PSI1'].values.reshape(KS, JS, IS, order='C')
    PSIY = df['PSI2'].values.reshape(KS, JS, IS, order='C')
    PSIZ = df['PSI3'].values.reshape(KS, JS, IS, order='C')
    DEFORMX = df['DEFORMX'].values.reshape(KS, JS, IS, order='C')
    DEFORMY = df['DEFORMY'].values.reshape(KS, JS, IS, order='C')
    DEFORMZ = df['DEFORMZ'].values.reshape(KS, JS, IS, order='C')

    ALL_PM_XS.append(XS)
    ALL_PM_YS.append(YS)
    ALL_PM_ZS.append(ZS)
    ALL_PM_UXS.append(UXS)
    ALL_PM_UYS.append(UYS)
    ALL_PM_UZS.append(UZS)
    ALL_PM_QXS.append(QXS)
    ALL_PM_QYS.append(QYS)
    ALL_PM_QZS.append(QZS)
    ALL_PM_PSIX.append(PSIX)
    ALL_PM_PSIY.append(PSIY)
    ALL_PM_PSIZ.append(PSIZ)
    ALL_PM_DEFORMX.append(DEFORMX)
    ALL_PM_DEFORMY.append(DEFORMY)
    ALL_PM_DEFORMZ.append(DEFORMZ)
    
print(len(ALL_XS))
print(len(ALL_PM_XS))
x_min = np.min([np.min(x) for x in ALL_XS])
x_max = np.max([np.max(x) for x in ALL_XS])
y_min = np.min([np.min(y) for y in ALL_YS])
y_max = np.max([np.max(y) for y in ALL_YS])
z_min = np.min([np.min(z) for z in ALL_ZS])
z_max = np.max([np.max(z) for z in ALL_ZS])
print(x_min, x_max)
print(y_min, y_max)
print(z_min, z_max)

start_frame = 1
total_frames = len(ALL_XS) - start_frame

# Create the figure with 4 subplots
fig = plt.figure(figsize=(20, 10))
ax_particles: Axes3D = fig.add_subplot(221, projection='3d')
ax_grid_pm: Axes3D = fig.add_subplot(222, projection='3d')
ax_vorticity1: Axes = fig.add_subplot(223)
ax_vorticity2: Axes = fig.add_subplot(224)

# Set labels and views for 3D plots
for ax in [ax_particles, ax_grid_pm]:
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev=30, azim=30)
    ax.set_aspect('equal')
    ax.set_xlim(x_min, x_max)

# Get the shapes of the arrays
shape = ALL_PM_XS[start_frame].shape

# Calculate the middle indices for each dimension
mid_x = shape[2] // 2
mid_y = shape[1] // 2
mid_z = shape[0] // 2

# Create proper slices (remember, the order is Z, Y, X in the array)
slice_xy = (mid_z, slice(None), slice(None))  # X-Y plane at middle Z
slice_xz = (slice(None), mid_y, slice(None))  # X-Z plane at middle Y

# Initial scatter plots for particles and mesh points
particle_velocities = np.sqrt(ALL_UXS[start_frame]**2 + ALL_UYS[start_frame]**2 + ALL_UZS[start_frame]**2)
particles = ax_particles.scatter(ALL_XS[start_frame], ALL_YS[start_frame], ALL_ZS[start_frame], 
                                c=particle_velocities, s=10, marker='o')
fig.colorbar(particles, ax=ax_particles)

mesh_velocities = np.sqrt(ALL_PM_UXS[start_frame]**2 + ALL_PM_UYS[start_frame]**2 + ALL_PM_UZS[start_frame]**2)
mesh_points = ax_grid_pm.scatter(ALL_PM_XS[start_frame], ALL_PM_YS[start_frame], ALL_PM_ZS[start_frame],
                                c=mesh_velocities, s=0.5, marker='x')
fig.colorbar(mesh_points, ax=ax_grid_pm)

# Initial vorticity plots
vorticity = np.sqrt(ALL_PM_QXS[start_frame]**2 + ALL_PM_QYS[start_frame]**2 + ALL_PM_QZS[start_frame]**2)

# X-Y plane contour
x_xy = ALL_PM_XS[start_frame][slice_xy].T
y_xy = ALL_PM_YS[start_frame][slice_xy].T
vorticity_xy = vorticity[slice_xy].T
contour1 = ax_vorticity1.contourf(x_xy, y_xy, vorticity_xy, cmap='viridis')
fig.colorbar(contour1, ax=ax_vorticity1)

# X-Z plane contour
x_xz = ALL_PM_XS[start_frame][slice_xz].T
z_xz = ALL_PM_ZS[start_frame][slice_xz].T
vorticity_xz = vorticity[slice_xz].T
contour2 = ax_vorticity2.contourf(x_xz, z_xz, vorticity_xz, cmap='viridis')
fig.colorbar(contour2, ax=ax_vorticity2)

# Set labels for all plots
ax_particles.set_title('Particles')
ax_grid_pm.set_title('Mesh Points')
ax_vorticity1.set_title(f'Vorticity (X-Y plane) Z={ALL_PM_ZS[start_frame][mid_z, 0, 0]:.2f}')
ax_vorticity2.set_title(f'Vorticity (X-Z plane) Y={ALL_PM_YS[start_frame][0, mid_y, 0]:.2f}')

ax_vorticity1.set_xlabel('X')
ax_vorticity1.set_ylabel('Y')
ax_vorticity2.set_xlabel('X')
ax_vorticity2.set_ylabel('Z')

def update(frame):
    frame = frame + start_frame
    print(f"Frame {frame} / {total_frames}")
    global contour1, contour2, particles, mesh_points

    # Update particle positions
    particle_velocities = np.sqrt(ALL_UXS[frame]**2 + ALL_UYS[frame]**2 + ALL_UZS[frame]**2)
    particles.set_offsets(np.c_[ALL_XS[frame], ALL_YS[frame]])
    particles.set_3d_properties(ALL_ZS[frame], 'z')
    particles.set_array(particle_velocities)
    particles.set_clim(particle_velocities.min(), particle_velocities.max())
    ax_particles.autoscale_view()

    # Calculate new limits for particles
    x_min, x_max = np.min(ALL_XS[frame]), np.max(ALL_XS[frame])
    y_min, y_max = np.min(ALL_YS[frame]), np.max(ALL_YS[frame])
    z_min, z_max = np.min(ALL_ZS[frame]), np.max(ALL_ZS[frame])
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    # Calculate the range and center for each dimension
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    z_center = (z_max + z_min) / 2

    # Expand the range to occupy 75% of the view
    expansion_factor = 1 / 0.75
    x_range *= expansion_factor
    y_range *= expansion_factor
    z_range *= expansion_factor

    # Set new limits
    ax_particles.set_xlim(x_center - x_range/2, x_center + x_range/2)
    ax_particles.set_ylim(y_center - y_range/2, y_center + y_range/2)
    ax_particles.set_zlim(z_center - z_range/2, z_center + z_range/2)

    # Update mesh point positions
    mesh_velocities = np.sqrt(ALL_PM_UXS[frame]**2 + ALL_PM_UYS[frame]**2 + ALL_PM_UZS[frame]**2)
    mesh_points.set_offsets(np.c_[ALL_PM_XS[frame].ravel(), ALL_PM_YS[frame].ravel()])
    mesh_points.set_3d_properties(ALL_PM_ZS[frame].ravel(), 'z')
    mesh_points.set_array(mesh_velocities.ravel())
    mesh_points.set_clim(mesh_velocities.min(), mesh_velocities.max())
    ax_grid_pm.autoscale_view()

    # Calculate new limits for mesh points (similar to particles)
    x_min, x_max = np.min(ALL_PM_XS[frame]), np.max(ALL_PM_XS[frame])
    y_min, y_max = np.min(ALL_PM_YS[frame]), np.max(ALL_PM_YS[frame])
    z_min, z_max = np.min(ALL_PM_ZS[frame]), np.max(ALL_PM_ZS[frame])

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    z_center = (z_max + z_min) / 2

    x_range *= expansion_factor
    y_range *= expansion_factor
    z_range *= expansion_factor

    ax_grid_pm.set_xlim(x_center - x_range/2, x_center + x_range/2)
    ax_grid_pm.set_ylim(y_center - y_range/2, y_center + y_range/2)
    ax_grid_pm.set_zlim(z_center - z_range/2, z_center + z_range/2)

    vorticity = np.sqrt(ALL_PM_QXS[frame]**2 + ALL_PM_QYS[frame]**2 + ALL_PM_QZS[frame]**2)    
    # X-Y plane contour
    mid_z = ALL_PM_ZS[frame].shape[0] // 2
    slice_xy = (mid_z, slice(None), slice(None))
    x_xy = ALL_PM_XS[frame][slice_xy].T
    y_xy = ALL_PM_YS[frame][slice_xy].T
    vorticity_xy = vorticity[slice_xy].T
    for c in contour1.collections:
        c.remove()
    contour1 = ax_vorticity1.contourf(x_xy, y_xy, vorticity_xy, cmap='viridis')
    # Update the boundary of the contour plot
    ax_vorticity1.set_xlim(x_xy.min(), x_xy.max())
    ax_vorticity1.set_ylim(y_xy.min(), y_xy.max())
    # Update title
    ax_vorticity1.set_title(f'Vorticity (X-Y plane) Z={ALL_PM_ZS[frame][mid_z, 0, 0]:.2f}')

    # X-Z plane contour
    mid_y = ALL_PM_YS[frame].shape[1] // 2
    slice_xz = (slice(None), mid_y, slice(None))
    x_xz = ALL_PM_XS[frame][slice_xz].T
    z_xz = ALL_PM_ZS[frame][slice_xz].T
    vorticity_xz = vorticity[slice_xz].T
    
    for c in contour2.collections:
        c.remove()
    contour2 = ax_vorticity2.contourf(x_xz, z_xz, vorticity_xz, cmap='viridis')

    # Update the boundary of the contour plot
    ax_vorticity2.set_xlim(x_xz.min(), x_xz.max())
    ax_vorticity2.set_ylim(z_xz.min(), z_xz.max())
    # Update title
    ax_vorticity1.set_title(f'Vorticity (X-Y plane) Z={ALL_PM_ZS[frame][mid_z, 0, 0]:.2f}')

    # Update figure title
    fig.suptitle(f"Frame {frame} / {total_frames + start_frame - 1}, Realtime = {frame * 0.1:.1f}s")

    return particles, mesh_points, contour1, contour2

fig.tight_layout()
ani = FuncAnimation(fig, update, frames=total_frames-1, blit=False, repeat=False, interval=100)
ani.event_source.stop()
# Save the animation
ani.save(f"{folder}/anim.gif")