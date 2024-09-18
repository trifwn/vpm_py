import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
import scipy.integrate
from tqdm import tqdm
import scipy

from vpm_py.files import process_particle_ouput_file as process_file_par
from vpm_py.files import process_pm_output_file as process_file_pm

# User input for plot selection
print("Select the plots you want to animate:")
print("1. Particles (scatter)")
print("2. Mesh (scatter)")
print("3. Mesh slice (contours)")
print("Enter the numbers of your choices separated by commas (e.g., 1,2,3):")
user_choice = input().split(',')
plot_choices = [int(choice.strip()) for choice in user_choice]

# User input for configuration
# print("Enter the folder name where the data is stored:")
# folder = input()

# Configuration
folder = 'results/'
# Get the files in the folder
files = os.listdir(folder)
files_vr = [f for f in files if f.endswith('particles.dat')]
files_vr.sort()
files_sol = [f for f in files if f.endswith('pm_output.dat')]
files_sol.sort()
# Print the number of files found for each type
print(f"Found {len(files_vr)} particle files and {len(files_sol)} mesh files.")
print("Enter the starting frame:")
start_frame = int(input())
print("Enter the total number of frames:")
total_frames = int(input())
# Select the Quantity to plot:
print("Select the Quantity to plot (Default if Vorticity Magnitude):")
print("1. - Velocity x")
print("2. - Velocity y")
print("3. - Velocity z")
print("4. - Velocity magnitude")
print("5. - Vorticity x")
print("6. - Vorticity y")
print("7. - Vorticity z")
print("8. - Vorticity magnitude")
print("Enter the number of your choice:")
quantity_to_plot = input()
quantity_to_plot = int(quantity_to_plot) if quantity_to_plot else None

name = "Vorticity magnitude"
quantity = "QMAG"
if quantity_to_plot == 1:
    name = "Velocity x"
    quantity = "UXS"
elif quantity_to_plot == 2:
    name = "Velocity y"
    quantity = "UYS"
elif quantity_to_plot == 3:
    name = "Velocity z"
    quantity = "UZS"
elif quantity_to_plot == 4:
    name = "Velocity magnitude"
    quantity = "UMAG"
elif quantity_to_plot == 5:
    name = "Vorticity x"
    quantity = "QXS"
elif quantity_to_plot == 6:
    name = "Vorticity y"
    quantity = "QYS"
elif quantity_to_plot == 7:
    name = "Vorticity z"
    quantity = "QZS"
elif quantity_to_plot == 8:
    name = "Vorticity magnitude"
    quantity = "QMAG"

# User input for threshold
print("Enter the particle/mesh strength threshold (leave empty for no filtering):")
strength_threshold = input()
strength_threshold = float(strength_threshold) if strength_threshold else None

# Create the figure based on user selection
num_plots = len(plot_choices)
if 3 in plot_choices:
    num_plots += 1
fig = plt.figure( figsize=(15, 10)) 
axes: list[Axes | Axes3D] = []
rows = 1 if num_plots == 1 else 2
columns = num_plots // 2 + num_plots % 2

plots = []
plot_index = 0
if 1 in plot_choices:
    axes.append(fig.add_subplot(columns, rows, plot_index+1, projection='3d'))
    axes[plot_index].set_title(f"{name} on particles")
    plots.append(axes[plot_index].scatter([], [], [], c=[], s=10, marker='o'))
    # Add colorbar
    cbar_par = fig.colorbar(plots[plot_index], ax=axes[plot_index])
    plot_index += 1
if 2 in plot_choices:
    axes.append(fig.add_subplot(columns, rows, plot_index+1, projection='3d'))
    plots.append(axes[plot_index].scatter([], [], [], c=[], s=5, marker='x'))
    # Add title
    axes[plot_index].set_title(f"{name} on mesh nodes")
    # Add colorbar
    cbar_mesh = fig.colorbar(plots[plot_index], ax=axes[plot_index])
    plot_index += 1
if 3 in plot_choices:
    dummy_data = np.zeros((2, 2))
    axes.append(fig.add_subplot(columns, rows, plot_index+1))
    plots.append(axes[plot_index].contourf(dummy_data, dummy_data, dummy_data))
    # Add colorbar
    divider = make_axes_locatable(axes[plot_index])
    cax0 = divider.append_axes('right', size='5%', pad=0.05)
    cbar_slice_xy = fig.colorbar(plots[plot_index], cax=cax0, orientation='vertical')

    plot_index += 1
    axes.append(fig.add_subplot(columns, rows, plot_index+1))
    plots.append(axes[plot_index].contourf(dummy_data, dummy_data, dummy_data))
    # Add colorbar
    divider = make_axes_locatable(axes[plot_index])
    cax1 = divider.append_axes('right', size='5%', pad=0.05)
    cbar_slice_xz = fig.colorbar(plots[plot_index], cax=cax1, orientation='vertical')

    plot_index += 1


# Set labels and views for 3D plots
for ax in axes:
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    if isinstance(ax, Axes3D):
        ax.set_zlabel('Z')
        ax.view_init(elev=30, azim=30)

# Create a progress bar
pbar = tqdm(total=total_frames, desc="Creating animation")
fig.suptitle("Frame 0 / 0, Realtime = 0.0s")
# Make the plots as tight as possible
fig.tight_layout()
fig.subplots_adjust(top=0.9, hspace=0.3, wspace=0.3, bottom=0.1)
plt.ion()
fig.show()

def update(frame):
    frame = frame + start_frame
    pbar.update(1)
    pbar.set_description(f"Processing frame {frame}/{total_frames}")

    plot_index = 0
    if 1 in plot_choices:
        (
            particle_positions, particle_velocities, particle_strengths, particle_deformations
        )= process_file_par(files_vr[frame], folder)
        particle_data = {
            "XS": particle_positions[0,:],
            "YS": particle_positions[1,:],
            "ZS": particle_positions[2,:],
            "UXS": particle_velocities[0,:],
            "UYS": particle_velocities[1,:],
            "UZS": particle_velocities[2,:],
            "QXS": particle_strengths[0,:],
            "QYS": particle_strengths[1,:],
            "QZS": particle_strengths[2,:],
            "UMAG": np.sqrt(particle_velocities[0,:]**2 + particle_velocities[1,:]**2 + particle_velocities[2,:]**2),
            "QMAG": np.sqrt(particle_strengths[0,:]**2 + particle_strengths[1,:]**2 + particle_strengths[2,:]**2)
        }
        particle_quantity = particle_data[quantity]

        # Apply particle strength thresholding
        if strength_threshold is not None:
            mask = particle_quantity >= strength_threshold
            particle_subset_q = particle_quantity[mask]
            particle_subset_x = particle_data["XS"]
            particle_subset_y = particle_data["YS"]
            particle_subset_z = particle_data["ZS"]
        else:
            particle_subset_q = particle_quantity
            particle_subset_x = particle_data["XS"]
            particle_subset_y = particle_data["YS"]
            particle_subset_z = particle_data["ZS"]
    
        plots[plot_index]._offsets3d = (particle_subset_x, particle_subset_y, particle_subset_z)
        plots[plot_index].set_array(particle_subset_q)

        # Adjust the limits of the plot
        plots[plot_index].set_clim(particle_quantity.min(), particle_quantity.max())
        axes[plot_index].set_xlim(1.5* particle_data["XS"].min(), 1.5 * particle_data["XS"].max())
        axes[plot_index].set_ylim(1.5* particle_data["YS"].min(), 1.5 * particle_data["YS"].max())
        axes[plot_index].set_zlim(1.5* particle_data["ZS"].min(), 1.5 * particle_data["ZS"].max())
        # Set equal aspect ratio
        axes[plot_index].set_box_aspect([1.5, 1.5, 1.5])
        plot_index += 1

    if 2 in plot_choices or 3 in plot_choices:
        mesh_data = process_file_pm(files_sol[frame], folder)
        dx = (mesh_data["XS"][0, 0, 1] - mesh_data["XS"][0, 0, 0])
        dy = (mesh_data["YS"][0, 1, 0] - mesh_data["YS"][0, 0, 0])
        dz = (mesh_data["ZS"][1, 0, 0] - mesh_data["ZS"][0, 0, 0])
        # Find the volume of a cell
        dV = dx * dy * dz
        mesh_quantity: np.ndarray = mesh_data[quantity]                        

        if 2 in plot_choices:
            # Apply mesh strength thresholding
            if strength_threshold is not None:
                if quantity in ["QXS", "QYS", "QZS", "QMAG"]:
                    strength_threshold_pm = strength_threshold * dV
                else:
                    strength_threshold_pm = strength_threshold
                mask = mesh_quantity >= strength_threshold_pm
                mesh_subset_q = mesh_quantity[mask]
                mesh_subset_xs = mesh_data["XS"][mask]
                mesh_subset_ys = mesh_data["YS"][mask]
                mesh_subset_zs = mesh_data["ZS"][mask]
            else:
                mesh_subset_q = mesh_quantity.ravel()
                mesh_subset_xs = mesh_data["XS"].ravel()
                mesh_subset_ys = mesh_data["YS"].ravel()
                mesh_subset_zs = mesh_data["ZS"].ravel()
            
            plots[plot_index]._offsets3d = (mesh_subset_xs, mesh_subset_ys, mesh_subset_zs)
            plots[plot_index].set_array(mesh_subset_q)
            plots[plot_index].set_clim(mesh_subset_q.min(), mesh_subset_q.max())
            # Make the markers smaller according to the mesh strength
            plots[plot_index].set_sizes(10 * mesh_subset_q / mesh_subset_q.max())

            # Adjust the limits of the plot
            axes[plot_index].set_xlim(1.5 * mesh_data["XS"].min(), 1.5 * mesh_data["XS"].max())
            axes[plot_index].set_ylim(1.5 * mesh_data["YS"].min(), 1.5 * mesh_data["YS"].max())
            axes[plot_index].set_box_aspect([1.5, 1.5, 1.5])
            axes[plot_index].set_zlim(1.5 * mesh_data["ZS"].min(), 1.5 * mesh_data["ZS"].max())

            # Set equal aspect ratio
            plot_index += 1

        if 3 in plot_choices: 
            axes[plot_index].clear()
            axes[plot_index].clear()

            # Find Z-plane with maximum vorticity (integrate over X and Y)
            int_y = scipy.integrate.trapezoid(mesh_quantity, dx=dx, axis=1)
            int_xy = scipy.integrate.trapezoid(int_y, dx=dy, axis=0)
            max_z_index = np.argmax(int_xy)

            slice_xy = (max_z_index, slice(None), slice(None))
            x_xy, y_xy = mesh_data["XS"][slice_xy].T, mesh_data["YS"][slice_xy].T
            quantity_xy = mesh_quantity[slice_xy].T
            
            plots[plot_index] = axes[plot_index].contourf(x_xy, y_xy, quantity_xy, cmap='viridis')
            axes[plot_index].set_xlim(x_xy.min(), x_xy.max())
            axes[plot_index].set_ylim(y_xy.min(), y_xy.max())
            axes[plot_index].set_title(f'{name} (X-Y plane) Z={mesh_data["ZS"][max_z_index, 0, 0]:.2f}')
            axes[plot_index].set_xlabel('X')
            axes[plot_index].set_ylabel('Y')
            # Add colorbar
            cax0.clear()
            fig.colorbar(plots[plot_index], cax=cax0, orientation='vertical')
            
            plot_index += 1

            # Find mid Y-plane 
            mid_y = JS // 2
            slice_xz = (slice(None), mid_y, slice(None))
            x_xz, z_xz = mesh_data["XS"][slice_xz].T, mesh_data["ZS"][slice_xz].T
            quantity_xz = mesh_quantity[slice_xz].T
            
            plots[plot_index] = axes[plot_index].contourf(x_xz, z_xz, quantity_xz, cmap='viridis')
            axes[plot_index].set_xlim(x_xz.min(), x_xz.max())
            axes[plot_index].set_ylim(z_xz.min(), z_xz.max())
            axes[plot_index].set_title(f'{name} (X-Z plane) Y={mesh_data["YS"][0, mid_y, 0]:.2f}')
            axes[plot_index].set_xlabel('X')
            axes[plot_index].set_ylabel('Z')
            # Add colorbar
            cax1.clear()
            fig.colorbar(plots[plot_index], cax=cax1, orientation='vertical')

            plot_index += 1

    fig.suptitle(f"Frame {frame} / {total_frames + start_frame - 1}, Realtime = {frame * 0.1:.1f}s")
    renderer = fig.canvas.get_renderer()
    fig.draw(renderer=renderer)
    fig.canvas.draw()
    plt.pause(0.05)
    return plots

ani = FuncAnimation(fig, update, frames=total_frames-1, blit=False, repeat=False, interval=100)
ani.event_source.stop()
ani.save(f"{folder}/animation_{quantity}.mp4", writer='ffmpeg', fps=10, dpi=300)