import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PathCollection
from matplotlib.contour import QuadContourSet
from matplotlib.axes import Axes
import numpy as np
from vpm_py.vpm_io import print_green, print_IMPORTANT
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
from matplotlib.animation import FuncAnimation
import os
from tqdm import tqdm
from vpm_py.files import process_particle_ouput_file, process_pm_output_file

class Visualizer:
    """
    Class to plot the live position of particles in a 3D domain.
    """

    def __init__(
        self,
        plot_particles: bool = True,
        plot_particle_mesh: bool = False,
        plot_particle_mesh_slices: bool = False, 
        quantity_to_plot: str = "Particle Strengths"
    ) -> None:
        """
        Initialize the particle plot.
        """
        print_IMPORTANT("Initializing the plot")
        self.plot_options = []
        if plot_particles:
            self.plot_options.append("particles")
        if plot_particle_mesh:
            self.plot_options.append("mesh")
        if plot_particle_mesh_slices:
            self.plot_options.append("slice")
        
        self.fig = plt.figure(figsize=(15, 10)) 
        self.axes: dict[str, Axes | Axes3D] = {}
        self.plots: dict[str, PathCollection | QuadContourSet] = {}
        self.quantity_to_plot = quantity_to_plot
        self.quantity_dict = {
            "Particle Velx": "UXS",
            "Particle Vely": "UYS",
            "Particle Velz": "UZS",
            "Particle Velocities": "UMAG",
            "Particle Strength X": "QXS",
            "Particle Strength Y": "QYS",
            "Particle Strength Z": "QZS",
            "Particle Strengths": "QMAG",
        }
        self.quantity_key = self.quantity_dict[quantity_to_plot]
        self.strength_threshold = None
        
        self._setup_subplots()
        
        self.fig.suptitle("Frame 0 / 0, Realtime = 0.0s")
        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.9, hspace=0.3, wspace=0.3, bottom=0.1)
        plt.ion()
        self.fig.show()

    def _setup_subplots(self):
        plot_count = len(self.plot_options)
        if 'slice' in self.plot_options:
            plot_count += 1
        self.rows = 1 if plot_count == 1 else 2
        self.columns = plot_count // 2 + plot_count % 2

        idx = 1
        if 'particles' in self.plot_options:
            self._setup_particle_subplot(idx)
            idx += 1

        if 'mesh' in self.plot_options:
            self._setup_mesh_subplot(idx)
            idx += 1

        if 'slice' in self.plot_options:
            self._setup_slice_subplots(idx)

    def _setup_particle_subplot(self, idx):
        ax: Axes3D = self.fig.add_subplot(self.columns, self.rows, idx, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'{self.quantity_to_plot} on particles')
        self.axes['particles'] = ax
        self.plots['particles'] = ax.scatter([], [], [], c=[], s=10, cmap='viridis')
        self.plots['particles_cbar'] = plt.colorbar(self.plots['particles'], ax=ax)

    def _setup_mesh_subplot(self, idx):
        ax: Axes3D = self.fig.add_subplot(self.columns, self.rows, idx, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'{self.quantity_to_plot} on Particle Mesh')
        self.axes['mesh'] = ax
        self.plots['mesh'] = ax.scatter([], [], [], c=[], s=0.2, cmap='viridis')
        self.plots['mesh_cbar'] = self.fig.colorbar(self.plots['mesh'], ax=ax)

    def _setup_slice_subplots(self, idx):
        dummy_data = np.zeros((2,2))
        
        # XY slice
        ax_xy: Axes = self.fig.add_subplot(self.columns, self.rows, idx)
        ax_xy.set_xlabel('X')
        ax_xy.set_ylabel('Y')
        self.axes['slice_xy'] = ax_xy
        self.plots['slice_xy'] = ax_xy.contourf(dummy_data, dummy_data, dummy_data, cmap='viridis')
        divider = make_axes_locatable(ax_xy)
        self.axes['slice_xy_cbar'] = divider.append_axes('right', size='5%', pad=0.05)
        self.plots['slice_xy_cbar'] = self.fig.colorbar(self.plots['slice_xy'], cax=self.axes['slice_xy_cbar'], orientation='vertical')
        idx += 1

        # XZ slice
        ax_xz: Axes = self.fig.add_subplot(self.columns, self.rows, idx)
        self.axes['slice_xz'] = ax_xz
        self.plots['slice_xz'] = ax_xz.contourf(dummy_data, dummy_data, dummy_data, cmap='viridis')
        ax_xz.set_xlabel('X')
        ax_xz.set_ylabel('Z')
        divider = make_axes_locatable(ax_xz)
        self.axes['slice_xz_cbar'] = divider.append_axes('right', size='5%', pad=0.05)
        self.plots['slice_xz_cbar'] = self.fig.colorbar(self.plots['slice_xz'], cax=self.axes['slice_xz_cbar'], orientation='vertical')

    def results_animation_loop(
        self, 
        frame: int, 
        total_frames: int,
        folder: str,
        files_vr: list[str],
        files_sol: list[str],
        pbar: tqdm | None = None, 
    ):
        if 'particles' in self.plot_options:
            if pbar:
                pbar.set_description(f"Processing Particle Output File for frame: {frame}/{total_frames}")
            (
                particle_positions, particle_velocities, particle_strengths, particle_deformations
            ) = process_particle_ouput_file(files_vr[frame], folder)
            if pbar:
                pbar.set_description(f"Updating Particle Plot for frame: {frame}/{total_frames}")
            self._update_particle_plot(
                frame,
                particle_positions,
                particle_strengths,
                particle_velocities,
                particle_deformations,
            )
        
        if 'mesh' in self.plot_options or 'slice' in self.plot_options:
            if pbar:
                pbar.set_description(f"Processing PM Output File for frame: {frame}/{total_frames}")
            pm_data = process_pm_output_file(files_sol[frame], folder)
            if 'mesh' in self.plot_options:
                if pbar:
                    pbar.set_description(f"Updating Mesh Plot for frame: {frame}/{total_frames}")
                dx, dy, dz = self._calculate_grid_spacing(pm_data)
                dV = dx * dy * dz
                self._update_mesh_plot(pm_data) 
            if 'slice' in self.plot_options:
                if pbar:
                    pbar.set_description(f"Updating Slice Plots for frame: {frame}/{total_frames}")
                self._update_slice_plots(pm_data)
        self._update_figure_title(frame)
        self._render_plot()
        if pbar:
            pbar.update(1)

    def animate_result_folder(
        self, 
        folder: str, 
        total_frames: int,
        start_frame: int = 0
    ):
        files = os.listdir(folder)
        files_vr = sorted([f for f in files if f.endswith('particles.dat')])
        files_sol = sorted([f for f in files if f.endswith('pm_output.dat')])

        pbar = tqdm(total=total_frames, desc="Creating animation")
        ani = FuncAnimation(
            self.fig, 
            self.results_animation_loop, 
            frames=total_frames, blit=False, repeat=False, interval=100,                
            fargs=(total_frames, folder, files_vr, files_sol, pbar)
        )
        ani.event_source.stop()
        ani.save(f"{folder}/animation_{self.quantity_to_plot}.mp4", writer='ffmpeg', fps=10, dpi=300)

    def update_all_plots(
        self,
        iteration: int,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
        pm_positions: np.ndarray,
        pm_strengths: np.ndarray,
        pm_velocities: np.ndarray,
        pm_deformations: np.ndarray,
    ):
        self._update_particle_plot(
            iteration,
            particle_positions,
            particle_strengths,
            particle_velocities,
            particle_deformations,
        )

        mesh_data = self._prepare_mesh_data(pm_positions, pm_velocities, pm_strengths)

        if 'mesh' in self.plot_options or 'slice' in self.plot_options:

            if 'mesh' in self.plot_options:
                self._update_mesh_plot(mesh_data)

            if 'slice' in self.plot_options:
                self._update_slice_plots(mesh_data, mesh_quantity)
        self._update_figure_title(iteration)

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def update_particle_plot(
        self,
        iteration: int,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        self._update_particle_plot(
            iteration,
            particle_positions,
            particle_strengths,
            particle_velocities,
            particle_deformations,
        )
        self._update_figure_title(iteration)
        self._render_plot()

    def _update_particle_plot(
        self,
        iteration: int,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        axes = self.axes    
        plots = self.plots
        quantity_key = self.quantity_key
        plot_options = self.plot_options
        strength_threshold = self.strength_threshold

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
            "UMAG": np.sqrt(
                particle_velocities[0,:]**2 + particle_velocities[1,:]**2 + particle_velocities[2,:]**2
            ),
            "QMAG": np.sqrt(
                particle_strengths[0,:]**2 + particle_strengths[1,:]**2 + particle_strengths[2,:]**2
            ),
        }
        # Plot particles if selected
        if 'particles' in plot_options:
            particle_quantity = particle_data[quantity_key]

            if strength_threshold:
                mask = particle_quantity >= strength_threshold
                particle_subset_q = particle_quantity[mask]
                particle_subset_x = particle_data["XS"][mask]
                particle_subset_y = particle_data["YS"][mask]
                particle_subset_z = particle_data["ZS"][mask]
            else:
                particle_subset_q = particle_quantity
                particle_subset_x = particle_data["XS"]
                particle_subset_y = particle_data["YS"]
                particle_subset_z = particle_data["ZS"]
        
            plots['particles']._offsets3d = (particle_subset_x, particle_subset_y, particle_subset_z)
            plots['particles'].set_array(particle_subset_q)
            plots['particles'].set_clim(particle_quantity.min(), particle_quantity.max())

            # Adjust plot limits
            axes['particles'].set_xlim(1.5 * particle_data["XS"].min(), 1.5 * particle_data["XS"].max())
            axes['particles'].set_ylim(1.5 * particle_data["YS"].min(), 1.5 * particle_data["YS"].max())
            axes['particles'].set_box_aspect([1.5, 1.5, 1.5])
            axes['particles'].set_zlim(1.5 * particle_data["ZS"].min(), 1.5 * particle_data["ZS"].max())

    def _prepare_mesh_data(self, pm_positions, pm_velocities, pm_strengths):
        mesh_data = {
            "XS": pm_positions[0, :, :, :],
            "YS": pm_positions[1, :, :, :],
            "ZS": pm_positions[2, :, :, :],
            "UXS": pm_velocities[0, :, :, :],
            "UYS": pm_velocities[1, :, :, :],
            "UZS": pm_velocities[2, :, :, :],
            "QXS": pm_strengths[0, :, :, :],
            "QYS": pm_strengths[1, :, :, :],
            "QZS": pm_strengths[2, :, :, :],
        }
        mesh_data["UMAG"] = np.sqrt(
            mesh_data["UXS"]**2 + mesh_data["UYS"]**2 + mesh_data["UZS"]**2
        )
        mesh_data["QMAG"] = np.sqrt(
            mesh_data["QXS"]**2 + mesh_data["QYS"]**2 + mesh_data["QZS"]**2
        )
        return mesh_data

    def _calculate_grid_spacing(self, mesh_data):
        dx = (mesh_data["XS"][0, 0, 1] - mesh_data["XS"][0, 0, 0])
        dy = (mesh_data["YS"][0, 1, 0] - mesh_data["YS"][0, 0, 0])
        dz = (mesh_data["ZS"][1, 0, 0] - mesh_data["ZS"][0, 0, 0])
        return dx, dy, dz

    def _update_mesh_plot(self, mesh_data):
        dx, dy, dz = self._calculate_grid_spacing(mesh_data)
        dV = dx * dy * dz
        mesh_quantity = mesh_data[self.quantity_key]
        if self.strength_threshold:
            if self.quantity_key in ["QXS", "QYS", "QZS", "QMAG"]:
                mask = mesh_quantity >= self.strength_threshold * dV
            else:
                mask = mesh_quantity >= self.strength_threshold
            mask = self._get_strength_mask(mesh_quantity, dV)
            mesh_subset_q = mesh_quantity[mask]
            mesh_subset_x = mesh_data["XS"][mask]
            mesh_subset_y = mesh_data["YS"][mask]
            mesh_subset_z = mesh_data["ZS"][mask]
        else:
            mesh_subset_q = mesh_quantity.ravel()
            mesh_subset_x = mesh_data["XS"].ravel()
            mesh_subset_y = mesh_data["YS"].ravel()
            mesh_subset_z = mesh_data["ZS"].ravel()
        
        self.plots['mesh']._offsets3d = (mesh_subset_x, mesh_subset_y, mesh_subset_z)
        self.plots['mesh'].set_array(mesh_subset_q)
        self.plots['mesh'].set_clim(mesh_quantity.min(), mesh_quantity.max())

        self._adjust_mesh_plot_limits(mesh_data)

    def _adjust_mesh_plot_limits(self, mesh_data):
        self.axes['mesh'].set_xlim(1.5 * mesh_data["XS"].min(), 1.5 * mesh_data["XS"].max())
        self.axes['mesh'].set_ylim(1.5 * mesh_data["YS"].min(), 1.5 * mesh_data["YS"].max())
        self.axes['mesh'].set_box_aspect([1.5, 1.5, 1.5])
        self.axes['mesh'].set_zlim(1.5 * mesh_data["ZS"].min(), 1.5 * mesh_data["ZS"].max())

    def _update_slice_plots(self, mesh_data):
        mesh_quantity = mesh_data[self.quantity_key]
        int_y = scipy.integrate.trapezoid(mesh_quantity, axis=1)
        int_xy = scipy.integrate.trapezoid(int_y, axis=0)
        max_z_index = np.argmax(int_xy)
        mid_y = mesh_data["YS"].shape[1] // 2

        self._update_xy_slice(mesh_data, mesh_quantity, max_z_index)
        self._update_xz_slice(mesh_data, mesh_quantity, mid_y)
        self._update_slice_planes(mesh_data, max_z_index, mid_y)

    def _update_xy_slice(self, mesh_data, mesh_quantity, max_z_index):
        slice_xy = (max_z_index, slice(None), slice(None))
        x_xy, y_xy = mesh_data["XS"][slice_xy].T, mesh_data["YS"][slice_xy].T
        quantity_xy = mesh_quantity[slice_xy].T

        self.plots['slice_xy'] = self.axes['slice_xy'].contourf(x_xy, y_xy, quantity_xy, cmap='viridis')
        self.plots['slice_xy'].set_clim(mesh_quantity.min(), mesh_quantity.max())
        self.axes['slice_xy'].set_xlim(y_xy.min(), y_xy.max())
        self.axes['slice_xy'].set_ylim(x_xy.min(), x_xy.max())
        self.axes['slice_xy'].set_title(f'{self.quantity_key} (X-Y plane) Z={mesh_data["ZS"][max_z_index, 0, 0]:.2f}')
        self.axes['slice_xy_cbar'].clear()
        self.fig.colorbar(self.plots['slice_xy'], cax=self.axes['slice_xy_cbar'])

    def _update_xz_slice(self, mesh_data, mesh_quantity, mid_y):
        slice_xz = (slice(None), mid_y, slice(None))
        x_xz, z_xz = mesh_data["XS"][slice_xz].T, mesh_data["ZS"][slice_xz].T
        quantity_xz = mesh_quantity[slice_xz].T

        self.plots['slice_xz'] = self.axes['slice_xz'].contourf(x_xz, z_xz, quantity_xz, cmap='viridis')
        self.plots['slice_xz'].set_clim(mesh_quantity.min(), mesh_quantity.max())
        self.axes['slice_xz'].set_xlim(x_xz.min(), x_xz.max())
        self.axes['slice_xz'].set_ylim(z_xz.min(), z_xz.max())
        self.axes['slice_xz'].set_title(f'{self.quantity_key} (X-Z plane) Y={mesh_data["YS"][0, mid_y, 0]:.2f}')
        self.axes['slice_xz_cbar'].clear()
        self.fig.colorbar(self.plots['slice_xz'], cax=self.axes['slice_xz_cbar'])

    def _update_slice_planes(self, mesh_data, max_z_index, mid_y):
        slice_plane_ax = self._get_slice_plane_axis()
        if slice_plane_ax:
            self._remove_existing_slice_planes()
            self._add_xy_slice_plane(slice_plane_ax, mesh_data, max_z_index)
            self._add_xz_slice_plane(slice_plane_ax, mesh_data, mid_y)

    def _get_slice_plane_axis(self):
        if 'particles' in self.plot_options:
            return self.axes['particles']
        elif 'mesh' in self.plot_options:
            return self.axes['mesh']
        return None

    def _remove_existing_slice_planes(self):
        for plane in ['particles_xy_plane', 'particles_xz_plane']:
            if plane in self.plots:
                self.plots[plane].remove()
                del self.plots[plane]

    def _add_xy_slice_plane(self, ax, mesh_data, max_z_index):
        min_x, max_x = ax.get_xlim()
        min_y, max_y = ax.get_ylim()
        x_plane, y_plane = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100))
        z_plane = np.ones_like(x_plane) * mesh_data["ZS"][max_z_index, 0, 0]
        self.plots['particles_xy_plane'] = ax.plot_surface(x_plane, y_plane, z_plane, alpha=0.1, rstride=1, cstride=1, color='red')

    def _add_xz_slice_plane(self, ax, mesh_data, mid_y):
        min_x, max_x = ax.get_xlim()
        min_z, max_z = ax.get_zlim()
        x_plane, z_plane = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_z, max_z, 100))
        y_plane = np.ones_like(x_plane) * mesh_data["YS"][0, mid_y, 0]
        self.plots['particles_xz_plane'] = ax.plot_surface(x_plane, y_plane, z_plane, alpha=0.1, rstride=1, cstride=1, color='blue')

    def _update_figure_title(self, iteration, total_iterations = None, t = None):
        string = f"Frame {iteration}"
        if total_iterations:
            string += f" / {total_iterations}"
        if t:
            string += f", Realtime = {t:.2f}s"
        self.fig.suptitle(string)

    def _render_plot(self):
        renderer = self.fig.canvas.get_renderer()
        self.fig.draw(renderer=renderer)
        self.fig.canvas.draw()
        plt.pause(0.05)