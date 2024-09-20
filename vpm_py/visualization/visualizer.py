import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
import numpy as np
import os
from tqdm import tqdm

from vpm_py.vpm_io import print_IMPORTANT
from vpm_py.files import process_particle_ouput_file, process_pm_output_file

from . import ResultPlot
from . import QuantityOfInterest
from . import Plotter
from . import DataAdapter

class Visualizer:
    def __init__(
        self,
        figure_size: tuple[int, int] = (10, 10),
        plot_options: list[ResultPlot] = [],
    ) -> None:
        """
        Initialize the particle plot.
        """
        print_IMPORTANT("Initializing the plot")
        self.plot_options: list[ResultPlot] = plot_options
        self.num_plots = len(self.plot_options) 
        self.fig = plt.figure(figsize=figure_size) 
        self.plotters: dict[str, Plotter] = {}
        self.data_adapters: dict[str, DataAdapter] = {}

        self.mesh_quantities: dict[str, QuantityOfInterest] = {}
        self.particle_quantities: dict[str, QuantityOfInterest] = {}

        self.has_mesh = False
        self.has_particles = False
        self._setup_subplots()

    def _setup_subplots(self):
        self.rows = 1 if self.num_plots <= 2 else 2
        self.columns = self.num_plots // self.rows + self.num_plots % 2

        gs = gridspec.GridSpec( self.rows, self.columns)

        for i, plot_option in enumerate(self.plot_options):
            row_num = i // self.columns
            col_num = i % self.columns

            if plot_option.get_plot_dimensions() == 2:
                ax = self.fig.add_subplot(gs[row_num, col_num])
            else:
                ax = self.fig.add_subplot(gs[row_num, col_num], projection='3d')
            
            name = f'{plot_option.type}_{plot_option.quantity}_{plot_option.component}_{i}'
            self.plotters[name] = plot_option.get_plotter(self.fig, ax)
            self.data_adapters[name] = plot_option.get_data_adapter()
            if plot_option.type == 'particle':
                self.has_particles = True
                self.particle_quantities[name] = plot_option.get_quantity()
            elif plot_option.type == 'mesh':
                self.has_mesh = True
                self.mesh_quantities[name] = plot_option.get_quantity()
    
        self.gridspec = gs
        self.fig.suptitle("Frame 0, Realtime = 0.0s")
        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.9, hspace=0.3, wspace=0.3, bottom=0.1)
        plt.ion()
        self.fig.show()

    def load_results_from_disk(
        self, 
        frame: int, 
        total_frames: int,
        file_vr: str,
        file_pm: str,
        time: float | None = None
    ):
        if self.has_particles:
            (
                particle_positions, particle_velocities, particle_strengths, particle_deformations
            ) = process_particle_ouput_file(file_vr)
            self._update_particle_plots(
                particle_positions,
                particle_strengths,
                particle_velocities,
                particle_deformations,
            )
        
        if self.has_mesh:
            (
                mesh_positions, mesh_velocities, mesh_strengths, mesh_deformations
            ) = process_pm_output_file(file_pm)
            self._update_mesh_plots(
                mesh_positions,
                mesh_strengths,
                mesh_velocities,
                mesh_deformations,
            )
        if time:
            self._update_figure_title(frame, total_frames, time)
        else:
            self._update_figure_title(frame, total_frames)
        self._render_plot()

    def animate_result_folder(
        self, 
        folder: str, 
        total_frames: int,
        start_frame: int = 0,
        filename: str | None = None,
        dt: float | None = None,
        format: str = 'mp4'
    ):
        files = os.listdir(folder)
        files_vr = sorted([f for f in files if f.endswith('particles.dat')])
        files_sol = sorted([f for f in files if f.endswith('pm_output.dat')])

        pbar = tqdm(total=total_frames, desc="Creating animation")
        def animation_loop(frame):
            pbar.update(1)
            pbar.set_postfix_str(f"Frame {frame+1}/{total_frames}")
            pbar.set_description(f"Drawing frame {frame+1}/{total_frames}")

            file_vr = files_vr[frame + start_frame]
            file_vr = os.path.join(folder, file_vr)
            file_pm = files_sol[frame + start_frame]
            file_pm = os.path.join(folder, file_pm)
            time = (frame + start_frame) * dt if dt else None
            self.load_results_from_disk(frame, total_frames, file_vr, file_pm, time)

        ani = FuncAnimation(
            self.fig, 
            animation_loop, 
            frames=total_frames, blit=False, repeat=False, interval=100,                
        )
        ani.event_source.stop()
        if filename:
            ani_name = f"{folder}/{filename}.{format}"
        else:
            ani_name = f"{folder}/animation_"
            for qoi in self.particle_quantities.values():
                ani_name += f"{qoi.name}_"
            for qoi in self.mesh_quantities.values():
                ani_name += f"{qoi.name}_"
            ani_name += f".{format}"
        ani.save(ani_name, writer='ffmpeg', fps=10, dpi=300)

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

        if self.has_particles:
            self._update_particle_plots(
                particle_positions,
                particle_strengths,
                particle_velocities,
                particle_deformations,
            )

        if self.has_mesh:
            self._update_mesh_plots(
                pm_positions,
                pm_strengths,
                pm_velocities,
                pm_deformations,
            )
        self._update_figure_title(iteration)
        self._render_plot()

    def _update_particle_plots(
        self,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        for key, qoi in self.particle_quantities.items():
            plot_data, data, info = self.data_adapters[key].transform(
                particle_positions, particle_strengths, particle_velocities, particle_deformations
            )
            title = f"Particle {qoi.quantity_name.capitalize()}"
            if qoi.component:
                title += f" {qoi.component.capitalize()}"
            self.plotters[key].update(plot_data, title)
            self.plotters[key].side_effect(plot_data, data, info)

    def _update_mesh_plots(
        self, 
        pm_positions: np.ndarray, 
        pm_strengths: np.ndarray, 
        pm_velocities: np.ndarray, 
        pm_deformations: np.ndarray
    ):
        for key, qoi in self.mesh_quantities.items():
            plot_data, data, info = self.data_adapters[key].transform(
                pm_positions, pm_strengths, pm_velocities, pm_deformations
            )
            title = f"Mesh {qoi.quantity_name.capitalize()}"
            if qoi.component:
                title += f" {qoi.component.capitalize()}"
            self.plotters[key].update(plot_data, title)
            self.plotters[key].side_effect(plot_data, data, info)

    def update_particle_plots(
        self,
        iteration: int,
        particle_positions: np.ndarray,
        particle_strengths: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        self._update_particle_plots(
            particle_positions,
            particle_strengths,
            particle_velocities,
            particle_deformations,
        )
        self._update_figure_title(iteration)
        self._render_plot() 

    def update_mesh_plots(
        self,
        iteration: int,
        pm_positions: np.ndarray,
        pm_strengths: np.ndarray,
        pm_velocities: np.ndarray,
        pm_deformations: np.ndarray,
    ):
        self._update_mesh_plots(
            pm_positions,
            pm_strengths,
            pm_velocities,
            pm_deformations,
        )
        self._update_figure_title(iteration)
        self._render_plot()

    def _update_figure_title(self, iteration, total_iterations = None, t = None):
        string = f"Frame {iteration}"
        if total_iterations:
            string += f" / {total_iterations}"
        if t:
            string += f", Realtime = {t:.2f}s"
        self.fig.suptitle(string)

    def _render_plot(self):
        self.fig.canvas.flush_events()
        # renderer = self.fig.canvas.get_renderer()
        # self.fig.draw(renderer)
        self.fig.canvas.draw()
        plt.pause(0.05)