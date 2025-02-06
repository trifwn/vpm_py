import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button
from matplotlib.artist import Artist
import numpy as np
import os
from tqdm import tqdm
import re
import time

from vpm_py.console_io import print_IMPORTANT
from vpm_py.file_io import process_particle_file, process_pm_file

from . import ResultPlot
from . import QuantityOfInterest
from . import Plotter
from . import DataAdapter
from . import OptimizedFFMpegWriter


class Visualizer:
    def __init__(
        self,
        figure_size: tuple[int, int] = (10, 10),
        plot_options: list[ResultPlot] = [],
        blitting: bool = False,
    ) -> None:
        """
        Initialize the particle plot.
        """
        print_IMPORTANT("Initializing the plot")
        self.plot_options: list[ResultPlot] = plot_options
        self.num_plots = len(self.plot_options) 
        print(f"Number of plots: {self.num_plots}")
        self.fig = plt.figure(figsize=figure_size) 
        self.plotters: dict[str, Plotter] = {}
        self.data_adapters: dict[str, DataAdapter] = {}

        self.mesh_quantities: dict[str, QuantityOfInterest] = {}
        self.particle_quantities: dict[str, QuantityOfInterest] = {}
        self.current_frame = 0

        self.has_mesh = False
        self.has_particles = False
        
        self._background = None
        self._setup_subplots()
        self.title = self.fig.suptitle("")

        self.blitting = blitting
        if blitting:
            self.cid = self.fig.canvas.mpl_connect("draw_event", self.on_draw)

    def _setup_subplots(self) -> None:
        self.rows = 1 if self.num_plots <= 2 else 2
        # We need to make sure that the number of columns is such that the number of plots is 
        # evenly distributed across the rows
        self.columns = int(np.ceil(self.num_plots / self.rows)) 
        print(f"Rows: {self.rows}, Columns: {self.columns}")
        gs = gridspec.GridSpec( self.rows, self.columns)

        for i, plot_option in enumerate(self.plot_options):
            row_num = i // self.columns
            col_num = i % self.columns

            if plot_option.get_plot_dimensions() == 2:
                ax = self.fig.add_subplot(gs[row_num, col_num], rasterized=True)
            else:
                ax = self.fig.add_subplot(gs[row_num, col_num], projection='3d', rasterized=True)
            
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
        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.9, hspace=0.3, wspace=0.3, bottom=0.1)

        # Set a box for the description on the bottom right corner 
        # self.description = self.fig.text(
        #     0.9, 0.1,
        #     '', 
        #     fontsize=10, 
        #     color='black', 
        #     bbox=dict(facecolor='navajowhite', alpha=0.5)
        # )


        self.fig.show()

    def set_problem_info(
        self,
        num_particles: int = 0,
        grid_size: tuple[int, int, int] = (0, 0, 0),
        dpm: tuple[float, float, float] = (0, 0, 0),
        dt: float | None = None,
    ):
        return
        s = f"PM Gird size: {grid_size} = {np.prod(grid_size)} cells\n"
        s += f"Number of particles: {num_particles}\n"
        s += f"Dpm = {dpm}\n"
        s += f"Dt  = {dt:.2f}" if dt else ""
        self.description.set_text(s)

    def get_artists(self) -> list[Artist]:
        all_artists = []
        for plotter in self.plotters.values():
            for artist in plotter.get_artists():
                all_artists.append(artist)
        return all_artists

    def add_folder(
        self,
        folder: str,
        particle_filename_pattern: str = r'particles.*h5',
        mesh_filename_pattern: str = r'particle_mesh.*h5',
        show_animation_buttons: bool = True,
        show_folder_name: bool = True,
    ):
        """
        Adds a folder to the visualizer to animate the results contained in it.
        Adds buttons to the plot to control the animation.

        Args:
            folder (str): Folder containing the results.
            particle_filename_pattern (str, optional): Regex pattern to match particle files. Defaults to r'particles.*\.h5'.
            mesh_filename_pattern (str, optional): Regex pattern to match mesh files. Defaults to r'pm_output.*\.h5'.
            show_animation_buttons (bool, optional): Whether to show the animation buttons. Defaults to True.
            show_folder_name (bool, optional): Whether to show the folder name in the title. Defaults to True.
        """
        # Compile the regex patterns for particles and mesh filenames
        particle_regex = re.compile(particle_filename_pattern)
        mesh_regex = re.compile(mesh_filename_pattern)
        
        # List all files in the folder
        files = os.listdir(folder)
        
        # Filter and sort files based on the regex pattern for particle and mesh files
        files_vr = sorted([f for f in files if particle_regex.search(f)])
        files_sol = sorted([f for f in files if mesh_regex.search(f)])
        print(f"Adding folder: {folder}")
        print(f"Number of particle files: {len(files_vr)}")
        print(f"Number of mesh files: {len(files_sol)}")

        # Get the folder name
        if show_folder_name:
            folder_name = os.path.basename(folder)
        else:
            folder_name = None
        
        # Add the animation buttons to the plot and a callback to control the animation
        if show_animation_buttons:
            self._add_animation_controls(
                folder,
                files_vr,
                files_sol,
                folder_name if show_folder_name else None,
            )
    
    def _add_animation_controls(
        self,
        folder: str,
        files_vr: list[str],
        files_sol: list[str],
        folder_name: str | None,
    ):
        """
        Add animation buttons to the plot to control the animation.

        Args:
            folder (str): Folder containing the results.
            files_vr (list[str]): List of particle files.
            files_sol (list[str]): List of mesh files.
            folder_name (str | None): Folder name to show in the title.
        """
        def update_frame():
            print("update frame")
            frame = self.current_frame
            # Load the results from disk
            file_vr = files_vr[frame]
            file_vr = os.path.join(folder, file_vr)
            file_pm = files_sol[frame]
            file_pm = os.path.join(folder, file_pm)
            self.load_results_from_disk(frame, len(files_vr), file_vr, file_pm)
            # Update the title with the current frame
            self._update_figure_title(f"Frame {frame}, Realtime = {frame* 0.1:.1f}s")
            self.fig.canvas.draw()

        def next_frame(event):
            print("next frame")
            if self.current_frame < len(files_vr) - 1:
                self.current_frame += 1
            else:
                self.current_frame = 0  # Loop back to the first frame
            update_frame()

        def prev_frame(event):
            print("prev frame")
            if self.current_frame > 0:
                self.current_frame -= 1
            else:
                self.current_frame = len(files_vr) - 1  # Loop back to the last frame
            update_frame()

        def on_key_press(event):
            """Handle keyboard events to control frames."""
            print(f"Key pressed: {event.key}")
            if event.key == 'right':  # Move forward
                self.on_next(None)
            elif event.key == 'left':  # Move backward
                self.on_prev(None)

        # Create the buttons for next/previous frame. Buttons should be placed at the bottom of the plot
        # On the left and be relatively small
        axprev = plt.axes((0.1, 0.01, 0.1, 0.05))
        axnext = plt.axes((0.8, 0.01, 0.1, 0.05))
        bnext = Button(axnext, '+')
        bprev = Button(axprev, '-')
        # Connect the buttons to their callback functions
        bnext.on_clicked(next_frame)
        bprev.on_clicked(prev_frame)
        # Add buttons to self.fig
        self.fig.canvas.draw()
        self.fig.canvas.mpl_connect('key_press_event', on_key_press)

        # Set the initial title of the plot
        if folder_name:
            self._update_figure_title(f"{folder_name} - Frame {self.current_frame}, Realtime = 0.0s")
        else:
            self._update_figure_title(f"Frame {self.current_frame}, Realtime = 0.0s")

    def load_results_from_disk(
        self, 
        frame: int, 
        total_frames: int,
        file_particles: str,
        file_particle_mesh: str,
        time: float | None = None,
        render: bool = True
    ):
        if self.has_particles:
            (
                particle_positions, particle_velocities, particle_charges, particle_deformations
            ) = process_particle_file(file_particles)
            self._update_particle_plots(
                particle_positions,
                particle_charges,
                particle_velocities,
                particle_deformations,
            )
        
        if self.has_mesh:
            (
                neq, mesh_positions, mesh_velocities, mesh_charges, mesh_vortex_stretching, 
                mesh_solutions, mesh_pressure, mesh_q_pressure, mesh_u_pressure
            ) = process_pm_file(file_particle_mesh)
            self._update_mesh_plots(
                mesh_positions,
                mesh_charges,
                mesh_velocities,
                mesh_vortex_stretching,
                mesh_solutions,
                mesh_pressure,
                mesh_q_pressure,
                mesh_u_pressure,
                neq
            )
        title = f"Frame {frame}, Realtime = {frame* 0.1:.1f}s"
        if time:
            title += f", Time = {time:.2f}s"
        self._update_figure_title(title)
        
        if render:
            print(f"Rendering plot for frame {frame + 1}/{total_frames}")
            self._render_plot()

    def setup_animation_writer(
        self, 
        filename: str, 
        fps: int = 10,
        codec: str = 'libx264',
        bitrate: int = 1800,
        dpi: int = 200
    ):
        """
        Set up the writer for the animation.
        """
        self.writer = OptimizedFFMpegWriter(
            fps=fps,
            codec=codec,
            bitrate=bitrate,
            extra_args=['-preset', 'fast', '-pix_fmt', 'yuv420p']
        )
        self.writer.setup(self.fig, outfile = filename, dpi = dpi)
        self.fig.set_dpi(dpi)
    
    def grab_frame(self):
        """
        Capture the current frame of the animation.
        """
        s_time = time.time()
        self.writer.grab_frame()
        e_time = time.time()
        # Log frame time for monitoring
        print(f"\tFrame grabbed in {e_time - s_time:.2f}s")

    
    def finish_animation(self):
        """
        Finish the animation and save the file.
        """
        self.writer.finish()

    def animate_result_folder(
        self, 
        folder: str, 
        total_frames: int,
        start_frame: int = 0,
        save_filename: str | None = None,
        particle_filename_pattern: str = r'particles.*\.h5',
        mesh_filename_pattern: str = r'pm_output.*\.h5',
        dt: float | None = None,
        format: str = 'mp4'
    ):
        """
        Animate the results from a specified folder containing particle and mesh files.
        Parameters:
        -----------
        folder : str
            The directory containing the result files.
        total_frames : int
            The total number of frames to animate.
        start_frame : int, optional
            The starting frame for the animation (default is 0).
        save_filename : str or None, optional
            The filename to save the animation. If None, a default name is generated (default is None).
        particle_filename_pattern : str, optional
            The regex pattern to match particle files (default is r'particles.*\.h5').
        mesh_filename_pattern : str, optional
            The regex pattern to match mesh files (default is r'pm_output.*\.h5').
        dt : float or None, optional
            The time step between frames. If None, time is not considered (default is None).
        format : str, optional
            The format to save the animation (default is 'mp4').
        Returns:
        --------
        None
        """
        # Compile the regex patterns for particles and mesh filenames
        particle_regex = re.compile(particle_filename_pattern)
        mesh_regex = re.compile(mesh_filename_pattern)
        
        # List all files in the folder
        files = os.listdir(folder)
        
        # Filter and sort files based on the regex pattern for particle and mesh files
        files_vr = sorted([f for f in files if particle_regex.search(f)])
        files_sol = sorted([f for f in files if mesh_regex.search(f)])

        pbar = tqdm(total=total_frames, desc="Creating animation")
        def animation_loop(frame):
            pbar.update(1)
            pbar.set_postfix_str(f"Frame {frame +1}/{total_frames - start_frame}")
            pbar.set_description(f"Drawing frame {frame +1}/{total_frames - start_frame}")

            file_vr = files_vr[frame + start_frame]
            file_vr = os.path.join(folder, file_vr)
            file_pm = files_sol[frame + start_frame]
            file_pm = os.path.join(folder, file_pm)
            time = (frame + start_frame) * dt if dt else None
            self.load_results_from_disk(frame + start_frame, total_frames, file_vr, file_pm, time)

        ani = FuncAnimation(
            self.fig, 
            animation_loop, 
            frames=total_frames, blit=False, repeat=False, interval=100,                
        )
        ani.event_source.stop()
        if save_filename:
            ani_name = f"{folder}/{save_filename}.{format}"
        else:
            ani_name = f"{folder}/animation_"
            for qoi in self.particle_quantities.values():
                ani_name += f"{qoi.name}_"
            for qoi in self.mesh_quantities.values():
                ani_name += f"{qoi.name}_"
            ani_name += f".{format}"
        ani.save(ani_name, writer='ffmpeg', fps=10, dpi=200)

    def update_all_plots(
        self,
        title: str,
        particle_positions: np.ndarray,
        particle_charges: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
        pm_positions: np.ndarray,
        pm_charges: np.ndarray,
        pm_velocities: np.ndarray,
        pm_vortex_stretching: np.ndarray,
        pm_solutions: np.ndarray | None = None,
        pm_pressure: np.ndarray | None = None,
        pm_q_pressure: np.ndarray | None = None,
        pm_u_pressure: np.ndarray | None = None,
        neq: int = 3,
    ):
        if self._background:
            s_time = time.time()
            self.fig.canvas.restore_region(self._background)
            print(f"\tBackground restored in {time.time() - s_time:.2f}s")
        
        # Restore background
        if self.has_particles:
            s_time = time.time()
            self._update_particle_plots(
                particle_positions,
                particle_charges,
                particle_velocities,
                particle_deformations,
            )
            print(f"\tParticle plots updated in {time.time() - s_time:.2f}s")

        if self.has_mesh:
            s_time = time.time()
            self._update_mesh_plots(
                pm_positions,
                pm_charges,
                pm_velocities,
                pm_vortex_stretching,
                pm_solutions,
                pm_pressure,
                pm_q_pressure,
                pm_u_pressure,
                neq
            )
            print(f"\tMesh plots updated in {time.time() - s_time:.2f}s")

        s_time = time.time()
        self._update_figure_title(title)
        print(f"\tTitle updated in {time.time() - s_time:.2f}s")
        
        self._render_plot()

    def _update_particle_plots(
        self,
        particle_positions: np.ndarray,
        particle_charges: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        for key, qoi in self.particle_quantities.items():
            plot_data, data, info = self.data_adapters[key].transform(
                particle_positions, particle_charges, particle_velocities, particle_deformations
            )
            title = f"Particle {qoi.quantity_name.capitalize()}"
            if qoi.component:
                title += f" {qoi.component.capitalize()}"
            self.plotters[key].update(plot_data, title)
            self.plotters[key].side_effect(plot_data, data, info)

    def _update_mesh_plots(
        self, 
        pm_positions: np.ndarray, 
        pm_charges: np.ndarray, 
        pm_velocities: np.ndarray, 
        pm_vortex_stretching: np.ndarray | None,
        pm_solutions: np.ndarray | None = None,
        pm_pressure: np.ndarray | None = None,
        pm_q_pressure: np.ndarray | None = None,
        pm_u_pressure: np.ndarray | None = None,
        neq: int = 3,
    ):
        for key, qoi in self.mesh_quantities.items():
            plot_data, data, info = self.data_adapters[key].transform(
                neq, pm_positions, pm_charges, pm_velocities, pm_vortex_stretching, 
                pm_solutions, pm_pressure, pm_q_pressure, pm_u_pressure
            )
            title = f"Mesh {qoi.quantity_name.capitalize()}"
            if qoi.component:
                title += f" {qoi.component.capitalize()}"
            self.plotters[key].update(plot_data, title)
            self.plotters[key].side_effect(plot_data, data, info)

    def update_particle_plots(
        self,
        title: str,
        particle_positions: np.ndarray,
        particle_charges: np.ndarray,
        particle_velocities: np.ndarray,
        particle_deformations: np.ndarray,
    ):
        self._update_particle_plots(
            particle_positions,
            particle_charges,
            particle_velocities,
            particle_deformations,
        )
        self._update_figure_title(title)
        self._render_plot() 

    def update_mesh_plots(
        self,
        title: str,
        pm_positions: np.ndarray,
        pm_charges: np.ndarray,
        pm_velocities: np.ndarray,
        pm_vortex_stretching: np.ndarray | None = None,
        pm_solutions: np.ndarray | None = None,
        pm_pressure: np.ndarray | None = None,
        pm_q_pressure: np.ndarray | None = None,
        pm_u_pressure: np.ndarray | None = None,
        neq: int = 3,
    ):
        self._update_mesh_plots(
            pm_positions = pm_positions,
            pm_charges = pm_charges,
            pm_velocities = pm_velocities,
            pm_vortex_stretching = pm_vortex_stretching,
            pm_solutions = pm_solutions,
            pm_pressure = pm_pressure,
            pm_q_pressure = pm_q_pressure,
            pm_u_pressure = pm_u_pressure,
            neq = neq,
        )
        self._update_figure_title(title)
        self._render_plot()

    def _update_figure_title(self, title:str)-> None:
        self.title.set_text(title)

    def _draw_animated(self):
        """Draw all of the animated artists."""
        fig = self.fig
        for artist in self.get_artists():
            try:
                # Check if artist is animated
                if not artist.get_animated():
                    artist.set_animated(True)
                fig.draw_artist(artist)
            except (ValueError, AttributeError):
                print(f"Error drawing artist: {artist}")

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        canvas = self.fig.canvas
        if event is not None:
            if event.canvas != canvas:
                raise RuntimeError
        self._background = canvas.copy_from_bbox(canvas.figure.bbox)
        self._draw_animated()

    def _render_plot(self):
        fig = self.fig
        canvas = fig.canvas

        if self.blitting:
            if self._background is None:
                timer = time.time()
                self.on_draw(None)
                print(f"Background created in {time.time() - timer:.2f}s")

            else:
                timer = time.time()
                canvas.restore_region(self._background)
                print(f"Background restored in {time.time() - timer:.2f}s")

                timer = time.time()
                self._draw_animated()
                print(f"\tArtists drawn in {time.time() - timer:.2f}s")

                timer = time.time()
                canvas.blit(fig.bbox)
                print(f"\tPlot blitted in {time.time() - timer:.2f}s")
        else:
            timer = time.time()
            canvas.draw_idle()
            print(f"\tPlot drawn in {time.time() - timer:.2f}s")

        timer = time.time()
        canvas.flush_events()  # Ensures the GUI updates immediately without blocking.
        print(f"\tEvents flushed in {time.time() - timer:.2f}s")

    def __del__(self):
        plt.close(self.fig)
        if hasattr(self, 'writer'):
            self.finish_animation()