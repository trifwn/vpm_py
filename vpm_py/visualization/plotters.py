import numpy as np
from abc import ABC, abstractmethod
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.artist import Artist
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from typing import Any, TYPE_CHECKING
from matplotlib.colorbar import Colorbar
from .annotate3D import Annotation3D

if TYPE_CHECKING:
    from vpm_py.visualization import ResultPlot

class Plotter(ABC):
    """
    Abstract base class for creating different types of plotters.
    """

    def __init__(
        self, 
        fig: Figure, 
        ax: Axes | Axes3D, 
        options: dict[str, Any] = {},
    ):
        """
        Initialize the Plotter.

        Parameters:
        fig (Figure): The matplotlib figure object.
        ax (Axes | Axes3D): The matplotlib axes object.
        options (dict, optional): Additional options for the plot. Defaults to an empty dictionary.
                    Each plotter class can have its own set of options. Uknown options are ignored.
        """
        self.fig = fig
        self.ax = ax
        self.plot = None
        self.colorbar = None
        self.is_3d = isinstance(ax, Axes3D)
        self.options = options
        self._title = self.ax.set_title('')
        self.setup_plot()

    @abstractmethod
    def update(self, *args, **kwargs):
        """
        Abstract method to update the plot with new data.
        """
        pass

    @abstractmethod
    def side_effect(self, *args, **kwargs):
        """
        Abstract method to add side effects to the plot.
        """
        pass

    @abstractmethod
    def setup_plot(self) -> None:
        """
        Abstract method to set up the initial plot.
        """
        pass

    @abstractmethod
    def get_artists(self) -> list[Artist]:
        """
        Get the artists for the plot.
        """
        return [self._title]

    def update_colorbar(self)-> None:
        """
        Update the colorbar based on the current plot.
        """
        if self.colorbar is None:
            self.colorbar: Colorbar = self.fig.colorbar(self.plot, ax=self.ax)
            self.colorbar.ax.tick_params(length=0, labelsize=8)  # Reduce tick size

        else:
            self.colorbar.update_normal(self.plot)

    def set_limits(self, x_lim, y_lim, z_lim=None):
        """
        Set the limits for the axes.

        Parameters:
        x_lim (tuple): The limits for the x-axis.
        y_lim (tuple): The limits for the y-axis.
        z_lim (tuple, optional): The limits for the z-axis (for 3D plots).
        """
        self.ax.set_xlim(x_lim)
        self.ax.set_ylim(y_lim)
        if z_lim and self.is_3d:
            self.ax.set_zlim(z_lim)

    def set_labels(self, x_label, y_label, z_label=None):
        """
        Set the labels for the axes.

        Parameters:
        x_label (str): The label for the x-axis.
        y_label (str): The label for the y-axis.
        z_label (str, optional): The label for the z-axis (for 3D plots).
        """
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        if z_label and self.is_3d:
            self.ax.set_zlabel(z_label)
    
    def set_title(self, title):
        """
        Set the title for the plot.

        Parameters:
        title (str): The title of the plot.
        """
        self._title.set_text(title)
        self._title.update({'fontsize': 10})

class ParticlePlotter(Plotter):
    """
    Plotter for particle data.
    """

    def setup_plot(self) -> None:
        """
        Set up the initial particle plot.
        """
        if 'cmap' in self.options.keys():
            self.cmap = self.options['cmap']
        else:
            self.cmap = 'viridis'
        if 'norm' in self.options.keys():
            self.norm = self.options['norm']
        else:
            self.norm = None
        if 's' in self.options.keys():
            self.s = self.options['s']
        else:
            self.s = 'auto'
        self.base_s = 10 if self.s == 'auto' else self.s

        if self.is_3d:
            self.plot = self.ax.scatter([], [], [], c=[], s=self.base_s, cmap=self.cmap, norm=self.norm, rasterized=True)
            self.set_labels('X', 'Y', 'Z')
        else:
            self.plot = self.ax.scatter([], [], c=[], s=self.base_s, cmap=self.cmap, norm=self.norm, rasterized=True)
            self.set_labels('X', 'Y')
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)
    
    def side_effect(self, *args, **kwargs):
        pass

    def get_artists(self):
        return super().get_artists() + [self.plot , self.colorbar]

    def update(self, plot_data, title=None):
        """
        Update the particle plot with new data.

        Parameters:
        data (dict): The data used to update the plot.
        title (str, optional): The title of the plot.
        """
        x = plot_data['positions']['x']
        y = plot_data['positions']['y']
        c = plot_data['colors']

        if self.is_3d:
            z = plot_data['positions']['z']
            self.plot._offsets3d = (x, y, z)
            # Set axis limits with equal scaling
            x_range = x.max() - x.min()
            y_range = y.max() - y.min()
            z_range = z.max() - z.min()
            max_range = max(x_range, y_range, z_range)

            x_mid = (x.max() + x.min()) * 0.5
            y_mid = (y.max() + y.min()) * 0.5
            z_mid = (z.max() + z.min()) * 0.5

            self.ax.set_xlim(x_mid - max_range/2, x_mid + max_range/2)
            self.ax.set_ylim(y_mid - max_range/2, y_mid + max_range/2)
            self.ax.set_zlim(z_mid - max_range/2, z_mid + max_range/2)

            # Set equal aspect ratio for 3D plots
            self.ax.set_box_aspect([1, 1, 1])
        else:
            self.plot.set_offsets(np.c_[x, y])
            x_range = x.max() - x.min()
            y_range = y.max() - y.min()
            max_range = max(x_range, y_range)

            x_mid = (x.max() + x.min()) * 0.5
            y_mid = (y.max() + y.min()) * 0.5

            self.ax.set_xlim(x_mid - max_range/2, x_mid + max_range/2)
            self.ax.set_ylim(y_mid - max_range/2, y_mid + max_range/2)

        self.plot.set_array(c)
        self.plot.set_clim(c.min(), c.max())
        if self.s == "auto":
            c = np.maximum(c, 1e-10)
            s = np.abs(c.ravel()/c.max())
            s = self.base_s * (1 + np.log10(s))
            # Set s to a minimum value 
            s = np.maximum(s, 0.01)
            self.plot.set_sizes(s.ravel())
        if title:
            self.set_title(title)
        self.update_colorbar()

class MeshPlotter(Plotter):
    """
    Plotter for mesh data.
    """

    def setup_plot(self) -> None:
        """
        Set up the initial mesh plot.
        """
        if 'cmap' in self.options.keys():
            cmap = self.options['cmap']
        else:
            cmap = 'viridis'
        if 'norm' in self.options.keys():
            norm = self.options['norm']
        else:
            norm = None
        if 's' in self.options.keys():
            self.s = self.options['s'] 

        self.base_s = 5
        if self.is_3d:
            self.plot = self.ax.scatter([], [], [], c=[], s=self.base_s, cmap=cmap, norm=norm, rasterized=True)
            self.set_labels('X', 'Y', 'Z')
        else:
            self.plot = self.ax.scatter([], [], c=[], s=self.base_s, cmap=cmap, norm=norm, rasterized=True)
            self.set_labels('X', 'Y')
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)
    
    def side_effect(self, *args, **kwargs):
        pass

    def get_artists(self):
        return super().get_artists() + [self.plot , self.colorbar]

    def update(self, plot_data, title=None):
        """
        Update the mesh plot with new data.

        Parameters:
        plot_data (dict): The data used to update the plot.
        title (str, optional): The title of the plot.
        """
        x = plot_data['positions']['x']
        y = plot_data['positions']['y']
        c = plot_data['colors']
        if self.is_3d:
            z = plot_data['positions']['z']
            self.plot._offsets3d = (x.ravel(), y.ravel(), z.ravel())
            self.ax.set_zlim(1.1 * z.min(), 1.1 * z.max())
        else:
            self.plot.set_offsets(np.c_[x.ravel(), y.ravel()])
        # Set the color and size of the points
        self.plot.set_array(c.ravel())
        self.plot.set_clim(c.min(), c.max())
        self.ax.set_xlim(1.1 * x.min(), 1.1 * x.max())
        self.ax.set_ylim(1.1 * y.min(), 1.1 * y.max())

        if self.s == "auto":
            c = np.maximum(c, 1e-10)
            s = np.abs(c.ravel()/c.max())
            s = self.base_s * (1 + np.log10(s))
            # Set s to a minimum value 
            s = np.maximum(s, 0.01)
            self.plot.set_sizes(s.ravel())
        
        if title:
            self.set_title(title)
        self.update_colorbar()
    
class SlicePlotter(Plotter):
    """
    Plotter for slice data.
    """

    def update_colorbar(self) -> None:
        """
        Update the colorbar for the slice plot.
        """
        if self.colorbar is None:
            divider = make_axes_locatable(self.ax)
            self.cax = divider.append_axes('right', size='5%', pad=0.05)
            self.colorbar = self.fig.colorbar(self.plot, cax=self.cax)
        else:
            self.colorbar.update_normal(self.plot)
    
    def setup_plot(self)-> None:
        """
        Set up the initial slice plot.
        """
        if 'cmap' in self.options.keys():
            self.cmap = self.options['cmap']
        else:
            self.cmap = 'viridis'

        if 'norm' in self.options.keys():
            self.norm = self.options['norm']
        else:
            self.norm = None

        if 'levels' in self.options.keys():
            self.levels = self.options['levels']
        else:
            self.levels = 100
        
        if 'add_slice_plane' in self.options.keys():
            self.slice_plane_plot: ResultPlot = self.options['add_slice_plane']
            self.slice_plane_ax = None
            self.slice_plane_surf = None
            self.annotation = None
        else:
            self.slice_plane_plot = None
            self.slice_plane_ax = None
            self.slice_plane_surf = None
            self.annotation = None

        self.plot = self.ax.pcolormesh(np.zeros((2,2)), np.zeros((2,2)), np.zeros((2,2)), cmap=self.cmap, norm=self.norm, shading='auto')
        self.ax.set_aspect('equal', 'box')
        self.set_labels('X', 'Y')
        self.update_colorbar() 

    def side_effect(self, plot_data, data, info):
        """
        Side effect method to add the slice plane to another plot. 
        Further side effects can be added here.
        """
        if self.slice_plane_ax:
            ax = self.slice_plane_ax

            # Add the slice plane to the other plot
            min_x, max_x = self.slice_plane_ax.get_xlim()
            min_y, max_y = self.slice_plane_ax.get_ylim()
            min_z, max_z = self.slice_plane_ax.get_zlim()

            if info['plane'] == 'X':
                y = data['position']['x']
                z = data['position']['y']
                x = data['position']['z']
                color = 'red'
            elif info['plane'] == 'Y':
                x = data['position']['x']
                z = data['position']['y']
                y = data['position']['z']
                color = 'green'
            elif info['plane'] == 'Z':
                x = data['position']['x']
                y = data['position']['y']
                z = data['position']['z']
                color = 'blue'
            else:
                print('Invalid plane.')
                return 

            if not self.annotation:
                self.annotation = Annotation3D(
                    s='Slice Plane', 
                    xyz=(min_x, min_y, min_z), 
                    xytext=(0, 0),  # Adjust text placement (negative x moves left)
                    textcoords='offset points',  # Use offset points for relative positioning
                    fontsize=10,
                    ha='right',  # Align text to the right, so it's outside the plot
                    va='center'
                )  
                # self.slice_plane_ax.add_artist(self.annotation)

            if np.all(x[0,0] == x):
                Y, Z = np.meshgrid(np.linspace(min_y, max_y, 3), np.linspace(min_z, max_z, 3))
                X = np.ones_like(Y) * x[0,0]
                # Name the axes appropriately
                self.ax.set_xlabel('Y')
                self.ax.set_ylabel('Z')
                title = f' Slice at X = {x[0,0]:.4f}' 
                # self.annotation.update(xyz=(x[0,0], min_y, min_z), s='X={:.2f}'.format(x[0,0]))
            elif np.all(y[0,0] == y):
                X, Z = np.meshgrid(np.linspace(min_x, max_x, 3), np.linspace(min_z, max_z, 3))
                Y = np.ones_like(X) * y[0,0]
                # Name the axes appropriately
                self.ax.set_xlabel('X')
                self.ax.set_ylabel('Z')
                title = f' Slice at Y = {y[0,0]:.4f}' 
                # self.annotation.update(xyz=(min_x, y[0,0], min_z), s='Y={:.2f}'.format(y[0,0]))
            elif np.all(z[0,0] == z):
                X, Y = np.meshgrid(np.linspace(min_x, max_x, 3), np.linspace(min_y, max_y, 3))
                Z = np.ones_like(X) * z[0,0]
                # Name the axes appropriately
                self.ax.set_xlabel('X')
                self.ax.set_ylabel('Y')
                title = f' Slice at Z = {z[0,0]:.4f}' 
                # self.annotation.update(xyz=(max_x, min_y, z[0,0]), s='Z={:.2f}'.format(z[0,0]))
            else:
                print('No fixed position found.')
            
            if self.slice_plane_surf:
                self.slice_plane_surf.remove()
                self.slice_plane_surf = ax.plot_surface(X,Y,Z, alpha=0.3, rstride=1, cstride=1, color=color)
            else:
                self.slice_plane_surf = ax.plot_surface(X,Y,Z, alpha=0.3, rstride=1, cstride=1, color=color)
            self.set_title(self.ax.get_title() + title)
        else:
            if self.slice_plane_plot:
                try:
                    plotter = self.slice_plane_plot.plotter
                    self.slice_plane_ax = plotter.ax
                    self.side_effect(plot_data, data, info)
                except AttributeError:
                    # print('No slice plot found. Cannot add side effect.')
                    return
            else:
                # print('No slice plot found.')
                return
        
    def get_artists(self):
        return super().get_artists() + [self.plot, self.colorbar, self.slice_plane_surf, self.annotation]
    
    def update(self, plot_data, title=None):
        """
        Update the slice plot with new data efficiently using pcolormesh.

        Parameters:
        plot_data (dict): The data used to update the plot.
        title (str, optional): The title of the plot.
        """
        x = plot_data['positions']['x']
        y = plot_data['positions']['y']
        z = plot_data['colors']

        # Check if the dimensions of x, y, or z have changed
        if hasattr(self, 'plot'):
            # Compare current dimensions with previous dimensions
            prev_x_shape = self.plot.get_array().shape if hasattr(self.plot, 'get_array') else None
            # print(f"Previous shape: {prev_x_shape} vs Current shape: {x.shape}")
            if prev_x_shape is not None and (x.shape != prev_x_shape or y.shape != prev_x_shape or z.shape != prev_x_shape):
                print('Grid Dimensions have changed, removing colormesh plot and creating a new one.')
                # Dimensions have changed, remove the old plot and create a new one
                self.plot.remove()
                self.plot = None

        # If the plot hasn't been initialized or dimensions have changed, create a new plot
        if not hasattr(self, 'plot') or self.plot is None:
            print(f"Creating new pcolormesh plot with dimensions: {x.shape} for {title}")
            self.plot = self.ax.pcolormesh(x, y, z, cmap=self.cmap, norm=self.norm, shading='auto')
            self.ax.set_aspect('equal', 'box')
            self.update_colorbar()
        else:
            # Update the data array and color limits
            self.plot.set_array(z.ravel())
            
        self.plot.set_clim(z.min(), z.max())
        # Adjust axes limits if the grid has changed
        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max()) 

        self.update_colorbar()
        # Update colorbar and title
        if title:
            self.set_title(title)