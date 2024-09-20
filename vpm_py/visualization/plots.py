from vpm_py.visualization import DataAdapterFactory, Filter, MeshPlotter, ParticlePlotter, SliceFilter_3D, SlicePlotter
from vpm_py.visualization.quantities import MeshQuantityOfInterest, ParticleQuantityOfInterest

from typing import Literal


class ResultPlot:
    """
    A class to represent a result plot for different types of data (particle or mesh).

    Attributes:
    -----------
    type : Literal['particle', 'mesh']
        The type of the plot, either 'particle' or 'mesh'.
    quantity : str
        The quantity to be plotted.
    component : str
        The component of the quantity to be plotted.
    filters : list[Filter], optional
        A list of filters to be applied to the data (default is an empty list).
    options : dict, optional
        Additional options for the plot (default is an empty dictionary).
    plotter : Plotter
        The plotter object to be used for plotting. Initialized when get_plotter() is called.
    data_adapter : DataAdapter
        The data adapter object to be used for data retrieval. Initialized when get_data_adapter() is called.
    quantity_of_interest : ParticleQuantityOfInterest or MeshQuantityOfInterest
        The quantity of interest object to be used for data retrieval. Initialized when get_quantity() is called.

    Methods:
    --------
    get_plotter(fig, ax):
        Returns the appropriate plotter object based on the type and filters.
    get_plot_dimensions():
        Returns the number of dimensions for the plot (2 or 3).
    get_data_adapter():
        Returns the appropriate data adapter based on the type.
    get_quantity():
        Returns the quantity of interest object based on the type.
    """
    def __init__(
        self,
        type: Literal['particle', 'mesh'],
        quantity: str,
        component: str,
        filters: list[Filter] = [],
        options: dict = {}
    ) -> None:
        """
        Initializes a new instance of the class.
        Args:
            type (Literal['particle', 'mesh']): The type of the plot, either 'particle' or 'mesh'.
            quantity (str): The quantity to be plotted.
            component (str): The component of the quantity to be plotted.
            filters (list[Filter], optional): A list of filters to apply to the data. Defaults to an empty list.
            options (dict, optional): Additional options for the plot. Defaults to an empty dictionary.
        Returns:
            None
        """

        self.type = type
        self.quantity = quantity
        self.component = component
        self.filters = filters
        self.options = options
        self.plotter = None
        self.data_adapter = None
        self.quantity_of_interest = None

    def get_plotter(self,fig, ax):
        """
        Returns an appropriate plotter instance based on the type of the object.

        Parameters:
        fig (matplotlib.figure.Figure): The figure object to plot on.
        ax (matplotlib.axes.Axes): The axes object to plot on.

        Returns:
        Plotter: An instance of a plotter class (ParticlePlotter, SlicePlotter, or MeshPlotter).

        Raises:
        ValueError: If the type of the object is not 'particle' or 'mesh'.
        """
        if self.type == 'particle':
            plotter =  ParticlePlotter(fig, ax, self.options)
        elif self.type == 'mesh':
            if any(isinstance(f, SliceFilter_3D) for f in self.filters):
                plotter = SlicePlotter(fig, ax, self.options)
            else:
                plotter = MeshPlotter(fig, ax, self.options)
        self.plotter = plotter
        return self.plotter

    def get_plot_dimensions(self):
        """
        Determine the plot dimensions based on the type and filters.

        Returns:
            int: The number of dimensions for the plot. Returns 3 if the type is 'particle'.
                 If the type is 'mesh', returns 2 if any filter is an instance of SliceFilter_3D,
                 otherwise returns 3.
        """
        if self.type == 'particle':
            return 3
        elif self.type == 'mesh':
            # Check if we have slice filters
            if any(isinstance(f, SliceFilter_3D) for f in self.filters):
                return 2
            else:
                return 3

    def get_data_adapter(self):
        """
        Retrieves the appropriate data adapter based on the type of data.

        Returns:
            DataAdapter: An instance of a data adapter for either 'particle' or 'mesh' data.

        Raises:
            ValueError: If the type is neither 'particle' nor 'mesh'.
        """
        if self.type == 'particle':
            data_adapter = DataAdapterFactory.create_adapter(
                'particles',
                ParticleQuantityOfInterest.create_quantity_of_interest(
                    self.quantity, 3, self.component
                ),
                self.filters
            )
        elif self.type == 'mesh':
            data_adapter = DataAdapterFactory.create_adapter(
                'mesh',
                MeshQuantityOfInterest.create_quantity_of_interest(
                    self.quantity, 3, self.component
                ),
                self.filters
            )
        self.data_adapter = data_adapter
        return self.data_adapter

    def get_quantity(self):
        """
        Retrieve the quantity of interest based on the type of object.

        Returns:
            ParticleQuantityOfInterest or MeshQuantityOfInterest: The quantity of interest object.

        Raises:
            ValueError: If the type is neither 'particle' nor 'mesh'.
        """
        if self.type == 'particle':
            qoi = ParticleQuantityOfInterest.create_quantity_of_interest(self.quantity, 3, self.component)
        elif self.type == 'mesh':
            qoi = MeshQuantityOfInterest.create_quantity_of_interest(self.quantity, 3, self.component)
        self.quantity_of_interest = qoi
        return self.quantity_of_interest

# Standard Plot Options
# Particle Options
particle_velocity_x_plot_option = ResultPlot('particle', 'velocity', 'x')
particle_velocity_y_plot_option = ResultPlot('particle', 'velocity', 'y')
particle_velocity_z_plot_option = ResultPlot('particle', 'velocity', 'z')
particle_vorticity_x_plot_option = ResultPlot('particle', 'strength', 'x')
particle_vorticity_y_plot_option = ResultPlot('particle', 'strength', 'y')
particle_vorticity_z_plot_option = ResultPlot('particle', 'strength', 'z')
particle_deformation_x_plot_option = ResultPlot('particle', 'deformation', 'x')
particle_deformation_y_plot_option = ResultPlot('particle', 'deformation', 'y')
particle_deformation_z_plot_option = ResultPlot('particle', 'deformation', 'z')
particle_velocity_magnitude_plot_option = ResultPlot('particle', 'velocity', 'magnitude')
particle_vorticity_magnitude_plot_option = ResultPlot('particle', 'strength', 'magnitude')
particle_deformation_magnitude_plot_option = ResultPlot('particle', 'deformation', 'magnitude')
# Mesh Options
mesh_velocity_x_plot_option = ResultPlot('mesh', 'velocity', 'x')
mesh_velocity_y_plot_option = ResultPlot('mesh', 'velocity', 'y')
mesh_velocity_z_plot_option = ResultPlot('mesh', 'velocity', 'z')
mesh_vorticity_x_plot_option = ResultPlot('mesh', 'strength', 'x')
mesh_vorticity_y_plot_option = ResultPlot('mesh', 'strength', 'y')
mesh_vorticity_z_plot_option = ResultPlot('mesh', 'strength', 'z')
mesh_deformation_x_plot_option = ResultPlot('mesh', 'deformation', 'x')
mesh_deformation_y_plot_option = ResultPlot('mesh', 'deformation', 'y')
mesh_deformation_z_plot_option = ResultPlot('mesh', 'deformation', 'z')
mesh_velocity_magnitude_plot_option = ResultPlot('mesh', 'velocity', 'magnitude')
mesh_vorticity_magnitude_plot_option = ResultPlot('mesh', 'strength', 'magnitude')
mesh_deformation_magnitude_plot_option = ResultPlot('mesh', 'deformation', 'magnitude')