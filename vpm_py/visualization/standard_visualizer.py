from . import Visualizer, ResultPlot, SliceFilter_3D, SliceStrategy

class StandardVisualizer:
    def __new__(
        self,
        plot_particles: tuple[str, str] | None = None,
        plot_mesh: tuple[str, str] | None = None,
        plot_slices: tuple[str, str] | None = None,
        figure_size=(12, 10),
    ):
        
        plot_options = []
        if plot_particles:
            particle_option = ResultPlot(
                "particle",
                quantity= plot_particles[0],
                component= plot_particles[1],
                filters=[],
            )
            plot_options.append(particle_option)

        if plot_mesh:
            mesh_option = ResultPlot(
                "mesh",
                quantity=plot_mesh[0],
                component=plot_mesh[1],
                filters=[],
                options={
                    "s": "auto",
                },
            )
            plot_options.append(mesh_option)
        if plot_slices:
            z_slice = SliceFilter_3D(plane="Z",strategy=SliceStrategy.MAX_INTEGRAL)
            y_slice = SliceFilter_3D(plane="Y",strategy=SliceStrategy.MAX_INTEGRAL)
            plot_mesh_quantity_z = ResultPlot(
                "mesh",
                quantity=plot_slices[0],
                component=plot_slices[1],
                filters=[z_slice],
                options={
                    "add_slice_plane": plot_particles,
                },
            )
            plot_mesh_quantity_y = ResultPlot(
                "mesh",
                quantity=plot_slices[0],
                component=plot_slices[1],
                filters=[y_slice],
                options={
                    "add_slice_plane": plot_particles,
                },
            )
            plot_options.append(plot_mesh_quantity_z)
            plot_options.append(plot_mesh_quantity_y)


        return Visualizer(
            figure_size=figure_size,
            plot_options= plot_options
        )