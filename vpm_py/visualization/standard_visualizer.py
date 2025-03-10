from . import Visualizer, ResultPlot, SliceFilter_3D, SliceStrategy, ValueSelector, MeshQuantityOfInterest, PositionFilter

class StandardVisualizer(Visualizer):
    def __new__(
        self,
        plot_particles: tuple[str, str] | None = None,
        plot_mesh: tuple[str, str] | None = None,
        plot_slices: list[tuple[str, str]] | None = None,
        figure_size=(12, 10),
    ):
        
        plot_options = []
        if plot_particles:
            particle_option = ResultPlot(
                "particle",
                quantity= plot_particles[0],
                component= plot_particles[1],
                filters=[
                    # ValueSelector('top_num',  50000),
                    # PositionFilter('greater', axis= 1, position = 0.0),
                ],
                options={
                    "s": "auto",
                },
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
            for plot_slice in plot_slices:
                qoi = MeshQuantityOfInterest.create_quantity_of_interest("charge", 3, "magnitude")

                z_slice = SliceFilter_3D(plane="Z",strategy=SliceStrategy.MAX_INTEGRAL, filter_quantity= qoi )
                y_slice = SliceFilter_3D(plane="Y",strategy=SliceStrategy.POSITION, filter_quantity= qoi, value= 0.0)
                if plot_particles:
                    options = {
                        "add_slice_plane": particle_option,
                    }
                else:
                    options = {}

                plot_mesh_quantity_z = ResultPlot(
                    "mesh",
                    quantity=plot_slice[0],
                    component=plot_slice[1],
                    filters=[z_slice],
                    options= options,
                )
                plot_mesh_quantity_y = ResultPlot(
                    "mesh",
                    quantity=plot_slice[0],
                    component=plot_slice[1],
                    filters=[y_slice],
                    options= options,
                )
                plot_options.append(plot_mesh_quantity_z)
                plot_options.append(plot_mesh_quantity_y)

        return Visualizer(
            figure_size=figure_size,
            plot_options= plot_options
        )