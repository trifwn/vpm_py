import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from vpm_py.vpm_io import print_green, print_IMPORTANT

class Particle3DPlot:
    """
    Class to plot the live position of particles in a 3D domain.
    """

    def __init__(self) -> None:
        """
        Initialize the particle plot.
        """
        print_IMPORTANT("Creating 3D plot")
        self.fig = plt.figure()
        self.ax: Axes3D = self.fig.add_subplot(111, projection="3d")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")
        self.ax.set_title("Particles in 3D")

        self.sc = self.ax.scatter([], [], [], c=[], marker="o")
        # Color bar
        self.fig.colorbar(self.sc, ax=self.ax, label="Strenghts")
        plt.ion()
        self.fig.show()

    def update(self, x, y, z, c):
        """
        Update the plot with new particle positions.

        Args:
        x: x-coordinates of the particles.
        y: y-coordinates of the particles.
        z: z-coordinates of the particles.
        c: color of the particles.
        """
        self.fig.canvas.flush_events()
        print_green(f"\tUpdating plot")
        print_green(f"\tNumber of particles: {len(x)}")
        self.sc.set_offsets(np.c_[x, y])
        self.sc.set_3d_properties(z, 'z')
        self.sc.set_array(c)
        # Relimit the axes
        self.ax.set_xlim(min(x), max(x))
        self.ax.set_ylim(min(y), max(y))
        self.ax.set_zlim(min(z), max(z))
        # Chage the colorbar limits
        self.sc.set_clim(min(c), max(c))
        
        # Redraw the plot
        self.ax.relim()
        self.ax.autoscale_view()
        # self.ax.autoscale()
        self.fig.canvas.draw()
        plt.pause(0.01)

    def close(self):
        """
        Close the plot.
        """
        plt.close()


