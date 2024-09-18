import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from abc import ABC, abstractmethod

class Plotter(ABC):
    def __init__(self, fig, ax, title):
        self.fig = fig
        self.ax = ax
        self.title = title
        self.plot = None
        self.colorbar = None

    @abstractmethod
    def setup_plot(self):
        pass

    @abstractmethod
    def update_plot(self, data):
        pass

    def set_labels(self, xlabel, ylabel, zlabel=None):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        if zlabel:
            self.ax.set_zlabel(zlabel)

    def set_title(self, title):
        self.ax.set_title(title)

class ParticlePlotter(Plotter):
    def setup_plot(self):
        self.plot = self.ax.scatter([], [], [], c=[], s=0.2, cmap='viridis')
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)
        self.set_labels('X', 'Y', 'Z')
        self.set_title(self.title)

    def update_plot(self, data):
        x, y, z = data['positions'].T
        c = data['colors']
        self.plot._offsets3d = (x, y, z)
        self.plot.set_array(c)
        self.plot.set_clim(c.min(), c.max())
        self.ax.set_xlim(1.1 * x.min(), 1.1 * x.max())
        self.ax.set_ylim(1.1 * y.min(), 1.1 * y.max())
        self.ax.set_zlim(1.1 * z.min(), 1.1 * z.max())

class MeshPlotter(Plotter):
    def setup_plot(self):
        self.plot = self.ax.scatter([], [], [], c=[], s=0.2, cmap='viridis')
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)
        self.set_labels('X', 'Y', 'Z')
        self.set_title(self.title)

    def update_plot(self, data):
        x, y, z = data['positions']
        c = data['colors']
        self.plot._offsets3d = (x.ravel(), y.ravel(), z.ravel())
        self.plot.set_array(c.ravel())
        self.plot.set_clim(c.min(), c.max())
        self.ax.set_xlim(1.1 * x.min(), 1.1 * x.max())
        self.ax.set_ylim(1.1 * y.min(), 1.1 * y.max())
        self.ax.set_zlim(1.1 * z.min(), 1.1 * z.max())

class SlicePlotter(Plotter):
    def setup_plot(self):
        self.plot = self.ax.contourf(np.zeros((2,2)), np.zeros((2,2)), np.zeros((2,2)), cmap='viridis')
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)
        self.set_labels('X', 'Y')
        self.set_title(self.title)

    def update_plot(self, data):
        x, y = data['positions']
        z = data['colors']
        self.ax.clear()
        self.plot = self.ax.contourf(x, y, z, cmap='viridis')
        self.ax.set_title(self.title)
        self.set_labels('X', 'Y')
        if self.colorbar:
            self.colorbar.remove()
        self.colorbar = self.fig.colorbar(self.plot, ax=self.ax)