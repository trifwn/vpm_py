import numpy as np
from abc import ABC, abstractmethod

class DataAdapter(ABC):
    @abstractmethod
    def get_particle_data(self):
        pass

    @abstractmethod
    def get_mesh_data(self):
        pass

    @abstractmethod
    def get_slice_data(self):
        pass

class VPMDataAdapter(DataAdapter):
    def __init__(self, particle_data, mesh_data):
        self.particle_data = particle_data
        self.mesh_data = mesh_data

    def get_particle_data(self):
        return {
            'positions': np.array([self.particle_data['XS'], 
                                   self.particle_data['YS'], 
                                   self.particle_data['ZS']]).T,
            'colors': self.particle_data[self.quantity_key]
        }

    def get_mesh_data(self):
        return {
            'positions': np.array([self.mesh_data['XS'], 
                                   self.mesh_data['YS'], 
                                   self.mesh_data['ZS']]),
            'colors': self.mesh_data[self.quantity_key]
        }

    def get_slice_data(self, axis='z'):
        if axis == 'z':
            max_z_index = np.argmax(np.sum(self.mesh_data[self.quantity_key], axis=(0,1)))
            return {
                'positions': (self.mesh_data['XS'][:,:,max_z_index], 
                              self.mesh_data['YS'][:,:,max_z_index]),
                'colors': self.mesh_data[self.quantity_key][:,:,max_z_index]
            }
        elif axis == 'y':
            mid_y = self.mesh_data['YS'].shape[1] // 2
            return {
                'positions': (self.mesh_data['XS'][:,mid_y,:], 
                              self.mesh_data['ZS'][:,mid_y,:]),
                'colors': self.mesh_data[self.quantity_key][:,mid_y,:]
            }
        else:
            raise ValueError("Invalid axis. Choose 'y' or 'z'.")

    @property
    def quantity_key(self):
        return 'QMAG'  # This could be made configurable