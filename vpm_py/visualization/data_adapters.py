import numpy as np
from abc import ABC, abstractmethod
from .filters import Filter 
from .quantities import QuantityOfInterest, Quantity

class DataAdapter(ABC):
    def __init__(
        self,
        quantity_of_interest: QuantityOfInterest,
        filters: list[Filter ] = []
    ):
        self.quantity_of_interest = quantity_of_interest
        self.filters = filters
    
    def apply_filters(self, data):
        filter_info = {}
        for filt in self.filters:
            data, info = filt.apply(self.quantity_of_interest, data)
            filter_info.update(info)
        return data , filter_info
    
    def transform(self, *args, **kwargs):
        data = self._preprocess(*args, **kwargs)
        filtered_data, filter_info = self.apply_filters(data)
        plot_data = self._transform(filtered_data)
        return plot_data, filtered_data, filter_info

    @abstractmethod
    def _preprocess(self, *args, **kwargs):
        pass

    @abstractmethod
    def _transform(self, data):
        pass

class ParticleDataAdapter(DataAdapter):
    def _preprocess(self, partcile_positions, particle_charges, particle_velocities, particle_deformations): 
        positions = {
            'x': partcile_positions[0,:],
            'y': partcile_positions[1,:],
            'z': partcile_positions[2,:]
        }
        charges = {
            'x': particle_charges[0,:],
            'y': particle_charges[1,:],
            'z': particle_charges[2,:]
        }
        velocities = {
            'x': particle_velocities[0,:],
            'y': particle_velocities[1,:],
            'z': particle_velocities[2,:]
        }
        deformations = {
            'x': particle_deformations[0,:],
            'y': particle_deformations[1,:],
            'z': particle_deformations[2,:]
        }
        data = {
            'position': positions,
            'charge': charges,
            'velocity': velocities,
            'deformation': deformations
        }
        return data

    def _transform(self, data):
        positions = data['position']
        colors = self.quantity_of_interest.get_quantity(data)
        return {
            'positions': positions,
            'colors': colors
        }

class MeshDataAdapter(DataAdapter):
    def _preprocess(
        self, 
        neq, 
        pm_positions, 
        pm_charges, 
        pm_velocities, 
        pm_deformations = None, 
        pm_solutions = None
    ):
        positions = {
            'x': pm_positions[0,:, :, :],
            'y': pm_positions[1,:, :, :],
            'z': pm_positions[2,:, :, :]
        }
        if pm_velocities is not None:
            velocities = {
                'x': pm_velocities[0,:, :, :],
                'y': pm_velocities[1,:, :, :],
                'z': pm_velocities[2,:, :, :]
            }
        else:
            velocities = None

        if pm_deformations is not None:
            deformations = {
                'x': pm_deformations[0,:, :, :],
                'y': pm_deformations[1,:, :, :],
                'z': pm_deformations[2,:, :, :]
            }
        else:
            deformations = None
        
        if pm_charges is not None:
            charges = {
                'x': pm_charges[0,:, : , :],
                'y': pm_charges[1,:, : , :],
                'z': pm_charges[2,:, : , :]
            }
        else:
            charges = None

        if pm_solutions is not None:
            solution = {
                'x': pm_solutions[0],
                'y': pm_solutions[1],
                'z': pm_solutions[2]
            }
        else:
            solution = None
            
        data = {
            'position': positions,
            'charge': charges,
            'velocity': velocities,
            'deformation': deformations,
            "solution": solution
        }
        return data
    
    def _transform(self, data):
        positions = data['position']
        colors = self.quantity_of_interest.get_quantity(data)
        return {
            'positions': positions,
            'colors': colors
        }
class DataAdapterFactory:
    @staticmethod
    def create_adapter(
        plot_type: str, 
        quantity_of_interest: QuantityOfInterest, 
        filters: list[Filter] = []
    ):
        if plot_type == 'particles':
            return ParticleDataAdapter(quantity_of_interest, filters)
        elif plot_type == 'mesh':
            return MeshDataAdapter(quantity_of_interest, filters)
        else:
            raise ValueError(f"Invalid plot type: {plot_type}")