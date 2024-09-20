import numpy as np

class Quantity:
    def __init__(self, name: str, num_dimensions: int, num_components: int):
        """
        Description of Quantity. For example, a quantity could be a vector field with 3 components 
        (e.g. velocity field) or a scalar field with 1 component (e.g. temperature field).

        Args:
            name (str): Name of the quantity.
            num_dimensions (int): Number of dimensions of the quantity (e.g. 2D or 3D).
            num_components (int): Number of components of the quantity.
        """
        self.quantity_name = name
        self.num_dimensions = num_dimensions
        self.num_components = num_components
    
    @property
    def name(self):
        return self.quantity_name

class QuantityOfInterest(Quantity):
    def __init__(self, name: str, num_dimensions: int, num_components: int, component: str | None):
        """
        Represents a specific quantity of interest for visualization.

        Args:
            name (str): Name of the quantity.
            num_dimensions (int): Number of dimensions of the quantity.
            num_components (int): Number of components of the quantity.
            component (str): Component of interest ('x', 'y', 'z', 'magnitude', or None for scalar fields).
        """
        super().__init__(name, num_dimensions, num_components)
        if num_components > 1 and component is None:
            raise ValueError("Component must be specified for vector fields.")
        self.field_type = 'vector' if num_components > 1 else 'scalar'
        self.component = component
    
    @property
    def name(self):
        return f"{self.quantity_name}_{self.component}" if self.component is not None else self.quantity_name

    def get_quantity(self, data: dict[str, np.ndarray]) -> np.ndarray:
        """
        Extracts the quantity of interest from the input data.

        Args:
            data (dict): Dictionary containing the field data.

        Returns:
            numpy.ndarray: The extracted quantity of interest.
        """
        if self.field_type == 'vector':
            if self.component in ['x', 'y', 'z']:
                return data[self.quantity_name][self.component]
            elif self.component == 'magnitude':
                data_np = np.array([data[self.quantity_name][c] for c in ['x', 'y', 'z']])
                return np.linalg.norm(data_np, axis=0)
            else:
                raise ValueError(f"Invalid component: {self.component}")
        elif self.field_type == 'scalar':
            return data[self.quantity_name]
        else: 
            raise ValueError(f"Invalid field type or component: {self.field_type}, {self.component}")

class ParticleQuantityOfInterest:
    @staticmethod
    def create_quantity_of_interest(name: str, num_components: int, component: str | None):
        return QuantityOfInterest(name, 1, num_components, component)
    
class MeshQuantityOfInterest:
    @staticmethod
    def create_quantity_of_interest(name: str, num_components: int, component: str):
        return QuantityOfInterest(name, 3, num_components, component)
    