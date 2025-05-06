import numpy as np
from abc import ABC, abstractmethod
from . import QuantityOfInterest, Slicer, SliceStrategy, Plane
from typing import Any, Literal

class Filter(ABC):
    """
    Abstract base class for filters used in data visualization.

    Methods
    -------
    apply(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]
        Abstract method to apply the filter to the given data. Must be implemented by subclasses.

    Parameters
    ----------
    data : dict[str, np.ndarray]
        A dictionary where keys are strings representing data labels and values are numpy arrays containing the data.

    Returns
    -------
    dict[str, np.ndarray]
        A dictionary with the same structure as the input, containing the filtered data.
    """
    @abstractmethod
    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray]
    ) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
        pass

class SliceFilter_3D(Filter):
    """
    A filter that slices a given quantity of interest on a specified plane using a defined slicing strategy.

    Attributes:
        quantity_of_interest (QuantityOfInterest): The quantity of interest to be sliced.
        slicer (Slicer): The slicer object that performs the slicing operation.

    Methods:
        __init__(quantity_of_interest, plane, slice_strategy):
            Initializes the SliceFilter with the given quantity of interest, plane, and slicing strategy.
        
        apply(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
            Applies the slice filter to the provided data and returns the sliced data.
            Args:
                data (dict[str, np.ndarray]): The input data containing positions and quantities.
            Returns:
                dict[str, np.ndarray]: The sliced data containing positions and colors.
    """
    def __init__(
        self,
        plane: Plane | str,
        strategy: SliceStrategy,
        filter_quantity: QuantityOfInterest | None = None,
        value: float | None = None,
    ):
        self.slicer = Slicer(strategy, plane, value)
        self.filter_quantity = filter_quantity

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        if self.filter_quantity:
            mesh_qoi = self.filter_quantity.get_quantity(data)
        else:
            mesh_qoi = quantity_of_interest.get_quantity(data)

        # Create slices
        slices = self.slicer.calculate_slicer(mesh_qoi, data['position'])

        # Extract slice data
        sliced_data = {}
        for key, slice in data.items():
            if slice is None:
                continue

            if key == 'position':
                continue
            if isinstance(slice, np.ndarray):
                sliced_data[key] = slice[slices]
            elif isinstance(slice, dict):
                sliced_data[key] = {k: v[slices] for k, v in slice.items()}
            elif key == 'solution':
                pass
            else:
                raise TypeError(f"Unsupported data type: {type(slice)}")
        # Reorder the missing axis out of the data. We want to put
        # The missing axis on the data['positions']['z'] axis
        # The functioning axis should be on the 'x' and 'y' axes
        x = data['position']['x']
        y = data['position']['y']
        z = data['position']['z']
        sliced_data['position'] = {}
        if self.slicer.plane == Plane.X.value:
            sliced_data['position']['x'] = y[slices]
            sliced_data['position']['y'] = z[slices]
            sliced_data['position']['z'] = x[slices]
        elif self.slicer.plane == Plane.Y.value:
            sliced_data['position']['x'] = x[slices]
            sliced_data['position']['y'] = z[slices]
            sliced_data['position']['z'] = y[slices]
        elif self.slicer.plane == Plane.Z.value:
            sliced_data['position']['x'] = x[slices]
            sliced_data['position']['y'] = y[slices]
            sliced_data['position']['z'] = z[slices]
        else:
            raise ValueError(f"Invalid plane: {self.slicer.plane}")
        
        return sliced_data, {'plane': self.slicer.plane, 'strategy': self.slicer.slice_strategy}

class ValueFilter(Filter):
    """
    A filter that selects data based on a specified value and tolerance for a given quantity of interest.

    Attributes:
        quantity_of_interest (QuantityOfInterest): The quantity of interest to filter on.
        tolerance (float): The tolerance within which the value should fall.
        value (float): The target value for the quantity of interest. Defaults to 0.

    Methods:
        apply(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
            Applies the filter to the provided data and returns the masked data.
    """
    def __init__(
        self,
        tolerance: float,
        value: float = 0,
        type = Literal['greater', 'less', 'equal'],
    ):
        self.value = value
        self.tolerance = tolerance
        self.type = type

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        quantity = quantity_of_interest.get_quantity(data)

        if self.type == 'greater':
            mask = (quantity - self.value) > self.tolerance
        elif self.type == 'less':
            mask = (quantity - self.value) < self.tolerance
        elif self.type == 'equal':
            mask = np.abs(quantity - self.value) < self.tolerance
        else:
            raise ValueError(f"Invalid type: {self.type}")
            
        
        masked_data = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                masked_data[key] = value[mask]
            elif isinstance(value, dict):
                masked_data[key] = {k: v[mask] for k, v in value.items()}
            else:
                raise TypeError(f"Unsupported data type: {type(value)}")
        
        return masked_data, {'value': self.value, 'tolerance': self.tolerance}

class PositionFilter(Filter):
    """
    A filter that selects data based on a specified position and tolerance for a given quantity of interest.

    Attributes:
        quantity_of_interest (QuantityOfInterest): The quantity of interest to filter on.
        tolerance (float): The tolerance within which the position should fall.
        position (np.ndarray): The target position for the quantity of interest.

    Methods:
        apply(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
            Applies the filter to the provided data and returns the masked data.
    """
    def __init__(
        self,
        type: Literal['greater', 'less', 'equal', 'within'],
        axis: int,
        position: float | np.ndarray,
        tolerance: float = 0.,
    ):
        self.position = position
        if type == 'greater':
            self.tolerance = np.inf
        elif type == 'less':
            self.tolerance = -np.inf
        elif type == 'equal':
            self.tolerance = 0
        elif type == 'within':
            self.tolerance = tolerance

        self.axis = axis

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        positions = data['position']
        positions = np.array([positions['x'], positions['y'], positions['z']])
        positions = positions[self.axis, :]
        if self.tolerance == np.inf:
            mask = positions > self.position
        elif self.tolerance == -np.inf:
            mask = positions < self.position
        elif self.tolerance == 0:
            mask = positions == self.position
        else:
            mask = np.abs(positions - self.position) < self.tolerance

        masked_data = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                masked_data[key] = value[mask]
            elif isinstance(value, dict):
                masked_data[key] = {k: v[mask] for k, v in value.items()}
            else:
                raise TypeError(f"Unsupported data type: {type(value)}")
        
        return masked_data, {'position': self.position, 'tolerance': self.tolerance}

class ValueSelector(Filter):
    """
    A filter that selects data based on a specified value and tolerance for a given quantity of interest.

    Attributes:
        quantity_of_interest (QuantityOfInterest): The quantity of interest to filter on.
        tolerance (float): The tolerance within which the value should fall.
        value (float): The target value for the quantity of interest. Defaults to 0.

    Methods:
        apply(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
            Applies the filter to the provided data and returns the masked data.
    """
    def __init__(
        self,
        type = Literal['top_num', 'bottom_num', 'top_percentage', 'bottom_percentage'],
        value: float = 0,
    ):
        self.value = value
        self.type = type

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        quantity = quantity_of_interest.get_quantity(data)

        # Top and bottom are used to indicate top % and bottom % of the data
        if self.type == 'top_percentage':
            mask = quantity > np.percentile(quantity, 100 - self.value)
        elif self.type == 'bottom_percentage':
            mask = quantity < np.percentile(quantity, self.value)
        # Max and min are used to indicate top and bottom number of data points
        elif self.type == 'top_num':
            val = min(self.value, len(quantity)- 1)
            mask = quantity > np.sort(quantity, axis= None)[-val]
            
        
        masked_data = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                masked_data[key] = value[mask]
            elif isinstance(value, dict):
                masked_data[key] = {k: v[mask] for k, v in value.items()}
            else:
                raise TypeError(f"Unsupported data type: {type(value)}")
        
        return masked_data, {'value': self.value, 'type': self.type}