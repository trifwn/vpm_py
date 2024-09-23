import numpy as np
from abc import ABC, abstractmethod
from . import QuantityOfInterest, Slicer, SliceStrategy, Plane
from typing import Any

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
    ):
        self.slicer = Slicer(strategy, plane)

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        mesh_qoi = quantity_of_interest.get_quantity(data)

        # Create slices
        slices = self.slicer.calculate_slicer(mesh_qoi)

        # Extract slice data
        sliced_data = {}
        for key, slice in data.items():
            if key == 'position':
                continue
            if isinstance(slice, np.ndarray):
                sliced_data[key] = slice[slices]
            elif isinstance(slice, dict):
                sliced_data[key] = {k: v[slices] for k, v in slice.items()}
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
    ):
        self.value = value
        self.tolerance = tolerance

    def apply(
        self, 
        quantity_of_interest: QuantityOfInterest,
        data: dict[str, np.ndarray],
    ):
        quantity = quantity_of_interest.get_quantity(data)
        mask = (quantity - self.value) < self.tolerance
        
        masked_data = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                masked_data[key] = value[mask]
            elif isinstance(value, dict):
                masked_data[key] = {k: v[mask] for k, v in value.items()}
            else:
                raise TypeError(f"Unsupported data type: {type(value)}")
        
        return masked_data, {'value': self.value, 'tolerance': self.tolerance}