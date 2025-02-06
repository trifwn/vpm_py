import numpy as np
from enum import Enum
from scipy.integrate import simpson as simps

class SliceStrategy(Enum):
    MIDDLE = "middle" 
    MAX_INTEGRAL = "max_integral"
    POSITION = "position"

    def __eq__(self, other):
        if isinstance(other, str):
            return self.value == other
        elif isinstance(other, SliceStrategy):
            return self.value == other.value
        else:
            return False

class Plane(Enum):
    X = "X" 
    Y = "Y" 
    Z = "Z"

    def __eq__(self, other):
        if isinstance(other, str):
            return self.value == other
        elif isinstance(other, Plane):
            return self.value == other.value
        else:
            return False

class Slicer:
    def __init__(
        self,
        slice_strategy: SliceStrategy | str,
        plane: Plane | str,
        value: float | None = None,
    ):
        if isinstance(slice_strategy, str):
            slice_strategy = SliceStrategy(slice_strategy)
        self.slice_strategy = slice_strategy
        if isinstance(plane, str):
            plane = Plane(plane)
        self.plane = plane
        self.value = value

    def calculate_slicer(
        self, 
        mesh_quantity: np.ndarray,
        positions: dict[str, np.ndarray], 
    ) -> tuple[slice | int, slice | int, slice | int]:
        if self.slice_strategy == SliceStrategy.MIDDLE:
            slice_idx = self.calculate_middle_slice(mesh_quantity)
        elif self.slice_strategy == SliceStrategy.MAX_INTEGRAL:
            slice_idx = self.calculate_max_integral_slice(mesh_quantity)
        elif self.slice_strategy == SliceStrategy.POSITION:
            slice_idx = self.calculate_idx_from_position(positions)
        else:
            raise ValueError(f"Invalid slice strategy: {self.slice_strategy}")
        
        if self.plane == Plane.X:
            return (slice_idx, slice(None), slice(None))
        elif self.plane == Plane.Y:
            return (slice(None), slice_idx, slice(None))
        elif self.plane == Plane.Z:
            return (slice(None), slice(None), slice_idx)
        else:
            raise ValueError(f"Invalid plane: {self.plane}")

    def calculate_middle_slice(
        self, 
        mesh_quantity: np.ndarray,
    ) -> int:
        if self.plane == Plane.X:
            return mesh_quantity.shape[0] // 2
        elif self.plane == Plane.Y:
            return mesh_quantity.shape[1] // 2
        elif self.plane == Plane.Z:
            return mesh_quantity.shape[2] // 2
        else:
            raise ValueError(f"Invalid plane: {self.plane}")

    def calculate_max_integral_slice(
        self, 
        mesh_quantity: np.ndarray,
    ) -> int:
        if self.plane == Plane.X:
            integral = simps( simps( np.abs(mesh_quantity[:, :, :]), axis=2), axis=1)
        elif self.plane == Plane.Y:
            integral = simps( simps( np.abs(mesh_quantity[:, :, :]), axis=2), axis=0)
        elif self.plane == Plane.Z:
            integral = simps( simps( np.abs(mesh_quantity[:, :, :]), axis=1), axis=0)
        else:
            raise ValueError(f"Invalid plane: {self.plane}")
        return int(np.argmax(integral))

    def calculate_idx_from_position(
        self, 
        positions: dict[str, np.ndarray],
    ) -> int:
        """
        Calculate the slice index from the given position
        """ 
        if self.plane == Plane.X:
            positions = positions['x']
            slice_idx = np.argmin(np.abs(positions[:, 0, 0] - self.value))
        elif self.plane == Plane.Y:
            positions = positions['y']
            slice_idx = np.argmin(np.abs(positions[0, :, 0] - self.value))
        elif self.plane == Plane.Z:
            positions = positions['z']
            slice_idx = np.argmin(np.abs(positions[0, 0, :] - self.value))
        else:
            raise ValueError(f"Invalid plane: {self.plane}")
        return int(slice_idx)
