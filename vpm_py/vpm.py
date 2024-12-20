import os
from ctypes import (_Pointer, byref,  c_double, c_int)

import numpy as np

# Local imports
from . import ParticleMesh, Particles
from .arrays import F_Array, F_Array_Struct
from .console_io import print_blue, print_green, print_IMPORTANT, print_red
from .utils import divide_processors
from .vpm_dtypes import dp_array_to_pointer, pointer_to_dp_array
from .vpm_lib import VPM_Lib
from .visualization import Visualizer


class VPM(object):
    """
    Interface to the VPM Fortran routines.
    """

    ### Initialization ###
    def __init__(
        self,
        number_of_equations: int = 3,
        number_of_processors: int = 1,
        rank: int = 0,
        verbocity: int = 1,
        dx_particle_mesh: float = 0.2,
        dy_particle_mesh: float = 0.2,
        dz_particle_mesh: float = 0.2,
        case_folder: str | None = None,
    ):
        lib = VPM_Lib()
        self._lib = lib._lib_vpm
        self.vpm_lib = lib

        self.rank = rank
        self.num_processors = number_of_processors
        
        # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
        NBI, NBJ, NBK = divide_processors(number_of_processors)
        if self.rank == 0:
            print(f"Number of processors: {self.num_processors}")
            print(f"NBI: {NBI}, NBJ: {NBJ}, NBK: {NBK}")
        
        self.dpm = np.array([dx_particle_mesh, dy_particle_mesh, dz_particle_mesh])
        self.NBI = NBI
        self.NBJ = NBJ
        self.NBK = NBK
        self.num_equations = number_of_equations

        self.set_verbosity(verbocity)
        if case_folder is not None:
            self.set_case_folder(case_folder)
        # Initialize the VPM
        self.initialize(self.dpm[0], self.dpm[1], self.dpm[2], NBI, NBJ, NBK)
        self.particle_mesh = ParticleMesh()
        self.particles = Particles(number_equations= self.num_equations)

        # Visualization
        self.visualizer: Visualizer | None = None
        self.has_animation_writer = False

    def initialize(
        self, 
        DXpm: float, 
        DYpm: float, 
        DZpm: float, 
        NBI: int, 
        NBJ: int, 
        NBK: int, 
        ncoarse = 8, 
        projection_type: int = 4, 
        boundary_condition: int =2, 
        variable_volume: bool = True, 
        remesh: bool = True, 
        use_tree: bool = True, 
        ilevmax: int = 4,
        OMPTHREADS: int = 1,
        is_domain_fixed: bool = False, 
        IPMWRITE: int = 0, 
        IPMWSTART: int = -1, 
        IPMWSTEPS: int = -1,
    ):
        """Wrapper function for calling Fortran subroutine `init`.

        Args:
            DXpm (float): Grid spacing in x-direction.
            DYpm (float): Grid spacing in y-direction.
            DZpm (float): Grid spacing in z-direction.
            NBI (int): Number of cells in x-direction.
            NBJ (int): Number of cells in y-direction.
            NBK (int): Number of cells in z-direction.
            ncoarse (int, optional): Number of coarse cells. Defaults to 8.
            projection_type (int, optional): Projection type. Defaults to 4.
            boundary_condition (int, optional): Boundary condition. Defaults to 2.
            variable_volume (bool, optional): Variable volume. Defaults to True.
            remesh (bool, optional): Remesh. Defaults to True.
            use_tree (bool, optional): Use tree. Defaults to True.
            ilevmax (int, optional): Maximum level. Defaults to 4.
            OMPTHREADS (int, optional): Number of threads. Defaults to 1.
            is_domain_fixed (bool, optional): Is box fixed. Defaults to True.
            IPMWRITE (int, optional): Write IPM. Defaults to 1.
            IPMWSTART (int, optional): Write IPM start. Defaults to 1.
            IPMWSTEPS (int, optional): Write IPM steps. Defaults to 1.
        """
        self.dpm = np.array([DXpm, DYpm, DZpm])
        self._lib.init(
            byref(c_double(DXpm)), byref(c_double(DYpm)), byref(c_double(DZpm)), byref(c_int(projection_type)),
            byref(c_int(boundary_condition)), byref(c_int(variable_volume)), byref(c_int(ncoarse)),
            byref(c_int(NBI)), byref(c_int(NBJ)), byref(c_int(NBK)), byref(c_int(remesh)), 
            byref(c_int(use_tree)), byref(c_int(ilevmax)), byref(c_int(OMPTHREADS)), byref(c_int(is_domain_fixed)),
            byref(c_int(IPMWRITE)), byref(c_int(IPMWSTART)), byref(c_int(IPMWSTEPS)),
        )
        if(self.rank == 0):
            print_green(f"Finished initializing VPM {self.rank}:")
            # Print the arguments passed
            print(f"\tDXpm= {DXpm}")
            print(f"\tDYpm= {DYpm}")
            print(f"\tDZpm= {DZpm}")
            print(f"\tinterf_iproj= {projection_type}")
            print(f"\tibctyp= {boundary_condition}")
            print(f"\tIDVPM= {variable_volume}")
            print(f"\tncoarse= {ncoarse}")
            print(f"\tNBI= {NBI}")
            print(f"\tNBJ= {NBJ}")
            print(f"\tNBK= {NBK}")
            print(f"\tNREMESH= {remesh}")
            print(f"\tiyntree= {use_tree}")
            print(f"\tilevmax= {ilevmax}")
            print(f"\tOMPTHREADS= {OMPTHREADS}")
            print(f"\tidefine= {is_domain_fixed}")
            print(f"\tIPMWRITE= {IPMWRITE}")
            print(f"\tIPMWSTART= {IPMWSTART}")
            print(f"\tIPMWSTEPS= {IPMWSTEPS}")

    def attach_visualizer(self, plotter: Visualizer):
        if not isinstance(plotter, Visualizer):
            raise ValueError("Invalid plotter type")
        self.visualizer = plotter

    def setup_animation_writer(
            self, 
            filename: str, 
            fps: int = 30, 
            codec: str = 'libx264', 
            bitrate: int = 1800,
            dpi: int = 100
    ):
        if self.visualizer is None:
            raise ValueError("No visualizer attached")
        self.visualizer.setup_animation_writer(
            filename = filename, 
            fps = fps, 
            codec =codec, 
            bitrate= bitrate,
            dpi= dpi
        )
        self.has_animation_writer = True

    def update_plot(self, title: int):
        if self.visualizer is None:
            return
        
        self.visualizer.update_particle_plots(
                title = title,
                particle_positions= self.particles.XP,
                particle_charges= self.particles.QP,
                particle_velocities= self.particles.UP,
                particle_deformations= self.particles.GP
            )
        self.visualizer.update_mesh_plots(
                title = title,
                pm_positions= self.particle_mesh.grid_positions,
                pm_velocities= self.particle_mesh.U,
                pm_charges= self.particle_mesh.RHS,
                pm_deformations= self.particle_mesh.deformation
            )
        if self.has_animation_writer:
            self.visualizer.grab_frame()

    ### VPM functions ###
    def vpm(
        self, 
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        mode: int,
        timestep: int,
        viscosity: float,
        num_particles: int | None = None,
    ):
        """
        Wrapper function for calling Fortran subroutine `vpm`.

        Args:
            num_equations (int): Number of equations to model
            particle_positions (np.ndarray): Particle positions array of shape (3, NVR_in)
            particle_charges (np.ndarray): Particle charges array of shape (num_equations + 1, NVR_in)
            mode (int): 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse, 6 - correct vorticity
            timestep (int): Timestep
            viscosity (float): Viscosity term for the diffusion equation
            num_particles (int, optional): Number of points to treat as particles. When None defaults
                        to all the points in the arrays. If different from None, the rest of the points
                        will be treated as source terms. Defaults to None.
        """
        # Check if the arrays are F_Arrays
        if isinstance(particle_positions, F_Array):
            particle_positions = particle_positions.data
        if isinstance(particle_charges, F_Array):
            particle_charges = particle_charges.data

        # Check that the arrays have the same number of particles
        NVR_size = particle_positions.shape[1]
        if not (particle_charges.shape[1] == NVR_size):
            raise ValueError("Number of particles in particle_charges does not match particle_positions")
        
        # Create F_Arrays to store the results
        particle_velocities = np.zeros_like(particle_positions, dtype=np.float64, order='K')  
        particle_deformations = np.zeros_like(particle_positions, dtype=np.float64, order='K')

        # Get the pointers to arrays for the particles
        XP_ptr = dp_array_to_pointer(particle_positions, copy = True)
        UP_ptr = dp_array_to_pointer(particle_velocities, copy = True)
        QP_ptr = dp_array_to_pointer(particle_charges, copy = True)
        GP_ptr = dp_array_to_pointer(particle_deformations, copy = True)

        # Get the pointers to arrays for the grid values
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size= 1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size= 1)

        if num_particles is None:
            num_particles = NVR_size

        self._lib.vpm(
            XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(num_equations)),byref(c_int(mode)), 
            byref(RHS_pm_ptr), byref(Velocity_ptr),
            byref(c_int(timestep)), byref(c_double(viscosity)),byref(c_int(NVR_size))
        )
        # store the results of the particles
        self.particles.particle_positions = pointer_to_dp_array(XP_ptr, particle_positions.shape)
        self.particles.particle_charges = pointer_to_dp_array(QP_ptr, particle_charges.shape)
        self.particles.particle_velocities = pointer_to_dp_array(UP_ptr, particle_velocities.shape)
        self.particles.particle_deformations = pointer_to_dp_array(GP_ptr, particle_deformations.shape)
                
        # store the results of the particle mesh
        neq = self.num_equations

        self.particle_mesh.number_equations = neq
        if not Velocity_ptr.is_null():
            U_arr = F_Array.from_ctype(Velocity_ptr, ownership=True, name = "U")
            self.particle_mesh.U = U_arr.transfer_data_ownership()
        
        if not RHS_pm_ptr.is_null():
            RHS_arr = F_Array.from_ctype(RHS_pm_ptr, ownership=True, name = "RHS")
            self.particle_mesh.RHS = RHS_arr.transfer_data_ownership()

    def vpm_project_solve(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        timestep: int,
        num_particles: int | None = None,
    ):
        """Project and solve VPM equations."""
        # Convert arrays
        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        
        # Validate
        NVR_size = self._validate_particle_counts(positions, charges)
        if num_particles is None:
            num_particles = NVR_size

        # Create pointers
        XP_ptr, QP_ptr = self._create_array_pointers(positions, charges)
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        # Call Fortran
        self._lib.vpm_project_solve(
            byref(c_int(timestep)),
            XP_ptr, QP_ptr,
            byref(c_int(num_particles)),
            byref(c_int(NVR_size)),
            byref(c_int(num_equations)),
            byref(RHS_pm_ptr),
        )

        # Store results
        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape)
        )
        self._store_mesh_results(RHS_pm_ptr)

    def vpm_define(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        timestep: int = 0,
        num_particles: int | None = None,
    ):
        """Define VPM problem."""
        if particle_positions is None or particle_charges is None:
            particle_positions = self.particles.particle_positions
            particle_charges = self.particles.particle_charges

        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        
        NVR_size = self._validate_particle_counts(positions, charges)
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr = self._create_array_pointers(positions, charges)

        self._lib.vpm_define(
            byref(c_int(timestep)), XP_ptr, QP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations))
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape)
        )

    def vpm_solve_velocity(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        timestep: int,
        num_particles: int | None = None,
    ):
        """Solve velocity field."""
        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        
        NVR_size = self._validate_particle_counts(positions, charges)
        velocities, deformations = self._initialize_particle_arrays(positions.shape)
        
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self._create_array_pointers(
            positions, charges, velocities, deformations
        )
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_solve_velocity(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr), byref(Velocity_ptr)
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape),
            pointer_to_dp_array(UP_ptr, velocities.shape),
            pointer_to_dp_array(GP_ptr, deformations.shape)
        )
        self._store_mesh_results(RHS_pm_ptr, Velocity_ptr)

    def vpm_solve_velocity_deformation(
        self,
        timestep: int = 0,
        num_equations: int | None = None,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        num_particles: int | None = None,
    ):
        """Solve velocity and deformation fields."""
        if particle_positions is None or particle_charges is None:
            particle_positions = self.particles.particle_positions
            particle_charges = self.particles.particle_charges
        
        if num_equations is None:
            num_equations = self.num_equations

        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        
        NVR_size = self._validate_particle_counts(positions, charges)
        velocities, deformations = self._initialize_particle_arrays(positions.shape)
        
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self._create_array_pointers(
            positions, charges, velocities, deformations
        )
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Deform_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_solve_velocity_deformation(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr),
            byref(Velocity_ptr), byref(Deform_ptr)
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape),
            pointer_to_dp_array(UP_ptr, velocities.shape),
            pointer_to_dp_array(GP_ptr, deformations.shape)
        )
        self._store_mesh_results(RHS_pm_ptr, Velocity_ptr, Deform_ptr)

    def vpm_interpolate(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        particle_velocities: np.ndarray | F_Array,
        particle_deformations: np.ndarray | F_Array,
        timestep: int,
        num_particles: int | None = None,
    ):
        """Interpolate particle fields."""
        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        velocities = self._convert_to_numpy(particle_velocities)
        deformations = self._convert_to_numpy(particle_deformations)
        
        NVR_size = self._validate_particle_counts(positions, charges, velocities, deformations)
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self._create_array_pointers(
            positions, charges, velocities, deformations
        )
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_interpolate(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr)
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape),
            pointer_to_dp_array(UP_ptr, velocities.shape),
            pointer_to_dp_array(GP_ptr, deformations.shape)
        )
        self._store_mesh_results(RHS_pm_ptr)

    def vpm_diffuse(
        self,
        num_equations: int,
        viscosity: float,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        num_particles: int | None = None,
    ):
        """Apply diffusion to particles."""
        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)

        (particle_velocities, particle_deformations) = self._initialize_particle_arrays(positions.shape)
        velocities = self._convert_to_numpy(particle_velocities)
        deformations = self._convert_to_numpy(particle_deformations)
        
        NVR_size = self._validate_particle_counts(positions, charges, velocities, deformations)
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self._create_array_pointers(
            positions, charges, velocities, deformations
        )
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_diffuse(
            byref(c_double(viscosity)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr)
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape),
            pointer_to_dp_array(UP_ptr, velocities.shape),
            pointer_to_dp_array(GP_ptr, deformations.shape)
        )
        self._store_mesh_results(RHS_pm_ptr)

    def vpm_correct_vorticity(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        num_particles: int | None = None,
    ):
        """Correct vorticity field."""
        positions = self._convert_to_numpy(particle_positions)
        charges = self._convert_to_numpy(particle_charges)
        
        NVR_size = self._validate_particle_counts(positions, charges)
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr  = self._create_array_pointers(
            positions, charges
        )

        self._lib.vpm_correct_vorticity(
            XP_ptr, QP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations))
        )

        self._store_particle_results(
            pointer_to_dp_array(XP_ptr, positions.shape),
            pointer_to_dp_array(QP_ptr, charges.shape),
        )

    def remesh_particles(self, project_particles: bool, particles_per_cell: int= 1, cut_off: float = 1e-9):
        """Remesh particles in 3D

        Args:
            project_particles (bool): Whether to project the particles or use RHS_pm 
            particles_per_cell (int): Number of particles per cell
            cut_off (float): Cut off value for the remeshing
        """
        NVR = self.particles.NVR
        neq = self.num_equations
        XP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        UP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        QP_struct = F_Array_Struct.null(ndims=2, total_size=(neq + 1)*NVR)
        GP_struct = F_Array_Struct.null(ndims=2, total_size=3*NVR)
        NVR = c_int(NVR)
        self._lib.remesh_particles_3d(
            byref(c_int(project_particles)), byref(c_int(particles_per_cell)),
            byref(XP_struct), byref(QP_struct),
            byref(UP_struct), byref(GP_struct),
            byref(NVR), byref(c_double(cut_off))
        )
        XP_arr = F_Array.from_ctype(XP_struct, ownership=True, name = "XP_remesh")
        QP_arr = F_Array.from_ctype(QP_struct, ownership=True, name = "QP_remesh")
        UP_arr = F_Array.from_ctype(UP_struct, ownership=True, name = "UP_remesh")
        GP_arr = F_Array.from_ctype(GP_struct, ownership=True, name = "GP_remesh")

        # store the results
        self.particles.particle_positions = XP_arr.transfer_data_ownership()
        self.particles.particle_velocities = UP_arr.transfer_data_ownership()
        self.particles.particle_charges = QP_arr.transfer_data_ownership()
        self.particles.particle_deformations = GP_arr.transfer_data_ownership()
        return XP_arr, QP_arr

    def set_verbosity(self, verbosity: int):
        """
            Set the verbosity of the VPM
        """
        self._lib.set_verbosity(byref(c_int(verbosity)))
    
    def set_case_folder(self, case_folder: str):
        """Set the case folder for VPM."""
        folder_b = case_folder.encode('utf-8')
        os.makedirs(case_folder, exist_ok=True)
        self._lib.set_case_folder(folder_b)
        
    def finalize(self):
        self._lib.finalize()
        if self.rank == 0:
            print_green(f"Finalized VPM. {self.rank}")

    def __del__(self):
        self.finalize()
        os.remove(self.vpm_lib._lib_vpm_path)

    ### Helper functions ###
    def _initialize_particle_arrays(self, shape: tuple) -> tuple[np.ndarray, np.ndarray]:
        """Initialize velocity and deformation arrays."""
        return (
            np.zeros(shape, dtype=np.float64, order='F'),
            np.zeros(shape, dtype=np.float64, order='F')
        )

    def _convert_to_numpy(self, array: np.ndarray | F_Array) -> np.ndarray:
        """Convert F_Array to numpy array if needed."""
        if isinstance(array, F_Array):
            return array.data
        
        if not isinstance(array, np.ndarray):
            raise ValueError("Invalid array type")
        
        if not array.flags['F_CONTIGUOUS']:
            return np.asfortranarray(array, dtype=np.float64)
        else:
            return array

    def _validate_particle_counts(self, positions: np.ndarray, *arrays: np.ndarray) -> int:
        """Validate particle counts match across arrays and return NVR_size."""
        NVR_size = positions.shape[1]
        if any(arr.shape[1] != NVR_size for arr in arrays):
            raise ValueError("Mismatched particle counts between arrays")
        return NVR_size

    def _create_array_pointers(self, *arrays: np.ndarray) -> list[_Pointer]:
        """Create pointers for arrays with copy."""
        return [dp_array_to_pointer(arr, copy=True) for arr in arrays]

    def _store_particle_results(self, positions, charges, velocities=None, deformations=None):
        """Store particle results back to class attributes."""
        self.particles.particle_positions = positions
        self.particles.particle_charges = charges
        if velocities is not None:
            self.particles.particle_velocities = velocities
        if deformations is not None:
            self.particles.particle_deformations = deformations

    def _store_mesh_results(self, RHS_ptr, Velocity_ptr=None, Deform_ptr=None):
        """Store mesh results if pointers are not null."""
        if not RHS_ptr.is_null():
            self.particle_mesh.RHS = F_Array.from_ctype(RHS_ptr, ownership=True, name="RHS").transfer_data_ownership()
        if Velocity_ptr and not Velocity_ptr.is_null():
            self.particle_mesh.U = F_Array.from_ctype(Velocity_ptr, ownership=True, name="U").transfer_data_ownership()
        if Deform_ptr and not Deform_ptr.is_null():
            self.particle_mesh.deformation = F_Array.from_ctype(Deform_ptr, ownership=True, name="Deform").transfer_data_ownership()

