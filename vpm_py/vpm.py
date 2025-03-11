import os
import numpy as np
from ctypes import (byref,  c_double, c_int)
import gc

# Local imports
from . import ParticleMesh, Particles
from .arrays import F_Array, F_Array_Struct
from .console_io import print_green, print_IMPORTANT
from .utils import divide_processors
from .vpm_dtypes import pointer_to_dp_array
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
        self.verbosity = verbocity

        if case_folder is not None:
            self.set_case_folder(case_folder)
        # Initialize the VPM
        self.initialize(self.dpm[0], self.dpm[1], self.dpm[2], NBI, NBJ, NBK)
        self.particle_mesh = ParticleMesh()
        self.particles = Particles(number_equations= self.num_equations)

        # Visualization
        self.visualizer: Visualizer | None = None
        self.has_animation_writer = False

    @property
    def verbosity(self):
        return self._verbosity

    @verbosity.setter
    def verbosity(self, value: int):
        self._verbosity = value
        self.set_verbosity(value)

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
        constant_particle_volume: bool = True, 
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
            constant_particle_volume (bool, optional): Constant particle volume. Defaults to True.
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
            byref(c_int(boundary_condition)), byref(c_int(constant_particle_volume)), byref(c_int(ncoarse)),
            byref(c_int(NBI)), byref(c_int(NBJ)), byref(c_int(NBK)), byref(c_int(remesh)), 
            byref(c_int(use_tree)), byref(c_int(ilevmax)), byref(c_int(OMPTHREADS)), byref(c_int(is_domain_fixed)),
            byref(c_int(IPMWRITE)), byref(c_int(IPMWSTART)), byref(c_int(IPMWSTEPS)),
        )
        if(self.rank == 0) and (self.verbosity > 1):
            print_green(f"Finished initializing VPM {self.rank}:")
            # Print the arguments passed
            print(f"\tDXpm= {DXpm}")
            print(f"\tDYpm= {DYpm}")
            print(f"\tDZpm= {DZpm}")
            print(f"\tinterf_iproj= {projection_type}")
            print(f"\tibctyp= {boundary_condition}")
            print(f"\tIDVPM= {constant_particle_volume}")
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
        # Check if filename exists
        if os.path.exists(filename):
            # Then append a number to the filename
            i = 1
            while os.path.exists(f"{filename}_{i}.mp4"):
                i += 1
            filename = f"{filename}_{i}.mp4"
            
        self.visualizer.setup_animation_writer(
            filename = filename, 
            fps = fps, 
            codec =codec, 
            bitrate= bitrate,
            dpi= dpi
        )
        self.has_animation_writer = True

    def update_plot(self, title: str, dt = None):
        if self.visualizer is None:
            return
        
        self.visualizer.update_all_plots(
                title = title,
                particle_positions= self.particles.particle_positions,
                particle_charges= self.particles.particle_charges,
                particle_velocities= self.particles.particle_velocities,
                particle_deformations= self.particles.particle_deformations,
                pm_positions= self.particle_mesh.grid_positions,
                pm_velocities= self.particle_mesh.velocity,
                pm_charges= self.particle_mesh.RHS,
                pm_vortex_stretching= self.particle_mesh.deformation,
                pm_pressure= self.particle_mesh.pressure,
                pm_q_pressure= self.particle_mesh.q_pressure,
                pm_u_pressure= self.particle_mesh.u_pressure,
            )
        self.visualizer.set_problem_info(
            num_particles = self.particles.NVR,
            grid_size = self.particle_mesh.grid_size,
            dpm= self.particle_mesh.grid_spacing,
            dt = dt,
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
        # Get the pointers to arrays for the particles
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        XP_ptr, UP_ptr, QP_ptr, GP_ptr = self.particles.get_array_pointers(return_all=True) 
        NVR_size = self.particles.NVR
        if num_particles is None:
            num_particles = NVR_size

        # Get the pointers to arrays for the grid values
        # Create null pointers for the arrays using ctypes
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm(
            XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(num_equations)),byref(c_int(mode)), 
            byref(RHS_pm_ptr), byref(Velocity_ptr),
            byref(c_int(timestep)), byref(c_double(viscosity)),byref(c_int(NVR_size))
        )

        # store the results of the particles
        NVR = self.particles.NVR
        self.particles.store_particles(
            positions= pointer_to_dp_array(XP_ptr, (3, NVR)), 
            charges = pointer_to_dp_array(QP_ptr, (self.num_equations + 1, NVR)), 
            velocities = pointer_to_dp_array(UP_ptr, (3, NVR)), 
            deformations = pointer_to_dp_array(GP_ptr, (3, NVR)), 
        )
        self.particle_mesh.store_mesh_results(
            rhs         = self.dereference_F_Array(RHS_pm_ptr),
            velocity    = self.dereference_F_Array(Velocity_ptr),
        )
        # Call garbage collector
        gc.collect()

    def vpm_project_solve(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        timestep: int,
        num_particles: int | None = None,
    ):
        """Project and solve VPM equations."""
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        # Create pointers
        XP_ptr, QP_ptr = self.particles.get_array_pointers()
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
        self.particle_mesh.store_mesh_results(rhs = self.dereference_F_Array(RHS_pm_ptr))
        # Call garbage collector
        gc.collect()

    def vpm_define(
        self,
        num_equations: int | None = None,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        timestep: int = 0,
        num_particles: int | None = None,
    ):
        """Define VPM problem."""
        if num_equations is None:
            num_equations = self.num_equations

        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.particle_positions.shape[1] 
        self.particles.validate_stored_particle_counts(NVR_size) 
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr = self.particles.get_array_pointers()

        self._lib.vpm_define(
            byref(c_int(timestep)), XP_ptr, QP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations))
        )

        NVR = self.particles.NVR
        self.particles.store_particles(
            positions= pointer_to_dp_array(XP_ptr, (3, NVR)), 
            charges= pointer_to_dp_array(QP_ptr, (self.num_equations + 1, NVR)), 
        )
        # Call garbage collector
        gc.collect()

    def vpm_solve_velocity(
        self,
        num_equations: int,
        particle_positions: np.ndarray | F_Array,
        particle_charges: np.ndarray | F_Array,
        timestep: int,
        num_particles: int | None = None,
    ):
        """Solve velocity field."""
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self.particles.get_array_pointers(return_all = True)
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_solve_velocity(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr), byref(Velocity_ptr)
        )

        NVR = self.particles.NVR
        self.particles.store_particles(
            positions = pointer_to_dp_array(XP_ptr, (3, NVR)), 
            charges = pointer_to_dp_array(QP_ptr, (self.num_equations + 1, NVR)), 
            velocities = pointer_to_dp_array(UP_ptr, (3, NVR)), 
        )
        self.particle_mesh.store_mesh_results(
            rhs         = self.dereference_F_Array(RHS_pm_ptr),
            velocity    = self.dereference_F_Array(Velocity_ptr),
        )
        # Call garbage collector
        gc.collect()

    def vpm_solve_velocity_deformation(
        self,
        timestep: int = 0,
        num_equations: int | None = None,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        num_particles: int | None = None,
    ):
        """Solve velocity and deformation fields."""
        if num_equations is None:
            num_equations = self.num_equations
        
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self.particles.get_array_pointers(return_all = True)
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Velocity_ptr = F_Array_Struct.null(ndims=4, total_size=1)
        Deform_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_solve_velocity_deformation(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)), byref(RHS_pm_ptr),
            byref(Velocity_ptr), byref(Deform_ptr)
        )

        NVR = self.particles.NVR
        self.particles.store_particles(
            positions = pointer_to_dp_array(XP_ptr, (3, NVR)), 
            charges = pointer_to_dp_array(QP_ptr, (self.num_equations + 1, NVR)), 
            velocities = pointer_to_dp_array(UP_ptr, (3, NVR)), 
            deformations = pointer_to_dp_array(GP_ptr, (3, NVR)), 
        )
        # self._store_mesh_results(RHS_pm_ptr, Velocity_ptr, Deform_ptr)
        self.particle_mesh.store_mesh_results(
            rhs         = self.dereference_F_Array(RHS_pm_ptr),
            velocity    = self.dereference_F_Array(Velocity_ptr),
            deformation = self.dereference_F_Array(Deform_ptr),
        )
        # Call garbage collector
        gc.collect()

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
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges

        if particle_velocities is not None:
            self.particles.particle_velocities = particle_velocities
        
        if particle_deformations is not None:
            self.particles.particle_deformations = particle_deformations
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self.particles.get_array_pointers(return_all = True)
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_interpolate(
            byref(c_int(timestep)), XP_ptr, QP_ptr, UP_ptr, GP_ptr,
            byref(c_int(num_particles)), byref(c_int(NVR_size)),
            byref(c_int(num_equations)),
            byref(RHS_pm_ptr)
        )

        NVR = self.particles.NVR
        self.particles.store_particles(
            positions = pointer_to_dp_array(XP_ptr, (3, NVR)), 
            charges = pointer_to_dp_array(QP_ptr, (self.num_equations + 1, NVR)), 
            velocities = pointer_to_dp_array(UP_ptr, (3, NVR)), 
            deformations = pointer_to_dp_array(GP_ptr, (3, NVR)), 
        )
        self.particle_mesh.store_mesh_results(rhs = self.dereference_F_Array(RHS_pm_ptr))
        # Call garbage collector
        gc.collect()

    def vpm_diffuse(
        self,
        viscosity: float,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        num_particles: int | None = None,
    ):
        """Apply diffusion to particles."""
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions

        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR 
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        XP_ptr, QP_ptr, UP_ptr, GP_ptr = self.particles.get_array_pointers(return_all = True)
        RHS_pm_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_diffuse(
            byref(c_double(viscosity)), XP_ptr, QP_ptr, UP_ptr, GP_ptr, 
            byref(c_int(num_particles)),
            byref(c_int(NVR_size)),
            byref(c_int(self.num_equations)),
            byref(RHS_pm_ptr)
        )

        NVR = self.particles.NVR
        self.particles.store_particles(
            deformations = pointer_to_dp_array(GP_ptr, (3, NVR)), 
        )
        self.particle_mesh.store_mesh_results(rhs = self.dereference_F_Array(RHS_pm_ptr))
        # Call garbage collector
        gc.collect()

    def vpm_correct_vorticity(
        self,
        particle_positions: np.ndarray | F_Array | None = None,
        particle_charges: np.ndarray | F_Array | None = None,
        num_particles: int | None = None,
    ):
        """Correct vorticity field."""
        if particle_positions is not None:
            self.particles.particle_positions = particle_positions
        if particle_charges is not None:
            self.particles.particle_charges = particle_charges
        
        NVR_size = self.particles.NVR
        if num_particles is None:
            num_particles = NVR_size
        self.particles.validate_stored_particle_counts(NVR_size)

        XP_ptr, QP_ptr = self.particles.get_array_pointers()
        num_equations = self.num_equations

        NVR_size = c_int(NVR_size)
        self._lib.vpm_correct_vorticity(
            XP_ptr, QP_ptr,
            byref(c_int(num_particles)),
            byref(c_int(num_equations)),
            byref(NVR_size),
        )
        NVR = self.particles.NVR
        self.particles.store_particles(
            charges = pointer_to_dp_array(QP_ptr, (num_equations + 1, NVR))
        )
        # Call garbage collector
        gc.collect()
    
    def vpm_solve_pressure( 
        self,
        velocity: F_Array | np.ndarray | None = None,
        vorticity: F_Array | np.ndarray | None = None,
        density: float = 1.225,
        viscosity: int = 1,
    ):
        """Solve pressure field."""
        if velocity is None:
            if self.particle_mesh.velocity is not None:
                velocity = self.particle_mesh.velocity[:, :, :, :]
            else:
                velocity = self.particle_mesh.velocity
        if vorticity is None:
            if self.particle_mesh.RHS is not None:
                vorticity = self.particle_mesh.RHS[:3, :, :, :]
            else:
                vorticity = self.particle_mesh.RHS

        # Assert that the arrays are of the same shape
        if velocity.shape != vorticity.shape:
            print_IMPORTANT(f"Velocity shape: {velocity.shape}, Vorticity shape: {vorticity.shape}")
            print_IMPORTANT("Velocity and vorticity arrays must be of the same shape")
            return
        
        # If numpy arrays are passed, convert them to F_Array
        if isinstance(velocity, np.ndarray) or velocity is None:
            velocity = F_Array.from_ndarray(velocity)
        if isinstance(vorticity, np.ndarray) or vorticity is None:
            vorticity = F_Array.from_ndarray(vorticity)

        if not isinstance(velocity, F_Array) or not isinstance(vorticity, F_Array):
            raise ValueError("Invalid array type")

        Velocity_ptr  = velocity.to_ctype()
        Vorticity_ptr = vorticity.to_ctype()
        Pressures_ptr = F_Array_Struct.null(ndims=4, total_size=1)

        self._lib.vpm_solve_pressure(
            Vorticity_ptr, Velocity_ptr, Pressures_ptr, byref(c_double(density)), byref(c_double(viscosity))
        )

        self.particle_mesh.store_mesh_results(pressures= self.dereference_F_Array(Pressures_ptr))
        # Call garbage collector
        gc.collect()


    def remesh_particles(
        self, 
        project_particles: bool, 
        particles_per_cell: int= 1, 
        cut_off: float = 1e-9
    ):
        """Remesh particles in 3D

        Args:
            project_particles (bool): Whether to project the particles or use RHS_pm 
            particles_per_cell (int): Number of particles per cell
            cut_off (float): Cut off value for the remeshing
        """
        NVR = self.particles.NVR

        self.particles.validate_stored_particle_counts()
        XP_struct = F_Array.from_ndarray(self.particles.particle_positions).to_ctype()
        UP_struct = F_Array.from_ndarray(self.particles.particle_velocities).to_ctype()
        QP_struct = F_Array.from_ndarray(self.particles.particle_charges).to_ctype()
        GP_struct = F_Array.from_ndarray(self.particles.particle_deformations).to_ctype()
        NVR = c_int(NVR)
        self._lib.remesh_particles_3d(
            byref(c_int(project_particles)), byref(c_int(particles_per_cell)),
            byref(XP_struct), byref(QP_struct),
            byref(UP_struct), byref(GP_struct),
            byref(NVR), 
            byref(c_double(cut_off))
        )
        NVR = NVR.value

        assert(NVR == self.particles.NVR)
        # store the results
        XP_arr = F_Array.from_ctype(XP_struct)
        QP_arr = F_Array.from_ctype(QP_struct)
        UP_arr = F_Array.from_ctype(UP_struct)
        GP_arr = F_Array.from_ctype(GP_struct)

        self.particles.store_particles(
            positions = XP_arr.to_numpy(copy= True) if XP_arr.total_size > 0 else None, 
            charges = QP_arr.to_numpy(copy= True) if QP_arr.total_size > 0 else None,
            velocities =UP_arr.to_numpy(copy= True) if UP_arr.total_size > 0 else None,
            deformations = GP_arr.to_numpy(copy= True) if GP_arr.total_size > 0 else None,
        )
        # Call garbage collector
        gc.collect()
        return self.particles.particle_positions, self.particles.particle_charges 

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
    def dereference_F_Array(self, F_Array_ptr: F_Array_Struct) -> F_Array | None:
        """Dereference F_Array_Struct to F_Array."""
        if not F_Array_ptr.is_null() and F_Array_ptr.total_size > 0:
            return F_Array.from_ctype(F_Array_ptr)
        return None
