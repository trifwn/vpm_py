module file_io
    use vpm_types, only: dp
    use console_io, only: vpm_print, nocolor, dummy_string
    use h5fortran

    implicit none
    integer, parameter :: MAX_STRING_LENGTH = 256

    ! FILES
    character(len=MAX_STRING_LENGTH), save :: case_folder = "vpm_case/"

    ! Particle output
    character(len=MAX_STRING_LENGTH), save :: particle_folder = "results/"
    character(len=MAX_STRING_LENGTH), save :: particle_output_file = "particles"

    ! Mesh output
    character(len=MAX_STRING_LENGTH), save :: mesh_folder = "results/"
    character(len=MAX_STRING_LENGTH), save :: mesh_output_file = "particle_mesh"
    character(len=MAX_STRING_LENGTH), save :: field_output_file = "field"

    ! Stat and Info output
    character(len=MAX_STRING_LENGTH), save :: yaps_time_file = "yaps_time.dat"
    character(len=MAX_STRING_LENGTH), save :: solve_stats_file = "solve_information.csv"
    public :: particle_output_file, mesh_output_file, case_folder

contains
    subroutine write_pm_solution_dat(NTIME, comp_grid, neqpm, RHS, SOL, velocity, deform_pm)
        use vpm_types, only: cartesian_grid

        implicit none
        type(cartesian_grid), intent(in) :: comp_grid
        integer, intent(in)              :: NTIME
        integer, intent(in)              :: neqpm
        real(dp), intent(in)             :: RHS(neqpm, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3)), &
                                            SOL(neqpm, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3))
        real(dp), intent(in)             :: velocity(3, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3))
        real(dp), intent(in), optional   :: deform_pm(3, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3))

        character*50                     :: filout
        integer                          :: i, j, k
        logical                          :: exist_flag
        real(dp)                         :: cellx, celly, cellz, velocx, velocy, velocz, &
                                            wmegax, wmegay, wmegaz, psi_1, psi_2, psi_3, &
                                            defx, defy, defz

        integer                          :: NXs, NYs, NZs, NXf, NYf, NZf
        real(dp)                         :: XMIN_pm, YMIN_pm, ZMIN_pm, XMAX_pm, YMAX_pm, ZMAX_pm
        real(dp)                         :: DXpm, DYpm, DZpm

        write (filout, '(a,i5.5,a)') trim(case_folder)//trim(mesh_folder), &
            NTIME, trim(mesh_output_file)//'.dat'

        write (dummy_string, "(A)") achar(9)//'Writing PM solution to file: '//trim(filout)
        call vpm_print(dummy_string, nocolor, 2)

        NXs = comp_grid%NN_bl(1)
        NYs = comp_grid%NN_bl(2)
        NZs = comp_grid%NN_bl(3)
        NXf = comp_grid%NN_bl(4)
        NYf = comp_grid%NN_bl(5)
        NZf = comp_grid%NN_bl(6)

        XMIN_pm = comp_grid%Xbound(1)
        YMIN_pm = comp_grid%Xbound(2)
        ZMIN_pm = comp_grid%Xbound(3)
        XMAX_pm = comp_grid%Xbound(4)
        YMAX_pm = comp_grid%Xbound(5)
        ZMAX_pm = comp_grid%Xbound(6)

        DXpm = comp_grid%Dpm(1)
        DYpm = comp_grid%Dpm(2)
        DZpm = comp_grid%Dpm(3)

        ! INQUIRE if file exists
        inquire (file=trim(filout), exist=exist_flag)
        ! if file exists open it with overwrite mode
        if (exist_flag) then
            open (1, file=trim(filout), status='replace')
        else
            open (1, file=trim(filout))
        end if

        if (present(deform_pm)) then
            write (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" '// &
                '"VORTX" "VORTY" "VORTZ" "PSI1" "PSI2" '// &
                '"PSI3" "DEFORMX" "DEFORMY" "DEFORMZ"'
        else
            write (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ" "PSI1" "PSI2" "PSI3" '
        end if

        write (1, "(a,I5.5,a,I5.5,a,I5.5)") 'ZONE I=', NXf - NXs + 1, &
            ' J=', NYf - NYs + 1, &
            ' K=', NZf - NZs + 1, ' F=POINT'

        do k = NZs, NZf
            do j = NYs, NYf
                do i = NXs, NXf
                    ! Structured grid coordinates
                    cellx = XMIN_pm + (I - 1)*DXpm
                    celly = YMIN_pm + (J - 1)*DYpm
                    cellz = ZMIN_pm + (K - 1)*DZpm

                    ! Calculated Velocity
                    velocx = velocity(1, i, j, k)
                    velocy = velocity(1, i, j, k)
                    velocz = velocity(1, i, j, k)

                    ! Deformation of vorticity
                    if (present(deform_pm)) then
                        defx = deform_pm(1, i, j, k)
                        defy = deform_pm(2, i, j, k)
                        defz = deform_pm(3, i, j, k)
                    else
                        defx = 0.d0
                        defy = 0.d0
                        defz = 0.d0
                    end if

                    ! RHS
                    wmegax = -RHS(1, I, J, K)
                    wmegay = -RHS(2, I, J, K)
                    wmegaz = -RHS(3, I, J, K)

                    ! Solution
                    psi_1 = SOL(1, I - NXs + 1, J - NYs + 1, K - NZs + 1)
                    psi_2 = SOL(2, I - NXs + 1, J - NYs + 1, K - NZs + 1)
                    psi_3 = SOL(3, I - NXs + 1, J - NYs + 1, K - NZs + 1)

                write (1, '(15(E14.7,1x))') cellx, celly, cellz,      &
                                            velocx, velocy, velocz,   &
                                            wmegax, wmegay, wmegaz,   &
                                            psi_1, psi_2, psi_3,      &
                        defx, defy, defz
                end do
            end do
        end do
        close (1)
    end subroutine write_pm_solution_dat

    subroutine write_particles_dat(NTIME, XPR, UPR, QPR, GPR, neq, NVR, NVR_size)
        implicit none
        integer, intent(in) :: NTIME, NVR, neq, NVR_size
        real(dp), intent(in):: XPR(3, NVR_size), QPR(neq + 1, NVR_size), UPR(3, NVR), GPR(3, NVR)
        integer ::i
        logical :: exist_flag
        character*80 :: filout1

        write (filout1, '(a,i5.5,a)') trim(case_folder)//trim(particle_folder), &
            NTIME, trim(particle_output_file)//'.dat'
        write (dummy_string, "(A)") achar(9)//'Writing particles to file: '//trim(filout1)
        call vpm_print(dummy_string, nocolor, 2)

        ! INQUIRE if file exists
        inquire (file=trim(filout1), exist=exist_flag)
        ! if file exists open it with overwrite mode
        if (exist_flag) then
            open (10, file=trim(filout1), status='replace')
        else
            open (10, file=trim(filout1))
        end if
        write (10, "(A)") 'VARIABLES = "i" "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
        write (10, "(a,I5.5,a)") 'ZONE I=', NVR, ' F=POINT'
        do i = 1, NVR
            write (10, '(i5.5,1x,12(E15.7,1x))') i, &
                XPR(1, i), XPR(2, i), XPR(3, i), &
                UPR(1, i), UPR(2, i), UPR(3, i), &
                QPR(:, i)
        end do

        ! call system('~/bin/preplot '&
        !    //filout1//' >/dev/null')
        ! call system('rm '//filout1)
        close (10)
    end subroutine write_particles_dat

    !-----------------------------------------------------------------
    !> \brief Writes particle data to an HDF5 file.
    !>
    !> This subroutine writes particle data including position, velocity,
    !> charge, and other properties to an HDF5 file. The data is written
    !> for a specific time step and can be optionally compressed.
    !>
    !> \param[in] NTIME        Integer representing the current time step.
    !> \param[in] XPR          Real array containing particle positions.
    !> \param[in] UPR          Real array containing particle velocities.
    !> \param[in] QPR          Real array containing particle charges.
    !> \param[in] GPR          Real array containing additional particle properties.
    !> \param[in] neq          Integer representing the number of equations.
    !> \param[in] NVR          Integer array containing the number of variables.
    !> \param[in] NVR_size     Integer representing the size of the NVR array.
    !> \param[in] compression  Integer representing the compression level.
    !-----------------------------------------------------------------
    subroutine write_particles_hdf5(NTIME, XPR, UPR, QPR, GPR, neq, NVR, NVR_size, compression)
        use h5fortran, only: hdf5_file
        implicit none

        integer, intent(in)              :: NTIME, NVR, neq, NVR_size
        real(dp), intent(in)             :: XPR(3, NVR_size), QPR(neq+1, NVR_size), &
                                            UPR(3, NVR),      GPR(3, NVR)
        integer, optional                :: compression
        integer                          :: comp_level = 4
        type(hdf5_file)                  :: h5f
        character(len=MAX_STRING_LENGTH) :: filout1

        ! Create the filename for the output
        write (filout1, '(A,I5.5,A)') trim(case_folder)//trim(particle_folder), &
            NTIME, trim(particle_output_file)//".h5"
        write (dummy_string, "(A)") achar(9)//'Writing particles to HDF5 file: '//trim(filout1)
        call vpm_print(dummy_string, nocolor, 2)
        ! Open the HDF5 file for writing
        if (present(compression)) then
            comp_level = compression
        end if
        call h5f%open(trim(filout1), action='w', comp_lvl=comp_level)

        ! Write attributes to the HDF5 file
        call h5f%writeattr('/', 'NTIME', NTIME)
        call h5f%writeattr('/', 'NVR', NVR)
        call h5f%writeattr('/', 'neq', neq)
        call h5f%writeattr('/', 'NVR_size', NVR_size)

        ! Write XPR (Particle positions) to the HDF5 file
        call h5f%write('/XPR', XPR)

        ! Write UPR (Particle velocities)
        call h5f%write('/UPR', UPR)

        ! Write QPR (Vorticity or other quantities)
        call h5f%write('/QPR', QPR)

        ! Write GPR (Other relevant particle data)
        call h5f%write('/GPR', GPR)

        ! Close the HDF5 file to ensure all data is written to disk
        call h5f%close()

    end subroutine write_particles_hdf5

    !----------------------------------------------------------------------
    !> \brief Writes the PM solution to an HDF5 file.
    !>
    !> This subroutine handles the output of the PM solution data into an
    !> HDF5 file format. It includes various parameters such as time step,
    !> node information, right-hand side values, solution vectors, velocity
    !> components, deformation components, and compression data.
    !>
    !> \param[in] NTIME        Time step.
    !> \param[in] NN_in        Number of nodes.
    !> \param[in] NNbl_in      Number of nodes excluding dumb cells 
    !> \param[in] neqpm        Number of equations in the PM solver.
    !> \param[in] RHS          Right-hand side of PM equation (Charges).
    !> \param[in] SOL          Solution vector.
    !> \param[in] velocity     Velocity components.
    !> \param[in] deformation  Deformation components.
    !> \param[in] compression  Compression level.
    !----------------------------------------------------------------------
    subroutine write_pm_solution_hdf5(                                  &
        NTIME, NN_in, NNbl_in, neqpm, RHS, SOL, velocity,               &
        deformation, compression                 &
    )
        use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm
        use vpm_vars, only: timestep_info
        use h5fortran

        integer, intent(in)              :: NTIME
        integer, intent(in)              :: NN_in(3), NNbl_in(6)
        integer, intent(in)              :: neqpm
        real(dp), intent(in)             :: RHS(neqpm, NN_in(1), NN_in(2), NN_in(3)), &
                                            SOL(neqpm, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in)             :: velocity(3, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in), optional   :: deformation(3, NN_in(1), NN_in(2), NN_in(3))
        integer, optional                :: compression 
        integer                          :: comp_level = 3
        
        integer                          :: NXs, NYs, NZs, NXf, NYf, NZf
        character(len=MAX_STRING_LENGTH) :: filout
        integer                          :: i, j, k
        real(dp), dimension(NN_in(1), NN_in(2), NN_in(3)) :: X, Y, Z

        type(hdf5_file)                  :: h5f

        ! Construct file name
        write(filout, '(a,i5.5,a)') trim(case_folder)//trim(mesh_folder), &
                NTIME, trim(mesh_output_file) // ".h5"
        ! Informative print statement
        write(dummy_string, "(A)") achar(9)//'Writing PM solution to HDF5 file: '//trim(filout)
        call vpm_print(dummy_string, nocolor, 2)

        NXs = NNbl_in(1)
        NYs = NNbl_in(2)
        NZs = NNbl_in(3)
        NXf = NNbl_in(4)
        NYf = NNbl_in(5)
        NZf = NNbl_in(6)

        if (present(compression)) then
            comp_level = compression
        end if
        ! Open HDF5 file
        call h5f%open(trim(filout), action='w', comp_lvl= comp_level)
        ! Write Attributes
        call h5f%writeattr('/','NTIME', NTIME)
        call h5f%writeattr('/','neqpm', neqpm)
        call h5f%writeattr('/','NN_in', NN_in)
        call h5f%writeattr('/','NNbl_in', NNbl_in)

        ! Create datasets for each variable
        call h5f%create('/X', H5T_NATIVE_DOUBLE, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
        call h5f%create('/Y', H5T_NATIVE_DOUBLE, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
        call h5f%create('/Z', H5T_NATIVE_DOUBLE, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])

        ! call h5f%create('/VEL', H5T_NATIVE_DOUBLE, dset_dims=[3, NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])

        ! Iterate over the 3D grid and write data to HDF5 file
        do k = NZs, NZf
            do j = NYs, NYf
                do i = NXs, NXf
                ! Structured grid coordinates
                X(i - NXs + 1, j - NYs + 1, k - NZs + 1) = XMIN_pm + (i - 1) * DXpm
                Y(i - NXs + 1, j - NYs + 1, k - NZs + 1) = YMIN_pm + (j - 1) * DYpm
                Z(i - NXs + 1, j - NYs + 1, k - NZs + 1) = ZMIN_pm + (k - 1) * DZpm
                end do
            end do
        end do

                ! Write the structured grid coordinates to the HDF5 file
        call h5f%write('/X', X)
        call h5f%write('/Y', Y)
        call h5f%write('/Z', Z)
        
        ! Write the vorticity field to the HDF5 file
        call h5f%write('/RHS', RHS(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf))

        ! Write the solution field to the HDF5 file
        call h5f%write('/SOL', SOL(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf))
        
        ! Write the velocity field to the HDF5 file
        call h5f%write('/VEL', velocity(:, NXs:NXf, NYs:NYf, NZs:NZf))
        
        ! Write the deformation field to the HDF5 file
        if (present(deformation)) then
            call h5f%create('/VORTEXSTRETCH', H5T_NATIVE_DOUBLE, dset_dims=[3, NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
            call h5f%write('/VORTEXSTRETCH', deformation(1:3, NXs:NXf, NYs:NYf, NZs:NZf)) 
        end if
        
        ! Close HDF5 file
        call h5f%close()

        ! Add timestep metadata to the HDF5 file
        call write_timestep_metadata_hdf5(filout, timestep_info) 
        
    end subroutine write_pm_solution_hdf5

    subroutine write_pressure_hdf5(                                     &
        NTIME, NN_in, NNbl_in, pressure, compression                    &
    )
        use h5fortran
        integer, intent(in)                 :: NTIME
        integer, intent(in)                 :: NN_in(3), NNbl_in(6)
        real(dp), intent(in)                :: pressure(2, NN_in(1), NN_in(2), NN_in(3))
        integer, optional                   :: compression 
        integer                             :: comp_level = 3
        integer                             :: NXs, NYs, NZs, NXf, NYf, NZf
        character(len=MAX_STRING_LENGTH)    :: filout
        type(hdf5_file)                     :: h5f

        write(filout, '(a,i5.5,a)') trim(case_folder)//trim(mesh_folder), &
                NTIME, trim(mesh_output_file) // ".h5"
        
        write(dummy_string, "(A)") achar(9)//'Writing Pressure to HDF5 file: '//trim(filout)
        call vpm_print(dummy_string, nocolor, 2)

        NXs = NNbl_in(1)
        NYs = NNbl_in(2)
        NZs = NNbl_in(3)
        NXf = NNbl_in(4)
        NYf = NNbl_in(5)
        NZf = NNbl_in(6)

        if (present(compression)) then
            comp_level = compression
        end if
        ! Open HDF5 file
        call h5f%open(trim(filout), action='rw', comp_lvl= comp_level)  
        ! Write the velocity field to the HDF5 file
        call h5f%write('/Q', pressure(1, NXs:NXf, NYs:NYf, NZs:NZf))
        call h5f%write('/P', pressure(2, NXs:NXf, NYs:NYf, NZs:NZf))
        ! Close HDF5 file
        call h5f%close()
    end subroutine write_pressure_hdf5

    subroutine write_timestep_metadata_hdf5(filename , timestep)
        use vpm_types, only: timestepInformation
        implicit none
        type(timestepInformation), intent(in) :: timestep
        character(len=256), intent(in)        :: filename
        type(hdf5_file)    :: h5f

        print *, "Writing timestep metadata to HDF5 file: ", trim(filename)
        ! Determine the file status and write mode
        call h5f%open(trim(filename), action='rw')

        ! Write Attributes
        ! call h5f%writeattr('/', 'NTIME', timestep%n)
        ! call h5f%writeattr('/', 'dt', timestep_info%dt)
        ! call h5f%writeattr('/', 't', timestep_info%t)
        call h5f%writeattr('/', 'NVR', timestep%NVR)
        call h5f%writeattr('/', 'NVR_size', timestep%NVR_size)
        call h5f%writeattr('/', 'NN', timestep%NN)
        call h5f%writeattr('/', 'Xbound', timestep%Xbound)
        call h5f%writeattr('/', 'Dpm', timestep%Dpm)
        call h5f%writeattr('/', 'solver', timestep%solver)

        call h5f%writeattr('/', 'mean_div_w', timestep%mean_div_w)
        call h5f%writeattr('/', 'max_div_w', timestep%max_div_w)
        call h5f%writeattr('/', 'min_div_w', timestep%min_div_w)
        call h5f%writeattr('/', 'mean_div_u', timestep%mean_div_u)
        call h5f%writeattr('/', 'max_div_u', timestep%max_div_u)
        call h5f%writeattr('/', 'min_div_u', timestep%min_div_u)

        ! Particle Mesh Attributes
        call h5f%writeattr('/', 'total_pm_momentum_x', timestep%total_momentum_x_pm)
        call h5f%writeattr('/', 'total_pm_momentum_y', timestep%total_momentum_y_pm)
        call h5f%writeattr('/', 'total_pm_momentum_z', timestep%total_momentum_z_pm)
        call h5f%writeattr('/', 'total_pm_kinetic_energy', timestep%total_kinetic_energy_pm)
        call h5f%writeattr('/', 'total_pm_vorticity', timestep%total_vorticity_pm)
        call h5f%writeattr('/', 'total_pm_enstrophy', timestep%total_enstrophy_pm)

        ! Particle Attributes
        call h5f%writeattr('/', 'total_particle_momentum_x', timestep%total_momentum_x_particles)
        call h5f%writeattr('/', 'total_particle_momentum_y', timestep%total_momentum_y_particles)
        call h5f%writeattr('/', 'total_particle_momentum_z', timestep%total_momentum_z_particles)
        call h5f%writeattr('/', 'total_particle_kinetic_energy', timestep%total_kinetic_energy_particles)
        call h5f%writeattr('/', 'total_particle_vorticity', timestep%total_vorticity_particles)
        call h5f%writeattr('/', 'total_particle_enstrophy', timestep%total_enstrophy_particles)

        ! Close HDF5 file
        call h5f%close()
        print *, "Timestep metadata written to HDF5 file: ", trim(filename)
    end subroutine write_timestep_metadata_hdf5

    subroutine write_timestep_information_dat(timestep_info, NTIME_pm)
        use vpm_types, only: timestepInformation
        implicit none
        type(timestepInformation), intent(in) :: timestep_info
        integer                               :: NTIME_pm
        character(len=256) :: filename
        logical :: file_exists
        real(dp)           :: mean_div_w, max_div_w, min_div_w
        real(dp)           :: mean_div_u, max_div_u, min_div_u
        real(dp)           :: total_kinetic_energy, total_vorticity, total_enstrophy
        real(dp)           :: total_momentum_x, total_momentum_y, total_momentum_z

        mean_div_w = timestep_info%mean_div_w
        max_div_w  = timestep_info%max_div_w
        min_div_w  = timestep_info%min_div_w
        
        mean_div_u = timestep_info%mean_div_u
        max_div_u  = timestep_info%max_div_u
        min_div_u  = timestep_info%min_div_u
        
        total_momentum_x = timestep_info%total_momentum_x_pm 
        total_momentum_y = timestep_info%total_momentum_y_pm 
        total_momentum_z = timestep_info%total_momentum_z_pm 
        total_kinetic_energy = timestep_info%total_kinetic_energy_pm
        total_vorticity = timestep_info%total_vorticity_pm
        total_enstrophy = timestep_info%total_enstrophy_pm

        ! Construct the filename
        write (filename, "(A, A)") trim(case_folder), trim(solve_stats_file)
        ! Determine the file status and write mode
        inquire(file=filename, exist=file_exists)

        if (NTIME_pm == 1 .or. NTIME_pm == 0 .or. .not. file_exists) then
            open (unit=10, file=filename, status='replace', action='write')
            ! Write CSV header
            write (10, '(A)') "Iteration, Total_Enstrophy, Total_Vorticity,"              &
                            //"Total_Kinetic_Energy, Momentum_X, Momentum_Y, Momentum_Z," &
                            //"MEAN_DIV_W,MAX_DIV_W, MEAN_DIV_VEL, MAX_DIV_VEL"
        else
            open (unit=10, file=filename, status='old', position='append', action='write')
        end if

        write (10, '(I0,16(",",ES14.7))') NTIME_pm,                                 &
            total_enstrophy, total_vorticity, total_kinetic_energy,                 &                  
            total_momentum_x, total_momentum_y, total_momentum_z,                   &
            mean_div_w, max_div_w, mean_div_u, max_div_u
        close (10)
    end subroutine write_timestep_information_dat

    subroutine write_field_hdf5(neq, comp_grid, field, field_name, compression)
        use h5fortran
        use vpm_types, only: cartesian_grid
        implicit none

        integer, intent(in)              :: neq
        type(cartesian_grid), intent(in) :: comp_grid
        real(dp), intent(in)             :: field(neq, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3))
        character(len=MAX_STRING_LENGTH) :: field_name
        integer, optional                :: compression
        integer                          :: comp_level = 0
        character(len=MAX_STRING_LENGTH) :: filout
        real(dp), dimension(3)           :: origin, spacing
        type(hdf5_file)                  :: h5f

        origin = [0.0_dp, 0.0_dp, 0.0_dp]
        spacing = [comp_grid%Dpm(1), comp_grid%Dpm(2), comp_grid%Dpm(3)]
        ! Construct file name
        write (filout, '(a)') trim(case_folder)//trim(mesh_folder)//trim(field_name)//".h5"
        ! Informative print statement
        write (dummy_string, "(A)") achar(9)//'Writing field to HDF5 file: '//trim(filout)
        call vpm_print(dummy_string, nocolor, 1)

        if (present(compression)) then
            comp_level = compression
        end if
        ! Open HDF5 file
        call h5f%open(trim(filout), action='w', comp_lvl=comp_level)
        ! Write Attributes
        call h5f%writeattr('/', 'dims', comp_grid%NN)
        call h5f%writeattr('/', 'origin', origin)
        call h5f%writeattr('/', 'spacing', spacing)
        call h5f%writeattr('/', 'NN_bl', comp_grid%NN_bl)

        ! Create datasets for each variable
        call h5f%create('/FIELD', H5T_NATIVE_DOUBLE, dset_dims=[neq, comp_grid%NN(1), comp_grid%NN(2), comp_grid%NN(3)])

        ! Write the field to the HDF5 file
        call h5f%write('/FIELD', field)
        call h5f%writeattr('/FIELD', 'AttributeType', 'Vector')
        call h5f%writeattr('/FIELD', 'Center', 'Node')

        ! Close HDF5 file
        call h5f%close()

    end subroutine write_field_hdf5

end module file_io