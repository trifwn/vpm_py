module file_io
    use base_types, only: dp
    use console_io, only: vpm_print, nocolor, dummy_string

    implicit none
    integer, parameter :: MAX_STRING_LENGTH= 256

    ! FILES
    character (len=MAX_STRING_LENGTH), save :: case_folder = "vpm_case/"

    ! Particle output
    character (len=MAX_STRING_LENGTH), save :: particle_folder = "results/"
    character (len=MAX_STRING_LENGTH), save :: particle_output_file = "particles"

    ! Mesh output
    character (len=MAX_STRING_LENGTH), save :: mesh_folder = "results/"
    character (len=MAX_STRING_LENGTH), save :: mesh_output_file = "particle_mesh"

    ! Stat and Info output
    character (len=MAX_STRING_LENGTH), save :: yaps_time_file = "yaps_time.dat"
    character (len=MAX_STRING_LENGTH), save :: solve_stats_file = "solve_information.csv"
    public :: particle_output_file, mesh_output_file, case_folder

contains

    subroutine write_pm_solution(NTIME, NN_in, NNbl_in, neqpm, RHS, SOL, velocity, deform_pm)
        use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm
        
        integer, intent(in)              :: NTIME
        integer, intent(in)              :: NN_in(3), NNbl_in(6)
        integer, intent(in)              :: neqpm
        real(dp), intent(in)             :: RHS(neqpm, NN_in(1), NN_in(2), NN_in(3)), &
                                            SOL(neqpm, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in)             :: velocity(3, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in), optional   :: deform_pm(3, NN_in(1), NN_in(2), NN_in(3))                                            
        integer                          :: NXs, NYs, NZs, NXf, NYf, NZf

        character*50                     :: filout
        integer                          :: i, j, k, size_eq
        logical                          :: exist_flag
        real(dp)                         :: cellx, celly, cellz, velocx, velocy, velocz, &
                                            wmegax, wmegay, wmegaz, psi_1, psi_2, psi_3, &
                                            defx, defy, defz

        
        write (filout, '(a,i5.5,a)') trim(case_folder)//trim(mesh_folder), &
                NTIME, trim(mesh_output_file) // '.dat'

        write (dummy_string, "(A)") achar(9)//'Writing PM solution to file: '//trim(filout)
        call vpm_print(dummy_string, nocolor, 2)

        NXs = NNbl_in(1)
        NYs = NNbl_in(2)
        NZs = NNbl_in(3)
        NXf = NNbl_in(4)
        NYf = NNbl_in(5)
        NZf = NNbl_in(6)

        ! INQUIRE if file exists
        inquire (file=trim(filout), exist=exist_flag)
        ! if file exists open it with overwrite mode
        if (exist_flag) then
            open (1, file=trim(filout), status='replace')
        else
            open (1, file=trim(filout))
        end if

        if (present(deform_pm)) then
            write (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" '//&
                            '"VORTX" "VORTY" "VORTZ" "PSI1" "PSI2" '//&
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
                cellx =  XMIN_pm + (I - 1)*DXpm
                celly =  YMIN_pm + (J - 1)*DYpm
                cellz =  ZMIN_pm + (K - 1)*DZpm

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
    end subroutine write_pm_solution

    subroutine write_particles(NTIME, XPR, UPR, QPR, GPR, neq, NVR, NVR_size)
        Implicit None
        integer, intent(in) :: NTIME, NVR, neq, NVR_size
        real(dp), intent(in):: XPR(3, NVR_size), QPR(neq +1, NVR_size), UPR(3, NVR), GPR(3, NVR)
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
        close(10)
    end subroutine write_particles

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
        integer                          :: comp_level = 1
        type(hdf5_file)                  :: h5f
        character(len=MAX_STRING_LENGTH) :: filout1

        ! Create the filename for the output
        write(filout1, '(A,I5.5,A)') trim(case_folder)//trim(particle_folder) , &
                NTIME, trim(particle_output_file) // ".h5" 
        write (dummy_string, "(A)") achar(9)//'Writing particles to HDF5 file: '//trim(filout1)
        call vpm_print(dummy_string, nocolor, 2)
        ! Open the HDF5 file for writing
        if (present(compression)) then
            comp_level = compression
        end if
        call h5f%open(trim(filout1), action='w', comp_lvl = comp_level)

        ! Write attributes to the HDF5 file
        call h5f%writeattr('/','NTIME', NTIME)
        call h5f%writeattr('/','NVR', NVR)
        call h5f%writeattr('/','neq', neq)
        call h5f%writeattr('/','NVR_size', NVR_size)

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
        use h5fortran

        integer, intent(in)              :: NTIME
        integer, intent(in)              :: NN_in(3), NNbl_in(6)
        integer, intent(in)              :: neqpm
        real(dp), intent(in)             :: RHS(neqpm, NN_in(1), NN_in(2), NN_in(3)), &
                                            SOL(neqpm, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in)             :: velocity(3, NN_in(1), NN_in(2), NN_in(3))
        real(dp), intent(in), optional   :: deformation(3, NN_in(1), NN_in(2), NN_in(3))
        integer, optional                :: compression 
        integer                          :: comp_level = 4
        
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
        call h5f%create('/X', H5T_NATIVE_REAL, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
        call h5f%create('/Y', H5T_NATIVE_REAL, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
        call h5f%create('/Z', H5T_NATIVE_REAL, dset_dims=[NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])

        call h5f%create('/VELOCITY', H5T_NATIVE_REAL, dset_dims=[3, NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])

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
        call h5f%write('/VEL', velocity(1:3, NXs:NXf, NYs:NYf, NZs:NZf))
        
        ! Write the deformation field to the HDF5 file
        if (present(deformation)) then
            call h5f%create('/DEFORM', H5T_NATIVE_REAL, dset_dims=[3, NXf - NXs + 1, NYf - NYs + 1, NZf - NZs + 1])
            call h5f%write('/DEFORM', deformation(1:3, NXs:NXf, NYs:NYf, NZs:NZf)) 
        end if

        ! Close HDF5 file
        call h5f%close()
        
    end subroutine write_pm_solution_hdf5

end module file_io