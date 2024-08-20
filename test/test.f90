module test_mod
   double precision, allocatable, target:: XPR(:, :), QPR(:, :), UPR(:, :), GPR(:, :), &
                                           XPO(:, :), QPO(:, :), &
                                           XP_in(:, :), QP_in(:, :)
   double precision, pointer:: velx(:, :, :), vely(:, :, :), velz(:, :, :)
   double precision, pointer:: RHS_pm_in(:, :, :, :), RHS_pm_out(:, :, :, :)
   integer, allocatable    :: qflag(:)
   integer :: NVR_ext, NVR_ext_init
end module test_mod

Program test_pm
   use pmgrid, only:    XMIN_pm, NXs_bl, DXpm, DYpm, DZpm, EPSVOL
   use vpm_vars, only:  mrem, interf_iproj, ncoarse, ncell_rem, iynslice, &
                        IPMWRITE, idefine, IPMWSTART, IPMWSTEPS
   use vpm_size, only:  st, et, NREMESH, iyntree, ilevmax, ibctyp, NBI, &
                        NBJ, NBK
   use test_mod, only:  XPR, QPR, UPR, GPR, NVR_ext, &
                        QPO, XPO, Qflag, &
                        QP_in, XP_in, RHS_pm_in, &
                        velx, vely, velz
   use pmeshpar, only:  IDVPM
   use parvar, only:    QP, XP, UP, NVR
   use openmpth, only:  OMPTHREADS
   use vpm_lib, only:   vpm, remesh_particles_3d, write_particles, write_pm_solution
   use MPI

   Implicit None
   double precision              :: Vref, NI_in, DT_in, RMETM, FACDEF, T, &
                                    XMIN, XMAX, UINF(3)
   integer                       :: NVR_MAX
   integer                       :: my_rank, np, ierr, i, neq, j, TMAX
   ! integer                       :: debug_switch
   ! external                      :: sleep

   call MPI_INIT(ierr)
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   ! ! DEBUG SWITCH USED TO ATTACH TO DEBUGGER
   ! debug_switch = 0
   ! ! While loop for debugging
   ! do while (debug_switch .eq. 0)
   !    if (my_rank .eq. 0) then
   !       print *, 'Debugging'
   !    end if
   !    call sleep(1)
   !    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !    call MPI_BCAST(debug_switch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   ! end do
   ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   ! call MPI_BCAST(debug_switch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   print *, 'Processor ', my_rank, ' of ', np, ' started'
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! READ SETTINGS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   open (1, file='pm.inp')
   read (1, *) DXpm, DYpm, DZpm     ! CELL SIZES
   read (1, *) interf_iproj         ! INTERFACE PROJECTION FUNCTION 1 , 2, 3 , 4
   !NOTE PEZEI NA EINAI GIA DOMAIN CIZE
   read (1, *) ibctyp               ! 1 - PERIODIC, 2 - INFLOW, 3 - OUTFLOW, 4 - WALL
   read (1, *) IDVPM                ! Variable/Constant Volume(0,1)
   read (1, *) EPSVOL               ! CUT-OFF VOLUME (?)
   read (1, *) ncoarse              ! NUMBER OF FINE CELLS PER COARSE CELL per dir
   read (1, *) NBI, NBJ, NBK        !  NBI x NBJ x NBK = NUM OF PROCESSORS (NP)
   read (1, *) NREMESH, ncell_rem   ! 0: NO REMESHING, 1: REMESHING, ncell_rem: PARTICLE PER CELL
   read (1, *) iyntree, ilevmax     ! 1: TREE 0: NO TREE, 3: NUMB OF SUBDIVISION (2^3)
   read (1, *) OMPTHREADS           ! 1 - OPENMP THREADS
   read (1, *) idefine              ! 0: FREE GRID, 1: FIXED GRID
   read (1, *) iynslice             !
   read (1, *) IPMWRITE             ! SAVING PARAMETER
   if (IPMWRITE .gt. 10) stop       !maximume value of writes equal to 10
   if (IPMWRITE .GT. 0) then
      do i = 1, IPMWRITE            !max value 10
         read (1, *) IPMWSTART(i), IPMWSTEPS(i) ! START AND FREQ OF WRITING
      end do
   end if
   close (1)
   if (my_rank .eq. 0) then
      print *, 'Inputs read:'
      print *, achar(9), 'DXpm=', DXpm
      print *, achar(9), 'DYpm=', DYpm
      print *, achar(9), 'DZpm=', DZpm
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   NI_in = -0.1         ! Viscosity
   DT_in = 0.05         ! =dx/U
   NVR_MAX = 5000        ! MAX NUMBER OF PARTICLES
   neq = 3              ! NUMBER OF EQUATIONS
   mrem = 1             ! REMESHING FACTOR

   UINF = 0;            ! INFLOW VELOCITY
   UINF(1) = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !--- INPUT OF PARTICLES
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (my_rank .eq. 0) then
      print *, 'Reading Particles:'
      open (11, file='particles.bin', form='unformatted', access='sequential', status='old') ! read particles
      print *, 'particles.bin'

      read (11) NVR_ext
      read (11) Vref
      print *, achar(9), 'NVR_ext=', NVR_ext
      print *, achar(9), 'Vref=', Vref

      ! Allocate XPR and QPR arrays based on NVR_ext
      allocate (XPR(3, NVR_ext))
      allocate (QPR(4, NVR_ext))

      ! Read XPR and QPR arrays
      do i = 1, NVR_ext
         read (11) XPR(:, i)
         read (11) QPR(:, i)
         write (*, *) achar(9), 'Particle: i', i
         write (*, *) achar(9), 'XPR=', XPR(:, i)
         write (*, *) achar(9), 'QPR=', QPR(:, i)
      end do
      close (11)

      QPR(1:3, :) = -QPR(1:3, :)*Vref
      QPR(4, :) = QPR(4, :)*Vref
      QPR(neq + 1, :) = Vref
      RMETM = 0.001
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !--- INITIALIZATION VPM
   call vpm(XPR, QPR, UPR, GPR, NVR_ext, neq, 0, RHS_pm_in, velx, vely, velz, 0, NI_in, NVR_MAX)
   !--- END INITIALIZATION VPM

   !------------ Remeshing ----------------
   ! We remesh the particles in order to properly distribute them in the domain
   if (my_rank .eq. 0) then
      st = MPI_WTIME()
      write (*, *) achar(27)//'[1;31mRemeshing '//achar(27)//'[0m'
   end if
   call remesh_particles_3d(1, XPR, QPR, NVR_EXT)
   if (my_rank .eq. 0) then
      et = MPI_WTIME()
      write (*, *) achar(9), 'Remeshing took', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !------------ End Remeshing ----------------

   ! ------- Allocate memory for all particles
   if (my_rank .eq. 0) then
      allocate (QPO(1:neq + 1, NVR_ext))     ! OUTPUT STRENGTH
      allocate (XPO(1:3, NVR_ext))           ! OUTPUT POSITION
      
      allocate (QP_in(1:neq + 1, NVR_ext))   ! INPUT STRENGTH
      allocate (XP_in(1:3, NVR_ext))         ! INPUT POSITION

      allocate (Qflag(NVR_ext))
      Qflag = 0 ! 

      QPO = QPR
      XPO = XPR
      
      XMIN = XMIN_pm + (NXs_bl(1) + 4 - 1)*DXpm
      XMAX = maxval(XPO(1, :))
      XPO(1, :) = XPO(1, :) - (XMAX - XMIN)
      
      QP_in = QPO
      XP_in = XPO
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ------- End Allocate memory for all particles

   ! Initialize VPM to call define_size and definepm
   call vpm(XPR, QPR, UPR, GPR, NVR_ext, neq, 0, RHS_pm_in, velx, vely, velz, 0, NI_in, NVR_MAX)

   !--- MAIN LOOP
   T = 0
   TMAX = 100
   do i = 1, TMAX
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_rank .eq. 0) then
         write (*, *)
         write (*, *) "---------------------------------"
         write (*, *) achar(27)//'[1;92mITERATION= ', i, ' of', TMAX, achar(27)//'[0m'
         write (*, *) achar(27)//'[1;92mT=', DT_in*i, achar(27)//'[0m'
         write (*, *) achar(27)//'[1;92mDT=', DT_in, achar(27)//'[0m'
         write (*, *) "---------------------------------"
      end if
      !get velocities and deformations
      T = DT_in

      !--- ALLOCATIONS FOR ALL PARTICLES and sources
      if (my_rank .eq. 0) then
         if (allocated(UPR)) deallocate(UPR)
         if (allocated(GPR)) deallocate(GPR)
         ! Only particles
         allocate (UPR(3, NVR_EXT))
         allocate (GPR(3, NVR_EXT))
         UPR = 0; GPR = 0
      end if

      !--- VPM GETS VELOCITIES AND DEFORMATIONS FROM THE PM SOLUTION
      call vpm(XPR, QPR, UPR, GPR, NVR_EXT, neq, 2, RHS_pm_in, velx, vely, velz, i, NI_in, NVR_MAX)
      !--- END VPM GETS VELOCITIES AND DEFORMATIONS FROM THE PM SOLUTION

      ! WRITE PARTICLES
      if (my_rank .eq. 0) then
         write (*, *) achar(27)//'[1;31mWriting Particles '//achar(27)//'[0m'
         st = MPI_WTIME()
         call write_particles(i, XP, UP, QP, NVR)
         et = MPI_WTIME()
         write (*, *) achar(9), 'VPM: Writing Particles', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if
      ! END WRITE PARTICLES

      ! CONVECTING PARTICLES
      if (my_rank .eq. 0) then
         write (*, *) achar(27)//'[1;31mCONVECTING '//achar(27)//'[0m'
         write (*, *) achar(9), "Max Particle Velocity is now:", maxval(abs(UPR))

         st = MPI_WTIME()

         ! MOVE PARTICLES
         !!$omp parallel private(j,FACDEF)
         !!$omp do
         do j = 1, NVR
            XPR(1:3, j) = XPR(1:3, j) + (UPR(1:3, j) + UINF(1:3))*DT_in
            FACDEF = 1.
            !OMET   = sqrt ( QPR(1,j)**2 + QPR(2,j)**2 + QPR(3,j)**2 )
            !OG     = sqrt ( GPR(1,j)**2 + GPR(2,j)**2 + GPR(3,j)**2 )
            ! if (OG.ne.0.) then                                                      !RMETM
            !    if (OMET.gt.0.001)  FACDEF = OMET*MIN(RMETM,DT_in*OG/OMET)/OG/DT_in  !RMETM
            ! endif                                                                   !RMETM

            ! DEFORMATION OF PARTICLES
            QPR(1:3, j) = QPR(1:3, j) - FACDEF*GPR(1:3, j)*DT_in
            !minus beacuse defromation is negative

            !if(mod(j,10).eq.0)write(28,*) GPR(1:3,j)
            !if(mod(j,10).eq.0)write(29,*) UPR(1:3,j)
         end do
         !!$omp enddo
         !!$omp end parallel
         et = MPI_WTIME()
         write (*, *) achar(9), 'Convecting Particles:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         !    call move_par_out(DT_in)
         st = MPI_WTIME()
         ! PROBABLY MOVES PARTICLES IN AND OUT OF THE DOMAIN
         write (*, *) achar(27)//'[1;31mMoving Particles in/out '//achar(27)//'[0m'
         call find_par_in(T, UINF(1))
         call find_par_out
         et = MPI_WTIME()
         write (*, *) achar(9), 'Moving Particles IN/OUT', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if
      ! END CONVECTING PARTICLES

      !--- VPM INITIALIZATION AND REMESHING
      ! Used to redefine the sizes
      call vpm(XPR, QPR, UPR, GPR, NVR_EXT, neq, 0, RHS_pm_in, velx, vely, velz, i, NI_in, NVR_MAX)

      if (my_rank .eq. 0) then
         st = MPI_WTIME()
         write (*, *) achar(27)//'[1;31mRemeshing '//achar(27)//'[0m'
      end if
      
      ! Remesh the particles
      ! call remesh_particles_3d(1, XPR, QPR, NVR_EXT)
      if (my_rank.eq.0) then
         write (*, *) achar(27)//'[1;31m AFTER Remeshing '//achar(27)//'[0m'
         write(*,*) achar(9), 'NVR_EXT', NVR_EXT
         write(*,*) achar(9), 'NVR', NVR
      end if
      ! ! BCAST NVR_EXT
      call MPI_BCAST(NVR_EXT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! !--- END VPM INITIALIZATION AND REMESHING


      !--- VPM DIFFUSION
      !call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,5,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
      !if (my_rank.eq.0) then
      !   write(*,*) maxval(GPR(:,:))
      !    do j= 1,NVR_ext
      !        QPR(1:3,j) = QPR(1:3,j)  -GPR(1:3,j) * DT_in
      !    enddo
      !endif
      !get velocities and deformation
   end do
   !--- END MAIN LOOP
   
   call MPI_FINALIZE(ierr)
end Program test_pm

Subroutine find_par_out
   use pmgrid, only: DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm, &
                     NXf_bl, NXs_bl, NYf_bl, NYs_bl, NZf_bl, NZs_bl, &
                     NXpm, NYpm, NZpm
   use vpm_vars, only: interf_iproj, neqpm, mrem
   use test_mod, only: NVR_ext, XPR, QPR

   Implicit None
   integer  :: neq
   integer  :: i, NVR_in, NVR_out, NVR_out_max
   integer, allocatable :: NVR_projout(:)
   double precision             :: XMAX, XMIN, YMAX, YMIN, ZMIN, ZMAX, EPSX, EPSY, EPSZ
   double precision             :: xc, yc, zc, dx, dy, dz
   double precision, allocatable :: XP_out(:, :), QP_out(:, :), XP_tmp(:, :), QP_tmp(:, :)
   ! character *50                ::filout

   !an eps for arithmetic reasons
   EPSX = 0.01*DXpm
   EPSY = 0.01*DYpm
   EPSZ = 0.01*DZpm
   XMAX = XMIN_pm + (NXf_bl(1) - (interf_iproj/2 + 1) - 1)*DXpm - EPSX
   XMIN = XMIN_pm + (NXs_bl(1) + (interf_iproj/2 + 1) - 1)*DXpm + EPSX
   YMAX = YMIN_pm + (NYf_bl(1) - (interf_iproj/2 + 1) - 1)*DYpm - EPSY
   YMIN = YMIN_pm + (NYs_bl(1) + (interf_iproj/2 + 1) - 1)*DYpm + EPSY
   ZMAX = ZMIN_pm + (NZf_bl(1) - (interf_iproj/2 + 1) - 1)*DZpm - EPSZ
   ZMIN = ZMIN_pm + (NZs_bl(1) + (interf_iproj/2 + 1) - 1)*DZpm + EPSZ
   xc = 0.5d0*(XMAX + XMIN)
   yc = 0.5d0*(YMAX + YMIN)
   zc = 0.5d0*(ZMAX + ZMIN)
   dx = XMAX - XMIN
   dy = YMAX - YMIN
   dz = ZMAX - ZMIN
   NVR_out = 0

   ! do i=1,NVR_ext
   !    if (XPR(1,i).lt.XMIN.or.XPR(1,i).gt.XMAX) then
   !         XPR(1,i) = XPR(1,i)  - int( (XPR(1,i) - xc)/(0.5d0*dx))*dx
   !       !NVR_out = NVR_out + 1
   !       !XP(1,i) = XPR(1,i)  -2*(XPR(1,i)-xc)
   !    endif
   !    if (XPR(2,i).lt.YMIN.or.XPR(2,i).gt.YMAX) then
   !        XPR(2,i) = XPR(2,i)  - int( (XPR(2,i) - yc)/(0.5d0*dy))*dy
   !       !XPR(2,i) = XPR(2,i)  -2*(XPR(2,i)-yc)
   !       !NVR_out = NVR_out + 1
   !    endif
   !    if (XPR(3,i).lt.ZMIN.or.XPR(3,i).gt.ZMAX) then
   !        XPR(3,i) = XPR(3,i)  - int( (XPR(3,i) - zc)/(0.5d0*dz))*dz
   !       !XPR(3,i) = XPR(3,i)  -2*(XPR(3,i)-zc)
   !      !NVR_out = NVR_out + 1
   !    endif
   ! enddo

   neq = neqpm + 1
   NVR_out_max = (2*NXpm*NYpm + 2*NYpm*NZpm + 2*NZpm*NXpm)*3*mrem**2
   allocate (XP_out(1:3, NVR_out_max), QP_out(1:neq, NVR_out_max), NVR_projout(NVR_out_max))
   allocate (XP_tmp(1:3, NVR_ext), QP_tmp(1:neq, NVR_ext))
   NVR_projout = 2!interf_iproj
   NVR_out = 0
   NVR_in = 0
   do i = 1, NVR_ext
      if (XPR(2, i) .lt. YMIN .or. XPR(2, i) .gt. YMAX &
            .or. XPR(3, i) .lt. ZMIN .or. XPR(3, i) .gt. ZMAX) then
         NVR_out = NVR_out + 1
         XP_out(1:3, NVR_out) = XPR(1:3, i)
         QP_out(1:neq, NVR_out) = QPR(1:neq, i)
      else
         NVR_in = NVR_in + 1
         XP_tmp(1:3, NVR_in) = XPR(1:3, i)
         QP_tmp(1:neq, NVR_in) = QPR(1:neq, i)
      end if
   end do
   !if (NVR_out.eq.0) RHS_pm_out=0.d0
   write (*, *) achar(9), 'Particles out', NVR_out
   deallocate (XPR, QPR)
   allocate (XPR(3, NVR_in), QPR(neq, NVR_in))
   XPR(1:3, 1:NVR_in) = XP_tmp(1:3, 1:NVR_in); QPR(1:neq, 1:NVR_in) = QP_tmp(1:neq, 1:NVR_in)
   NVR_ext = NVR_in
   !deallocate(XP_tmp,QP_tmp)
   !allocate(ieq(neqpm+1),QINF(neqpm+1))
   !QINF=0.d0
   !do i=1,neqpm+1
   !   ieq(i)=i
   !enddo
   !call project_particles_3D(RHS_pm_out,QP_out,XP_out,NVR_projout,NVR_out,neqpm+1,ieq,neqpm+1,QINF,NVR_out_max)
   !RHS_pm_out(neqpm+1,:,:,:)=DVpm
   !write(*,*) maxval(abs(QP_out(1,:))),maxval(abs(QPR(1,:)))
   !call project_vol3d(RHS_pm_out,neqpm+1,ieq,neqpm+1,1)
   !write(*,*) 'out',NVR_out,NVR_in,maxval(abs(RHS_pm_out(1,:,:,:)))
   deallocate (XP_out, QP_out)
End Subroutine find_par_out

Subroutine find_par_in(T_in, U)
   use vpm_vars, only: neqpm, mrem
   use test_mod, only: Qflag, QPO, QP_in, XPO, XP_in, NVR_ext, XPR, QPR, NVR_ext_init
   use pmgrid, only: XMIN_pm, NXs_bl, DXpm, NXpm, NYpm, NZpm

   !
   ! This subroutine moves the particles in and out of the domain
   ! It also resets the inflow particles
   !

   Implicit None
   double precision, intent(in):: T_in, U
   double precision           :: XO, XMIN
   integer                    :: NVR_in, neq
   integer  :: i, NVR_in_max
   double precision, allocatable :: XP_tmp(:, :), QP_tmp(:, :)

   ! Find the inflow particles
   XMIN = XMIN_pm + (NXs_bl(1) + 4 - 1)*DXpm
   XO = XMIN
   
   ! Q flag is used to determine if a particle is inflow or not
   if (minval(Qflag) .eq. 1) then
      Qflag = 0
      QPO = QP_in; XPO = XP_in
      write (155, *) 'Reseting inflow Particles'
   end if

   ! Real number of equations
   neq = neqpm + 1

   ! Max number of particles that can enter the domain
   NVR_in_max = (2*NXpm*NYpm + 2*NYpm*NZpm + 2*NZpm*NXpm)*3*mrem**2
   allocate (XP_tmp(3, NVR_ext + NVR_in_max), QP_tmp(neq, NVR_ext + NVR_in_max))
   
   XP_tmp(1:3, 1:NVR_ext) = XPR(1:3, 1:NVR_ext)
   QP_tmp(1:neq, 1:NVR_ext) = QPR(1:neq, 1:NVR_ext)
   NVR_in = 0
   do i = 1, NVR_ext_init
      XPO(1, i) = XPO(1, i) + U*T_in
      if (XPO(1, i) .gt. XO .and. qflag(i) .eq. 0) then
         NVR_in = NVR_in + 1
         Qflag(i) = 1
         XP_tmp(1:3, NVR_ext + NVR_in) = XPO(1:3, i)
         QP_tmp(1:neq, NVR_ext + NVR_in) = QPO(1:neq, i)
      end if
   end do

   deallocate (XPR, QPR)
   NVR_ext = NVR_ext + NVR_in
   allocate (XPR(1:3, NVR_ext), QPR(1:neq, NVR_ext))
   XPR = XP_tmp
   QPR = QP_tmp
   deallocate (XP_tmp, QP_tmp)
   write (*, *) achar(9), 'Particles in', NVR_in

   ! ADD SOURCES
   ! allocate (XP_tmp(3, NVR_sources + NVR_in_max))
   ! allocate (QP_tmp(neq, NVR_sources + NVR_in_max))
   ! XP_tmp(1:3, 1:NVR_sources) = XSOUR(1:3, 1:NVR_sources)
   ! QP_tmp(1:neq, 1:NVR_sources) = QSOUR(1:neq, 1:NVR_sources)
   NVR_in = 0
   ! do i = 1, NVR_sources_init
   !    XPO(1, NVR_ext_init + i) = XPO(1, NVR_ext_init + i) + U*T_in
   !    if (XPO(1, NVR_ext_init + i) .gt. XO .and. qflag(NVR_ext_init + i) .eq. 0) then
   !       NVR_in = NVR_in + 1
   !       Qflag(NVR_ext_init + i) = 1
   !       XP_tmp(1:3, NVR_sources + NVR_in) = XPO(1:3, NVR_ext_init + i)
   !       QP_tmp(1:neq, NVR_sources + NVR_in) = QPO(1:neq, NVR_ext_init + i)
   !    end if
   ! end do
   ! deallocate (XSOUR, QSOUR)
   ! NVR_sources = NVR_sources + NVR_in
   ! allocate (XSOUR(1:3, NVR_sources), QSOUR(1:neq, NVR_sources))
   ! XSOUR = XP_tmp
   ! QSOUR = QP_tmp
   ! deallocate (XP_tmp, QP_tmp)

   write (*, *) achar(9), 'Sources in', NVR_in
End Subroutine find_par_in

Subroutine move_par_out(DT)
   use vpm_vars, only: interf_iproj
   use test_mod, only: NVR_ext, XPR, UPR
   use pmgrid, only: DXpm, DYpm, DZpm, XMIN_pm, NXf_bl, NXs_bl, &
                     YMIN_pm, NYf_bl, NYs_bl, NZf_bl, NZs_bl, ZMIN_pm

   ! Subroutine to move particles out of the domain
   ! It moves the particles out of the domain if they are outside the domain
   
   Implicit None
   double precision, intent(in)::DT
   integer  :: i
   double precision             :: XMAX, XMIN, YMAX, YMIN, ZMIN, ZMAX, EPSX, EPSY, EPSZ
   double precision             :: xc, yc, zc, dx, dy, dz
   ! character *50                ::filout

   !an eps for arithmetic reasons
   EPSX = 0.01*DXpm
   EPSY = 0.01*DYpm
   EPSZ = 0.01*DZpm
   XMAX = XMIN_pm + (NXf_bl(1) - (interf_iproj/2 + 1) - 1)*DXpm - EPSX
   XMIN = XMIN_pm + (NXs_bl(1) + (interf_iproj/2 + 1) - 1)*DXpm + EPSX
   YMAX = YMIN_pm + (NYf_bl(1) - (interf_iproj/2 + 1) - 1)*DYpm - EPSY
   YMIN = YMIN_pm + (NYs_bl(1) + (interf_iproj/2 + 1) - 1)*DYpm + EPSY
   ZMAX = ZMIN_pm + (NZf_bl(1) - (interf_iproj/2 + 1) - 1)*DZpm - EPSZ
   ZMIN = ZMIN_pm + (NZs_bl(1) + (interf_iproj/2 + 1) - 1)*DZpm + EPSZ
   xc = 0.5d0*(XMAX + XMIN)
   yc = 0.5d0*(YMAX + YMIN)
   zc = 0.5d0*(ZMAX + ZMIN)
   dx = XMAX - XMIN
   dy = YMAX - YMIN
   dz = ZMAX - ZMIN

   do i = 1, NVR_ext
      if (XPR(2, i) .lt. YMIN .or. XPR(2, i) .gt. YMAX) then
         XPR(2, i) = XPR(2, i) - UPR(2, i)*DT

      end if
      if (XPR(3, i) .lt. ZMIN .or. XPR(3, i) .gt. ZMAX) then
         XPR(3, i) = XPR(3, i) - UPR(3, i)*DT
      end if
   end do

   return
End Subroutine move_par_out

! Subroutine create_sources(UINF)
!    use vpm_size, only: Dpm
!    use vpm_vars, only: interf_iproj, neqpm
!    use test_mod, only: NVR_sources, XSOUR, QSOUR
!    use pmgrid, only: NXf_bl, NXs_bl, NYf_bl, NYs_bl, NZf_bl, NZs_bl, &
!                      DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm

!    Implicit None
!    double precision, intent(in) :: UINF(3)
!    double precision :: X(8), Y(8), Z(8)
!    integer :: nxfin, nxstart, nyfin, nystart, nzfin, nzstart, NXpm1, NYpm1, NZpm1
!    integer :: i, j, k, npar
!    nxfin = NXf_bl(1) - interf_iproj/2
!    nxstart = NXs_bl(1) + interf_iproj/2
!    nyfin = NYf_bl(1) - interf_iproj/2
!    nystart = NYs_bl(1) + interf_iproj/2
!    nzfin = NZf_bl(1) - interf_iproj/2
!    nzstart = NZs_bl(1) + interf_iproj/2
!    NXpm1 = nxfin - nxstart
!    NYpm1 = nyfin - nystart
!    NZpm1 = nzfin - nzstart
!    NVR_sources = NXpm1*NYpm1*2 + NXpm1*NZpm1*2
!    Dpm(1) = DXpm; Dpm(2) = DYpm; Dpm(3) = DZpm
!    allocate (XSOUR(3, NVR_sources), QSOUR(neqpm + 1, NVR_sources))
!    XSOUR = 0
!    QSOUR = 0
!    npar = 0
!    j = nystart
!    do k = nzstart, nzfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)
!          Z(2) = ZMIN_pm + Dpm(3)*(k)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = Y(1)
!          XSOUR(3, npar) = 0.5d0*(Z(1) + Z(2))

!          QSOUR(3, npar) = -UINF(1)*Dpm(1)*Dpm(3)
!       end do
!    end do

!    j = nyfin
!    do k = nzstart, nzfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)
!          Z(2) = ZMIN_pm + Dpm(3)*(k)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = Y(1)
!          XSOUR(3, npar) = 0.5d0*(Z(1) + Z(2))

!          QSOUR(3, npar) = UINF(1)*Dpm(1)*Dpm(3)
!       end do
!    end do

!    k = nzstart
!    do j = nystart, nyfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)
!          Y(2) = YMIN_pm + Dpm(2)*(j)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = 0.5d0*(Y(1) + Y(2))
!          XSOUR(3, npar) = Z(1)

!          QSOUR(2, npar) = UINF(1)*Dpm(1)*Dpm(2)
!       end do
!    end do

!    k = nzfin
!    do j = nystart, nyfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)
!          Y(2) = YMIN_pm + Dpm(2)*(j)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = 0.5d0*(Y(1) + Y(2))
!          XSOUR(3, npar) = Z(1)

!          QSOUR(2, npar) = -UINF(1)*Dpm(1)*Dpm(2)
!       end do
!    end do

!    QSOUR = -QSOUR
! End Subroutine create_sources
