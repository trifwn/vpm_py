!--------------------------------------------------------------------------------
!>@file
!>@brief Poisson solver library.This source code file(pmlib.f90) contains the library
!!for solving poisson problem on a structured grid with constant DX,DY,DZ
!--------------------------------------------------------------------------------

!>@brief This module defines the variables that will be used internally in the library
!>       All variables are private
module pmlib
   use base_types, only : dp
   use constants, only : pi, pi2, pi4
   implicit none
   real(dp), save                :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
   real(dp), save                :: DXpm, DYpm, DZpm

   integer, save                 :: NVR, NXpm, NYpm, NZpm, ND
   integer, save                 :: NXs_bl, NYs_bl, NXf_bl, NYf_bl, NZs_bl, NZf_bl, NBlocks

   integer, save                 :: nbound, levmax
   integer, save                 :: itree, ibctyp 

   !Here pointers are defined which will be assigned in the external data to save up space
   real(dp), pointer             :: SOL_pm(:, :, :, :), RHS_pm(:, :, :, :), QP(:, :), XP(:, :)

   real(dp), allocatable         :: SOL_0_pm(:,:,:,:), source_bound(:,:),x_s(:,:),y_s(:,:),z_s(:,:),d_s(:)
   real(dp), allocatable, save   :: source_bound_lev(:, :, :), xs_lev(:, :), ys_lev(:, :), zs_lev(:, :), ds_lev(:, :)
   integer,  allocatable, save   :: nbound_lev(:), ilev_t(:, :, :)

   private ::XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm, DXpm, DYpm, DZpm
   private ::NVR, NXpm, NYpm, NZPm, ND
   private ::NXs_bl, NYs_bl, NXf_bl, NYf_bl, NZs_bl, NZf_bl, NBlocks
   private ::nbound, ilev_t
   private ::SOL_pm, RHS_pm,QP, XP,SOL_0_pm
   private ::source_bound, x_s, y_s, z_s, d_s
   private ::source_bound_lev, xs_lev, ys_lev, zs_lev, ds_lev

! pinfdomain.f90
   interface infdomain
      module subroutine infdomain(neqs, neqf)
         Implicit None
         integer, intent(in) :: neqs, neqf
      end subroutine infdomain
   end interface
   interface build_level_nbound
      module subroutine build_level_nbound(NXs, NXf, NYs, NYf, neqs, neqf)
         Implicit None
         integer, intent(in)     :: NXs, NXf, NYs, NYf, neqs, neqf
      end subroutine build_level_nbound
   end interface
   interface infdomain_3D
      module subroutine infdomain_3D(neqs, neqf)
         Implicit None
         integer, intent(in) :: neqs, neqf
      end subroutine infdomain_3D
   end interface

   interface build_level_nbound_3d
      module subroutine build_level_nbound_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
         Implicit None
         integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
      end subroutine build_level_nbound_3d
   end interface

! pmsolve.f90
   interface solve_eq
      module subroutine solve_eq(NXs, NXf, NYs, NYf, neq)
         Implicit None
         integer, intent(in)  :: NXs, NXf, NYs, NYf, neq
      end subroutine solve_eq
   end interface

   interface solve_eq_0
      module subroutine solve_eq_0(NXs, NXf, NYs, NYf, neq)
         Implicit None
         integer, intent(in)   :: NXs, NXf, NYs, NYf, neq
      end subroutine solve_eq_0
   end interface

   interface solve_eq_3D
      module subroutine solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
         Implicit None
         integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
      end subroutine solve_eq_3D
   end interface

   interface solve_eq_0_3D
      module subroutine solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
         Implicit None
         integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
      end subroutine solve_eq_0_3D
   end interface

! pmbound.f90
   interface Bounds2D
      module subroutine Bounds2d(itype, NXs, NXf, NYs, NYf, neqs, neqf)
         Implicit None
         integer, intent(in)  :: itype, NXs, NXf, NYs, NYf, neqs, neqf
      end subroutine Bounds2D
   end interface

   interface Bounds2d_lev
      module subroutine Bounds2d_lev(itype, NXs, NXf, NYs, NYf, neqs, neqf)
         Implicit None
         integer, intent(in):: itype, NXs, NXf, NYs, NYf, neqs, neqf
      end subroutine Bounds2d_lev
   end interface

   interface Bounds3d
      module subroutine Bounds3d(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
         Implicit None
         integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
      end subroutine Bounds3d
   end interface

   interface Bounds3d_lev
      module subroutine Bounds3d_lev(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
         Implicit None
         integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
      end subroutine Bounds3d_lev
   end interface

contains
   !--------------------------------------------------------------------------------
   !>@function
   ! subroutine  pmesh
   !>
   !>@author Papis
   !>
   !>@brief
   !>This is the main subroutine of the Poisson solver library
   !>The input of the subroutine is DSOL_pm,DRHS_pm,DQP,DXP and d velocities.These variables are assigned the
   !>specified pointers.What is also needed is
   !> Dpm,NN,NN_bl,Nblocks,Xbound which define the grid
   !> Dpm(3)   is DX,DY,DZ
   !> NN(3)    is the size of extend domain which coincides with the size of DSOL_pm,DRHS_pm
   !> NN_bl(6) is the size of the original grid (smaller domain) in which the poisson problem will be solved
   !>          (NN_bl(1:3),the starting nodes of the grid(X,Y,Z) with respect to the extended domain)
   !>          (NN_bl(4:6),the last nodes of the grid(X,Y,Z)  with respect to the extended domain)
   !>Nblocks is deprecated
   !>Xbound(6) Xmin,Ymin,Zmin,Xmax,Ymax,Zmax of the extended domain
   !>ibctyp (1 for bio savart law) (2 for infinite domain boundary conditions)
   !>neqs,neqf is an option if we want to solve more than one equation.neqs,neqf shoud coresspond to
   !>DSOL_pm(:,:,:,neqs,neqf)...
   !>iynbc is in case we want to add externally bc's for our equation.0 means normal solve (which means the
   !>solution at boundaries is from ibctyp.1 means keeps the solution at the boundaries at what is defined
   !> externally)
   !>IMPORTANT NOTE: The following library assumes that there are two domains.
   !>                -->the original domain and
   !>                -->an extended domain(that's why NN differ from NN_bl)
   !>                The solution is found on the extended domain.
   !>                IF you don't want an extended domain then set NN_bl(1)=1,NN_bl(4)=NN(1)..NN_bl(2)..
   !>                EXTERNALLY
   !REVISION HISTORY
   !> TODO_dd
   !>
   !>@param [in]
   !>@param [out]
   !>--------------------------------------------------------------------------------
   subroutine pmesh(DSOL_pm, DRHS_pm, DQP, DXP, Xbound, Dpm, NN, NN_bl, ND_in, Nblocks_in, ibctyp_in, &
                    neqs, neqf, iynbc, NVR_in, itree_in, levmax_in)
      ! use parvar, only : XP, QP
      use MPI
      Implicit None
      real(dp), intent(inout), target  :: DSOL_pm(:, :, :, :), DRHS_pm(:, :, :, :), DQP(:, :), DXP(:, :)
      integer, intent(in)              :: ibctyp_in, itree_in, NVR_in, levmax_in, iynbc, neqs, neqf
      real(dp), intent(in)             :: Xbound(6), Dpm(3)
      integer, intent(in)              :: NN_bl(6), NN(3), ND_in, Nblocks_in
      ! Local variables
      integer                          ::  nb, NXs, NYs, NXf, NYf, NZs, NZf, neq

      ibctyp = ibctyp_in
      itree = itree_in
      NVR = NVR_in
      levmax = levmax_in
      ND = ND_in
      NBlocks = Nblocks_in

      ! integer                        :: ierr, my_rank, np, rank
      ! real(dp)                       :: XPM, YPM, velx, vely
      ! real(dp)                       :: xi, yi, ksi1, ksi2, th1, th2, w1, w2
      ! real(dp)                       :: R, DX, DY, GreenF, nv

      !-->Pass the external data to the namespace of the library
      XMIN_pm = Xbound(1); YMIN_pm = Xbound(2); ZMIN_pm = Xbound(3)
      XMAX_pm = Xbound(4); YMAX_pm = Xbound(5); ZMAX_pm = Xbound(6)
      DXpm = Dpm(1); DYpm = Dpm(2); DZpm = Dpm(3)
      !NXs_bl,NXf_bl refer to the starting ending node of the domain we want to solve
      NXpm = NN(1); NYpm = NN(2); NZpm = NN(3)
      NXs_bl = NN_bl(1); NYs_bl = NN_bl(2); NZs_bl = NN_bl(3)
      NXf_bl = NN_bl(4); NYf_bl = NN_bl(5); NZf_bl = NN_bl(6)
      !Assign the pointers to the external data
      SOL_pm => DSOL_pm; RHS_pm => DRHS_pm; !QP => DQP; XP => DXP
      !-->
      ! call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      ! call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
      ! do rank = 0, np - 1
      !    if (rank.eq.my_rank) then
      !       write (*,*) 'Rank:', rank
      !       write (*,*) '------------------------'
      !       write (*,*) maxval(abs(RHS_pm))
      !       write (*,*) '------------------------'
      !    end if
      !    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! enddo

      !iynbc 1 normal poisson solve.(Calculation of bc's is done here)
      if (iynbc .eq. 1) then
         !NXs,NXf is the points we want to calculate boundary conditions
         if (ND .eq. 2) then
            if (ibctyp .eq. 1) then
               !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
               NXs = 1              ! NXs_bl
               NXf = NXpm           !  NXf_bl
               NYs = 1              !NYs_bl
               NYf = NYpm           !NYf_bl
               call Bounds2D(ibctyp, NXs, NXf, NYs, NYf, neqs, neqf)
            else
               !Infinite domain boundary conditions(asume zero bc' at the boundary
               call infdomain(neqs, neqf)
            end if
         else if (ND .eq. 3) then
            if (ibctyp .eq. 1) then
               !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
               NXs = 1           !NXs_bl
               NXf = NXpm        !NXf_bl
               NYs = 1           !NYs_bl
               NYf = NYpm        !NYf_bl
               NZs = 1           !NZs_bl
               NZf = NZpm        !NZf_bl
               call Bounds3D(ibctyp, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
            else
               !Infinite domain boundary conditions(asume zero bc' at the boundary
               call infdomain_3D(neqs, neqf)
            end if
         end if

         !Boundary calculations have finished now we have the boundary conditions at NXs,NXf...
         if (ND .eq. 2) then
            nb = 1
            NXs = 1           !NXs_bl(1)
            NXf = NXpm        !NXf_bl
            NYs = 1           !NYs_bl
            NYf = NYpm        !NYf_bl
            !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
            do neq = neqs, neqf
               call solve_eq(NXs, NXf, NYs, NYf, neq)
            end do
         else
            nb = 1
            NXs = 1         !NXs_bl
            NXf = NXpm      !NXf_bl
            NYs = 1         !NYs_bl
            NYf = NYpm      !NYf_bl
            NZs = 1         !NZs_bl
            NZf = NZpm
            !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
            do neq = neqs, neqf
               call solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
            end do
         end if
      else
         !In case iynbc = 0 then the boundary calculations are asummed ok from the external data and thus
         !going fo the final poisson solve.
         !IMPORTANT:: because the original idea was to have a solution at the smaller domain (extension of
         !the domain used for solving the small).The poisson solver is now solved for the small domain,since
         !boundary conditions are assumed ok.
         if (ND .eq. 2) then
            nb = 1
            NXs = NXs_bl
            NXf = NXf_bl
            NYs = NYs_bl
            NYf = NYf_bl
            !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
            do neq = neqs, neqf
               call solve_eq(NXs, NXf, NYs, NYf, neq)
            end do
         else
            nb = 1
            NXs = NXs_bl
            NXf = NXf_bl
            NYs = NYs_bl
            NYf = NYf_bl
            NZs = NZs_bl
            NZf = NZf_bl
            !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
            do neq = neqs, neqf
               call solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
            end do
         end if
      end if

      ! do rank = 0, np - 1
      !    if (rank.eq.my_rank) then
      !       write (*,*) 'Rank:', rank
      !       write (*,*) '------------------------'
      !       write (*,*) maxval(abs(SOL_pm))
      !       write (*,*) '------------------------'
      !    end if
      !    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! enddo

      nullify (SOL_pm, RHS_pm, QP, XP)
   
      
   end subroutine pmesh
   ! contains
   !    include 'pmsolve.f90'
   !    include 'pmbound.f90'
   !    include 'pinfdomain.f90'

!--------------------------------------------------------------------------------
!>@function
! subroutine  definepm
!>
!>@author Papis
!>
!>@brief
!>definepm  defines the characteristics of the poisson solver grid
!> the Input  is Dpm(3) (DX,DY,DZ)
!>               ndum   (the number of Dpm's which the domain will be extended see above)
!>               Xbound(6) Xmin,Ymin,Zmin,Xmax,Ymax,Zmax of the original domain
!> the Output is Xbound(6) but of the extended domain
!>               NN(3) (NX,NY,NZ) the nodes of the extended domain (1-->NXpm)
!>               NN_bl(6) the starting and ending nodes of the original domain(defined by xbound at input)
!>               (NXs,NYs,NZs,NXf,NYf,NZf)
!>
!REVISION HISTORY
!> TODO_dd
!>
!>@param [in]
!>@param [out]
!--------------------------------------------------------------------------------
   subroutine definepm(itype, Xbound, Dpm, ND, ndum, nsize, NN, NN_bl)
      Implicit None
      integer, intent(in)     :: itype, ND, nsize(3)
      integer, intent(in)     :: ndum
      real(dp), intent(inout) :: Xbound(6)
      real(dp), intent(inout) :: Dpm(3)
      integer, intent(out)    :: NN(3), NN_bl(6)
      integer                 :: ndum_new(3), nn1, nn2
      ! real(dp)              :: Xbound_old(6)

      !-> Define Pmesh X,Y,Z min/max boundaries
      if (ND .eq. 2) then
         Xbound(6) = 0
         Xbound(3) = 0
      end if

      ndum_new = ndum
      
      select case(itype)
      case(1)
         !Itype 1:
         ! extends the domain by ndum_new cells
         Xbound(1) = Xbound(1) - ((ndum_new(1))*Dpm(1))
         Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))

         Xbound(2) = Xbound(2) - ((ndum_new(2))*Dpm(2))
         Xbound(5) = Xbound(5) + ((ndum_new(2))*Dpm(2))
         
         !Find number of nodes with the Dpm given from input
         NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/(Dpm(2)))) + 1

         ! Do the same for the 3d direction
         if (ND .eq. 3) then
            Xbound(3) = Xbound(3) - ((ndum_new(3))*Dpm(3))
            Xbound(6) = Xbound(6) + ((ndum_new(3))*Dpm(3))
            NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/(Dpm(3)))) + 1
         end if
      case(2)
         !Itype 2:
         ! extends the domain by ndum_new cells and changes Dpm so that the number
         !of cells are divided exactly by nsize
         Xbound(1) = Xbound(1) - ((ndum_new(1))*Dpm(1))
         Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))

         Xbound(2) = Xbound(2) - ((ndum_new(2))*Dpm(2))
         Xbound(5) = Xbound(5) + ((ndum_new(2))*Dpm(2))

         if (ND .eq. 3) then
            Xbound(3) = Xbound(3) - ((ndum_new(3))*Dpm(3))
            Xbound(6) = Xbound(6) + ((ndum_new(3))*Dpm(3))
         end if

         !Find number of nodes with the Dpm given from input
         NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/(Dpm(2)))) + 1
         if (ND .eq. 3) NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/(Dpm(3)))) + 1
         NN(1) = NN(1) + nsize(1) - mod(NN(1) - 1, nsize(1))
         NN(2) = NN(2) + nsize(2) - mod(NN(2) - 1, nsize(2))
         if (ND .eq. 3) NN(3) = NN(3) + nsize(3) - mod(NN(3) - 1, nsize(3))
         Dpm(1) = (abs(Xbound(4) - Xbound(1))/(NN(1) - 1))
         Dpm(2) = (abs(Xbound(5) - Xbound(2))/(NN(2) - 1))
         if (ND .eq. 3) Dpm(3) = (abs(Xbound(6) - Xbound(3))/(NN(3) - 1))
         ! write(*,*) 'New Dpm(1),Dpm(2),Dpm(3)'
         ! write(*,*) Dpm(1),Dpm(2),Dpm(3)
         ! write(*,*) NN
      case(4)
         !Itype 4:
         ! extends the domain by ndum_new cells 
         ! and adds dummy cells at both directions
         Xbound(1) = Xbound(1) - ((ndum_new(1))*Dpm(1))
         Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))

         NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         ndum_new(1) = nsize(1) - mod(NN(1) - 1, nsize(1))
         !if (mod(ndum_new(1),2).ne.0) then
         !   write(*,*) 'error sizes',ndum_new,nsize(1),NN(1)
         ! !   stop
         !endif
         !if (mod(ndum_new(2),2).ne.0) then
         !   write(*,*) 'error sizes',ndum_new,nsize(2),NN(2)
         ! !   stop
         !endif
         !if (mod(ndum_new(3),2).ne.0.and.ND.eq.3) then
         !   write(*,*) 'error sizes',ndum_new,nsize(3),NN(3)
         ! !   stop
         !endif
         Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))
         ndum_new(1) = ndum_new(1) + ndum
         NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
      case(3)
         !Itype 3: 
         ! extends the domain by ndum_new cells 
         ! and adds dummy cells at both directions so that the total cells are divided by nsize
         Xbound(1) = Xbound(1) - ((ndum_new(1))*Dpm(1))
         Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))

         Xbound(2) = Xbound(2) - ((ndum_new(2))*Dpm(2))
         Xbound(5) = Xbound(5) + ((ndum_new(2))*Dpm(2))
         NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/(Dpm(2)))) + 1

         if (ND .eq. 3) then
            Xbound(3) = Xbound(3) - ((ndum_new(3))*Dpm(3))
            Xbound(6) = Xbound(6) + ((ndum_new(3))*Dpm(3))
            NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/(Dpm(3)))) + 1
         end if
         
         ndum_new(1) = nsize(1) - mod(NN(1) - 1, nsize(1)) ! X direction
         ndum_new(2) = nsize(2) - mod(NN(2) - 1, nsize(2)) ! Y direction
         ndum_new(3) = nsize(3) - mod(NN(3) - 1, nsize(3)) ! Z direction

         if (mod(ndum_new(1), 2) .eq. 0) then
            ndum_new(1) = ndum_new(1)/2
            Xbound(1) = Xbound(1) - ((ndum_new(1))*Dpm(1))
            Xbound(4) = Xbound(4) + ((ndum_new(1))*Dpm(1))
            ndum_new(1) = ndum_new(1) + ndum
            NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         else
            nn1 = mod(ndum_new(1), 2)
            nn2 = int(ndum_new(1))/int(2)
            Xbound(1) = Xbound(1) - ((nn1 + nn2)*Dpm(1))
            Xbound(4) = Xbound(4) + ((nn2)*Dpm(1))
            NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/(Dpm(1)))) + 1
         end if

         if (mod(ndum_new(2), 2) .eq. 0) then
            ndum_new(2) = ndum_new(2)/2
            Xbound(2) = Xbound(2) - ((ndum_new(2))*Dpm(2))
            Xbound(5) = Xbound(5) + ((ndum_new(2))*Dpm(2))
            ndum_new(2) = ndum_new(2) + ndum !add the initial dummy cells
            NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/(Dpm(2)))) + 1
         else
            nn1 = mod(ndum_new(2), 2)
            nn2 = int(ndum_new(2))/int(2)
            Xbound(2) = Xbound(2) - ((nn1 + nn2)*Dpm(2))
            Xbound(5) = Xbound(5) + ((nn2)*Dpm(2))
            NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/(Dpm(2)))) + 1
         end if

         if (ND .eq. 3) then
            if (mod(ndum_new(3), 2) .eq. 0) then
               ndum_new(3) = ndum_new(3)/2
               Xbound(3) = Xbound(3) - ((ndum_new(3))*Dpm(3))
               Xbound(6) = Xbound(6) + ((ndum_new(3))*Dpm(3))
               ndum_new(3) = ndum_new(3) + ndum !add the initial dummy cells
               NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/(Dpm(3)))) + 1
            else
               nn1 = mod(ndum_new(3), 2)
               nn2 = int(ndum_new(3))/int(2)
               Xbound(3) = Xbound(3) - ((nn1 + nn2)*Dpm(3))
               Xbound(6) = Xbound(6) + ((nn2)*Dpm(3))
               NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/(Dpm(3)))) + 1
            end if
         endif   
      end select
      
      !2d do not have Z nodes
      if (ND .eq. 2) NN(3) = 1
      
      !Define the nodes which the original domain starts (corresponding to Xbound_in)
      !NN_bl(1) = ndum_new(1)+ ndum + 1
      !NN_bl(1) = ndum_new(1)+ ndum + 1
      NN_bl(1) = ndum + 1
      NN_bl(4) = NN(1) - ndum
      
      !NN_bl(2) = ndum_new(2)+ ndum + 1
      !NN_bl(5) = NN(2) - ndum
      NN_bl(2) = ndum + 1
      NN_bl(5) = NN(2) - ndum
      
      !2d do not have Z nodes
      if (ND .eq. 2) then
         NN(3) = 1
         NN_bl(3) = 1
         NN_bl(6) = 1
      else if (ND .eq. 3) then
         !NN_bl(3) = ndum_new(3) + ndum + 1
         !NN_bl(6) = NN(3) - ndum
         NN_bl(3) = ndum + 1
         NN_bl(6) = NN(3) - ndum
      end if
   end subroutine definepm

end module pmlib
