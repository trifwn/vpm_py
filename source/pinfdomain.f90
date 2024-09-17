submodule (pmlib) pinfdomain
   implicit none
contains
   !--------------------------------------------------------------------------------
   !>@file
   !>@brief Find the boundary conditions for the poisson solver using infinite domain assumption
   !------------------------------------------------------------------------!
   !->subroutine infdomain                                                  !
   !  This subroutine calculates Boundary conditions using the infinite     !
   !  Domain Poisson Poblem.                                                !
   ! - First grad^2(phi) = 0 is solved(no boundary conditions)              !
   ! - From this solution th(phi)/thn is constructed at the boundaries      !
   !   which gives the boundary conditions                                  !
   !------------------------------------------------------------------------!
   !----------------------------------------------------------------------------------!
   module subroutine infdomain(neqs, neqf)
      use MPI
      Implicit None
      integer, intent(in)  :: neqs, neqf
      integer              :: NXs, NXf, NYs, NYf, nn, ndum, neq_siz
      integer              :: my_rank, ierr
      integer              :: nworkb, neq, i, j , k 
      real(dp)             :: et, st

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

      !--Solve PM with zero boundary conditions
      neq_siz = neqf
      allocate (SOL_0_pm(neq_siz, NXpm, NYpm, NZpm))
      allocate (source_bound(neq_siz, NXpm*2 + NYpm*2))
      allocate (x_s(2, NXpm*2 + NYpm*2))
      allocate (y_s(2, NXpm*2 + NYpm*2))
      allocate (d_s((NXpm*2 + NYpm*2)))

      !these variables are used in case we want a higher order source definition(linear for example)
      d_s = 0.d0
      x_s = 0.d0
      y_s = 0.d0
      SOL_0_pm = 0.d0
      !sources will reside 2 dummy cells away from the original domain
      source_bound = 0
      !We will solve a poisson problem with zero bc's at the point of calculation of the sources
      NXs = 1     !NXs_bl - ndum
      NXf = NXpm  !NXf_bl + ndum
      NYs = 1     !NYs_bl - ndum
      NYf = NYpm
      !-->Solve for zero boundary conditions
      do neq = neqs, neqf
         call solve_eq_0(NXs, NXf, NYs, NYf, neq)
      end do
      !-->Calculate normal gradient (which defines the sources)
      nbound = 0
      call calc_normalderiv(NXs, NXf, NYs, NYf, neqs, neqf)
      deallocate (SOL_0_pm)
      if (my_rank .eq. 0) st = MPI_WTIME()
      if (itree .eq. 1) then
         call build_level_nbound(NXs, NXf, NYs, NYf, neqs, neqf)
      end if
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (199, *) 'tree struc', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if
      !-->SOL_pm_0 not needed for anything else

      !-->Using the sources calculated above define the correct bc's to the extended domain(BS law)
      NXs = 1        
      NXf = NXpm     
      NYs = 1        
      NYf = NYpm     
      if (my_rank .eq. 0) st = MPI_WTIME()
      if (itree .eq. 1) then
         call Bounds2d_lev(ibctyp, NXs, NXf, NYs, NYf, neqs, neqf)
         deallocate (source_bound_lev, ds_lev, xs_lev, ys_lev, nbound_lev, ilev_t)
         deallocate (source_bound, d_s, x_s, y_s)
      else
         call Bounds2D(ibctyp, NXs, NXf, NYs, NYf, neqs, neqf)
         deallocate (source_bound, d_s, x_s, y_s)
      end if
      if (my_rank .eq. 0) et = MPI_WTIME()
      if (my_rank .eq. 0) write (199, *) 'Bounds', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'

   end subroutine infdomain

   !------------------------------------------------------------------------!
   !->subroutine infdomain                                                  !
   !  This subroutine calculates Boundary conditions using the infinite     !
   !  Domain Poisson Poblem.                                                !
   ! - First grad^2(phi) = 0 is solved(no boundary conditions)              !
   ! - From this solution th(phi)/thn is constructed at the boundaries      !
   !   which gives the boundary conditions                                  !
   !------------------------------------------------------------------------!
   module subroutine infdomain_3d(neqs, neqf)
      use MPI
      Implicit None
      integer, intent(in)  :: neqs, neqf
      integer              :: NXs, NXf, NYs, NYf, NZs, NZf, nn, ndum, neq, nworkb, neq_siz, nb

      integer              :: my_rank, ierr
      real(dp)             :: et, st

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      nworkb = 2*NXpm*NYpm + 2*NXpm*NZpm + 2*NZPm*Nypm
      neq_siz = neqf
      allocate (SOL_0_pm(neq_siz, NXpm, NYpm, NZpm))
      allocate (source_bound(neq_siz, nworkb))
      allocate (x_s(4,nworkb))
      allocate (y_s(4,nworkb))
      allocate (z_s(4,nworkb))
      allocate (d_s(nworkb))
      d_s = 0.d0
      x_s = 0.d0
      y_s = 0.d0
      z_s = 0.d0
      SOL_0_pm = 0.d0
      !sources will reside 2 dummy cells away from the original domain
      ndum = 2**levmax
      nb = 1

      !We will solve a poisson problem with zero bc's at the point of calculation of the sources
      NXs = 1         !NXs_bl(nb) - ndum
      NXf = NXpm      !NXf_bl(nb) + ndum
      NYs = 1         !NYs_bl(nb) - ndum
      NYf = NYpm      !NYf_bl(nb) + ndum
      NZs = 1         !NZs_bl(nb) - ndum
      NZf = NZpm      !NZf_bl(nb) + ndum
      !-->Solve for zero boundary conditions
      do neq = neqs, neqf
         call solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
      end do
      !-->Calculate normal gradient
      nbound = 0
      call calc_normalderiv_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)

      deallocate (SOL_0_pm)
      if (my_rank .eq. 0) st = MPI_WTIME()
      if (itree .eq. 1) then
         call build_level_nbound_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
      end if
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (199, *) 'tree struc', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if
      !-->SOL_pm_0 not needed for anything else

      !-->Using the sources calculated above define the correct bc's to the extended domain(BS law)
      NXs = 1     !NXs_bl
      NXf = NXpm  !NXf_bl
      NYs = 1     !NYs_bl
      NYf = NYpm  !NYf_bl
      NZs = 1     !NZs_bl
      NZf = NZpm  !NZf_bl
      if (my_rank .eq. 0) st = MPI_WTIME()
      if (itree .eq. 1) then
         call Bounds3d_lev(ibctyp, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
         deallocate (source_bound_lev, ds_lev, xs_lev, ys_lev, zs_lev, nbound_lev, ilev_t)
         deallocate (source_bound, d_s, x_s, y_s, z_s)
      else
         call Bounds3d(ibctyp, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
         deallocate (source_bound, d_s, x_s, y_s, z_s)
      end if
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, *) "Calculated boundary conditions with:"
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,I1)") achar(9)//"itree = ", itree
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,I1)") achar(9)//"ibctyp = ", ibctyp
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,I4,A,I4,A,I4)") achar(9)//"For a domain:", NXpm, "x", NYpm, "x", NZpm
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,I3,A,F8.3,A)") achar(9)//"Bounds finished in:", int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, yellow, 2)
         write (199, *) 'Bounds', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if

   end subroutine infdomain_3d
   !-------------------------------------------------------------------------!
   !-> subroutine calc_normalderiv                                           !
   !   This subroutine calculated normal derivate at the boundary of the     !
   !   domain.This together with the appropriate Green function calculate    !
   !   the boundary conditions at the domain boundaries.                     !
   !   For the calculation of the derivative  a fourth order one sided       !
   !   difference approximation is used:                                     !
   !   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
   !   Note that one [X,Y,Z]MIN boundaries its is a forwardd difference (and !
   !   thus with the difference scheme is with a minus but because the       !
   !   the sources definitions points outward the minus sign is cancelled    !
   !   and so the scheme is the same for MIN and MAX  boundaries             !
   !-------------------------------------------------------------------------!
   subroutine calc_normalderiv(NXs, NXf, NYs, NYf, neqs, neqf)
      Implicit None
      integer, intent(in)  :: NXs, NXf, NYs, NYf, neqs, neqf
      integer              :: i, j, j1, i1, neq
      real(dp)             :: a1, a2, a3, a4, a5, psi1, psi2

      !Sources are defined in NXs,NXf,Nys,NYf using SOL_pm0(zero bc solution)
      a1 = 25.d0/12.d0
      a2 = -4.d0
      a3 = 3.d0
      a4 = -4.d0/3.d0
      a5 = 0.25d0
      i = Nxs
      do j = NYs, NYf - 1
         nbound = nbound + 1
         j1 = j + 1

         do neq = neqs, neqf
            psi1 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, 1) + a2*SOL_0_pm(neq, i + 1, j, 1) + a3*SOL_0_pm(neq, i + 2, j, 1) + &
                           a4*SOL_0_pm(neq, i + 3, j, 1) + a5*SOL_0_pm(neq, i + 4, j, 1))
            psi2 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, 1) + a2*SOL_0_pm(neq, i + 1, j1, 1) + a3*SOL_0_pm(neq, i + 2, j1, 1) + &
                           a4*SOL_0_pm(neq, i + 3, j1, 1) + a5*SOL_0_pm(neq, i + 4, j1, 1))
            source_bound(neq, nbound) = 0.5d0*(psi1 + psi2)
         end do

         x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
         x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
         y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm!+ 0.5d0 * DYpm
         y_s(2, nbound) = YMIN_pm + (j)*DYpm!+ 0.5d0 * DYpm
         d_s(nbound) = DYpm

      end do

      !---XMAX BOUNDARY----
      i = NXf
      do j = NYs, NYf - 1
         nbound = nbound + 1
         j1 = j + 1

         do neq = neqs, neqf
            psi1 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, 1) + a2*SOL_0_pm(neq, i - 1, j, 1) + a3*SOL_0_pm(neq, i - 2, j, 1) + &
                           a4*SOL_0_pm(neq, i - 3, j, 1) + a5*SOL_0_pm(neq, i - 4, j, 1))
            psi2 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, 1) + a2*SOL_0_pm(neq, i - 1, j1, 1) + a3*SOL_0_pm(neq, i - 2, j1, 1) + &
                           a4*SOL_0_pm(neq, i - 3, j1, 1) + a5*SOL_0_pm(neq, i - 4, j1, 1))
            source_bound(neq, nbound) = 0.5d0*(psi1 + psi2)
         end do

         x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
         x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
         y_s(1, nbound) = YMIN_pm + (j)*DYpm
         y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
         d_s(nbound) = DYpm

      end do

      j1 = 0

      !---YMIN BOUNDARY----(fn points at the y direction)
      j = NYs
      do i = NXs, NXf - 1
         nbound = nbound + 1
         i1 = i + 1

         do neq = neqs, neqf
            psi1 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, 1) + a2*SOL_0_pm(neq, i, j + 1, 1) + a3*SOL_0_pm(neq, i, j + 2, 1) + &
                           a4*SOL_0_pm(neq, i, j + 3, 1) + a5*SOL_0_pm(neq, i, j + 4, 1))
            psi2 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, 1) + a2*SOL_0_pm(neq, i1, j + 1, 1) + a3*SOL_0_pm(neq, i1, j + 2, 1) + &
                           a4*SOL_0_pm(neq, i1, j + 3, 1) + a5*SOL_0_pm(neq, i1, j + 4, 1))
            source_bound(neq, nbound) = 0.5d0*(psi1 + psi2)
         end do

         x_s(1, nbound) = XMIN_pm + (i)*DXpm!+ 0.5d0 * DXpm
         x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm!+ 0.5d0 * DXpm
         y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm
         y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
         d_s(nbound) = DXpm

      end do

      !---YMAX BOUNDARY----
      j = NYf
      do i = NXs, NXf - 1
         nbound = nbound + 1
         i1 = i + 1

         do neq = neqs, neqf
            psi1 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, 1) + a2*SOL_0_pm(neq, i, j - 1, 1) + a3*SOL_0_pm(neq, i, j - 2, 1) + &
                           a4*SOL_0_pm(neq, i, j - 3, 1) + a5*SOL_0_pm(neq, i, j - 4, 1))
            psi2 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, 1) + a2*SOL_0_pm(neq, i1, j - 1, 1) + a3*SOL_0_pm(neq, i1, j - 2, 1) + &
                           a4*SOL_0_pm(neq, i1, j - 3, 1) + a5*SOL_0_pm(neq, i1, j - 4, 1))
            source_bound(neq, nbound) = 0.5d0*(psi1 + psi2)
         end do

         x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm! + 0.5d0 * DXpm
         x_s(2, nbound) = XMIN_pm + (i)*DXpm! + 0.5d0 * DXpm
         y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm
         y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
         d_s(nbound) = DXpm

      end do

   end subroutine calc_normalderiv

   !-------------------------------------------------------------------------!
   !-> subroutine calc_normalderiv                                           !
   !   This subroutine calculated normal derivate at the boundary of the     !
   !   domain.This together with the appropriate Green function calculate    !
   !   the boundary conditions at the domain boundaries.                     !
   !   For the calculation of the derivative  a fourth order one sided       !
   !   difference approximation is used:                                     !
   !   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
   !   Note that one [X,Y,Z]MIN boundaries its is a forwardd difference (and !
   !   thus with the difference scheme is with a minus but because the       !
   !   the sources definitions points outward the minus sign is cancelled    !
   !   and so the scheme is the same for MIN and MAX  boundaries             !
   !-------------------------------------------------------------------------!
   subroutine calc_normalderiv_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
      Implicit None
      integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
      integer              :: i, j, k, j1, i1, k1, neq
      real(dp)             :: a1, a2, a3, a4, a5, psi1, psi2, psi3, psi4

      a1 = 25.d0/12.d0
      a2 = -4.d0
      a3 = 3.d0
      a4 = -4.d0/3.d0
      a5 = 0.25d0
      !Sources are defined AT CELL CENTERS
      !---XMIN BOUNDARY----(fn points at the x direction)
      i = Nxs
      do k = NZs, NZf - 1
         do j = NYs, NYf - 1
            nbound = nbound + 1
            j1 = j + 1
            k1 = k + 1

            do neq = neqs, neqf
               psi1 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i + 1, j, k) + a3*SOL_0_pm(neq, i + 2, j, k) + &
                              a4*SOL_0_pm(neq, i + 3, j, k) + a5*SOL_0_pm(neq, i + 4, j, k))
               psi2 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, k) + a2*SOL_0_pm(neq, i + 1, j1, k) + a3*SOL_0_pm(neq, i + 2, j1, k) + &
                              a4*SOL_0_pm(neq, i + 3, j1, k) + a5*SOL_0_pm(neq, i + 4, j1, k))
               psi3 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, k1) + a2*SOL_0_pm(neq, i + 1, j1, k1) + a3*SOL_0_pm(neq, i + 2, j1, k1) + &
                              a4*SOL_0_pm(neq, i + 3, j1, k1) + a5*SOL_0_pm(neq, i + 4, j1, k1))
               psi4 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, k1) + a2*SOL_0_pm(neq, i + 1, j, k1) + a3*SOL_0_pm(neq, i + 2, j, k1) + &
                              a4*SOL_0_pm(neq, i + 3, j, k1) + a5*SOL_0_pm(neq, i + 4, j, k1))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(4, nbound) = XMIN_pm + (i - 1)*DXpm

            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(3, nbound) = YMIN_pm + (j)*DYpm
            y_s(4, nbound) = YMIN_pm + (j)*DYpm

            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k)*DZpm
            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DYpm*DZpm

         end do
      end do
      !---XMAX BOUNDARY----
      i = NXf
      do k = NZs, NZf - 1
         do j = NYs, NYf - 1
            nbound = nbound + 1
            j1 = j + 1
            k1 = k + 1

            do neq = neqs, neqf
               psi1 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i - 1, j, k) + a3*SOL_0_pm(neq, i - 2, j, k) + &
                              a4*SOL_0_pm(neq, i - 3, j, k) + a5*SOL_0_pm(neq, i - 4, j, k))
               psi2 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, k) + a2*SOL_0_pm(neq, i - 1, j1, k) + a3*SOL_0_pm(neq, i - 2, j1, k) + &
                              a4*SOL_0_pm(neq, i - 3, j1, k) + a5*SOL_0_pm(neq, i - 4, j1, k))
               psi3 = 1./DXpm*(a1*SOL_0_pm(neq, i, j1, k1) + a2*SOL_0_pm(neq, i - 1, j1, k1) + a3*SOL_0_pm(neq, i - 2, j1, k1) + &
                              a4*SOL_0_pm(neq, i - 3, j1, k1) + a5*SOL_0_pm(neq, i - 4, j1, k1))
               psi4 = 1./DXpm*(a1*SOL_0_pm(neq, i, j, k1) + a2*SOL_0_pm(neq, i - 1, j, k1) + a3*SOL_0_pm(neq, i - 2, j, k1) + &
                              a4*SOL_0_pm(neq, i - 3, j, k1) + a5*SOL_0_pm(neq, i - 4, j, k1))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(4, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm

            y_s(4, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(3, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(2, nbound) = YMIN_pm + (j)*DYpm
            y_s(1, nbound) = YMIN_pm + (j)*DYpm

            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k)*DZpm
            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DYpm*DZpm

         end do
      end do
      j1 = 0

      !---YMIN BOUNDARY----(fn points at the y direction)
      j = NYs
      do k = NZs, NZf - 1
         do i = NXs, NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            k1 = k + 1

            do neq = neqs, neqf
               psi1 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i, j + 1, k) + a3*SOL_0_pm(neq, i, j + 2, k) + &
                              a4*SOL_0_pm(neq, i, j + 3, k) + a5*SOL_0_pm(neq, i, j + 4, k))
               psi2 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, k) + a2*SOL_0_pm(neq, i1, j + 1, k) + a3*SOL_0_pm(neq, i1, j + 2, k) + &
                              a4*SOL_0_pm(neq, i1, j + 3, k) + a5*SOL_0_pm(neq, i1, j + 4, k))
               psi3 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, k1) + a2*SOL_0_pm(neq, i1, j + 1, k1) + a3*SOL_0_pm(neq, i1, j + 2, k1) + &
                              a4*SOL_0_pm(neq, i1, j + 3, k1) + a5*SOL_0_pm(neq, i1, j + 4, k1))
               psi4 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, k1) + a2*SOL_0_pm(neq, i, j + 1, k1) + a3*SOL_0_pm(neq, i, j + 2, k1) + &
                              a4*SOL_0_pm(neq, i, j + 3, k1) + a5*SOL_0_pm(neq, i, j + 4, k1))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i)*DXpm
            x_s(4, nbound) = XMIN_pm + (i)*DXpm

            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(3, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(4, nbound) = YMIN_pm + (j - 1)*DYpm

            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k)*DZpm
            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DXpm*DZpm

         end do
      end do

      !---YMAX BOUNDARY----
      j = NYf
      do k = NZs, NZf - 1
         do i = NXs, NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            k1 = k + 1

            do neq = neqs, neqf
               psi1 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i, j - 1, k) +  &
                              a3*SOL_0_pm(neq, i, j - 2, k) + &
                              a4*SOL_0_pm(neq, i, j - 3, k) + a5*SOL_0_pm(neq, i, j - 4, k))
               psi2 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, k) + a2*SOL_0_pm(neq, i1, j - 1, k) +  &
                              a3*SOL_0_pm(neq, i1, j - 2, k) + &
                              a4*SOL_0_pm(neq, i1, j - 3, k) + a5*SOL_0_pm(neq, i1, j - 4, k))
               psi3 = 1./DYpm*(a1*SOL_0_pm(neq, i1, j, k1) + a2*SOL_0_pm(neq, i1, j - 1, k1) + &
                              a3*SOL_0_pm(neq, i1, j - 2, k1) + &
                              a4*SOL_0_pm(neq, i1, j - 3, k1) + a5*SOL_0_pm(neq, i1, j - 4, k1))
               psi4 = 1./DYpm*(a1*SOL_0_pm(neq, i, j, k1) + a2*SOL_0_pm(neq, i, j - 1, k1) + &
                              a3*SOL_0_pm(neq, i, j - 2, k1) + &
                              a4*SOL_0_pm(neq, i, j - 3, k1) + a5*SOL_0_pm(neq, i, j - 4, k1))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(4, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i)*DXpm
            x_s(1, nbound) = XMIN_pm + (i)*DXpm

            y_s(4, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(3, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(2, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm

            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k)*DZpm
            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DXpm*DZpm

         end do
      end do
      !---ZMIN BOUNDARY----
      k1 = 0
      k = NZs
      do j = NYs, NYf - 1
         do i = NXs, NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            j1 = j + 1

            do neq = neqs, neqf
               psi1 = 1./DZpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i, j, k + 1) + &
                              a3*SOL_0_pm(neq, i, j, k + 2) + &
                              a4*SOL_0_pm(neq, i, j, k + 3) + a5*SOL_0_pm(neq, i, j, k + 4))
               psi2 = 1./DZpm*(a1*SOL_0_pm(neq, i1, j, k) + a2*SOL_0_pm(neq, i1, j, k + 1) + &
                              a3*SOL_0_pm(neq, i1, j, k + 2) + &
                              a4*SOL_0_pm(neq, i1, j, k + 3) + a5*SOL_0_pm(neq, i1, j, k + 4))
               psi3 = 1./DZpm*(a1*SOL_0_pm(neq, i1, j1, k) + a2*SOL_0_pm(neq, i1, j1, k + 1) + &
                              a3*SOL_0_pm(neq, i1, j1, k + 2) + &
                              a4*SOL_0_pm(neq, i1, j1, k + 3) + a5*SOL_0_pm(neq, i1, j1, k + 4))
               psi4 = 1./DZpm*(a1*SOL_0_pm(neq, i, j1, k) + a2*SOL_0_pm(neq, i, j1, k + 1) + &
                              a3*SOL_0_pm(neq, i, j1, k + 2) + &
                              a4*SOL_0_pm(neq, i, j1, k + 3) + a5*SOL_0_pm(neq, i, j1, k + 4))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i)*DXpm
            x_s(4, nbound) = XMIN_pm + (i)*DXpm

            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(2, nbound) = YMIN_pm + (j)*DYpm
            y_s(3, nbound) = YMIN_pm + (j)*DYpm
            y_s(4, nbound) = YMIN_pm + (j - 1)*DYpm

            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DXpm*DYpm

            x_s(1, nbound) = XMIN_pm + (i - 1)*DXpm
            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm!+ 0.5d0 * DYpm
            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm!+ 0.5d0 * DYpm
         end do
      end do
      !---ZMAX BOUNDARY
      k = NZf
      do j = NYs, NYf - 1
         do i = NXs, NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            j1 = j + 1

            do neq = neqs, neqf
               psi1 = 1./DZpm*(a1*SOL_0_pm(neq, i, j, k) + a2*SOL_0_pm(neq, i, j, k - 1) + &
                              a3*SOL_0_pm(neq, i, j, k - 2) + &
                              a4*SOL_0_pm(neq, i, j, k - 3) + a5*SOL_0_pm(neq, i, j, k - 4))
               psi2 = 1./DZpm*(a1*SOL_0_pm(neq, i1, j, k) + a2*SOL_0_pm(neq, i1, j, k - 1) + &
                              a3*SOL_0_pm(neq, i1, j, k - 2) + &
                              a4*SOL_0_pm(neq, i1, j, k - 3) + a5*SOL_0_pm(neq, i1, j, k - 4))
               psi3 = 1./DZpm*(a1*SOL_0_pm(neq, i1, j1, k) + a2*SOL_0_pm(neq, i1, j1, k - 1) + &
                              a3*SOL_0_pm(neq, i1, j1, k - 2) + &
                              a4*SOL_0_pm(neq, i1, j1, k - 3) + a5*SOL_0_pm(neq, i1, j1, k - 4))
               psi4 = 1./DZpm*(a1*SOL_0_pm(neq, i, j1, k) + a2*SOL_0_pm(neq, i, j1, k - 1) + &
                              a3*SOL_0_pm(neq, i, j1, k - 2) + &
                              a4*SOL_0_pm(neq, i, j1, k - 3) + a5*SOL_0_pm(neq, i, j1, k - 4))
               source_bound(neq, nbound) = 0.25d0*(psi1 + psi2 + psi3 + psi4)
            end do

            x_s(4, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(3, nbound) = XMIN_pm + (i - 1)*DXpm
            x_s(2, nbound) = XMIN_pm + (i)*DXpm
            x_s(1, nbound) = XMIN_pm + (i)*DXpm

            y_s(4, nbound) = YMIN_pm + (j - 1)*DYpm
            y_s(3, nbound) = YMIN_pm + (j)*DYpm
            y_s(2, nbound) = YMIN_pm + (j)*DYpm
            y_s(1, nbound) = YMIN_pm + (j - 1)*DYpm

            z_s(4, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(3, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(2, nbound) = ZMIN_pm + (k - 1)*DZpm
            z_s(1, nbound) = ZMIN_pm + (k - 1)*DZpm

            d_s(nbound) = DXpm*DYpm

         end do
      end do
   end subroutine calc_normalderiv_3d

   !This subroutine builds the nbounds_lev matrrix using the values calulated at the finer level.

   module subroutine build_level_nbound(NXs, NXf, NYs, NYf, neqs, neqf)
      Implicit None
      integer, intent(in)     :: NXs, NXf, NYs, NYf, neqs, neqf
      integer                 :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin,&
                                 ires, leafacc, lf
      integer                 :: ncountlev(0:levmax), nj, nn, npre
      integer, allocatable    :: nn_lev(:, :)
      real(dp)                :: x, y, s, source(neqf), xc, yc, sc, sourcec(neqf)
      real(dp), allocatable   :: xs_tmp(:, :), ds_tmp(:), s_tmp(:, :)
      integer                 :: i, j, k, neq 

      leafmax = 1
      leafmax = 2
      allocate (source_bound_lev(nbound, neqs:neqf, 0:levmax)); source_bound_lev = 0.d0
      allocate (xs_lev(nbound, 0:levmax)); xs_lev = 0.d0
      allocate (ys_lev(nbound, 0:levmax)); ys_lev = 0.d0
      allocate (ds_lev(nbound, 0:levmax)); ds_lev = 0.d0
      allocate (ilev_t(nbound, 0:levmax, 6)); ilev_t = 0
      allocate (nbound_lev(0:levmax)); nbound_lev = 0
      allocate (nn_lev(3, 0:levmax)); nn_lev = 0.d0
      allocate (xs_tmp(nbound, 3), ds_tmp(nbound), s_tmp(nbound, 1:neqf))

      xs_tmp(1:nbound, 1) = 0.5d0*(x_s(1, :) + x_s(2, :))
      xs_tmp(1:nbound, 2) = 0.5d0*(y_s(1, :) + y_s(2, :))
      ds_tmp(1:nbound) = d_s(:)
      do neq = neqs, neqf
         s_tmp(1:nbound, neq) = source_bound(neq, 1:nbound)
      end do

      nbound_lev(0) = nbound

      ilev_t = 0
      icount = 0
      leafcount = 0
      leafacc = 0
      istep = 2**levmax
      nbound_lev = 0; ncountlev = 0
      i = Nxs

      !FINER LEVEL 0 -COARSER LEVEL levmax

      !The tree level is built.It supports levels that are not exactly divided.In this case a dummy cell is added
      !in the coarser level and when you encounter that in the tree routine you automatically go to the next level.
      !The tree is defined at cells not at nodes.
      !the basic structure is the ilev_t matrix
      !ilev_t(1:2_ are the corresponding finer level cells for each coarse cell
      !ilev_t(5) are how many cells are there(normally 2,dummy 1)
      !dummy can have more than 1 cell because each coarser level is define as groups of two
      !ilev_t(6) >0 if normal cell,-1 if dummy cell in i direction,-2 in dummy cell in j direction
      !store  how many cells in each level in each direction for later use in building the tree

      do lev = 0, levmax
         istep = 2**lev
         nn_lev(1, lev) = int((NXf - NXs)/istep)
         nn_lev(2, lev) = int((NYf - NYs)/istep)
      !if not divided exactly dummy cell
         if (mod(NXf - NXs, istep) .ne. 0) nn_lev(1, lev) = nn_lev(1, lev) + 1
         if (mod(NYf - NYs, istep) .ne. 0) nn_lev(2, lev) = nn_lev(2, lev) + 1
      end do

      !----------Start for each plane

      do lev = 0, levmax
         istep = 2**lev
         icount = nbound_lev(lev)
         leafcount = nbound_lev(0)
         do j = NYs, NYf - istep, istep
            !the levels are build based on the finer level
            xc = 0; yc = 0; sc = 0; sourcec = 0
            do k = 1, istep
               leafcount = leafcount + 1
               x = xs_tmp(leafcount, 1)
               y = xs_tmp(leafcount, 2)
               s = ds_tmp(leafcount)
               source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)
               xc = x + xc; yc = y + yc; sc = s + sc; sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
            end do
            icount = icount + 1
            xs_lev(icount, lev) = xc/float(istep)
            ys_lev(icount, lev) = yc/float(istep)
            ds_lev(icount, lev) = sc
            source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
            npre = nn_lev(2, lev)
            !for levels greater than 0 find the child cells
            if (lev .gt. 0) then
               nn = nn_lev(2, lev - 1)
               if (j .eq. NYs .and. icount - 1 .eq. nbound_lev(lev)) then
                  nj = nbound_lev(lev - 1) + 1
               !for the first cell of the plane
               else if (j .eq. NYs .and. icount - 1 .gt. nbound_lev(lev)) then
                  !for the first cell of the line
                  nj = ilev_t(icount - 1, lev, 2) + 1
                  !if the previous cell is dummy
                  if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 1) + 1
               else
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if

               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
            end if
            ilev_t(icount, lev, 6) = npre
         end do
         if (mod(NYf - Nys, istep) .ne. 0 .or. NYf - NYs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here'
               STOP
            end if
            !for dummy cells in the ith direction
            icount = icount + 1
            nn = nn_lev(2, lev - 1)
            nj = ilev_t(icount - 1, lev, 2) + 1
            ilev_t(icount, lev, 1) = nj
            ilev_t(icount, lev, 5) = 1
            ilev_t(icount, lev, 6) = -1
         end if
         ncountlev(lev) = icount
      end do
      nbound_lev = ncountlev

      i = Nxf
      do lev = 0, levmax
         istep = 2**lev
         icount = nbound_lev(lev)
         leafcount = nbound_lev(0)
         do j = NYs, NYf - istep, istep
            !for all leaves in
            xc = 0; yc = 0; sc = 0; sourcec = 0
            do k = 1, istep
               leafcount = leafcount + 1
               x = xs_tmp(leafcount, 1)
               y = xs_tmp(leafcount, 2)
               s = ds_tmp(leafcount)
               source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)
               xc = x + xc; yc = y + yc; sc = s + sc; sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
            end do
            icount = icount + 1
            xs_lev(icount, lev) = xc/float(istep)
            ys_lev(icount, lev) = yc/float(istep)
            ds_lev(icount, lev) = sc
            source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
            npre = nn_lev(2, lev)
            if (lev .gt. 0) then
            !---         if(icount.eq.0)(2,lev-1
               nn = nn_lev(2, lev - 1)
               if (j .eq. NYs .and. icount - 1 .eq. nbound_lev(lev)) then
                  nj = nbound_lev(lev - 1) + 1
               else if (j .eq. NYs .and. icount - 1 .gt. nbound_lev(lev)) then
                  nj = ilev_t(icount - 1, lev, 2) + 1
                  if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 1) + 1
               else
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if

               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
            end if
            ilev_t(icount, lev, 6) = npre
         end do
         if (mod(NYf - Nys, istep) .ne. 0 .or. NYf - NYs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here'
               STOP
            end if
            icount = icount + 1
            nn = nn_lev(2, lev - 1)
            nj = ilev_t(icount - 1, lev, 2) + 1
            ilev_t(icount, lev, 1) = nj
            ilev_t(icount, lev, 5) = 1
            ilev_t(icount, lev, 6) = -1
         end if
         ncountlev(lev) = icount
      end do
      nbound_lev = ncountlev

      j = NYs
      do lev = 0, levmax
         istep = 2**lev
         icount = nbound_lev(lev)
         leafcount = nbound_lev(0)
         do i = NXs, NXf - istep, istep
            !for all leaves in
            xc = 0; yc = 0; sc = 0; sourcec = 0
            do k = 1, istep
               leafcount = leafcount + 1
               x = xs_tmp(leafcount, 1)
               y = xs_tmp(leafcount, 2)
               s = ds_tmp(leafcount)
               source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)
               xc = x + xc; yc = y + yc; sc = s + sc; sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
            end do
            icount = icount + 1
            xs_lev(icount, lev) = xc/float(istep)
            ys_lev(icount, lev) = yc/float(istep)
            ds_lev(icount, lev) = sc
            source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
            npre = nn_lev(1, lev)
            if (lev .gt. 0) then
            !---         if(icount.eq.0)(2,lev-1
               nn = nn_lev(1, lev - 1)
               if (i .eq. NXs .and. icount - 1 .eq. nbound_lev(lev)) then
                  nj = nbound_lev(lev - 1) + 1
               else if (i .eq. NXs .and. icount - 1 .gt. nbound_lev(lev)) then
                  nj = ilev_t(icount - 1, lev, 2) + 1
                  if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 1) + 1
               else
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if

               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
            end if
            ilev_t(icount, lev, 6) = npre
         end do
         if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here'
               STOP
            end if
            icount = icount + 1
            nn = nn_lev(1, lev - 1)
            nj = ilev_t(icount - 1, lev, 2) + 1
            ilev_t(icount, lev, 1) = nj
            ilev_t(icount, lev, 5) = 1
            ilev_t(icount, lev, 6) = -1
         end if
         ncountlev(lev) = icount
      end do
      nbound_lev = ncountlev

      j = NYf
      do lev = 0, levmax
         istep = 2**lev
         icount = nbound_lev(lev)
         leafcount = nbound_lev(0)
         do i = NXs, NXf - istep, istep
            !for all leaves in
            xc = 0; yc = 0; sc = 0; sourcec = 0
            do k = 1, istep
               leafcount = leafcount + 1
               x = xs_tmp(leafcount, 1)
               y = xs_tmp(leafcount, 2)
               s = ds_tmp(leafcount)
               source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)
               xc = x + xc; yc = y + yc; sc = s + sc; sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
            end do
            icount = icount + 1
            xs_lev(icount, lev) = xc/float(istep)
            ys_lev(icount, lev) = yc/float(istep)
            ds_lev(icount, lev) = sc
            source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
            npre = nn_lev(1, lev)
            if (lev .gt. 0) then
            !---         if(icount.eq.0)(2,lev-1
               nn = nn_lev(1, lev - 1)
               if (i .eq. NXs .and. icount - 1 .eq. nbound_lev(lev)) then
                  nj = nbound_lev(lev - 1) + 1
               else if (i .eq. NXs .and. icount - 1 .gt. nbound_lev(lev)) then
                  nj = ilev_t(icount - 1, lev, 2) + 1
                  if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 1) + 1
               else
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if

               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
            end if
            ilev_t(icount, lev, 6) = npre
         end do
         if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here'
               STOP
            end if
            icount = icount + 1
            nn = nn_lev(1, lev - 1)
            nj = ilev_t(icount - 1, lev, 2) + 1
            ilev_t(icount, lev, 1) = nj
            ilev_t(icount, lev, 5) = 1
            ilev_t(icount, lev, 6) = -1
         end if
         ncountlev(lev) = icount
      end do
      nbound_lev = ncountlev

      if (nbound_lev(0) .ne. nbound) then
         write (*, *) 'something is very wrong with levels.lev0 bounds should be nbounds'
         STOP

      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      !   do lev=0,levmax
      !      write(filout,'(a,i2.2)') 'lev',lev
      !      open(27,file=filout)
      !      do k=1,nbound_lev(lev)
      !         if(ilev_t(k,lev).eq.1) write(27,*)xs_lev(k,lev,2),ys_lev(k,lev,2)
      !      enddo
      !   enddo
      !stop
   end subroutine build_level_nbound

   !-------------------------------------------------------------------------!
   !-> subroutine calc_normalderiv                                           !
   !   This subroutine calculated normal derivate at the boundary of the     !
   !   domain.This together with the appropriate Green function calculate    !
   !   the boundary conditions at the domain boundaries.                     !
   !   For the calculation of the derivative  a fourth order one sided       !
   !   difference approximation is used:                                     !
   !   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
   !-------------------------------------------------------------------------!
   module subroutine build_level_nbound_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
      use mpi
      Implicit None
      integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
      integer                 :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin, ires, leafacc, lf
      integer                 :: nleaflev, nleafroot, ncountlev(0:levmax), istepj, n1s, n1f, n2s, n2f, nj, npre
      integer                 :: m, l, il, im, mm, nn, iresroot, NNX, NNY, NNZ, NNX0, NNY0, NNZ0
      real(dp)                :: x, y, z, s, source(neqf), xc, yc, zc, sc, sourcec(neqf)
      real(dp), allocatable   :: xs_tmp(:, :), ds_tmp(:), s_tmp(:, :)
      integer, allocatable    :: nn_lev(:, :)
      integer                 :: my_rank, ierr, i, j, k, neq

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      leafmax = 2
      istep = 2**levmax

      allocate (nbound_lev(0:levmax)); nbound_lev = 0
      nbound_lev(0) = nbound
      allocate (source_bound_lev(nbound_lev(0), neqs:neqf, 0:levmax)); source_bound_lev = 0.d0
      allocate (xs_lev(nbound_lev(0), 0:levmax)); xs_lev = 0.d0
      allocate (ys_lev(nbound_lev(0), 0:levmax)); ys_lev = 0.d0
      allocate (zs_lev(nbound_lev(0), 0:levmax)); zs_lev = 0.d0
      allocate (ds_lev(nbound_lev(0), 0:levmax)); ds_lev = 0.d0
      allocate (ilev_t(nbound_lev(0), 0:levmax, 6)); ilev_t = 0
      allocate (nn_lev(3, 0:levmax)); nn_lev = 0.d0
      allocate (xs_tmp(nbound, 3), ds_tmp(nbound), s_tmp(nbound, 1:neqf))

      xs_tmp(1:nbound, 1) = 0.25d0*(x_s(1, 1:nbound) + x_s(2, 1:nbound) + x_s(3, 1:nbound) + x_s(4, 1:nbound))
      xs_tmp(1:nbound, 2) = 0.25d0*(y_s(1, 1:nbound) + y_s(2, 1:nbound) + y_s(3, 1:nbound) + y_s(4, 1:nbound))
      xs_tmp(1:nbound, 3) = 0.25d0*(z_s(1, 1:nbound) + z_s(2, 1:nbound) + z_s(3, 1:nbound) + z_s(4, 1:nbound))
      ds_tmp(1:nbound) = d_s(1:nbound)
      do neq = neqs, neqf
         s_tmp(1:nbound, neq) = source_bound(neq, 1:nbound)
      end do

      !FINER LEVEL 0 -COARSER LEVEL levmax

      !The tree level is built.It supports levels that are not exactly divided.In this case a dummy cell is added
      !in the coarser level and when you encounter that in the tree routine you automatically go to the next level.
      !The tree is defined at cells not at nodes.
      !the basic structure is the ilev_t matrix
      !ilev_t(1:4) are the corresponding finer level cells for each coarse cell
      !dummy can have more than 2 cell because each coarser level is define as groups of four
      !ilev_t(5) are how many cells are there(normally 4,dummy 2)
      !ilev_t(6) >0 if normal cell,-1 if dummy cell in i direction,-2 in dummy cell in j direction
      NNZ = NZf - NZs
      NNY = NYf - NYs
      NNX = NXf - NXs
      icount = 0
      leafcount = 0
      leafacc = 0
      ncountlev = 0

      !store  how many cells in each level in each direction for later use in building the tree
      do lev = 0, levmax
         istep = 2**lev
         nn_lev(1, lev) = int((NXf - NXs)/istep)
         nn_lev(2, lev) = int((NYf - NYs)/istep)
         nn_lev(3, lev) = int((NZf - NZs)/istep)
         !if not divided exactly dummy cell
         if (mod(NXf - NXs, istep) .ne. 0) nn_lev(1, lev) = nn_lev(1, lev) + 1
         if (mod(NYf - NYs, istep) .ne. 0) nn_lev(2, lev) = nn_lev(2, lev) + 1
         if (mod(NZf - NZs, istep) .ne. 0) nn_lev(3, lev) = nn_lev(3, lev) + 1
      end do
      !----------Start for each plane
      nbound_lev = 0
      i = Nxs
      do lev = 0, levmax
         istep = 2**lev
         !n1s,n2s to define the same if's in each plane
         n2s = NZs; n2f = NZf - istep
         n1s = NYs; n1f = NYf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do k = n2s, n2f, istep
            do j = n1s, n1f, istep
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               !the levels are build based on the finer level
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (k - NZs)*NNY + (j - NYs) + (mm - 1)*NNY + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(2, lev)
               !for levels greater than 0 find the child cells
               if (lev .gt. 0) then
                  nn = nn_lev(2, lev - 1)
                  !for the first cell of the plane
                  if (j .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                     !for the first cell of the line
                  else if (j .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     !if the previous cell is dummy
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if

                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4
               end if
               ilev_t(icount, lev, 6) = npre
            end do!j
            if (mod(NYf - NYs, istep) .ne. 0 .or. NYf - NYs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 1'
                  write (*, *) istep, NYf - NYs
                  STOP
               end if
               !for dummy cells in the ith direction
               icount = icount + 1
               nn = nn_lev(2, lev - 1)
               nj = ilev_t(icount - 1, lev, 2) + 1
               !in the ith direction the two cells are nj,nj+nn
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NZf - NZs, istep) .ne. 0 .or. NZf - NZs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here 2'
               write (*, *) istep, NZf - NZs
               STOP
            end if
            !for dummy cells in the jth direction a loop is needed
            do j = n1s, n1f, istep
               !last corner not taken twice
               if (j .eq. n1f .and. mod(NYf - NYs, istep) .ne. 0) cycle
               icount = icount + 1
               if (j .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               !in the jth direction the two cells are nj,nj+1
               nn = nn_lev(2, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      !---XMAX BOUNDARY----
      i = Nxf
      do lev = 0, levmax
         istep = 2**lev
         !n1s,n2s to define the same if's in each plane
         n2s = NZs; n2f = NZf - istep
         n1s = NYs; n1f = NYf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do k = n2s, n2f, istep
            do j = n1s, n1f, istep
               !the levels are build based on the finer level
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (k - NZs)*NNY + (j - NYs) + (mm - 1)*NNY + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(2, lev)
               !for levels greater than 0 find the child cells
               if (lev .gt. 0) then
                  nn = nn_lev(2, lev - 1)
                  !for the first cell of the plane
                  if (j .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                  else if (j .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     !for the first cell of the line
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     !if the previous cell is dummy
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if

                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4

               end if
               ilev_t(icount, lev, 6) = npre
            end do!j
            if (mod(NYf - NYs, istep) .ne. 0 .or. NYf - NYs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 3'
                  write (*, *) istep, NYf - NYs
                  STOP
               end if
               icount = icount + 1
               nj = ilev_t(icount - 1, lev, 2) + 1
               nn = nn_lev(2, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NZf - NZs, istep) .ne. 0 .or. NZf - NZs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here 4'
               write (*, *) istep, NZf - NZs
               STOP
            end if
            do j = n1s, n1f, istep
               if (j .eq. n1f .and. mod(NYf - NYs, istep) .ne. 0) cycle
               icount = icount + 1
               if (j .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               nn = nn_lev(2, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      !---YMIN BOUNDARY----
      j = NYs
      do lev = 0, levmax
         istep = 2**lev
         n2s = NZs; n2f = NZf - istep
         n1s = NXs; n1f = NXf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do k = n2s, n2f, istep
            do i = n1s, n1f, istep
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (k - NZs)*NNX + (i - NXs) + (mm - 1)*NNX + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(1, lev)
               if (lev .gt. 0) then

                  nn = nn_lev(1, lev - 1)
                  if (i .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                  else if (i .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if
                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4
               end if
               ilev_t(icount, lev, 6) = npre
               !if (mod(NXf-NXs,istep).ne.0) ilev_t(icount,lev)=ilev_t(icount,lev) + 1
            end do!j
            if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 5'
                  write (*, *) istep, NXf - NXs
                  STOP
               end if
               icount = icount + 1
               nj = ilev_t(icount - 1, lev, 2) + 1
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NZf - NZs, istep) .ne. 0 .or. NZf - NZs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here 6'
               write (*, *) istep, NZf - NZs
               STOP
            end if
            do i = n1s, n1f, istep
               if (i .eq. n1f .and. mod(NXf - NXs, istep) .ne. 0) cycle
               icount = icount + 1
               if (i .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      !---YMAX BOUNDARY----
      j = NYf
      do lev = 0, levmax
         istep = 2**lev
         n2s = NZs; n2f = NZf - istep
         n1s = NXs; n1f = NXf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do k = n2s, n2f, istep
            do i = n1s, n1f, istep
               !we use l-1 and m-1 because the 1 is already asigned at il,im
               !sources at points
               !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
               !sources at cells
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (k - NZs)*NNX + (i - NXs) + (mm - 1)*NNX + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(1, lev)
               if (lev .gt. 0) then
                  nn = nn_lev(1, lev - 1)
                  if (i .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                  else if (i .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if
                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4
               end if
               ilev_t(icount, lev, 6) = npre
               !if (mod(NXf-NXs,istep).ne.0) ilev_t(icount,lev)=ilev_t(icount,lev) + 1
            end do!j
            if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 7'
                  write (*, *) istep, NXf - NXs
                  STOP
               end if
               icount = icount + 1
               nj = ilev_t(icount - 1, lev, 2) + 1
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NZf - NZs, istep) .ne. 0 .or. NZf - NZs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here 8'
               write (*, *) istep, NZf - NZs
               STOP
            end if
            do i = n1s, n1f, istep
               if (i .eq. n1f .and. mod(NXf - NXs, istep) .ne. 0) cycle
               icount = icount + 1
               if (i .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      !---ZMIN BOUNDARY----
      k = NZs
      do lev = 0, levmax
         istep = 2**lev
         n2s = NYs; n2f = NYf - istep
         n1s = NXs; n1f = NXf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do j = n2s, n2f, istep
            do i = n1s, n1f, istep
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (j - NYs)*NNX + (i - NXs) + (mm - 1)*NNX + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(1, lev)
               if (lev .gt. 0) then
                  nn = nn_lev(1, lev - 1)
                  if (i .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                  else if (i .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if
                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4
               end if
               ilev_t(icount, lev, 6) = npre
               !  if (mod(NXf-NXs,istep).ne.0) ilev_t(icount,lev)=ilev_t(icount,lev) + 1
            end do!j
            if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 9'
                  write (*, *) istep, NXf - NXs
                  STOP
               end if
               icount = icount + 1
               nj = ilev_t(icount - 1, lev, 2) + 1
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NYf - NYs, istep) .ne. 0 .or. NYf - NYs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here 10'
               write (*, *) istep, NYf - NYs
               STOP
            end if
            do i = n1s, n1f, istep
               if (i .eq. n1f .and. mod(NXf - NXs, istep) .ne. 0) cycle
               icount = icount + 1
               if (i .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      !---ZMAX BOUNDARY----
      k = NZf
      do lev = 0, levmax
         istep = 2**lev
         n2s = NYs; n2f = NYf - istep
         n1s = NXs; n1f = NXf - istep
         icount = nbound_lev(lev)
         lf = nbound_lev(0)
         do j = n2s, n2f, istep
            do i = n1s, n1f, istep
               xc = 0; yc = 0; zc = 0; sc = 0; sourcec = 0
               do mm = 1, istep
                  do nn = 1, istep
                     leafcount = lf + (j - NYs)*NNX + (i - NXs) + (mm - 1)*NNX + nn
                     x = xs_tmp(leafcount, 1)
                     y = xs_tmp(leafcount, 2)
                     z = xs_tmp(leafcount, 3)
                     s = ds_tmp(leafcount)
                     source(neqs:neqf) = s_tmp(leafcount, neqs:neqf)

                     xc = x + xc; yc = y + yc; zc = z + zc; sc = s + sc
                     sourcec(neqs:neqf) = source(neqs:neqf)*s + sourcec(neqs:neqf)
                  end do
               end do
               icount = icount + 1
               xs_lev(icount, lev) = xc/float(istep**2)
               ys_lev(icount, lev) = yc/float(istep**2)
               zs_lev(icount, lev) = zc/float(istep**2)
               ds_lev(icount, lev) = sc
               source_bound_lev(icount, neqs:neqf, lev) = sourcec(neqs:neqf)/sc
               npre = nn_lev(1, lev)
               if (lev .gt. 0) then
                  nn = nn_lev(1, lev - 1)
                  if (i .eq. n1s .and. icount - 1 .eq. nbound_lev(lev)) then
                     nj = nbound_lev(lev - 1) + 1
                  else if (i .eq. n1s .and. icount - 1 .gt. nbound_lev(lev)) then
                     nj = ilev_t(icount - 1, lev, 3) + 1
                     if (ilev_t(icount - 1, lev, 6) .lt. 0) nj = ilev_t(icount - 1, lev, 2) + 1
                  else
                     nj = ilev_t(icount - 1, lev, 2) + 1
                  end if
                  nn = nn_lev(1, lev - 1)
                  ilev_t(icount, lev, 1) = nj
                  ilev_t(icount, lev, 2) = nj + 1
                  ilev_t(icount, lev, 3) = nj + nn + 1
                  ilev_t(icount, lev, 4) = nj + nn
                  ilev_t(icount, lev, 5) = 4
               end if
               ilev_t(icount, lev, 6) = npre
               !       if (mod(NXf-NXs,istep).ne.0) ilev_t(icount,lev)=ilev_t(icount,lev) + 1
            end do!j
            if (mod(NXf - NXs, istep) .ne. 0 .or. NXf - NXs .le. istep) then
               if (lev .eq. 0) then
                  write (*, *) 'something is very wrong with levels.lev0 should not be in here 11'
                  write (*, *) istep, NXf - NXs
                  STOP
               end if
               icount = icount + 1
               nj = ilev_t(icount - 1, lev, 2) + 1
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + nn
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -1
            end if
         end do!k
         if (mod(NYf - NYs, istep) .ne. 0 .or. NYf - NYs .le. istep) then
            if (lev .eq. 0) then
               write (*, *) 'something is very wrong with levels.lev0 should not be in here12'
               write (*, *) istep, NYf - NYs
               STOP
            end if
            do i = n1s, n1f, istep
               if (i .eq. n1f .and. mod(NXf - NXs, istep) .ne. 0) cycle
               icount = icount + 1
               if (i .eq. n1s .and. ilev_t(icount - 1, lev, 6) .gt. 0) then
                  nj = ilev_t(icount - 1, lev, 3) + 1
               else!even if ilev_t(icount-1,lev,6).lt.0 I want ilev_t(icount-1,lev,2)
                  nj = ilev_t(icount - 1, lev, 2) + 1
               end if
               nn = nn_lev(1, lev - 1)
               ilev_t(icount, lev, 1) = nj
               ilev_t(icount, lev, 2) = nj + 1
               ilev_t(icount, lev, 5) = 2
               ilev_t(icount, lev, 6) = -2
            end do
         end if
         ncountlev(lev) = icount
      end do!lev
      nbound_lev = ncountlev

      if (nbound_lev(0) .ne. nbound) then
         write (*, *) 'nbound_lev problem', nbound_lev(0), icount, nbound
         stop
      end if

      !do lev=0,levmax
      !   open(13,file='level'//char(48+lev)//'_'//char(48+my_rank))
      !   do icount=1,nbound_lev(lev)
      !      write(13,'(15(e28.17,1x))') xs_lev(icount,lev,2),ys_lev(icount,lev,2),zs_lev(icount,lev,2),&
      !                                 ds_lev(icount,lev,2),source_bound_lev(icount,neqs:neqf,lev,2)

      !   enddo
      !   close(13)
      !enddo

      deallocate (xs_tmp, ds_tmp, s_tmp)
   end subroutine build_level_nbound_3d

end submodule pinfdomain