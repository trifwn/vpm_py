!Subroutine remesh_particles
!This Subroutine remeshes particles  on the pm grid (4 particles per cell)
!---------------------------------------------------------------------------------------------------

Subroutine remesh_particles_3d(iflag, XP_out, QP_out, NVR_out)
   use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, &
                     DZpm, YMAX_pm, XMAX_pm, ZMAX_pm, RHS_pm, NXpm, &
                     NYpm, NZpm, NXs_bl, NYs_bl, NZs_bl, NXf_bl, &
                     NYf_bl, NZf_bl, DVpm, EPSVOL
   use projlib, only: projlibinit, project_particles_3D, project_vol3d
   use vpm_vars, only: mrem, NEQPM, NVR_p, XP_scatt, QP_scatt, NVR_projscatt, INTERF_IPROJ, V_REF, NCELL_REM, NVR_SIZE
   use pmeshpar, only: NPAR_CELL, IDVPM, ND

   use parvar, only: NVR, XP, QP, GP, UP
   ! use openmpth, only: OMPTHREADS
   use MPI

   Implicit None

   ! PARAMETERS
   integer, intent(in)  :: iflag
   double precision, intent(out),allocatable, target :: XP_out(:, :), QP_out(:, :)
   integer, intent(out) :: NVR_out

   ! LOCAL VARIABLES
   double precision, dimension(8)   :: X, Y, Z
   double precision    :: XMIN_vr, YMIN_vr, DXvr, DYvr, DZvr, ZMIN_vr
   double precision, allocatable:: XC(:), YC(:), ZC(:)
   double precision, allocatable, target::XP_tmp(:, :), QP_tmp(:, :)
   integer             :: i, j, k, NXpm1, NYpm1, NZpm1, ncell, npar, ndumc
   integer             :: nxstart, nxfin, nystart, nyfin, nzstart, nzfin, ndum_rem, nnod, nc
   integer, allocatable :: ieq(:)
   double precision    :: Xbound(6), Dpm(3), wmag
   double precision, allocatable :: QINF(:)
   integer    :: NVR_OLD 
   integer             :: my_rank, ierr, np, NN(3), NN_bl(6)

   ! DEPRECATED
   ! double precision,intent(inout):: XP_in(:,:),QP_in(:,:)
   ! double precision :: Vol, ANG, dens1, dens2, Mach,
   ! integer   :: nv, inode, jnode, knode, itype
   ! integer   :: iis, jjs, kks, iif, jjf, kkf
   ! integer   :: iis2, jjs2, kks2, iif2, jjf2, kkf2
   ! double precision    :: fx, fy, f
   ! double precision    :: w1, w2, r1, r2, core, radi, th, xx, yy
   ! double precision, allocatable :: rhsper(:, :, :, :)

   NVR_OLD = NVR

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   XMIN_vr = XMIN_pm
   YMIN_vr = YMIN_pm
   ZMIN_vr = ZMIN_pm
   DXvr = DXpm
   DYvr = DYpm
   DZvr = DZpm
   Ndumc = 1
   npar_cell = 1

   Dpm(1) = DXpm; Dpm(2) = DYpm; Dpm(3) = DZpm
   Xbound(1) = XMIN_pm;
   Xbound(2) = YMIN_pm; 
   Xbound(3) = ZMIN_pm

   Xbound(4) = XMAX_pm;
   Xbound(5) = YMAX_pm; 
   Xbound(6) = ZMAX_pm
   !->PM grid is orthogonal (Volume of PM cell

   !--------------------------------------------------------------------------------!
   !-->The loops starts from 2 because we need cells that DO NOT contain particles  !
   !-->Total Number of Cells                                                        !
   !--------------------------------------------------------------------------------!
   NN(1) = NXpm; NN(2) = NYpm; NN(3) = NZpm
   NN_bl(1) = NXs_bl(1); NN_bl(2) = NYs_bl(1); NN_bl(3) = NZs_bl(1)
   NN_bl(4) = NXf_bl(1); NN_bl(5) = NYf_bl(1); NN_bl(6) = NZf_bl(1)

   NN = NN*mrem
   NN_bl = NN_bl*mrem
   Dpm(1) = (Xbound(4) - Xbound(1))/(NN(1) - 1)
   Dpm(2) = (Xbound(5) - Xbound(2))/(NN(2) - 1)
   Dpm(3) = (Xbound(6) - Xbound(3))/(NN(3) - 1)

   DVpm = Dpm(1)*Dpm(2)*Dpm(3)
   ! Here we project the particles on the PM grid
   if (iflag .eq. 1) then
      if (allocated(RHS_pm)) then
         deallocate (RHS_pm)
         allocate (RHS_pm(neqpm + 1, NN(1), NN(2), NN(3)))
      else
         allocate (RHS_pm(neqpm + 1, NN(1), NN(2), NN(3)))
      end if
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      NVR_p = NVR/np
      if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)
      allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm + 1, NVR_p), NVR_projscatt(NVR_p))
      NVR_projscatt = interf_iproj
      call particles_scat
      ! PARTICLES ARE SCATTERED ON PROCESSORS

      call projlibinit(Xbound, Dpm, NN, NN_bl, EPSVOL, IDVPM, ND)
      allocate (ieq(neqpm + 1), QINF(neqpm + 1))
      QINF = 0.d0
      do i = 1, neqpm + 1
         ieq(i) = i
      end do

      ! PROJECT PARTICLES ON PM GRID
      call project_particles_3D(RHS_pm, QP_scatt, XP_scatt, NVR_projscatt, NVR_p, neqpm + 1, &
                                ieq, neqpm + 1, QINF, NVR_p)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call proj_gath(NN)
      ! RHS IS GATHERED BACK ON PROCESSOR 0
      deallocate (XP_scatt, QP_scatt, NVR_projscatt)

      if (my_rank .eq. 0) then
         RHS_pm(neqpm + 1, :, :, :) = DVpm
         call project_vol3d(RHS_pm, neqpm + 1, ieq, neqpm + 1, IDVPM)
         !call hill_assign(NN,NN_bl,Xbound,Dpm,RHS_pm,neqpm+1)
      end if
      deallocate (ieq, QINF)
   end if

   ! HERE WE REMESH PARTICLES
   if (my_rank .eq. 0) then
      ncell = ncell_rem
      ndum_rem = 2
      nnod = 1 !if ncell gt 1 particles IN cell else in nodes
      if (ncell .eq. 1) nnod = 0
      allocate (XC(ncell), YC(ncell), ZC(ncell))
      if (ncell .eq. 1) then
         nxfin = NN_bl(4) - interf_iproj/2
         nxstart = NN_bl(1) + interf_iproj/2
         nyfin = NN_bl(5) - interf_iproj/2
         nystart = NN_bl(2) + interf_iproj/2
         nzfin = NN_bl(6) - interf_iproj/2
         nzstart = NN_bl(3) + interf_iproj/2
      else
         nxfin = NN_bl(4) - 1
         nxstart = NN_bl(1)
         nyfin = NN_bl(5) - 1
         nystart = NN_bl(2)
         nzfin = NN_bl(6) - 1
         nzstart = NN_bl(3)
      end if

      NXpm1 = nxfin - nxstart + 1
      NYpm1 = nyfin - nystart + 1
      NZpm1 = nzfin - nzstart + 1
      NVR = NXpm1*NYpm1*NZpm1*ncell
      
      print *, achar(9), achar(9), 'Allocating XP_tmp and QP_tmp with NVR:', NVR
      allocate (XP_tmp(3, NVR), QP_tmp(neqpm + 1, NVR))
      XP => XP_tmp
      QP => QP_tmp
      XP = 0
      QP = 0
      npar = 0
      V_ref = 1.d0/float(ncell)*DVpm
      ! !$omp parallel private(i,j,k,npar,X,Y,Z) num_threads(OMPTHREADS)
      ! !$omp do
      do k = nzstart, nzfin
         do j = nystart, nyfin
            do i = nxstart, nxfin
               ! !-> Get PM cell nodes (orthogonal structured grid
               X(1) = XMIN_pm + Dpm(1)*(float(i) - 1.d0)
               X(2) = XMIN_pm + Dpm(1)*(float(i))
               X(3) = XMIN_pm + Dpm(1)*(float(i))
               X(4) = XMIN_pm + Dpm(1)*(float(i) - 1.d0)
               X(5) = XMIN_pm + Dpm(1)*(float(i) - 1.d0)
               X(6) = XMIN_pm + Dpm(1)*(float(i))
               X(7) = XMIN_pm + Dpm(1)*(float(i))
               X(8) = XMIN_pm + Dpm(1)*(float(i) - 1.d0)

               Y(1) = YMIN_pm + Dpm(2)*(j - 1)
               Y(2) = YMIN_pm + Dpm(2)*(j - 1)
               Y(3) = YMIN_pm + Dpm(2)*(j)
               Y(4) = YMIN_pm + Dpm(2)*(j)
               Y(5) = YMIN_pm + Dpm(2)*(j - 1)
               Y(6) = YMIN_pm + Dpm(2)*(j - 1)
               Y(7) = YMIN_pm + Dpm(2)*(j)
               Y(8) = YMIN_pm + Dpm(2)*(j)

               Z(1) = ZMIN_pm + Dpm(3)*(k - 1)
               Z(2) = ZMIN_pm + Dpm(3)*(k - 1)
               Z(3) = ZMIN_pm + Dpm(3)*(k - 1)
               Z(4) = ZMIN_pm + Dpm(3)*(k - 1)
               Z(5) = ZMIN_pm + Dpm(3)*(k)
               Z(6) = ZMIN_pm + Dpm(3)*(k)
               Z(7) = ZMIN_pm + Dpm(3)*(k)
               Z(8) = ZMIN_pm + Dpm(3)*(k)

               !npar = ((k-nzstart)*NXpm1*NYpm1 +(j-nystart)*NXpm1 + i-nxstart)*ncell

               if (ncell .gt. 1) then
                  YC = cell3d_interp_euler(Y, ncell, 2)
                  ZC = cell3d_interp_euler(Z, ncell, 2)
                  XC = cell3d_interp_euler(X, ncell, 2)
                  do nc = 1, ncell
                     npar = npar + 1
                     XP(1, npar) = XC(nc)
                     XP(2, npar) = YC(nc)
                     XP(3, npar) = ZC(nc)
                     QP(neqpm + 1, npar) = 1.d0/float(ncell)*DVpm
                  end do

               else
                  wmag = sqrt(RHS_pm(1, i, j, k)**2 + RHS_pm(2, i, j, k)**2 + RHS_pm(3, i, j, k)**2)
                  if (wmag .lt. 1e-09) cycle
                  npar = npar + 1
                  XP(1, npar) = X(1)
                  XP(2, npar) = Y(1)
                  XP(3, npar) = Z(1)

                  QP(1:neqpm, npar) = RHS_pm(1:neqpm, i, j, k)*DVpm
                  QP(neqpm + 1, npar) = DVpm
               end if
            end do
         end do
      end do
      !!$omp enddo
      !!$omp endparallel
      NVR = npar
      NVR_size = NVR
      NVR_out = NVR

      allocate (XP_out(3,npar), QP_out(neqpm + 1, npar))
      XP_out = XP(:, 1:npar)
      QP_out = QP(:, 1:npar)
      nullify (XP, QP)
      XP => XP_out
      QP => QP_out
      if (ncell .gt. 1) call back_to_particles_3D_rem(RHS_pm, XP, QP, Xbound, Dpm, NN, NVR, 4)
      if (iflag .eq. 0) deallocate (RHS_pm)

      write (*, *) achar(9), 'After remesh'
      write (*, *) achar(9), achar(9), 'Number of particles before', NVR_OLD
      write (*, *) achar(9), achar(9), 'Number of particles after', NVR
      write (*, *) achar(9), achar(9), 'Volume of a cell', DVpm
      write (*, *) achar(9), achar(9), 'Number of cells', NXpm, NYpm, NZpm
      write (*, *) achar(9), achar(9), 'Size of XP', size(XP,1) , size(XP,2)
      write (*, *) achar(9), achar(9), 'Size of QP', size(QP,1) , size(QP,2)
      write (*, *) achar(9), achar(9), 'Maximal value of QPR', maxval(QP(neqpm, :))
   end if
   !call back_to_particles_2D(4)

   ! open(1,file='vr.dat')
   !  WRITE(1,*)'VARIABLES = "X" "Y" "Z"'
   !  do  i=1, NVR
   !       write(1,'(2(F20.10,1x))') XP(i,1), XP(i,2),XP(i,3)
   !  enddo
    !!-----FOR PLOTTING PURPOSES ONLY
   !     call system('preplot vr.dat>/dev/null')
   !     call system('rm vr.dat')
    !!----FOR PLOTTING PURPOSES ONLY
   ! close(1)
   if (allocated(RHS_pm)) deallocate (RHS_pm)
End Subroutine remesh_particles_3d

!---------------------------------------------------------------------------!
!-> Subroutine back_to_particles                                            !
!   This subroutine interpolates PM grid values back to particles at the    !
!   positions they ARE.Convections takes place afterwards.                  !
!   Input :                                                                 !
!          itype (1,2) defines what value to interpolate to the particles   !
!---------------------------------------------------------------------------!
Subroutine back_to_particles_3D_rem(RHS_pm, XP, QP, Xbound, Dpm, NN, NVR, iproj)
   use openmpth
   use projlib, only: projection_fun
   Implicit None
   integer, intent(in)             :: NN(3), NVR, iproj
   ! double precision, intent(in), dimension(:, :, :, :) :: RHS_pm
   double precision, intent(in)    :: RHS_pm(4, NN(1), NN(2), NN(3))
   ! f2py depend(NN) :: RHS_pm(4, NN(1), NN(2), NN(3))
   double precision, intent(inout) :: QP(4, NVR), XP(3, NVR)

   double precision, intent(in)    :: Xbound(6), Dpm(3)

   double precision :: fx, fy, fz, f, x, y, z
   integer          :: inode, jnode, knode, i, j, k, nv, ips, ipf

   if (iproj .eq. 2) then
      ips = 0
      ipf = 1
   else if (iproj .eq. 3) then
      ips = 1
      ipf = 2
   else if (iproj .eq. 4) then
      ips = 1
      ipf = 2
   end if
   QP(:, 1:3) = 0
   do nv = 1, NVR
      !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
      inode = int((XP(1, nv) - XBound(1))/Dpm(1)) + 1
      jnode = int((XP(2, nv) - XBound(2))/Dpm(2)) + 1
      knode = int((XP(3, nv) - XBound(3))/Dpm(3)) + 1

      !--We search the 4 nodes close to the particles
      do k = knode - ips, knode + ipf
         do j = jnode - ips, jnode + ipf
            do i = inode - ips, inode + ipf
               x = (XP(1, nv) - XBound(1) - (i - 1)*Dpm(1))/Dpm(1)
               fx = projection_fun(iproj, x)

               y = (XP(2, nv) - XBound(2) - (j - 1)*Dpm(2))/Dpm(2)
               fy = projection_fun(iproj, y)

               z = (XP(3, nv) - XBound(3) - (k - 1)*Dpm(3))/Dpm(3)
               fz = projection_fun(iproj, z)

               f = fx*fy*fz
               QP(1:3, nv) = QP(1:3, nv) + f*RHS_pm(1:3, i, j, k)

            end do
         end do
      end do
      QP(1:3, nv) = QP(1:3, nv)*QP(4, nv)
   end do

End Subroutine back_to_particles_3D_rem

!--------------------------------------------------------------------------------
!>@function
! Subroutine    cell3d_interp_euler
!>
!>@author Papis
!>
!>@brief
!>Subroutine cell3d_interp_euler creates 4 or more particles per cell using ksi ita
!>coordinates
!REVISION HISTORY
!> 17/7/2013 - Initial Version
!> TODO_dd
!>
!>@param [in]  F is the value at the global coordinates
!>@param [out] FC is the value at global coordinates of the interpolated value
!--------------------------------------------------------------------------------
function cell3d_interp_euler(F, N, M) result(FC)
   use iso_fortran_env
   implicit none

   integer, parameter :: dp = real64

   double precision, dimension(8), intent(in) :: F
   integer, intent(in) :: N, M
   double precision, dimension(N) :: FC

   ! LOCAL VARIABLES
   double precision :: KSIC(8), HTAC(8), ZETAC(8)
   double precision, dimension(N) :: KSI, HTA, ZETA
   integer :: i, j
   double precision :: addit
   integer :: NTEMP

   NTEMP = N
   !-->Define KSI,HTA corners
   KSIC = [-1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0, -1.d0]
   HTAC = [-1.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0]
   ZETAC = [-1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0, 1.d0, 1.d0]

   FC = 0.0_dp
   write (*, *) 'NTEMP', NTEMP, N, M
   call get_ksi_ita_pos_3d(N, M, KSIC, HTAC, ZETAC, KSI, HTA, ZETA)

   do i = 1, 8 !cell nodes
      do j = 1, NTEMP
         addit = F(i)
         addit = addit*(1.0_dp + KSI(j)*KSIC(i))
         addit = addit*(1.0_dp + HTA(j)*HTAC(i))
         addit = addit*(1.0_dp + ZETA(j)*ZETAC(i))
         FC(j) = FC(j) + addit
      end do
   end do

   FC = 0.125d0*FC ! 1/8
end function cell3d_interp_euler

!--------------------------------------------------------------------------------
!>@function
! Subroutine    get_ksi_ita_pos
!>
!>@author Papis
!>
!>@brief
!>Subroutine get_ksi_ita_pos  depending on the defined number of particles(must be perfect square)
!REVISION HISTORY
!> 22/7/2013 - Initial Version
!> TODO_dd
!>
!>@param [in]  N is the number of particles
!>@param [in]  KSIC(4),HTAC(4) is the corner coordinates in the KSI,HTA
!>@param [out] KSI(2*N),HTA(2*N) local position
!--------------------------------------------------------------------------------
Subroutine get_ksi_ita_pos_3d(N, M, KSIC, HTAC, ZETAC, KSI, HTA, ZETA)
   Implicit None

   integer, intent(in)           :: M
   integer, intent(in)           :: N
   ! N is the number of particles to remesh
   ! M is 2
   double precision, intent(in)  :: KSIC(8), HTAC(8), ZETAC(8)
   double precision, intent(out) :: KSI(N), HTA(N), ZETA(N)
   double precision              :: DKSI, DHTA, DZETA
   integer                       :: i, j, k, nod

   !--> find position minus ksi minus ita quadrant and then by symmerty * 4
   KSI = 0.d0
   HTA = 0.d0

   DKSI = dabs(2.d0/float(M))
   DHTA = dabs(2.d0/float(M))
   DZETA = dabs(2.d0/float(M))
   write (*, *) 'SIZE of KSI ', size(KSI)
   do k = 1, M
      do j = 1, M
         do i = 1, M
            nod = (k - 1)*M*M + (j - 1)*M + i
            write (*, *) 'nod', nod
            KSI(nod) = KSIC(1) + (i - 1./2.)*DKSI
            HTA(nod) = HTAC(1) + (j - 1./2.)*DHTA
            ZETA(nod) = ZETAC(1) + (k - 1./2.)*DZETA
         end do
      end do
   end do

End Subroutine get_ksi_ita_pos_3d
