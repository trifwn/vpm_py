
Subroutine rhsbcast(RHS_pm, NN, neq)
   use MPI
   use mpi_matrices, only: mpimat4
   Implicit None

   integer, intent(in)                                                  :: NN(3), neq
   ! double precision, dimension(:,:,:,:), intent(inout)                  :: RHS_pm
   double precision, dimension(neq, NN(1), NN(2), NN(3)), intent(inout)    :: RHS_pm
   !f2py depend(neq, NN) :: RHS_pm(neq, NN(1), NN(2), NN(3))
   integer                                                              :: my_rank, np, ierr, mat4

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   !---------------------------------------------
   call mpimat4(mat4, neq, NN(1), NN(2), NN(3))
   call MPI_BCAST(RHS_pm, 1, mat4, 0, MPI_COMM_WORLD, ierr)
   call MPI_TYPE_FREE(mat4, ierr)
   !-----------------------------------
End Subroutine rhsbcast

Subroutine rhsscat(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, RHS_pm_bl)
   use vpm_vars, only: neqpm
   use pmgrid, only: RHS_pm
   use MPI
   Implicit None
   integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), nb_i, nb_j, NN_bl(6), NN_tmp(3)
   ! double precision, dimension(:,:,:,:)  :: RHS_pm_bl
   double precision, intent(out) ::RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   !f2py depend(neqpm, NN_tmp) :: RHS_pm(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   integer :: my_rank, ierr
   integer :: ixs, jxs, ixf, jxf, nb, NXs, NXf, NYs, NYf, NN(3)

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

   RHS_pm_bl = 0
   nb = my_rank + 1
   NN(1:3) = NNbl(1:3, nb)
   NXs = NNbl_bl(1, nb)
   NXf = NNbl_bl(4, nb)

   NYs = NNbl_bl(2, nb)
   NYf = NNbl_bl(5, nb)
   ixs = (nb_i - 1)*(NXf - NXs) + NN_bl(1)
   jxs = (nb_j - 1)*(NYf - NYs) + NN_bl(2)

   ixf = ixs + (NXf - NXs + 1) - 1
   jxf = jxs + (NYf - NYs + 1) - 1

   RHS_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, 1) = RHS_pm(1:neqpm, ixs:ixf, jxs:jxf, 1)
   if (nb_i .gt. 1) RHS_pm_bl(:, NXs, :, :) = 0.d0
   if (nb_j .gt. 1) RHS_pm_bl(:, :, NYs, :) = 0.d0

End Subroutine rhsscat

Subroutine solget(BLOCKS, NBI, NBJ, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl)
   ! use pmgrid
   use vpm_vars, only: neqpm
   use pmeshpar, only: SOL_pm
   use mpi_matrices, only: mpimat4
   use MPI
   Implicit None
   integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), NBI, NBJ, NN_bl(6), NN_tmp(3)
   ! double precision, dimension(:,:,:,:), intent(in)  :: SOL_pm_bl
   double precision, intent(in)  :: SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   !f2py depend(neqpm, NN_tmp) :: SOL_pm(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   double precision, allocatable :: SOL_pm_tmp(:, :, :, :)
   integer :: my_rank, ierr, source, dest, status(MPI_STATUS_SIZE), mat4
   integer :: ixs, jxs, ixf, jxf, nb, NXs, NXf, NYs, NYf, j, NN_block(3), i, nbs

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

   nb = my_rank + 1
   !NN(1:3)=NN_tmp(1:3)
   !do nbs=1,BLOCKS
   !  write(*,*) my_rank,NN_tmp
   !enddo
   if (my_rank .eq. 0) then
      j = 1
      i = 1
      NXs = NNbl_bl(1, 1)
      NXf = NNbl_bl(4, 1)

      NYs = NNbl_bl(2, 1)
      NYf = NNbl_bl(5, 1)
      ixs = (i - 1)*(NXf - NXs) + NN_bl(1)
      jxs = (j - 1)*(NYf - NYs) + NN_bl(2)

      ixf = ixs + (NXf - NXs + 1) - 1
      jxf = jxs + (NYf - NYs + 1) - 1

      SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, 1) = SOL_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, 1)
      !-->Assign
      do j = 1, NBJ
         do i = 1, NBI
            nbs = (j - 1)*NBI + i
            if (nbs .eq. 1) cycle
            NN_block(1:3) = NNbl(1:3, nbs)
            allocate (SOL_pm_tmp(neqpm, NN_block(1), NN_block(2), NN_block(3)))
            NXs = NNbl_bl(1, nbs)
            NXf = NNbl_bl(4, nbs)

            NYs = NNbl_bl(2, nbs)
            NYf = NNbl_bl(5, nbs)
            ixs = (i - 1)*(NXf - NXs) + NN_bl(1)
            jxs = (j - 1)*(NYf - NYs) + NN_bl(2)

            ixf = ixs + (NXf - NXs + 1) - 1
            jxf = jxs + (NYf - NYs + 1) - 1
            call mpimat4(mat4, neqpm, NN_block(1), NN_block(2), NN_block(3))
            source = nbs - 1
            call MPI_RECV(SOL_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
            SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, 1) = SOL_pm_tmp(1:neqpm, NXs:NXf, NYs:NYf, 1)
            call MPI_TYPE_FREE(mat4, ierr)
            deallocate (SOL_pm_tmp)
         end do
      end do
   else
      dest = 0
      call mpimat4(mat4, neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      call MPI_SEND(SOL_pm_bl, 1, mat4, dest, 1, MPI_COMM_WORLD, ierr)
      call MPI_TYPE_FREE(mat4, ierr)
   end if

End Subroutine solget

Subroutine rhsscat_3d(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, nb_k, RHS_pm_bl)
   ! use pmeshpar
   use vpm_vars, only: neqpm
   use pmgrid, only: RHS_pm
   use MPI
   Implicit None
   integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), nb_i, nb_j, nb_k, NN_bl(6), NN_tmp(3)
   ! double precision, dimension(:,:,:,:), intent(out)  :: RHS_pm_bl
   double precision, intent(out) ::RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   !f2py depend(neqpm, NN_tmp) :: RHS_pm(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   integer :: my_rank, ierr
   integer :: ixs, jxs, kxs, ixf, jxf, kxf, nb, NXs, NXf, NYs, NYf, NZs, NZf, NN(3)

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

   RHS_pm_bl = 0
   nb = my_rank + 1
   NN(1:3) = NNbl(1:3, nb)

   NXs = NNbl_bl(1, nb) ! DUMMY START
   NXf = NNbl_bl(4, nb) ! DUMMY FINISH

   NYs = NNbl_bl(2, nb)
   NYf = NNbl_bl(5, nb)

   NZs = NNbl_bl(3, nb)
   NZf = NNbl_bl(6, nb)
   ixs = (nb_i - 1)*(NXf - NXs) + NN_bl(1)
   jxs = (nb_j - 1)*(NYf - NYs) + NN_bl(2)
   kxs = (nb_k - 1)*(NZf - NZs) + NN_bl(3)

   ixf = ixs + (NXf - NXs + 1) - 1
   jxf = jxs + (NYf - NYs + 1) - 1
   kxf = kxs + (NZf - NZs + 1) - 1

   RHS_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf) = RHS_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf)
   if (nb_i .gt. 1) RHS_pm_bl(:, NXs, :, :) = 0.d0
   if (nb_j .gt. 1) RHS_pm_bl(:, :, NYs, :) = 0.d0
   if (nb_k .gt. 1) RHS_pm_bl(:, :, :, NZs) = 0.d0

End Subroutine rhsscat_3d

Subroutine solget_3d(BLOCKS, NBI, NBJ, NBK, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl)
   ! use pmeshpar
   use vpm_vars, only: neqpm
   use pmeshpar, only: SOL_pm
   use mpi_matrices, only: mpimat4
   use MPI
   Implicit None
   integer, intent(in)           :: BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), NBI, NBJ, NBK, NN_bl(6), NN_tmp(3)
   ! double precision, dimension(:,:,:,:), intent(in)  :: SOL_pm_bl
   double precision, intent(in)  :: SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   ! f2py depend(neqpm, NN_tmp) :: SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
   double precision, allocatable :: SOL_pm_tmp(:, :, :, :)
   integer :: my_rank, ierr, source, dest, status(MPI_STATUS_SIZE), mat4
   integer :: ixs, jxs, kxs, ixf, jxf, kxf, nb, NXs, NXf, NYs, NYf, NZs, NZf, j, k, NN_block(3), i, nbs

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

   if (my_rank .eq. 0) write (*, *) achar(9), achar(27)//'[1;34m', 'Assembling Solution', achar(27)//'[0m'

   nb = my_rank + 1
   !NN(1:3)=NN_tmp(1:3)
   !do nbs=1,BLOCKS
   !  write(*,*) my_rank,NN_tmp
   !enddo
   if (my_rank .eq. 0) then
      j = 1
      i = 1
      k = 1
      NXs = NNbl_bl(1, nb)
      NXf = NNbl_bl(4, nb)

      NYs = NNbl_bl(2, nb)
      NYf = NNbl_bl(5, nb)

      NZs = NNbl_bl(3, nb)
      NZf = NNbl_bl(6, nb)

      ixs = (i - 1)*(NXf - NXs) + NN_bl(1)
      jxs = (j - 1)*(NYf - NYs) + NN_bl(2)
      kxs = (k - 1)*(NZf - NZs) + NN_bl(3)

      ixf = ixs + (NXf - NXs + 1) - 1
      jxf = jxs + (NYf - NYs + 1) - 1
      kxf = kxs + (NZf - NZs + 1) - 1

      SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf) = SOL_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf)
      !write(*,*) maxval(SOL_pm_bl),minval(SOL_pm_bl)
      !-->Assign
      do k = 1, NBK
         do j = 1, NBJ
            do i = 1, NBI
               nbs = (k - 1)*NBI*NBJ + (j - 1)*NBI + i
               if (nbs .eq. 1) cycle
               NN_block(1:3) = NNbl(1:3, nbs)
               allocate (SOL_pm_tmp(neqpm, NN_block(1), NN_block(2), NN_block(3)))
               NXs = NNbl_bl(1, nbs)
               NXf = NNbl_bl(4, nbs)

               NYs = NNbl_bl(2, nbs)
               NYf = NNbl_bl(5, nbs)

               NZs = NNbl_bl(3, nbs)
               NZf = NNbl_bl(6, nbs)

               ixs = (i - 1)*(NXf - NXs) + NN_bl(1)
               jxs = (j - 1)*(NYf - NYs) + NN_bl(2)
               kxs = (k - 1)*(NZf - NZs) + NN_bl(3)

               ixs = (i - 1)*(NXf - NXs) + NN_bl(1)
               jxs = (j - 1)*(NYf - NYs) + NN_bl(2)
               kxs = (k - 1)*(NZf - NZs) + NN_bl(3)

               ixf = ixs + (NXf - NXs + 1) - 1
               jxf = jxs + (NYf - NYs + 1) - 1
               kxf = kxs + (NZf - NZs + 1) - 1

               call mpimat4(mat4, neqpm, NN_block(1), NN_block(2), NN_block(3))
               source = nbs - 1
               call MPI_RECV(SOL_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
               ! write (*, *) achar(9), achar(9), 'RECEIVED TEMPORARY SOLUTION OF SOL_PM_TMP from ', source
               ! write (*, *) achar(9), achar(9), 'max = ', maxval(SOL_pm_tmp), 'min = ', minval(SOL_pm_tmp)
               SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf) = SOL_pm_tmp(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf)
               call MPI_TYPE_FREE(mat4, ierr)
               deallocate (SOL_pm_tmp)
            end do
         end do
      end do
   else
      dest = 0
      call mpimat4(mat4, neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      call MPI_SEND(SOL_pm_bl, 1, mat4, dest, 1, MPI_COMM_WORLD, ierr)
      call MPI_TYPE_FREE(mat4, ierr)
   end if

End Subroutine solget_3d

Subroutine velbcast_3d
   use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm, NXpm, NYpm, NZpm
   use mpi_matrices, only: mpimat3_pm
   use MPI
   Implicit None
   integer :: my_rank, np, ierr, mat3

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   !---------------------------------------------
   call mpimat3_pm(mat3, NXpm, NYpm, NZpm)
   call MPI_BCAST(velvrx_pm, 1, mat3, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(velvry_pm, 1, mat3, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(velvrz_pm, 1, mat3, 0, MPI_COMM_WORLD, ierr)
   call MPI_TYPE_FREE(mat3, ierr)
   !--------------------------------------------
End Subroutine velbcast_3d

Subroutine particles_scat
   use vpm_vars, only: NVR_p, neqpm, XP_scatt, QP_scatt, NVR_size
   use parvar, only: XP, QP, NVR
   use mpi_matrices, only: mpimat2_pm
   use MPI

   Implicit None
   integer :: my_rank, np, ierr, i
   integer :: dest, NVR_pr, NVR_r, mat2
   integer :: status(MPI_STATUS_SIZE)
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   !---------------------------------------------
   if (my_rank .eq. 0) then
      XP_scatt(1:3, 1:NVR_p) = XP(1:3, 1:NVR_p)
      QP_scatt(1:neqpm + 1, 1:NVR_p) = QP(1:neqpm + 1, 1:NVR_p)
      NVR_pr = NVR_p
      NVR_r = NVR/np
      if (NVR_r .gt. 0) then
         do i = 2, np
            dest = i - 1
            call mpimat2_pm(mat2, 3, NVR_size, 3, NVR_r, NVR_pr)!NVR_pr because counting starts from 0
            call MPI_SEND(XP, 1, mat2, dest, 1, MPI_COMM_WORLD, ierr)
            call MPI_TYPE_FREE(mat2, ierr)

            call mpimat2_pm(mat2, neqpm + 1, NVR_size, neqpm + 1, NVR_r, NVR_pr)
            call MPI_SEND(QP, 1, mat2, dest, 1, MPI_COMM_WORLD, ierr)
            call MPI_TYPE_FREE(mat2, ierr)
            NVR_pr = NVR_pr + NVR_r
         end do
      end if
   else
      if (NVR_p .gt. 0) then
         call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
         call MPI_RECV(XP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_TYPE_FREE(mat2, ierr)

         call mpimat2_pm(mat2, neqpm + 1, NVR_p, neqpm + 1, NVR_p, 0)
         call MPI_RECV(QP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_TYPE_FREE(mat2, ierr)
      end if
   end if
End Subroutine particles_scat

Subroutine particles_gath
   use vpm_vars, only: NVR_p, neqpm, XP_scatt, QP_scatt, UP_scatt, GP_scatt
   use parvar, only: XP, QP, UP, GP, NVR
   use mpi_matrices, only: mpimat2_pm
   use MPI

   Implicit None
   integer :: my_rank, np, ierr, i
   integer :: dest, NVR_pr, NVR_r, mat2
   integer :: status(MPI_STATUS_SIZE)
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   !---------------------------------------------
   if (my_rank .eq. 0) then
      XP(1:3, 1:NVR_p) = XP_scatt(1:3, 1:NVR_p)
      QP(1:neqpm + 1, 1:NVR_p) = QP_scatt(1:neqpm + 1, 1:NVR_p)
      UP(1:3, 1:NVR_p) = UP_scatt(1:3, 1:NVR_p)
      GP(1:3, 1:NVR_p) = GP_scatt(1:3, 1:NVR_p)
      deallocate (QP_scatt, XP_scatt, UP_scatt, GP_scatt)
      NVR_pr = NVR_p
      NVR_r = NVR/np
      allocate (XP_scatt(3, NVR_r), QP_scatt(neqpm + 1, NVR_r), UP_scatt(3, NVR_r), GP_scatt(3, NVR_r))
      if (NVR_r .gt. 0) then
         do i = 2, np
            dest = i - 1
            call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
            call MPI_RECV(XP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_FREE(mat2, ierr)

            call mpimat2_pm(mat2, neqpm + 1, NVR_r, neqpm + 1, NVR_r, 0)
            call MPI_RECV(QP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_FREE(mat2, ierr)

            call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
            call MPI_RECV(UP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_FREE(mat2, ierr)

            call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
            call MPI_RECV(GP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_FREE(mat2, ierr)

            XP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = XP_scatt(1:3, 1:NVR_r)
            QP(1:neqpm + 1, NVR_pr + 1:NVR_pr + NVR_r) = QP_scatt(1:neqpm + 1, 1:NVR_r)
            UP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = UP_scatt(1:3, 1:NVR_r)
            GP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = GP_scatt(1:3, 1:NVR_r)

            NVR_pr = NVR_pr + NVR_r
         end do
      end if
   else
      if (NVR_p .gt. 0) then
         call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
         call MPI_SEND(XP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_FREE(mat2, ierr)

         call mpimat2_pm(mat2, neqpm + 1, NVR_p, neqpm + 1, NVR_p, 0)
         call MPI_SEND(QP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_FREE(mat2, ierr)

         call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
         call MPI_SEND(UP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_FREE(mat2, ierr)

         call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
         call MPI_SEND(GP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_FREE(mat2, ierr)
      end if
   end if

   !write(filname,'(i1)') my_rank+1
   !open(15,file=filname)
   !write(15,*) 'VARIABLES="X" "Y" "Z"'
   !do i = 1,NVR_p
   !   write(15,'(7(e28.17,1x))') XP_scatt(1:3,i)!,QP_scatt(1:neqpm+1,i)
   !enddo
End Subroutine particles_gath

Subroutine proj_gath(NN)
   use vpm_vars, only: neqpm
   use pmgrid, only: RHS_pm
   ! use pmeshpar
   ! use parvar
   use mpi_matrices, only: mpimat4
   use MPI
   Implicit None
   integer, intent(in) :: NN(3)
   integer :: my_rank, np, ierr, i
   integer :: dest, source, mat4
   integer :: status(MPI_STATUS_SIZE)
   double precision, allocatable:: RHS_pm_tmp(:, :, :, :)
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   if (my_rank .eq. 0) then
      allocate (RHS_pm_tmp(neqpm + 1, NN(1), NN(2), NN(3)))
      do i = 2, np
         source = i - 1
         call mpimat4(mat4, neqpm + 1, NN(1), NN(2), NN(3))
         call MPI_RECV(RHS_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_TYPE_FREE(mat4, ierr)
         RHS_pm = RHS_pm + RHS_pm_tmp
      end do
      deallocate (RHS_pm_tmp)
   else
      dest = 0
      call mpimat4(mat4, neqpm + 1, NN(1), NN(2), NN(3))
      call MPI_SEND(RHS_pm, 1, mat4, dest, 1, MPI_COMM_WORLD, ierr)
      call MPI_TYPE_FREE(mat4, ierr)
   end if
End Subroutine proj_gath

Subroutine proj_gath_new(NN)
   use vpm_vars, only: interf_iproj, neqpm, XP_scatt, neqpm
   use pmgrid, only: RHS_pm, XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm
   use pmeshpar, only: ND
   ! use parvar
   use mpi_matrices, only: mpimat4
   use MPI

   Implicit None
   integer, intent(in) :: NN(3)
   integer :: my_rank, np, ierr, NN_proj(6), nn1, nn2, nn3, NN_tmp(6)
   integer :: source, mat4
   integer :: status(MPI_STATUS_SIZE)
   integer :: imax, imin, jmax, jmin, kmax, kmin
   double precision, allocatable:: RHS_pmtmp(:, :, :, :)
   double precision            :: xpmax, xpmin, ypmax, ypmin, zpmax, zpmin

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   if (ND .eq. 2) then
      xpmax = maxval(XP_scatt(1, :)); xpmin = minval(XP_scatt(1, :))
      ypmax = maxval(XP_scatt(2, :)); ypmin = minval(XP_scatt(2, :))
      imax = int((xpmax - XMIN_pm)/DXpm) + 1; imin = int((xpmin - XMIN_pm)/DXpm) + 1
      jmax = int((ypmax - YMIN_pm)/DYpm) + 1; jmin = int((ypmin - YMIN_pm)/DYpm) + 1
      NN_proj(1) = max(imin - interf_iproj, 1); NN_proj(2) = max(jmin - interf_iproj, 1); NN_proj(3) = 1       !call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
      NN_proj(4) = min(imax + interf_iproj, NN(1)); NN_proj(5) = min(jmax + interf_iproj, NN(2)); NN_proj(6) = 1       !call  projlib
   else
      xpmax = maxval(XP_scatt(1, :)); xpmin = minval(XP_scatt(1, :))
      ypmax = maxval(XP_scatt(2, :)); ypmin = minval(XP_scatt(2, :))
      zpmax = maxval(XP_scatt(3, :)); zpmin = minval(XP_scatt(3, :))
      imax = int((xpmax - XMIN_pm)/DXpm) + 1; imin = int((xpmin - XMIN_pm)/DXpm) + 1
      jmax = int((ypmax - YMIN_pm)/DYpm) + 1; jmin = int((ypmin - YMIN_pm)/DYpm) + 1
      kmax = int((zpmax - ZMIN_pm)/DZpm) + 1; kmin = int((zpmin - ZMIN_pm)/DZpm) + 1
      NN_proj(1) = max(imin - interf_iproj, 1); NN_proj(2) = max(jmin - interf_iproj, 1); NN_proj(3) = max(kmin - interf_iproj, 1)     !call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
      NN_proj(4) = min(imax + interf_iproj, NN(1)); NN_proj(5) = min(jmax + interf_iproj, NN(2)); NN_proj(6) = min(kmax + interf_iproj, NN(3))
   end if
   if (my_rank .ne. 0) then
      call MPI_SEND(NN_proj, 6, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
      nn1 = NN_proj(4) - NN_proj(1) + 1
      nn2 = NN_proj(5) - NN_proj(2) + 1
      nn3 = NN_proj(6) - NN_proj(3) + 1
      allocate (RHS_pmtmp(neqpm + 1, nn1, nn2, nn3))
      RHS_pmtmp(1:neqpm + 1, 1:nn1, 1:nn2, 1:nn3) = &
         RHS_pm(1:neqpm + 1, NN_proj(1):NN_proj(4), NN_proj(2):NN_proj(5), NN_proj(3):NN_proj(6))
      call mpimat4(mat4, neqpm + 1, nn1, nn2, nn3)
      call MPI_SEND(RHS_pmtmp, 1, mat4, 0, 1, MPI_COMM_WORLD, ierr)
      call MPI_TYPE_FREE(mat4, ierr)
      deallocate (RHS_pmtmp)
   else
      do source = 1, np - 1
         call MPI_RECV(NN_tmp, 6, MPI_INTEGER, source, 1, MPI_COMM_WORLD, status, ierr)
         nn1 = NN_tmp(4) - NN_tmp(1) + 1
         nn2 = NN_tmp(5) - NN_tmp(2) + 1
         nn3 = NN_tmp(6) - NN_tmp(3) + 1
         allocate (RHS_pmtmp(neqpm + 1, nn1, nn2, nn3))
         call mpimat4(mat4, neqpm + 1, nn1, nn2, nn3)
         call MPI_RECV(RHS_pmtmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
         RHS_pm(1:neqpm + 1, NN_tmp(1):NN_tmp(4), NN_tmp(2):NN_tmp(5), NN_tmp(3):NN_tmp(6)) = &
            RHS_pm(1:neqpm + 1, NN_tmp(1):NN_tmp(4), NN_tmp(2):NN_tmp(5), NN_tmp(3):NN_tmp(6)) + &
            RHS_pmtmp(1:neqpm + 1, 1:nn1, 1:nn2, 1:nn3)
         deallocate (RHS_pmtmp)
         call MPI_TYPE_FREE(mat4, ierr)
      end do
   end if

End Subroutine proj_gath_new

