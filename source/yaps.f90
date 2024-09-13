!This library solves the poisson problem using domain decomposition method
Module yaps
   use base_types, only: dp
   use mpi_matrices!, only: mpimat4, MPIMAT5_OLD, mpimat2, mpimat5, mpimat2_pm_
   use MPI

   private

   real(dp), pointer       :: SOL_pm_bl(:, :, :, :), RHS_pm_bl(:, :, :, :), QP(:, :), XP(:, :)
   real(dp), allocatable   :: SOL_pm_coarse(:, :, :, :), RHS_pm_coarse(:, :, :, :), SOL_pm_sample(:, :, :, :, :)
   real(dp), allocatable   :: SOL_pm_sumsample(:, :, :, :)
   real(dp), allocatable   :: BBound(:, :, :)
   integer, allocatable    :: NN_coarse_map(:, :), map_nodes(:, :, :, :), nnb(:)

   integer                 :: NB, i , j, k
   real(dp)                :: DXpm_c, DYpm_c, DZpm_c
   integer                 :: NXpm_c, NYpm_c, NZPm_c

   real(dp)                :: fx, fy, fz, f, xc, yc, zc, X(3), starttime, endtime, total_starttime
   integer                 :: ic, jc, kc, inode, jnode, knode
   
   integer                 :: NXs, NYs, NZs, NXf, NYf, NZf
   integer                 :: NN_map(6), node, neq, isizex, isizey, isizez, nbc
   integer                 :: iproj
   integer                 :: status(MPI_STATUS_SIZE), source, ierr, my_rank, np, mat4, mat5, dest


   public                  :: yaps2d, yaps3d

   interface yaps2d
      module subroutine yaps2d(DSOL_pm, DRHS_pm, Xbound_bl, Xbound_coarse, Dpm_fine, Dpm_coarse, NNbl, NNbl_bl, &
            NN_coarse, NN_bl_coarse, ND, BLOCKS, ibctyp, neqs, neqf, nc, NBI, NBJ, nb_i, nb_j, ireturn, &
            iyntree, ilevmax, npmsize)
         integer, intent(in)             :: ibctyp, neqs, neqf, nc, ireturn, iyntree, ilevmax, npmsize
         integer, intent(in)             :: ND, BLOCKS, NNbl(3, BLOCKS), NNBl_bl(6, BLOCKS)
         integer, intent(in)             :: NN_coarse(3), NN_bl_coarse(6), nb_i, nb_j, NBI, NBJ

         real(dp), intent(in)            :: Xbound_bl(6, BLOCKS), Xbound_coarse(6)
         real(dp), intent(in)            :: Dpm_fine(3), Dpm_coarse(3)
         real(dp), intent(inout), target :: DSOL_pm(:, :, :, :), DRHS_pm(:, :, :, :)
      end subroutine yaps2d
   end interface yaps2d

   interface yaps3d
      module subroutine yaps3d(DSOL_pm, DRHS_pm, Xbound_bl, Xbound_coarse, Dpm_fine, Dpm_coarse, NNbl, NNbl_bl, &
               NN_coarse, NN_bl_coarse, ND, BLOCKS, ibctyp, neqs, neqf, nc, NBI, NBJ, NBK, nb_i, &
               nb_j, nb_k, ireturn, iyntree, ilevmax, npmsize)
         integer, intent(in)              :: ibctyp, neqs, neqf, nc, ireturn, iyntree, ilevmax, npmsize
         integer, intent(in)              :: ND, BLOCKS, NNbl(3, BLOCKS), NNBl_bl(6, BLOCKS)
         integer, intent(in)              :: NN_coarse(3), NN_bl_coarse(6), nb_i, nb_j, nb_k, NBI, NBJ, NBK

         real(dp), intent(in)             :: Xbound_bl(6, BLOCKS), Xbound_coarse(6)
         real(dp), intent(in)             :: Dpm_fine(3), Dpm_coarse(3)
         real(dp), intent(inout), target  :: DSOL_pm(:, :, :, :), DRHS_pm(:, :, :, :)
      end subroutine yaps3d
   end interface yaps3d
End Module yaps
