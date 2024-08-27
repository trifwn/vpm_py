!This library solves the poisson problem using domain decomposition method
Module yapslib
   use base_types, only: dp
   use projlib
   ! use pmlib
   use mpi_matrices!, only: mpimat4, MPIMAT5_OLD, mpimat2, mpimat5, mpimat2_pm_
   use MPI

   private

   real(dp), allocatable      :: SOL_pm_coarse(:, :, :, :), RHS_pm_coarse(:, :, :, :), SOL_pm_sample(:, :, :, :, :)
   real(dp), allocatable      :: SOL_pm_sumsample(:, :, :, :)

   real(dp)             :: DXpm, DYpm, DZpm, DXpm_c, DYpm_c, DZpm_c
   integer                      :: NXpm, NYpm, NZpm, NXpm_c, NYpm_c, NZPm_c
   integer                      :: ibctyp, idvpm, epsvol, ND, iproj, ndumcell, npmsize
   real(dp)             :: XMIN_pm, YMIN_pm, ZMIN_pm, XMAX_pm, YMAX_pm, ZMAX_pm
   real(dp)             :: XMIN_pm_c, YMIN_pm_c, ZMIN_pm_c, XMAX_pm_c, YMAX_pm_c, ZMAX_pm_c
   real(dp)             :: MACH
   real(dp), allocatable :: Xbound_bl(:, :)
   integer                      :: BLOCKS, NXB, NYB, NBI, NBJ, NBK, NB, i, j, k, NXbl, NYbl, NN(3), NN_bl(6)
   integer, allocatable          :: NNbl(:, :), NNbl_bl(:, :), NN_coarse_map(:, :), map_nodes(:, :, :, :), nnb(:)
   real(dp)             :: projection_fun, fx, fy, fz, f, xc, yc, zc, X(3), addsol, starttime, endtime
   integer                      :: ic, jc, kc, inode, jnode, knode
   integer                      :: i_nb, j_nb, k_nb
   integer                      :: NXs, NYs, NZs, NXf, NYf, NZf, ib, jb, kb, ibj, jbj, kbj, ixs, ixf, jxs, jxf, izs, izf
   integer                      :: nc, NN_map(6), isubtrackt, node, neq, isizex, isizey, isizez, nbc

   real(dp), pointer                 :: SOL_pm_bl(:, :, :, :), RHS_pm_bl(:, :, :, :), QP(:, :), XP(:, :)

   real(dp), allocatable             :: BBound(:, :, :)

   integer  :: status(MPI_STATUS_SIZE), source, ierr, my_rank, np, mat4, mat5, dest
   character*25                :: outfil1, outfil2

   public  :: yaps2d, yaps3d
contains
   include 'yaps2d.f90'
   include 'yaps3d.f90'
End Module yapslib
