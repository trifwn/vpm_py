Module vpm_size
   double precision, save              :: Xbound(6), Dpm(3), Xbound0(6), Dpm0(3)
   integer, save                       :: NN_bl(6), NN(3), NN0_bl(6), NN0(3)
   integer, save                       :: NXs0_bl(10), NYs0_bl(10), NXf0_bl(10), NYf0_bl(10), NZs0_bl(10), NZf0_bl(10)

   double precision, allocatable, save :: Xbound_bl(:, :)
   integer, allocatable, save          :: NNbl_bl(:, :), NNbl(:, :)

   double precision, save              :: Xbound_tmp(6), Xbound_coarse(6), Dpm_coarse(3)
   integer, save                       :: NN_tmp(3), NN_bl_tmp(6), NN_coarse(3), NN_bl_coarse(6)
   integer, save                       :: nb_i, nb_j, nb_k, NBB, NXbl, NYbl, NZbl, BLOCKS, NXB, NYB, NZB, ndumcell_coarse, ndumcell_bl
   double precision                    :: starttime, endtime, st, et, ct
   integer, save                       :: iynbc, iret, NBI, NBJ, NBK, NVR_out_thres, NREMESH, ntorder, &
                                          iyntree, ilevmax, itree, nsize_out(3), ibctyp, NWRITE

   public :: print_vpm_size_info, get_NN_bl, get_NN, get_Xbound

contains
   Subroutine get_NN_bl(NN_bl_out)
      integer :: NN_bl_out(6)
      NN_bl_out = NN_bl
   end Subroutine get_NN_bl

   Subroutine get_NN(NN_out)
      integer :: NN_out(3)
      NN_out = NN
   end Subroutine get_NN
   
   Subroutine get_Xbound(Xbound_out)
      double precision :: Xbound_out(6)
      Xbound_out = Xbound
   end Subroutine get_Xbound

   subroutine print_vpm_size_info()
      use io, only: dp_1d_array_info, i_1d_array_info, dp_2d_alloc_info, i_2d_alloc_info, &
                     dp_1d_alloc_info
      print *, "VPM_SIZE INFO"
      print *, "============"
      print *, ""
      call dp_1d_array_info("Xbound", Xbound, 6)
      call dp_1d_array_info("Dpm", Dpm, 3)
      call dp_1d_array_info("Xbound0", Xbound0, 6)
      call dp_1d_array_info("Dpm0", Dpm0, 3)
      call i_1d_array_info("NN_bl", NN_bl, 6)
      call i_1d_array_info("NN", NN, 3)
      call i_1d_array_info("NN0_bl", NN0_bl, 6)
      call i_1d_array_info("NN0", NN0, 3)
      call i_1d_array_info("NXs0_bl", NXs0_bl, 10)
      call i_1d_array_info("NYs0_bl", NYs0_bl, 10)
      call i_1d_array_info("NXf0_bl", NXf0_bl, 10)
      call i_1d_array_info("NYf0_bl", NYf0_bl, 10)
      call i_1d_array_info("NZs0_bl", NZs0_bl, 10)
      call i_1d_array_info("NZf0_bl", NZf0_bl, 10)
      call dp_2d_alloc_info("Xbound_bl", Xbound_bl)
      call i_2d_alloc_info("NNbl_bl", NNbl_bl)
      call i_2d_alloc_info("NNbl", NNbl)
      call dp_1d_array_info("Xbound_tmp", Xbound_tmp, 6)
      call dp_1d_array_info("Xbound_coarse", Xbound_coarse, 6)
      call dp_1d_array_info("Dpm_coarse", Dpm_coarse, 3)
      call i_1d_array_info("NN_tmp", NN_tmp, 3)
      call i_1d_array_info("NN_bl_tmp", NN_bl_tmp, 6)
      call i_1d_array_info("NN_coarse", NN_coarse, 3)
      call i_1d_array_info("NN_bl_coarse", NN_bl_coarse, 6)
      print *, achar(9)//"nb_i", " = ", nb_i
      print  *, achar(9)//"nb_j", " = ", nb_j
      print  *, achar(9)//"nb_k", " = ", nb_k
      print  *, achar(9)//"NBB", " = ", NBB
      print  *, achar(9)//"NXbl", " = ", NXbl
      print  *, achar(9)//"NYbl", " = ", NYbl
      print  *, achar(9)//"NZbl", " = ", NZbl
      print  *, achar(9)//"BLOCKS", " = ", BLOCKS
      print  *, achar(9)//"NXB", " = ", NXB
      print  *, achar(9)//"NYB", " = ", NYB
      print  *, achar(9)//"NZB", " = ", NZB
      print  *, achar(9)//"ndumcell_coarse", " = ", ndumcell_coarse
      print  *, achar(9)//"ndumcell_bl", " = ", ndumcell_bl
      print  *, achar(9)//"iynbc", " = ", iynbc
      print  *, achar(9)//"iret", " = ", iret
      print  *, achar(9)//"NBI", " = ", NBI
      print  *, achar(9)//"NBJ", " = ", NBJ
      print  *, achar(9)//"NBK", " = ", NBK
      print  *, achar(9)//"NVR_out_thres", " = ", NVR_out_thres
      print  *, achar(9)//"NREMESH", " = ", NREMESH
      print  *, achar(9)//"ntorder", " = ", ntorder
      print  *, achar(9)//"iyntree", " = ", iyntree
      print  *, achar(9)//"ilevmax", " = ", ilevmax
      print  *, achar(9)//"itree", " = ", itree
      print  *, achar(9)//"nsize_out", " = ", nsize_out
      print  *, achar(9)//"ibctyp", " = ", ibctyp
      print  *, achar(9)//"NWRITE", " = ", NWRITE
   end subroutine print_vpm_size_info


End Module vpm_size
