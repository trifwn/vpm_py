Module vpm_size
   use base_types, only: dp      
   use io, only: dp_1d_array_info, i_1d_array_info, dp_2d_alloc_info, i_2d_alloc_info, &
                     dp_1d_alloc_info

   ! Fine grid
   real(dp), save              :: Xbound(6), Dpm(3)
   integer, save               :: NN_bl(6), NN(3)

   ! Block grid
   real(dp), allocatable, save :: Xbound_block(:, :)
   integer, allocatable, save  :: NN_block(:, :), NN_bl_block(:, :) 

   ! Coarse grid
   real(dp), save              :: Xbound_coarse(6), Dpm_coarse(3)
   integer, save               :: NN_coarse(3), NN_bl_coarse(6)
   
   integer, save               :: NN_tmp(3)
   integer, save               :: nb_i, nb_j, nb_k, BLOCKS, NXB, NYB, NZB, ndumcell_coarse, ndumcell_bl
   real(dp)                    :: starttime, endtime, st, et, ct
   integer, save               :: iynbc, iret, NBI, NBJ, NBK, NREMESH, ntorder, &
                                  iyntree, ilevmax, itree, nsize_out(3), ibctyp, NWRITE

   public :: print_vpm_size_info, get_NN_bl, get_NN, get_Xbound

contains
   subroutine get_NN(NN_out) bind(C, name='get_NN') 
      use iso_c_binding
      implicit none
      integer(c_int), dimension(3) :: NN_out

      NN_out = NN
   End subroutine get_NN

   subroutine get_NN_bl(NN_bl_out) bind(C, name='get_NN_bl') 
      use iso_c_binding
      implicit none
      integer(c_int), dimension(6) :: NN_bl_out
      NN_bl_out = NN_bl
   End subroutine get_NN_bl

   subroutine get_Xbound(Xbound_out) bind(C, name='get_Xbound') 
      use iso_c_binding
      implicit none
      real(c_double), dimension(6) :: Xbound_out      
      Xbound_out = Xbound
   End subroutine get_Xbound

   subroutine print_vpm_size_info()

      print *, "VPM_SIZE INFO"
      print *, "============"
      print *, ""
      call dp_1d_array_info("Xbound", Xbound, 6)
      call dp_1d_array_info("Dpm", Dpm, 3)
      call i_1d_array_info("NN_bl", NN_bl, 6)
      call i_1d_array_info("NN", NN, 3)
      call dp_2d_alloc_info("Xbound_bl", Xbound_block)
      call i_2d_alloc_info("NNbl_bl", NN_bl_block)
      call i_2d_alloc_info("NNbl", NN_block)
      call dp_1d_array_info("Xbound_coarse", Xbound_coarse, 6)
      call dp_1d_array_info("Dpm_coarse", Dpm_coarse, 3)
      call i_1d_array_info("NN_tmp", NN_tmp, 3)
      call i_1d_array_info("NN_coarse", NN_coarse, 3)
      call i_1d_array_info("NN_bl_coarse", NN_bl_coarse, 6)
      print *, achar(9)//"nb_i", " = ", nb_i
      print  *, achar(9)//"nb_j", " = ", nb_j
      print  *, achar(9)//"nb_k", " = ", nb_k
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
