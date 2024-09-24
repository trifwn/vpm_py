Module vpm_size
   use base_types, only: dp      
   use console_io, only: dp_1d_array_info, i_1d_array_info, dp_2d_alloc_info, i_2d_alloc_info, &
                     dp_1d_alloc_info

   ! Fine grid
   real(dp), save              :: Xbound(6), Dpm(3)
   integer, save               :: NN_bl(6), NN(3)

   ! Block grid 
   real(dp), allocatable       :: SOL_pm_bl(:, :, :, :), RHS_pm_bl(:, :, :, :)
   !Blcok info 1st dimension is the block number, after that the same as fine grid
   real(dp), allocatable, save :: Xbound_block(:, :)
   integer, allocatable, save  :: NN_block(:, :), NN_bl_block(:, :) 
   integer                     :: my_block

   ! Coarse grid
   real(dp), save              :: Xbound_coarse(6), Dpm_coarse(3)
   integer, save               :: NN_coarse(3), NN_bl_coarse(6)
   
   integer, save               :: NN_tmp(3)
   integer, save               :: nb_i, nb_j, nb_k, BLOCKS, NXB, NYB, NZB, ndumcell_coarse, ndumcell_bl
   real(dp)                    :: st, et
   integer, save               :: iynbc, iret, NBI, NBJ, NBK, nremesh, &
                                  iyntree, ilevmax, itree, ibctyp

   public :: print_vpm_size_info, get_NN_bl, get_NN, get_Xbound

contains
   !> Defines Coarse and Fine GRID
   subroutine define_sizes
      use vpm_vars, only:  interf_iproj, idefine, NTIME_PM
      use pmgrid, only:    XMIN_pm, YMIN_pm, ZMIN_pm, &
                           XMAX_pm, YMAX_pm, ZMAX_pm, &
                           NXs_coarse_bl, NXf_coarse_bl, NYs_coarse_bl, &
                           NYf_coarse_bl, NZs_coarse_bl, NZf_coarse_bl, &
                           NXpm_coarse, NYpm_coarse, NZpm_coarse, DVPM, &
                           DXpm, DYpm, DZpm, ND, ndumcell, ncoarse
      use pmlib, only: definepm
      use parvar, only: XP
      use console_io, only: dummy_string, vpm_print, nocolor, blue, yellow, red, tab_level
      use MPI

      Implicit None
      integer    :: nsiz(3), nsiz_bl(3)
      integer    :: i, j, k, np, my_rank, ierr, nb, istep, lev
      real(dp)   :: Xbound_tmp(6)
      integer    :: redifine_pm
      integer    :: NN_bl_tmp(6), NXbl, NYbl, NZbl 

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      if ((NTIME_pm .eq. 0) .or. idefine .eq. 0) then
         redifine_pm = 1
      else
         redifine_pm = 0
      end if

      if(my_rank.eq.0) then
         st = MPI_WTIME()
         write (dummy_string, *) 'Defining Sizes'
         call vpm_print(dummy_string, blue, 2)
         write (dummy_string, '(A, I5, A, I1)') achar(9)//'NTIME_PM=', NTIME_pm, &
         achar(9)//'Redefine=', redifine_pm
         call vpm_print(dummy_string,nocolor, 2)
      end if

      tab_level = tab_level + 1

      if (my_rank .eq. 0) then
         print *, "minval X" , minval(XP(1, :))
         if (redifine_pm.eq.1) then
            call get_domain_bounds_from_particles
            write (dummy_string, '(A)') 'The computational domain bounds are recalculated from the particle positions'
            call vpm_print(dummy_string,red, 2)
            ! Write the min and max values of the particle positions
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                     achar(9)//'Particle XMIN=', minval(XP(1, :)), &
                     achar(9)//'Particle XMAX=', maxval(XP(1, :))
            call vpm_print(dummy_string,yellow, 2)
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                     achar(9)//'Particle YMIN=', minval(XP(2, :)), &
                     achar(9)//'Particle YMAX=', maxval(XP(2, :))
            call vpm_print(dummy_string,yellow, 2)
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                     achar(9)//'Particle ZMIN=', minval(XP(3, :)), &
                     achar(9)//'Particle ZMAX=', maxval(XP(3, :))
            call vpm_print(dummy_string,yellow, 2)
         end if

         write (dummy_string, '(A)') 'The computational domain bounds are:'
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, F10.5, A, F8.5, A, F10.5)') &
                  achar(9)//'XMIN=', XMIN_pm, &
                  achar(9)//'YMIN=', YMIN_pm, &
                  achar(9)//'ZMIN=', ZMIN_pm
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                  achar(9)//'XMAX=', XMAX_pm, &
                  achar(9)//'YMAX=', YMAX_pm, &
                  achar(9)//'ZMAX=', ZMAX_pm
         call vpm_print(dummy_string,nocolor, 2)
      end if

      ! BRADCAST THE MIN AND MAX OF THE COMPUTATIONAL DOMAIN
      call MPI_BCAST(XMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_BCAST(XMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      
      ! CREATE THE BOUNDING BOX
      Xbound(1) = XMIN_pm
      Xbound(2) = YMIN_pm
      Xbound(3) = ZMIN_pm
      Xbound(4) = XMAX_pm 
      Xbound(5) = YMAX_pm 
      Xbound(6) = ZMAX_pm
      
      Dpm(1) = DXpm
      Dpm(2) = DYpm
      Dpm(3) = DZpm
      
      nsiz(1) = NBI*ncoarse
      nsiz(2) = NBJ*ncoarse
      nsiz(3) = NBK*ncoarse

      if(my_rank.eq.0) then
         write (dummy_string, '(A)') 'The extended fine domain is defined (no dummy cells)'
         call vpm_print(dummy_string,red, 2)
         write (dummy_string, '(A)') achar(9)//'The number of cells (nodes-1) must be divisible by the processor subdivision:'
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//achar(9)//'X-dir: ', nsiz(1), &
                  achar(9)//achar(9)//'Y-dir: ', nsiz(2), &
                  achar(9)//achar(9)//'Z-dir: ', nsiz(3)
         call vpm_print(dummy_string,nocolor, 2)
      endif
      
      ! First Change Dpm so that the numbers of cells divides (that would be case 2)
      ! Here we change the number of cells
      ! by nsize i.e with NBI,NBJ,NBK, ncoarse,levmax depending on the criterion
      ! thats why ndumcell=0
      ndumcell = 0
      if (redifine_pm .eq. 1) then
         call definepm(3, Xbound, Dpm, ND, ndumcell, nsiz, NN, NN_bl)
      endif

      ! THE new extended domain is defined
      XMIN_pm = Xbound(1)
      YMIN_pm = Xbound(2)
      ZMIN_pm = Xbound(3)
      XMAX_pm = Xbound(4)
      YMAX_PM = Xbound(5)
      ZMAX_pm = Xbound(6)

      NXpm_coarse = NN(1)
      NYpm_coarse = NN(2)
      NZpm_coarse = NN(3)

      NXs_coarse_bl = NN_bl(1)
      NYs_coarse_bl = NN_bl(2)
      NZs_coarse_bl = NN_bl(3)
      NXf_coarse_bl = NN_bl(4)
      NYf_coarse_bl = NN_bl(5)
      NZf_coarse_bl = NN_bl(6)

      DXpm = Dpm(1)
      DYpm = Dpm(2)
      DZpm = Dpm(3)

      DVpm = DXpm*DYpm
      if (ND .eq. 3) then
         DVpm = DVpm*DZpm
      end if

      if (my_rank .eq. 0) then
         write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                  achar(9)//'XMIN='//achar(9), XMIN_pm, &
                  achar(9)//'YMIN='//achar(9), YMIN_pm, &
                  achar(9)//'ZMIN='//achar(9), ZMIN_pm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F10.5, A, F10.5, A, F8.5)') &
                  achar(9)//'XMAX='//achar(9), XMAX_pm, &
                  achar(9)//'YMAX='//achar(9), YMAX_pm, &
                  achar(9)//'ZMAX='//achar(9), ZMAX_pm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                  achar(9)//'DX='//achar(9), DXpm, &
                  achar(9)//'DY='//achar(9), DYpm, &
                  achar(9)//'DZ='//achar(9), DZpm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'Nodes X= ', NXpm_coarse, &
                  achar(9)//achar(9)//'Nodes Y= ', NYpm_coarse, &
                  achar(9)//achar(9)//'Nodes Z= ', NZpm_coarse
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, "(A)") "The indexes of the coarse grid that do not include the dummy cells are:"
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'NXs_coarse=', NXs_coarse_bl, &
                  achar(9)//'NYs_coarse=', NYs_coarse_bl, &
                  achar(9)//'NZs_coarse=', NZs_coarse_bl
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'NXf_coarse=', NXf_coarse_bl, &
                  achar(9)//'NYf_coarse=', NYf_coarse_bl, &
                  achar(9)//'NZf_coarse=', NZf_coarse_bl
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5)') achar(9) // 'DV=', DVpm
         call vpm_print(dummy_string,nocolor,2)
      end if

      ! define block grids
      ! so they are divided by ncoarse and ilevmax
      ! so as to have the coarse information at the boundaries exactly.
      BLOCKS = np
      if (.not. allocated(Xbound_block)) then
         allocate (Xbound_block(6, BLOCKS), NN_block(3, BLOCKS), NN_bl_block(6, BLOCKS))
      end if
      
      ! DEVIDE THE DOMAIN INTO EQUAL BLOCKS 
      NXB = int(nint(((Xbound(4) - Xbound(1))/Dpm(1)))) ! NUMBER OF CELLS IN A BLOCK IN X DIRECTION 
      NYB = int(nint(((Xbound(5) - Xbound(2))/Dpm(2)))) ! NUMBER OF CELLS IN A BLOCK IN Y DIRECTION
      NZB = int(nint(((Xbound(6) - Xbound(3))/Dpm(3)))) ! NUMBER OF CELLS IN A BLOCK IN Z DIRECTION
      NXbl = NXB / NBI
      NYbl = NYB / NBJ
      NZbl = NZB / NBK
      
      ndumcell_bl = ncoarse
      nsiz_bl = ncoarse!ndumcell_coarse!*2*2**ilevmax

      if (my_rank.eq.0) then
         write (dummy_string, '(A)') 'The fine block domains are defined for each processor'
         call vpm_print(dummy_string,red, 2)
         write(dummy_string, '(A)') achar(9) // 'The number of cells in each direction for the block grid'
         call vpm_print(dummy_string,nocolor,2)
         write(dummy_string, '(A, I5, A, I5, A, I5)') achar(9) // &
                  achar(9)//'Total Num Cells:'//                   &
                  achar(9)//'X='//achar(9), NXB,                  &
                  achar(9)//'Y='//achar(9), NYB,                  &
                  achar(9)//'Z='//achar(9), NZB       
         call vpm_print(dummy_string,nocolor,2)                          
         write(dummy_string, '(A, I5, A, I5, A, I5)') achar(9) // &
                  achar(9)// 'Cells per Block:'//                  &
                  achar(9)//'X='//achar(9), NXbl,                 &
                  achar(9)//'Y='//achar(9), NYbl,                 &
                  achar(9)//'Z='//achar(9), NZbl
         call vpm_print(dummy_string,nocolor,2)
         write(dummy_string, '(A, I5)') achar(9)//achar(9)//'Block dummy cells:', ndumcell_bl
         call vpm_print(dummy_string,nocolor,2)
      end if

      nb_i = -1
      nb_j = -1
      nb_k = -1
      do k = 1, NBK
         do j = 1, NBJ
            do i = 1, NBI
               nb = (k - 1)*NBJ*NBI + (j - 1)*NBI + i
               Xbound_block(1, nb) = Xbound(1) + (i - 1)*(NXbl)*Dpm(1)
               Xbound_block(4, nb) = Xbound(1) + (i)*(NXbl)*Dpm(1)
               
               Xbound_block(2, nb) = Xbound(2) + (j - 1)*(NYbl)*Dpm(2)
               Xbound_block(5, nb) = Xbound(2) + (j)*(NYbl)*Dpm(2)
               
               Xbound_block(3, nb) = Xbound(3) + (k - 1)*(NZbl)*Dpm(3)
               Xbound_block(6, nb) = Xbound(3) + (k)*(NZbl)*Dpm(3)
               Xbound_tmp(1:6) = Xbound_block(1:6, nb)
               call definepm(1, Xbound_tmp, Dpm, ND, ndumcell_bl, nsiz_bl, NN_tmp, NN_bl_tmp)
               
               Xbound_block(1:6, nb)   = Xbound_tmp(1:6)  ! Block domain boundaries
               NN_block(1:3, nb)        = NN_tmp(1:3)      ! Block number of nodes
               NN_bl_block(1:6, nb)     = NN_bl_tmp(1:6)   ! Indices of the block domain that do not include the dummy cells
               
               if (my_rank .eq. 0) then
                  write (dummy_string, '(A,I3,A,3I3,A)') achar(9)//'Block ', nb, " = (", i, j, k, ")"
                  call vpm_print(dummy_string,blue,2)
                  write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') achar(9)//&
                           achar(9)//'XMIN= '//achar(9), Xbound_block(1, nb), &
                           achar(9)//'YMIN= '//achar(9), Xbound_block(2, nb), &
                           achar(9)//'ZMIN= '//achar(9), Xbound_block(3, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') achar(9)// &
                           achar(9)//'XMAX= '//achar(9), Xbound_block(4, nb), &
                           achar(9)//'YMAX= '//achar(9), Xbound_block(5, nb), &
                           achar(9)//'ZMAX= '//achar(9), Xbound_block(6, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, '(A, I5, A, I5, A, I5)') achar(9)//achar(9)//&
                           "Cells X:"//achar(9), NN_block(1, nb),      &
                           achar(9)//"Y:"//achar(9), NN_block(2, nb),      &
                           achar(9)//"Z:"//achar(9), NN_block(3, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, "(A)") achar(9) // achar(9) // &
                           "Real Cells (not dummy):"
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, '(A, I5, A, I5, A, I5)') achar(9)//achar(9)//&
                           achar(9)//'Xs:'//achar(9), NN_bl_block(1, nb), &
                           achar(9)//'Ys:'//achar(9), NN_bl_block(2, nb), &
                           achar(9)//'Zs:'//achar(9), NN_bl_block(3, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, '(A, I5, A, I5, A, I5)') achar(9)//achar(9)//&
                           achar(9)//'Xf:'//achar(9), NN_bl_block(4, nb), &
                           achar(9)//'Yf:'//achar(9), NN_bl_block(5, nb), &
                           achar(9)//'Zf:'//achar(9), NN_bl_block(6, nb)
                  call vpm_print(dummy_string,nocolor,2)
               end if

               if (nb .eq. my_rank + 1) then
                  nb_i = i
                  nb_j = j
                  nb_k = k
               end if
            end do
         end do
      end do
      if ((nb_i .eq. -1) .or. (nb_j .eq. -1) .or. (nb_k .eq. -1)) then
         print *, "Processor with rank ", my_rank, " does not have a block"
         stop
      end if

      !define coarse grid must cover block grids
      Xbound(1) = XMIN_pm!minval(Xbound_bl(1,:))
      Xbound(2) = YMIN_pm!minval(Xbound_bl(2,:))
      Xbound(3) = ZMIN_pm!minval(Xbound_bl(3,:))
      Xbound(4) = XMAX_pm!maxval(Xbound_bl(4,:))
      Xbound(5) = YMAX_pm!maxval(Xbound_bl(5,:))
      Xbound(6) = ZMAX_pm!maxval(Xbound_bl(6,:))
      Xbound_coarse = Xbound
      Dpm_coarse = ncoarse*Dpm
      ndumcell_coarse = 4!2**ilevmax
      nsiz_bl = 2**ilevmax
      call definepm(1, Xbound_coarse, Dpm_coarse, ND, ndumcell_coarse, nsiz_bl, NN_coarse, NN_bl_coarse)

      if (my_rank .eq. 0) then
         write (dummy_string, '(A)') 'The extended coarse domain is redefined (with dummy cells)'
         call vpm_print(dummy_string,red, 2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMIN='//achar(9), Xbound_coarse(1), &
                  achar(9)//'YMIN='//achar(9), Xbound_coarse(2), &
                  achar(9)//'ZMIN='//achar(9), Xbound_coarse(3)
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMAX='//achar(9), Xbound_coarse(4), &
                  achar(9)//'YMAX='//achar(9), Xbound_coarse(5), &
                  achar(9)//'ZMAX='//achar(9), Xbound_coarse(6)
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'DX='//achar(9), DXpm, &
                  achar(9)//'DY='//achar(9), DYpm, &
                  achar(9)//'DZ='//achar(9), DZpm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'Nodes X= ',           NN_coarse(1), &
                  achar(9)//achar(9)//'Nodes Y= ', NN_coarse(2), &
                  achar(9)//achar(9)//'Nodes Z= ', NN_coarse(3)
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, "(A)") "The indexes of the coarse grid that do not include the dummy cells are:"
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'NXs_coarse=', NN_bl_coarse(1), &
                  achar(9)//'NYs_coarse=', NN_bl_coarse(2), &
                  achar(9)//'NZs_coarse=', NN_bl_coarse(3)
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'NXf_coarse=', NN_bl_coarse(4), &
                  achar(9)//'NYf_coarse=', NN_bl_coarse(5), &
                  achar(9)//'NZf_coarse=', NN_bl_coarse(6)
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5)') achar(9) // 'DV=', DVpm
         call vpm_print(dummy_string,nocolor,2)
      end if
      
      ! CALCULATE THE NUMBER OF LEVELS
      if (my_rank .eq. 0) then
         do lev = 0, ilevmax
            istep = 2**lev
            !!!!!!!!if not divided exactly dummy cell
            if (ND .eq. 2) then
               if (int((NN_bl_block(4, 1) - NN_bl_block(1, 1))/istep) .eq. 0 .or. &
               int((NN_bl_block(5, 1) - NN_bl_block(2, 1))/istep) .eq. 0) then
                  ilevmax = lev - 1
                  print *, 'Changing number of levels', ilevmax
                  exit
               end if
            else
               if (  int((NN_bl_block(4, 1) - NN_bl_block(1, 1))/istep) .eq. 0 .or. &
                     int((NN_bl_block(5, 1) - NN_bl_block(2, 1))/istep) .eq. 0 .or. &
                     int((NN_bl_block(6, 1) - NN_bl_block(3, 1))/istep) .eq. 0) then
                  ilevmax = lev - 1
                  print *, 'Changing number of levels', ilevmax
                  exit
               end if
            end if
         end do
      end if
      call MPI_BCAST(ilevmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      tab_level = tab_level - 1
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, *) 'Defining Sizes finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         call vpm_print(dummy_string, yellow, 2)
      end if      
   end subroutine define_sizes

   subroutine get_domain_bounds_from_particles
      use parvar, only: NVR, XP
      use vpm_vars, only: interf_iproj
      use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm, XMAX_pm, &
                        YMAX_pm, ZMAX_pm
      real(dp) :: XMIN_pm_old, YMIN_pm_old, ZMIN_pm_old, XMAX_pm_old, YMAX_pm_old, ZMAX_pm_old, &
                  X_mean, Y_mean, Z_mean
      logical  :: bounds_changed

      XMIN_pm_old = XMIN_pm
      YMIN_pm_old = YMIN_pm
      ZMIN_pm_old = ZMIN_pm
      XMAX_pm_old = XMAX_pm
      YMAX_pm_old = YMAX_pm
      ZMAX_pm_old = ZMAX_pm

      ! Get mean X, Y, Z
      X_mean = sum(XP(1, 1:NVR))/NVR
      Y_mean = sum(XP(2, 1:NVR))/NVR
      Z_mean = sum(XP(3, 1:NVR))/NVR

      XMIN_pm = minval(XP(1, 1:NVR)) - interf_iproj*DXpm
      YMIN_pm = minval(XP(2, 1:NVR)) - interf_iproj*DYpm
      ZMIN_pm = minval(XP(3, 1:NVR)) - interf_iproj*DZpm

      XMAX_pm = maxval(XP(1, 1:NVR)) + interf_iproj*DXpm
      YMAX_pm = maxval(XP(2, 1:NVR)) + interf_iproj*DYpm
      ZMAX_pm = maxval(XP(3, 1:NVR)) + interf_iproj*DZpm

      ! Normalize the domain bounds to be multiple of DX, DY, DZ
      XMIN_pm = X_mean - ceiling((X_mean - XMIN_pm)/DXpm)*DXpm
      YMIN_pm = Y_mean - ceiling((Y_mean - YMIN_pm)/DYpm)*DYpm
      ZMIN_pm = Z_mean - ceiling((Z_mean - ZMIN_pm)/DZpm)*DZpm
      
      XMAX_pm = X_mean + ceiling((XMAX_pm - X_mean)/DXpm)*DXpm
      YMAX_pm = Y_mean + ceiling((YMAX_pm - Y_mean)/DYpm)*DYpm
      ZMAX_pm = Z_mean + ceiling((ZMAX_pm - Z_mean)/DZpm)*DZpm

      ! Check if the domain bounds have changed
      if (XMIN_pm .ne. XMIN_pm_old .or. YMIN_pm .ne. YMIN_pm_old .or. ZMIN_pm .ne. ZMIN_pm_old .or. &
          XMAX_pm .ne. XMAX_pm_old .or. YMAX_pm .ne. YMAX_pm_old .or. ZMAX_pm .ne. ZMAX_pm_old) then
         bounds_changed =  .TRUE.
      else
         bounds_changed = .FALSE.
      end if
   end subroutine get_domain_bounds_from_particles

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
      print  *, achar(9)//"NREMESH", " = ", nremesh
      print  *, achar(9)//"iyntree", " = ", iyntree
      print  *, achar(9)//"ilevmax", " = ", ilevmax
      print  *, achar(9)//"itree", " = ", itree
      print  *, achar(9)//"ibctyp", " = ", ibctyp
   end subroutine print_vpm_size_info
End Module vpm_size
