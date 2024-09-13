submodule(vpm_lib) vpm_gcalc
   implicit none

contains
   module subroutine calc_velocity_serial_3d(idcalc)
      ! use vpm_vars
      use pmeshpar, only: SOL_pm
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm, RHS_pm,     &
                        deformx_pm, deformy_pm, deformz_pm,          &
                        DXpm, DYpm, DZpm,                            &
                        NXs_coarse_bl, NYs_coarse_bl, NZs_coarse_bl, &
                        NXf_coarse_bl, NYf_coarse_bl, NZf_coarse_bl
      ! use parvar
      ! use openmpth
      Implicit None
      integer, intent(in) :: idcalc
      real(dp) ::  dpsidx(3), dpsidy(3), dpsidz(3)
      ! real(dp) ::  dphidx, dphidy, dphidz,
      real(dp) ::  wdudx, wdvdy, wdwdz
      real(dp) ::  upi, umi, upj, umj, upk, umk
      real(dp) ::  vpi, vmi, vpj, vmj, vpk, vmk
      real(dp) ::  wpi, wmi, wpj, wmj, wpk, wmk
      real(dp) ::  DXpm2, DYpm2, DZpm2
      integer  :: i, j, k

      if (idcalc .ge. 0) then
         DXpm2 = 2*DXpm
         DYpm2 = 2*DYpm
         DZpm2 = 2*DZpm
         !$omp parallel  &
         !$omp private(i,j,k,dpsidx,dpsidy,dpsidz)  &
         !$omp shared(velvrx_pm,velvry_pm,velvrz_pm,SOL_pm)
         !!$omp          num_threads(OMPTHREADS)
         !$omp do
         do k = NZs_coarse_bl + 1, NZf_coarse_bl - 1
            do j = NYs_coarse_bl + 1, NYf_coarse_bl - 1
               do i = NXs_coarse_bl + 1, NXf_coarse_bl - 1

                  !--> dpsi(x,y,z) / d(xyz)
                  dpsidx(1:3) = (SOL_pm(1:3, i + 1, j, k) - SOL_pm(1:3, i - 1, j, k))/DXpm2
                  dpsidy(1:3) = (SOL_pm(1:3, i, j + 1, k) - SOL_pm(1:3, i, j - 1, k))/DYpm2
                  dpsidz(1:3) = (SOL_pm(1:3, i, j, k + 1) - SOL_pm(1:3, i, j, k - 1))/DZpm2
                  ! U = grad x psi
                  velvrx_pm(i, j, k) = velvrx_pm(i, j, k) + dpsidy(3) - dpsidz(2)
                  velvry_pm(i, j, k) = velvry_pm(i, j, k) - (dpsidx(3) - dpsidz(1))
                  velvrz_pm(i, j, k) = velvrz_pm(i, j, k) + dpsidx(2) - dpsidy(1)
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel

         if (neqpm .eq. 4) then
            !$omp parallel &
            !$omp private(i,j,k,dpsidx,dpsidy,dpsidz) &
            !$omp shared(velvrx_pm,velvry_pm,velvrz_pm,SOL_pm)
            !!$omp num_threads(OMPTHREADS)
            !$omp do
            do k = NZs_coarse_bl + 1, NZf_coarse_bl - 1
               do j = NYs_coarse_bl + 1, NYf_coarse_bl - 1
                  do i = NXs_coarse_bl + 1, NXf_coarse_bl - 1

                     !--> dpsi(x,y,z)d(xyz)
                     dpsidx(1) = (SOL_pm(4, i + 1, j, k) - SOL_pm(4, i - 1, j, k))/DXpm2
                     dpsidy(1) = (SOL_pm(4, i, j + 1, k) - SOL_pm(4, i, j - 1, k))/DYpm2
                     dpsidz(1) = (SOL_pm(4, i, j, k + 1) - SOL_pm(4, i, j, k - 1))/DZpm2
                     ! U = grad x psi
                     velvrx_pm(i, j, k) = velvrx_pm(i, j, k) + dpsidx(1)
                     velvry_pm(i, j, k) = velvry_pm(i, j, k) + dpsidy(1)
                     velvrz_pm(i, j, k) = velvrz_pm(i, j, k) + dpsidz(1)
                  end do
               end do
            end do
            !$omp end do
            !$omp end parallel
         end if
      end if!

      if (idcalc .eq. 0) return

      !Sol of vorticity is no longer need thus we use it for storing deformation OLD NOW WE HAVE VAR
      ! SOL_pm = 0.d0
      !REMEMBER VORTICITY CARRIED IS -OMEGA and the quantity transfered is -OMEGA thus
      !deformation = - omega*\grad u
      deformx_pm = 0.d0
      deformy_pm = 0.d0
      deformz_pm = 0.d0

      !$omp parallel private(upi,upj,upk,vpi,vpj,vpk,wpi,wpj,wpk,umi,umj,umk,vmi,vmj,vmk,&
      !$omp                  wmi,wmj,wmk,          &
      !$omp                  i,j,k,                &  
      !$omp                  wdudx,wdvdy,wdwdz )   & 
      !$omp shared(velvrx_pm,velvry_pm,velvrz_pm,RHS_pm,SOL_pm,deformx_pm,deformy_pm,deformz_pm)
      !!$omp   num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_coarse_bl + 2, NZf_coarse_bl - 2
         do j = NYs_coarse_bl + 2, NYf_coarse_bl - 2
            do i = NXs_coarse_bl + 2, NXf_coarse_bl - 2
               ! velxp = velvrx_pm(i + 1, j, k)
               ! velxm = velvrx_pm(i - 1, j, k)

               ! velyp = velvry_pm(i, j + 1, k)
               ! velym = velvry_pm(i, j - 1, k)

               ! velzp = velvrz_pm(i, j, k + 1)
               ! velzm = velvrz_pm(i, j, k - 1)

               upi = velvrx_pm(i + 1, j, k)
               upj = velvrx_pm(i, j + 1, k)
               upk = velvrx_pm(i, j, k + 1)

               vpi = velvry_pm(i + 1, j, k)
               vpj = velvry_pm(i, j + 1, k)
               vpk = velvry_pm(i, j, k + 1)

               wpi = velvrz_pm(i + 1, j, k)
               wpj = velvrz_pm(i, j + 1, k)
               wpk = velvrz_pm(i, j, k + 1)

               umi = velvrx_pm(i - 1, j, k)
               umj = velvrx_pm(i, j - 1, k)
               umk = velvrx_pm(i, j, k - 1)

               vmi = velvry_pm(i - 1, j, k)
               vmj = velvry_pm(i, j - 1, k)
               vmk = velvry_pm(i, j, k - 1)

               wmi = velvrz_pm(i - 1, j, k)
               wmj = velvrz_pm(i, j - 1, k)
               wmk = velvrz_pm(i, j, k - 1)
               !DEFORMATION WITH A MINUS BECAUSE WE HAVE STORED MINUS VORTICITY
               wdudx = -(RHS_pm(1, i + 1, j, k)*upi - (RHS_pm(1, i - 1, j, k))*umi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*upj - (RHS_pm(2, i, j - 1, k))*umj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*upk - (RHS_pm(3, i, j, k - 1))*umk)/DZpm2

               deformx_pm(i, j, k) = wdudx + wdvdy + wdwdz

               ! Wy * (thu/thx + thv/thy + thw/thz)
               wdudx = -(RHS_pm(1, i + 1, j, k)*vpi - (RHS_pm(1, i - 1, j, k))*vmi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*vpj - (RHS_pm(2, i, j - 1, k))*vmj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*vpk - (RHS_pm(3, i, j, k - 1))*vmk)/DZpm2

               deformy_pm(i, j, k) = wdudx + wdvdy + wdwdz

               ! Wy * (thu/thx + thv/thy + thw/thz)
               wdudx = -(RHS_pm(1, i + 1, j, k)*wpi - (RHS_pm(1, i - 1, j, k))*wmi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*wpj - (RHS_pm(2, i, j - 1, k))*wmj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*wpk - (RHS_pm(3, i, j, k - 1))*wmk)/DZpm2

               deformz_pm(i, j, k) = wdudx + wdvdy + wdwdz
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel
   end subroutine calc_velocity_serial_3d

   module subroutine diffuse_vort_3d
      use vpm_vars
      use pmeshpar
      use parvar
      use pmgrid
      use openmpth
      Implicit None
      real(dp) ::  dwxdx, dwydy, dwzdz, VIS
      real(dp) ::  DXpm2, DYpm2, DZpm2 
      integer  :: i, j, k

      DXpm2 = DXpm**2
      DYpm2 = DYpm**2
      DZpm2 = DZpm**2
      deformx_pm = 0.d0
      deformy_pm = 0.d0
      deformz_pm = 0.d0
      !$omp parallel private(i,j,k,dwxdx,dwydy,dwzdz,VIS) shared(RHS_pm,deformx_pm,deformy_pm,deformz_pm)
      !!$omp num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_coarse_bl + 1, NZf_coarse_bl - 1
         do j = NYs_coarse_bl + 1, NYf_coarse_bl - 1
            do i = NXs_coarse_bl + 1, NXf_coarse_bl - 1
               if (neqpm .eq. 3) then
                  VIS = NI
               else
                  VIS = RHS_pm(4, i, j, k) + NI
               end if
               !--> Remember that RHS = -w
               dwxdx = (RHS_pm(1, i + 1, j, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i - 1, j, k))/DXpm2
               dwydy = (RHS_pm(1, i, j + 1, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j - 1, k))/DYpm2
               dwzdz = (RHS_pm(1, i, j, k + 1) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformx_pm(i, j, k) = -VIS*(dwxdx + dwydy + dwzdz) ! because RHS=-w

               dwxdx = (RHS_pm(2, i + 1, j, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i - 1, j, k))/DXpm2
               dwydy = (RHS_pm(2, i, j + 1, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j - 1, k))/DYpm2
               dwzdz = (RHS_pm(2, i, j, k + 1) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformy_pm(i, j, k) = -VIS*(dwxdx + dwydy + dwzdz) ! because RHS=-w

               dwxdx = (RHS_pm(3, i + 1, j, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i - 1, j, k))/DXpm2
               dwydy = (RHS_pm(3, i, j + 1, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j - 1, k))/DYpm2
               dwzdz = (RHS_pm(3, i, j, k + 1) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformz_pm( i, j, k) = -VIS*(dwxdx + dwydy + dwzdz) ! because RHS=-w
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel

      !!$omp parallel private(i,j,k) num_threads(OMPTHREADS)
      !!$omp do
      !do k = NZs_bl + 1, NZf_bl- 1
      !    do j = NYs_bl + 1, NYf_bl(1 )- 1
      !       do i = NXs_bl + 1, NXf_bl - 1

      !            !--> Remember that RHS = -w
      !            RHS_pm(1, i, j, k) = RHS_pm(1, i, j, k) - NI * deformx_pm(i, j, k)
      !            RHS_pm(2, i, j, k) = RHS_pm(2, i, j, k) - NI * deformy_pm(i, j, k)
      !            RHS_pm(3, i, j, k) = RHS_pm(3, i, j, k) - NI * deformz_pm(i, j, k)
      !        enddo
      !    enddo
      !enddo
      !!$omp enddo
      !!$omp endparallel
   end subroutine diffuse_vort_3d

   module subroutine calc_antidiffusion
      use vpm_vars
      use pmeshpar
      use parvar
      use pmgrid
      use openmpth
      Implicit None
      real(dp)                ::  dwxdx, dwydy, dwzdz, Ct
      real(dp)                ::  DXpm2, DYpm2, DZpm2
      integer                 :: i, j, k
      real(dp), allocatable   :: laplvort(:, :, :, :)

      allocate (laplvort(3, NXpm_coarse, NYpm_coarse, NZpm_coarse))
      laplvort = 0.d0
      DXpm2 = DXpm**2
      DYpm2 = DYpm**2
      DZpm2 = DZpm**2
      Sol_pm = 0.d0
      Ct = 6.8*DXpm**2/4
      do k = NZs_coarse_bl + 1, NZf_coarse_bl - 1
         do j = NYs_coarse_bl + 1, NYf_coarse_bl - 1
            do i = NXs_coarse_bl + 1, NXf_coarse_bl - 1
               !--> Remember that RHS = -w
               dwxdx = -(RHS_pm(1, i + 1, j, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i - 1, j, k))/DXpm2
               dwydy = -(RHS_pm(1, i, j + 1, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j - 1, k))/DYpm2
               dwzdz = -(RHS_pm(1, i, j, k + 1) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j, k - 1))/DZpm2
               ! U = grad x psi
               laplvort(1, i, j, k) = dwxdx + dwydy + dwzdz

               dwxdx = -(RHS_pm(2, i + 1, j, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i - 1, j, k))/DXpm2
               dwydy = -(RHS_pm(2, i, j + 1, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j - 1, k))/DYpm2
               dwzdz = -(RHS_pm(2, i, j, k + 1) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j, k - 1))/DZpm2
               ! U = grad x psi
               laplvort(2, i, j, k) = dwxdx + dwydy + dwzdz

               dwxdx = -(RHS_pm(3, i + 1, j, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i - 1, j, k))/DXpm2
               dwydy = -(RHS_pm(3, i, j + 1, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j - 1, k))/DYpm2
               dwzdz = -(RHS_pm(3, i, j, k + 1) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm2
               ! U = grad x psi
               laplvort(3, i, j, k) = dwxdx + dwydy + dwzdz
            end do
         end do
      end do

      do k = NZs_coarse_bl + 1, NZf_coarse_bl - 1
         do j = NYs_coarse_bl + 1, NYf_coarse_bl - 1
            do i = NXs_coarse_bl + 1, NXf_coarse_bl - 1
               !Minus because of (-w) has been included in laplvort
               dwxdx = (laplvort(1, i + 1, j, k) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i - 1, j, k))/DXpm2
               dwydy = (laplvort(1, i, j + 1, k) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i, j - 1, k))/DYpm2
               dwzdz = (laplvort(1, i, j, k + 1) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformx_pm(i, j, k) = deformx_pm(i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplvort(2, i + 1, j, k) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i - 1, j, k))/DXpm2
               dwydy = (laplvort(2, i, j + 1, k) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i, j - 1, k))/DYpm2
               dwzdz = (laplvort(2, i, j, k + 1) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformy_pm(i, j, k) = deformy_pm(i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplvort(3, i + 1, j, k) - 2*laplvort(3, i, j, k) &
                        + laplvort(3, i - 1, j, k))/DXpm2
               dwydy = (laplvort(3, i, j + 1, k) - 2*laplvort(3, i, j, k) &
                        + laplvort(3, i, j - 1, k))/DYpm2
               dwzdz = (laplvort(3, i, j, k + 1) - 2*laplvort(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm2
               ! U = grad x psi
               deformz_pm(i, j, k) = deformz_pm(i, j, k) + Ct*dwxdx + dwydy + dwzdz
            end do
         end do
      end do

   end subroutine calc_antidiffusion

end submodule vpm_gcalc