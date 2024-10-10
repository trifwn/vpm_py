module vpm_gcalc
   use base_types, only: dp
   use serial_vector_field_operators
   implicit none

contains
   !> \brief Calculates the velocity in a 3D space for a given calculation ID.
   !>
   !> This subroutine performs the velocity calculation in a serial manner for a 
   !> three-dimensional space. It takes an identifier for the calculation as input 
   !> and processes the velocity based on the provided ID. All calculations are performed
   !> using central difference schemes for the gradient approximation.
   !>
   !> \param[in] idcalc An integer representing the calculation ID. If zero, then only velocity 
   !>             calculation is performed. If greater than zero, then the velocity calculation
   !>             is performed along with the deformation calculation. If less than zero, then
   !>             only the deformation calculation is performed.
   subroutine calc_velocity_serial_3d(idcalc)
      use MPI
      use pmgrid, only: velocity_pm, deform_pm, RHS_pm,        &
                        DXpm, DYpm, DZpm,                      &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        SOL_pm, NXpm_fine, NYpm_fine, NZpm_fine
      use vpm_vars, only: neqpm
      use console_io, only: vpm_print, blue, yellow, dummy_string
      use base_types, only: dp
      Implicit None
      integer, intent(in) :: idcalc
      real(dp) ::  dphidx, dphidy, dphidz
      ! real(dp) ::  dphidx, dphidy, dphidz,
      real(dp) ::  wdudx, wdvdy, wdwdz
      real(dp) ::  upi, umi, upj, umj, upk, umk
      real(dp) ::  vpi, vmi, vpj, vmj, vpk, vmk
      real(dp) ::  wpi, wmi, wpj, wmj, wpk, wmk
      real(dp) ::  DXpm2, DYpm2, DZpm2
      real(dp), allocatable :: deformation(:,:,:,:)
      real(dp), allocatable :: jac_u(:,:,:,:,:)
      real(dp)              :: total_error, max_error
      integer  :: i, j, k
      real(dp) :: st, et
      logical :: caluclate_velocity, calculate_deformation

      caluclate_velocity = idcalc >= 0
      calculate_deformation = idcalc /= 0

      st = MPI_WTIME()
      write (*, *) ""
      if (idcalc .eq. 0) then
         write (dummy_string, "(A)") "Calculating Velocities on PM using FD"
      else if (idcalc .gt. 0) then
         write (dummy_string, "(A)") "Calculating Velocities and Deformations on PM using FD"
      else
         write (dummy_string, "(A)") 'Calculating Deformations on PM using FD'
      end if
      call vpm_print(dummy_string, blue, 1)

      DXpm2 = 2*DXpm
      DYpm2 = 2*DYpm
      DZpm2 = 2*DZpm
      if (caluclate_velocity) then
         ! Velocity due to vorticity
         velocity_pm = curl(SOL_pm(1:3,:,:,:), DXpm, DYpm, DZpm)
         
         ! Velocity due to irrotational part
         if (neqpm .eq. 4) then
            velocity_pm = velocity_pm + gradient(SOL_pm(4,:,:,:), DXpm, DYpm, DZpm) 
         endif
      end if
      
      total_error = 0.d0
      max_error = 0.d0
      if (calculate_deformation) then
      ! if (.true.) then
         !REMEMBER VORTICITY CARRIED IS -OMEGA and the quantity transfered is -OMEGA thus
         !deformation = - (omega*\nabla) u =  \nabla(w x u) - u * (\grad w) = 
         ! We will calculate the deformation by using \grad(w x u) 
         deform_pm = 0.d0
         jac_u = jacobian(velocity_pm, DXpm, DYpm, DZpm)
         ! Jac_u is a 5D array with the following dimensions
         ! jac_u(3, 3, NXpm_fine, NYpm_fine, NZpm_fine)
         ! jac_u(1, 1, i, j, k) = du/dx
         ! jac_u(1, 2, i, j, k) = du/dy
         ! jac_u(1, 3, i, j, k) = du/dz
         ! jac_u(2, 1, i, j, k) = dv/dx
         ! ... etc
         allocate(deformation(3, NXpm_fine, NYpm_fine, NZpm_fine))
         deformation = 0.d0
         !$omp parallel private(upi,upj,upk,vpi,vpj,vpk,wpi,wpj,wpk,umi,umj,umk,vmi,vmj,vmk,&
         !$omp                  wmi,wmj,wmk,          &
         !$omp                  i,j,k,                &  
         !$omp                  wdudx,wdvdy,wdwdz )   & 
         !$omp shared(velocity_pm,RHS_pm,SOL_pm,deform_pm, &
         !$omp        deformation, total_error, max_error)
         !!$omp   num_threads(OMPTHREADS)
         !$omp do
         do k = NZs_fine_bl + 2, NZf_fine_bl - 2
            do j = NYs_fine_bl + 2, NYf_fine_bl - 2
               do i = NXs_fine_bl + 2, NXf_fine_bl - 2
                  ! deform_x = du/dx * w1 + dv/dx * w2 + dw/dx * w3
                  ! deform_y = du/dy * w1 + dv/dy * w2 + dw/dy * w3
                  ! deform_z = du/dz * w1 + dv/dz * w2 + dw/dz * w3
                  deformation(1,i,j,k) = dot_product(jac_u(1, : , i, j, k), RHS_pm(1:3, i, j, k))
                  deformation(2,i,j,k) = dot_product(jac_u(2, : , i, j, k), RHS_pm(1:3, i, j, k))
                  deformation(3,i,j,k) = dot_product(jac_u(3, : , i, j, k), RHS_pm(1:3, i, j, k))

                  upi = velocity_pm(1, i + 1, j, k)
                  upj = velocity_pm(1, i, j + 1, k)
                  upk = velocity_pm(1, i, j, k + 1)

                  vpi = velocity_pm(2, i + 1, j, k)
                  vpj = velocity_pm(2, i, j + 1, k)
                  vpk = velocity_pm(2, i, j, k + 1)

                  wpi = velocity_pm(3, i + 1, j, k)
                  wpj = velocity_pm(3, i, j + 1, k)
                  wpk = velocity_pm(3, i, j, k + 1)

                  umi = velocity_pm(1, i - 1, j, k)
                  umj = velocity_pm(1, i, j - 1, k)
                  umk = velocity_pm(1, i, j, k - 1)

                  vmi = velocity_pm(2, i - 1, j, k)
                  vmj = velocity_pm(2, i, j - 1, k)
                  vmk = velocity_pm(2, i, j, k - 1)

                  wmi = velocity_pm(3, i - 1, j, k)
                  wmj = velocity_pm(3, i, j - 1, k)
                  wmk = velocity_pm(3, i, j, k - 1)
                  !DEFORMATION WITH A MINUS BECAUSE WE HAVE STORED MINUS VORTICITY
                  wdudx = -(RHS_pm(1, i + 1, j, k)*upi - (RHS_pm(1, i - 1, j, k))*umi)/DXpm2
                  wdvdy = -(RHS_pm(2, i, j + 1, k)*upj - (RHS_pm(2, i, j - 1, k))*umj)/DYpm2
                  wdwdz = -(RHS_pm(3, i, j, k + 1)*upk - (RHS_pm(3, i, j, k - 1))*umk)/DZpm2

                  deform_pm(1,i, j, k) = wdudx + wdvdy + wdwdz

                  ! Wy * (thu/thx + thv/thy + thw/thz)
                  wdudx = -(RHS_pm(1, i + 1, j, k)*vpi - (RHS_pm(1, i - 1, j, k))*vmi)/DXpm2
                  wdvdy = -(RHS_pm(2, i, j + 1, k)*vpj - (RHS_pm(2, i, j - 1, k))*vmj)/DYpm2
                  wdwdz = -(RHS_pm(3, i, j, k + 1)*vpk - (RHS_pm(3, i, j, k - 1))*vmk)/DZpm2

                  deform_pm(2,i, j, k) = wdudx + wdvdy + wdwdz

                  ! Wy * (thu/thx + thv/thy + thw/thz)
                  wdudx = -(RHS_pm(1, i + 1, j, k)*wpi - (RHS_pm(1, i - 1, j, k))*wmi)/DXpm2
                  wdvdy = -(RHS_pm(2, i, j + 1, k)*wpj - (RHS_pm(2, i, j - 1, k))*wmj)/DYpm2
                  wdwdz = -(RHS_pm(3, i, j, k + 1)*wpk - (RHS_pm(3, i, j, k - 1))*wmk)/DZpm2

                  deform_pm(3,i, j, k) = wdudx + wdvdy + wdwdz

                  total_error = total_error + abs(deform_pm(1,i,j,k) - deformation(1, i,j,k)) + &
                                              abs(deform_pm(2,i,j,k) - deformation(2, i,j,k)) + &
                                              abs(deform_pm(3,i,j,k) - deformation(3, i,j,k))
                  
                  if (max_error .lt. abs(deform_pm(1,i, j, k) - deformation(1, i, j, k))) then
                     max_error = abs(deform_pm(1,i, j, k) - deformation(1, i, j, k))
                  end if

                  if (max_error .lt. abs(deform_pm(2,i, j, k) - deformation(2, i, j, k))) then
                     max_error = abs(deform_pm(2,i, j, k) - deformation(2, i, j, k))
                  end if

                  if (max_error .lt. abs(deform_pm(3,i, j, k) - deformation(3, i, j, k))) then
                     max_error = abs(deform_pm(3,i, j, k) - deformation(3, i, j, k))
                  end if
               end do
            end do
         end do
         !$omp enddo
         !$omp endparallel
         print *, "Deformation  Max error", total_error, "Max:", max_error
         print *, "Max val deformation      ", maxval(abs(deformation)), "Mean val deformation",  &
                                         sum(abs(deformation))/(NXpm_fine*NYpm_fine*NZpm_fine)
         print *, "MAX val deformation (old)", maxval(abs(deform_pm)), "Mean val deformation", &
                                     sum(abs(deform_pm))/(NXpm_fine*NYpm_fine*NZpm_fine)
      endif 

      ! Open 3 files for writing the deformation fields
      ! The first one will consist of the [x,y,z] components of the deformation field stored in the deformation array
      ! The second one will use the deform_pm arrays
      ! The third one will write the difference between the two deformation fields

      ! open(10, file='deformation_field.dat', status='unknown')
      ! open(11, file='deformation_field_pm.dat', status='unknown')
      ! open(12, file='deformation_field_diff.dat', status='unknown')

      ! do k = NZs_fine_bl + 1, NZf_fine_bl - 1
      !    do j = NYs_fine_bl + 1, NYf_fine_bl - 1
      !       do i = NXs_fine_bl + 1, NXf_fine_bl - 1
      !          write(10, *) i,j,k, deformation(1, i, j, k), deformation(2, i, j, k), deformation(3, i, j, k)
      !          write(11, *) i,j,k, deform_pm(1,i, j, k), deform_pm(2,i, j, k), deform_pm(3,i, j, k)
      !          write(12, *) i,j,k, deformation(1, i, j, k) - deform_pm(1,i, j, k), &
      !                       deformation(2, i, j, k) - deform_pm(2,i, j, k), &
      !                       deformation(3, i, j, k) - deform_pm(3,i, j, k)
      !       end do
      !    end do
      ! end do

      ! close(10)
      ! close(11)
      ! close(12)

      et = MPI_WTIME()
      write (dummy_string, "(A,I5,A,F8.2,A)") &
         achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
      call vpm_print(dummy_string, yellow, 1)
   end subroutine calc_velocity_serial_3d

   module subroutine diffuse_vort_3d
      use vpm_vars, only: OMPTHREADS, neqpm, NI
      use pmgrid, only: RHS_pm, deform_pm, &
                        DXpm, DYpm, DZpm,                            &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
      Implicit None
      real(dp) ::  dwxdx, dwydy, dwzdz, viscosity
      real(dp) ::  DXpm_sq, DYpm_sq, DZpm_sq 
      integer  :: i, j, k

      DXpm_sq = DXpm**2
      DYpm_sq = DYpm**2
      DZpm_sq = DZpm**2
      deform_pm = 0.d0
      !$omp parallel private(i,j,k,dwxdx,dwydy,dwzdz,viscosity) shared(RHS_pm,deform_pm)
      !!$omp num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 1, NZf_fine_bl - 1
         do j = NYs_fine_bl + 1, NYf_fine_bl - 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
               if (neqpm .eq. 3) then
                  viscosity = NI
               else
                  viscosity = RHS_pm(4, i, j, k) + NI
               end if
               !--> Remember that RHS = -w
               dwxdx = (RHS_pm(1, i + 1, j, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i - 1, j, k))/DXpm_sq
               dwydy = (RHS_pm(1, i, j + 1, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j - 1, k))/DYpm_sq
               dwzdz = (RHS_pm(1, i, j, k + 1) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(1,i, j, k) = -viscosity*(dwxdx + dwydy + dwzdz) ! because RHS=-w

               dwxdx = (RHS_pm(2, i + 1, j, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i - 1, j, k))/DXpm_sq
               dwydy = (RHS_pm(2, i, j + 1, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j - 1, k))/DYpm_sq
               dwzdz = (RHS_pm(2, i, j, k + 1) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(2,i, j, k) = -viscosity*(dwxdx + dwydy + dwzdz) ! because RHS=-w

               dwxdx = (RHS_pm(3, i + 1, j, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i - 1, j, k))/DXpm_sq
               dwydy = (RHS_pm(3, i, j + 1, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j - 1, k))/DYpm_sq
               dwzdz = (RHS_pm(3, i, j, k + 1) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(3, i, j, k) = -viscosity*(dwxdx + dwydy + dwzdz) ! because RHS=-w
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel

      !$omp parallel private(i,j,k) num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 1, NZf_fine_bl- 1
         do j = NYs_fine_bl + 1, NYf_fine_bl- 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
                 !--> Remember that RHS = -w
                 RHS_pm(1, i, j, k) = RHS_pm(1, i, j, k) - NI * deform_pm(1,i, j, k)
                 RHS_pm(2, i, j, k) = RHS_pm(2, i, j, k) - NI * deform_pm(2,i, j, k)
                 RHS_pm(3, i, j, k) = RHS_pm(3, i, j, k) - NI * deform_pm(3,i, j, k)
             enddo
         enddo
      enddo
      !$omp enddo
      !$omp endparallel
   end subroutine diffuse_vort_3d

   module subroutine calc_antidiffusion
      use pmgrid, only: RHS_pm, SOL_pm,                              &
                        deform_pm,          &
                        DXpm, DYpm, DZpm,                            &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        NXpm_fine, NYpm_fine, NZpm_fine
      Implicit None
      real(dp)                ::  dwxdx, dwydy, dwzdz, Ct
      real(dp)                ::  DXpm_sq, DYpm_sq, DZpm_sq
      integer                 :: i, j, k
      real(dp), allocatable   :: laplvort(:, :, :, :)

      allocate (laplvort(3, NXpm_fine, NYpm_fine, NZpm_fine))
      laplvort = 0.d0
      DXpm_sq = DXpm**2
      DYpm_sq = DYpm**2
      DZpm_sq = DZpm**2
      Sol_pm = 0.d0
      Ct = 6.8*DXpm**2/4
      do k = NZs_fine_bl + 1, NZf_fine_bl - 1
         do j = NYs_fine_bl + 1, NYf_fine_bl - 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
               !--> Remember that RHS = -w
               dwxdx = -(RHS_pm(1, i + 1, j, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i - 1, j, k))/DXpm_sq
               dwydy = -(RHS_pm(1, i, j + 1, k) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j - 1, k))/DYpm_sq
               dwzdz = -(RHS_pm(1, i, j, k + 1) - 2*RHS_pm(1, i, j, k) &
                        + RHS_pm(1, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               laplvort(1, i, j, k) = dwxdx + dwydy + dwzdz

               dwxdx = -(RHS_pm(2, i + 1, j, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i - 1, j, k))/DXpm_sq
               dwydy = -(RHS_pm(2, i, j + 1, k) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j - 1, k))/DYpm_sq
               dwzdz = -(RHS_pm(2, i, j, k + 1) - 2*RHS_pm(2, i, j, k) &
                        + RHS_pm(2, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               laplvort(2, i, j, k) = dwxdx + dwydy + dwzdz

               dwxdx = -(RHS_pm(3, i + 1, j, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i - 1, j, k))/DXpm_sq
               dwydy = -(RHS_pm(3, i, j + 1, k) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j - 1, k))/DYpm_sq
               dwzdz = -(RHS_pm(3, i, j, k + 1) - 2*RHS_pm(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               laplvort(3, i, j, k) = dwxdx + dwydy + dwzdz
            end do
         end do
      end do

      do k = NZs_fine_bl + 1, NZf_fine_bl - 1
         do j = NYs_fine_bl + 1, NYf_fine_bl - 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
               !Minus because of (-w) has been included in laplvort
               dwxdx = (laplvort(1, i + 1, j, k) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i - 1, j, k))/DXpm_sq
               dwydy = (laplvort(1, i, j + 1, k) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i, j - 1, k))/DYpm_sq
               dwzdz = (laplvort(1, i, j, k + 1) - 2*laplvort(1, i, j, k) &
                        + laplvort(1, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(1,i, j, k) = deform_pm(1,i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplvort(2, i + 1, j, k) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i - 1, j, k))/DXpm_sq
               dwydy = (laplvort(2, i, j + 1, k) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i, j - 1, k))/DYpm_sq
               dwzdz = (laplvort(2, i, j, k + 1) - 2*laplvort(2, i, j, k) &
                        + laplvort(2, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(2,i, j, k) = deform_pm(2,i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplvort(3, i + 1, j, k) - 2*laplvort(3, i, j, k) &
                        + laplvort(3, i - 1, j, k))/DXpm_sq
               dwydy = (laplvort(3, i, j + 1, k) - 2*laplvort(3, i, j, k) &
                        + laplvort(3, i, j - 1, k))/DYpm_sq
               dwzdz = (laplvort(3, i, j, k + 1) - 2*laplvort(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(3,i, j, k) = deform_pm(3,i, j, k) + Ct*dwxdx + dwydy + dwzdz
            end do
         end do
      end do

      deallocate(laplvort)
   end subroutine calc_antidiffusion

end module vpm_gcalc