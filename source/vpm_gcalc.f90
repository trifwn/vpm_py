module vpm_gcalc
   use base_types, only: dp
   use serial_vector_field_operators
   implicit none
contains
   subroutine calc_velocity
      use MPI
      use pmgrid, only: velocity_pm, DXpm, DYpm, DZpm, SOL_pm
      use vpm_vars, only: neqpm
      use console_io, only: vpm_print, blue, yellow, dummy_string
      use base_types, only: dp
      implicit none
      real(dp) ::  st, et
      st = MPI_WTIME()
      write (dummy_string, "(A)") "Calculating Velocities on PM using FD"
      call vpm_print(dummy_string, blue, 1)

      ! Velocity due to vorticity
      velocity_pm = curl(SOL_pm(1:3,:,:,:), DXpm, DYpm, DZpm)
      ! Velocity due to irrotational part
      if (neqpm .eq. 4) then
         velocity_pm = velocity_pm + gradient(SOL_pm(4,:,:,:), DXpm, DYpm, DZpm) 
      endif

      et = MPI_WTIME()
      write (dummy_string, "(A,I5,A,F8.2,A)") &
      achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
      call vpm_print(dummy_string, yellow, 1)
   end subroutine calc_velocity

   subroutine calc_vortex_stretching_conservative(velocity, deformation)
      use MPI
      use pmgrid, only: deform_pm, RHS_pm, DXpm, DYpm, DZpm,   &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        SOL_pm, NXpm_fine, NYpm_fine, NZpm_fine
      use vpm_vars, only: neqpm
      use console_io, only: vpm_print, blue, yellow, dummy_string
      use base_types, only: dp
      implicit none
      real(dp) ::  DXpm2, DYpm2, DZpm2
      real(dp) ::  st, et
      real(dp), intent(in) :: velocity(:,:,:,:)
      real(dp), allocatable, intent(out):: deformation(:,:,:,:)
      real(dp) ::  upi, upj, upk, vpi, vpj, vpk, wpi, wpj, wpk, &
                   umi, umj, umk, vmi, vmj, vmk, wmi, wmj, wmk, &
                   wdudx, wdvdy, wdwdz
      integer  :: i, j, k

      st = MPI_WTIME()
      write (dummy_string, "(A)") 'Calculating Deformations on PM using FD'
      call vpm_print(dummy_string, blue, 1)

      if (allocated(deformation)) deallocate(deformation)
      allocate(deformation, mold=velocity)

      DXpm2 = 2*DXpm
      DYpm2 = 2*DYpm
      DZpm2 = 2*DZpm
      deform_pm = 0.d0
      !$omp parallel private(upi,upj,upk,vpi,vpj,vpk,wpi,wpj,wpk,umi,umj,umk,vmi,vmj,vmk,&
      !$omp                  wmi,wmj,wmk,          &
      !$omp                  i,j,k,                &  
      !$omp                  wdudx,wdvdy,wdwdz )   & 
      !$omp shared(velocity,RHS_pm,SOL_pm,deform_pm, deformation, total_error, max_error)
      !$omp num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 2, NZf_fine_bl - 2
         do j = NYs_fine_bl + 2, NYf_fine_bl - 2
            do i = NXs_fine_bl + 2, NXf_fine_bl - 2
               !REMEMBER VORTICITY CARRIED IS -OMEGA and the quantity transfered is -OMEGA thus
               ! deformation = - (omega*\nabla) u =  \div(w (x) u) - u * (\grad w) = 
               ! We will calculate the deformation by using \div(w (x) u) as it is a conservative form
               upi = velocity(1, i + 1, j, k)
               upj = velocity(1, i, j + 1, k)
               upk = velocity(1, i, j, k + 1)

               vpi = velocity(2, i + 1, j, k)
               vpj = velocity(2, i, j + 1, k)
               vpk = velocity(2, i, j, k + 1)

               wpi = velocity(3, i + 1, j, k)
               wpj = velocity(3, i, j + 1, k)
               wpk = velocity(3, i, j, k + 1)

               umi = velocity(1, i - 1, j, k)
               umj = velocity(1, i, j - 1, k)
               umk = velocity(1, i, j, k - 1)

               vmi = velocity(2, i - 1, j, k)
               vmj = velocity(2, i, j - 1, k)
               vmk = velocity(2, i, j, k - 1)

               wmi = velocity(3, i - 1, j, k)
               wmj = velocity(3, i, j - 1, k)
               wmk = velocity(3, i, j, k - 1)
               !DEFORMATION WITH A MINUS BECAUSE WE HAVE STORED MINUS VORTICITY
               wdudx = -(RHS_pm(1, i + 1, j, k)*upi - (RHS_pm(1, i - 1, j, k))*umi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*upj - (RHS_pm(2, i, j - 1, k))*umj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*upk - (RHS_pm(3, i, j, k - 1))*umk)/DZpm2
               deform_pm(1,i, j, k) = wdudx + wdvdy + wdwdz

               wdudx = -(RHS_pm(1, i + 1, j, k)*vpi - (RHS_pm(1, i - 1, j, k))*vmi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*vpj - (RHS_pm(2, i, j - 1, k))*vmj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*vpk - (RHS_pm(3, i, j, k - 1))*vmk)/DZpm2
               deform_pm(2,i, j, k) = wdudx + wdvdy + wdwdz

               wdudx = -(RHS_pm(1, i + 1, j, k)*wpi - (RHS_pm(1, i - 1, j, k))*wmi)/DXpm2
               wdvdy = -(RHS_pm(2, i, j + 1, k)*wpj - (RHS_pm(2, i, j - 1, k))*wmj)/DYpm2
               wdwdz = -(RHS_pm(3, i, j, k + 1)*wpk - (RHS_pm(3, i, j, k - 1))*wmk)/DZpm2
               deform_pm(3,i, j, k) = wdudx + wdvdy + wdwdz
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel
      
      et = MPI_WTIME()
      write (dummy_string, "(A,I5,A,F8.2,A)") &
      achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
      call vpm_print(dummy_string, yellow, 1)
   end subroutine calc_vortex_stretching_conservative

   subroutine calc_vortex_stretching_definition(velocity, deformation)
      use MPI
      use pmgrid, only: RHS_pm, DXpm, DYpm, DZpm,              &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        SOL_pm, NXpm_fine, NYpm_fine, NZpm_fine
      use vpm_vars, only: neqpm
      use console_io, only: vpm_print, blue, yellow, dummy_string
      use base_types, only: dp
      implicit none
      real(dp) ::  DXpm2, DYpm2, DZpm2
      real(dp) ::  st, et
      real(dp), intent(in) :: velocity(:,:,:,:)
      real(dp), allocatable, intent(out):: deformation(:,:,:,:)
      real(dp), allocatable :: jac_u(:,:,:,:,:)
      real(dp) ::  upi, upj, upk, vpi, vpj, vpk, wpi, wpj, wpk, &
                   umi, umj, umk, vmi, vmj, vmk, wmi, wmj, wmk, &
                   wdudx, wdvdy, wdwdz
      integer  :: i, j, k
      if (allocated(deformation)) deallocate(deformation)
      allocate(deformation(3, NXpm_fine, NYpm_fine, NZpm_fine))
      ! We will calculate the deformation using the exact expression deformation = - (omega*\nabla) u
      
      jac_u = jacobian(velocity, DXpm, DYpm, DZpm)
      deformation = 0.d0
      !$omp parallel private(i,j,k ) 
      !$omp shared(velocity_pm,RHS_pm,SOL_pm, deformation)
      !$omp num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 2, NZf_fine_bl - 2
         do j = NYs_fine_bl + 2, NYf_fine_bl - 2
            do i = NXs_fine_bl + 2, NXf_fine_bl - 2
               deformation(1,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(1, :, i, j, k))
               deformation(2,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(2, :, i, j, k))
               deformation(3,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(3, :, i, j, k))
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel
      deallocate(jac_u)
   end subroutine calc_vortex_stretching_definition

   subroutine calc_vortex_stretching_conservative2(velocity, deformation)
      use MPI
      use pmgrid, only: deform_pm, RHS_pm, DXpm, DYpm, DZpm,   &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        SOL_pm, NXpm_fine, NYpm_fine, NZpm_fine
      use vpm_vars, only: neqpm
      use console_io, only: vpm_print, blue, yellow, dummy_string
      use base_types, only: dp
      implicit none
      real(dp) ::  DXpm2, DYpm2, DZpm2
      real(dp) ::  st, et
      real(dp), intent(in) :: velocity(:,:,:,:)
      real(dp), allocatable, intent(out):: deformation(:,:,:,:)
      real(dp), allocatable :: jac_u(:,:,:,:,:)
      real(dp), allocatable :: div_w(:,:,:)
      real(dp) ::  upi, upj, upk, vpi, vpj, vpk, wpi, wpj, wpk, &
                   umi, umj, umk, vmi, vmj, vmk, wmi, wmj, wmk, &
                   wdudx, wdvdy, wdwdz
      integer  :: i, j, k

      if (allocated(deformation)) deallocate(deformation)
      allocate(deformation(3, NXpm_fine, NYpm_fine, NZpm_fine))
      deformation = 0.d0

      ! We will calculate the deformation using the exact expression
      ! deformation = - (omega*\nabla) u = \div(w (x) u) - u * (\grad w) =
      ! We will calculate the deformation by using \div(w (x) u) as it is a conservative form
      ! \div(w (x) u) = \div(w) * u + w \div(u)
      jac_u = jacobian(velocity, DXpm, DYpm, DZpm)
      div_w = divergence(-RHS_PM(1:3,:,:,:), DXpm, DYpm, DZpm)
      
      !$omp parallel private(i,j,k ) 
      !$omp shared(velocity_pm,RHS_pm,SOL_pm, deformation)
      !$omp num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 2, NZf_fine_bl - 2
         do j = NYs_fine_bl + 2, NYf_fine_bl - 2
            do i = NXs_fine_bl + 2, NXf_fine_bl - 2
               deformation(1,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(1, :, i, j, k))
               deformation(2,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(2, :, i, j, k))
               deformation(3,i,j,k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(3, :, i, j, k))
   
               deformation(1:3,i,j,k) = deformation(1:3,i,j,k) + velocity(1:3, i, j, k) * div_w(i, j, k)
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel
      deallocate(jac_u)
      deallocate(div_w)
   end subroutine calc_vortex_stretching_conservative2

   module subroutine diffuse_vort_3d
      use MPI
      use vpm_vars, only: OMPTHREADS, neqpm, NI
      use pmgrid, only: RHS_pm, deform_pm, &
                        DXpm, DYpm, DZpm,                            &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
      implicit none
      real(dp) :: viscosity
      real(dp), allocatable :: laplace_vort(:,:,:,:)
      real(dp) ::  st, et, total_error
      integer  :: i, j, k
      laplace_vort = laplacian(RHS_pm(1:3,:,:,:), DXpm, DYpm, DZpm) 
      !$omp parallel private(i,j,k) num_threads(OMPTHREADS)
      !$omp do
      do k = NZs_fine_bl + 1, NZf_fine_bl- 1
         do j = NYs_fine_bl + 1, NYf_fine_bl- 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
                 !--> Remember that RHS = -w
               if (neqpm .eq. 3) then
                  viscosity = NI
               else
                  viscosity = RHS_pm(4, i, j, k) + NI
               end if
               laplace_vort(1:3, i, j, k) = viscosity*laplace_vort(1:3, i, j, k)
               RHS_pm(1, i, j, k) = RHS_pm(1, i, j, k) - laplace_vort(1, i, j, k) 
               RHS_pm(2, i, j, k) = RHS_pm(2, i, j, k) - laplace_vort(2, i, j, k) 
               RHS_pm(3, i, j, k) = RHS_pm(3, i, j, k) - laplace_vort(3, i, j, k) 
            enddo
         enddo
      enddo
      !$omp enddo
      !$omp endparallel
      deallocate(laplace_vort)
   end subroutine diffuse_vort_3d

   module subroutine calc_antidiffusion
      use pmgrid, only: RHS_pm, SOL_pm, deform_pm,             &
                        DXpm, DYpm, DZpm,                      &
                        NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl, &
                        NXpm_fine, NYpm_fine, NZpm_fine
      implicit none
      real(dp), allocatable   :: laplace_vort(:, :, :, :)
      real(dp)                :: dwxdx, dwydy, dwzdz, Ct
      real(dp)                :: DXpm_sq, DYpm_sq, DZpm_sq
      integer                 :: i, j, k

      DXpm_sq = DXpm**2
      DYpm_sq = DYpm**2
      DZpm_sq = DZpm**2
      Ct = 6.8*DXpm**2/4

      laplace_vort = -laplacian(RHS_pm(1:3,:,:,:), DXpm, DYpm, DZpm)
      !$omp parallel private(i,j,k, dwxdx, dwydy, dwzdz) num_threads(OMPTHREADS)
      !$omp shared(laplace_vort, deform_pm)
      !$omp do
      do k = NZs_fine_bl + 1, NZf_fine_bl - 1
         do j = NYs_fine_bl + 1, NYf_fine_bl - 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
               !Minus because of (-w) has been included in laplvort
               dwxdx = (laplace_vort(1, i + 1, j, k) - 2*laplace_vort(1, i, j, k) &
                        + laplace_vort(1, i - 1, j, k))/DXpm_sq
               dwydy = (laplace_vort(1, i, j + 1, k) - 2*laplace_vort(1, i, j, k) &
                        + laplace_vort(1, i, j - 1, k))/DYpm_sq
               dwzdz = (laplace_vort(1, i, j, k + 1) - 2*laplace_vort(1, i, j, k) &
                        + laplace_vort(1, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(1,i, j, k) = deform_pm(1,i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplace_vort(2, i + 1, j, k) - 2*laplace_vort(2, i, j, k) &
                        + laplace_vort(2, i - 1, j, k))/DXpm_sq
               dwydy = (laplace_vort(2, i, j + 1, k) - 2*laplace_vort(2, i, j, k) &
                        + laplace_vort(2, i, j - 1, k))/DYpm_sq
               dwzdz = (laplace_vort(2, i, j, k + 1) - 2*laplace_vort(2, i, j, k) &
                        + laplace_vort(2, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(2,i, j, k) = deform_pm(2,i, j, k) + Ct*dwxdx + dwydy + dwzdz

               dwxdx = (laplace_vort(3, i + 1, j, k) - 2*laplace_vort(3, i, j, k) &
                        + laplace_vort(3, i - 1, j, k))/DXpm_sq
               dwydy = (laplace_vort(3, i, j + 1, k) - 2*laplace_vort(3, i, j, k) &
                        + laplace_vort(3, i, j - 1, k))/DYpm_sq
               dwzdz = (laplace_vort(3, i, j, k + 1) - 2*laplace_vort(3, i, j, k) &
                        + RHS_pm(3, i, j, k - 1))/DZpm_sq
               ! U = grad x psi
               deform_pm(3,i, j, k) = deform_pm(3,i, j, k) + Ct*dwxdx + dwydy + dwzdz
            end do
         end do
      end do
      !$omp enddo
      !$omp endparallel
      deallocate(laplace_vort)
   end subroutine calc_antidiffusion

   subroutine write_deformation_fields(deform1, deform2)
      use pmgrid, only: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, &
                        NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
      implicit none
      real(dp), intent(in) :: deform1(:,:,:,:), deform2(:,:,:,:)
      integer :: i, j, k
      ! Writes 3 files
      ! 1) will consist of the [x,y,z] components of the field stored in deform1
      ! 2) will consist of the [x,y,z] components of the field stored in deform2
      ! 3) will consist of the absolute difference between the two deformation fields
      open(10, file='deformation_field.dat', status='unknown')
      open(11, file='deformation_field_pm.dat', status='unknown')
      open(12, file='deformation_field_diff.dat', status='unknown')
      do k = NZs_fine_bl + 1, NZf_fine_bl - 1
         do j = NYs_fine_bl + 1, NYf_fine_bl - 1
            do i = NXs_fine_bl + 1, NXf_fine_bl - 1
               write(10, *) i,j,k, deform1(1, i, j, k), deform1(2, i, j, k), deform1(3, i, j, k)
               write(11, *) i,j,k, deform2(1, i, j, k), deform2(2, i, j, k), deform2(3, i, j, k)
               write(12, *) i,j,k, abs(deform1(1, i, j, k) - deform2(1,i, j, k)), &
                                   abs(deform1(2, i, j, k) - deform2(2,i, j, k)), &
                                   abs(deform1(3, i, j, k) - deform2(3,i, j, k))
            end do
         end do
      end do
      close(10)
      close(11)
      close(12)
   end subroutine write_deformation_fields
end module vpm_gcalc