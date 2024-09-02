submodule(pmlib) pmsolve
   implicit none
contains

   !-----------------------------------------------------------------------!
   !-> subroutine solve_phiz                                                !
   !   This subroutines calls the fft library to solve for Phi poisson     !
   !   in all the points of Particle mesh.Dirichlet Boundary Cond. are used!
   !-----------------------------------------------------------------------!
   module subroutine solve_eq(NXs, NXf, NYs, NYf, neq)
      use MKL_POISSON

      Implicit None
      integer, intent(in)              :: NXs, NXf, NYs, NYf, neq
      Integer                          :: i, j, NX, NY
      integer                          :: ipar(128), stat
      integer                          :: NN, nod
      real(dp), allocatable            :: dpar(:)
      real(dp)                         :: XminCalc, XmaxCalc, YminCalc, YmaxCalc
      real(dp), allocatable            :: f(:), bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:)
      type(DFTI_DESCRIPTOR), pointer   :: xhandle

      !--> Assignment of Boundary Values
      ipar = 0
      XminCalc = XMIN_pm + (NXs - 1)*DXpm
      XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

      YminCalc = YMIN_pm + (NYs - 1)*DYpm
      YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

      NX = NXf - NXs + 1
      NY = NYf - NYs + 1
      !-->Set Right hand Side term (except for boundary conditions)

      NN = NX*NY
      allocate (f(NN))
      do j = 1, NY
         do i = 1, NX
            nod = (j - 1)*NX + i
            f(nod) = -RHS_pm(neq, NXs + i - 1, NYs + j - 1, 1)
         end do
      end do

      !-->Set Boundary Conditions
      !---> XMIN,XMAX

      !    Psiz_pm2(1:5,:)  = 0.d0
      !    Psiz_pm2(NX:NX-5,:) = 0.d0
      !    Psiz_pm2(:,1:5)  = 0.d0
      !    Psiz_pm2(:,NY:NY-5) = 0.d0

      NN = NY
      allocate (bd_ax(NN), bd_bx(NN))
      do j = 1, NY
         bd_ax(j) = SOL_pm(neq, NXs, j + NYs - 1, 1)
         bd_bx(j) = SOL_pm(neq, NXf, j + NYs - 1, 1)
      end do
      !---> YMIN,YMAX
      NN = NX
      allocate (bd_ay(NN), bd_by(NN))

      do i = 1, NX
         bd_ay(i) = SOL_pm(neq, i + NXs - 1, NYs, 1)
         bd_by(i) = SOL_pm(neq, i + NXs - 1, NYf, 1)
      end do

      !--SOLVE PHI ON PMESH
      allocate (dpar(int(5*(NX - 1)/2) + 7))
      call d_init_Helmholtz_2D(XminCalc, XmaxCalc, YminCalc, YmaxCalc, NX - 1, NY - 1, 'DDDD', 0.d0, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_init_Helmholtz_3D'
         stop
      endif

      call d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_commit_Helmholtz_3D'
         stop
      endif

      call d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_Helmholtz_3D'
         stop
      endif
      
      call free_Helmholtz_2D(xhandle, ipar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: free_Helmholtz_3D'
         stop
      endif

      do j = 1, NY
         do i = 1, NX
            nod = (j - 1)*NX + i
            SOL_pm(neq, NXs + i - 1, NYs + j - 1, 1) = f(nod)
         end do
      end do

   End subroutine solve_eq!_i

   module subroutine solve_eq_0(NXs, NXf, NYs, NYf, neq)
      use MKL_POISSON
      Implicit None
      integer, intent(in)     :: NXs, NXf, NYs, NYf, neq
      Integer                 :: i, j, NX, NY
      integer                 :: ipar(128), stat
      integer                 :: NN, nod
      real(dp), allocatable   :: dpar(:)
      real(dp)                :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc
      real(dp), allocatable   :: f(:), bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:)

      type(DFTI_DESCRIPTOR), pointer    :: xhandle

      !--> Assignment of Boundary Values
      ipar = 0
      XminCalc = XMIN_pm + (NXs - 1)*DXpm
      XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

      YminCalc = YMIN_pm + (NYs - 1)*DYpm
      YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

      NX = NXf - NXs + 1
      NY = NYf - NYs + 1
      !-->Set Right hand Side term (except for boundary conditions)

      NN = NX*NY
      allocate (f(NN))
      do j = 1, NY
         do i = 1, NX
            nod = (j - 1)*NX + i
            f(nod) = -RHS_pm(neq, NXs + i - 1, NYs + j - 1, 1)
         end do
      end do

      !-->Set Boundary Conditions
      !---> XMIN,XMAX

      !    Psiz_pm2(1:5,:)  = 0.d0
      !    Psiz_pm2(NX:NX-5,:) = 0.d0
      !    Psiz_pm2(:,1:5)  = 0.d0
      !    Psiz_pm2(:,NY:NY-5) = 0.d0

      NN = NY
      allocate (bd_ax(NN), bd_bx(NN))
      bd_ax = 0.d0; bd_bx = 0.d0; !  SOL_pm(NXs,j + NYs -1,nb,neq)
      !---> YMIN,YMAX
      NN = NX
      allocate (bd_ay(NN), bd_by(NN))
      bd_ay = 0.d0; bd_by = 0.d0! SOL_pm(i + NXs-1,NYf ,nb,neq)

      !--SOLVE PHI ON PMESH
      allocate (dpar(int(5*(NX - 1)/2) + 7))
      call d_init_Helmholtz_2D(XminCalc, XmaxCalc, YminCalc, YmaxCalc, NX - 1, NY - 1, 'DDDD', 0.d0, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_init_Helmholtz_3D'
         stop
      endif
      
      call d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_commit_Helmholtz_3D'
         stop
      endif

      call d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_Helmholtz_3D'
         stop
      endif

      call free_Helmholtz_2D(xhandle, ipar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: free_Helmholtz_3D'
         stop
      endif

      !call hwscrt(XminCalc,XmaxCalc,NX-1,1,0,0,YminCalc,YmaxCalc,NY-1,1,0,0,  &
      !                   0,SOL_pm2,NX,pertrb,INFO,WORK)

      !  SOL_pm(NXs:NXf,NYs:NYf,nb,neq)=SOL_pm2(1:NX,1:NY)
      do j = 1, NY
         do i = 1, NX
            nod = (j - 1)*NX + i
            SOL_0_pm(neq, NXs + i - 1, NYs + j - 1, 1) = f(nod)
         end do
      end do

   End subroutine solve_eq_0

   module subroutine solve_eq_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
      use MKL_POISSON
      Implicit None
      integer, intent(in) :: NXs, NXf, NYs, NYf, NZs, NZf, neq
      integer            :: i, j, k, NX, NY, NZ
      real(dp)   :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, ZminCalc, ZmaxCalc
      integer              :: ipar(128), stat
      integer              :: dim, dim1, dim2, dim3, nod
      real(dp), allocatable ::dpar(:)
      real(dp), allocatable, dimension(:):: f
      real(dp), allocatable, dimension(:):: bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz
      type(DFTI_DESCRIPTOR), pointer   :: xhandle, yhandle

      !--> Assignment of Boundary Values
      ipar = 0
      XminCalc = XMIN_pm + (NXs - 1)*DXpm
      XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

      YminCalc = YMIN_pm + (NYs - 1)*DYpm
      YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

      ZminCalc = ZMIN_pm + (NZs - 1)*DZpm
      ZmaxCalc = ZMIN_pm + (NZf - 1)*DZpm

      NX = NXf - NXs + 1!-2
      NY = NYf - NYs + 1!-2
      NZ = NZf - NZs + 1!-2

      dim = NX*NY*NZ
      dim1 = NY*NZ
      dim2 = NX*NY
      dim3 = NX*NZ

      allocate (bd_ax(dim1), bd_bx(dim1))
      allocate (bd_az(dim2), bd_bz(dim2))
      allocate (bd_ay(dim3), bd_by(dim3))
      allocate (f(dim))

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               nod = (k - 1)*NX*NY + (j - 1)*NX + i
               f(nod) = -RHS_pm(neq, NXs + i - 1, NYs + j - 1, NZs + k - 1)
            end do
         end do
      end do

      do k = 1, NZ
         do j = 1, NY
            nod = (k - 1)*NY + j
            bd_ax(nod) = SOL_pm(neq, NXs, j + NYs - 1, k + NZs - 1)
            bd_bx(nod) = SOL_pm(neq, NXf, j + NYs - 1, k + NZs - 1)
         end do
      end do

      do j = 1, NY
         do i = 1, NX
            nod = (j - 1)*NX + i
            bd_az(nod) = SOL_pm(neq, i + NXs - 1, j + NYs - 1, NZs)
            bd_bz(nod) = SOL_pm(neq, i + NXs - 1, j + NYs - 1, NZf)
         end do
      end do

      !---> YMIN,YMAX
      do k = 1, NZ
         do i = 1, NX
            nod = (k - 1)*NX + i
            bd_ay(nod) = SOL_pm(neq, i + NXs - 1, NYs, k + NZs - 1)
            bd_by(nod) = SOL_pm(neq, i + NXs - 1, NYf, k + NZs - 1)
         end do
      end do

      allocate (dpar(int(5*(NX - 1 + NY - 1)/2) + 9))
      call d_init_Helmholtz_3D(XminCalc,XmaxCalc,YminCalc,YmaxCalc,ZminCalc,ZmaxCalc,NX-1,NY-1,NZ-1,'DDDDDD',0.d0,ipar,dpar,stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_init_Helmholtz_3D'
         stop
      endif
      
      call d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle, yhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_commit_Helmholtz_3D'
         stop
      endif

      call d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle, yhandle, ipar, dpar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: d_Helmholtz_3D'
         stop
      endif

      call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
      if (stat .ne. 0) then
         write(*,*) 'Error in solve_eq_0_3d: free_Helmholtz_3D'
         stop
      endif


      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               nod = (k - 1)*NX*NY + (j - 1)*NX + i
               SOL_pm(neq, NXs + i - 1, NYs + j - 1, NZs + k - 1) = f(nod)
            end do
         end do
      end do

   End subroutine solve_eq_3d

   module subroutine solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
      use MKL_POISSON

      Implicit None
      integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
      integer              :: i, j, k, NX, NY, NZ
      real(dp)     :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, ZminCalc, ZmaxCalc
      integer              :: ipar(128), stat
      integer              :: NN, nod, NNX, NNY, NNZ
      ! real(dp), allocatable:: f(:), bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:), bd_az(:), bd_bz(:)
      real(dp)     :: dpar((5*(NXf - NXs + NYf - NYs)/2) + 9)
      real(dp)     ::  f((NXf - NXs + 1)*(NYf - NYs + 1)*(NZf - NZs + 1))
      real(dp)     ::  bd_ax((NYf - NYs + 1)*(NZf - NZs + 1)), bd_bx((NYf - NYs + 1)*(NZf - NZs + 1))
      real(dp)     ::  bd_ay((NXf - NXs + 1)*(NZf - NZs + 1)), bd_by((NXf - NXs + 1)*(NZf - NZs + 1))
      real(dp)     ::  bd_az((NXf - NXs + 1)*(NYf - NYs + 1)), bd_bz((NXf - NXs + 1)*(NYf - NYs + 1))
      type(DFTI_DESCRIPTOR), pointer    :: xhandle, yhandle
      character(6) BCtype


      ! Setting the type of the boundary conditions on each surface of the parallelepiped domain:
      ! On the boundary laying on the plane x=0(=ax) Dirichlet boundary condition will be used
      ! On the boundary laying on the plane x=1(=bx) Dirichlet boundary condition will be used
      ! On the boundary laying on the plane y=0(=ay) Dirichlet boundary condition will be used
      ! On the boundary laying on the plane y=1(=by) Dirichlet boundary condition will be used
      ! On the boundary laying on the plane z=0(=az) Dirichlet boundary condition will be used
      ! On the boundary laying on the plane z=1(=bz) Dirichlet boundary condition will be used
      BCtype = 'DDDDDD'

      !--> Assignment of Boundary Values
      ipar = 0
      XminCalc = XMIN_pm + (NXs - 1)*DXpm
      XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

      YminCalc = YMIN_pm + (NYs - 1)*DYpm
      YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

      ZminCalc = ZMIN_pm + (NZs - 1)*DZpm
      ZmaxCalc = ZMIN_pm + (NZf - 1)*DZpm

      NX = NXf - NXs + 1!-2
      NY = NYf - NYs + 1!-2
      NZ = NZf - NZs + 1!-2

      NN = NX*NY*NZ
      NNX = NY*NZ
      NNY = NX*NZ
      NNZ = NX*NY

      ! DYNAMIC Allocation produced errors on certain systems. Current impl is reverted to static allocation
      ! allocate (f(NN))
      ! allocate(bd_ax(NNX))
      ! allocate(bd_bx(NNX))
      ! allocate(bd_ay(NNY))
      ! allocate(bd_by(NNY))
      ! allocate(bd_az(NNZ))
      ! allocate(bd_bz(NNZ))

      !-->Set Boundary Conditions
      !---> XMIN,XMAX
      bd_ax = 0
      bd_bx = 0
      
      !---> YMIN,YMAX
      bd_ay = 0
      bd_by = 0 

      !---> ZMIN,ZMAX
      bd_az = 0
      bd_bz = 0; 

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               nod = (k - 1)*NX*NY + (j - 1)*NX + i
               f(nod) = -RHS_pm(neq, NXs + i - 1, NYs + j - 1, NZs + k - 1)
               ! CHECK F BOUNDS
               if (nod .gt. NN) then
                  stop
               end if
               ! CHECK RHS BOUNDS
               if (                                   &
                  (NXS+i-1 .gt. size(RHS_pm,2)) .or.  &
                  (NYS+j-1 .gt. size(RHS_pm,3)) .or.  &
                  (NZS+k-1 .gt. size(RHS_pm,4))       &
               ) then
                  stop
               end if
            end do
         end do
      end do
      
      ! Initializing ipar array to make it free from garbage
      do i=1,128
         ipar(i)=0
      enddo

      ! Computing the approximate solution of 3D Laplace problem
      ! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz should not be changed
      ! between the Commit step and the subsequent call to the Solver routine!
      ! Otherwise the results may be wrong.
      call d_init_Helmholtz_3D(&
         XminCalc,XmaxCalc,YminCalc,YmaxCalc,ZminCalc,ZmaxCalc,NX-1,NY-1,NZ-1,BCtype,0.d0,ipar,dpar,stat &
      )
      if (stat .ne. 0) then 
         write(*,*) 'Error in solve_eq_0_3d: d_init_Helmholtz_3D'
         stop
      enendif

      ! Initializing complex data structures of Poisson Library for 3D Laplace Solver
      ! NOTE: Right-hand side f may be altered after the Commit step. If you want to keep it,
      ! you should save it in another memory location!
      call d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle, yhandle, ipar, dpar, stat)
      if (stat .ne. 0) then 
         write(*,*) 'Error in solve_eq_0_3d: d_commit_Helmholtz_3D'
         stop
      endif

      ! Computing the approximate solution of 3D Laplace problem
      ! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz should not be changed
      ! between the Commit step and the subsequent call to the Solver routine!
      ! Otherwise the results may be wrong.
      call d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle, yhandle, ipar, dpar, stat)
      if (stat .ne. 0) then 
         write(*,*) 'Error in solve_eq_0_3d: d_Helmholtz_3D'
         stop
      endif

      ! Cleaning the memory used by xhandle and yhandle
      call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
      if (stat .ne. 0) then 
         write(*,*) 'Error in solve_eq_0_3d: free_Helmholtz_3D'
         stop
      endif

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               nod = (k - 1)*NX*NY + (j - 1)*NX + i
               SOL_0_pm(neq, NXs + i - 1, NYs + j - 1, NZs + k - 1) = f(nod)
            end do
         end do
      end do
   End subroutine solve_eq_0_3d
End Submodule pmsolve