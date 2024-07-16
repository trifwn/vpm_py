!    Subroutine writesolXavatar(itypesliceavatar)

!       ! AVATAR T2.4 specific: rotor + turbulent particles inflow [PM]

!       use vpm_vars
!       use pmeshpar
!       use parvar
!       use pmgrid
!       use MPI

!       character*50      :: filout
!       integer, intent(IN) :: itypesliceavatar
!       integer           :: i, j, k
!       double precision  :: XPM, YPM, ZPM, velocx, velocy, velocz, POSX(15)
!       integer, dimension(15) :: NX_AVA
!       integer               :: NXPOS_AVA_512_128_128, ii
!       integer               :: iynsliceavatar

!       if (iynsliceavatar .eq. 1) then
!          NXPOS_AVA_512_128_128 = 10
!          POSX(1) = -350
!          POSX(2) = -250
!          POSX(3) = -150
!          POSX(4) = -50
!          POSX(5) = 50
!          POSX(6) = 150
!          POSX(7) = 250
!          POSX(8) = 350
!          POSX(9) = 450

!          NX_AVA(1) = 2
!          NX_AVA(2:10) = int((POSX(1:9) - XMIN_pm)/DXpm) + 1
!       else
!          NXPOS_AVA_512_128_128 = 10
!          POSX(1) = -204
!          POSX(2) = -102
!          POSX(3) = -50
!          POSX(4) = -4
!          POSX(5) = 4
!          POSX(6) = 50
!          POSX(7) = 102
!          POSX(8) = 204
!          POSX(9) = 306

!          NX_AVA(1) = 2
!          NX_AVA(2:10) = int((POSX(1:9) - 7.1 - XMIN_pm)/DXpm) + 1
!       end if
!       write (filout, '(i5.5,a)') NTIME_pm, 'solX.bin'
!       open (1, file=filout, form='unformatted')

!       !    WRITE(1,'(a100)')'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'! "Vol" "PSIX" "PSIY" "PSIZ"'
!       do ii = 1, NXPOS_AVA_512_128_128
!          i = NX_AVA(ii) !/2  PM: DX=DY=DZ=8m
!          !WRITE(1,'(a11,i3,a8,i4,a7,i4,a7,i4,a11)') 'ZONE T= "',i,'", I=',1,&
!          !                                     ', J=',NYf_bl(1)-NYs_bl(1)+1, &
!          !                                     ', K=',NZf_bl(1)-NZs_bl(1)+1,  ', F=POINT'
!          write (1) NYf_bl(1) - NYs_bl(1) + 1, NZf_bl(1) - NZs_bl(1) + 1
!          do k = NZs_bl(1), NZf_bl(1)
!          do j = NYs_bl(1), NYf_bl(1)
!             XPM = XMIN_pm + (I - 1)*DXpm
!             YPM = YMIN_pm + (J - 1)*DYpm
!             ZPM = ZMIN_pm + (K - 1)*DZpm
!             velocx = VelvrX_pm(i, j, k)
!             velocy = VelvrY_pm(i, j, k)
!             velocz = VelvrZ_pm(i, j, k)
!             WRITE (1) XPM, YPM, ZPM, velocx, velocy, velocz, -RHS_pm(1, I, J, K), &
!                -RHS_pm(2, I, J, K), &
!                -RHS_pm(3, I, J, K)!,RHS_pm(4,I,J,K),SOL_pm(1,I,J,K),SOL_pm(2,I,J,K), SOL_pm(3,I,J,K)
!          end do !j
!          end do !k
!       end do !ii
!       close (1)

!       return

!    End Subroutine writesolXavatar

!    Subroutine writeline
!       use vpm_vars
!       use pmeshpar
!       use parvar
!       use pmgrid
!       use MPI
      
!       character*50        :: filout
!       integer           :: i, j, k, jmat(9), kmat(9), NNJ, NNK, il
!       double precision  :: XPM, YPM, ZPM, velocx, velocy, velocz
      
!       NNJ = NYf_bl(1) - NYs_bl(1) + 1
!       NNK = NZf_bl(1) - NZs_bl(1) + 1

!       kmat(1) = 0.5*NNK
!       kmat(2) = 0.25*NNK
!       kmat(3) = 0.75*NNK
!       kmat(4) = 0.25*NNK
!       kmat(5) = 0.75*NNK
!       kmat(6) = 0.25*NNK
!       kmat(7) = 0.75*NNK
!       kmat(8) = 0.5*NNK
!       kmat(9) = 0.5*NNK
      
!       jmat(1) = 0.5*NNJ;  
!       jmat(2) = 0.25*NNJ; 
!       jmat(3) = 0.25*NNJ; 
!       jmat(4) = 0.75*NNJ; 
!       jmat(5) = 0.75*NNJ; 
!       jmat(6) = 0.5*NNJ;  
!       jmat(7) = 0.5*NNJ;  
!       jmat(8) = 0.25*NNJ; 
!       jmat(9) = 0.75*NNJ; 
!       do i = 1, 9
!          j = jmat(i); k = kmat(i)
!          write (filout, '(i2.2,a)') i, 'hist.bin'
!          open (1, file=filout, access='APPEND', form='unformatted')
!          if (NTIME_pm .eq. 1) then
!             rewind (1)
!             write (1) XMIN_pm, YMIN_pm, ZMIN_pm
!             write (1) DXpm, DYpm, DZpm
!             write (1) NXs_bl(1), NXf_bl(1), j, k
!          end if
!          WRITE (1) NTIME_pm, (velvrx_pm(il, j, k), il=NXs_bl(1), NXf_bl(1)), &
!             (velvry_pm(il, j, k), il=NXs_bl(1), NXf_bl(1)), &
!             (velvrz_pm(il, j, k), il=NXs_bl(1), NXf_bl(1))

!          close (1)
!       end do

!    End Subroutine writeline