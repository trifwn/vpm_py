# VPM API

# Modules in project
- api
- fish
- pmmeshpar
- pmgrid
- parvar
- pmlib
- projlib
- vpm_vars
- vpm_size
- vpm_lib
- yapslib

# TODO:
- API TO READ AND WRITE COMMON BLOCKS MAINLY FROM MAIN_PM and VPM



# Functions in vpm.f90

- define_sizes
- project_particles
- diffuse_vort_3d
- rhsbcast
- rhsscat
- pmesh_solve
- write_pm_solution
- calc_velocity_serial_3d
- back_to_particles_3d
- back_to_particles_par
- pmesh
- yaps3d
- solget_3d
- velbcast_3d
- particles_scat
- particles_gath
- definepm


# TEST RECREATION
The main program is the test.f90

- Reads Particles.bin and sources.bin
```fortran
    open(1,file='particles.bin',form='unformatted') ! read particles
    read(1) NVR_ext
    read(1) Vref
    allocate(XPR(3,NVR_ext),QPR(neq+1,NVR_ext))
    QPR=0;XPR=0
    write(*,*) 'NVR=',NVR_ext,Vref
    do i=1,NVR_ext
       read(1) XPR(1,i),XPR(2,i),XPR(3,i),QPR(1,i),QPR(2,i),QPR(3,i)
    enddo

    QPR(1:3,:)   = -QPR(1:3,:)*Vref 
    QPR(4,:)     = QPR(4,:)*Vref
    QPR(neq+1,:) = Vref
    RMETM=0.001
    NVR_sources=0
   !open(1,file='sources.bin',form='unformatted')
   !read(1) NVR_sources
   allocate(XSOUR(3,NVR_sources),QSOUR(neq+1,NVR_sources))
   !write(*,*) 'NVR_sources=',NVR_sources
   !do i=1,NVR_sources
   !   read(1) XSOUR(1,i),XSOUR(2,i),XSOUR(3,i),QSOUR(1:4,i)
   !enddo
```

- Calls vpm with whattodo = 0

- Commented create_sources

- Calls remesh_particles_3d and 



# Modules:

```fortran
module pmeshpar
    double precision , save              :: PI, PI2, PI4,DT

    double precision, allocatable,save   :: velx_pm(:,:,:), vely_pm(:,:,:),velz_pm(:,:,:),qx_pm(:,:,:), qy_pm(:,:,:)
    double precision, allocatable,save   :: velphix_pm(:,:,:),velphiy_pm(:,:,:),velphiz_pm(:,:,:)
    double precision,allocatable, save   :: SOL_pm(:,:,:,:),SOL_0_pm(:,:,:,:)
    double precision,allocatable,save    :: source_bound(:,:),x_s(:,:),y_s(:,:),z_s(:,:),d_s(:),cos_s(:),sin_s(:)
    double precision,allocatable,save    :: source_bound_lev(:,:,:),xs_lev(:,:,:),ys_lev(:,:,:),zs_lev(:,:,:),ds_lev(:,:,:)
    integer         ,allocatable,save    :: nbound_lev(:)
    double precision,allocatable,save    :: Psiz_pm_0(:,:,:),Psiz_pm_f(:,:,:)
    double precision, allocatable,save   :: Cont_pm(:,:)
    integer                              :: levmax,npar_cell,ND

    integer, save                        :: nbound,ndumcell,NVR_CFD_sa,IDVPM

end module pmeshpar

module pmgrid
    double precision,allocatable,target,save        :: velvrx_pm(:,:,:),velvry_pm(:,:,:),velvrz_pm(:,:,:)
    double precision,allocatable, target,save       :: RHS_pm(:,:,:,:)

    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DD, DXpm, DYpm, DZpm,DXpm2,DYpm2,DZpm2,DVpm,EPSVOL
    integer ,save                        :: NXpm, NYpm, NZpm, NXpm_par,NYpm_par,NZpm_par    
    integer ,save                        :: NXs_bl(10),NYs_bl(10),NXf_bl(10),NYf_bl(10),NZs_bl(10),NZf_bl(10),NBlocks

end module pmgrid

module parvar
    integer                             ::NVR
    double precision, pointer , save    :: XP(:,:),QP(:,:)
    !     QP
    !   -->1 Vorticity X
    !   -->2 Vorticity Y
    !   -->3 Vorticity Z
    !   -->4 Dilatation
    !   -->5 Pseudopressure
    !   -->6 Mass
    !   -->7 Volume
    double precision, pointer , save    :: UP(:,:),GP(:,:)
    integer,allocatable,save            :: NVR_projtype(:)

    double precision,allocatable        :: XP_CFD_Sa(:,:),UT_CFD_Sa(:),UN_CFD_Sa(:),&
                                        DS_CFD_Sa(:)  ,Vo_CFD_Sa(:)
end module parvar

MODULE MKL_DFT_TYPE

end MODULE MKL_DFT_TYPE

module pmlib
    use mkl_dfti
    double precision , save              :: PI, PI2, PI4
    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DXpm, DYpm, DZpm,DXpm2,DYpm2,DZpm2

    integer ,save                        :: NVR, NXpm, NYpm, NZpm,ND
    integer ,save                        :: NXs_bl(10),NYs_bl(10),NXf_bl(10),NYf_bl(10),NZs_bl(10),NZf_bl(10),NBlocks
 
    integer, save                        :: nbound,levmax
!Here pointers are defined which will be assigned in the external data to save up space
    double precision,pointer             :: SOL_pm(:,:,:,:), RHS_pm(:,:,:,:),QP(:,:),XP(:,:)
    double precision,allocatable         :: SOL_0_pm(:,:,:,:), source_bound(:,:),x_s(:,:),y_s(:,:),z_s(:,:),d_s(:),cos_s(:),sin_s(:)
    double precision,allocatable,save    :: source_bound_lev(:,:,:),xs_lev(:,:),ys_lev(:,:),zs_lev(:,:),ds_lev(:,:)
    integer,allocatable,save             :: nbound_lev(:),ilev_t(:,:,:)

    private ::PI,PI2,PI4,XMIN_pm,XMAX_pm,YMIN_pm,YMAX_pm,ZMIN_pm,ZMAX_pm,DXpm,DYpm,DZpm,NVR,NXpm,NYpm,NZPm,ND
    private ::NXs_bl,NYs_bl,NXf_bl,NYf_bl,NZs_bl,NZf_bl,NBlocks,DXpm2,DYpm2,DZpm2
    private ::SOL_pm,RHS_pm,SOL_0_pm,QP,XP
    private ::source_bound,x_s,y_s,z_s,d_s,cos_s,sin_s
    private ::source_bound_lev,xs_lev,ys_lev,zs_lev,ds_lev
    private ::nbound,ilev_t

    contains 
    ...
end module pmlib

Module projlib

    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DXpm, DYpm, DZpm,DVpm
    double precision, save               :: EPSVOL
    integer ,save                        :: NXpm, NYpm, NZpm,NXs,NXf,NYs,NYf,NZs,NZf
    integer ,save                        :: IDVPM,ND



    private :: XMIN_pm,XMAX_pm,YMIN_pm,YMAX_pm,ZMIN_pm,ZMAX_pm,DXpm,DYpm,DZpm,NXpm,NYpm,NZpm,ND,EPSVOL,DVpm,IDVPM
    private :: NXs,NXf,NYs,NYf,NZs,NZf
end Module projlib

Module vpm_vars
   double precision, allocatable       :: XP_scatt(:, :), QP_scatt(:, :), UP_scatt(:, :), GP_scatt(:, :)
   double precision                    :: DT_c, V_ref, NI
   integer, allocatable                :: NVR_projscatt(:)
   integer                             :: interf_iproj, ncell_rem, iynhj

   integer                             :: ncoarse, nparcell1d
   integer                             :: neqpm, NVR_p, NVR_size, iwrite, NTIME_pm

   integer, save                       :: IPMWRITE, mrem, idefine, iynslice
   integer, save                       :: IPMWSTART(10), IPMWSTEPS(10)

End Module vpm_vars

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

End Module vpm_size

Module openmpth
   integer                       ::OMPTHREADS
End Module openmpth

Module vpm_lib
   use vpm_vars
   use vpm_size
   use openmpth

   !private ::  starttime,endtime,st,et,ct
   !private ::  nb_i,nb_j,nb_k,NBB,NXbl,NYbl,NZbl,BLOCKS,NXB,NYB,NZB,ndumcell_coarse ,ndumcell_bl
   !private :: II,iynbc,iret,NBI,NBJ,NBK,NVR_out_thres,NREMESH,ntorder,&
   !                                   iyntree,ilevmax,itree,nsize_out,ibctyp
end Module vpm_lib

!This library solves the poisson problem using domain decomposition method

Module yapslib
    use projlib
    use pmlib
    use MPI

    private

    double precision,allocatable      :: SOL_pm_coarse(:,:,:,:),RHS_pm_coarse(:,:,:,:),SOL_pm_sample(:,:,:,:,:)
    double precision,allocatable      :: SOL_pm_sumsample(:,:,:,:)

    double precision             :: DXpm,DYpm,DZpm,DXpm_c,DYpm_c,DZpm_c
    integer                      :: NXpm,NYpm,NZpm,NXpm_c,NYpm_c,NZPm_c
    integer                      :: ibctyp,idvpm,epsvol,ND,iproj,ndumcell,npmsize
    double precision             :: XMIN_pm,YMIN_pm,ZMIN_pm,XMAX_pm,YMAX_pm,ZMAX_pm
    double precision             :: XMIN_pm_c,YMIN_pm_c,ZMIN_pm_c,XMAX_pm_c,YMAX_pm_c,ZMAX_pm_c
    double precision             :: MACH
    double precision,allocatable :: Xbound_bl(:,:)
    integer                      :: BLOCKS,NXB,NYB,NBI,NBJ,NBK,NB,i,j,k,NXbl,NYbl,NN(3),NN_bl(6)
    integer,allocatable          :: NNbl(:,:),NNbl_bl(:,:),NN_coarse_map(:,:),map_nodes(:,:,:,:),nnb(:)
    double precision             :: projection_fun,fx,fy,fz,f,xc,yc,zc,X(3),addsol,starttime,endtime
    integer                      :: ic,jc,kc,inode,jnode,knode
    integer                      :: i_nb,j_nb,k_nb
    integer                      :: NXs,NYs,NZs,NXf,NYf,NZf,ib,jb,kb,ibj,jbj,kbj,ixs,ixf,jxs,jxf,izs,izf
    integer                      :: nc,NN_map(6),isubtrackt,node,neq,isizex,isizey,isizez,nbc


    double precision,pointer                 :: SOL_pm_bl(:,:,:,:), RHS_pm_bl(:,:,:,:),QP(:,:),XP(:,:)

    double precision,allocatable             :: BBound(:,:,:)


    integer  :: status(MPI_STATUS_SIZE),source,ierr,my_rank,np,mat4,mat5,dest
    character *25                :: outfil1,outfil2

    public  :: yaps2d,yaps3d
end Module yapslib

Module test_mod
  double precision,allocatable,target:: XPR(:,:),QPR(:,:),UPR(:,:),GPR(:,:),XPO(:,:),QPO(:,:),&
                                        XP_in(:,:),QP_in(:,:),XSOUR(:,:),QSOUR(:,:),XP_all(:,:),QP_all(:,:)
  double precision,pointer:: velx(:,:,:), vely(:,:,:), velz(:,:,:)
  double precision,pointer:: RHS_pm_in(:,:,:,:),RHS_pm_out(:,:,:,:)
  integer ,allocatable    :: qflag(:)
  integer :: NVR_ext,NVR_sources,NVR_sources_init,NVR_ext_init
end Module test_mod

```