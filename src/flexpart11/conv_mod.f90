!*******************************************************************************
!   Include file for convection
!   This file contains a global common block used by convect
!   and other subroutines
!   Author: P. Ferstl
!
!   Feb 2001
!
!   Changes
!       2021 L. Bakels:
!          - Array operations in convect subroutine
!          - OpenMP parallelisation in convmix and redist
!          - Moved all subroutines related to the convection to this module
!*******************************************************************************

module conv_mod

  use com_mod, only: lconvection, numthreads,nxmaxn, nymaxn,numbnests
  use windfields_mod, only: metdata_format,akz,bkz,akm,bkm,nuvz, &
    uvheight,ps,tt2,td2,tth,qvh,pplev,tt,qv,nx,ny,tt2n,td2n,psn,tthn,qvhn, &
    yln,yrn,xln,xrn,yresoln,xresoln,nxn,nyn,dxn,dyn, &
    nxmax,nymax,nconvlevmax,na,nuvzmax
  use sort_mod
  implicit none

  !integer,parameter :: nconvlevmax = nuvzmax-1, &
  !                     na = nconvlevmax+1
  !these parameters are defined in par_mod now!
  real,allocatable,dimension(:,:) :: &
    pconv,                           & ! 
    phconv,                          & !
    dpr,                             & !
    pconv_hpa,                       & !
    phconv_hpa,                      & !
    sub,                             & ! subsidence
    tconv,                           & !
    qconv,                           & !
    qsconv                            
  real,allocatable,dimension(:,:,:) :: & !
    fmass,                           & !
    fmassfrac,                       & !
    cbaseflux                         
  real,allocatable,dimension(:,:,:,:) :: &
    cbasefluxn
  ! integer,dimension(na) :: &
  !   NENT
  ! real,dimension(na,na) :: &
  !   MENT,QENT,ELIJ,SIJ
  ! real,dimension(na) ::    &
  !   fup,fdown,M,MP,TVP,TV, &
  !   WATER,QP,EP,TH,WT,     &
  !   EVAP,CLW,SIGP,TP,CPN,  &
  !   LV,LVCP,H,HP,GZ,HM
  real,allocatable,dimension(:,:) ::    &
    uvzlev,wsub
  real :: psconv,tt2conv,td2conv
  
  integer :: nconvlev ! global variable set at start in FLEXPART.f90
  integer :: nconvtop ! needs to be private, set in convect

  save :: uvzlev

! $OMP THREADPRIVATE( fmass, sub, fmassfrac, &
! $OMP pconv, phconv, dpr, pconv_hpa, phconv_hpa, &
! $OMP tconv, qconv, qsconv, psconv, tt2conv, td2conv, &
! $OMP nconvtop,uvzlev,wsub,cbaseflux,cbasefluxn)

!$OMP THREADPRIVATE( psconv, tt2conv, td2conv,nconvtop)
! , &
! !$OMP fup,fdown,MENT,NENT,M,MP,QENT,ELIJ,SIJ,TVP,TV, &
! !$OMP WATER,QP,EP,TH,WT,EVAP,CLW,SIGP,TP,CPN,LV,LVCP, &
! !$OMP H,HP,GZ,HM)

contains

subroutine alloc_convect
  implicit none
  integer :: stat
  if (.not.lconvection.eq.1) return
  ! nconvlevmax=nuvzmax-1
  ! na=nconvlevmax+1
  allocate( pconv(nconvlevmax,numthreads),phconv(na,numthreads), &
    dpr(nconvlevmax,numthreads), &
    pconv_hpa(nconvlevmax,numthreads),phconv_hpa(na,numthreads), &
    fmass(nconvlevmax,nconvlevmax,numthreads),        &
    sub(nconvlevmax,numthreads), &
    fmassfrac(nconvlevmax,nconvlevmax,numthreads),   &
    cbaseflux(0:nxmax-1,0:nymax-1,numthreads),                        &
    cbasefluxn(0:nxmaxn-1,0:nymaxn-1,numbnests,numthreads),            &
    tconv(na,numthreads),qconv(na,numthreads),qsconv(na,numthreads),stat=stat)
  if (stat.ne.0) error stop "Could not allocate convection arrays"

  allocate( uvzlev(nuvzmax,numthreads),wsub(nuvzmax,numthreads),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uvzlev or wsub"

  ! allocate(FUP(NA),FT(NA),FQ(NA),FDOWN(NA),NENT(NA),       &
  !   M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA),      &
  !   SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA),                   &
  !   QP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA),          &
  !   SIGP(NA),TP(NA),CPN(NA),                               &
  !   LV(NA),LVCP(NA),H(NA),HP(NA),GZ(NA),HM(NA))

end subroutine alloc_convect

subroutine dealloc_convect
  implicit none
  if (.not.lconvection.eq.1) return
  deallocate(pconv,phconv,dpr,pconv_hpa,phconv_hpa,sub,      &
    tconv,qconv,qsconv,fmass,fmassfrac,cbaseflux,cbasefluxn)
  deallocate(uvzlev,wsub)
  ! deallocate(fup,fdown,ment,ft,fq,M,MP,QENT,ELIJ,SIJ,TVP,TV, &
  !   WATER,QP,EP,TH,WT,EVAP,CLW,SIGP,TP,CPN,LV,LVCP, &
  !   H,HP,GZ,HM)
end subroutine dealloc_convect

subroutine set_conv_top()
! Determine the uppermost level for which the convection scheme shall be applied
! by assuming that there is no convection above 50 hPa (for standard SLP)
!*****************************************************************************  

  integer :: i
  real :: pint

  do i=1,nuvz-2
    pint=akz(i)+bkz(i)*101325.
    if (pint.lt.5000.) exit
  end do
  nconvlev=i
  if (nconvlev.gt.nconvlevmax-1) then
    nconvlev=nconvlevmax-1
    write(*,*) 'INFORMATION: Convection only calculated up to ', &
         akz(nconvlev)+bkz(nconvlev)*1013.25,' hPa'
  endif  

end subroutine set_conv_top

subroutine convmix(itime)
  !                     i
  !**************************************************************
  !handles all the calculations related to convective mixing
  !Petra Seibert, Bernd C. Krueger, Feb 2001
  !nested grids included, Bernd C. Krueger, May 2001
  !
  !Changes by Caroline Forster, April 2004 - February 2005:
  !  convmix called every lsynctime seconds
  !CHANGES by A. Stohl:
  !  various run-time optimizations - February 2005
  ! CHANGES by C. Forster, November 2005, NCEP GFS version
  !      in the ECMWF version convection is calculated on the
  !      original eta-levels
  !      in the GFS version convection is calculated on the
  !      FLEXPART levels
  !
  !   Unified ECMWF and GFS builds                                             
  !   Marian Harustak, 12.5.2017                                              
  !     - Merged convmix and convmix_gfs into one routine using if-then           
  !       for meteo-type dependent code                                        
  !**************************************************************
  use omp_lib
  use flux_mod
  use par_mod
  use com_mod
  use class_gribfile_mod
  use particle_mod

  implicit none

  integer :: igr,igrold, ipart, itime, ix, i, j, ik, inest
  integer :: ipconv,ithread,stat,countconv
  integer :: jy, kpart, ktop, ngrid,kz
  integer,allocatable :: igrid(:), ipoint(:), igridn(:,:)

  ! itime [s]                 current time
  ! igrid(maxpart)            horizontal grid position of each particle
  ! igridn(maxpart,maxnests)  dto. for nested grids
  ! ipoint(maxpart)           pointer to access particles according to grid position

  logical :: lconv,lcalcflux
  real :: x, y, xtn,ytn, ztold, delt
  real :: dt1,dt2,dtt
  integer :: mind1,mind2
  ! dt1,dt2,dtt,mind1,mind2       variables used for time interpolation
  integer :: itage,nage,inage

  ! OMP changes
  integer :: cnt,kk
  integer,allocatable,dimension(:) :: frst,kkcnt
  double precision :: tmarray(2)

  integer :: alivepart
  real:: eps
  eps=nxmax/3.e5
  ! Calculate auxiliary variables for time interpolation
  !*****************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)
  mind1=memind(1)
  mind2=memind(2)
  delt=real(abs(lsynctime))

  lconv = .false.

  ! if no particles are present return after initialization
  !********************************************************
  call get_alivepart_num(alivepart)
  if (alivepart.le.0 ) return

  !call get_totalpart_num(totpart)
  allocate( igrid(alivepart),stat=stat)
  if (stat.ne.0) error stop "Could not allocate igrid"
  allocate( ipoint(alivepart),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ipoint"
  allocate( igridn(alivepart,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate igridn"

  ! Assign igrid and igridn, which are pseudo grid numbers indicating particles
  ! that are outside the part of the grid under consideration
  ! (e.g. particles near the poles or particles in other nests).
  ! Do this for all nests but use only the innermost nest; for all others
  ! igrid shall be -1
  ! Also, initialize index vector ipoint
  !************************************************************************
  igrid(:) = -1
  do j=numbnests,1,-1
    igridn(:,j)=-1
  end do

!$OMP PARALLEL PRIVATE(ipart, i, j, x, y, ngrid, xtn, ytn, ix, jy)
!$OMP DO
  do i=1,alivepart
    ! do not consider particles that are not (yet) part of simulation
    ipart=count%ialive(i)

    ipoint(i)=ipart

    x = real(part(ipart)%xlon)
    y = real(part(ipart)%ylat)

  ! Determine which nesting level to be used
  !**********************************************************

    ngrid=0
    if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
      do j=numbnests,1,-1
        ! Temporary fix for nested layer edges: replaced eps with dxn and dyn (LB)
        if ( x.gt.xln(j)+dxn(j) .and. x.lt.xrn(j)-dxn(j) .and. &
             y.gt.yln(j)+dyn(j) .and. y.lt.yrn(j)-dyn(j) ) then
          ngrid=j
          exit
        endif
      end do
    else
      do j=numbnests,1,-1
        if ( x.gt.xln(j) .and. x.lt.xrn(j) .and. &
             y.gt.yln(j) .and. y.lt.yrn(j) ) then
          ngrid=j
          exit
        endif
      end do
    endif
 ! 23   continue

  ! Determine nested grid coordinates
  !**********************************

    if (ngrid.gt.0) then
  ! nested grids
      xtn=(x-xln(ngrid))*xresoln(ngrid)
      ytn=(y-yln(ngrid))*yresoln(ngrid)
      ix=nint(xtn)
      jy=nint(ytn)
      ! igridn(ipart,ngrid) = 1 + jy*nxn(ngrid) + ix
      igridn(i,ngrid) = 1 + ix*nyn(ngrid) + jy
    else if(ngrid.eq.0) then
  ! mother grid
      ix=nint(x)
      jy=nint(y)
      !igrid(ipart) = 1 + jy*nx + ix
      igrid(i) = 1 + ix*ny + jy
    endif
  end do
!$OMP END DO
!$OMP END PARALLEL

  ! sumall = 0.
  ! sumconv = 0.

  !*****************************************************************************
  ! 1. Now, do everything for the mother domain and, later, for all of the nested domains
  ! While all particles have to be considered for redistribution, the Emanuel convection
  ! scheme only needs to be called once for every grid column where particles are present.
  ! Therefore, particles are sorted according to their grid position. Whenever a new grid
  ! cell is encountered by looping through the sorted particles, the convection scheme is called.
  !*****************************************************************************

  ! sort particles according to horizontal position and calculate index vector IPOINT

  call sort2(alivepart,igrid,ipoint)

  ! Now visit all grid columns where particles are present
  ! by going through the sorted particles

  !LB changes following the CTM version
  allocate( frst(nx*(ny+1)+1),stat=stat)
  if (stat.ne.0) error stop "Could not allocate frst"
  frst(1) = 1
  cnt = 2
  igrold = igrid(1)
  ! Looping over all particles and counting how many in each igrid reside.
  ! This is saved in frst. The number of consecutive particles in igrid is saved in frst(i)
  do kpart=1,alivepart
    if (igrold.ne.igrid(kpart)) then
      frst(cnt) = kpart
      igrold=igrid(kpart)
      cnt=cnt+1
    endif
  end do 
  frst(cnt) = alivepart+1

  allocate(kkcnt(cnt-1))
  countconv=0
  do kk=1,cnt-1
    if (igrid(frst(kk)).eq.-1) cycle ! Only consider grids that have particles inside
    countconv=countconv+1
    kkcnt(countconv)=kk
  end do

!$OMP PARALLEL PRIVATE(kk,jy,ix,tmarray,j,kz,ktop,lconv,kpart,ipart,&
!$OMP ztold,nage,ipconv,itage,ithread,lcalcflux)

#if (defined _OPENMP)
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts at 1
#else
    ithread = 1
#endif

!$OMP DO SCHEDULE(dynamic)
  do ik=1,countconv
    
    !if (igrid(frst(kk)).eq.-1) cycle
    kk=kkcnt(ik)
    ! Find horizontal location of grid column
    ix = (igrid(frst(kk))-1)/ny
    jy = igrid(frst(kk)) - ix*ny - 1
    ! jy = (igrid(frst(kk))-1)/nx
    ! ix = igrid(frst(kk)) - jy*nx - 1

  ! Interpolate all meteorological data needed for the convection scheme
    psconv=(ps(ix,jy,1,mind1)*dt2+ps(ix,jy,1,mind2)*dt1)*dtt
    tt2conv=(tt2(ix,jy,1,mind1)*dt2+tt2(ix,jy,1,mind2)*dt1)*dtt
    td2conv=(td2(ix,jy,1,mind1)*dt2+td2(ix,jy,1,mind2)*dt1)*dtt

    if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
      do kz=1,nuvz-1           !bugfix
        tconv(kz,ithread)=(tth(ix,jy,kz+1,mind1)*dt2+ &
             tth(ix,jy,kz+1,mind2)*dt1)*dtt
        qconv(kz,ithread)=(qvh(ix,jy,kz+1,mind1)*dt2+ &
             qvh(ix,jy,kz+1,mind2)*dt1)*dtt
      end do
    else
      do kz=1,nuvz-1           !bugfix
        pconv(kz,ithread)=(pplev(ix,jy,kz,mind1)*dt2+ &
            pplev(ix,jy,kz,mind2)*dt1)*dtt
        tconv(kz,ithread)=(tt(ix,jy,kz,mind1)*dt2+ &
            tt(ix,jy,kz,mind2)*dt1)*dtt
        qconv(kz,ithread)=(qv(ix,jy,kz,mind1)*dt2+ &
            qv(ix,jy,kz,mind2)*dt1)*dtt
      end do
    end if

  ! Calculate translocation matrix
    call calcmatrix(lconv,delt,cbaseflux(ix,jy,ithread),ithread)
    
  ! treat particle only if column has convection
    if (lconv .eqv. .true.) then
      ktop = 0
  ! assign new vertical position to particle
      do kpart=frst(kk), frst(kk+1)-1
        ipart = ipoint(kpart)
        ztold=real(part(ipart)%z)
        call redist(itime,ipart,ktop,ipconv,ithread)
  !    if (ipconv.le.0) sumconv = sumconv+1

  ! Calculate the gross fluxes across layer interfaces
  !***************************************************

        if (iflux.eq.1) then
          lcalcflux=.true.
          if (lagespectra.eq.1) then
            nage=0
            itage=abs(itime-part(ipart)%tstart)
            do inage=1,nageclass
              if ((itage.lt.lage(inage)).and.(part(ipart)%alive)) exit
              nage=inage
            end do
            if (nage.eq.nageclass) lcalcflux=.false.
            nage=nage+1
          else
            nage=1
          endif

          if (lcalcflux) &
            call calcfluxes(itime,nage,ipart,real(part(ipart)%xlon), &
               real(part(ipart)%ylat),ztold,ithread)
        endif
      enddo

    endif   !(lconv .eqv. .true)
  end do
!$OMP END DO
!$OMP END PARALLEL

  deallocate(frst)

  ! OpenMP Reduction for dynamically allocated arrays. This is done manually since this
  ! is not yet supported in most OpenMP versions
  !************************************************************************************
  if (iflux.eq.1) then
    do ithread=1,numthreads
      flux(:,:,:,:,:,:,:)=flux(:,:,:,:,:,:,:)+flux_omp(:,:,:,:,:,:,:,ithread)
    end do
  endif

  !*****************************************************************************
  ! 2. Nested domains
  !*****************************************************************************

  ! sort particles according to horizontal position and calculate index vector IPOINT
  do inest=1,numbnests
    do i=1,alivepart
      ipoint(i)=count%ialive(i)
      igrid(i) = igridn(i,inest)
    enddo
    call sort2(alivepart,igrid,ipoint)

  ! Now visit all grid columns where particles are present
  ! by going through the sorted particles
!$OMP PARALLEL PRIVATE (igrold,kpart,ipart,igr,jy,ix,kz,lconv, &
!$OMP ktop,ztold,nage,ipconv,itage,lcalcflux,ithread)
    igrold = -1
#if (defined _OPENMP)
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts at 1
#else
    ithread = 1
#endif
!$OMP DO SCHEDULE(dynamic)
    do kpart=1,alivepart
      igr = igrid(kpart)
      if (igr .eq. -1) cycle
      ipart = ipoint(kpart)
      ! sumall = sumall + 1
      if (igr .ne. igrold) then
  ! we are in a new grid column
        jy = (igr-1)/nxn(inest)
        ix = igr - jy*nxn(inest) - 1

  ! Interpolate all meteorological data needed for the convection scheme
        psconv=(psn(ix,jy,1,mind1,inest)*dt2+ &
             psn(ix,jy,1,mind2,inest)*dt1)*dtt
        tt2conv=(tt2n(ix,jy,1,mind1,inest)*dt2+ &
             tt2n(ix,jy,1,mind2,inest)*dt1)*dtt
        td2conv=(td2n(ix,jy,1,mind1,inest)*dt2+ &
             td2n(ix,jy,1,mind2,inest)*dt1)*dtt
!!$        do kz=1,nconvlev+1    !old
        do kz=1,nuvz-1           !bugfix
          tconv(kz,ithread)=(tthn(ix,jy,kz+1,mind1,inest)*dt2+ &
               tthn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
          qconv(kz,ithread)=(qvhn(ix,jy,kz+1,mind1,inest)*dt2+ &
               qvhn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
        end do

  ! calculate translocation matrix
  !*******************************
        call calcmatrix(lconv,delt,cbasefluxn(ix,jy,inest,ithread),ithread)
        igrold = igr
        ktop = 0
      endif

  ! treat particle only if column has convection
      if (lconv .eqv. .true.) then
  ! assign new vertical position to particle
        ztold=real(part(ipart)%z)
        call redist(itime,ipart,ktop,ipconv,ithread)
  !      if (ipconv.le.0) sumconv = sumconv+1

  ! Calculate the gross fluxes across layer interfaces
  !***************************************************
        if (iflux.eq.1) then
          lcalcflux=.true.
          if (lagespectra.eq.1) then
            nage=0
            itage=abs(itime-part(ipart)%tstart)
            do inage=1,nageclass
              if ((itage.lt.lage(inage)).and.(part(ipart)%alive)) exit
              nage=inage
            end do
            if (nage.eq.nageclass) lcalcflux=.false.
            nage=nage+1
          else
            nage=1
          endif

          if (lcalcflux) &
            call calcfluxes(itime,nage,ipart,real(part(ipart)%xlon), &
               real(part(ipart)%ylat),ztold,ithread)
        endif
      endif !(lconv .eqv. .true.)

    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
  ! OpenMP Reduction for dynamically allocated arrays. This is done manually since this
  ! is not yet supported in most OpenMP versions
  !************************************************************************************
  if (iflux.eq.1) then
    do ithread=1,numthreads
      flux(:,:,:,:,:,:,:)=flux(:,:,:,:,:,:,:)+flux_omp(:,:,:,:,:,:,:,ithread)
    end do
  endif
  !--------------------------------------------------------------------------
  ! write(*,*)'############################################'
  ! write(*,*)'TIME=',&
  !    &  itime
  ! write(*,*)'fraction of particles under convection',&
  !    &  sumconv/(sumall+0.001)
  ! write(*,*)'total number of particles',&
  !    &  sumall
  ! write(*,*)'number of particles under convection',&
  !    &  sumconv
  ! write(*,*)'############################################'

  deallocate( igrid )
  deallocate( ipoint )
  deallocate( igridn )

  return
end subroutine convmix

subroutine calcmatrix(lconv,delt,cbmf,ithread)
  !                        o    i    o
  !*****************************************************************************
  !                                                                            *
  !  This subroutine calculates the matrix describing convective               *
  !  redistribution of mass in a grid column, using the subroutine             *
  !  convect43c.f provided by Kerry Emanuel.                                   *
  !                                                                            *
  !  Petra Seibert, Bernd C. Krueger, 2000-2001                                *
  !                                                                            *
  !*****************************************************************************
  ! Changes:                                                                   *
  !  changed by C. Forster, November 2003 - February 2004                      *
  !  array fmassfrac(nconvlevmax,nconvlevmax) represents                       *
  !  the convective redistribution matrix for the particles                    *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Merged calcmatrix and calcmatrix_gfs into one routine using if-then  *
  !       for meteo-type dependent code                                        *
  !*****************************************************************************
  !                                                                            *
  !  lconv        indicates whether there is convection in this cell, or not   *
  !  delt         time step for convection [s]                                 *
  !  cbmf         cloud base mass flux                                         *
  !  metdata_format format of metdata (ecmwf/gfs)                              *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use class_gribfile_mod
  use qvsat_mod

  implicit none

  integer,intent(in) :: ithread !OMP thread starting at 1
  real :: rlevmass,summe

  integer :: iflag, k, kk

  !1-d variables for convection
  !variables for redistribution matrix
  real :: cbmfold, precip, qprime
  real :: tprime, wd
  real :: delt,cbmf
  logical :: lconv

  lconv = .false.


  ! calculate pressure at eta levels for use in convect
  ! and assign temp & spec. hum. to 1D workspace
  ! -------------------------------------------------------

  ! pconv(1) is the pressure at the first level above ground
  ! phconv(k) is the pressure between levels k-1 and k
  ! dpr(k) is the pressure difference "around" tconv(k)
  ! phconv(kmax) must also be defined 1/2 level above pconv(kmax)
  ! Therefore, we define k = kuvz-1 and let kuvz start from 2
  ! top layer cannot be used for convection because p at top of this layer is
  ! not given


  phconv(1,ithread) = psconv
  ! Emanuel subroutine needs pressure in hPa, therefore convert all pressures
  ! do kuvz = 2,nuvz
  !   k = kuvz-1
  !   if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
  !     pconv(k) = (akz(kuvz) + bkz(kuvz)*psconv)
  !     phconv(kuvz) = (akm(kuvz) + bkm(kuvz)*psconv)
  !   else
  !     phconv(kuvz) =  0.5*(pconv(kuvz)+pconv(k))
  !   endif
  !   dpr(k) = phconv(k) - phconv(kuvz)
  !   qsconv(k) = f_qvsat( pconv(k), tconv(k) )

  ! initialize mass fractions
    ! do kk=1,nconvlev
    !   fmassfrac(k,kk)=0.
    ! end do
  ! end do
  ! LB 04.05.2021, replace above with array operations
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    pconv(1:nuvz-1,ithread) = (akz(2:nuvz) + bkz(2:nuvz)*psconv)
    phconv(2:nuvz,ithread) = (akm(2:nuvz) + bkm(2:nuvz)*psconv)
  else
    phconv(2:nuvz,ithread) =  0.5*(pconv(2:nuvz,ithread)+pconv(1:nuvz-1,ithread))
  endif
  dpr(1:nuvz-1,ithread) = phconv(1:nuvz-1,ithread) - phconv(2:nuvz,ithread)
  do k = 1,nuvz-1
    qsconv(k,ithread) = f_qvsat( pconv(k,ithread), tconv(k,ithread) )
  end do
  fmassfrac(1:nuvz-1,1:nconvlev,ithread)=0.
  ! LB end

  !note that Emanuel says it is important
  !a. to set this =0. every grid point
  !b. to keep this value in the calling programme in the iteration

  ! CALL CONVECTION
  !******************

  cbmfold = cbmf
  ! Convert pressures to hPa, as required by Emanuel scheme
  !********************************************************
  !!$    do k=1,nconvlev     !old
  ! do k=1,nconvlev+1      !bugfix
  !   pconv_hpa(k)=pconv(k)/100.
  !   phconv_hpa(k)=phconv(k)/100.
  ! end do
  ! phconv_hpa(nconvlev+1)=phconv(nconvlev+1)/100.
  ! LB 04.05.2021, replace above with array operations
  pconv_hpa(1:nconvlev+1,ithread)=pconv(1:nconvlev+1,ithread)/100.
  phconv_hpa(1:nconvlev+1,ithread)=phconv(1:nconvlev+1,ithread)/100.
  ! LB end
    
  call convect(nconvlevmax, nconvlev, delt, iflag, &
       precip, wd, tprime, qprime, cbmf, ithread)

  ! do not update fmassfrac and cloudbase massflux
  ! if no convection takes place or
  ! if a CFL criterion is violated in convect43c.f
  if (iflag .ne. 1 .and. iflag .ne. 4) then
   cbmf=cbmfold
   return
  endif

  ! do not update fmassfrac and cloudbase massflux
  ! if the old and the new cloud base mass
  ! fluxes are zero
  if (cbmf.le.0..and.cbmfold.le.0.) then
   cbmf=cbmfold
   return
  endif

  ! Update fmassfrac
  ! account for mass displaced from level k to level k

  lconv = .true.
  do k=1,nconvtop
    rlevmass = dpr(k,ithread)/ga
    summe = 0.
    do kk=1,nconvtop
      fmassfrac(k,kk,ithread) = delt*fmass(k,kk,ithread)
      summe = summe + fmassfrac(k,kk,ithread)
    end do
    fmassfrac(k,k,ithread)=fmassfrac(k,k,ithread) + rlevmass - summe
  end do
  ! LB 04.05.2021, replace above with array operations (not the problem)
  ! fmassfrac(1:nconvtop,1:nconvtop) = delt*fmass(1:nconvtop,1:nconvtop)
  ! do k=1, nconvtop
  !     fmassfrac(k, k) = fmassfrac(k, k) + dpr(k)/ga - sum(fmassfrac(k, 1:nconvtop))
  ! end do
  ! LB end
end subroutine calcmatrix

subroutine redist(itime,ipart,ktop,ipconv,ithread)

  !**************************************************************************
  ! Do the redistribution of particles due to convection
  ! This subroutine is called for each particle which is assigned
  ! a new vertical position randomly, based on the convective redistribution
  ! matrix
  !**************************************************************************

  ! Petra Seibert, Feb 2001, Apr 2001, May 2001, Jan 2002, Nov 2002 and
  ! Andreas Frank, Nov 2002

  ! Caroline Forster:  November 2004 - February 2005

  use par_mod
  use com_mod
  use random_mod
  use omp_lib
  use interpol_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use particle_mod
  use qvsat_mod

  implicit none

  integer,intent(in) :: ithread ! OMP thread starting at 1
  real,parameter :: const=r_air/ga
  integer :: ipart, ktop,ipconv,itime
  integer :: k, kz, levnew, levold

  real :: totlevmass, wsubpart
  real :: temp_levold,temp_levold1
  real :: sub_levold,sub_levold1
  real :: rn, dlevfrac
  real :: ztold,ffraction, dlogp
  real :: dz, dz1, dz2
#ifndef ETA
  real :: tv, tv1, tv2, tvold, pold, pint
#endif

  ! ipart   ... number of particle to be treated

  ipconv=1

  ! !  determine vertical grid position of particle in the eta system
  ! !****************************************************************
#ifdef ETA
  ztold = real(part(abs(ipart))%zeta)
  ! find old particle grid position
  levold = nconvtop
  do kz = 2, nconvtop
    if (wheight(kz) .le. ztold ) then
      levold = kz-1
      exit
    endif
  end do
#else

  ! determine height of the eta half-levels (uvzlev)
  ! do that only once for each grid column
  ! i.e. when ktop.eq.1
  !**************************************************************

  if (ktop .le. 1) then

    tvold=tt2conv*(1.+0.378*ew(td2conv,psconv)/psconv)
    pold=psconv
    uvzlev(1,ithread)=0.

    pint = phconv(2,ithread)
  !  determine next virtual temperatures
    tv1 = tconv(1,ithread)*(1.+0.608*qconv(1,ithread))
    tv2 = tconv(2,ithread)*(1.+0.608*qconv(2,ithread))
  !  interpolate virtual temperature to half-level
    tv = tv1 + (tv2-tv1)*(pconv(1,ithread)-phconv(2,ithread))/ &
      (pconv(1,ithread)-pconv(2,ithread))
    tv = tv1 + (tv2-tv1)*(pconv(1,ithread)-phconv(2,ithread))/ &
      (pconv(1,ithread)-pconv(2,ithread))
    if (abs(tv-tvold).gt.0.2) then
      uvzlev(2,ithread) = uvzlev(1,ithread) + &
           const*log(pold/pint)* &
           (tv-tvold)/log(tv/tvold)
    else
      uvzlev(2,ithread) = uvzlev(1,ithread)+ &
           const*log(pold/pint)*tv
    endif
    tvold=tv
    tv1=tv2
    pold=pint

  ! integrate profile (calculation of height agl of eta layers) as required
    do kz = 3, nconvtop+1
  !    note that variables defined in calcmatrix.f (pconv,tconv,qconv)
  !    start at the first real ECMWF model level whereas kz and
  !    thus uvzlev(kz) starts at the surface. uvzlev is defined at the
  !    half-levels (between the tconv, qconv etc. values !)
  !    Thus, uvzlev(kz) is the lower boundary of the tconv(kz) cell.
      pint = phconv(kz,ithread)
  !    determine next virtual temperatures
      tv2 = tconv(kz,ithread)*(1.+0.608*qconv(kz,ithread))
  !    interpolate virtual temperature to half-level
      tv = tv1 + (tv2-tv1)*(pconv(kz-1,ithread)-phconv(kz,ithread))/ &
           (pconv(kz-1,ithread)-pconv(kz,ithread))
      tv = tv1 + (tv2-tv1)*(pconv(kz-1,ithread)-phconv(kz,ithread))/ &
           (pconv(kz-1,ithread)-pconv(kz,ithread))
      if (abs(tv-tvold).gt.0.2) then
        uvzlev(kz,ithread) = uvzlev(kz-1,ithread) + &
             const*log(pold/pint)* &
             (tv-tvold)/log(tv/tvold)
      else
        uvzlev(kz,ithread) = uvzlev(kz-1,ithread)+ &
             const*log(pold/pint)*tv
      endif
      tvold=tv
      tv1=tv2
      pold=pint


    end do

    ktop = 2

  endif
  
  ztold = real(part(abs(ipart))%z)
  ! find old particle grid position
  levold = nconvtop
  do kz = 2, nconvtop
    if (uvzlev(kz,ithread) .ge. ztold ) then
      levold = kz-1
      exit
    endif
  end do
#endif

  ! If the particle is above the potentially convective domain, it will be skipped
  if (levold.ne.nconvtop) then

    ! now redistribute particles
    !****************************

    !  Choose a random number and find corresponding level of destination
    !  Random numbers to be evenly distributed in [0,1]

    rn = ran3(iseed2(ithread-1),ithread-1)

    ! initialize levnew

    levnew = levold

    ffraction = 0.
    totlevmass=dpr(levold,ithread)/ga
    loop1: do k = 1,nconvtop
    ! for backward runs use the transposed matrix
     if (ldirect.eq.1) then
       ffraction=ffraction+fmassfrac(levold,k,ithread) / totlevmass
     else
       ffraction=ffraction+fmassfrac(k,levold,ithread) / totlevmass
     endif
     if (rn.le.ffraction) then
       levnew=k
    ! avoid division by zero or a too small number
    ! if division by zero or a too small number happens the
    ! particle is assigned to the center of the grid cell
       if (ffraction.gt.1.e-20) then
        if (ldirect.eq.1) then
          dlevfrac = (ffraction-rn) / fmassfrac(levold,k,ithread) * totlevmass
        else
          dlevfrac = (ffraction-rn) / fmassfrac(k,levold,ithread) * totlevmass
        endif
       else
         dlevfrac = 0.5
       endif
       exit loop1
     endif
    end do loop1

    ! now assign new position to particle
#ifdef ETA
    if ((levnew.le.nconvtop).and.(levnew.ne.levold)) then
      dlogp = (1.-dlevfrac) * (wheight(levnew+1)-wheight(levnew))
      call set_zeta(ipart,wheight(levnew)+dlogp)
      if (part(abs(ipart))%zeta.ge.1.) call set_zeta(ipart,1.-(part(abs(ipart))%zeta-1.))
      if (part(abs(ipart))%zeta.eq.1.) call update_zeta(ipart,-1.e-4)
      if (ipconv.gt.0) ipconv=-1
    endif
#else
    if ((levnew.le.nconvtop).and.(levnew.ne.levold)) then
      dlogp = (1.-dlevfrac)* (log(phconv(levnew+1,ithread))-log(phconv(levnew,ithread)))
      pint = log(phconv(levnew,ithread))+dlogp
      dz1 = pint - log(phconv(levnew,ithread))
      dz2 = log(phconv(levnew+1,ithread)) - pint
      dz = dz1 + dz2
      call set_z(ipart,(uvzlev(levnew,ithread)*dz2+uvzlev(levnew+1,ithread)*dz1)/dz)
      if (part(abs(ipart))%z.lt.0.) call set_z(ipart,-1.*part(abs(ipart))%z)
      if (ipconv.gt.0) ipconv=-1
    endif
#endif
    ! displace particle according to compensating subsidence
    ! this is done to those particles, that were not redistributed
    ! by the matrix
    !**************************************************************

    if ((levnew.le.nconvtop).and.(levnew.eq.levold)) then

      ! determine compensating vertical velocity at the levels
      ! above and below the particel position
      ! increase compensating subsidence by the fraction that
      ! is displaced by convection to this level

      if (levold.gt.1) then
        temp_levold = tconv(levold-1,ithread) + &
            (tconv(levold,ithread)-tconv(levold-1,ithread)) &
            *(pconv(levold-1,ithread)-phconv(levold,ithread))/ &
            (pconv(levold-1,ithread)-pconv(levold,ithread))
        ! LB: the units seem to not add up correctly, but adding lsynctime gives incorrect mixing
        ! in the lowest km and too many right above the ground
        ! sub_levold = sub(levold,ithread)/(1.-ga*sub(levold,ithread)*lsynctime/dpr(levold,ithread))
        sub_levold = sub(levold,ithread)/(1.-ga*sub(levold,ithread)/dpr(levold,ithread))
        wsub(levold,ithread)=-1.*sub_levold*r_air*temp_levold/(phconv(levold,ithread))
      else
        wsub(levold,ithread)=0.
      endif

      temp_levold1 = tconv(levold,ithread) + &
          (tconv(levold+1,ithread)-tconv(levold,ithread)) &
          *(pconv(levold,ithread)-phconv(levold+1,ithread))/ &
          (pconv(levold,ithread)-pconv(levold+1,ithread))
      !sub_levold1 = sub(levold+1,ithread)/(1.-ga*sub(levold+1,ithread)*lsynctime/dpr(levold+1,ithread))
      sub_levold1 = sub(levold+1,ithread)/(1.-ga*sub(levold+1,ithread)/dpr(levold+1,ithread))
      wsub(levold+1,ithread)=-1.*sub_levold1*r_air*temp_levold1/ &
          (phconv(levold+1,ithread))

      ! interpolate wsub to the vertical particle position
#ifdef ETA
      ztold = real(part(abs(ipart))%zeta)
      dz1 = ztold - wheight(levold)
      dz2 = wheight(levold+1) - ztold
      dz = dz1 + dz2

      ! Convert z(eta) to z(m) in order to add subsidence
      call update_zeta_to_z(itime, ipart)
      ! call zeta_to_z(itime,part(abs(ipart))%xlon,part(abs(ipart))%ylat, &
      !   part(abs(ipart))%zeta,part(abs(ipart))%z)

      wsubpart = (dz2*wsub(levold,ithread)+dz1*wsub(levold+1,ithread))/dz

      call update_z(ipart,wsubpart*real(lsynctime))

      if (part(abs(ipart))%z.lt.0.) call set_z(ipart,-1.*part(abs(ipart))%z)

      ! Convert new z(m) back to z(eta)
      call update_z_to_zeta(itime, ipart)
          
#else
      ztold = real(part(abs(ipart))%z)
      dz1 = ztold - uvzlev(levold,ithread)
      dz2 = uvzlev(levold+1,ithread) - ztold
      dz = dz1 + dz2

      wsubpart = (dz2*wsub(levold,ithread)+dz1*wsub(levold+1,ithread))/dz

      call update_z(ipart,wsubpart*real(lsynctime))

      if (part(abs(ipart))%z.lt.0.) call set_z(ipart,-1.*part(abs(ipart))%z)

#endif
    endif      !(levnew.le.nconvtop.and.levnew.eq.levold)
  endif
  ! Maximum altitude .5 meter below uppermost model level
  !*******************************************************

#ifdef ETA
  if (part(abs(ipart))%zeta .lt. uvheight(nz)) call set_zeta(ipart,uvheight(nz)+1.e-4)
  if (part(abs(ipart))%zeta.ge.1.) call set_zeta(ipart,1.-(part(abs(ipart))%zeta-1.))
  if (part(abs(ipart))%zeta.eq.1.) call update_zeta(ipart,-1.e-4)
#else
  if (part(abs(ipart))%z .gt. height(nz)-0.5) call set_z(ipart,height(nz)-0.5)
#endif

end subroutine redist

!**************************************************************************
!****                       SUBROUTINE CONVECT                        *****
!****                          VERSION 4.3c                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************
!
  SUBROUTINE CONVECT &
         (ND,  NL,   DELT, IFLAG, &
         PRECIP, WD,   TPRIME, QPRIME, CBMF, ithread  )
  !
  !-cv *************************************************************************
  !-cv C. Forster, November 2003 - May 2004:
  !-cv
  !-cv The subroutine has been downloaded from Kerry Emanuel's homepage,
  !-cv where further infos on the convection scheme can be found
  !-cv http://www-paoc.mit.edu/~emanuel/home.html
  !-cv
  !-cv The following changes have been made to integrate this subroutine
  !-cv into FLEXPART
  !-cv
  !-cv Putting most of the variables in a new common block
  !-cv renaming eps to eps0 because there is some eps already in includepar
  !-cv
  !-cv removing the arrays U,V,TRA and related arrays
  !-cv
  !-cv renaming the original arrays T,Q,QS,P,PH to
  !-cv TCONV,QCONV,QSCONV,PCONV_HPA,PHCONV_HPA
  !-cv
  !-cv Initialization of variables has been put into parameter statements
  !-cv instead of assignment of values at each call, in order to save 
  !-cv computation time.
  !***************************************************************************
  !
  !-----------------------------------------------------------------------------
  !    *** On input:      ***
  !
  !T:   Array of absolute temperature (K) of dimension ND, with first
  !      index corresponding to lowest model level. Note that this array
  !      will be altered by the subroutine if dry convective adjustment
  !      occurs and if IPBL is not equal to 0.
  !
  !Q:   Array of specific humidity (gm/gm) of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !QS:  Array of saturation specific humidity of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
  !       index corresponding with the lowest model level. Defined at
  !       same levels as T. Note that this array will be altered if
  !       dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !V:   Same as U but for meridional velocity.
  !
  !TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
  !       where NTRA is the number of different tracers. If no
  !       convective tracer transport is needed, define a dummy
  !       input array of dimension (ND,1). Tracers are defined at
  !       same vertical levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !P:   Array of pressure (mb) of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T.
  !
  !PH:  Array of pressure (mb) of dimension ND+1, with first index
  !       corresponding to lowest level. These pressures are defined at
  !       levels intermediate between those of P, T, Q and QS. The first
  !       value of PH should be greater than (i.e. at a lower level than)
  !       the first value of the array P.
  !
  !ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
  !
  !NL:  The maximum number of levels to which convection can
  !       penetrate, plus 1.
  !       NL MUST be less than or equal to ND-1.
  !
  !NTRA:The number of different tracers. If no tracer transport
  !       is needed, set this equal to 1. (On most compilers, setting
  !       NTRA to 0 will bypass tracer calculation, saving some CPU.)
  !
  !DELT: The model time step (sec) between calls to CONVECT
  !
  !----------------------------------------------------------------------------
  !    ***   On Output:         ***
  !
  !IFLAG: An output integer whose value denotes the following:
  !
  !           VALUE                        INTERPRETATION
  !           -----                        --------------
  !             0               No moist convection; atmosphere is not
  !                             unstable, or surface temperature is less
  !                             than 250 K or surface specific humidity
  !                             is non-positive.
  !
  !             1               Moist convection occurs.
  !
  !             2               No moist convection: lifted condensation
  !                             level is above the 200 mb level.
  !
  !             3               No moist convection: cloud base is higher
  !                             then the level NL-1.
  !
  !             4               Moist convection occurs, but a CFL condition
  !                             on the subsidence warming is violated. This
  !                             does not cause the scheme to terminate.
  !
  !FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
  !        grid levels as T, Q, QS and P.
  !
  !FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
  !        defined at same grid levels as T, Q, QS and P.
  !
  !FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
  !        defined at same grid levels as T.
  !
  !FV:   Same as FU, but for forcing of meridional velocity.
  !
  !FTRA: Array of forcing of tracer content, in tracer mixing ratio per
  !        second, defined at same levels as T. Dimensioned (ND,NTRA).
  !
  !PRECIP: Scalar convective precipitation rate (mm/day).
  !
  !WD:    A convective downdraft velocity scale. For use in surface
  !        flux parameterizations. See convect.ps file for details.
  !
  !TPRIME: A convective downdraft temperature perturbation scale (K).
  !         For use in surface flux parameterizations. See convect.ps
  !         file for details.
  !
  !QPRIME: A convective downdraft specific humidity
  !         perturbation scale (gm/gm).
  !         For use in surface flux parameterizations. See convect.ps
  !         file for details.
  !
  !CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
  !         BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
  !         ITS NEXT CALL. That is, the value of CBMF must be "remembered"
  !         by the calling program between calls to CONVECT.
  !
  !-----------------------------------------------------------------------------
  !
  !    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
  !    ***                OR EQUAL TO  ND + 1                    ***
  !
  !
  use par_mod

  implicit none
  !
  !-cv====>Begin Module CONVECT    File convect.f      Undeclared variables
  !
  !Argument variables
  !
  integer, intent(in) :: ithread ! OMP thread starting at 1
  integer :: iflag, nd, nl
  !
  real :: cbmf, delt, precip, qprime, tprime, wd
  !
  !Local variables
  !
  integer :: i, icb, ihmin, inb, inb1, j, jtt, k
  integer :: nk
  !
  real :: ad, afac, ahmax, ahmin, alt, altem
  real :: am, amp1, anum, asij, awat, b6, bf2, bsum, by
  real :: byp, c6, cape, capem, cbmfold, chi, coeff
  real :: cpinv, cwat, damps
  real :: defrac, dei, delm, delp, delt0, delti, denom, dhdp
  real :: dpinv, dtma, dtmin, dtpbl, elacrit, ents
  real :: epmax, fac, fqold, frac, ftold
  real :: plcl, qp1, qsm, qstm, qti, rat
  real :: revap, rh, scrit, sigt, sjmax
  real :: sjmin, smid, smin, stemp, tca
  real :: tvaplcl, tvpplcl, tvx, tvy, wdtrain

  !integer jc,jn
  !real alvnew,a2,ahm,alv,rm,sum,qnew,dphinv,tc,thbar,tnew,x
  !REAL TOLD(NA)

  real :: FUP(NA),FDOWN(NA),FT(NA),FQ(NA)
  !
  !-cv====>End Module   CONVECT    File convect.f

  INTEGER :: NENT(NA)
  REAL :: M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
  REAL :: SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
  REAL :: QP(NA),EP(NA),WT(NA),EVAP(NA),CLW(NA)
  REAL :: SIGP(NA),TP(NA),CPN(NA)
  REAL :: LV(NA),LVCP(NA),H(NA),HP(NA),GZ(NA),HM(NA)
  !
  ! -----------------------------------------------------------------------
  !
  !   ***                     Specify Switches                         ***
  !
  !   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
  !   ***    Any other value results in dry adiabatic adjustment       ***
  !   ***     (Zero value recommended for use in models with           ***
  !   ***                   boundary layer schemes)                    ***
  !
  !   ***   MINORIG: Lowest level from which convection may originate  ***
  !   ***     (Should be first model level at which T is defined       ***
  !   ***      for models using bulk PBL schemes; otherwise, it should ***
  !   ***      be the first model level at which T is defined above    ***
  !   ***                      the surface layer)                      ***
  !
  INTEGER,PARAMETER :: IPBL=0
  INTEGER,PARAMETER :: MINORIG=1
  !
  !------------------------------------------------------------------------------
  !
  !   ***                    SPECIFY PARAMETERS                        ***
  !
  !   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
  !   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
  !   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
  !   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
  !   ***               BETWEEN 0 C AND TLCRIT)                        ***
  !   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
  !   ***                       FORMULATION                            ***
  !   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
  !   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
  !   ***                        OF CLOUD                              ***
  !   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
  !   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
  !   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  !   ***                          OF RAIN                             ***
  !   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  !   ***                          OF SNOW                             ***
  !   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
  !   ***                         TRANSPORT                            ***
  !   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
  !   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
  !   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
  !   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
  !   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
  !   ***                   (DAMP MUST BE LESS THAN 1)                 ***
  !
  REAL,PARAMETER :: ELCRIT=.0011
  REAL,PARAMETER :: TLCRIT=-55.0
  REAL,PARAMETER :: ENTP=1.5
  REAL,PARAMETER :: SIGD=0.05
  REAL,PARAMETER :: SIGS=0.12
  REAL,PARAMETER :: OMTRAIN=50.0
  REAL,PARAMETER :: OMTSNOW=5.5
  REAL,PARAMETER :: COEFFR=1.0
  REAL,PARAMETER :: COEFFS=0.8
  REAL,PARAMETER :: CU=0.7
  REAL,PARAMETER :: BETA=10.0
  REAL,PARAMETER :: DTMAX=0.9
  REAL,PARAMETER :: ALPHA=0.025  !original 0.2
  REAL,PARAMETER :: DAMP=0.1
  !
  !   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
  !   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
  !   ***             THESE SHOULD BE CONSISTENT WITH             ***
  !   ***              THOSE USED IN CALLING PROGRAM              ***
  !   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
  !
  REAL,PARAMETER :: CPD=1005.7
  REAL,PARAMETER :: CPV=1870.0
  REAL,PARAMETER :: CL=2500.0
  REAL,PARAMETER :: RV=461.5
  REAL,PARAMETER :: RD=287.04
  REAL,PARAMETER :: LV0=2.501E6
  REAL,PARAMETER :: G=9.81
  REAL,PARAMETER :: ROWL=1000.0
  !
  REAL,PARAMETER :: CPVMCL=CL-CPV
  REAL,PARAMETER :: EPS0=RD/RV
  REAL,PARAMETER :: EPSI=1./EPS0
  REAL,PARAMETER :: GINV=1.0/G
  REAL,PARAMETER :: EPSILON=1.e-20

  ! EPSILON IS A SMALL NUMBER USED TO EXCLUDE MASS FLUXES OF ZERO
  !
  DELTI=1.0/DELT
  !
  !      ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
  !

  FT(:NL+1)=0.0
  FQ(:NL+1)=0.0
  FDOWN(:NL+1)=0.0
  sub(:NL+1,ithread)=0.0
  FUP(:NL+1)=0.0
  M(:NL+1)=0.0
  MP(:NL+1)=0.0
  fmass(:NL+1,:NL+1,ithread)=0.0
  MENT(:NL+1,:NL+1)=0.0

  PRECIP=0.0
  WD=0.0
  TPRIME=0.0
  QPRIME=0.0
  IFLAG=0
  !
  !  IF(IPBL.NE.0)THEN
  !
  !***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
  !
  !  JC=0
  !  DO 30 I=NL-1,1,-1
  !   JN=0
  !    SUM=TH(I)*(1.+qconv(I)*EPSI-qconv(I))
  !   DO 10 J=I+1,NL
  !    SUM=SUM+TH(J)*(1.+qconv(J)*EPSI-qconv(J))
  !    THBAR=SUM/REAL(J+1-I)
  !    IF((TH(J)*(1.+qconv(J)*EPSI-qconv(J))).LT.THBAR)JN=J
  !  10    CONTINUE
  !   IF(I.EQ.1)JN=MAX(JN,2)
  !   IF(JN.EQ.0)GOTO 30
  !  12    CONTINUE
  !   AHM=0.0
  !   RM=0.0
  !   DO 15 J=I,JN
  !    AHM=AHM+(CPD*(1.-qconv(J))+qconv(J)*CPV)*tconv(J)*
  !    +   (phconv_hpa(J)-phconv_hpa(J+1))
  !    RM=RM+qconv(J)*(phconv_hpa(J)-phconv_hpa(J+1))
  !  15    CONTINUE
  !   DPHINV=1./(phconv_hpa(I)-phconv_hpa(JN+1))
  !   RM=RM*DPHINV
  !   A2=0.0
  !   DO 20 J=I,JN
  !    qconv(J)=RM
  !    RDCP=(RD*(1.-qconv(J))+qconv(J)*RV)/
  !    1     (CPD*(1.-qconv(J))+qconv(J)*CPV)
  !    X=(0.001*pconv_hpa(J))**RDCP
  !    TOLD(J)=tconv(J)
  !    tconv(J)=X
  !    A2=A2+(CPD*(1.-qconv(J))+qconv(J)*CPV)*X*
  !    1    (phconv_hpa(J)-phconv_hpa(J+1))
  !  20    CONTINUE
  !   DO 25 J=I,JN
  !    TH(J)=AHM/A2
  !    tconv(J)=tconv(J)*TH(J)
  !    TC=TOLD(J)-273.15
  !    ALV=LV0-CPVMCL*TC
  !    qsconv(J)=qsconv(J)+qsconv(J)*(1.+qsconv(J)*(EPSI-1.))*ALV*
  !    1    (tconv(J)- TOLD(J))/(RV*TOLD(J)*TOLD(J))
  ! if (qslev(j) .lt. 0.) then
  !   write(*,*) 'qslev.lt.0 ',j,qslev
  ! endif
  !  25    CONTINUE
  !   IF((TH(JN+1)*(1.+qconv(JN+1)*EPSI-qconv(JN+1))).LT.
  !    1    (TH(JN)*(1.+qconv(JN)*EPSI-qconv(JN))))THEN
  !    JN=JN+1
  !    GOTO 12
  !   END IF
  !   IF(I.EQ.1)JC=JN
  !  30   CONTINUE
  !
  !   ***   Remove any supersaturation that results from adjustment ***
  !
  !IF(JC.GT.1)THEN
  ! DO 38 J=1,JC
  !    IF(qsconv(J).LT.qconv(J))THEN
  !     ALV=LV0-CPVMCL*(tconv(J)-273.15)
  !     TNEW=tconv(J)+ALV*(qconv(J)-qsconv(J))/(CPD*(1.-qconv(J))+
  !    1      CL*qconv(J)+qsconv(J)*(CPV-CL+ALV*ALV/(RV*tconv(J)*tconv(J))))
  !     ALVNEW=LV0-CPVMCL*(TNEW-273.15)
  !     QNEW=(ALV*qconv(J)-(TNEW-tconv(J))*(CPD*(1.-qconv(J))
  !    1     +CL*qconv(J)))/ALVNEW
  !     PRECIP=PRECIP+24.*3600.*1.0E5*(phconv_hpa(J)-phconv_hpa(J+1))*
  !    1      (qconv(J)-QNEW)/(G*DELT*ROWL)
  !     tconv(J)=TNEW
  !     qconv(J)=QNEW
  !     qsconv(J)=QNEW
  !    END IF
  !  38  CONTINUE
  !END IF
  !
  !END IF
  !
  !  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
  !
  GZ(1)=0.0
  CPN(1)=CPD*(1.-qconv(1,ithread))+qconv(1,ithread)*CPV
  H(1)=tconv(1,ithread)*CPN(1)
  LV(1)=LV0-CPVMCL*(tconv(1,ithread)-273.15)
  HM(1)=LV(1)*qconv(1,ithread)
  TV(1)=tconv(1,ithread)*(1.+qconv(1,ithread)*EPSI-qconv(1,ithread))
  AHMIN=1.0E12
  IHMIN=NL

  DO I=2,NL+1
    TVX=tconv(I,ithread)*(1.+qconv(I,ithread)*EPSI-qconv(I,ithread))
    TVY=tconv(I-1,ithread)*(1.+qconv(I-1,ithread)*EPSI-qconv(I-1,ithread))
    GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(pconv_hpa(I-1,ithread)-pconv_hpa(I,ithread))/ &
         phconv_hpa(I,ithread)
    CPN(I)=CPD*(1.-qconv(I,ithread))+CPV*qconv(I,ithread)
    H(I)=tconv(I,ithread)*CPN(I)+GZ(I)
    LV(I)=LV0-CPVMCL*(tconv(I,ithread)-273.15)
    HM(I)=(CPD*(1.-qconv(I,ithread))+CL*qconv(I,ithread))*(tconv(I,ithread)- &
      tconv(1,ithread)) + LV(I)*qconv(I,ithread)+GZ(I)
    TV(I)=tconv(I,ithread)*(1.+qconv(I,ithread)*EPSI-qconv(I,ithread))
!
!  ***  Find level of minimum moist static energy    ***
!
    IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
      AHMIN=HM(I)
      IHMIN=I
    END IF
  END DO
  IHMIN=MIN(IHMIN, NL-1)
  !
  !  ***     Find that model level below the level of minimum moist       ***
  !  ***  static energy that has the maximum value of moist static energy ***
  !
  AHMAX=0.0
  !  ***  bug fixed: need to assign an initial value to NK
  !  HSO, 05.08.2009
  NK=MINORIG
  DO I=MINORIG,IHMIN
    IF(HM(I).GT.AHMAX)THEN
      NK=I
      AHMAX=HM(I)
    END IF
  END DO
  ! LB 04.05.2021, replace above with array operations (maxloc not working)
  ! NK=MINORIG+maxloc(HM(MINORIG:IHMIN))-1

  !
  !  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
  !  ***                          ARE REASONABLE                         ***
  !  ***      Skip convection if HM increases monotonically upward       ***
  !
  IF(tconv(NK,ithread).LT.250.0.OR.qconv(NK,ithread).LE.0.0.OR.IHMIN.EQ.(NL-1)) THEN
    IFLAG=0
    CBMF=0.0
    RETURN
  END IF
  !
  !   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
  !   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
  !
  RH=qconv(NK,ithread)/qsconv(NK,ithread)
  CHI=tconv(NK,ithread)/(1669.0-122.0*RH-tconv(NK,ithread))
  PLCL=pconv_hpa(NK,ithread)*(RH**CHI)
  IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
    IFLAG=2
    CBMF=0.0
    RETURN
  END IF
  !
  !   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
  !
  ICB=NL-1
  DO I=NK+1,NL
    IF(pconv_hpa(I,ithread).LT.PLCL)THEN
      ICB=MIN(ICB,I)
    END IF
  END DO
  IF(ICB.GE.(NL-1))THEN
    IFLAG=3
    CBMF=0.0
    RETURN
  END IF
  !
  !   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
  !
  !   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
  !   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
  !   ***                   LIQUID WATER CONTENT                             ***
  !
  CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,1,ithread)
  TVP(NK:ICB)=TVP(NK:ICB)-TP(NK:ICB)*qconv(NK,ithread)
  !
  !   ***  If there was no convection at last time step and parcel    ***
  !   ***       is stable at ICB then skip rest of calculation        ***
  !
  IF(CBMF.EQ.0.0.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
    IFLAG=0
    RETURN
  END IF
  !
  !   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
  !
  IF(IFLAG.NE.4)IFLAG=1
  !
  !   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
  !
  CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,2,ithread)
  !
  !   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
  !   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
  !   ***      THESE MAY BE FUNCTIONS OF TP(I), pconv_hpa(I) AND CLW(I)     ***
  !
  EP(1:NK)=0.0
  SIGP(1:NL)=SIGS

  DO I=NK+1,NL
    TCA=TP(I)-273.15
    IF(TCA.GE.0.0)THEN
      ELACRIT=ELCRIT
    ELSE
      ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
    END IF
    ELACRIT=MAX(ELACRIT,0.0)
    EPMAX=0.999
    EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
    EP(I)=MAX(EP(I),0.0)
    EP(I)=MIN(EP(I),EPMAX)
    SIGP(I)=SIGS
  END DO
  ! LB 04.05.2021, replace above with array operations 
  ! (this makes it less readable, and not any faster)
  ! PROBLEM 1 is within the statement below
  ! EPMAX=0.999
  ! where ((TP(NK+1:NL)-273.15).ge.0.0)
  !   EP(NK+1:NL)=EPMAX*(1.0-max(ELCRIT, 0.0)/MAX(CLW(NK+1:NL),1.0E-8))
  ! elsewhere
  !   EP(NK+1:NL)=EPMAX*(1.0-max(ELCRIT*(1.0-TCA/TLCRIT), 0.0)/MAX(CLW(NK+1:NL),1.0E-8))
  ! end where
  ! where (EP(NK+1:NL).lt.0.0)
  !   EP(NK+1:NL)=0.0
  ! elsewhere (EP(NK+1:NL).gt.EPMAX)
  !   EP(NK+1:NL)=EPMAX
  ! end where

  !
  !   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
  !   ***                    VIRTUAL TEMPERATURE                    ***
  ! !
  TVP(ICB+1:NL)=TVP(ICB+1:NL)-TP(ICB+1:NL)*qconv(NK,ithread)
  TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
  !
  !   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***

  HP(:NL+1)=H(:NL+1)
  NENT(:NL+1)=0
  WATER(:NL+1)=0.0
  EVAP(:NL+1)=0.0
  WT(:NL+1)=OMTSNOW
  LVCP(:NL+1)=LV(:NL+1)/CPN(:NL+1)
  ELIJ(:NL+1,:NL+1)=0.0
  SIJ(:NL+1,:NL+1)=0.0
  DO I=1,NL+1
    QENT(I,:NL+1)=qconv(:NL+1,ithread)
  END DO  
  QP(1)=qconv(1,ithread)
  QP(2:NL+1)=qconv(:NL,ithread)

  !
  !  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
  !  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
  !  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
  !
  CAPE=0.0
  CAPEM=0.0
  INB=ICB+1
  INB1=INB
  BYP=0.0
  DO I=ICB+1,NL-1
    BY=(TVP(I)-TV(I))*(phconv_hpa(I,ithread)-phconv_hpa(I+1,ithread))/ &
      pconv_hpa(I,ithread)
    CAPE=CAPE+BY
    IF(BY.GE.0.0)INB1=I+1
    IF(CAPE.GT.0.0)THEN
      INB=I+1
      BYP=(TVP(I+1)-TV(I+1))*(phconv_hpa(I+1,ithread)-phconv_hpa(I+2,ithread))/ &
           pconv_hpa(I+1,ithread)
      CAPEM=CAPE
    END IF
  END DO
  INB=MAX(INB,INB1)
  CAPE=CAPEM+BYP
  DEFRAC=CAPEM-CAPE
  DEFRAC=MAX(DEFRAC,0.001)
  FRAC=-CAPE/DEFRAC
  FRAC=MIN(FRAC,1.0)
  FRAC=MAX(FRAC,0.0)
  !
  !   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
  !
  HP(ICB:INB)=H(NK)+(LV(ICB:INB)+(CPD-CPV)*tconv(ICB:INB,ithread))* &
    EP(ICB:INB)*CLW(ICB:INB)
  !
  !   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
  !   ***                   AT EACH MODEL LEVEL                       ***
  !
  
  !
  !   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
  !   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
  !
  TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(pconv_hpa(ICB-1,ithread)-PLCL)/ &
       (CPN(ICB-1)*pconv_hpa(ICB-1,ithread))
  TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-pconv_hpa(ICB,ithread))/ &
       (pconv_hpa(ICB,ithread)-pconv_hpa(ICB+1,ithread))
  DTPBL=0.0

  DTPBL=sum((TVP(NK:ICB-1)-TV(NK:ICB-1))*(phconv_hpa(NK:ICB-1,ithread)- &
    phconv_hpa(NK+1:ICB,ithread)))/ &
    (phconv_hpa(NK,ithread)-phconv_hpa(ICB,ithread))
  DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
  DTMA=DTMIN
  !
  !   ***  ADJUST CLOUD BASE MASS FLUX   ***
  !
  CBMFOLD=CBMF
  ! *** C. Forster: adjustment of CBMF is not allowed to depend on FLEXPART timestep
  DELT0=DELT/3.
  DAMPS=DAMP*DELT/DELT0
  CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA
  CBMF=MAX(CBMF,0.0)
  !
  !   *** If cloud base mass flux is zero, skip rest of calculation  ***
  !
  IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN
    RETURN
  END IF

  !
  !   ***   CALCULATE RATES OF MIXING,  M(I)   ***
  M(ICB)=0.0
  M(ICB+1:INB1)=ABS(TV(ICB+1:INB1)-TVP(ICB+1:INB1))+ &
        ENTP*0.02*(phconv_hpa(ICB+1:INB1,ithread)-phconv_hpa(ICB+2:INB1+1,ithread))
  M(INB1:INB)=ABS(TV(INB1)-TVP(INB1))+ &
        ENTP*0.02*(phconv_hpa(INB1,ithread)-phconv_hpa(INB1+1,ithread))
  M(ICB+1:INB)=CBMF*M(ICB+1:INB)/sum(M(ICB+1:INB))

  !
  !   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
  !   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
  !   ***                        FRACTION (SIJ)                          ***
  !
  DO I=ICB+1,INB
    QTI=qconv(NK,ithread)-EP(I)*CLW(I)
    DO J=ICB,INB
      BF2=1.+LV(J)*LV(J)*qsconv(J,ithread)/(RV*tconv(J,ithread)*tconv(J,ithread)*CPD)
      ANUM=H(J)-HP(I)+(CPV-CPD)*tconv(J,ithread)*(QTI-qconv(J,ithread))
      DENOM=H(I)-HP(I)+(CPD-CPV)*(qconv(I,ithread)-QTI)*tconv(J,ithread)
      DEI=DENOM
      IF(ABS(DEI).LT.0.01)DEI=0.01
      SIJ(I,J)=ANUM/DEI
      SIJ(I,I)=1.0
      ALTEM=SIJ(I,J)*qconv(I,ithread)+(1.-SIJ(I,J))*QTI-qsconv(J,ithread)
      ALTEM=ALTEM/BF2
      CWAT=CLW(J)*(1.-EP(J))
      STEMP=SIJ(I,J)
      IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR. &
            ALTEM.GT.CWAT).AND.J.GT.I)THEN
        ANUM=ANUM-LV(J)*(QTI-qsconv(J,ithread)-CWAT*BF2)
        DENOM=DENOM+LV(J)*(qconv(I,ithread)-QTI)
        IF(ABS(DENOM).LT.0.01)DENOM=0.01
        SIJ(I,J)=ANUM/DENOM
        ALTEM=SIJ(I,J)*qconv(I,ithread)+(1.-SIJ(I,J))*QTI-qsconv(J,ithread)
        ALTEM=ALTEM-(BF2-1.)*CWAT
      END IF
      IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
        QENT(I,J)=SIJ(I,J)*qconv(I,ithread)+(1.-SIJ(I,J))*QTI
        ELIJ(I,J)=ALTEM
        ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
        MENT(I,J)=M(I)/(1.-SIJ(I,J))
        NENT(I)=NENT(I)+1
      END IF
      SIJ(I,J)=MAX(0.0,SIJ(I,J))
      SIJ(I,J)=MIN(1.0,SIJ(I,J))
    END DO
  !
  !   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
  !   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
  !
    IF(NENT(I).EQ.0)THEN
      MENT(I,I)=M(I)
      QENT(I,I)=qconv(NK,ithread)-EP(I)*CLW(I)
      ELIJ(I,I)=CLW(I)
      SIJ(I,I)=1.0
    END IF
  END DO
  SIJ(INB,INB)=1.0
  ! LB 04.05.2021, Attempt to array the loop above: PROBLEM 2 is here
  ! DO J=ICB,INB
  !   BF2=1.+LV(J)*LV(J)*qsconv(J)/(RV*tconv(J)*tconv(J)*CPD)
  !   CWAT=CLW(J)*(1.-EP(J))
  !   DO I=ICB+1,INB
  !     QTI=qconv(NK)-EP(I)*CLW(I)
  !     ANUM=H(J)-HP(I)+(CPV-CPD)*tconv(J)*(QTI-qconv(J))
  !     DENOM=H(I)-HP(I)+(CPD-CPV)*(qconv(I)-QTI)*tconv(J)
  !     DEI=DENOM
  !     IF(I.EQ.J)THEN
  !       SIJ(I,I)=1.0
  !     ELSE IF(ABS(DENOM).LT.0.01)THEN
  !       SIJ(I,J)=ANUM/0.01
  !     ELSE
  !       SIJ(I,J)=ANUM/DENOM
  !     END IF
  !     ALTEM=(SIJ(I,J)*qconv(I)+(1.-SIJ(I,J))*QTI-qsconv(J))/BF2
  !     IF((SIJ(I,J).LT.0.0.OR.SIJ(I,J).GT.1.0.OR. &
  !           ALTEM.GT.CWAT).AND.J.GT.I)THEN
  !       ANUM=ANUM-LV(J)*(QTI-qsconv(J)-CWAT*BF2)
  !       DENOM=DENOM+LV(J)*(qconv(I)-QTI)
  !       IF(ABS(DENOM).LT.0.01)DENOM=0.01
  !       SIJ(I,J)=ANUM/DENOM
  !       ALTEM=SIJ(I,J)*qconv(I)+(1.-SIJ(I,J))*QTI-qsconv(J)
  !       ALTEM=ALTEM-(BF2-1.)*CWAT
  !     END IF
  !     IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
  !       QENT(I,J)=SIJ(I,J)*qconv(I)+(1.-SIJ(I,J))*QTI
  !       ELIJ(I,J)=ALTEM
  !       ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
  !       MENT(I,J)=M(I)/(1.-SIJ(I,J))
  !       NENT(I)=NENT(I)+1
  !     END IF
  !     SIJ(I,J)=MAX(0.0,SIJ(I,J))
  !     SIJ(I,J)=MIN(1.0,SIJ(I,J))
  !   END DO
  ! END DO
  ! !
  ! !   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
  ! !   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
  ! !
  ! do I=ICB+1,INB
  !   IF(NENT(I).EQ.0)THEN
  !     MENT(I,I)=M(I)
  !     QENT(I,I)=qconv(NK)-EP(I)*CLW(I)
  !     ELIJ(I,I)=CLW(I)
  !     SIJ(I,I)=1.0
  !   END IF
  ! END DO
  ! SIJ(INB,INB)=1.0


  !
  !   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
  !   ***              PROBABILITIES OF MIXING                     ***
  !
  ! LB 04.05.2021, depending on how often NENT.ne.0, reversing the loop could
  ! speed it up...
  DO I=ICB+1,INB
    IF(NENT(I).NE.0)THEN
      QP1=qconv(NK,ithread)-EP(I)*CLW(I)
      ANUM=H(I)-HP(I)-LV(I)*(QP1-qsconv(I,ithread))
      DENOM=H(I)-HP(I)+LV(I)*(qconv(I,ithread)-QP1)
      IF(ABS(DENOM).LT.0.01)DENOM=0.01
      SCRIT=ANUM/DENOM
      ALT=QP1-qsconv(I,ithread)+SCRIT*(qconv(I,ithread)-QP1)
      IF(ALT.LT.0.0)SCRIT=1.0
      SCRIT=MAX(SCRIT,0.0)
      ASIJ=0.0
      SMIN=1.0
      DO J=ICB,INB
        IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
          IF(J.GT.I)THEN
            SMID=MIN(SIJ(I,J),SCRIT)
            SJMAX=SMID
            SJMIN=SMID
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
              SMIN=SMID
              SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
              SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
              SJMIN=MIN(SJMIN,SCRIT)
            END IF
          ELSE
            SJMAX=MAX(SIJ(I,J+1),SCRIT)
            SMID=MAX(SIJ(I,J),SCRIT)
            SJMIN=0.0
            IF(J.GT.1)SJMIN=SIJ(I,J-1)
            SJMIN=MAX(SJMIN,SCRIT)
          END IF
          DELP=ABS(SJMAX-SMID)
          DELM=ABS(SJMIN-SMID)
          ASIJ=ASIJ+(DELP+DELM)*(phconv_hpa(J,ithread)-phconv_hpa(J+1,ithread))
          MENT(I,J)=MENT(I,J)*(DELP+DELM)* &
               (phconv_hpa(J,ithread)-phconv_hpa(J+1,ithread))
        END IF
      END DO
      ASIJ=MAX(1.0E-21,ASIJ)
      ASIJ=1.0/ASIJ
      DO J=ICB,INB
        MENT(I,J)=MENT(I,J)*ASIJ
      END DO
      BSUM=0.0
      DO J=ICB,INB
        BSUM=BSUM+MENT(I,J)
      END DO
      IF(BSUM.LT.1.0E-18)THEN
        NENT(I)=0
        MENT(I,I)=M(I)
        QENT(I,I)=qconv(NK,ithread)-EP(I)*CLW(I)
        ELIJ(I,I)=CLW(I)
        SIJ(I,I)=1.0
      END IF
    END IF
  END DO

  !
  !   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
  !   ***             DOWNDRAFT CALCULATION                      ***
  !
  if (EP(INB).ge.0.0001) then
    !
    !   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
    !   ***                AND CONDENSED WATER FLUX                    ***
    !
    JTT=2
    !
    !    ***                    BEGIN DOWNDRAFT LOOP                    ***
    !
    DO I=INB,1,-1
    !
    !    ***              CALCULATE DETRAINED PRECIPITATION             ***
    !
      WDTRAIN=G*EP(I)*M(I)*CLW(I)
      IF(I.GT.1)THEN
        DO J=1,I-1
          AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
          AWAT=MAX(0.0,AWAT)
          WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
        END DO
      END IF
    !
    !    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
    !    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
    !
    !
    !  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
    !
      COEFF=COEFFS
      WT(I)=OMTSNOW
    !
    !  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
    !
      IF(tconv(I,ithread).GT.273.0)THEN
        COEFF=COEFFR
        WT(I)=OMTRAIN
      END IF
      QSM=0.5*(qconv(I,ithread)+QP(I+1))
      AFAC=COEFF*phconv_hpa(I,ithread)*(qsconv(I,ithread)-QSM)/ &
           (1.0E4+2.0E3*phconv_hpa(I,ithread)*qsconv(I,ithread))
      AFAC=MAX(AFAC,0.0)
      SIGT=SIGP(I)
      SIGT=MAX(0.0,SIGT)
      SIGT=MIN(1.0,SIGT)
      B6=100.*(phconv_hpa(I,ithread)-phconv_hpa(I+1,ithread))*SIGT*AFAC/WT(I)
      C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
      REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
      EVAP(I)=SIGT*AFAC*REVAP
      WATER(I)=REVAP*REVAP
    !
    !    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
    !    ***              HYDROSTATIC APPROXIMATION                 ***
    !
      if (.not. I.eq.1) then
        DHDP=(H(I)-H(I-1))/(pconv_hpa(I-1,ithread)-pconv_hpa(I,ithread))
        DHDP=MAX(DHDP,10.0)
        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
        MP(I)=MAX(MP(I),0.0)
      !
      !   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
      !
        FAC=20.0/(phconv_hpa(I-1,ithread)-phconv_hpa(I,ithread))
        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
      !
      !    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
      !    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
      !
        IF(pconv_hpa(I,ithread).GT.(0.949*pconv_hpa(1,ithread)))THEN
          JTT=MAX(JTT,I)
          MP(I)=MP(JTT)*(pconv_hpa(1,ithread)-pconv_hpa(I,ithread))/(pconv_hpa(1,ithread)- &
              pconv_hpa(JTT,ithread))
        END IF
      endif
    !
    !    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
    !
      if (.not. I.eq.INB) then
        IF(I.EQ.1)THEN
          QSTM=qsconv(1,ithread)
        ELSE
          QSTM=qsconv(I-1,ithread)
        END IF
        IF(MP(I).GT.MP(I+1))THEN
          RAT=MP(I+1)/MP(I)
          QP(I)=QP(I+1)*RAT+qconv(I,ithread)*(1.0-RAT)+100.*GINV* &
               SIGD*(phconv_hpa(I,ithread)-phconv_hpa(I+1,ithread))*(EVAP(I)/MP(I))
         ELSE
          IF(MP(I+1).GT.0.0)THEN
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+tconv(I+1,ithread)*(CL-CPD))+ &
              CPD*(tconv(I+1,ithread)-tconv(I,ithread)))/ &
              (LV(I)+tconv(I,ithread)*(CL-CPD))
          END IF
        END IF
        QP(I)=MIN(QP(I),QSTM)
        QP(I)=MAX(QP(I),0.0)
      endif
    END DO
    !
    !   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
    !
    PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
    !
  endif ! Downdraft calculation
  !
  !   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
  !   ***                    WATER VAPOR FLUCTUATIONS                      ***
  !
  WD=BETA*ABS(MP(ICB))*0.01*RD*tconv(ICB,ithread)/(SIGD*pconv_hpa(ICB,ithread))
  QPRIME=0.5*(QP(1)-qconv(1,ithread))
  TPRIME=LV0*QPRIME/CPD
  !
  !   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
  !   ***                      AND MIXING RATIO                        ***
  !

  DPINV=0.01/(phconv_hpa(1,ithread)-phconv_hpa(2,ithread))
  AM=0.0
  IF(NK.EQ.1)THEN
    AM = sum(M(2:INB))
  END IF
  ! save saturated upward mass flux for first level
  FUP(1)=AM
  IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4
  FT(1)=FT(1)+G*DPINV*AM*(tconv(2,ithread)-tconv(1,ithread)+(GZ(2)-GZ(1))/CPN(1))
  FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
  FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(tconv(2,ithread)- &
       tconv(1,ithread))*DPINV/CPN(1)
  FQ(1)=FQ(1)+G*MP(2)*(QP(2)-qconv(1,ithread))* &
       DPINV+SIGD*EVAP(1)
  FQ(1)=FQ(1)+G*AM*(qconv(2,ithread)-qconv(1,ithread))*DPINV

  FQ(1)=FQ(1)+G*DPINV*sum(MENT(2:INB,1)*(QENT(2:INB,1)-qconv(1,ithread)))
  !
  !   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
  !   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
  !
  !   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
  !   ***                      THROUGH EACH LEVEL                          ***
  !
  DO I=2,INB
    DPINV=0.01/(phconv_hpa(I,ithread)-phconv_hpa(I+1,ithread))
    CPINV=1.0/CPN(I)
    AMP1=0.0
    AD=0.0
    IF(I.GE.NK)THEN
      AMP1 = sum(M(I+1:INB+1))
    END IF
    AMP1 = AMP1 + sum(MENT(1:I,I+1:INB+1))
  ! save saturated upward mass flux
    FUP(I)=AMP1
    IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4

    AD = sum(MENT(I:INB,1:I-1))
  ! save saturated downward mass flux
    FDOWN(I)=AD
    FT(I)=FT(I)+G*DPINV*(AMP1*(tconv(I+1,ithread)-tconv(I,ithread)+(GZ(I+1)-GZ(I))* &
         CPINV)-AD*(tconv(I,ithread)-tconv(I-1,ithread)+(GZ(I)-GZ(I-1))*CPINV)) &
         -SIGD*LVCP(I)*EVAP(I)
    FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+ &
         tconv(I,ithread)*(CPV-CPD)*(qconv(I,ithread)-QENT(I,I)))*CPINV
    FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)* &
         (tconv(I+1,ithread)-tconv(I,ithread))*DPINV*CPINV
    FQ(I)=FQ(I)+G*DPINV*(AMP1*(qconv(I+1,ithread)-qconv(I,ithread))- &
         AD*(qconv(I,ithread)-qconv(I-1,ithread)))
    DO K=1,I-1
      AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
      AWAT=MAX(AWAT,0.0)
      FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-qconv(I,ithread))
    END DO

    FQ(I)=FQ(I)+G*DPINV*sum(MENT(I:INB,I)*(QENT(I:INB,I)-qconv(I,ithread)))
    FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)* &
         (QP(I+1)-qconv(I,ithread))-MP(I)*(QP(I)-qconv(I-1,ithread)))*DPINV
  END DO
  !
  !   *** Adjust tendencies at top of convection layer to reflect  ***
  !   ***       actual position of the level zero CAPE             ***
  !
  FQOLD=FQ(INB)
  FQ(INB)=FQ(INB)*(1.-FRAC)
  FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((phconv_hpa(INB,ithread)- &
       phconv_hpa(INB+1,ithread))/ &
       (phconv_hpa(INB-1,ithread)-phconv_hpa(INB,ithread)))*LV(INB)/LV(INB-1)
  FTOLD=FT(INB)
  FT(INB)=FT(INB)*(1.-FRAC)
  FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((phconv_hpa(INB,ithread)- &
       phconv_hpa(INB+1,ithread))/ &
       (phconv_hpa(INB-1,ithread)-phconv_hpa(INB,ithread)))*CPN(INB)/CPN(INB-1)
!
!   ***   Very slightly adjust tendencies to force exact   ***
!   ***     enthalpy, momentum and tracer conservation     ***
!
  ENTS=0.0

  ENTS = sum((CPN(1:INB)*FT(1:INB)+LV(1:INB)*FQ(1:INB))* &
        (phconv_hpa(1:INB,ithread)-phconv_hpa(2:INB+1,ithread)))

  ENTS=ENTS/(phconv_hpa(1,ithread)-phconv_hpa(INB+1,ithread))

  FT(1:INB)=FT(1:INB) - ENTS/CPN(1:INB)

  ! ************************************************
  ! **** DETERMINE MASS DISPLACEMENT MATRIX
  ! ***** AND COMPENSATING SUBSIDENCE
  ! ************************************************

  ! mass displacement matrix due to saturated up-and downdrafts
  ! inside the cloud and determine compensating subsidence
  ! FUP(I) (saturated updrafts), FDOWN(I) (saturated downdrafts) are assumed to be
  ! balanced by  compensating subsidence (SUB(I))
  ! FDOWN(I) and SUB(I) defined positive downwards

  ! nconvtop IS THE TOP LEVEL AT WHICH CONVECTIVE MASS FLUXES ARE DIAGNOSED
  ! EPSILON IS A SMALL NUMBER

  fmass(NK, :INB+1,ithread) = fmass(NK,:INB+1,ithread)+M(:INB+1)
  fmass(:INB+1,:INB+1,ithread) = fmass(:INB+1,:INB+1,ithread)+MENT(:INB+1,:INB+1)
  sub(1,ithread) = 0.
  sub(2:INB+1,ithread) = FUP(1:INB) - FDOWN(2:INB+1)
  nconvtop=1
  do i=1,INB+1
    do j=1,INB+1
      if (fmass(j,i,ithread).gt.EPSILON) nconvtop=MAX(nconvtop,i,j)
    end do
  end do
  nconvtop=nconvtop+1
  RETURN
  !
END SUBROUTINE CONVECT
!
! ---------------------------------------------------------------------------
!
SUBROUTINE TLIFT(GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK, ithread)
  !
  !-cv
  use par_mod

  implicit none
  !-cv
  !====>Begin Module TLIFT      File convect.f      Undeclared variables
  !
  !Argument variables
  !
  integer, intent(in) :: ithread
  integer :: icb, kk, nd, nk, nl
  !
  !Local variables
  !
  integer :: i, j, nsb, nst
  !
  real :: ah0, ahg, alv, cpinv, cpp, denom
  real :: es, qg, rg, s, tc, tg
  !
  !====>End Module   TLIFT      File convect.f

  REAL :: GZ(ND),TPK(ND),CLW(ND)
  REAL :: TVP(ND)
  !
  !   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
  !
  REAL,PARAMETER :: CPD=1005.7
  REAL,PARAMETER :: CPV=1870.0
  REAL,PARAMETER :: CL=2500.0
  REAL,PARAMETER :: RV=461.5
  REAL,PARAMETER :: RD=287.04
  REAL,PARAMETER :: LV0=2.501E6
  !
  REAL,PARAMETER :: CPVMCL=CL-CPV
  REAL,PARAMETER :: EPS0=RD/RV
  REAL,PARAMETER :: EPSI=1./EPS0
  !
  !   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
  !
  AH0=(CPD*(1.-qconv(NK,ithread))+CL*qconv(NK,ithread))*tconv(NK,ithread)+ &
    qconv(NK,ithread)*(LV0-CPVMCL*(tconv(NK,ithread)-273.15))+GZ(NK)
  CPP=CPD*(1.-qconv(NK,ithread))+qconv(NK,ithread)*CPV
  CPINV=1./CPP
  !
  IF(KK.EQ.1)THEN
  !
  !   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
  !
    CLW(1:ICB-1) = 0.0
    TPK(NK:ICB-1)=tconv(NK,ithread)-(GZ(NK:ICB-1)-GZ(NK))*CPINV
    TVP(NK:ICB-1)=TPK(NK:ICB-1)*(1.+qconv(NK,ithread)*EPSI)
  END IF
  !
  !    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
  !
  NST=ICB
  NSB=ICB
  IF(KK.EQ.2)THEN
    NST=NL
    NSB=ICB+1
  END IF
  DO I=NSB,NST
    TG=tconv(I,ithread)
    QG=qsconv(I,ithread)
    ALV=LV0-CPVMCL*(tconv(I,ithread)-273.15)
    DO J=1,2
      S=CPD+ALV*ALV*QG/(RV*tconv(I,ithread)*tconv(I,ithread))
      S=1./S
      AHG=CPD*TG+(CL-CPD)*qconv(NK,ithread)*tconv(I,ithread)+ALV*QG+GZ(I)
      TG=TG+S*(AH0-AHG)
      TG=MAX(TG,35.0)
      TC=TG-273.15
      DENOM=243.5+TC
      IF(TC.GE.0.0)THEN
        ES=6.112*EXP(17.67*TC/DENOM)
      ELSE
        ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
      END IF
      QG=EPS0*ES/(pconv_hpa(I,ithread)-ES*(1.-EPS0))
    END DO
    ALV=LV0-CPVMCL*(tconv(I,ithread)-273.15)
    TPK(I)=(AH0-(CL-CPD)*qconv(NK,ithread)*tconv(I,ithread)-GZ(I)-ALV*QG)/CPD
    CLW(I)=qconv(NK,ithread)-QG
    CLW(I)=MAX(0.0,CLW(I))
    RG=QG/(1.-qconv(NK,ithread))
    TVP(I)=TPK(I)*(1.+RG*EPSI)
  END DO
  RETURN
END SUBROUTINE TLIFT

end module conv_mod
