!*****************************************************************************
!                                                                            *
! L. Bakels 2021: This module contains dry deposition related subroutines    *
!                                                                            *
! To do: dry deposition mass loss substraction of individual particles       *
!        should be moved from to timemanager_mod.f90 to here                 *
!                                                                            *
!                                                                            *
!*****************************************************************************

module drydepo_mod
  use par_mod
  use com_mod
  use unc_mod
  use windfields_mod
  use erf_mod
  
  implicit none

  real,allocatable,dimension(:,:,:) :: xlanduse 
   ! area fractions in percent [0-1]
  real,allocatable,dimension(:,:,:,:) :: xlandusen
   ! nested area fractions in percent [0-1]
  real,allocatable,dimension(:,:,:,:) :: vdep ! deposition velocity [m/s]

  ! roughness length
  real,allocatable,dimension(:,:) :: z0_drydep
  ! roughtness lenght nested area
  real,allocatable,dimension(:,:,:) :: z0_drydepn

contains

subroutine alloc_drydepo

  implicit none
  integer:: stat

  if (.not. drydep) return
  write(*,*) 'allocate drydepo fields'
  allocate(xlanduse(0:nxmax-1,0:nymax-1,numclass),      &
           vdep(0:nxmax-1,0:nymax-1,maxspec,numwfmem), &
           z0_drydep(0:nxmax-1,0:nymax-1), stat=stat)
  if (stat.ne.0) error stop "Could not allocate drydepo fields"
  if (numbnests.ge.1) then
    allocate(xlandusen(0:nxmaxn-1,0:nymaxn-1,numclass,numbnests), &
           vdepn(0:nxmaxn-1,0:nymaxn-1,maxspec,numwfmem,numbnests), &
           z0_drydepn(0:nxmaxn-1,0:nymaxn-1,numbnests), stat=stat)
    if (stat.ne.0) error stop "Could not allocate drydepo fields"
  endif  
           
end subroutine alloc_drydepo

subroutine dealloc_drydepo

  if (.not. drydep) return
  deallocate(xlanduse,vdep,z0_drydep)
  if (numbnests.ge.1) then
    deallocate(xlandusen,vdepn,z0_drydepn)
  endif

end subroutine dealloc_drydepo

subroutine assignland

  !*****************************************************************************
  !                                                                            *
  !     This routine assigns fractions of the 13 landuse classes to each ECMWF *
  !     grid point.                                                            *
  !     The landuse inventory of                                               *
  !                                                                            *
  ! Belward, A.S., Estes, J.E., and Kline, K.D., 1999,                         *
  ! The IGBP-DIS 1-Km Land-Cover Data Set DISCover:                            *
  ! A Project Overview: Photogrammetric Engineering and Remote Sensing ,       *
  ! v. 65, no. 9, p. 1013-1020                                                 *
  !                                                                            *
  !     if there are no data in the inventory                                  *
  !     the ECMWF land/sea mask is used to distinguish                         *
  !     between sea (-> ocean) and land (-> grasslands).                       *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     5 December 1996                                                        *
  !     8 February 1999 Additional use of nests, A. Stohl                      *
  !    29 December 2006 new landuse inventory, S. Eckhardt                     *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! xlanduse          fractions of numclass landuses for each model grid point *
  ! landinvent       landuse inventory (0.3 deg resolution)                    *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: ix,jy,k,l,li,nrefine,iix,jjy,stat
  integer,parameter :: lumaxx=1200,lumaxy=600
  integer,parameter :: xlon0lu=-180,ylat0lu=-90
  real,parameter :: dxlu=0.3
  real :: xlon,ylat,sumperc,p,xi,yj
  real,allocatable,dimension(:,:,:) :: xlandusep
  ! character*2 ck

  if (.not.DRYDEP) return
  
  allocate( xlandusep(lumaxx,lumaxy,numclass), stat=stat)
  if (stat.ne.0) error stop "Could not allocate xlandusep in assignland"

  do ix=1,lumaxx
    do jy=1,lumaxy
       do k=1,numclass
         xlandusep(ix,jy,k)=0.
       end do
       sumperc=0.
       do li=1,3
         sumperc=sumperc+landinvent(ix,jy,li+3)
       end do
       do li=1,3
         k=landinvent(ix,jy,li)
       if (sumperc.gt.0) then
         p=landinvent(ix,jy,li+3)/sumperc
       else
         p=0
       endif
! p has values between 0 and 1
         xlandusep(ix,jy,k)=p
       end do
    end do
  end do

  ! do 13 k=1,11
  ! write (ck,'(i2.2)') k
  ! open(4,file='xlandusetest'//ck,form='formatted')
  ! do 11 ix=1,lumaxx
  !11       write (4,*) (xlandusep(ix,jy,k),jy=1,lumaxy)
  !11       write (4,*) (landinvent(ix,jy,k),jy=1,lumaxy)
  !13     close(4)

  ! write (*,*) xlon0,ylat0,xlon0n(1),ylat0n(1),nxmin1,nymin1
  ! write (*,*) dx, dy, dxout, dyout, ylat0, xlon0
  nrefine=10
  do ix=0,nxmin1
    do jy=0,nymin1
      do k=1,numclass
        sumperc=0.
        xlanduse(ix,jy,k)=0.
      end do
        do iix=1, nrefine
          xlon=(ix+(iix-1)/real(nrefine))*dx+xlon0        ! longitude, should be between -180 and 179
          if (xlon.ge.(xlon0lu+lumaxx*dxlu))  then
               xlon=xlon-lumaxx*dxlu
          endif
          do jjy=1, nrefine
           ylat=(jy+(jjy-1)/real(nrefine))*dy+ylat0       ! and lat. of each gridpoint
           xi=int((xlon-xlon0lu)/dxlu)+1
           yj=int((ylat-ylat0lu)/dxlu)+1
           if (xi.gt.lumaxx) xi=xi-lumaxx
           if (yj.gt.lumaxy) yj=yj-lumaxy
           if (xi.lt.0) then
              write (*,*) 'problem with landuseinv sampling: ', &
                   xlon,xlon0lu,ix,iix,xlon0,dx,nxmax
              error stop
           endif
           do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)+xlandusep(int(xi),int(yj),k)
             sumperc=sumperc+xlanduse(ix,jy,k)  ! just for the check if landuseinv. is available
           end do
          end do
        end do
          if (sumperc.gt.0) then                       ! detailed landuse available
          sumperc=0.
          do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)/real(nrefine*nrefine)
            sumperc=sumperc+xlanduse(ix,jy,k)
          end do
  !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
          if (sumperc.lt.1-1E-5) then
            do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)/sumperc
            end do
          endif
          else
            if (lsm(ix,jy).lt.0.1) then           ! over sea  -> ocean
              xlanduse(ix,jy,3)=1.
            else                                  ! over land -> rangeland
              xlanduse(ix,jy,7)=1.
            endif
          endif


    end do
  end do

  !***********************************
  ! for test: write out xlanduse

  ! open(4,file='landusetest',form='formatted')
  ! do 56 k=1,13
  ! do 55 ix=0,nxmin1
  !55    write (4,*) (xlanduse(ix,jy,k),jy=0,nymin1)
  !56    continue
  ! close(4)
  ! write (*,*) 'landuse written'
  !stop
  !  open(4,file='landseatest'//ck,form='formatted')
  ! do 57 ix=0,nxmin1
  !57       write (4,*) (lsm(ix,jy),jy=0,nymin1)
  !  write (*,*) 'landseamask written'

  !****************************************
  ! Same as above, but for the nested grids
  !****************************************

  !************** TEST ********************
  ! dyn(1)=dyn(1)/40
  ! dxn(1)=dxn(1)/40
  ! xlon0n(1)=1
  ! ylat0n(1)=50
  !************** TEST ********************

  do l=1,numbnests
    do ix=0,nxn(l)-1
      do jy=0,nyn(l)-1
        do k=1,numclass
          sumperc=0.
          xlandusen(ix,jy,k,l)=0.
        end do
          do iix=1, nrefine
           xlon=(ix+(iix-1)/real(nrefine))*dxn(l)+xlon0n(l)
           do jjy=1, nrefine
             ylat=(jy+(jjy-1)/real(nrefine))*dyn(l)+ylat0n(l)
             xi=int((xlon-xlon0lu)/dxlu)+1
             yj=int((ylat-ylat0lu)/dxlu)+1
             if (xi.gt.lumaxx) xi=xi-lumaxx
             if (yj.gt.lumaxy) yj=yj-lumaxy
             do k=1,numclass
                xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)+ &
                     xlandusep(int(xi),int(yj),k)
               sumperc=sumperc+xlandusen(ix,jy,k,l)
             end do
           end do
          end do
          if (sumperc.gt.0) then                     ! detailed landuse available
          sumperc=0.
            do k=1,numclass
               xlandusen(ix,jy,k,l)= &
                    xlandusen(ix,jy,k,l)/real(nrefine*nrefine)
              sumperc=sumperc+xlandusen(ix,jy,k,l)
            end do
  !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
            if (sumperc.lt.1-1E-5) then
              do k=1,numclass
                xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)/sumperc
              end do
            endif
          else                                    ! check land/sea mask
            if (lsmn(ix,jy,l).lt.0.1) then   ! over sea  -> ocean
              xlandusen(ix,jy,3,l)=1.
            else                        ! over land -> grasslands
              xlandusen(ix,jy,7,l)=1.
            endif
          endif
      end do
    end do
  end do

  !***********************************
  ! for test: write out xlanduse

  ! do 66 k=1,11
  ! write (ck,'(i2.2)') k
  ! open(4,file='nlandusetest'//ck,form='formatted')
  ! do 65 ix=0,nxn(1)-1
  !65       write (4,*) (xlandusen(ix,jy,k,1),jy=0,nyn(1)-1)
  !66      close(4)

  ! write (*,*) 'landuse nested written'
end subroutine assignland

real function raerod (l,ust,z0)
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the aerodynamical resistance ra from ground up to href  *
  !                                                                            *
  !     AUTHOR: Matthias Langer, modified by Andreas Stohl (6 August 1993)     *
  !                                                                            *
  !     Literature:                                                            *
  !     [1]  Hicks/Baldocchi/Meyers/Hosker/Matt (1987), A Preliminary          *
  !             Multiple Resistance Routine for Deriving Dry Deposition        *
  !             Velocities from Measured Quantities.                           *
  !             Water, Air and Soil Pollution 36 (1987), pp.311-330.           *
  !     [2]  Scire/Yamartino/Carmichael/Chang (1989),                          *
  !             CALGRID: A Mesoscale Photochemical Grid Model.                 *
  !             Vol II: User's Guide. (Report No.A049-1, June, 1989)           *
  !                                                                            *
  !     Variable list:                                                         *
  !     L     = Monin-Obukhov-length [m]                                       *
  !     ust   = friction velocity [m/sec]                                      *
  !     z0    = surface roughness length [m]                                   *
  !     href  = reference height [m], for which deposition velocity is         *
  !             calculated                                                     *
  !                                                                            *
  !     Constants:                                                             *
  !     karman    = von Karman-constant (~0.4)                                 *
  !     ramin = minimum resistence of ra (1 s/m)                               *
  !                                                                            *
  !     Subprograms and functions:                                             *
  !     function psih (z/L)                                                    *
  !                                                                            *
  !*****************************************************************************

  use pbl_profile_mod, only: psih

  implicit none

  real :: l,ust,z0

  raerod=(alog(href/z0)-psih(href,l)+psih(z0,l))/(karman*ust)

end function raerod

subroutine drydepo_massloss(ipart,ks,ldeltat,drydepopart)
  use particle_mod

  implicit none

  integer,intent(in) ::  &
    ipart,               & ! particle index
    ks,                  & ! species index
    ldeltat                ! radioactive decay time
  real(dep_prec),intent(out) ::  &
    drydepopart            ! drydeposit for particle ipart
  real decfact             ! radioactive decay factor

  if (decay(ks).gt.0.) then             ! radioactive decay
    decfact=exp(-real(abs(lsynctime))*decay(ks))
  else
    decfact=1.
  endif
  drydepopart=mass(ipart,ks)*prob(ipart,ks)*decfact
  
  drydeposit(ipart,ks)=drydeposit(ipart,ks)+ &
    mass(ipart,ks)*prob(ipart,ks)*decfact

  mass(ipart,ks)=mass(ipart,ks)*(1.-prob(ipart,ks))*decfact

  if (decay(ks).gt.0.) then   ! correct for decay (see wetdepo)
    drydepopart=drydepopart*exp(real(abs(ldeltat))*decay(ks))
  endif

end subroutine drydepo_massloss

subroutine drydepokernel(nunc,deposit,x,y,nage,kp,thread)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition to the grid using a uniform kernel with  *
  !     bandwidths dx and dy.                                                  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************
  ! Changes:
  ! eso 10/2016: Added option to disregard kernel 
  ! 
  !*****************************************************************************

  implicit none

  integer,intent(in) :: thread
  real(dep_prec), dimension(maxspec) :: deposit
  real :: x,y,ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,nunc,nage,kp


  xl=(x*dx+xoutshift)/dxout
  yl=(y*dy+youtshift)/dyout
  ix=int(xl)
  jy=int(yl)
  ddx=xl-real(ix)                   ! distance to left cell border
  ddy=yl-real(jy)                   ! distance to lower cell border

  if (ddx.gt.0.5) then
    ixp=ix+1
    wx=1.5-ddx
  else
    ixp=ix-1
    wx=0.5+ddx
  endif

  if (ddy.gt.0.5) then
    jyp=jy+1
    wy=1.5-ddy
  else
    jyp=jy-1
    wy=0.5+ddy
  endif

  ! If no kernel is used, direct attribution to grid cell
  !******************************************************

  if (.not.lusekerneloutput) then
    do ks=1,nspec
      if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then
        if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
             (jy.le.numygrid-1)) then
#ifdef _OPENMP
          gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
               gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)
#else
          drygridunc(ix,jy,ks,kp,nunc,nage)= &
               drygridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
#endif
        end if
      end if
    end do
  else ! use kernel 


  ! Determine mass fractions for four grid points
  !**********************************************
    do ks=1,nspec

     if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then

        if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
          (jy.le.numygrid-1)) then
          w=wx*wy
#ifdef _OPENMP
          gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
             gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
          drygridunc(ix,jy,ks,kp,nunc,nage)= &
             drygridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
       endif

      if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
         (jyp.le.numygrid-1)) then
        w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
        gridunc_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)= &
             gridunc_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygridunc(ixp,jyp,ks,kp,nunc,nage)= &
             drygridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
         (jy.le.numygrid-1)) then
        w=(1.-wx)*wy
#ifdef _OPENMP
        gridunc_omp(ixp,jy,1,ks,kp,nunc,nage,thread)= &
             gridunc_omp(ixp,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygridunc(ixp,jy,ks,kp,nunc,nage)= &
             drygridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
         (jyp.le.numygrid-1)) then
        w=wx*(1.-wy)
#ifdef _OPENMP
        gridunc_omp(ix,jyp,1,ks,kp,nunc,nage,thread)= &
             gridunc_omp(ix,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygridunc(ix,jyp,ks,kp,nunc,nage)= &
             drygridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      endif ! deposit>0
    end do
  end if
end subroutine drydepokernel

subroutine drydepokernel_nest(nunc,deposit,x,y,nage,kp,thread)
  !                               i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     nested deposition fields using a uniform kernel with bandwidths        *
  !     dxoutn and dyoutn.                                                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !      2 September 2004: Adaptation from drydepokernel.                      *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer,intent(in) :: thread
  real(dep_prec), dimension(maxspec) :: deposit
  real :: x,y,ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,kp,nunc,nage



  xl=(x*dx+xoutshiftn)/dxoutn
  yl=(y*dy+youtshiftn)/dyoutn
  ix=int(xl)
  jy=int(yl)
  ddx=xl-real(ix)                   ! distance to left cell border
  ddy=yl-real(jy)                   ! distance to lower cell border

  if (ddx.gt.0.5) then
    ixp=ix+1
    wx=1.5-ddx
  else
    ixp=ix-1
    wx=0.5+ddx
  endif

  if (ddy.gt.0.5) then
    jyp=jy+1
    wy=1.5-ddy
  else
    jyp=jy-1
    wy=0.5+ddy
  endif


  ! Determine mass fractions for four grid points
  !**********************************************
  do ks=1,nspec

    if (DRYDEPSPEC(ks).and.(abs(deposit(ks)).gt.0)) then

      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
           (jy.le.numygridn-1)) then
        w=wx*wy
#ifdef _OPENMP
        griduncn_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
             griduncn_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygriduncn(ix,jy,ks,kp,nunc,nage)= &
             drygriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgridn-1).and. &
           (jyp.le.numygridn-1)) then
        w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
        griduncn_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)= &
             griduncn_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygriduncn(ixp,jyp,ks,kp,nunc,nage)= &
             drygriduncn(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgridn-1).and. &
           (jy.le.numygridn-1)) then
        w=(1.-wx)*wy
#ifdef _OPENMP
        griduncn_omp(ixp,jy,1,ks,kp,nunc,nage,thread)= &
             griduncn_omp(ixp,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygriduncn(ixp,jy,ks,kp,nunc,nage)= &
             drygriduncn(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

      if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgridn-1).and. &
           (jyp.le.numygridn-1)) then
        w=wx*(1.-wy)
#ifdef _OPENMP
        griduncn_omp(ix,jyp,1,ks,kp,nunc,nage,thread)= &
             griduncn_omp(ix,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
        drygriduncn(ix,jyp,ks,kp,nunc,nage)= &
             drygriduncn(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
      endif

    endif

  end do
end subroutine drydepokernel_nest

subroutine part0(dquer,dsigma,density,ni,fract,schmi,cun,vsh)
  !                  i      i       i   i   o     o    o   o
  !*****************************************************************************
  !                                                                            *
  !      Calculation of time independent factors of the dry deposition of      *
  !      particles:                                                            *
  !      Log-Normal-distribution of mass [dM/dlog(dp)], unimodal               *
  !                                                                            *
  !      AUTHOR: Matthias Langer, adapted by Andreas Stohl, 13 November 1993   *
  !                                                                            *
  !      Literature:                                                           *
  !      [1]  Scire/Yamartino/Carmichael/Chang (1989),                         *
  !             CALGRID: A Mesoscale Photochemical Grid Model.                 *
  !             Vol II: User's Guide. (Report No.A049-1, June, 1989)           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! alpha            help variable                                             *
  ! cun              'slip-flow' correction after Cunningham                   *
  ! d01 [um]         upper diameter                                            *
  ! d02 [um]         lower diameter                                            *
  ! dc [m2/s]        coefficient of Brownian diffusion                         *
  ! delta            distance given in standard deviation units                *
  ! density [kg/m3]  density of the particle                                   *
  ! dmean            geometric mean diameter of interval                       *
  ! dquer [um]       geometric mass mean particle diameter                     *
  ! dsigma           e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass   *
  !                  are between 0.1*dquer and 10*dquer                        *
  ! fract(ni)        mass fraction of each diameter interval                   *
  ! kn               Knudsen number                                            *
  ! ni               number of diameter intervals, for which deposition        *
  !                  is calculated                                             *
  ! schmidt          Schmidt number                                            *
  ! schmi            schmidt**2/3                                              *
  ! vsh [m/s]        gravitational settling velocity of the particle           *
  ! x01              normalized upper diameter                                 *
  ! x02              normalized lower diameter                                 *
  !                                                                            *
  ! Constants:                                                                 *
  ! g [m/s2]         Acceleration of gravity                                   *
  ! kb [J/K]         Stefan-Boltzmann constant                                 *
  ! lam [m]          mean free path of air molecules                           *
  ! myl [kg/m/s]     dynamical viscosity of air                                *
  ! nyl [m2/s]       kinematic viscosity of air                                *
  ! tr               reference temperature                                     *
  !                                                                            *
  ! Function:                                                                  *
  ! erf              calculates the integral of the Gauss function             *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real,parameter :: tr=293.15

  integer :: i,ni
  real :: dquer,dsigma,density,xdummy,d01,d02,delta,x01,x02
  real :: dmean,alpha,cun,dc,schmidt,kn,erf,fract_norm
  real,dimension(ni),intent(inout) :: fract,schmi,vsh
  real,parameter :: myl=1.81e-5,nyl=0.15e-4
  real,parameter :: lam=6.53e-8,kb=1.38e-23,eps=1.2e-38

  ! xdummy constant for all intervals
  !**********************************

  xdummy=sqrt(2.)*alog(dsigma)


  ! particles diameters are split up to ni intervals between
  ! dquer-3*dsigma and dquer+3*dsigma
  !*********************************************************
  ! Normalisation. Why was it not normalised?
  !******************************************
  x01=alog(dsigma**3)/xdummy
  x02=alog(dsigma**(-3))/xdummy
  fract_norm=0.5*(erf(x01)-erf(x02))

  delta=6./real(ni)

  d01=dquer*dsigma**(-3)
  do i=1,ni
    d02=d01
    d01=dquer*dsigma**(-3.+delta*real(i))
    x01=alog(d01/dquer)/xdummy
    x02=alog(d02/dquer)/xdummy
    !print*,'part0:: d02=' , d02 , 'd01=', d01

  ! Area under Gauss-function is calculated and gives mass fraction of interval
  !****************************************************************************

    fract(i)=0.5*(erf(x01)-erf(x02))/fract_norm
    !print*,'part0:: fract(',i,')', fract(i)
    !print*,'part0:: fract', fract(i), x01, x02, erf(x01), erf(x02)

  ! Geometric mean diameter of interval in [m]
  !*******************************************

    dmean=1.E-6*exp(0.5*alog(d01*d02))
    !print*,'part0:: dmean=', dmean

  ! Calculation of time independent parameters of each interval
  !************************************************************

    kn=2.*lam/dmean
    if ((-1.1/kn).le.log10(eps)*log(10.)) then
      alpha=1.257
    else
      alpha=1.257+0.4*exp(-1.1/kn)
    endif
    cun=1.+alpha*kn
    dc=kb*tr*cun/(3.*pi*myl*dmean)
    schmidt=nyl/dc
    schmi(i)=schmidt**(-2./3.)
    vsh(i)=ga*density*dmean*dmean*cun/(18.*myl)

    !print*,'part0:: vsh(',i,')', vsh(i)

  end do

  !stop 'part0' 
end subroutine part0

subroutine get_vdep_prob(itime,xt,yt,zt,tmpprob,ithread)
  !                    i     i  i  i  o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of the probability for dry deposition                         * 
  !                                                                            *
  !  Particle positions are read in - prob returned                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          time at which this subroutine is entered                *
  ! itimec [s]         actual time, which is incremented in this subroutine    *
  ! href [m]           height for which dry deposition velocity is calculated  *
  ! ldirect            1 forward, -1 backward                                  *
  ! ldt [s]            Time step for the next integration                      *
  ! lsynctime [s]      Synchronisation interval of FLEXPART                    *
  ! ngrid              index which grid is to be used                          *
  ! prob               probability of absorption due to dry deposition         *
  ! vdepo              Deposition velocities for all species                   *
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use interpol_mod

  implicit none

  real,intent(in) :: xt,yt,zt
  integer,intent(in) :: itime,ithread !ithread starting at 1
  real,intent(out) :: tmpprob(maxspec)
  integer :: ks,m,memindnext!nix,njy,
  real :: vdepo(maxspec),vdeptemp(2)
  real :: eps

  eps=nxmax/3.e5

  if (DRYDEP) then    ! reset probability for deposition
    do ks=1,nspec
      depoindicator(ks,ithread)=.true.
      tmpprob(ks)=0.
    end do
  endif


  ! Determine whether lat/long grid or polarstereographic projection
  ! is to be used
  ! Furthermore, determine which nesting level to be used
  !*****************************************************************
  call find_ngrid(xt,yt)

  !***************************
  ! Interpolate necessary data
  !***************************

  if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
    memindnext=1
  else
    memindnext=2
  endif

  ! Determine nested grid coordinates
  !**********************************
  call find_grid_indices(xt,yt)

  ! Determine probability of deposition
  !************************************

  if ((DRYDEP).and.(real(zt).lt.2.*href)) then
    do ks=1,nspec
      if (DRYDEPSPEC(ks)) then
        if (depoindicator(ks,ithread)) then
          if (ngrid.le.0) then
            do m=1,2
              call hor_interpol(vdep,vdeptemp(m),ks,memind(m),maxspec)
            end do
          else
            do m=1,2
              call hor_interpol_nest(vdepn,vdeptemp(m),ks,memind(m),maxspec)
            end do
          endif
          call temporal_interpolation(vdeptemp(1),vdeptemp(2),vdepo(ks))
        endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
        tmpprob(ks)=vdepo(ks)
  !               prob(ks)=vdepo(ks)/2./href
  ! instead of prob - return vdepo -> result kg/m2/s
      endif
    end do
  endif
end subroutine get_vdep_prob

subroutine drydepo_probability(ipart,dt,zts,vdepo,ithread)
  use par_mod
  use com_mod
  use interpol_mod
  use particle_mod

  implicit none

  integer,intent(in) :: ithread ! OMP thread starting at 1
  integer,intent(in) :: ipart ! particle index
  real,intent(out) :: vdepo(maxspec)  ! deposition velocities for all species
  real,intent(in) :: dt,zts             ! real(ldt), real(zt)
  integer :: ns,m                      ! loop variable over species
  real :: vdeptemp(2)

  if (zts.lt.2.*href) then
    do ns=1,nspec
      if (DRYDEPSPEC(ns)) then
        if (depoindicator(ns,ithread)) then
          if (ngrid.le.0) then
            do m=1,2
              call hor_interpol(vdep,vdeptemp(m),ns,memind(m),maxspec)
            end do
          else
            do m=1,2
              call hor_interpol_nest(vdepn,vdeptemp(m),ns,memind(m),maxspec)
            end do
          endif
          call temporal_interpolation(vdeptemp(1),vdeptemp(2),vdepo(ns))
        endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
        prob(ipart,ns)=1.+(prob(ipart,ns)-1.)*exp(-vdepo(ns)*abs(dt)/(2.*href))
        !if (pp.eq.535) write(*,*) 'advance1', ks,dtt,p1,vdep(ix,jy,ks,1)
      endif
    end do
  endif
end subroutine drydepo_probability

subroutine getvdep(n,ix,jy,ust,temp,pa,L,gr,rh,rr,snow,vdepo)
  !                   i i  i   i   i   i  i i  i  i    i o
  !*****************************************************************************
  !                                                                            *
  !  This routine calculates the dry deposition velocities.                    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 December 1996                                                       *
  !     Sabine Eckhardt, Jan 07                                                *
  !     if the latitude is negative: add half a year to the julian day         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! gr [W/m2]         global radiation                                         *
  ! L [m]             Obukhov length                                           *
  ! nyl               kinematic viscosity                                      *
  ! pa [Pa]           surface air pressure                                     *
  ! ra [s/m]          aerodynamic resistance                                   *
  ! raquer [s/m]      average aerodynamic resistance                           *
  ! rh [0-1]          relative humidity                                        *
  ! rhoa              density of the air                                       *
  ! rr [mm/h]         precipitation rate                                       *
  ! temp [K]          2m temperature                                           *
  ! tc [C]            2m temperature                                           *
  ! ust [m/s]         friction velocity                                        *
  ! snow [m of water equivalent] snow depth                                    *
  ! xlanduse          fractions of numclasS landuses for each model grid point *
  !                                                                            *
  !*****************************************************************************
  use date_mod
  
  implicit none

  integer :: yyyymmdd,hhmmss,yyyy,mmdd,n,lseason,i,j,ix,jy
  real :: vdepo(maxspec),vd,rb(maxspec),rc(maxspec),raquer,ylat
  real :: ra,ust,temp,tc,pa,L,gr,rh,rr,myl,nyl,rhoa,diffh2o,snow
  real :: slanduse(numclass)
  real,parameter :: eps=1.e-5
  real(kind=dp) :: jul

  ! Calculate month and determine the seasonal category
  !****************************************************

  jul=bdate+real(wftime(n),kind=dp)/86400._dp

  ylat=jy*dy+ylat0
  if (ylat.lt.0) then
      jul=jul+365.*0.5
  endif


  call caldate(jul,yyyymmdd,hhmmss)
  yyyy=yyyymmdd/10000
  mmdd=yyyymmdd-10000*yyyy

  if ((ylat.gt.-20).and.(ylat.lt.20)) then
     mmdd=600 ! summer
  endif

  if ((mmdd.ge.1201).or.(mmdd.le.301)) then
    lseason=4
  else if ((mmdd.ge.1101).or.(mmdd.le.331)) then
    lseason=3
  else if ((mmdd.ge.401).and.(mmdd.le.515)) then
    lseason=5
  else if ((mmdd.ge.516).and.(mmdd.le.915)) then
    lseason=1
  else
    lseason=2
  endif

  ! Calculate diffusivity of water vapor
  !************************************
  diffh2o=2.11e-5*(temp/273.15)**1.94*(101325/pa)

  ! Conversion of temperature from K to C
  !**************************************

  tc=temp-273.15

  ! Calculate dynamic viscosity
  !****************************

  ! Why is this different from the viscosity funtion???

  if (tc.lt.0) then
    myl=(1.718+0.0049*tc-1.2e-05*tc**2)*1.e-05
  else
    myl=(1.718+0.0049*tc)*1.e-05
  endif

  ! Calculate kinematic viscosity
  !******************************

  rhoa=pa/(287.*temp)
  nyl=myl/rhoa


  ! 0. Set all deposition velocities zero
  !**************************************

  do i=1,nspec
    vdepo(i)=0.
  end do


  ! 1. Compute surface layer resistances rb
  !****************************************

  call getrb(nspec,ust,nyl,diffh2o,reldiff,rb)

  ! change for snow
  do j=1,numclass
    if (snow.gt.0.001) then ! 10 mm
       if (j.eq.12) then
         slanduse(j)=1.
       else
         slanduse(j)=0.
       endif
    else
       slanduse(j)=xlanduse(ix,jy,j)
    endif
  end do

  raquer=0.
  do j=1,numclass            ! loop over all landuse classes

    if (slanduse(j).gt.eps)  then

  ! 2. Calculate aerodynamic resistance ra
  !***************************************
      if (DRYDEP.and.(j.eq.7)) then
        ra=raerod(L,ust,z0_drydep(ix,jy))
      else
        ra=raerod(L,ust,z0(j))
      endif
      raquer=raquer+ra*slanduse(j)

  ! 3. Calculate surface resistance for gases
  !******************************************

      call getrc(nspec,lseason,j,tc,gr,rh,rr,rc)

  ! 4. Calculate deposition velocities for gases and ...
  ! 5. ... sum deposition velocities for all landuse classes
  !*********************************************************

      do i=1,nspec
        if (reldiff(i).gt.0.) then
          if ((ra+rb(i)+rc(i)).gt.0.) then
            vd=1./(ra+rb(i)+rc(i))
          else
            vd=9.999
          endif
          vdepo(i)=vdepo(i)+vd*slanduse(j)
        endif
      end do
    endif
  end do


  ! 6. Calculate deposition velocities for particles
  !*************************************************

  call partdep(nspec,density,fract,schmi,vset,raquer,ust,nyl, &
    rhoa,vdepo)

  !if (debug_mode) then
  !  print*,'getvdep:188: vdepo=', vdepo
    !stop
  !endif

  ! 7. If no detailed parameterization available, take constant deposition
  !    velocity if that is available
  !***********************************************************************

  do i=1,nspec
    if ((reldiff(i).lt.0.).and.(density(i).lt.0.).and. &
         (dryvel(i).gt.0.)) then
      vdepo(i)=dryvel(i)
    endif
  end do
end subroutine getvdep

subroutine getvdep_nest(n,ix,jy,ust,temp,pa, &
       L,gr,rh,rr,snow,vdepo,lnest)
  !                   i i  i   i   i   i  i i  i  i    i o i
  !*****************************************************************************
  !                                                                            *
  !  This routine calculates the dry deposition velocities.                    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 December 1996                                                       *
  !     Sabine Eckhardt, Jan 07                                                *
  !     if the latitude is negative: add half a year to the julian day         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! gr [W/m2]         global radiation                                         *
  ! L [m]             Obukhov length                                           *
  ! nyl               kinematic viscosity                                      *
  ! pa [Pa]           surface air pressure                                     *
  ! ra [s/m]          aerodynamic resistance                                   *
  ! raquer [s/m]      average aerodynamic resistance                           *
  ! rh [0-1]          relative humidity                                        *
  ! rhoa              density of the air                                       *
  ! rr [mm/h]         precipitation rate                                       *
  ! temp [K]          2m temperature                                           *
  ! tc [C]            2m temperature                                           *
  ! ust [m/s]         friction velocity                                        *
  ! snow [m of water equivalent] snow depth                                    *
  ! xlanduse          fractions of numclasS landuses for each model grid point *
  !                                                                            *
  !*****************************************************************************
  use date_mod

  implicit none

  integer :: yyyymmdd,hhmmss,yyyy,mmdd,n,lseason,i,j,ix,jy,lnest
  real :: vdepo(maxspec),vd,rb(maxspec),rc(maxspec),raquer,ylat
  real :: ra,ust,temp,tc,pa,L,gr,rh,rr,myl,nyl,rhoa,diffh2o,snow
  real :: slanduse(numclass)
  real,parameter :: eps=1.e-5
  real(kind=dp) :: jul

  ! Calculate month and determine the seasonal category
  !****************************************************

  jul=bdate+real(wftime(n),kind=dp)/86400._dp

  ylat=jy*dy+ylat0
  if (ylat.lt.0) then
      jul=jul+365.*0.5
  endif


  call caldate(jul,yyyymmdd,hhmmss)
  yyyy=yyyymmdd/10000
  mmdd=yyyymmdd-10000*yyyy

  if ((ylat.gt.-20).and.(ylat.lt.20)) then
     mmdd=600 ! summer
  endif

  if ((mmdd.ge.1201).or.(mmdd.le.301)) then
    lseason=4
  else if ((mmdd.ge.1101).or.(mmdd.le.331)) then
    lseason=3
  else if ((mmdd.ge.401).and.(mmdd.le.515)) then
    lseason=5
  else if ((mmdd.ge.516).and.(mmdd.le.915)) then
    lseason=1
  else
    lseason=2
  endif

  ! Calculate diffusivity of water vapor
  !************************************
  diffh2o=2.11e-5*(temp/273.15)**1.94*(101325/pa)

  ! Conversion of temperature from K to C
  !**************************************

  tc=temp-273.15

  ! Calculate dynamic viscosity
  !****************************

  if (tc.lt.0) then
    myl=(1.718+0.0049*tc-1.2e-05*tc**2)*1.e-05
  else
    myl=(1.718+0.0049*tc)*1.e-05
  endif

  ! Calculate kinematic viscosity
  !******************************

  rhoa=pa/(287.*temp)
  nyl=myl/rhoa


  ! 0. Set all deposition velocities zero
  !**************************************

  do i=1,nspec
    vdepo(i)=0.
  end do


  ! 1. Compute surface layer resistances rb
  !****************************************

  call getrb(nspec,ust,nyl,diffh2o,reldiff,rb)

  ! change for snow
  do j=1,numclass
    if (snow.gt.0.001) then ! 10 mm
       if (j.eq.12) then
         slanduse(j)=1.
       else
         slanduse(j)=0.
       endif
    else
       slanduse(j)=xlandusen(ix,jy,j,lnest)
    endif
  end do

  raquer=0.
  do j=1,numclass            ! loop over all landuse classes

    if (slanduse(j).gt.eps)  then

  ! 2. Calculate aerodynamic resistance ra
  !***************************************
      if (DRYDEP.and.(j.eq.7)) then
        ra=raerod(L,ust,z0_drydepn(ix,jy,lnest))
      else
        ra=raerod(L,ust,z0(j))
      endif
      raquer=raquer+ra*slanduse(j)

  ! 3. Calculate surface resistance for gases
  !******************************************

      call getrc(nspec,lseason,j,tc,gr,rh,rr,rc)

  ! 4. Calculate deposition velocities for gases and ...
  ! 5. ... sum deposition velocities for all landuse classes
  !*********************************************************

      do i=1,nspec
        if (reldiff(i).gt.0.) then
          if ((ra+rb(i)+rc(i)).gt.0.) then
            vd=1./(ra+rb(i)+rc(i))
  ! XXXXXXXXXXXXXXXXXXXXXXXXXX TEST
  !         vd=1./rc(i)
  ! XXXXXXXXXXXXXXXXXXXXXXXXXX TEST
          else
            vd=9.999
          endif
          vdepo(i)=vdepo(i)+vd*slanduse(j)
        endif
      end do
    endif
  end do


  ! 6. Calculate deposition velocities for particles
  !*************************************************

  call partdep(nspec,density,fract,schmi,vset,raquer,ust,nyl, &
    rhoa,vdepo)

  ! 7. If no detailed parameterization available, take constant deposition
  !    velocity if that is available
  !***********************************************************************

  do i=1,nspec
    if ((reldiff(i).lt.0.).and.(density(i).lt.0.).and. &
         (dryvel(i).gt.0.)) then
      vdepo(i)=dryvel(i)
    endif
  end do
end subroutine getvdep_nest

subroutine partdep(nc,density,fract,schmi,vset,ra,ustar,nyl,rhoa,vdep_tmp)
  !                   i     i      i     i    i   i    i    i  i, i, i/o
  !*****************************************************************************
  !                                                                            *
  !      Calculation of the dry deposition velocities of particles.            *
  !      This routine is based on Stokes' law for considering settling and     *
  !      assumes constant dynamic viscosity of the air.                        *
  !                                                                            *
  !     AUTHOR: Andreas Stohl, 12 November 1993                                *
  !                            Update: 20 December 1996                        *
  !                                                                            *
  !     Literature:                                                            *
  !     [1]  Hicks/Baldocchi/Meyers/Hosker/Matt (1987), A Preliminary          *
  !             Multiple Resistance Routine for Deriving Dry Deposition        *
  !             Velocities from Measured Quantities.                           *
  !             Water, Air and Soil Pollution 36 (1987), pp.311-330.           *
  !     [2]  Slinn (1982), Predictions for Particle Deposition to              *
  !             Vegetative Canopies. Atm.Env.16-7 (1982), pp.1785-1794.        *
  !     [3]  Slinn/Slinn (1980),  Predictions for Particle Deposition on       *
  !             Natural Waters. Atm.Env.14 (1980), pp.1013-1016.               *
  !     [4]  Scire/Yamartino/Carmichael/Chang (1989),                          *
  !             CALGRID: A Mesoscale Photochemical Grid Model.                 *
  !             Vol II: User's Guide. (Report No.A049-1, June, 1989)           *
  !     [5]  Langer M. (1992): Ein einfaches Modell zur Abschaetzung der       *
  !             Depositionsgeschwindigkeit von Teilchen und Gasen.             *
  !             Internal report.                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! alpha                help variable                                         *
  ! fract(nc,ni)         mass fraction of each diameter interval               *
  ! lpdep(nc)            1 for particle deposition, 0 else                     *
  ! nc                   actual number of chemical components                  *
  ! ni                   number of diameter intervals, for which vdepj is calc.*
  ! rdp [s/m]            deposition layer resistance                           *
  ! ra [s/m]             aerodynamical resistance                              *
  ! schmi(nc,ni)         Schmidt number**2/3 of each diameter interval         *
  ! stokes               Stokes number                                         *
  ! ustar [m/s]          friction velocity                                     *
  ! vdep_tmp(nc) [m/s]       deposition velocities of all components               *
  ! vdepj [m/s]          help, deposition velocity of 1 interval               *
  ! vset(nc,ni)          gravitational settling velocity of each interval      *
  !                                                                            *
  ! Constants:                                                                 *
  ! nc                   number of chemical species                            *
  ! ni                   number of diameter intervals, for which deposition    *
  !                      is calculated                                         *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real, intent(in) ::       &
    nyl,                    & ! kinematic viscosity
    rhoa,                   & ! air density
    ustar,                  & ! friction velocity
    ra,                     & ! aerodynamical resistance
    vset(maxspec,maxndia),  & ! gravitational settling velocity of each interval
    density(maxspec),       & ! density of the particle
    fract(maxspec,maxndia)    ! mass fraction of each diameter interval
  real, intent(inout) ::    &
    vdep_tmp(maxspec)
  real :: schmi(maxspec,maxndia)
  real :: stokes,vdepj,rdp,alpha
  real :: & ! Variables related to shape
    dfdr, alpha1, beta1, ks, kn, c_d, &
    settling, settling_old, reynolds, kn1

  real,parameter :: eps=1.e-5
  integer :: ic,j,nc,i


  do ic=1,nc                  ! loop over all species
    if (density(ic).gt.0.) then
      do j=1,ndia(ic)         ! loop over all diameter intervals
        if (ustar.gt.eps) then          
          if (ishape(ic).eq.0) then

            reynolds=dquer(ic)/1.e6*vset(ic,j)/nyl
            settling_old=-1.0*vset(ic,j)

            do i=1,20

              if (reynolds.le.0.02) then
                c_d=(24.0/reynolds)

              else ! Clif and Gauvin scheme is used
                c_d=(24.0/reynolds)*(1+0.15*(reynolds**0.687))+ &
                  0.42/(1.0+42500.0/(reynolds**1.16))
              endif


            ! Settling velocity of a particle is defined by the Newton's impact law:
              settling=-1.* &
                      sqrt(4.*ga*dquer(ic)/1.e6*density(ic)*cunningham(ic)/ &
                      (3.*c_d*rhoa))

              if (abs((settling-settling_old)/settling).lt.0.01) exit     
              reynolds=dquer(ic)/1.e6*abs(settling)/nyl
              settling_old=settling
            end do
      
            ! Stokes number for each diameter interval
            !*****************************************
            ! Use this stokes number for different shapes
            stokes=abs(settling)/ga*ustar*ustar/nyl
            alpha=-3./stokes

            ! Deposition layer resistance
            !****************************

            if (alpha.le.log10(eps)) then
              rdp=1./(schmi(ic,j)*ustar)
            else
              rdp=1./((schmi(ic,j)+10.**alpha)*ustar)
            endif

            vdepj=abs(settling)+1./(ra+rdp+ra*rdp*abs(settling))

          else ! Daria Tatsii: Drag coefficient scheme by Bagheri & Bonadonna 2016
               ! Settling velocities of other shapes

            reynolds=dquer(ic)/1.e6*vset(ic,j)/nyl
            settling_old=-1.0*vset(ic,j)

            ! Orientation of particles
            !*************************
            if (orient(ic).eq.0) then
              ! Horizontal orientation
              ks=ks2(ic)  ! B&B Figure 12 k_(s,max)
              kn=kn2(ic)
            else if (orient(ic).eq.1) then 
              ! Random orientation
              dfdr=density(ic)/rhoa

              alpha1=0.45+10.0/(exp(2.5*log10(dfdr))+30.0)
              beta1=1.-37.0/(exp(3.0*log10(dfdr))+100.0)
              ks=ks1(ic)
              kn=10.**(alpha1*(-log10(Fn(ic)))**beta1)
            else
              ! The average of random and horizontal orientation
              dfdr=density(ic)/rhoa

              alpha1=0.45+10.0/(exp(2.5*log10(dfdr))+30.0)
              beta1=1.-37.0/(exp(3.0*log10(dfdr))+100.0)
              kn1=10.**(alpha1*(-log10(Fn(ic)))**beta1)
              ks=(ks1(ic)+ks2(ic))*0.5
              kn=(kn1+kn2(ic))*0.5
            endif

            do i=1,20
              c_d=(24.*ks/reynolds)*(1.+0.125*((reynolds*kn/ks)**(2./3.)))+ &
                (0.46*kn/(1.+5330./(reynolds*kn/ks)))

              ! Settling velocity of a particle is defined by the Newton's impact law:
              settling=-1.* &
                      sqrt(4.*ga*dquer(ic)/1.e6*density(ic)*cunningham(ic)/ &
                      (3.*c_d*rhoa))

              if (abs((settling-settling_old)/settling).lt.0.01) exit

              reynolds=dquer(ic)/1.e6*abs(settling)/nyl
              settling_old=settling
            end do
            ! We assume aerodynamic resistance ra and quasi-laminar sub-layer resistance rdp
            ! Stokes number for each diameter interval
            !*****************************************
            ! Use this stokes number for different shapes
            stokes=abs(settling)/ga*ustar*ustar/nyl
            alpha=-3./stokes

            ! Deposition layer resistance
            !****************************

            if (alpha.le.log10(eps)) then
              rdp=1./(schmi(ic,j)*ustar)
            else
              rdp=1./((schmi(ic,j)+10.**alpha)*ustar)
            endif
            
            
            vdepj=abs(settling)+1./(ra+rdp+ra*rdp*abs(settling))

          endif

        else
          vdepj=vset(ic,j)
        endif

        ! deposition velocities of each interval are weighted with mass fraction
        !***********************************************************************

        vdep_tmp(ic)=vdep_tmp(ic)+vdepj*fract(ic,j)
          
      end do
    endif
  end do

end subroutine partdep

subroutine getrb(nc,ustar,nyl,diffh2o,reldiff,rb)
  !                 i    i    i     i       i    o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of the quasilaminar sublayer resistance to dry deposition.    *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 20 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! rb(ncmax)       sublayer resistance                                        *
  ! schmidt         Schmidt number                                             *
  ! ustar [m/s]     friction velocity                                          *
  ! diffh20 [m2/s]  diffusivity of water vapor in air                          *
  ! reldiff         diffusivity relative to H2O                                *
  !                                                                            *
  ! Constants:                                                                 *
  ! karman          von Karman constant                                        *
  ! pr              Prandtl number                                             *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real :: ustar,diffh2o,rb(maxspec),schmidt,nyl
  real :: reldiff(maxspec)
  integer :: ic,nc
  real,parameter :: pr=0.72

  do ic=1,nc
    if (reldiff(ic).gt.0.) then
      schmidt=nyl/diffh2o*reldiff(ic)
      rb(ic)=2.0*(schmidt/pr)**0.67/(karman*ustar)
    endif
  end do
end subroutine getrb

subroutine getrc(nc,i,j,t,gr,rh,rr,rc)
  !                 i  i i i i  i  i  o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of the surface resistance according to the procedure given    *
  !  in:                                                                       *
  !  Wesely (1989): Parameterization of surface resistances to gaseous         *
  !  dry deposition in regional-scale numerical models.                        *
  !  Atmos. Environ. 23, 1293-1304.                                            *
  !                                                                            *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 19 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! reldiff(maxspec)  diffusivity of H2O/diffusivity of component i            *
  ! gr [W/m2]       global radiation                                           *
  ! i               index of seasonal category                                 *
  ! j               index of landuse class                                     *
  ! ldep(maxspec)          1, if deposition shall be calculated for species i  *
  ! nc                   actual number of chemical components                  *
  ! rcl(maxspec,5,8) [s/m] Lower canopy resistance                             *
  ! rgs(maxspec,5,8) [s/m] Ground resistance                                   *
  ! rlu(maxspec,5,8) [s/m] Leaf cuticular resistance                           *
  ! rm(maxspec) [s/m]      Mesophyll resistance                                *
  ! t [C]           temperature                                                *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,j,ic,nc
  real :: gr,rh,rr,t,rs,rsm,corr,rluc,rclc,rgsc,rdc,rluo
  real :: rc(maxspec)


  ! Compute stomatal resistance
  !****************************
  ! Sabine Eckhardt, Dec 06: use 1E25 instead of 99999. for infinite res.

  if ((t.gt.0.).and.(t.lt.40.)) then
    rs=ri(i,j)*(1.+(200./(gr+0.1))**2)*(400./(t*(40.-t)))
  else
    rs=1.E25
  !  rs=99999.
  endif


  ! Correct stomatal resistance for effect of dew and rain
  !*******************************************************

  if ((rh.gt.0.9).or.(rr.gt.0.)) rs=rs*3.

  ! Compute the lower canopy resistance
  !************************************

  rdc=100.*(1.+1000./(gr+10.))


  corr=1000.*exp(-1.*t-4.)
  do ic=1,nc
    if (reldiff(ic).gt.0.) then

  ! Compute combined stomatal and mesophyll resistance
  !***************************************************

      rsm=rs*reldiff(ic)+rm(ic)

  ! Correct leaf cuticular, lower canopy and ground resistance
  !***********************************************************

      rluc=rlu(ic,i,j)+corr
      rclc=rcl(ic,i,j)+corr
      rgsc=rgs(ic,i,j)+corr

  ! Correct leaf cuticular resistance for effect of dew and rain
  !*************************************************************

      if (rr.gt.0.) then
        rluo=1./(1./1000.+1./(3.*rluc))
        rluc=1./(1./(3.*rluc)+1.e-7*henry(ic)+f0(ic)/rluo)
      else if (rh.gt.0.9) then
        rluo=1./(1./3000.+1./(3.*rluc))
        rluc=1./(1./(3.*rluc)+1.e-7*henry(ic)+f0(ic)/rluo)
      endif

  ! Combine resistances to give total resistance
  !*********************************************

      rc(ic)=1./(1./rsm+1./rluc+1./(rdc+rclc)+1./(rac(i,j)+rgsc))
  ! Sabine Eckhardt, Dec 06: avoid possible excessively high vdep
      if (rc(ic).lt.10.) rc(ic)=10.
    endif
  end do
end subroutine getrc

end module drydepo_mod
