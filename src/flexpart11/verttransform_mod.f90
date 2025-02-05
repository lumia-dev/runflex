! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  ! This module contains all subroutines computing the vertical coordinate     *
  ! transformation of the meteorological input data (L. Bakels 2021)           *
  !                                                                            *
  !*****************************************************************************

module verttransform_mod
  use par_mod
  use com_mod
  use qvsat_mod
  use cmapf_mod, only: cc2gll
  use windfields_mod

  implicit none

contains

subroutine verttransform_ecmwf(n,uuh,vvh,wwh,pvh)
  !                              i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !*****************************************************************************
  !  CHANGES                                                                   *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !   Variables tth and qvh (on eta coordinates) from common block             *
  !                                                                            *
  ! Sabine Eckhardt, March 2007                                                *
  ! added the variable cloud for use with scavenging - descr. in com_mod       *
  !                                                                            *
  ! Unified ECMWF and GFS builds                                               *
  ! Marian Harustak, 12.5.2017                                                 *
  !     - Renamed from verttransform to verttransform_ecmwf                    *
  !                                                                            *
  ! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
  !                                                                            *
  ! Lucie Bakels, 2022                                                         *
  !    - Separated the code into subroutines                                   *
  !    - In case of wind_coord_type='ETA': keep ECMWF vertical winds in eta    *
  !      coordinates                                                           *
  !    - OpenMP parallelisation                                                *
  !                                                                            *
  !  Petra Seibert, Anne Philipp, 2019-05-02: implement wetdepo quickfix       *
  !  Petra Seibert, Anne Tipka, 2020-11-19: reimplement in latest version      *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! Note PS, AT 2021-01-29: all these fields are 0:nxmax-1,0:nymax-1 !!        *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! icloudbot(0:nxmax,0:nymax,numwfmem) cloud bottom field for wet deposition  * 
  ! icloudtop(0:nxmax,0:nymax,numwfmem) cloud thickness for wet deposition    *
  ! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
  ! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
  ! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
  !                                          [deltaeta/s]                      *
  ! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
  ! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
  ! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in) :: n
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax) :: wwh

  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: rhoh
  real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: pinmconv
  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: prsh

  logical :: init = .true.

  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
  ! the vertical resolution in the z system is doubled. As reference point,*
  ! the lower left corner of the grid is used.                             *
  ! Unlike in the eta system, no difference between heights for u,v and    *
  ! heights for w exists.                                                  *
  !*************************************************************************


  !eso measure CPU time
  !  call mpif_mtime('verttransform',0)

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************
    call verttransform_init(n)

  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.
  endif


  ! Compute heights of eta levels and their respective pressure and density fields
  !*******************************************************************************
  call verttransform_ecmwf_heights(nxmin1,nymin1,tt2(0:nxmin1,0:nymin1,1,n), &
    td2(0:nxmin1,0:nymin1,1,n),ps(0:nxmin1,0:nymin1,1,n),qvh(0:nxmin1,0:nymin1,:,n), &
    tth(0:nxmin1,0:nymin1,:,n),prsh(0:nxmin1,0:nymin1,:), &
    rhoh(0:nxmin1,0:nymin1,:),pinmconv(0:nxmin1,0:nymin1,:), &
    etauvheight(0:nxmin1,0:nymin1,:,n),etawheight(0:nxmin1,0:nymin1,:,n))

  ! Transform the wind fields to the internal coordinate system and save the native ETA 
  ! fields when case wind_coord_type==ETA
  !*************************************************************
  call verttransform_ecmwf_windfields(n,nxmin1,nymin1, &
    uuh(0:nxmin1,0:nymin1,1:nuvzmax), &
    vvh(0:nxmin1,0:nymin1,1:nuvzmax),wwh(0:nxmin1,0:nymin1,1:nwzmax), &
    pvh(0:nxmin1,0:nymin1,1:nuvzmax),rhoh(0:nxmin1,0:nymin1,1:nuvzmax), &
    prsh(0:nxmin1,0:nymin1,1:nuvzmax),pinmconv(0:nxmin1,0:nymin1,1:nzmax))

  ! If north or south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************
  call verttransform_ecmwf_polar(n)

  ! Create cloud fields
  !*********************
#ifdef ETA
    call verttransform_ecmwf_cloud(lcw,lcwsum,nxmin1,nymin1, &
      ctwc(0:nxmin1,0:nymin1,n), &
      clwc(0:nxmin1,0:nymin1,:,n), &
      ciwc(0:nxmin1,0:nymin1,:,n), &
      icloudbot(0:nxmin1,0:nymin1,n), &
      icloudtop(0:nxmin1,0:nymin1,n), &
      lsprec(0:nxmin1,0:nymin1,1,:,n), &
      convprec(0:nxmin1,0:nymin1,1,:,n), &
      rhoeta(0:nxmin1,0:nymin1,:,n), &
      tteta(0:nxmin1,0:nymin1,:,n), &
      qv(0:nxmin1,0:nymin1,:,n), &
      etauvheight(0:nxmin1,0:nymin1,:,n), &
      etawheight(0:nxmin1,0:nymin1,:,n))
#else
    call verttransform_ecmwf_cloud(lcw,lcwsum,nxmin1,nymin1, &
      ctwc(0:nxmin1,0:nymin1,n),clwc(0:nxmin1,0:nymin1,:,n), &
      ciwc(0:nxmin1,0:nymin1,:,n), &
      icloudbot(0:nxmin1,0:nymin1,n), icloudtop(0:nxmin1,0:nymin1,n), &
      lsprec(0:nxmin1,0:nymin1,1,:,n),convprec(0:nxmin1,0:nymin1,1,:,n), &
      rho(0:nxmin1,0:nymin1,:,n), tt(0:nxmin1,0:nymin1,:,n), &
      qv(0:nxmin1,0:nymin1,:,n), etauvheight(0:nxmin1,0:nymin1,:,n), &
      etawheight(0:nxmin1,0:nymin1,:,n)) 
#endif
end subroutine verttransform_ecmwf

subroutine verttransform_nest(n,uuhn,vvhn,wwhn,pvhn)
  !                            i   i    i    i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !     It is similar to verttransform, but makes the transformations for      *
  !     the nested grids.                                                      *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !*****************************************************************************
  !  CHANGES                                                                   *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !  Bernd C. Krueger, Feb. 2001:                                              *
  !   Variables tthn and qvhn (on eta coordinates) from common block           *
  !                                                                            *
  ! Sabine Eckhardt, March 2007:                                               *
  ! added the variable cloud for use with scavenging - descr. in com_mod       *
  ! PS/AT 2018/-21: variable "cloud" is replaced by quickfix, see below        * 
  !                                                                            *
  ! ESO, 2016                                                                  *
  ! -note that divide-by-zero occurs when nxmaxn,nymaxn etc. are larger than   *
  !  the actual field dimensions                                               *
  !                                                                            *
  ! Don Morton, 2017-05-30:                                                    *
  !   modification of a bug in ew. Don Morton (CTBTO project)                  *
  !                                                                            *
  !  undocumented modifications by NILU for v10                                *
  !                                                                            *
  !  Petra Seibert, 2018-06-13:                                                *
  !   - put back SAVE attribute for INIT, just to be safe                      *
  !   - minor changes, most of them just cosmetics                             *
  !   for details see changelog.txt in branch unive                            *
  !                                                                            *
  !  Petra Seibert, Anne Philipp, 2019-05-02: implement wetdepo quickfix       *
  !  Petra Seibert, Anne Tipka, 2020-11-19: reimplement in latest version      *
  !                                                                            *
  ! ****************************************************************************
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! Note PS, AT 2021-01-29: all these fields are 0:nxmaxn-1,0:nymaxn-1 !!      *
  ! nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction    *
  ! icloudbot                       cloud bottom field for wet deposition      *
  ! icloudtop                      cloud thickness for wet deposition         *
  ! uun                             wind components in x-direction [m/s]       *
  ! vvn                             wind components in y-direction [m/s]       *
  ! wwn                             wind components in z-direction [deltaeta/s]*
  ! ttn                             temperature [K]                            *
  ! pvn                             potential vorticity (pvu)                  *
  ! psn                             surface pressure [Pa]                      *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numbnests) :: &
    uuhn,vvhn,pvhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nwzmax,numbnests) :: wwhn

  real,dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: rhohn,prshn
  real,dimension(0:nxmaxn-1,0:nymaxn-1,nzmax) :: pinmconv

  integer :: nxm1, nym1
  integer :: n,l


  ! Loop over all nests
  !********************

  do l=1,numbnests
    nxm1=nxn(l)-1
    nym1=nyn(l)-1
    if (nxm1.lt.1 .or. nym1.lt.1 ) cycle
    call verttransform_ecmwf_heights(nxm1,nym1, &
      tt2n(0:nxm1,0:nym1,1,n,l),td2n(0:nxm1,0:nym1,1,n,l),psn(0:nxm1,0:nym1,1,n,l), &
      qvhn(0:nxm1,0:nym1,:,n,l),tthn(0:nxm1,0:nym1,:,n,l),prshn(0:nxm1,0:nym1,:), &
      rhohn(0:nxm1,0:nym1,:),pinmconv(0:nxm1,0:nym1,:), &
      etauvheightn(0:nxm1,0:nym1,:,n,l),etawheightn(0:nxm1,0:nym1,:,n,l))

    call verttransform_ecmwf_windfields_nest(l,n,uuhn,vvhn,wwhn,pvhn, &
                                                        rhohn,prshn,pinmconv)

    ! Create cloud fields
    !*********************

#ifdef ETA
    call verttransform_ecmwf_cloud(lcw_nest(l),lcwsum_nest(l),nxm1,nym1,&
      ctwcn(0:nxm1,0:nym1,n,l), clwcn(0:nxm1,0:nym1,:,n,l), ciwcn(0:nxm1,0:nym1,:,n,l), &
      icloudbotn(0:nxm1,0:nym1,n,l), icloudtopn(0:nxm1,0:nym1,n,l), &
      lsprecn(0:nxm1,0:nym1,1,:,n,l),convprecn(0:nxm1,0:nym1,1,:,n,l), &
      rhoetan(0:nxm1,0:nym1,:,n,l),ttetan(0:nxm1,0:nym1,:,n,l), &
      qvn(0:nxm1,0:nym1,:,n,l), etauvheightn(0:nxm1,0:nym1,:,n,l), &
      etawheightn(0:nxm1,0:nym1,:,n,l))
#else
    call verttransform_ecmwf_cloud(lcw_nest(l),lcwsum_nest(l),nxm1,nym1,&
      ctwcn(0:nxm1,0:nym1,n,l), clwcn(0:nxm1,0:nym1,:,n,l), ciwcn(0:nxm1,0:nym1,:,n,l), &
      icloudbotn(0:nxm1,0:nym1,n,l), icloudtopn(0:nxm1,0:nym1,n,l), &
      lsprecn(0:nxm1,0:nym1,1,:,n,l),convprecn(0:nxm1,0:nym1,1,:,n,l), &
      rhon(0:nxm1,0:nym1,:,n,l),ttn(0:nxm1,0:nym1,:,n,l), &
      qvn(0:nxm1,0:nym1,:,n,l), etauvheightn(0:nxm1,0:nym1,:,n,l), &
      etawheightn(0:nxm1,0:nym1,:,n,l))
#endif
      
  end do ! end loop over nests
end subroutine verttransform_nest

subroutine verttransform_init(n)
  implicit none

  integer, intent(in) :: n
  real ::  tvold,pold,pint,tv
  integer :: ix,jy,kz,ixm,jym
  real,parameter :: const=r_air/ga

  if ((ipin.eq.1).or.(ipin.eq.4)) then
    call read_heightlevels(height,nmixz)
    return
  endif

  jym=nymin1
  ixm=0
  loop1: do jy=0,nymin1
    do ix=0,nxmin1
      if (ps(ix,jy,1,n).gt.100000.) then
        ixm=ix
        jym=jy
        exit loop1
      endif
    end do
  end do loop1

  tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
       ps(ixm,jym,1,n))
  pold=ps(ixm,jym,1,n)
  height(1)=0.

  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
    tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

    if (abs(tv-tvold).gt.0.2) then
      height(kz)= height(kz-1)+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
    else
      height(kz)=height(kz-1)+const*log(pold/pint)*tv
    endif

    tvold=tv
    pold=pint
  end do

  ! Determine highest levels that can be within PBL
  !************************************************

  do kz=1,nz
    if (height(kz).gt.hmixmax) then
      nmixz=kz
      exit
    endif
  end do

  if (loutrestart.ne.-1) then
    call output_heightlevels(height,nmixz)
  endif
end subroutine verttransform_init

subroutine output_heightlevels(height_tmp,nmixz_tmp)
  implicit none

  real,intent(in) :: height_tmp(nzmax)
  integer,intent(in) :: nmixz_tmp
  integer :: kz
  character(len=256) :: heightlevels_filename

  heightlevels_filename = path(2)(1:length(2))//'heightlevels.bin'

  write(*,*) 'Writing Initialised heightlevels to file:', trim(heightlevels_filename)
  
  open(unitheightlevels,file=trim(heightlevels_filename),form='unformatted')

  write(unitheightlevels) nmixz_tmp

  do kz=1,nz
    write(unitheightlevels) height_tmp(kz)
  end do
  close(unitheightlevels)
end subroutine output_heightlevels

subroutine read_heightlevels(height_tmp,nmixz_tmp)
  implicit none

  real,intent(out) :: height_tmp(nzmax)
  integer,intent(out) :: nmixz_tmp
  integer :: kz,ios
  character(len=256) :: heightlevels_filename

  heightlevels_filename = path(2)(1:length(2))//'heightlevels.bin'

  write(*,*) 'Reading heightlevels from file:', trim(heightlevels_filename)
  
  open(unitheightlevels,file=trim(heightlevels_filename),form='unformatted',err=9988)

  read(unitheightlevels,iostat=ios) nmixz_tmp

  do kz=1,nz
    read(unitheightlevels) height_tmp(kz)
  end do
  close(unitheightlevels)

  return

9988   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE             #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'heightlevels.bin'//'    #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS        #### '
  write(*,*) ' #### NAME DOES NOT EXISTS, REMOVE call read_heightlevels #### '
  write(*,*) ' #### FROM VERTTRANSFORM_MOD.                 #### '
end subroutine read_heightlevels

subroutine verttransform_ecmwf_windfields(n,nxlim,nylim,uuh,vvh,wwh,pvh,rhoh,prsh,pinmconv)
  implicit none

  integer,intent(in) :: n,nxlim,nylim
  real,intent(in),dimension(0:nxlim,0:nylim,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxlim,0:nylim,nwzmax) :: wwh
  real,intent(in),dimension(0:nxlim,0:nylim,nuvzmax) :: rhoh
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: pinmconv
  ! RLT added pressure
  real,intent(in),dimension(0:nxlim,0:nylim,nuvzmax) :: prsh

  !real,dimension(0:nxmax-1,0:nymax-1) ::  dpdeta

  real,dimension(0:nylim) :: cosf

  integer,dimension(0:nxlim,0:nylim,nzmax) :: idx,idxw

  integer :: ix,jy,kz,iz,ixp,jyp,ix1,jy1
  real :: dz1,dz2,dz,dpdeta
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2

  ! Copy fields for ETA coordinate interpolations
#ifdef ETA
!$OMP PARALLEL PRIVATE(ix,jy,kz)
!$OMP WORKSHARE
  uueta(0:nxlim,0:nylim,:,n) = uuh(0:nxlim,0:nylim,:)
  vveta(0:nxlim,0:nylim,:,n) = vvh(0:nxlim,0:nylim,:)
  tteta(0:nxlim,0:nylim,:,n) = tth(0:nxlim,0:nylim,:,n)
  qv(0:nxlim,0:nylim,:,n)    = qvh(0:nxlim,0:nylim,:,n)
  pveta(0:nxlim,0:nylim,:,n) = pvh(0:nxlim,0:nylim,:)
  rhoeta(0:nxlim,0:nylim,:,n) = rhoh(0:nxlim,0:nylim,:)
  prseta(0:nxlim,0:nylim,:,n) = prsh(0:nxlim,0:nylim,:)

  forall (ix=0:nxlim,jy=0:nylim,kz=2:nz-1)
    drhodzeta(ix,jy,kz,n)= &
         (rhoh(ix,jy,kz+1)-rhoh(ix,jy,kz-1))/ &
         (height(kz+1)-height(kz-1)) 
         ! Note that this is still in SI units and not in eta
  end forall
  drhodzeta(0:nxlim,0:nylim,1,n)=(rhoh(0:nxlim,0:nylim,2)-rhoh(0:nxlim,0:nylim,1))/(height(2)-height(1))
  drhodzeta(0:nxlim,0:nylim,nz,n)=drhodzeta(0:nxlim,0:nylim,nz-1,n)

  ! Convert w from Pa/s to eta/s, following FLEXTRA
  !************************************************
  ! z=1
  wweta(0:nxlim,0:nylim,1,n)=wwh(0:nxlim,0:nylim,1)/ &
    ((akm(2)-akm(1)+(bkm(2)-bkm(1))*ps(0:nxlim,0:nylim,1,n))/ &
        (wheight(2)-wheight(1)))
  ! z=nuvz-1
  wweta(0:nxlim,0:nylim,nuvz-1,n)=wwh(0:nxlim,0:nylim,nuvz-1)/ &
    ((akm(nuvz-1)-akm(nuvz-2)+(bkm(nuvz-1)-bkm(nuvz-2))*ps(0:nxlim,0:nylim,1,n))/ &
        (wheight(nuvz-1)-wheight(nuvz-2)))
  ! 1<z<nuvz-1
  forall (ix=0:nxlim,jy=0:nylim,kz=2:nuvz-2)
    wweta(ix,jy,kz,n)=wwh(ix,jy,kz)/ &
      ((akm(kz+1)-akm(kz-1)+(bkm(kz+1)-bkm(kz-1))*ps(ix,jy,1,n))/ &
        (wheight(kz+1)-wheight(kz-1)))
  end forall
!$OMP END WORKSHARE NOWAIT

  if (lcw) then
!$OMP DO
    do kz=1,nz
      clwc(0:nxlim,0:nylim,kz,n)=clwch(0:nxlim,0:nylim,kz,n)
      if (.not. lcwsum) ciwc(0:nxlim,0:nylim,kz,n)=ciwch(0:nxlim,0:nylim,kz,n)
    end do
!$OMP END DO
  endif
!$OMP END PARALLEL
#endif

!$OMP PARALLEL PRIVATE(jy,ix,iz,kz,dz1,dz2,dz,ix1,jy1,ixp,jyp,dzdx1,dzdx2,dzdx, &
!$OMP dzdy1,dzdy2,dzdy,dpdeta)
  ! Finding the index in eta levels (uv and w) that correspond to
  ! a certain height level in meters
  !**************************************************************
!$OMP WORKSHARE
  idx(0:nxlim,0:nylim,1:2)=2
  idxw(0:nxlim,0:nylim,1:2)=2
!$OMP END WORKSHARE

!$OMP DO
  do jy=0,nylim
    do iz=2,nz-1
      do ix=0,nxlim
        idx(ix,jy,iz)=idx(ix,jy,iz-1)
        idxw(ix,jy,iz)=idxw(ix,jy,iz-1)

        ! height in meters that corresponds to the eta w level
        innwz: do kz=idxw(ix,jy,iz),nuvz
          if ((idxw(ix,jy,iz).le.kz).and. &
            (height(iz).gt.etawheight(ix,jy,kz-1,n)).and. &
            (height(iz).le.etawheight(ix,jy,kz,n))) then
            idxw(ix,jy,iz)=kz
            exit innwz
          endif
        enddo innwz

        ! height in meters that corresponds to the eta uv level
        if(height(iz).gt.etauvheight(ix,jy,nuvz,n)) then
          cycle
        else
          innuvz: do kz=idx(ix,jy,iz),nuvz
            if ((idx(ix,jy,iz).le.kz).and. &
              (height(iz).gt.etauvheight(ix,jy,kz-1,n)).and. &
              (height(iz).le.etauvheight(ix,jy,kz,n))) then
              idx(ix,jy,iz)=kz
              exit innuvz
            endif
          enddo innuvz
        endif
      end do
    end do
  end do
!$OMP END DO NOWAIT

  ! Setting upper and lower levels
!$OMP WORKSHARE
  uu(0:nxlim,0:nylim,1,n)=uuh(0:nxlim,0:nylim,1)
  uu(0:nxlim,0:nylim,nz,n)=uuh(0:nxlim,0:nylim,nuvz)
  vv(0:nxlim,0:nylim,1,n)=vvh(0:nxlim,0:nylim,1)
  vv(0:nxlim,0:nylim,nz,n)=vvh(0:nxlim,0:nylim,nuvz)
  tt(0:nxlim,0:nylim,1,n)=tth(0:nxlim,0:nylim,1,n)
  tt(0:nxlim,0:nylim,nz,n)=tth(0:nxlim,0:nylim,nuvz,n)
  pv(0:nxlim,0:nylim,1,n)=pvh(0:nxlim,0:nylim,1)
  pv(0:nxlim,0:nylim,nz,n)=pvh(0:nxlim,0:nylim,nuvz)
#ifndef ETA
  qv(0:nxlim,0:nylim,1,n)=qvh(0:nxlim,0:nylim,1,n)
  qv(0:nxlim,0:nylim,nz,n)=qvh(0:nxlim,0:nylim,nuvz,n)
#endif
  rho(0:nxlim,0:nylim,1,n)=rhoh(0:nxlim,0:nylim,1)
  rho(0:nxlim,0:nylim,nz,n)=rhoh(0:nxlim,0:nylim,nuvz)
  ! RLT add pressure
  prs(0:nxlim,0:nylim,1,n)=prsh(0:nxlim,0:nylim,1)
  prs(0:nxlim,0:nylim,nz,n)=prsh(0:nxlim,0:nylim,nuvz)
  ! RLT
  ww(0:nxlim,0:nylim,1,n)=wwh(0:nxlim,0:nylim,1)*pinmconv(0:nxlim,0:nylim,1)
  ww(0:nxlim,0:nylim,nz,n)=wwh(0:nxlim,0:nylim,nwz)*pinmconv(0:nxlim,0:nylim,nz)
  forall (jy=0:nylim) 
    cosf(jy)=1./cos((real(jy)*dy+ylat0)*pi180) ! Needed in slope computations
  end forall
!$OMP END WORKSHARE NOWAIT

#ifndef ETA
  if (lcw) then !hg adding the cloud water 
!$OMP WORKSHARE
    clwc(0:nxlim,0:nylim,1,n)=clwch(0:nxlim,0:nylim,1,n)
    clwc(0:nxlim,0:nylim,nz,n)=clwch(0:nxlim,0:nylim,nuvz,n)
!$OMP END WORKSHARE NOWAIT
    if (.not. lcwsum) then
!$OMP WORKSHARE
      ciwc(0:nxlim,0:nylim,1,n)=ciwch(0:nxlim,0:nylim,1,n)
      ciwc(0:nxlim,0:nylim,nz,n)=ciwch(0:nxlim,0:nylim,nuvz,n) 
!$OMP END WORKSHARE NOWAIT
    endif
  endif
#endif

!$OMP BARRIER

!$OMP DO
  do iz=2,nz-1
    do jy=0,nylim
      do ix=0,nxlim
        ! Levels, where uv is given
        !*************************
        if (height(iz).gt.etauvheight(ix,jy,nuvz,n)) then
          uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
          vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
          tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
          pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
          rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
          prs(ix,jy,iz,n)=prs(ix,jy,nz,n)   ! RLT
#ifndef ETA
          qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
          !hg adding the cloud water
          if (lcw) then
            clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            if (.not.lcwsum) ciwc(ix,jy,iz,n)=ciwc(ix,jy,nz,n)
          end if
#endif
        else
          kz=idx(ix,jy,iz)
          dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
          dz2=etauvheight(ix,jy,kz,n)-height(iz)
          dz=dz1+dz2
          uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
          vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
          tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
               +tth(ix,jy,kz,n)*dz1)/dz
          pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
          rho(ix,jy,iz,n)=(rhoh(ix,jy,kz-1)*dz2+rhoh(ix,jy,kz)*dz1)/dz
          ! RLT add pressure
          prs(ix,jy,iz,n)=(prsh(ix,jy,kz-1)*dz2+prsh(ix,jy,kz)*dz1)/dz
#ifndef ETA
          qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2+qvh(ix,jy,kz,n)*dz1)/dz
          !hg adding the cloud water
          if  (lcw) then
            clwc(ix,jy,iz,n)= &
              (clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
            if (.not.lcwsum) ciwc(ix,jy,iz,n)= &
              (ciwch(ix,jy,kz-1,n)*dz2+ciwch(ix,jy,kz,n)*dz1)/dz
          end if
#endif
        endif
        ! Levels, where w is given
        !*************************
        kz=idxw(ix,jy,iz)

        dz1=height(iz)-etawheight(ix,jy,kz-1,n)
        dz2=etawheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(ix,jy,kz-1)*dz2 &
             +wwh(ix,jy,kz)*pinmconv(ix,jy,kz)*dz1)/dz

        !****************************************************************
        ! Compute slope of eta levels in windward direction and resulting
        ! vertical wind correction
        !****************************************************************
        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1
        if ((jy.eq.nylim).or.(jy.eq.0)) then
          cycle
        else if ((.not.xglobal).and.((ix.eq.nxlim).or.(ix.eq.0))) then
          cycle
        else if (ix.eq.nxlim) then
          ixp=0
        else if (ix.eq.0) then
          ix1=nxlim
        endif
        kz=idx(ix,jy,iz)
        dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
        dz2=etauvheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        
        dzdx1=(etauvheight(ixp,jy,kz-1,n)-etauvheight(ix1,jy,kz-1,n))*0.5
        dzdx2=(etauvheight(ixp,jy,kz,n)-etauvheight(ix1,jy,kz,n))*0.5

        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(etauvheight(ix,jyp,kz-1,n)-etauvheight(ix,jy1,kz-1,n))*0.5
        dzdy2=(etauvheight(ix,jyp,kz,n)-etauvheight(ix,jy1,kz,n))*0.5
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n) + dzdx*uu(ix,jy,iz,n)*dxconst*cosf(jy) &
                                      + dzdy*vv(ix,jy,iz,n)*dyconst
      end do
    end do
  end do
!$OMP END DO

!$OMP WORKSHARE
  ! Compute density gradients
  !**************************
  drhodz(0:nxlim,0:nylim,nz,n)=drhodz(0:nxlim,0:nylim,nz-1,n)
  drhodz(0:nxlim,0:nylim,1,n)=(rho(0:nxlim,0:nylim,2,n)-rho(0:nxlim,0:nylim,1,n))/ &
    (height(2)-height(1))
  forall (ix=0:nxlim,jy=0:nylim,iz=2:nz-1)
    drhodz(ix,jy,iz,n)=(rho(ix,jy,iz+1,n)-rho(ix,jy,iz-1,n))/ &
      (height(iz+1)-height(iz-1))
  end forall
!$OMP END WORKSHARE
!$OMP END PARALLEL

end subroutine verttransform_ecmwf_windfields

subroutine verttransform_ecmwf_polar(n)
  implicit none

  integer, intent(in) :: n

  integer :: ix,jy,iz
  real :: xlon,ylat,xlonr
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy

  if (nglobal) then
!$OMP PARALLEL PRIVATE(iz,jy,ix,xlon,ylat)
!$OMP DO
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
#ifdef ETA
          call cc2gll(northpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
               vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
               vvpoleta(ix,jy,iz,n))
#endif
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(iz,jy,ix,xlon,xlonr,ffpol,ddpol,uuaux,vvaux,uupolaux, &
!$OMP wdummy,vvpolaux)
!$OMP DO
    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+ &
           vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else if (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

      ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=90.
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)


      ! Fix: Set W at pole to the zonally averaged W of the next equator-
      ! ward parallel of latitude
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do

#ifdef ETA
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uueta(nx/2-1,nymin1,iz,n)**2+ &
           vveta(nx/2-1,nymin1,iz,n)**2)
      if (vveta(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uueta(nx/2-1,nymin1,iz,n)/ &
             vveta(nx/2-1,nymin1,iz,n))-xlonr
      else if (vveta(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uueta(nx/2-1,nymin1,iz,n)/ &
             vveta(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+wweta(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        wweta(ix,jy,iz,n)=wdummy
        uupoleta(ix,jy,iz,n)=uupolaux
        vvpoleta(ix,jy,iz,n)=vvpolaux
      end do

#endif
    end do
!$OMP END DO
!$OMP END PARALLEL


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    ! do iz=1,nz
    !   wdummy=0.
    !   jy=ny-2
    !   do ix=0,nxmin1
    !     wdummy=wdummy+ww(ix,jy,iz,n)
    !   end do
    !   wdummy=wdummy/real(nx)
    !   jy=nymin1
    !   do ix=0,nxmin1
    !     ww(ix,jy,iz,n)=wdummy
    !   end do
    ! end do

! #ifdef ETA
!     do iz=1,nz
!       wdummy=0.
!       jy=ny-2
!       do ix=0,nxmin1
!         wdummy=wdummy+wweta(ix,jy,iz,n)
!       end do
!       wdummy=wdummy/real(nx)
!       jy=nymin1
!       do ix=0,nxmin1
!         wweta(ix,jy,iz,n)=wdummy
!       end do
!     end do
! #endif

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
!$OMP PARALLEL PRIVATE(iz,jy,ix,xlon,ylat)
!$OMP DO
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
#ifdef ETA
          call cc2gll(southpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
               vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
               vvpoleta(ix,jy,iz,n))
#endif
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(iz,jy,ix,xlon,xlonr,ffpol,ddpol,uuaux,vvaux,uupolaux, &
!$OMP wdummy,vvpolaux)
!$OMP DO
    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+ &
           vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else if (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=-90.
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

    ! Fix: Set W at pole to the zonally averaged W of the next equator-
    ! ward parallel of latitude
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do

#ifdef ETA
  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uueta(nx/2-1,0,iz,n)**2+ &
           vveta(nx/2-1,0,iz,n)**2)
      if (vveta(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uueta(nx/2-1,0,iz,n)/ &
             vveta(nx/2-1,0,iz,n))+xlonr
      else if (vveta(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uueta(nx/2-1,0,iz,n)/ &
             vveta(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+wweta(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        wweta(ix,jy,iz,n)=wdummy
        uupoleta(ix,jy,iz,n)=uupolaux
        vvpoleta(ix,jy,iz,n)=vvpolaux
      end do
#endif
    end do
!$OMP END DO
!$OMP END PARALLEL
  endif
end subroutine verttransform_ecmwf_polar

subroutine verttransform_ecmwf_cloud(lcw_tmp,lcwsum_tmp,nxlim,nylim,&
  ctwc_tmp,clwc_tmp,ciwc_tmp,icloudbot_tmp,icloudtop_tmp,lsprec_tmp,convprec_tmp,rho_tmp, &
  tt_tmp,qv_tmp,uvzlev,wzlev)
  implicit none

  logical,intent(in) :: lcw_tmp,lcwsum_tmp
  integer, intent(in) :: nxlim,nylim
  real,intent(out) :: ctwc_tmp(0:nxlim,0:nylim)

  real,intent(inout) :: clwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: ciwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: lsprec_tmp(0:nxlim,0:nylim,numpf),convprec_tmp(0:nxlim,0:nylim,numpf)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: rho_tmp,tt_tmp,qv_tmp
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev,wzlev

  integer,intent(out) :: icloudbot_tmp(0:nxlim,0:nylim), icloudtop_tmp(0:nxlim,0:nylim)

  integer :: ix,jy

  ! converted parameters for eta coordinates:
  integer :: max_cloudthck_eta
  integer, dimension(2) :: conv_clrange_eta,highconvp_clrange_eta,lowconvp_clrange_eta
  
  ! AT, PS: for v11, we add back the quick fix to interpolate clouds in 
  !   interpol_rain developed by PS for v8 and extend it to using 
  !   cloud water fields (in apply_cloud_bounds)


! !$OMP PARALLEL PRIVATE(ix,jy,kz,k,lsp,convp,prec, &
! !$OMP max_cloudthck_eta,conv_clrange_eta,conv_clrange_eta,highconvp_clrange_eta, &
! $OMP highconvp_clrange_eta,lowconvp_clrange_eta,lowconvp_clrange_eta) REDUCTION(+:ctwc_tmp)
! !$OMP DO SCHEDULE(dynamic,max(1,nylim/50))
    
  do jy=0,nylim
    do ix=0,nxlim

#ifdef ETA
      call convert_cloud_params(ix,jy,nxlim,nylim,max_cloudthck_eta,conv_clrange_eta, &
        highconvp_clrange_eta,lowconvp_clrange_eta,uvzlev)
#endif

      icloudbot_tmp(ix,jy) = icmv 

      ! Find the bottom and top of clouds in grid cell ix, jy
      call identify_cloud(ix,jy,lcw_tmp,lcwsum_tmp,nxlim,nylim, &
        ctwc_tmp,clwc_tmp,ciwc_tmp,icloudbot_tmp,icloudtop_tmp,rho_tmp, &
        tt_tmp,qv_tmp,uvzlev,wzlev)

      ! Adjust clouds according to minimum thickness, height, lower level, etc.
      call apply_cloud_bounds(ix,jy,nxlim,nylim,lsprec_tmp,convprec_tmp,uvzlev, &
        icloudbot_tmp,icloudtop_tmp,max_cloudthck_eta,conv_clrange_eta, &
        highconvp_clrange_eta,lowconvp_clrange_eta)
    enddo ! ix loop
  enddo ! jy loop
end subroutine verttransform_ecmwf_cloud

subroutine convert_cloud_params(ix,jy,nxlim,nylim,max_cloudthck_eta,conv_clrange_eta, &
  highconvp_clrange_eta,lowconvp_clrange_eta,uvzlev)

  implicit none

  integer, intent(in) :: ix,jy,nxlim,nylim
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev

  ! converted parameters for eta coordinates:
  integer,intent(out) :: max_cloudthck_eta
  integer, dimension(2),intent(out) :: conv_clrange_eta,highconvp_clrange_eta, &
    lowconvp_clrange_eta
  integer :: i, kz


  ! Convert cloud parameters to eta coords.
  ! Reverse sign when using eta (eta:1-0, meter:0-max)
  max_cloudthck_eta=int(uvheight(nz)*eta_convert)
  conv_clrange_eta=int(uvheight(nz)*eta_convert)
  highconvp_clrange_eta=int(uvheight(nz)*eta_convert)
  lowconvp_clrange_eta=int(uvheight(nz)*eta_convert)
  do kz=1,nz
    if (uvzlev(ix,jy,kz).gt.max_cloudthck) then
      max_cloudthck_eta=int(uvheight(kz)*eta_convert)
      exit
    endif
  end do
  do i=1,2
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.conv_clrange(i)) then 
        conv_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.highconvp_clrange(i)) then 
        highconvp_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.lowconvp_clrange(i)) then 
        lowconvp_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
  end do
end subroutine convert_cloud_params

subroutine identify_cloud(ix,jy,lcw_tmp,lcwsum_tmp,nxlim,nylim, &
  ctwc_tmp,clwc_tmp,ciwc_tmp,icloudbot_tmp,icloudtop_tmp,rho_tmp, &
  tt_tmp,qv_tmp,uvzlev,wzlev)

  implicit none

  logical,intent(in) :: lcw_tmp,lcwsum_tmp
  integer, intent(in) :: ix,jy,nxlim,nylim
  real,intent(out) :: ctwc_tmp(0:nxlim,0:nylim)

  real,intent(inout) :: clwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: ciwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: rho_tmp,tt_tmp,qv_tmp
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev,wzlev

  integer,intent(out) :: icloudbot_tmp(0:nxlim,0:nylim), icloudtop_tmp(0:nxlim,0:nylim)

  integer :: kz
  real :: pressure,rh
  real :: clw

  !*******************************************************************************
  if (lcw_tmp) then ! identify clouds based on cloud water content
  !*******************************************************************************

    ctwc_tmp(ix,jy) = 0. ! initialise cloud total water content
    if (.not.lcwsum_tmp) clwc_tmp(ix,jy,:) = clwc_tmp(ix,jy,:) + ciwc_tmp(ix,jy,:)
    do kz = 1,nz-1 ! Changed order of loop to prevent ETA computation to be done unnecessarily

      ! vertically integrate cloud water and determine cloud bottom, top
      ! cloud water per cell in kg / m2
      ! calculate cloud water mass per area: kgCW/kgAIR * kgAIR/m3 * m = kgCW/m2
      
      ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3
#ifdef ETA
      clw = clwc_tmp(ix,jy,kz)*rho_tmp(ix,jy,kz)*(uvzlev(ix,jy,kz+1)-uvzlev(ix,jy,kz))
#else
      clw = clwc_tmp(ix,jy,kz)*rho_tmp(ix,jy,kz)*(height(kz+1)-height(kz)) 
#endif
      ! Add this layer to column cloud water [m3/m3]
      ctwc_tmp(ix,jy) = ctwc_tmp(ix,jy)+clw ! kg / m2 (or eta) in column

      if (clw .gt. 0.) then ! cloud layer - maybe use threshold?
#ifdef ETA
        if (icloudbot_tmp(ix,jy) .eq. icmv) & !cloud bottom set to first cloud instance
          icloudbot_tmp(ix,jy) = int(uvheight(kz)*eta_convert)
        icloudtop_tmp(ix,jy) = int(uvheight(kz)*eta_convert) !After the loop, icloudtop will be the top
#else
        if (icloudbot_tmp(ix,jy) .eq. icmv) &
            icloudbot_tmp(ix,jy) = (height(kz))
          icloudtop_tmp(ix,jy) = (height(kz))
#endif
      endif
    end do

  !**************************************************************************
  else       ! identify clouds using relative humidity
  !**************************************************************************
    do kz = 1,nz-1 ! Changed order of loop to prevent ETA computation to be done unnecessarily
      pressure=rho_tmp(ix,jy,kz)*r_air*tt_tmp(ix,jy,kz)
      rh=qv_tmp(ix,jy,kz)/f_qvsat(pressure,tt_tmp(ix,jy,kz))
      ! PS if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
      if (rh .ge. rhmin) then
#ifdef ETA
        if (icloudbot_tmp(ix,jy) .eq. icmv) then
          icloudbot_tmp(ix,jy)=int(uvheight(kz)*eta_convert)
        endif
        icloudtop_tmp(ix,jy)=int(uvheight(kz) *eta_convert)
#else
        if (icloudbot_tmp(ix,jy) .eq. icmv) then
          icloudbot_tmp(ix,jy)=(height(kz))
        endif
        icloudtop_tmp(ix,jy)=(height(kz))
#endif

      endif
    end do
!**************************************************************************
  endif ! lcw true/false
!**************************************************************************
end subroutine identify_cloud

subroutine apply_cloud_bounds(ix,jy,nxlim,nylim,lsprec_tmp,convprec_tmp,uvzlev, &
  icloudbot_tmp,icloudtop_tmp,max_cloudthck_eta,conv_clrange_eta, &
  highconvp_clrange_eta,lowconvp_clrange_eta)
  
  implicit none

  integer, intent(in) :: ix,jy,nxlim,nylim

  real,intent(in) :: lsprec_tmp(0:nxlim,0:nylim,numpf),convprec_tmp(0:nxlim,0:nylim,numpf)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev

  ! converted parameters for eta coordinates:
  integer,intent(in) :: max_cloudthck_eta
  integer,intent(in),dimension(2) :: conv_clrange_eta,highconvp_clrange_eta,lowconvp_clrange_eta

  integer,intent(inout) :: icloudbot_tmp(0:nxlim,0:nylim), icloudtop_tmp(0:nxlim,0:nylim)

  integer :: icloudtop_old, min_cloudtop
  integer :: kz
  real :: lsp,convp,prec
  logical lconvectprec


  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagn.


  ! memorise icloudtop
  icloudtop_old = icloudtop_tmp(ix,jy)
  ! top level, kz=nz-1
  ! limit cloud top to 19 km:
#ifdef ETA
  ! Enforce maximum cloudtop height in eta coords
  if (icloudbot_tmp(ix,jy) .eq. icmv) then !Can we skip to the next gridcell?
    icloudtop_tmp(ix,jy)=icmv
  else if (icloudtop_tmp(ix,jy) .lt. max_cloudthck_eta) then
    icloudtop_tmp(ix,jy) = max_cloudthck_eta
  endif

  ! To compute the minimum thickness in eta coordinates, this needs
  ! to be computed from the bottom level:
  ! 1) Convert icloudbot_tmp to meter coordinates
  ! 2) Add the min_cloudthck to the icloudbot_tmp(meter)
  ! 3) Convert this back to eta coordinates (min_cloudtop(eta))
  ! 4) Check if our cloudtop is higher (lower eta) than this minimum
  do kz=1,nz
    if (int(uvheight(kz)*eta_convert).le.icloudbot_tmp(ix,jy)) then
      min_cloudtop = uvzlev(ix,jy,kz)+min_cloudthck ! In meters
      exit
    endif
  end do
  do kz=1,nz
    if (int(uvzlev(ix,jy,kz)).gt.min_cloudtop) then ! back to eta
      min_cloudtop=int(uvheight(kz)*eta_convert)
      exit
    endif
  enddo
  ! PS  get rid of too thin clouds   
  if (icloudtop_tmp(ix,jy) .gt. min_cloudtop) then
    icloudbot_tmp(ix,jy)=icmv
    icloudtop_tmp(ix,jy)=icmv
  endif
  ! PS implement a rough fix for badly represented convection
  ! PS is based on looking at a limited set of comparison data
  lsp=  sum(   lsprec_tmp(ix,jy,:) )
  convp=sum( convprec_tmp(ix,jy,:) )
  prec=lsp+convp
  if (lsp.gt.convp) then !  prectype='lsp'
    lconvectprec = .false.
  else ! prectype='cp '
    lconvectprec = .true.
  endif
  if (lconvectprec .and. prec .gt. precmin .and.  &
    (icloudtop_old .gt. conv_clrange_eta(2) .or. &
      icloudbot_tmp(ix,jy) .lt. conv_clrange_eta(1)) ) then
    if (convp .lt. 0.1) then
      icloudbot_tmp(ix,jy) = lowconvp_clrange_eta(1)
      icloudtop_tmp(ix,jy) = lowconvp_clrange_eta(2)
    else
      icloudbot_tmp(ix,jy) = highconvp_clrange_eta(1)
      icloudtop_tmp(ix,jy) = highconvp_clrange_eta(2)
    endif
  endif

  !---------------------------------------------------------------------------------------
#else
  if (icloudbot_tmp(ix,jy) .eq. icmv) then ! if no bottom found, no top either
    icloudtop_tmp(ix,jy)=icmv
  else if (icloudtop_tmp(ix,jy) .gt. max_cloudthck) then ! max cloud height
    icloudtop_tmp(ix,jy) = max_cloudthck
  endif
 ! PS  get rid of too thin clouds
  if (icloudtop_tmp(ix,jy) .lt. icloudbot_tmp(ix,jy) + min_cloudthck) then
    icloudbot_tmp(ix,jy)=icmv
    icloudtop_tmp(ix,jy)=icmv
  endif
  ! PS implement a rough fix for badly represented convection
  ! PS is based on looking at a limited set of comparison data
  lsp=  sum(   lsprec_tmp(ix,jy,:) )
  convp=sum( convprec_tmp(ix,jy,:) )
  prec=lsp+convp
  if (lsp.gt.convp) then !  prectype='lsp'
    lconvectprec = .false.
  else ! prectype='cp '
    lconvectprec = .true.
  endif
  if (lconvectprec .and. prec .gt. precmin .and.  &
    (icloudtop_old .lt. conv_clrange(2) .or. &
      icloudbot_tmp(ix,jy) .gt. conv_clrange(1)) ) then
    if (convp .lt. 0.1) then
      icloudbot_tmp(ix,jy) = lowconvp_clrange(1)
      icloudtop_tmp(ix,jy) = lowconvp_clrange(2)
    else
      icloudbot_tmp(ix,jy) = highconvp_clrange(1)
      icloudtop_tmp(ix,jy) = highconvp_clrange(2)
    endif
  endif
#endif
end subroutine apply_cloud_bounds

subroutine verttransform_gfs(n,uuh,vvh,wwh,pvh)
  !                          i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !  CHANGES                                                                   *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !  Bernd C. Krueger, Feb. 2001:                                              *
  !   Variables tth and qvh (on eta coordinates) from common block             *
  !                                                                            *
  ! Sabine Eckhardt, March 2007:                                               *
  ! added the variable cloud for use with scavenging - descr. in com_mod       *
  ! PS/AT 2018/-21: variable "cloud" is replaced by quickfix, see below        * 
  !                                                                            *
  ! Unified ECMWF and GFS builds                                               *
  ! Marian Harustak, 12.5.2017                                                 *
  !     - Renamed from verttransform to verttransform_ecmwf                    *
  !                                                                            *
  !  undocumented modifications by NILU for v10                                *
  !                                                                            *
  !  Petra Seibert, 2018-06-13:                                                *
  !   - put back SAVE attribute for INIT, just to be safe                      *
  !   - minor changes, most of them just cosmetics                             *
  !   for details see changelog.txt in branch unive                            *
  !                                                                            *
  !  Petra Seibert, Anne Philipp, 2019-05-02: implement wetdepo quickfix       *
  !  Petra Seibert, Anne Tipka, 2020-11-19: reimplement in latest version      *
  !                                                                            *
  ! ****************************************************************************

  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! Note PS, AT 2021-01-29: all these fields are 0:nxmax-1,0:nymax-1 !!        *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! icloudbot(0:nxmax,0:nymax,numwfmem) cloud bottom field for wet deposition  * 
  ! icloudtop(0:nxmax,0:nymax,numwfmem) cloud thickness for wet deposition    *
  ! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
  ! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
  ! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
  !                                          [deltaeta/s]                      *
  ! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
  ! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
  ! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax)  :: wwh

  real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev
  real,dimension(nuvzmax) :: rhoh
  real,dimension(nwzmax) :: wzlev
  real,dimension(nzmax) :: pinmconv

! local array introduced in v10 by ?? to achieve better loop order (PS)
  integer,dimension(0:nxmax-1,0:nymax-1) :: idx

  integer :: icloudtop_old
  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym,kz_inv
  real :: clw ! cloud water in kg / m2 in a grid cell
  real :: pressure,rh,lsp,convp,prec

  real :: pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2,cosf
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  
  real,parameter :: const=r_air/ga
  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

  logical lconvectprec
  logical, save :: init = .true. !PS 2018-06-13: add back save attribute

  ! NCEP version
  integer :: llev, i

  integer :: rank_thread , nr_threads 


  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given.                                                *
  !*************************************************************************

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************
    call verttransform_init(n)
  
  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.

  endif


  ! Loop over the whole grid
  !*************************
!$OMP DO
  do jy=0,nymin1
    do ix=0,nxmin1


!   if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 1' 

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
        if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      enddo
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version


  ! compute height of pressure levels above ground
  !***********************************************

      tvold=tth(ix,jy,llev,n)*(1.+0.608*qvh(ix,jy,llev,n))
      pold=akz(llev)
      wzlev(llev)=0.
      uvwzlev(ix,jy,llev)=0.
      rhoh(llev)=pold/(r_air*tvold)

!    if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 2'

      do kz=llev+1,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
        rhoh(kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)* &
          (tv-tvold)/log(tv/tvold)
        else
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)*tv
        endif
        wzlev(kz)=uvwzlev(ix,jy,kz)

        tvold=tv
        pold=pint
      enddo

  ! pinmconv=(h2-h1)/(p2-p1)
! if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 3'

      pinmconv(llev)=(uvwzlev(ix,jy,llev+1)-uvwzlev(ix,jy,llev))/ &
           ((aknew(llev+1)+bknew(llev+1)*ps(ix,jy,1,n))- &
           (aknew(llev)+bknew(llev)*ps(ix,jy,1,n)))
      do kz=llev+1,nz-1
        pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps(ix,jy,1,n))- &
             (aknew(kz-1)+bknew(kz-1)*ps(ix,jy,1,n)))
      enddo
      pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps(ix,jy,1,n))- &
           (aknew(nz-1)+bknew(nz-1)*ps(ix,jy,1,n)))


  ! Levels, where u,v,t and q are given
  !************************************
! if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 4'


      uu(ix,jy,1,n)=uuh(ix,jy,llev)
      vv(ix,jy,1,n)=vvh(ix,jy,llev)
      tt(ix,jy,1,n)=tth(ix,jy,llev,n)
      qv(ix,jy,1,n)=qvh(ix,jy,llev,n)

      
  ! IP & SEC, 201812 add clouds
      if (lcw) then
         clwc(ix,jy,1,n)=clwch(ix,jy,llev,n)
      endif 
      pv(ix,jy,1,n)=pvh(ix,jy,llev)
      rho(ix,jy,1,n)=rhoh(llev)
      pplev(ix,jy,1,n)=akz(llev)
      uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
      vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
      tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
      qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
  ! IP & SEC, 201812 add clouds

   
      if (lcw) then
         clwc(ix,jy,nz,n)=clwch(ix,jy,nuvz,n)
      endif
      pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
      rho(ix,jy,nz,n)=rhoh(nuvz)
      pplev(ix,jy,nz,n)=akz(nuvz)
      kmin=llev+1

       
      do iz=2,nz-1
        do kz=kmin,nuvz
          ! print*, 'in loop 4.3.1', jy, ix, iz, kz  
          if (height(iz).gt.uvwzlev(ix,jy,nuvz)) then 
            uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
            vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
            tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
            qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  ! IP & SEC, 201812 add clouds
            if (lcw) then
               clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            endif
            pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
            rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
            pplev(ix,jy,iz,n)=pplev(ix,jy,nz,n)
            exit
          endif
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
           (height(iz).le.uvwzlev(ix,jy,kz))) then
            !real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
            vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
            tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2+tth(ix,jy,kz,n)*dz1)/dz
            qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2+qvh(ix,jy,kz,n)*dz1)/dz

  ! IP & SEC, 201812 add clouds
            if (lcw) then
              clwc(ix,jy,iz,n)= &
                (clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
            endif
            pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz

            rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz

            pplev(ix,jy,iz,n)=(akz(kz-1)*dz2+akz(kz)*dz1)/dz

          endif
        enddo
      enddo

 
  ! Interpolation of vertical motion (levels where w is given)
  !***********************************************************

      ww(ix,jy,1,n)=wwh(ix,jy,llev)*pinmconv(llev)
      ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
      kmin=llev+1
      do iz=2,nz
        do kz=kmin,nwz
          if ((height(iz).gt.wzlev(kz-1)).and. &
           (height(iz).le.wzlev(kz))) then
            dz1=height(iz)-wzlev(kz-1)
            dz2=wzlev(kz)-height(iz)
            dz=dz1+dz2
            ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(kz-1)*dz2 &
              +wwh(ix,jy,kz)*pinmconv(kz)*dz1)/dz
          endif
        enddo
      enddo

!if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 6'
  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/(height(2)-height(1))
      do kz=2,nz-1
        drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
          (height(kz+1)-height(kz-1))
      enddo
      drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    !if ((jy.eq.0).and.(ix.eq.0)) 


    enddo
  enddo
!$OMP END DO

  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************


!$OMP DO
  do jy=1,ny-2
    cosf=cos((real(jy)*dy+ylat0)*pi180)

    do ix=1,nx-2
   ! print*, '1] slope of eta levels jy, ix=',jy, ix 

   !if ((jy.le.2).and.(ix.eq.1)) print*, 'in eta loop 1'

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
       if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
      llev = llev+1

!if ((jy.le.2).and.(ix.eq.1)) print*, 'in eta loop 2'


      if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version

      kmin=llev+1

      !if ((jy.le.3).and.(ix.gt.250)) print*, 'in eta loop 3'

      do iz=2,nz-1

        ui=uu(ix,jy,iz,n)*dxconst/cosf
        vi=vv(ix,jy,iz,n)*dyconst

        klp=nz+1
        do kz=kmin,nz
 
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            exit
          endif

        enddo

        if (klp.eq.nz+1) then
          klp=nz
          kl=nz-1

         ! real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev
         
          dz1=uvwzlev(ix,jy,klp)-uvwzlev(ix,jy,kl)


          dz2=0.
        endif

        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))*0.5
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))*0.5
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))*0.5
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))*0.5
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi)

      enddo ! z

      !if ((jy.le.3).and.(ix.eq.1)) print*, 'in eta loop end z jy=',jy


    enddo


  enddo
!$OMP END DO


  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        enddo
      enddo
    enddo


    do iz=1,nz

      ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      elseif (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

      ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=90.
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)
      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      enddo
      
    enddo


    ! Fix: Set W (vertical wind) at pole to the zonally averaged W of the next 
    ! equator-ward parallel 


    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      enddo
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      enddo
    enddo

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
            vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        enddo
      enddo
    enddo

    do iz=1,nz

      ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))+xlonr
      elseif (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

      ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=-90.
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      enddo
      
    enddo


    ! Fix: Set W at pole to the zonally averaged W of the next equator-
    ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      enddo
    enddo
    
  endif


  ! PS, AT: for v10.5, we add back the quick fix to interpolate clouds in 
  !   interpol_rain.f90 developed by PS for v8 and extend it to using 
  !   cloud water fields

  !*******************************************************************************
  if (lcw) then ! identify clouds based on cloud water content
  !*******************************************************************************

    write(*,*) 'Global NCEP fields: using cloud water'
    
    ctwc(:,:,n)=0. ! initialise cloud total water content

    ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    if (.not. lcwsum) clwc(:,:,:,n) = clwc(:,:,:,n) + ciwc(:,:,:,n)

    do kz = 1,nz-1
      do jy=0,nymin1
        do ix=0,nxmin1
          if (kz .eq. 1) then
            icloudbot(ix,jy,n) = icmv
!!          icloudtop=icmv ! this is just a local variable
!           we will use icloudtop as workspace for cloud top
          endif

          ! vertically integrate cloud water and determine cloud bottom, top
          ! cloud water per cell in kg / m2
          ! calculate cloud water mass per area: kgCW/kgAIR * kgAIR/m3 * m = kgCW/m2

          clw = clwc(ix,jy,kz,n)*rho(ix,jy,kz,n)*(height(kz+1)-height(kz)) 
          ! Add this layer to column cloud water [m3/m3]
          ctwc(ix,jy,n) = ctwc(ix,jy,n)+clw ! kg / m2 in column

          if (clw .gt. 0.) then ! cloud layer - maybe use threshold?
            if (icloudbot(ix,jy,n) .eq. icmv) &
              icloudbot(ix,jy,n) = nint(height(kz))
            icloudtop(ix,jy,n) = nint(height(kz))
          endif

          if (kz .eq. nz-1) then ! top level
            ! memorise icloudtop
            icloudtop_old = icloudtop(ix,jy,n)
            ! limit cloud top to 19 km:
            if (icloudtop(ix,jy,n) .gt. 19000) icloudtop(ix,jy,n) = 19000 
            if (icloudbot(ix,jy,n) .eq. icmv) then
              icloudtop(ix,jy,n) = icmv
            endif

           ! PS  get rid of too thin clouds      
            if (icloudtop(ix,jy,n) .lt. 50) then
              icloudbot(ix,jy,n)=icmv
              icloudtop(ix,jy,n)=icmv
            endif
            
            ! PS implement a rough fix for badly represented convection
            ! PS is based on looking at a limited set of comparison data
            lsp=  sum(   lsprec(ix,jy,1,:,n) )
            convp=sum( convprec(ix,jy,1,:,n) )
            prec=lsp+convp
            if (lsp.gt.convp) then !  prectype='lsp'
              lconvectprec = .false.
            else ! prectype='cp '
              lconvectprec = .true.
            endif
            if (lconvectprec .and. prec .gt. precmin .and.  &
              (icloudtop_old .lt. 6000 .or. icloudbot(ix,jy,n) .gt. 3000) ) then
              if (convp .lt. 0.1) then
                icloudbot(ix,jy,n) = 500
                icloudtop(ix,jy,n) =         8000
              else
                icloudbot(ix,jy,n) = 0
                icloudtop(ix,jy,n) =      10000
              endif
            endif
          endif ! end top level

        enddo ! ix loop
      enddo ! jy loop
    enddo ! kz loop


!**************************************************************************
  else       ! identify clouds using relative humidity
!**************************************************************************
!   clouds occur where rh>90% (using rh_ice for T<-20 deg C)

    write(*,*) 'NCEP fields: using relative humidity for cloud &
        &identification'
    do kz = 1,nz-1
      do jy=0,nymin1
        do ix=0,nxmin1

!PS       note that original by Sabine Eckhart was 80%
!PS       however, for T<-20 C we consider saturation over ice
!PS       so I think 90% should be enough          
          if (kz .eq. 1) then
            icloudbot(ix,jy,n) = icmv
!!          icloudtop=icmv ! this is just a local variable
!           we will use icloudtop as workspace for cloud top
          endif
!98        do kz=1,nz
          pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
          rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
!PS            if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
          if (rh .ge. rhmin) then
            if (icloudbot(ix,jy,n) .eq. icmv) then
              icloudbot(ix,jy,n)=nint(height(kz))! use int to save memory
            endif
            icloudtop(ix,jy,n)=nint(height(kz)) ! use int to save memory
          endif
!          enddo

!PS/AT 2021: in this version, we skip the iteration with smaller rhmin
!PS try to get a cloud thicker than 50 m 
!PS if there is at least .01 mm/h  - changed to 0.002 and put into
!PS parameter precpmin        
!          if ((icloudbot(ix,jy,n) .eq. icmv .or. &
!            icloudtop-icloudbot(ix,jy,n) .lt. 50) .and. prec .gt. precmin) then
!              rhmin = rhmin - 0.05
!              if (rhmin .ge. 0.30) goto 98 ! give up for <= 25% rel.hum.
!          endif
!          if (icloudtop .ne. icmv) then
!            icloudtop(ix,jy,n) = icloudtop-icloudbot(ix,jy,n)
!          else
!            icloudtop(ix,jy,n) = icmv
!          endif
    
          if (kz .eq. nz-1) then ! top level
          
            ! memorise icloudtop
            icloudtop_old = icloudtop(ix,jy,n)
            ! limit cloud top to 19 km:
            if (icloudtop(ix,jy,n) .gt. 19000) icloudtop(ix,jy,n) = 19000 
            if (icloudbot(ix,jy,n) .ne. icmv) then
              icloudtop(ix,jy,n) = icloudtop(ix,jy,n)-icloudbot(ix,jy,n)
            else
              icloudtop(ix,jy,n) = icmv
            endif

           ! PS  get rid of too thin clouds      
            if (icloudtop(ix,jy,n) .lt. 50) then
              icloudbot(ix,jy,n)=icmv
              icloudtop(ix,jy,n)=icmv
            endif
            
            ! PS implement a rough fix for badly represented convection
            ! PS is based on looking at a limited set of comparison data
            lsp=  sum(   lsprec(ix,jy,1,:,n) )
            convp=sum( convprec(ix,jy,1,:,n) )
            prec=lsp+convp
            if (lsp.gt.convp) then !  prectype='lsp'
              lconvectprec = .false.
            else ! prectype='cp '
              lconvectprec = .true.
            endif
            if (lconvectprec .and. prec .gt. precmin .and.  &
              (icloudtop_old .lt. 6000 .or. icloudbot(ix,jy,n) .gt. 3000) ) then
              if (convp .lt. 0.1) then
                icloudbot(ix,jy,n) = 500
                icloudtop(ix,jy,n) = 8000
              endif
            else
              icloudbot(ix,jy,n) = 0
              icloudtop(ix,jy,n) = 10000
            endif
            
          endif ! end top level

        enddo ! ix loop
      enddo ! jy loop
    enddo ! kz loop

!**************************************************************************
  endif ! lcw true/false
!**************************************************************************

end subroutine verttransform_gfs

subroutine verttransform_ecmwf_heights(nxlim,nylim, &
  tt2_tmp,td2_tmp,ps_tmp,qvh_tmp,tth_tmp,prsh_tmp, &
  rhoh_tmp,pinmconv,uvzlev,wzlev)

  implicit none

  integer, intent(in) :: nxlim,nylim
  real,intent(in),dimension(0:nxlim,0:nylim) :: tt2_tmp,td2_tmp,ps_tmp
  real,intent(in),dimension(0:nxlim,0:nylim,nuvzmax) :: qvh_tmp,tth_tmp
  real,intent(out),dimension(0:nxlim,0:nylim,nuvzmax) :: rhoh_tmp,prsh_tmp
  real,intent(out),dimension(0:nxlim,0:nylim,nzmax) :: pinmconv
  real,intent(out),dimension(0:nxlim,0:nylim,nuvzmax) :: uvzlev,wzlev
  !real,dimension(0:nxlim,0:nylim) :: tvold,pold,pint,tv
  real :: tvold,pold,pint,tv
  real,parameter :: const=r_air/ga
  integer :: ix,jy,kz

  ! Loop over the whole grid
  !*************************
!$OMP PARALLEL PRIVATE(jy,ix,pint,tv,tvold,pold,kz)
!$OMP DO
  do jy=0,nylim
    do ix=0,nxlim
      tvold=tt2_tmp(ix,jy)*(1.+0.378*ew(td2_tmp(ix,jy),ps_tmp(ix,jy))/ &
           ps_tmp(ix,jy))
      pold=ps_tmp(ix,jy)
      uvzlev(ix,jy,1)=0.
      wzlev(ix,jy,1)=0.
      rhoh_tmp(ix,jy,1)=pold/(r_air*tvold)
      prsh_tmp(ix,jy,1)=ps_tmp(ix,jy)

      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*ps_tmp(ix,jy)
        prsh_tmp(ix,jy,kz)=pint
        tv=tth_tmp(ix,jy,kz)*(1.+0.608*qvh_tmp(ix,jy,kz))
        rhoh_tmp(ix,jy,kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvzlev(ix,jy,kz)=uvzlev(ix,jy,kz-1)+const* &
               log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          uvzlev(ix,jy,kz)=uvzlev(ix,jy,kz-1)+const* &
               log(pold/pint)*tv
        endif

        tvold=tv
        pold=pint

      end do

      do kz=2,nwz-1
        wzlev(ix,jy,kz)=(uvzlev(ix,jy,kz+1)+uvzlev(ix,jy,kz))*0.5
      end do
      wzlev(ix,jy,nwz)=wzlev(ix,jy,nwz-1)+ &
           uvzlev(ix,jy,nuvz)-uvzlev(ix,jy,nuvz-1)


      pinmconv(ix,jy,1)=(uvzlev(ix,jy,2))/ &
           ((aknew(2)+bknew(2)*ps_tmp(ix,jy))- &
           (aknew(1)+bknew(1)*ps_tmp(ix,jy)))
      do kz=2,nz-1
        pinmconv(ix,jy,kz)=(uvzlev(ix,jy,kz+1)-uvzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps_tmp(ix,jy))- &
             (aknew(kz-1)+bknew(kz-1)*ps_tmp(ix,jy)))
      end do
      pinmconv(ix,jy,nz)=(uvzlev(ix,jy,nz)-uvzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps_tmp(ix,jy))- &
           (aknew(nz-1)+bknew(nz-1)*ps_tmp(ix,jy)))
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  ! pold(:,:)=ps_tmp(:,:)
  ! uvzlev(:,:,1)=0.
  ! wzlev(:,:,1)=0.
  ! rhoh_tmp(:,:,1)=pold(:,:)/(r_air*tvold(:,:))
  ! prsh_tmp(:,:,1)=ps_tmp(:,:)

  ! Compute heights of eta levels
  !******************************

  ! do kz=2,nuvz
  !   pint(:,:)=akz(kz)+bkz(kz)*ps_tmp(:,:)
  !   prsh_tmp(:,:,kz)=pint(:,:)
  !   tv(:,:)=tth_tmp(:,:,kz)*(1.+0.608*qvh_tmp(:,:,kz))
  !   rhoh_tmp(:,:,kz)=pint(:,:)/(r_air*tv(:,:))

  !   where (abs(tv(:,:)-tvold(:,:)).gt.0.2) 
  !     uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*&
  !          &log(pold(:,:)/pint(:,:))* &
  !          (tv(:,:)-tvold(:,:))/&
  !          &log(tv(:,:)/tvold(:,:))
  !   elsewhere
  !     uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*&
  !          &log(pold(:,:)/pint(:,:))*tv(:,:)
  !   endwhere

  !   tvold(:,:)=tv(:,:)
  !   pold(:,:)=pint(:,:)

  ! end do

  ! do kz=2,nwz-1
  !   wzlev(:,:,kz)=(uvzlev(:,:,kz+1)+uvzlev(:,:,kz))/2.
  ! end do
  ! wzlev(:,:,nwz)=wzlev(:,:,nwz-1)+ &
  !      uvzlev(:,:,nuvz)-uvzlev(:,:,nuvz-1)


  ! pinmconv(:,:,1)=(uvzlev(:,:,2))/ &
  !      ((aknew(2)+bknew(2)*ps_tmp(:,:))- &
  !      (aknew(1)+bknew(1)*ps_tmp(:,:)))
  ! do kz=2,nz-1
  !   pinmconv(:,:,kz)=(uvzlev(:,:,kz+1)-uvzlev(:,:,kz-1))/ &
  !        ((aknew(kz+1)+bknew(kz+1)*ps_tmp(:,:))- &
  !        (aknew(kz-1)+bknew(kz-1)*ps_tmp(:,:)))
  ! end do
  ! pinmconv(:,:,nz)=(uvzlev(:,:,nz)-uvzlev(:,:,nz-1))/ &
  !      ((aknew(nz)+bknew(nz)*ps_tmp(:,:))- &
  !      (aknew(nz-1)+bknew(nz-1)*ps_tmp(:,:)))
end subroutine verttransform_ecmwf_heights

subroutine verttransform_ecmwf_windfields_nest(l,n, &
  uuhn,vvhn,wwhn,pvhn,rhohn,prshn,pinmconv)

  implicit none

  integer,intent(in) :: l,n
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numbnests) :: &
    uuhn,vvhn,pvhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nwzmax,numbnests) :: wwhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: rhohn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: prshn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nzmax) :: pinmconv
  real,dimension(0:nymaxn-1) :: cosf

  integer,dimension(0:nxmaxn-1,0:nymaxn-1) :: idx

  integer :: ix,jy,kz,iz,ix1,jy1,ixp,jyp

  real :: dz1,dz2,dz,dpdeta
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  integer :: nxm1, nym1

  nxm1=nxn(l)-1 
  nym1=nyn(l)-1 

  ! Levels, where u,v,t and q are given
  !************************************
!$OMP PARALLEL PRIVATE(jy,ix,kz,dz1,dz2,dz,ix1,jy1,ixp,jyp,dzdx1,dzdx2,dzdx, &
!$OMP dzdy1,dzdy2,dzdy,dpdeta)

!$OMP DO
  do jy=0,nym1
    do ix=0,nxm1
      uun(ix,jy,1,n,l)=uuhn(ix,jy,1,l)
      vvn(ix,jy,1,n,l)=vvhn(ix,jy,1,l)
      ttn(ix,jy,1,n,l)=tthn(ix,jy,1,n,l)
#ifndef ETA
      qvn(ix,jy,1,n,l)=qvhn(ix,jy,1,n,l)
#endif
      if (lcw_nest(l)) then
        clwcn(ix,jy,1,n,l)=clwchn(ix,jy,1,n,l)
        if (.not.lcwsum_nest(l)) ciwcn(ix,jy,1,n,l)=ciwchn(ix,jy,1,n,l)
      end if
      pvn(ix,jy,1,n,l)=pvhn(ix,jy,1,l)
      rhon(ix,jy,1,n,l)=rhohn(ix,jy,1)
      prsn(ix,jy,1,n,l)=prshn(ix,jy,1)

      uun(ix,jy,nz,n,l)=uuhn(ix,jy,nuvz,l)
      vvn(ix,jy,nz,n,l)=vvhn(ix,jy,nuvz,l)
      ttn(ix,jy,nz,n,l)=tthn(ix,jy,nuvz,n,l)
#ifndef ETA
      qvn(ix,jy,nz,n,l)=qvhn(ix,jy,nuvz,n,l)
      if (lcw_nest(l)) then
        clwcn(ix,jy,nz,n,l)=clwchn(ix,jy,nuvz,n,l)
        if (.not.lcwsum_nest(l)) ciwcn(ix,jy,nz,n,l)=ciwchn(ix,jy,nuvz,n,l)
      endif
#endif
      pvn(ix,jy,nz,n,l)=pvhn(ix,jy,nuvz,l)
      rhon(ix,jy,nz,n,l)=rhohn(ix,jy,nuvz)
      prsn(ix,jy,nz,n,l)=prshn(ix,jy,nuvz)

      idx(ix,jy)=2
    enddo
  enddo
!$OMP END DO

  do iz=2,nz-1
!$OMP DO SCHEDULE(dynamic)
    do jy=0,nym1
      do ix=0,nxm1
        if(height(iz).gt.etauvheightn(ix,jy,nuvz,n,l)) then
          uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
          vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
          ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
          pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
#ifndef ETA
          qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
          !hg adding the cloud water
          if (lcw_nest(l)) then
            clwcn(ix,jy,iz,n,l)=clwcn(ix,jy,nz,n,l)
            if (.not.lcwsum_nest(l)) ciwcn(ix,jy,iz,n,l)=ciwcn(ix,jy,nz,n,l)
          endif
#endif
          rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
          prsn(ix,jy,iz,n,l)=prsn(ix,jy,nz,n,l)
        else
          innuvz: do kz=idx(ix,jy),nuvz
            if ((idx(ix,jy).lt.kz).and. & 
              (height(iz).gt.etauvheightn(ix,jy,kz-1,n,l)).and. &
              (height(iz).le.etauvheightn(ix,jy,kz,n,l))) then
              idx(ix,jy)=kz
              exit innuvz
            endif
          enddo innuvz
        endif

        if(height(iz).le.etauvheightn(ix,jy,nuvz,n,l)) then
          kz=idx(ix,jy)
          dz1=height(iz)-etauvheightn(ix,jy,kz-1,n,l)
          dz2=etauvheightn(ix,jy,kz,n,l)-height(iz)
          dz=dz1+dz2
          uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+uuhn(ix,jy,kz,l)*dz1)/dz
          vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+vvhn(ix,jy,kz,l)*dz1)/dz
          ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2 &
               +tthn(ix,jy,kz,n,l)*dz1)/dz
          pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+pvhn(ix,jy,kz,l)*dz1)/dz
#ifndef ETA
          qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2 &
               +qvhn(ix,jy,kz,n,l)*dz1)/dz
          !hg adding the cloud water
          if (lcw_nest(l)) then
            clwcn(ix,jy,iz,n,l)=(clwchn(ix,jy,kz-1,n,l)*dz2+clwchn(ix,jy,kz,n,l)*dz1)/dz
            if (.not.lcwsum_nest(l)) ciwcn(ix,jy,iz,n,l) = &
              (ciwchn(ix,jy,kz-1,n,l)*dz2+ciwchn(ix,jy,kz,n,l)*dz1)/dz
          end if
#endif
          rhon(ix,jy,iz,n,l)=(rhohn(ix,jy,kz-1)*dz2+rhohn(ix,jy,kz)*dz1)/dz
          prsn(ix,jy,iz,n,l)=(prshn(ix,jy,kz-1)*dz2+prshn(ix,jy,kz)*dz1)/dz
        endif
      enddo
    enddo
!$OMP END DO
!$OMP BARRIER
  enddo

  ! Interpolation of vertical motion (levels where w is given)
  !***********************************************************

!$OMP DO
  do jy=0,nym1
    do ix=0,nxm1
      idx(ix,jy)=2
      wwn(ix,jy,1,n,l)=wwhn(ix,jy,1,l)*pinmconv(ix,jy,1)
      wwn(ix,jy,nz,n,l)=wwhn(ix,jy,nwz,l)*pinmconv(ix,jy,nz)
    enddo
  enddo
!$OMP END DO

  do iz=2,nz-1
!$OMP DO SCHEDULE(dynamic)
    do jy=0,nym1
      do ix=0,nxm1

        inn: do kz=idx(ix,jy),nwz
          if((idx(ix,jy).le.kz) .and. &
            (height(iz).gt.etawheightn(ix,jy,kz-1,n,l)).and. &
            (height(iz).le.etawheightn(ix,jy,kz,n,l))) then
            idx(ix,jy)=kz
            exit inn
          endif
        enddo inn

        kz=idx(ix,jy)
        dz1=height(iz)-etawheightn(ix,jy,kz-1,n,l)
        dz2=etawheightn(ix,jy,kz,n,l)-height(iz)
        dz=dz1+dz2
        wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(ix,jy,kz-1)*dz2 &
             +wwhn(ix,jy,kz,l)*pinmconv(ix,jy,kz)*dz1)/dz
        drhodzn(ix,jy,iz,n,l)=(rhon(ix,jy,iz+1,n,l)-rhon(ix,jy,iz-1,n,l))/ &
             (height(iz+1)-height(iz-1))
      enddo
    enddo
!$OMP END DO
!$OMP BARRIER
  enddo

  ! Compute density gradients at intermediate levels
  !*************************************************
!$OMP DO
  do jy=0,nym1
    do ix=0,nxm1
      drhodzn(ix,jy,nz,n,l)=drhodzn(ix,jy,nz-1,n,l)
      drhodzn(ix,jy,1,n,l)=(rhon(ix,jy,2,n,l)-rhon(ix,jy,1,n,l))/ &
           (height(2)-height(1))
    end do
  end do
!$OMP END DO

  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

!$OMP DO
  do jy=1,nyn(l)-2
    cosf(jy)=1./cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
    do ix=1,nxn(l)-2
      idx(ix,jy)=2
    enddo
  enddo
!$OMP END DO

  do iz=2,nz-1
!$OMP DO SCHEDULE(guided)
    do jy=1,nyn(l)-2
      do ix=1,nxn(l)-2

        inneta: do kz=idx(ix,jy),nz
          if((idx(ix,jy).le.kz).and. &
            (height(iz).gt.etauvheightn(ix,jy,kz-1,n,l)).and. &
            (height(iz).le.etauvheightn(ix,jy,kz,n,l))) then
            idx(ix,jy)=kz
            exit inneta
          endif
        enddo inneta

        kz=idx(ix,jy)
        dz1=height(iz)-etauvheightn(ix,jy,kz-1,n,l)
        dz2=etauvheightn(ix,jy,kz,n,l)-height(iz)
        dz=dz1+dz2
        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(etauvheightn(ixp,jy,kz-1,n,l)-etauvheightn(ix1,jy,kz-1,n,l))*0.5
        dzdx2=(etauvheightn(ixp,jy,kz,n,l)-etauvheightn(ix1,jy,kz,n,l))*0.5
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(etauvheightn(ix,jyp,kz-1,n,l)-etauvheightn(ix,jy1,kz-1,n,l))*0.5
        dzdy2=(etauvheightn(ix,jyp,kz,n,l)-etauvheightn(ix,jy1,kz,n,l))*0.5
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l) + &
          (dzdx*uun(ix,jy,iz,n,l)*dxconst*xresoln(l)*cosf(jy)+ &
          dzdy*vvn(ix,jy,iz,n,l)*dyconst*yresoln(l))

      enddo
    enddo
!$OMP END DO
!$OMP BARRIER
  enddo

  ! Keep original fields if wind_coord_type==ETA
#ifdef ETA
!$OMP DO

    do kz=1,nz
      do jy=0,nym1
        do ix=0,nxm1
          uuetan(ix,jy,kz,n,l) = uuhn(ix,jy,kz,l)
          vvetan(ix,jy,kz,n,l) = vvhn(ix,jy,kz,l)
          ttetan(ix,jy,kz,n,l) = tthn(ix,jy,kz,n,l)
          qvn(ix,jy,kz,n,l) = qvhn(ix,jy,kz,n,l)
          pvetan(ix,jy,kz,n,l) = pvhn(ix,jy,kz,l)
          rhoetan(ix,jy,kz,n,l) = rhohn(ix,jy,kz)
          prsetan(ix,jy,kz,n,l) = prshn(ix,jy,kz)
          ! eq A11 from Mid-latitude atmospheric dynamics by Jonathan E. Martin
          ! tvirtualn(ix,jy,kz,n,l)=ttetan(ix,jy,kz,n,l)* &  
          !   ((qvn(ix,jy,kz,n,l)+0.622)/(0.622*qvn(ix,jy,kz,n,l)+0.622))
          if ((kz.gt.1).and.(kz.lt.nz)) drhodzetan(ix,jy,kz,n,l)= &
            (rhohn(ix,jy,kz+1)-rhohn(ix,jy,kz-1))/(height(kz+1)-height(kz-1))
          if (lcw) then
            clwcn(ix,jy,kz,n,l)=clwchn(ix,jy,kz,n,l)
            if (.not.lcwsum_nest(l)) ciwcn(ix,jy,kz,n,l)=ciwchn(ix,jy,kz,n,l)
          endif
        end do
      end do
    end do
!$OMP END DO

!$OMP DO
    do jy=0,nym1
      do ix=0,nxm1
        drhodzetan(ix,jy,1,n,l)=(rhoetan(ix,jy,2,n,l)-rhoetan(ix,jy,1,n,l))/ &
             (height(2)-height(1))
        drhodzetan(ix,jy,nz,n,l)=drhodzetan(ix,jy,nz-1,n,l)
        ! tvirtualn(ix,jy,1,n,l)=tt2n(ix,jy,1,n,l)* &
        !   (1.+0.378*ew(td2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))/ps(ix,jy,1,n,l))
        
        ! Convert w from Pa/s to eta/s, following FLEXTRA
        !************************************************
        do kz=1,nuvz-1
          if (kz.eq.1) then
            dpdeta=(akm(kz+1)-akm(kz)+(bkm(kz+1)-bkm(kz))*ps(ix,jy,1,n))/ &
              (wheight(kz+1)-wheight(kz))
          else if (kz.eq.nuvz-1) then
            dpdeta=(akm(kz)-akm(kz-1)+(bkm(kz)-bkm(kz-1))*ps(ix,jy,1,n))/ &
              (wheight(kz)-wheight(kz-1))
          else
            dpdeta=(akm(kz+1)-akm(kz-1)+(bkm(kz+1)-bkm(kz-1))*ps(ix,jy,1,n))/ &
              (wheight(kz+1)-wheight(kz-1))
          endif
          wwetan(ix,jy,kz,n,l)=wwhn(ix,jy,kz,l)/dpdeta
        enddo
        wwetan(ix,jy,nuvz,n,l)=wwetan(ix,jy,nuvz-1,n,l)
      enddo
    enddo 
!$OMP END DO
#endif
!$OMP END PARALLEL
end subroutine verttransform_ecmwf_windfields_nest

end module verttransform_mod
