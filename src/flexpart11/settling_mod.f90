! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later
module settling_mod
  
  implicit none

  integer :: re_nsteps ! Number of Reynolds number table points
  ! real, allocatable ,dimension(:) :: re_lookup, cd_lookup


  private :: viscosity
  public :: get_settling
contains

! subroutine init_dragcoeff_lookup()
!   implicit none
!   !*****************************
!   ! Clift and Guavin 1971 model*
!   !*****************************
!   real :: re_min, re_max ! Range of Reynolds numbers
!   integer :: i

!   re_nsteps=1000

!   allocate (re_lookup(re_nsteps), cd_lookup(re_nsteps))

!   ! Reynolds numbers range between 0.2 and 8000. 
!   ! Below this range, the approx of re/24 is sufficient and above 8000, the 
!   ! full equation is executed during run-time.
!   ! logarithmic spacing with a 1000 steps gives a maximum deviation of drag_coeff of
!   ! 0.07% as compared to using the full equation (mean=0.0037% deviation)
!   !**********************************************************************
!   re_min=log10(0.02)
!   re_max=log10(8000.)
!   do i=1,re_nsteps
!     re_lookup(i) = 10.**((re_max-re_min) * real(i-1) / real(re_nsteps-1) + re_min)
!   end do

!   ! Computing drag coefficients
!   !****************************
!   cd_lookup=(24.0/re_lookup)*(1+0.15*(re_lookup**0.687))+ &
!       0.42/(1.0+42500.0/(re_lookup**1.16))

! end subroutine init_dragcoeff_lookup

! subroutine find_dragcoeff(reynolds, c_drag)
!   implicit none
!   real, intent(in) :: reynolds
!   real, intent(out) :: c_drag
!   integer :: i, i_re, minsteps
!   real :: dre,dre1,dre2
!   ! This lookup table is extremely slow. For now it will just return
!   ! the function instead.


!   ! If reynolds<0.2, approximation 24/reynolds is valid
!   ! 24/reynolds makes up >99% of all components
!   !****************************************************
!   if (reynolds.le.0.02) then
!     c_drag=(24.0/reynolds)

!   else if (reynolds.gt.8000) then ! Outside of lookup table range
!     c_drag=(24.0/reynolds)*(1+0.15*(reynolds**0.687))+ &
!       0.42/(1.0+42500.0/(reynolds**1.16))

!   else
!     ! Linear search for correct indices in lookup table
!     !**************************************************
!     i_re = re_nsteps
!     do i = 0, re_nsteps
!       if (re_lookup(i).ge.reynolds) then
!         i_re=i
!         exit
!       endif
!     end do

!     ! Linear interpolation (maybe also do logarithmic?)
!     !*********************
!     if (i_re.eq.re_nsteps) then
!       c_drag=cd_lookup(re_nsteps)
!     else if (i_re.eq.1) then
!       c_drag=cd_lookup(1)
!     else
!       dre=1./(cd_lookup(i_re+1)-cd_lookup(i_re))
!       dre1=(reynolds-re_lookup(i_re))*dre
!       dre2=(re_lookup(i_re+1)-reynolds)*dre
!       c_drag=dre1*cd_lookup(i_re+1)+dre2*cd_lookup(i_re)
!     endif
!   endif

! end subroutine

subroutine get_settling(xt,yt,zt,nsp,settling)
  !                          i   i  i  i   i     o
  !*****************************************************************************
  !                                                                            *
  !  This subroutine calculates particle settling velocity.                    *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    May 2010                                                                *
  !                                                                            *
  !  Improvement over traditional settling calculation in FLEXPART:            *
  !  generalize to higher Reynolds numbers and also take into account the      *
  !  temperature dependence of dynamic viscosity.                              *
  !                                                                            *
  !  Based on:                                                                 *
  !  Naeslund E., and Thaning, L. (1991): On the settling velocity in a        *
  !  nonstationary atmosphere, Aerosol Science and Technology 14, 247-256.     *
  !                                                                            *
  !  Changes                                                                   *
  !     Daria Tatsii 2022: implementation of shape factor according to         *
  !                        Bagheri & Bonadonna 2016                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zt           coordinates position for which wind data shall be cal-  *
  !                    culated                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  ! dfdr               fluid density/particle density                          *
  ! Veq [m^3]           equivalent volume of a sphere                          *
  ! dcyl [m]           diameter of a cylinder (fiber)                          * 
  ! f                  flatness parameters, S/I                               *
  ! e                  elongation parameters, I/L                              *
  ! Fs                 Stokes form factor, f e^1.3                             *
  ! Fn                 Newton's  form factor                                   *
  ! Ks                 Stokes' drag correction                                 * 
  ! vsp                help variable                                           *
  ! x                  aspect ratio of cylinder height to its diameter         *
  !                                                                            *
  ! Variables:                                                                 *
  ! c_d                drag coefficient                                        *
  ! settling [m/s]     settling velocity                                       *
   !*****************************************************************************

  use par_mod
  use com_mod
  use windfields_mod

  implicit none

  integer, intent(in) :: nsp
  real, intent(in) :: xt, yt, zt
  real, intent(out) :: settling
  integer :: indz

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: rho1(2),tt1(2),temperature,airdens,vis_dyn,vis_kin
  real :: settling_old,reynolds,c_d
  integer :: i,n,nix,njy,indzh

  ! Variables needed for drag coefficient calculation
  real :: dfdr,kn,ks,alpha1,beta1,kn1

  !*****************************************************************************
  ! 1. Interpolate temperature and density: nearest neighbor interpolation sufficient
  !*****************************************************************************

  nix=int(xt)
  njy=int(yt)

  ! Determine the level below the current position for u,v
  !*******************************************************
  indz=nz-1
  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      exit
    endif
  end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz


  ! Bilinear horizontal interpolation
  !**********************************

  ! Loop over 2 levels
  !*******************

  do n=1,2
    indzh=indz+n-1
    rho1(n)=rho(nix,njy,indzh,1)
    tt1(n)=tt(nix,njy,indzh,1)
  end do


  ! Linear vertical interpolation
  !******************************

  temperature=dz2*tt1(1)+dz1*tt1(2)
  airdens=dz2*rho1(1)+dz1*rho1(2)

  vis_dyn=viscosity(temperature)
  vis_kin=vis_dyn/airdens

  reynolds=dquer(nsp)/1.e6*abs(vsetaver(nsp))/vis_kin

  ! Iteration to determine both Reynolds number and settling velocity
  !******************************************************************

  settling_old=vsetaver(nsp)    ! initialize iteration with Stokes' law to define settling velocity of a sphere, constant viscosity estimate

  if (ishape(nsp).eq.0) then
    do i=1,20    ! do a few iterations

      ! if (reynolds.lt.1.917) then
      !   c_d=24./reynolds
      ! else if (reynolds.lt.500.) then
      !   c_d=18.5/(reynolds**0.6)
      ! else
      !   c_d=0.44
      ! endif

      ! Clift and Guavin 1971 model

      ! call find_dragcoeff(reynolds, c_d)

      ! If reynolds<0.2, approximation 24/reynolds is valid
      ! 24/reynolds makes up >99% of all components
      !****************************************************
      if (reynolds.le.0.02) then
        c_d=(24.0/reynolds)

      else
        c_d=(24.0/reynolds)*(1+0.15*(reynolds**0.687))+ &
          0.42/(1.0+42500.0/(reynolds**1.16))
      endif

      settling=-1.* &
           sqrt(4*ga*dquer(nsp)/1.e6*density(nsp)*cunningham(nsp)/ &
           (3.*c_d*airdens))

      if (abs((settling-settling_old)/settling).lt.0.01) exit  ! stop iteration

      reynolds=dquer(nsp)/1.e6*abs(settling)/vis_kin
      settling_old=settling
    end do

  else ! Drag coefficient scheme by Bagheri & Bonadonna, 2016 to define settling velocities of other shapes (by D.Tatsii)

    ! Orientation of particles
    !*************************
    if (orient(nsp).eq.0) then
      ! Horizontal orientation
      ks=ks2(nsp)  ! B&B Figure 12 k_(s,max)
      kn=kn2(nsp)
    else if (orient(nsp).eq.1) then 
      ! Random orientation
      dfdr=density(nsp)/airdens

      alpha1=0.45+10.0/(exp(2.5*log10(dfdr))+30.0)
      beta1=1.-37.0/(exp(3.0*log10(dfdr))+100.0)
      ks=ks1(nsp)
      kn=10.**(alpha1*(-log10(Fn(nsp)))**beta1)
    else
      ! The average of random and horizontal orientation
      dfdr=density(nsp)/airdens

      alpha1=0.45+10.0/(exp(2.5*log10(dfdr))+30.0)
      beta1=1.-37.0/(exp(3.0*log10(dfdr))+100.0)
      kn1=10.**(alpha1*(-log10(Fn(nsp)))**beta1)
      ks=(ks1(nsp)+ks2(nsp))*0.5
      kn=(kn1+kn2(nsp))*0.5
    endif

    do i=1,20
      c_d=(24.*ks/reynolds)*(1.+0.125*((reynolds*kn/ks)**(2./3.)))+ &
        (0.46*kn/(1.+5330./(reynolds*kn/ks)))

      settling=-1.* &
              sqrt(4.*ga*dquer(nsp)/1.e6*density(nsp)*cunningham(nsp)/ &
              (3.*c_d*airdens))

      if (abs((settling-settling_old)/settling).lt.0.01) exit

      reynolds=dquer(nsp)/1.e6*abs(settling)/vis_kin
      settling_old=settling
    end do
  endif
end subroutine get_settling

real function viscosity(t)
  ! Function calculates dynamic viscosity of air (kg/m/s) as function of
  ! temperature (K) using Sutherland's formula
  implicit none

  real :: t
  real,parameter :: c=120.,t_0=291.15,eta_0=1.827e-5

  viscosity=eta_0*(t_0+c)/(t+c)*(t/t_0)**1.5

  return

end function viscosity

end module settling_mod
