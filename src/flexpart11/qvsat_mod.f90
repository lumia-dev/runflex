! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!##################################################################
!##################################################################
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
module qvsat_mod
  
  implicit none
  private

  public :: f_qvsat, ew

contains

real function f_qvsat( p, t )

  !PURPOSE:
  !
  !Calculate the saturation specific humidity using enhanced Teten's
  !formula.
  !
  !AUTHOR: Yuhe Liu
  !01/08/1998
  !
  !MODIFICATION HISTORY:
  !
  !INPUT :
  !  p        Pressure (Pascal)
  !  t        Temperature (K)
  !OUTPUT:
  !  f_qvsat  Saturation water vapor specific humidity (kg/kg).
  !
  !Variable Declarations.
  !

  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: fespt

  real,parameter ::  rd     = 287.0 ! Gas constant for dry air  (m**2/(s**2*K))
  real,parameter ::  rv     = 461.0 ! Gas constant for water vapor  (m**2/(s**2*K)).
  real,parameter ::  rddrv  = rd/rv


  ! Change by A. Stohl to save computation time:
  ! IF ( t.ge.273.15 ) THEN     ! for water
  if ( t.ge.253.15 ) then      ! modification Petra Seibert
                               ! (supercooled water may be present)
    fespt=f_esl(p,t)
  else
    fespt=f_esi(p,t)
  endif

!!$  f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)      !old
  if (p-(1.0-rddrv)*fespt == 0.) then                     !bugfix
     f_qvsat = 1.
  else
     f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)
  end if

  return
end function f_qvsat


real function f_esl( p, t )
  ! Saturation water vapor pressure over liquid water
  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: f

  !#######################################################################
  !
  !Saturation specific humidity parameters used in enhanced Teten's
  !formula. (See A. Buck, JAM 1981)
  !
  !#######################################################################

  real,parameter ::  satfwa = 1.0007
  real,parameter ::  satfwb = 3.46e-8  ! for p in Pa

  real,parameter ::  satewa = 611.21   ! es in Pa
  real,parameter ::  satewb = 17.502
  real,parameter ::  satewc = 32.18

  real,parameter ::  satfia = 1.0003
  real,parameter ::  satfib = 4.18e-8  ! for p in Pa

  real,parameter ::  sateia = 611.15   ! es in Pa
  real,parameter ::  sateib = 22.452
  real,parameter ::  sateic = 0.6

  f = satfwa + satfwb * p
  f_esl = f * satewa * exp( satewb*(t-273.15)/(t-satewc) )

  return
end function f_esl

real function f_esi( p, t )
  ! Saturation water vapor pressure over ice (Pa)
  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: f

  !#######################################################################
  !
  !Saturation specific humidity parameters used in enhanced Teten's
  !formula. (See A. Buck, JAM 1981)
  !
  !#######################################################################
  !
  real,parameter ::  satfwa = 1.0007
  real,parameter ::  satfwb = 3.46e-8  ! for p in Pa

  real,parameter ::  satewa = 611.21   ! es in Pa
  real,parameter ::  satewb = 17.502
  real,parameter ::  satewc = 32.18

  real,parameter ::  satfia = 1.0003
  real,parameter ::  satfib = 4.18e-8  ! for p in Pa

  real,parameter ::  sateia = 611.15   ! es in Pa
  real,parameter ::  sateib = 22.452
  real,parameter ::  sateic = 0.6

  f = satfia + satfib * p
  f_esi = f * sateia * exp( sateib*(t-273.15)/(t-sateic) )

  return
end function f_esi

real function ew(x,p) ! p is not used

  !****************************************************************
  !SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
  !NACH DER GOFF-GRATCH-FORMEL.
  !****************************************************************

  implicit none

  real :: x, y, a, p , c, d

  ew=0.
  if(x.le.0.) error stop 'sorry: t not in [k]'
  ! Formula of Goff and Gratch (after Murray, 1966)
  ! if (x.lt.273.15) then
  ! ! Above ice
  !   a = 273.15/x
  !   y = -20.947031*a - 3.56654*log(a) - 2.01889049/a
  !   ew = 5.75185606E10*exp(y)
  ! else
  ! ! Above water
  !   a = 373.15/x 
  !   y = -18.1972839*a + 5.02808*log(a) - 70242.1852*exp(-26.1205253/a) + &
  !     58.0691913*exp(-8.03945282*a)
  !   ew = 7.95357242E10*exp(y)
  ! endif

  ! ! Formula of Magnus (after Murray, 1966)
  ! if (x.lt.273.15) then
  ! ! Above ice
  !   ew = 6.1078*exp(21.8745584*(x-273.15)/(x-7.66))
  ! else 
  ! ! Above water
  !   ew = 6.1078*exp(17.2693882*(x-273.15)/(x-35.86))
  ! endif

  ! Formula of Buck 1981
  ! ew = f_qvsat(p,x)

  ! ! Original
  y=373.16/x
  a=-7.90298*(y-1.)
  a=a+(5.02808*0.43429*alog(y))
  c=(1.-(1./y))*11.344
  c=-1.+(10.**c)
  c=-1.3816*c/(10.**7)
  d=(1.-y)*3.49149
  d=-1.+(10.**d)
  d=8.1328*d/(10.**3)
  y=a+c+d
  ew=101324.6*(10.**y)       ! Saettigungsdampfdruck in Pa

end function ew

end module qvsat_mod