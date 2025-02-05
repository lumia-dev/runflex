module pbl_profile_mod
  use par_mod
  use qvsat_mod

  implicit none

contains

function psih (z,l)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the stability correction term                           *
  !                                                                            *
  !     AUTHOR: Matthias Langer, adapted by Andreas Stohl (6 August 1993)      *
  !             Update: G. Wotawa, 11 October 1994                             *
  !                                                                            *
  !     Literature:                                                            *
  !     [1] C.A.Paulson (1970), A Mathematical Representation of Wind Speed    *
  !           and Temperature Profiles in the Unstable Atmospheric Surface     *
  !           Layer. J.Appl.Met.,Vol.9.(1970), pp.857-861.                     *
  !                                                                            *
  !     [2] A.C.M. Beljaars, A.A.M. Holtslag (1991), Flux Parameterization over*
  !           Land Surfaces for Atmospheric Models. J.Appl.Met. Vol. 30,pp 327-*
  !           341                                                              *
  !                                                                            *
  !     Variables:                                                             *
  !     L     = Monin-Obukhov-length [m]                                       *
  !     z     = height [m]                                                     *
  !     zeta  = auxiliary variable                                             *
  !                                                                            *
  !     Constants:                                                             *
  !     eps   = 1.2E-38, SUN-underflow: to avoid division by zero errors       *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real :: psih,x,z,zeta,l
  real,parameter :: a=1.,b=0.667,c=5.,d=0.35,eps=1.e-20

  if ((l.ge.0).and.(l.lt.eps)) then
    l=eps
  else if ((l.lt.0).and.(l.gt.(-1.*eps))) then
    l=-1.*eps
  endif

  if ((log10(z)-log10(abs(l))).lt.log10(eps)) then
    psih=0.
  else
    zeta=z/l
    if (zeta.gt.0.) then
      psih = - (1.+0.667*a*zeta)**(1.5) - b*(zeta-c/d)*exp(-d*zeta) &
           - b*c/d + 1.
    else
      x=(1.-16.*zeta)**(.25)
      psih=2.*log((1.+x*x)*0.5)
    end if
  end if

end function psih

real function psim(z,al)

  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION: CALCULATION OF THE STABILITY CORRECTION FUNCTION FOR   *
  !              MOMENTUM AS FUNCTION OF HEIGHT Z AND OBUKHOV SCALE     *
  !              HEIGHT L                                               *
  !                                                                     *
  !**********************************************************************

  implicit none

  real :: z,al,zeta,x,a1,a2

  zeta=z/al
  if(zeta.le.0.) then
  ! UNSTABLE CASE
    x=(1.-15.*zeta)**0.25
    a1=((1.+x)*0.5)**2
    a2=(1.+x**2)*0.5
    psim=log(a1*a2)-2.*atan(x)+pi*0.5
  else
  ! STABLE CASE
    psim=-4.7*zeta
  endif

end function psim

subroutine pbl_profile(ps,td2m,zml1,t2m,tml1,u10m,uml1,stress,hf)

  !********************************************************************
  !                                                                   *
  !                    G. WOTAWA, 1995-07-07                          *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! DESCRIPTION: CALCULATION OF FRICTION VELOCITY AND SURFACE SENS-   *
  !              IBLE HEAT FLUX USING THE PROFILE METHOD (BERKOVICZ   *
  !              AND PRAHM, 1982)                                     *
  !                                                                   *
  ! Output now is surface stress instead of ustar                     *
  !                                                                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! INPUT:                                                            *
  !                                                                   *
  !                                                                   *
  ! ps      surface pressure(Pa)                                      *
  ! td2m    two metre dew point(K)                                    *
  ! zml1    heigth of first model level (m)                           *
  ! t2m     two metre temperature (K)                                 *
  ! tml1    temperature first model level (K)                         *
  ! u10m    ten metre wind speed (ms-1)                               *
  ! uml1    wind speed first model level (ms-1)                       *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! OUTPUT:                                                           *
  !                                                                   *
  ! stress  surface stress (i.e., friction velocity (ms-1) squared    *
  !                         multiplied with air density)              *
  ! hf      surface sensible heat flux (Wm-2)                         *
  !                                                                   *
  !********************************************************************
  ! ustar   friction velocity (ms-1)                                  *
  ! maxiter maximum number of iterations                              *
  !********************************************************************

  implicit none

  integer :: iter
  real :: ps,td2m,rhoa,zml1,t2m,tml1,u10m,uml1,ustar,hf
  real :: al,alold,aldiff,tmean,crit
  real :: deltau,deltat,thetastar,e,tv,stress
  integer,parameter :: maxiter=10
  real,parameter    :: r1=0.74

  e=ew(td2m,ps)               ! vapor pressure
  tv=t2m*(1.+0.378*e/ps)   ! virtual temperature
  rhoa=ps/(r_air*tv)       ! air density

  deltau=uml1-u10m         !! Wind Speed difference between
                           !! Model level 1 and 10 m

  if(deltau.le.0.001) then    !! Monin-Obukhov Theory not
    al=9999.               !! applicable --> Set dummy values
    ustar=0.01
    stress=ustar*ustar*rhoa
    hf=0.0
    return
  endif
  deltat=tml1-t2m+0.0098*(zml1-2.)  !! Potential temperature difference
                                    !! between model level 1 and 10 m

  if(abs(deltat).le.0.03) then    !! Neutral conditions
    hf=0.0
    al=9999.
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    stress=ustar*ustar*rhoa
    return
  endif

  tmean=0.5*(t2m+tml1)
  crit=(0.0219*tmean*(zml1-2.0)*deltau**2)/ &
       (deltat*(zml1-10.0)**2)
  if((deltat.gt.0).and.(crit.le.1.)) then
                                    !! Successive approximation will
    al=50.                          !! not converge
    ustar=(vonkarman*deltau)/ &
         (log(zml1*0.1)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1*0.5)-psih(zml1,al)+psih(2.,al))
    hf=rhoa*cpa*ustar*thetastar
    stress=ustar*ustar*rhoa
    return
  endif

  al=9999.                 ! Start iteration assuming neutral conditions
  do iter=1,maxiter
    alold=al
    ustar=(vonkarman*deltau)/ &
         (log(zml1*0.1)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1*0.5)-psih(zml1,al)+psih(2.,al))
    al=(tmean*ustar**2)/(ga*vonkarman*thetastar)
    aldiff=abs((al-alold)/alold)
    if(aldiff.lt.0.01) exit  !! Successive approximation successful
  end do
  hf=rhoa*cpa*ustar*thetastar
  if(al.gt.9999.) al=9999.
  if(al.lt.-9999.) al=-9999.

  stress=ustar*ustar*rhoa
end subroutine pbl_profile

end module pbl_profile_mod