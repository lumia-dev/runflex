! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later
module cbl_mod
    
    implicit none

    private :: cuberoot

    public :: cbl,reinit_particle,init_cbl_vel
    
contains

subroutine cbl(wp,zp,wst,h,rhoa,rhograd,sigmaw,dsigmawdz,tlw,ptot,Q,phi,ath,bth,ol,flagrein)
!                i  i i  i   i  i    i     i       i         i   o   o   o   o    o  i    o
     
    use par_mod, only:pi
    use com_mod, only:ldirect
    
    implicit none

!*******************************************************************************
! CBL skewed vertical profiles and formulation of LHH 1996 with profile of w^3
! from LHB 2000 
! LHH formulation has been modified to account for variable density profiles 
! and backward in time or forward in time simulations 
! see Cassiani et al. BLM 2014 doi  for explanations and references  
!*******************************************************************************

    real :: usurad2,usurad2p,C0,costluar4,eps 
    parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=3,costluar4=0.66667,eps=0.000001)

       
    integer :: flagrein
    real :: wp,zp,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,rhoa,rhograd
    real ::fluarw,fluarw2
    real ::w3,w2
    real ::dw3,dw2
    real ::wb,wa
    real ::deltawa,deltawb
    real ::wold,wold2
    real ::pa,pb,alfa
    real ::Phi,Q,ptot  
    real :: timedir
    real ::ol,transition   
    

    real :: &
    erf, &
    aperfa, &
    aperfb, &
    ath, &
    bth

    real ::  &
    z, &    
    skew, &
    skew2, &
    radw2, &
    rluarw, &
    xluarw, &
    aluarw, &
    bluarw, &
    sigmawa, &
    sigmawb, &
    dskew, &
    dradw2, &
    dfluarw, &
    drluarw, &
    dxluarw, &
    daluarw, &
    dbluarw, &
    dsigmawa, &
    dsigmawb, &
    dwa, &
    dwb, &
    sigmawa2, &
    sigmawb2
    
    
    dens=rhoa      !put to 1 for test constant density simulation
    ddens=rhograd  !put to 0 for test constant density simulation
             
    
    timedir=ldirect !ldirect contains direction of time forward (1) or backward(-1)
    ! assegnazione z 
    z=(zp/h)
    
    ! stability transition function see Cassiani et al(2015) BLM
    transition=1.
    !if (ol.lt.-50) transition=((sin(((ol+100.)/100.)*pi))-1.)/2.
    if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))*0.5+0.5

    ! momento secondo 
    
    w2=(sigmaw*sigmaw)
    dw2=(2.*sigmaw*dsigmawdz)
    

    !=================== dissipazione =======================================
    
    alfa=2.*w2/(C0*tlw)

    !========================================================================
    
    wold=timedir*wp
    
    !=========================================================================
    !------------------------------ momento terzo ============================
                                
    w3=((1.2*z*((1.-z)**(3./2.)))+eps)*(wst**3)*transition
    dw3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))*(wst**3)*(1./h)*transition

!===========================================================================0
   
    skew=w3/(w2**1.5)
    skew2=skew*skew
    dskew=(dw3*w2**(1.5)-w3*1.5*w2**0.5*dw2)/w2**3
    radw2=w2**0.5
    dradw2=0.5*w2**(-0.5)*dw2
    !costluar4=0.66667  ! costante da LHH
    fluarw=costluar4*(cuberoot(skew))   !skew**(1./3.)
    fluarw2=fluarw*fluarw
    
    if (skew.ne.0) then
      
        dfluarw=costluar4*(1./3.)*cuberoot(skew**(-2.))*dskew
    
        rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
        xluarw=(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)    !----> r^1/2
        
        drluarw=( ((3.*(1.+fluarw2)**2*(2.*fluarw*dfluarw)*skew2)+ &
            (1.+fluarw2)**3*2.*skew*dskew) *(3.+fluarw2)**2.*fluarw2 - &
            (1.+fluarw2)**3*skew2* &
            ( (2.*(3.+fluarw2)*(2.*fluarw*dfluarw)*fluarw2) + &
            (3.+fluarw2)**2*2.*fluarw*dfluarw) )/ &
            (((3.+fluarw2)**2.*fluarw2)**2)
                                     
        dxluarw=( ((1.5*(1.+fluarw2)**0.5*(2.*fluarw*dfluarw)*skew)+ &
            (1.+fluarw2)**1.5*dskew) *(3.+fluarw2)*fluarw - &
            (1.+fluarw2)**1.5*skew* &
            (3.*dfluarw+3*fluarw2*dfluarw) )/ &
            (((3.+fluarw2)*fluarw)**2)
    
    else        
        dfluarw=0.
        rluarw=0.
        drluarw=0.
        xluarw=0.
        dxluarw=0.
    end if



   aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
   bluarw=1.-aluarw

   daluarw=-0.5*(  (dxluarw*(4.+rluarw)**0.5) - &
           (0.5*xluarw*(4.+rluarw)**(-0.5)*drluarw) ) &
           /(4.+rluarw) 
   dbluarw=-daluarw
   
   sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
   sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

   dsigmawa=dradw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5+ &
            radw2*( &
            (0.5*(bluarw/(aluarw*(1.+fluarw2)))**(-0.5))      * &
            ( &
            (dbluarw*(aluarw*(1.+fluarw2))- &
            bluarw*(daluarw*(1.+fluarw2)+aluarw*2.*fluarw*dfluarw)) &
            /((aluarw*(1.+fluarw2))**2) &
            ) &
            )
  dsigmawb=dradw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5+ &
           radw2*( &
           (0.5*(aluarw/(bluarw*(1.+fluarw2)))**(-0.5))      * &
           ( &
           (daluarw*(bluarw*(1.+fluarw2))- &
           aluarw*(dbluarw*(1.+fluarw2)+bluarw*2.*fluarw*dfluarw)) &
           /((bluarw*(1.+fluarw2))**2) &
           ) &
           )
    
           wa=(fluarw*sigmawa)
           wb=(fluarw*sigmawb)
           dwa=dfluarw*sigmawa+fluarw*dsigmawa
           dwb=dfluarw*sigmawb+fluarw*dsigmawb

           deltawa=wold-wa
           deltawb=wold+wb
           wold2=wold*wold
           sigmawa2=sigmawa*sigmawa
           sigmawb2=sigmawb*sigmawb
                            
           if (abs(deltawa).gt.6.*sigmawa.and.abs(deltawb).gt.6.*sigmawb) flagrein=1

           pa=(usurad2p*(1./sigmawa))*(exp(-(0.5*((deltawa/sigmawa)**2.))))
           pb=(usurad2p*(1./sigmawb))*(exp(-(0.5*((deltawb/sigmawb)**2.))))
                                    
           ptot=dens*aluarw*pa+dens*bluarw*pb
                               
           aperfa=deltawa*usurad2/sigmawa
           aperfb=deltawb*usurad2/sigmawb

           Phi=-0.5* &
               (aluarw*dens*dwa+dens*wa*daluarw+aluarw*wa*ddens)*erf(aperfa) &
               +sigmawa*(aluarw*dens*dsigmawa*(wold2/sigmawa2+1.)+ &
               sigmawa*dens*daluarw+sigmawa*ddens*aluarw+ &
               aluarw*wold*dens/sigmawa2*(sigmawa*dwa-wa*dsigmawa))*pa &
               +0.5* &
               (bluarw*dens*dwb+wb*dens*dbluarw+wb*bluarw*ddens)*erf(aperfb) &
               +sigmawb*(bluarw*dens*dsigmawb*(wold2/sigmawb2+1.)+ &
               sigmawb*dens*dbluarw+sigmawb*ddens*bluarw+ &
               bluarw*wold*dens/sigmawb2*(-sigmawb*dwb+wb*dsigmawb))*pb

               Q=timedir*((aluarw*dens*deltawa/sigmawa2)*pa+ &
               (bluarw*dens*deltawb/sigmawb2)*pb)

               ath=(1./ptot)*(-(C0/2.)*alfa*Q+phi)
               bth=sqrt(C0*alfa)
              !bth=sngl(sigmaw*sqrt(2.*tlw))

    return
end subroutine cbl

function cuberoot(x) result(y)
    
    implicit none

    real, intent(in) :: x
    real :: y 
    real, parameter :: third=0.333333333

    y=sign((abs(x))**third,x)
end function cuberoot

subroutine reinit_particle(zp,wst,h,sigmaw,wp,nrand,ol)
!                                      i   i  i   i  i    io  io    i 

!******************************************************************************
! CBL skewed vertical profiles and formulation of LHH 1996 with profile of w^3
! from lHB 2000 
! LHH formulation has been modified to account for variable density profiles 
! and backward in time or forward in time simulations
! This routine re-initialize particle velocity if a numerical instability 
! in the cbl scheme generated a NaN value 
! The particle velocity is extracted from the updraft and downdraft 
! distribution as required  
! The re-initialization si not perfect 
! See e.g. Cassiani et al(2015) BLM           
!******************************************************************************
  use par_mod, only:pi
  use com_mod, only:ldirect,rannumb

  implicit none


  real :: usurad2,usurad2p,C0,costluar4,eps 
  parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

  integer nrand
  real :: wp,zp,wst,h,sigmaw,dcas1 !,ran3,gasdev
  real :: w3,w2
  real ::  z, &    
       skew, &
       skew2, &
       radw2, &
       fluarw,fluarw2, &
       rluarw, &
       xluarw, &
       aluarw, &
       bluarw, &
       sigmawa, &
       sigmawb, &  
       wb,wa 
  real timedir
  real ol,transition

!---------------------------------------------------------------------------  
!timedir direction of time forward (1) or backward(-1)
  nrand=nrand+1
  dcas1=rannumb(nrand)
  timedir=ldirect 
  z=zp/h
  transition=1.

  if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))*0.5+0.5

  w2=sigmaw*sigmaw  
  w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3)*transition 
  skew=w3/(w2**1.5)
  skew2=skew*skew
  radw2=sqrt(w2) !sigmaw

  fluarw=costluar4*skew**0.333333333333333
  fluarw2=fluarw*fluarw
  rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
  xluarw=rluarw**0.5 !(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)    !----> r^1/2

  aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
  bluarw=1.-aluarw

  sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
  sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

  wa=(fluarw*sigmawa)
  wb=(fluarw*sigmawb)



  if ((sign(1.,wp)*timedir).gt.0) then !updraft   
100 wp=(dcas1*sigmawa+wa)
    if (wp.lt.0)  then
      nrand=nrand+1
      dcas1=rannumb(nrand)
      goto 100
    end if
    wp=wp*timedir
  else if ((sign(1.,wp)*timedir).lt.0) then !downdraft    
101 wp=(dcas1*sigmawb-wb)
    if (wp.gt.0)  then 
      nrand=nrand+1
      dcas1=rannumb(nrand)
      goto 101
    end if
    wp=wp*timedir
  end if

  return
end subroutine reinit_particle

subroutine init_cbl_vel(idum,zp,wst,h,sigmaw,wp,ol,ithread)
  !                              i/o   i  i   i  i     i  o   i  

  use par_mod, only:pi
  use com_mod, only:ldirect
  use random_mod, only: gasdev, ran3

  implicit none
  !===============================================================================
  ! CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3
  ! from LHB 2000
  ! LHH formulation has been modified to account for variable density profiles and
  ! backward in time or forward in time simulations
  ! see Cassiani et al. BLM 2014 doi  for explanations and references
  !===============================================================================

  integer,intent(in) :: ithread
  real :: usurad2,usurad2p,C0,costluar4,eps 
  parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

  
  real :: wp,zp,wst,h,sigmaw,dcas,dcas1!,ran3,gasdev
  real :: w3,w2
  real ::  z, &    
       skew, &
       skew2, &
       radw2, &
       fluarw,fluarw2, &
       rluarw, &
       xluarw, &
       aluarw, &
       bluarw, &
       sigmawa, &
       sigmawb, &  
       wb,wa 
  real timedir
  real ol, transition
  integer :: idum
  
  !---------------------------------------------------------------------------  
  timedir=ldirect !direction of time forward (1) or backward(-1)
  z=zp/h


  transition=1.
  if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))*0.5+0.5  !see also in cbl.f90

  w2=sigmaw*sigmaw
  w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3) *transition

  skew=w3/(w2**1.5)
  skew2=skew*skew

  radw2=sqrt(w2) !sigmaw

  fluarw=costluar4*skew**0.333333333333333
  fluarw2=fluarw*fluarw
  rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
  xluarw=rluarw**0.5 !----> r^1/2

  aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
  bluarw=1.-aluarw

  sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
  sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

  wa=(fluarw*sigmawa)
  wb=(fluarw*sigmawb)

  dcas=ran3(idum,ithread) 

  if (dcas.le.aluarw) then
    dcas1=gasdev(idum,ithread)
    wp=timedir*(dcas1*sigmawa+wa)
  else
    dcas1=gasdev(idum,ithread)
    wp=timedir*(dcas1*sigmawb-wb)
  end if

  return
end subroutine init_cbl_vel

end module cbl_mod
