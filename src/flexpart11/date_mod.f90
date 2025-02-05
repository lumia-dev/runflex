! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
!   L. Bakels 2022: this module contains all subroutines related to          *
!                   calculations between dates: caldate and juldate          *
!                                                                            *
!*****************************************************************************

module date_mod
  use par_mod, only: dp
  
  implicit none

contains

subroutine caldate(juliandate,yyyymmdd,hhmiss)
  !                      i       o       o
  !*****************************************************************************
  !                                                                            *
  !     Calculates the Gregorian date from the Julian date                     *
  !                                                                            *
  !     AUTHOR: Andreas Stohl (21 January 1994), adapted from Numerical Recipes*
  !                                                                            *
  !     PS 2020-07-27: add a check to avoid giving back 240000 for hhmiss      *
  !                                                                            *
  !                                                                            *
  !     Variables:                                                             *
  !     dd             Day                                                     *
  !     hh             Hour                                                    *
  !     hhmiss         Hour, Minute, Second                                    *
  !     ja,jb,jc,jd,je help variables                                          *
  !     jalpha         help variable                                           *
  !     juliandate        Julian Date                                          *
  !     julday         help variable                                           *
  !     mi             Minute                                                  *
  !     mm             Month                                                   *
  !     ss             Seconds                                                 *
  !     yyyy           Year                                                    *
  !     yyyymmdd       Year, Month, Day                                        *
  !                                                                            *
  !     Constants:                                                             *
  !     igreg          help constant                                           *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer           :: yyyymmdd,yyyy,mm,dd,hhmiss,hh,mi,ss
  integer           :: julday,ja,jb,jc,jd,je,jalpha
  real(kind=dp)     :: juliandate
  integer,parameter :: igreg=2299161

  julday=int(juliandate)
  ! PS check to avoid 240000 as hhmiss:  
  if ((juliandate-julday)*86400._dp .ge. 86399.5_dp) then
    juliandate = juliandate + juliandate-julday-86399.5_dp/86400._dp
    julday=int(juliandate)
  endif
  if(julday.ge.igreg)then
    jalpha=int(((julday-1867216)-0.25)/36524.25)
    ja=julday+1+jalpha-int(0.25*jalpha)
  else
    ja=julday
  endif
  jb=ja+1524
  jc=int(6680.+((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)
  dd=jb-jd-int(30.6001*je)
  mm=je-1
  if (mm.gt.12) mm=mm-12
  yyyy=jc-4715
  if (mm.gt.2) yyyy=yyyy-1
  if (yyyy.le.0) yyyy=yyyy-1

  yyyymmdd=10000*yyyy+100*mm+dd
  hh=int(24._dp*(juliandate-real(julday,kind=dp)))
  mi=int(1440._dp*(juliandate-real(julday,kind=dp))-60._dp*real(hh,kind=dp))
  ss=nint(86400._dp*(juliandate-real(julday,kind=dp))-3600._dp*real(hh,kind=dp)- &
       60._dp*real(mi,kind=dp))
  if (ss.eq.60) then  ! 60 seconds = 1 minute
    ss=0
    mi=mi+1
  endif
  if (mi.eq.60) then
    mi=0
    hh=hh+1
  endif
  hhmiss=10000*hh+100*mi+ss

end subroutine caldate

real(kind=dp) function juldate(yyyymmdd,hhmiss)

  !*****************************************************************************
  !                                                                            *
  !     Calculates the Julian date                                             *
  !                                                                            *
  !     AUTHOR: Andreas Stohl (15 October 1993)                                *
  !                                                                            *
  !     Variables:                                                             *
  !     dd             Day                                                     *
  !     hh             Hour                                                    *
  !     hhmiss         Hour, minute + second                                   *
  !     ja,jm,jy       help variables                                          *
  !     juldate        Julian Date                                             *
  !     julday         help variable                                           *
  !     mi             Minute                                                  *
  !     mm             Month                                                   *
  !     ss             Second                                                  *
  !     yyyy           Year                                                    *
  !     yyyymmddhh     Date and Time                                           *
  !                                                                            *
  !     Constants:                                                             *
  !     igreg          help constant                                           *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer           :: yyyymmdd,yyyy,mm,dd,hh,mi,ss,hhmiss
  integer           :: julday,jy,jm,ja
  integer,parameter :: igreg=15+31*(10+12*1582)
  !real(kind=dp)     :: juldate

  yyyy=yyyymmdd/10000
  mm=(yyyymmdd-10000*yyyy)/100
  dd=yyyymmdd-10000*yyyy-100*mm
  hh=hhmiss/10000
  mi=(hhmiss-10000*hh)/100
  ss=hhmiss-10000*hh-100*mi

  if (yyyy.eq.0) then
     error stop 'juldate: there is no year zero'
  end if
  if (yyyy.lt.0) yyyy=yyyy+1
  if (mm.gt.2) then
    jy=yyyy
    jm=mm+1
  else
    jy=yyyy-1
    jm=mm+13
  endif
  julday=int(365.25*jy)+int(30.6001*jm)+dd+1720995
  if (dd+31*(mm+12*yyyy).ge.igreg) then
    ja=int(0.01*jy)
    julday=julday+2-ja+int(0.25*ja)
  endif

  juldate=real(julday,kind=dp)   + real(hh,kind=dp)/24._dp + &
       real(mi,kind=dp)/1440._dp  + real(ss,kind=dp)/86400._dp

end function juldate

  !*****************************************************************************
  !                                                                            * 
  !    Calculates number of days in a month                                    * 
  !                                                                            * 
  !    Author: Rona Thompson (Sep 2023)                                        * 
  !                                                                            * 
  !    Variables:                                                              * 
  !    yyyymm       year and month                                             * 
  !    eomday       number of days in month (end of month day)                 * 
  !                                                                            * 
  !*****************************************************************************

  integer function calceomday(yyyymm)

    integer, intent(in) :: yyyymm
    integer :: yyyy,mm
    integer, dimension(12) :: leapdays,days
    integer :: eomday

    leapdays=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    days=(/31,28,31,30,31,30,31,31,30,31,30,31/)

    yyyy=floor(yyyymm/100.)
    mm=yyyymm-yyyy*100

    if((float(yyyy)/100.).eq.float(yyyy/100)) then
      if((float(yyyy)/400.).eq.float(yyyy/400)) then
        eomday=leapdays(mm)
      else
        eomday=days(mm)
      endif
    else
      if((float(yyyy)/4.).eq.float(yyyy/4)) then
        eomday=leapdays(mm)
      else
        eomday=days(mm)
      endif
    endif

    calceomday=eomday

  end function calceomday


end module date_mod
