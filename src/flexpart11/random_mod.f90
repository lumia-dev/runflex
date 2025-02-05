! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!  Taken from Press et al., Numerical Recipes

module random_mod
  
  implicit none

  integer, parameter :: ran1_ntab=32
  integer, allocatable :: ran1_iv(:,:), ran1_iy(:)

  integer, allocatable :: gasdev_iset(:)
  real, allocatable :: gasdev_gset(:)

  integer, allocatable :: ran3_iff(:)
  integer, allocatable :: ran3_inext(:),ran3_inextp(:)
  integer, allocatable :: ma(:,:)

  integer, allocatable :: iseed1(:), iseed2(:)

contains
  
  subroutine alloc_random(num_threads)

    implicit none

    integer :: num_threads, i,stat

    allocate(ran1_iv(ran1_ntab,0:num_threads-1),ran1_iy(0:num_threads-1), stat=stat)
    if (stat.ne.0) error stop "Could not allocate ran1_iv"
    allocate(gasdev_iset(0:num_threads-1),gasdev_gset(0:num_threads-1), stat=stat)
    if (stat.ne.0) error stop "Could not allocate gasdev_iset"
    allocate(ran3_iff(0:num_threads-1),ran3_inext(0:num_threads-1),ran3_inextp(0:num_threads-1), stat=stat)
    if (stat.ne.0) error stop "Could not allocate ran3_iff"
    allocate(ma(55,0:num_threads-1), stat=stat)
    if (stat.ne.0) error stop "Could not allocate ma"
    allocate(iseed1(0:num_threads-1),iseed2(0:num_threads-1), stat=stat)
    if (stat.ne.0) error stop "Could not allocate iseed1"

    do i=0,num_threads-1
      iseed1(i) = -7-i
      iseed2(i) = -88-i
    end do
    ran3_iff(0:num_threads-1)=0
    ran1_iv(:,0:num_threads-1)=0
    ran1_iy(0:num_threads-1)=0
    gasdev_iset(0:num_threads-1)=0
    gasdev_gset(0:num_threads-1)=0
  end subroutine alloc_random

  subroutine dealloc_random()

    deallocate(ran1_iv,ran1_iy)
    deallocate(gasdev_iset,gasdev_gset)
    deallocate(ran3_iff,ran3_inext,ran3_inextp)
    deallocate(ma)
    deallocate(iseed1,iseed2)
  end subroutine dealloc_random

  function ran1(idum,ithread)

    implicit none

    integer :: idum,ithread
    real    :: ran1
    integer,parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
    integer,parameter :: ndiv=1+(im-1)/ran1_ntab
    real,parameter    :: am=1./im, eps=1.2e-7, rnmx=1.-eps
    integer :: j, k

    if (idum.le.0.or.ran1_iy(ithread).eq.0) then
      idum=max(-idum,1)
      do j=ran1_ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        if (j.le.ran1_ntab) ran1_iv(j,ithread)=idum
      enddo
      ran1_iy(ithread)=ran1_iv(1,ithread)
    endif
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum=idum+im
    j=1+ran1_iy(ithread)/ndiv
    ran1_iy(ithread)=ran1_iv(j,ithread)
    ran1_iv(j,ithread)=idum
    ran1=min(am*ran1_iy(ithread),rnmx)
  end function ran1


  function gasdev(idum,ithread)

    implicit none

    integer :: idum,ithread
    real    :: gasdev, fac, r, v1, v2

    if (gasdev_iset(ithread).eq.0) then
1     v1=2.*ran3(idum,ithread)-1.
      v2=2.*ran3(idum,ithread)-1.
      r=v1**2+v2**2
      if(r.ge.1.0 .or. r.eq.0.0) go to 1
      fac=sqrt(-2.*log(r)/r)
      gasdev_gset(ithread)=v1*fac
      gasdev=v2*fac
      gasdev_iset(ithread)=1
    else
      gasdev=gasdev_gset(ithread)
      gasdev_iset(ithread)=0
    endif
  end function gasdev


  subroutine gasdev1(idum,random1,random2)

    implicit none

    integer :: idum
    real :: random1, random2, fac, v1, v2, r

1   v1=2.*ran3(idum,0)-1.
    v2=2.*ran3(idum,0)-1.
    r=v1**2+v2**2
    if(r.ge.1.0 .or. r.eq.0.0) go to 1
    fac=sqrt(-2.*log(r)/r)
    random1=v1*fac
    random2=v2*fac
! Limit the random numbers to lie within the interval -3 and +3
!**************************************************************
    if (random1.lt.-3.) random1=-3.
    if (random2.lt.-3.) random2=-3.
    if (random1.gt.3.) random1=3.
    if (random2.gt.3.) random2=3.
  end subroutine gasdev1


  function ran3(idum,ithread)

    implicit none

    integer :: idum,ithread
    real :: ran3

    integer,parameter :: mbig=1000000000, mseed=161803398, mz=0
    real,parameter    :: fac=1./mbig
    integer :: i,ii,k
    integer :: mj,mk

    if(idum.lt.0 .or. ran3_iff(ithread).eq.0)then
      ran3_iff(ithread)=1
      mj=mseed-iabs(idum)
      mj=mod(mj,mbig)
      ma(55,ithread)=mj
      mk=1
      do i=1,54
        ii=mod(21*i,55)
        ma(ii,ithread)=mk
        mk=mj-mk
        if(mk.lt.mz)mk=mk+mbig
        mj=ma(ii,ithread)
      end do
      do k=1,4
        do i=1,55
          ma(i,ithread)=ma(i,ithread)-ma(1+mod(i+30,55),ithread)
          if(ma(i,ithread).lt.mz) ma(i,ithread)=ma(i,ithread)+mbig
        end do
      end do
      ran3_inext(ithread)=0
      ran3_inextp(ithread)=31
      idum=1
    endif
    ran3_inext(ithread)=ran3_inext(ithread)+1
    if(ran3_inext(ithread).eq.56) ran3_inext(ithread)=1
    ran3_inextp(ithread)=ran3_inextp(ithread)+1
    if(ran3_inextp(ithread).eq.56) ran3_inextp(ithread)=1
    mj=ma(ran3_inext(ithread),ithread)-ma(ran3_inextp(ithread),ithread)
    if(mj.lt.mz)mj=mj+mbig
    ma(ran3_inext(ithread),ithread)=mj
    ran3=mj*fac
  end function ran3
!  (C) Copr. 1986-92 Numerical Recipes Software US.

end module random_mod
