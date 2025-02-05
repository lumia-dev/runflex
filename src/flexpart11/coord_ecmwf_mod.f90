! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
!  This module handles conversions between ECMWF eta coordinates and         *
!  internal meter coordinates                                                *
!                                                                            *
!     Author: L. Bakels                                                      *
!*****************************************************************************

module coord_ecmwf_mod

  use par_mod
  use com_mod
  use windfields_mod

contains
#ifdef ETA
subroutine update_zeta_to_z(itime, ipart)
  use particle_mod
  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index

  if (.not. part(ipart)%alive) return  
  if (part(ipart)%etaupdate) return

  call zeta_to_z(itime,part(ipart)%xlon,part(ipart)%ylat, &
    part(ipart)%zeta,part(ipart)%z)
  part(ipart)%etaupdate = .true.
  part(ipart)%meterupdate = .true.
end subroutine update_zeta_to_z

subroutine update_z_to_zeta(itime, ipart)
  use particle_mod
  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index

  if (.not. part(ipart)%alive) return
  if (part(ipart)%meterupdate) return

  call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat, &
      part(ipart)%z,part(ipart)%zeta)
  part(ipart)%etaupdate = .true.
  part(ipart)%meterupdate = .true.
end subroutine update_z_to_zeta

subroutine z_to_zeta(itime,xt,yt,zold,zteta)
  !*****************************************************************************
  ! Converting z from meter coordinates to eta using logarithmic vertical      *
  ! interpolation                                                              *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zold,zold    spatial positions of trajectory (meters)                *
  ! zteta              vertical position in eta coordinates (output)           *
  !                                                                            *
  ! etauvheight defined in windfields: half model heights for ETA coordinates  *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,m,k,n                         ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zold                            ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    zteta                           ! converted output z in meters
  real ::                         &
    frac,                         & ! fraction between z levels
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
  real ::                         &
    prx,pr1,pr2     ! pressure of encompassing levels

  if (.not. log_interpol) then
    call z_to_zeta_lin(itime,xt,yt,zold,zteta)
    return
  endif

  call find_ngrid(xt,yt)
  call find_grid_indices(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_vars(itime)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! First estimate the level it is at, to reduce computation time
  n=nz-3
  if (ngrid.le.0) then
    do i=2,nz-1
      if ((etauvheight(ix,jy,i,memind(1)).gt.real(zold)) .or. &
          (etauvheight(ixp,jy,i,memind(1)).gt.real(zold)) .or. & 
          (etauvheight(ix,jyp,i,memind(1)).gt.real(zold)) .or. &
          (etauvheight(ixp,jyp,i,memind(1)).gt.real(zold))) then
        n=i-2
        exit
      endif
    end do
  else
    do i=2,nz-1
      if ((etauvheightn(ix,jy,i,memind(1),ngrid).gt.real(zold)) .or. &
          (etauvheightn(ixp,jy,i,memind(1),ngrid).gt.real(zold)) .or. & 
          (etauvheightn(ix,jyp,i,memind(1),ngrid).gt.real(zold)) .or. &
          (etauvheightn(ixp,jyp,i,memind(1),ngrid).gt.real(zold))) then
        n=i-2
        exit
      endif
    end do
  endif
  n=max(n,2)

  ztemp1 = 0.
  do i=n,nz-1
    k=i
    if (ngrid.le.0) then
      do m=1,2
        call hor_interpol(etauvheight,ttemp1(m),i,memind(m),nzmax)
      end do
    else
      do m=1,2
        call hor_interpol_nest(etauvheightn,ttemp1(m),i,memind(m),nzmax)
      end do
    endif
    call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

    if (ztemp2.gt.real(zold)) then
      !frac = (real(zold)-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ztemp1=ztemp2
  end do

  if (k.lt.nz-1) then
    if (ngrid.le.0) then
      do m=1,2
        call hor_interpol(ps,psint1(m),1,memind(m),1)
      end do
    else
      do m=1,2
        call hor_interpol_nest(psn,psint1(m),1,memind(m),1)
      end do
    endif
    call temporal_interpolation(psint1(1),psint1(2),psint)  
    pr1=akz(k-1) + bkz(k-1)*psint
    pr2=akz(k) + bkz(k)*psint

    prx=pr1/exp(log(pr2/pr1)/(ztemp2-ztemp1)*ztemp1) * &
      exp(log(pr2/pr1)/(ztemp2-ztemp1)*real(zold))
    frac=(prx-pr1)/(pr2 - pr1)
  endif

  zteta=real(uvheight(k-1)*(1.-frac)+uvheight(k)*frac,kind=dp)
end subroutine z_to_zeta

subroutine zeta_to_z(itime,xt,yt,zteta,ztout)
  !*****************************************************************************
  ! Converting z from eta coordinates to meters using logarithmic              *
  ! vertical interpolation                                                     *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zteta        spatial position of trajectory                          *
  ! ztout              vertical postion in meter (output)                      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,k,m,ii                        ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zteta                           ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    ztout                           ! converted output z in meters
  real(kind=dp) ::                &
    frac,                         &  ! fraction between z levels
    prx,pr1,pr2                      ! pressure of encompassing levels
  real ::                         &
    ztemp1(2),                    & ! z positions of the two encompassing levels
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
 
  if (.not. log_interpol) then
    call zeta_to_z_lin(itime,xt,yt,zteta,ztout)
    return
  endif

  ! Convert eta z coordinate to meters
  !***********************************
  call find_ngrid(xt,yt)
  call find_grid_indices(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_vars(itime)

  k=nz-1
  frac=1.
  do i=2,nz-1
    k=i
    if (zteta.ge.real(uvheight(k),kind=dp)) then
      frac=(zteta-real(uvheight(k-1),kind=dp))/(real(uvheight(k)-uvheight(k-1),kind=dp))
      exit
    endif
  end do

  if (ngrid.le.0) then
    do m=1,2
      call hor_interpol(ps,psint1(m),1,memind(m),1)
    end do
  else
    do m=1,2
      call hor_interpol_nest(psn,psint1(m),1,memind(m),1)
    end do
  endif

  call temporal_interpolation(psint1(1),psint1(2),psint)  
  pr1=real(akz(k-1) + bkz(k-1)*psint,kind=dp)
  pr2=real(akz(k) + bkz(k)*psint,kind=dp)
  prx=pr1*(1.-frac) + pr2*frac
  
  if (ngrid.le.0) then
    do ii=1,2
      do m=1,2
        call hor_interpol(etauvheight,ttemp1(m),k+ii-2,memind(m),nzmax)
      end do
      call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1(ii))
    end do
  else
    do ii=1,2
      do m=1,2
        call hor_interpol_nest(etauvheightn, &
          ttemp1(m),k+ii-2,memind(m),nzmax)
      end do
      call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1(ii))
    end do
  endif

  if ((pr2.eq.0).or.(pr1.eq.0)) then
    ztout = real(ztemp1(1),kind=dp)*(1.-frac)+real(ztemp1(2),kind=dp)*frac
    return
  endif
  
  ztout = ztemp1(1) + (ztemp1(2)-ztemp1(1))/log(pr2/pr1)*log(prx/pr1)
end subroutine zeta_to_z

subroutine w_to_weta(itime,dt,xt,yt,z_old,zeta_old,w_in,weta_out)
  !*****************************************************************************
  ! Converting z from meter coordinates to eta using logarithmic vertical      *
  ! interpolation                                                              *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zold,zold    spatial positions of trajectory (meters)                *
  ! zteta              vertical position in eta coordinates (output)           *
  !                                                                            *
  ! etauvheight defined in windfields: half model heights for ETA coordinates  *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in)          ::  &
    itime                            ! time index
  real, intent(in)             ::  &
    dt                               ! time step
  real(kind=dp), intent(in)    ::  &
    xt,yt,z_old,zeta_old              ! particle position
  real, intent(in)             ::  &
    w_in                             ! w in meters/s
  real, intent(inout)          ::  &
    weta_out                         ! converted output w in meters to eta
  real(kind=dp)                ::  &
    znew

  call z_to_zeta(itime,xt,yt,z_old+real(w_in*dt,kind=dp),znew)

  weta_out=real(znew-zeta_old)/dt

end subroutine w_to_weta

subroutine z_to_zeta_lin(itime,xt,yt,zold,zteta)
  !*****************************************************************************
  ! Converting z from meter coordinates to eta using linear interpolation      *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zold,zold    spatial positions of trajectory (meters)                *
  ! zteta              vertical position in eta coordinates (output)           *
  !                                                                            *
  ! etauvheight defined in windfields: half model heights for ETA coordinates  *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,m,k,n                         ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zold                            ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    zteta                           ! converted output z in meters
  real ::                         &
    frac,                         & ! fraction between z levels
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp1(2)                       ! storing virtual temperature

  call find_ngrid(xt,yt)
  call find_grid_indices(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_vars(itime)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! First estimate the level it is at, to reduce computation time
  n=nz-3
  if (ngrid.le.0) then
    do i=2,nz-1
      if ((etauvheight(ix,jy,i,memind(1)).gt.real(zold)) .or. &
          (etauvheight(ixp,jy,i,memind(1)).gt.real(zold)) .or. &
          (etauvheight(ix,jyp,i,memind(1)).gt.real(zold)) .or. &
          (etauvheight(ixp,jyp,i,memind(1)).gt.real(zold))) then
        n=i-2
        exit
      endif
    end do
  else
    do i=2,nz-1
      if ((etauvheightn(ix,jy,i,memind(1),ngrid).gt.real(zold)) .or. &
          (etauvheightn(ixp,jy,i,memind(1),ngrid).gt.real(zold)) .or. &
          (etauvheightn(ix,jyp,i,memind(1),ngrid).gt.real(zold)) .or. &
          (etauvheightn(ixp,jyp,i,memind(1),ngrid).gt.real(zold))) then
        n=i-2
        exit
      endif
    end do
  endif
  n=max(n,2)

  ztemp1 = 0.
  frac = 0.
  k=n
  do i=n,nz-1
    k=i
    if (ngrid.le.0) then
      do m=1,2
        call hor_interpol(etauvheight,ttemp1(m),i,memind(m),nzmax)
      end do
    else
      do m=1,2
        call hor_interpol_nest(etauvheightn,ttemp1(m),i,memind(m),nzmax)
      end do
    endif    
    call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

    if (ztemp2.gt.real(zold)) then
      frac = (real(zold)-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ztemp1=ztemp2
  end do

  zteta=real(uvheight(k-1)*(1.-frac)+uvheight(k)*frac,kind=dp)
end subroutine z_to_zeta_lin

subroutine zeta_to_z_lin(itime,xt,yt,zteta,ztout)

  !*****************************************************************************
  ! Converting z from eta coordinates to meters using linear interpolation     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zteta        spatial position of trajectory                          *
  ! ztout              vertical postion in meter (output)                      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,k,m,ii                        ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zteta                           ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    ztout                           ! converted output z in meters
  real(kind=dp) ::                &
    frac                            ! fraction between z levels
  real ::                         &
    ztemp1(2),                    & ! z positions of the two encompassing levels
    ttemp1(2)                       ! storing virtual temperature


  ! Convert eta z coordinate to meters
  !***********************************
  call find_ngrid(xt,yt)
  call find_grid_indices(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_vars(itime)

  k=nz-1
  frac=1.
  do i=2,nz-1
    k=i
    if (zteta.ge.real(uvheight(k),kind=dp)) then
      frac=(zteta-real(uvheight(k-1),kind=dp))/(real(uvheight(k)-uvheight(k-1),kind=dp))
      exit
    endif
  end do
  
  if (ngrid.le.0) then
    do ii=1,2
      do m=1,2
        call hor_interpol(etauvheight,ttemp1(m),k+ii-2,memind(m),nzmax)
      end do
      call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1(ii))
    end do
  else
    do ii=1,2
      do m=1,2
        call hor_interpol_nest(etauvheightn,ttemp1(m),k+ii-2,memind(m),nzmax)
      end do
      call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1(ii))
    end do
  endif
  
  ztout = real(ztemp1(1),kind=dp)*(1.-frac)+real(ztemp1(2),kind=dp)*frac
end subroutine zeta_to_z_lin
#endif
end module coord_ecmwf_mod
