! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
!   L. Bakels 2022: This module contains the computation of particle         *
!                   trajectories                                             *
!                                                                            *
!*****************************************************************************

module advance_mod
  use point_mod
  use par_mod
  use com_mod
  use interpol_mod
  use cmapf_mod
  use random_mod, only: ran3,iseed1
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use particle_mod
  use turbulence_mod
  use settling_mod
  use windfields_mod, only: nxmax

  implicit none 
    real, parameter ::              &
      eps2=1.e-9,                   &
      eps3=tiny(1.0),               &
      eps_eta=1.e-4
    real ::                         &
      eps

  private :: adv_above_pbl, adv_in_pbl, petterssen_corr, update_xy, pushpartdown

contains
  
subroutine advance(itime,ipart,ithread)

  !*****************************************************************************
  !                                                                            *
  !  Calculation of turbulent particle trajectories utilizing a                *
  !  zero-acceleration scheme, which is corrected by a numerically more        *
  !  accurate Petterssen scheme whenever possible.                             *
  !                                                                            *
  !  Particle positions are read in, incremented, and returned to the calling  *
  !  program.                                                                  *
  !                                                                            *
  !  In different regions of the atmosphere (PBL vs. free troposphere),        *
  !  different parameters are needed for advection, parameterizing turbulent   *
  !  velocities, etc. For efficiency, different interpolation routines have    *
  !  been written for these different cases, with the disadvantage that there  *
  !  exist several routines doing almost the same. They all share the          *
  !  included file 'interpol_mod'. The following                               *
  !  interpolation routines are used:                                          *
  !                                                                            *
  !  interpol_all(_nest)     interpolates everything (called inside the PBL)  *
  !  interpol_misslev(_nest) if a particle moves vertically in the PBL,       *
  !                           additional parameters are interpolated if it     *
  !                           crosses a model level                            *
  !  interpol_wind(_nest)    interpolates the wind and determines the         *
  !                           standard deviation of the wind (called outside   *
  !                           PBL) also interpolates potential vorticity       *
  !  interpol_wind_short(_nest) only interpolates the wind (needed for the    *
  !                           Petterssen scheme)                               *
  !  interpol_vdep(_nest)    interpolates deposition velocities               *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Changes:                                                                  *
  !                                                                            *
  !  8 April 2000: Deep convection parameterization                            *
  !                                                                            *
  !  May 2002: Petterssen scheme introduced                                    *
  !                                                                            *
  !  2021, L. Bakels:                                                          *
  !         - Separated PBL and above PBL computations in different            *
  !           subroutines                                                      *
  !         - Moved all turbulence computations to turbulence_mod.f90          *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! icbt               1 if particle not transferred to forbidden state,       *
  !                    else -1                                                 *
  ! dawsave            accumulated displacement in along-wind direction        *
  ! dcwsave            accumulated displacement in cross-wind direction        *
  ! dxsave             accumulated displacement in longitude                   *
  ! dysave             accumulated displacement in latitude                    *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          time at which this subroutine is entered                *
  ! itimec [s]         actual time, which is incremented in this subroutine    *
  ! href [m]           height for which dry deposition velocity is calculated  *
  ! ladvance [s]       Total integration time period                           *
  ! ldirect            1 forward, -1 backward                                  *
  ! ldt [s]            Time step for the next integration                      *
  ! lsynctime [s]      Synchronisation interval of FLEXPART                    *
  ! ngrid              index which grid is to be used                          *
  ! nrand              index for a variable to be picked from rannumb          *
  ! nstop              if > 1 particle has left domain and must be stopped     *
  ! prob               probability of absorption due to dry deposition         *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! rhoa               air density                                             *
  ! rhograd            vertical gradient of the air density                    *
  ! up,vp,wp           random velocities due to turbulence (along wind, cross  *
  !                    wind, vertical wind                                     *
  ! usig,vsig,wsig     mesoscale wind fluctuations                             *
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  implicit none
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart,                        & ! particle index
    ithread                         ! OMP thread starting at 0
  integer ::                      &
    itimec,                       &
    i,j,                          & ! loop variables
    nrand,                        & ! random number used for turbulence
    memindnext,                   & ! seems useless
    ngr                             ! temporary new grid index of moved particle
  real ::                         &
    ux,vy,                        & ! random turbulent velocities above PBL
    tropop,                       & ! height of troposphere
    dxsave,dysave,                & ! accumulated displacement in long and lat
    dawsave,dcwsave,              & ! accumulated displacement in wind directions
    ztseta,wsigeta_tmp              ! real of eta z position
  logical ::                      &
    abovePBL
    ! flag will be set to 'true' if computation needs to be completed above PBL
#ifdef ETA
  real :: weta_settling             ! Settling velocity in eta coordinates
#endif

  eps=nxmax/3.e5

  part(ipart)%nstop=.false.
  do i=1,nmixz
    indzindicator(i,ithread+1)=.true.
  end do
  
  if (DRYDEP) then    ! reset probability for deposition
    depoindicator(:,ithread+1)=.true.
    prob(ipart,:)=0.
  endif
  
  if (lsettling) part(ipart)%settling=0.

  dxsave=0.           ! reset position displacements
  dysave=0.           ! due to mean wind
  dawsave=0.          ! and turbulent wind
  dcwsave=0.

  itimec=itime

  nrand=int(ran3(iseed1(ithread),ithread)*real(maxrand-1))+1

  ! Determine whether lat/long grid or polarstereographic projection
  ! is to be used
  ! Furthermore, determine which nesting level to be used
  !*****************************************************************
  call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)

  !***************************
  ! Interpolate necessary data
  !***************************

  if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
    memindnext=1
  else
    memindnext=2
  endif

  ! Convert z(eta) to z(m) for the turbulence scheme, w(m/s) 
  ! is computed in verttransform_ecmwf.f90
#ifdef ETA
  call update_zeta_to_z(itime,ipart)
  ztseta=real(part(ipart)%zeta)
#else
  ztseta=0.
#endif
  ! Determine nested grid coordinates
  ! Determine the lower left corner and its distance to the current position
  ! Calculate variables for time interpolation
  !*******************************************
  call init_interpol(itime, &
    real(part(ipart)%xlon),real(part(ipart)%ylat), &
    real(part(ipart)%z),ztseta)

  ! Compute maximum mixing height around particle position
  !*******************************************************

  ! Compute height of troposphere and PBL at x-y location of particle
  call interpol_htropo_hmix(tropop,h)
  zeta=real(part(ipart)%z)/h

  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************
  ! In the PBL we use meters instead of eta coordinates for vertical transport
  
  abovePBL=.true.
  if (zeta.le.1.) then
  
    abovePBL=.false.

    call adv_in_pbl(itime,itimec,&
      dxsave,dysave,dawsave,dcwsave,abovePBL,nrand,ipart,ithread)
#ifdef ETA
    if (lsettling) then
      call w_to_weta(itime,real(part(ipart)%idt),part(ipart)%xlon, &
        part(ipart)%ylat,part(ipart)%z,part(ipart)%zeta, &
        part(ipart)%settling,weta_settling)
      weta=weta+weta_settling
    endif
#endif
  endif 

  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************

  ! Interpolate the wind
  !*********************
  
  if (abovePBL) call adv_above_pbl(itime,itimec,dxsave,dysave, &
    ux,vy,tropop,nrand,ipart)
    ! Above PBL computation

  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************

  ! Mesoscale wind velocity fluctuations are obtained by scaling
  ! with the standard deviation of the grid-scale winds surrounding
  ! the particle location, multiplied by a factor fturbmeso.
  ! The autocorrelation time constant is taken as half the
  ! time interval between wind fields
  !****************************************************************

  if (lturbulence.eq.1) then 
    ! mesoscale turbulence is found to give issues, so turned off
    if (lmesoscale_turb) then
#ifdef ETA
      ztseta=real(part(ipart)%zeta)
      wsigeta_tmp=wsigeta
#else
      ztseta=0.
      wsigeta_tmp=0.
#endif
      call interpol_mesoscale( &
        real(part(ipart)%xlon),real(part(ipart)%ylat), &
        real(part(ipart)%z),ztseta)
      call turbulence_mesoscale(nrand,dxsave,dysave,ipart, &
        usig,vsig,wsig,wsigeta_tmp,eps_eta)
    endif

    !*************************************************************
    ! Transform along and cross wind components to xy coordinates,
    ! add them to u and v, transform u,v to grid units/second
    ! and calculate new position
    !*************************************************************

    call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
    dxsave=dxsave+ux   
     ! comment by MC: comment this line to stop particles horizontally for tests
    dysave=dysave+vy
  endif

  call update_xy(dxsave,dysave,ipart)
  if (part(ipart)%nstop) return

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************
  call pushpartdown(ipart)
  
  !************************************************************************
  ! Now we could finish, as this was done in FLEXPART versions up to 4.0.
  ! However, truncation errors of the advection can be significantly
  ! reduced by doing one iteration of the Petterssen scheme, if this is
  ! possible.
  ! Note that this is applied only to the grid-scale winds, not to
  ! the turbulent winds.
  !************************************************************************

  ! The Petterssen scheme can only applied with long time steps (only then u
  ! is the "old" wind as required by the scheme); otherwise do nothing
  !*************************************************************************

  if (part(ipart)%idt .ne. abs(lsynctime)) return

 ! The Petterssen scheme can only be applied if the ending time of the time step
 ! (itime+ldt*ldirect) is still between the two wind fields held in memory;
 ! otherwise do nothing
 !******************************************************************************

  if (abs(itime+part(ipart)%idt*ldirect).gt.abs(memtime(2))) return

  ! Apply it also only if starting and ending point of current time step are on
  ! the same grid; otherwise do nothing
  !*****************************************************************************
  ! ngr = ngrid 
  ! call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)

  if (nglobal .and. real(part(ipart)%ylat).gt.switchnorthg) then
    ngr=-1
  else if (sglobal.and. real(part(ipart)%ylat).lt.switchsouthg) then
    ngr=-2
  else
    ngr=0
    ! Temporary fix for nested layer edges: replaced eps with dxn and dyn (LB)
    do j=numbnests,1,-1
      if (real(part(ipart)%xlon).gt.xln(j)+dxn(j) .and. &
          real(part(ipart)%xlon).lt.xrn(j)-dxn(j) .and. &
          real(part(ipart)%ylat).gt.yln(j)+dyn(j) .and. &
          real(part(ipart)%ylat).lt.yrn(j)-dyn(j)) then
        ngr=j
        exit
      endif
    end do
  endif

  if (ngr.ne.ngrid) return

  call petterssen_corr(itime,ipart)

end subroutine advance

subroutine adv_above_pbl(itime,itimec,dxsave,dysave,ux,vy,tropop,nrand,ipart)

  implicit none
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  integer, intent(inout) ::       &
    itimec,                       & ! next timestep
    nrand                           ! random number used for turbulence
  real, intent(in) ::             &
    tropop                          ! height of troposphere
  real, intent(inout) ::          &
    ux,vy,                        & ! random turbulent velocities above PBL
    dxsave,dysave                   ! accumulated displacement in long and lat
  real ::                         &
    dt,                           & ! real(ldt)
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    wp                              ! random turbulence velocities
#ifdef ETA
  real :: weta_settling                   ! settling velocity in eta coordinates
#endif
  integer ::                      &
    insp,nsp                        ! loop variables for number of species

  zts=real(part(ipart)%z)
#ifdef ETA
  ztseta=real(part(ipart)%zeta)
#else
  ztseta=0.
#endif
  xts=real(part(ipart)%xlon)
  yts=real(part(ipart)%ylat)
  if (lsettling) part(ipart)%settling=0.

  call interpol_wind(itime,xts,yts,zts,ztseta)

  ! Compute everything for above the PBL

  ! Assume constant, uncorrelated, turbulent perturbations
  ! In the stratosphere, use a small vertical diffusivity d_strat,
  ! in the troposphere, use a larger horizontal diffusivity d_trop.
  ! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
  !******************************************************************

  part(ipart)%idt=abs(lsynctime-itimec+itime)
  dt=real(part(ipart)%idt)

  if (lturbulence.eq.1) then
    call turbulence_above_pbl(dt,nrand,ux,vy,wp,tropop,zts)
  else
    !sec switch off turbulence
    ux=0.0
    vy=0.0
    wp=0.0
  endif

  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases
  !*************************************************************************
  ! Does not work in eta coordinates yet
  if (mdomainfill.eq.0) then
    if (lsettling) then
      if ((ipin.ne.3).and.(ipin.ne.4)) then
        do insp=1,nspec
          nsp=insp
          if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
        end do
      else
        nsp=1
      endif
      ! LB change to eta coords?
      if (density(nsp).gt.0.) then
        call get_settling(xts,yts,zts,nsp,part(ipart)%settling)
#ifdef ETA
        call update_zeta_to_z(itime,ipart)
        if ((ldirect.eq.1).and.(part(ipart)%z+part(ipart)%settling*dt.lt.0.)) then
          part(ipart)%settling=-part(ipart)%z/dt
        endif
        call w_to_weta(itime,dt,part(ipart)%xlon,part(ipart)%ylat, &
          part(ipart)%z,part(ipart)%zeta,part(ipart)%settling,weta_settling)
        weta=weta+weta_settling
#else
        if ((ldirect.eq.1).and.(part(ipart)%z+part(ipart)%settling*dt.lt.0.)) then
          part(ipart)%settling=-part(ipart)%z/dt
        endif
        w=w+part(ipart)%settling
#endif
      end if
    endif
  end if

  ! Calculate position at time step itime+lsynctime
  !************************************************
  dxsave=dxsave+(u+ux)*dt
  dysave=dysave+(v+vy)*dt

#ifdef ETA
      if (wp.ne.0.) then
        call update_zeta_to_z(itime,ipart)
        call update_z(ipart,wp*dt*real(ldirect))
        if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))
          ! if particle below ground -> reflection
        call update_z_to_zeta(itime,ipart)
      endif
      call update_zeta(ipart,weta*dt*real(ldirect))
      if (part(ipart)%zeta.ge.1.) call set_zeta(ipart,1.-(part(ipart)%zeta-1.))
      if (part(ipart)%zeta.eq.1.) call update_zeta(ipart,-eps_eta)
#else
      call update_z(ipart,(w+wp)*dt*real(ldirect))
      if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))
#endif

end subroutine adv_above_pbl

subroutine adv_in_pbl(itime,itimec, dxsave,dysave,dawsave,dcwsave, abovePBL,  &
  nrand,ipart,ithread)

  use drydepo_mod, only: drydepo_probability

  implicit none

  logical, intent(inout) ::       &
    abovePBL                      
  ! flag will be set to 'true' if computation needs to be completed above PBL
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart,                        & ! particle index
    ithread                         ! number of the omp thread starting at 0
  real, intent(inout) ::          &
    dxsave,dysave,                & ! accumulated displacement in long and lat
    dawsave,dcwsave               ! accumulated displacement in wind directions
  integer, intent(inout) ::       &
    itimec,                       & ! next timestep
    nrand                           ! random number used for turbulence
  real ::                         &
    dt,                           & ! real(ldt)
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    rhoa,                         & ! air density, used in CBL
    rhograd                     ! vertical gradient of air density, used in CBL
  integer ::                      &
    loop,                         & ! loop variable for time in the PBL
    nsp,insp                        ! loop variable for species
  real :: vdepo(maxspec)  ! deposition velocities for all species

  eps=nxmax/3.e5
  if (lsettling) part(ipart)%settling=0.

  ! BEGIN TIME LOOP
  !================
  ! For wind_coord_type=ETA:
  ! Within this loop, only METER coordinates are used, and the new z value will
  ! be updated to ETA coordinates at the end
  !****************************************************************************
#ifdef ETA
  call update_zeta_to_z(itime,ipart)
  ztseta=real(part(ipart)%zeta)
#else
  ztseta=0.
#endif

  loop=0
  pbl_loop: do
  
    loop=loop+1
    if (method.eq.1) then
      part(ipart)%idt=min(part(ipart)%idt,abs(lsynctime-itimec+itime))
      itimec=itimec+part(ipart)%idt*ldirect
    else
      part(ipart)%idt=abs(lsynctime)
      itimec=itime+lsynctime
    endif
    dt=real(part(ipart)%idt)
    xts=real(part(ipart)%xlon)
    yts=real(part(ipart)%ylat)
    zts=real(part(ipart)%z)

    zeta=zts/h
    if (loop.eq.1) then ! Temporal interpolation only for the first iteration

      if (ngrid.le.0) then
        xts=real(part(ipart)%xlon)
        yts=real(part(ipart)%ylat)
        call interpol_pbl(itime,xts,yts,zts,ztseta,ithread+1)
      else
        call interpol_pbl(itime,xtn,ytn,zts,ztseta,ithread+1)
      endif

    else

      ! Determine the level below the current position for u,v,rho
      !***********************************************************
      call find_z_level_meters(zts)

      ! If one of the levels necessary is not yet available,
      ! calculate it
      !*****************************************************
      call interpol_pbl_misslev(ithread+1)

    endif

  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************

    call interpol_pbl_short(zts,rhoa,rhograd,ithread+1) ! Vertical interpolation

  ! Compute the turbulent disturbances
  ! Determine the sigmas and the timescales 
  !****************************************

    if (lturbulence.eq.1) then
      call turbulence_pbl(ipart,nrand,dt,zts,rhoa,rhograd,ithread) 
      ! Note: zts and nrand get updated

      ! Determine time step for next integration
      !*****************************************
      if (turbswitch) then
        part(ipart)%idt = int( &
          min( tlw, &
               h/max( 2.*abs(part(ipart)%turbvel%w*sigw), 1.e-5 ), &
               0.5/abs(dsigwdz) &
             ) *ctl)
      else
        part(ipart)%idt = int( &
          min( tlw, & 
               h/max( 2.*abs(part(ipart)%turbvel%w), 1.e-5) &
              ) *ctl)
      endif
    else
      part(ipart)%turbvel%u=0.
      part(ipart)%turbvel%v=0.
      part(ipart)%turbvel%w=0.
    endif

    part(ipart)%idt=max(part(ipart)%idt,mintime)


  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases, or if particle
  ! represents more than one species
  !*************************************************************************

    if (mdomainfill.eq.0) then
      if (lsettling) then
        if ((ipin.ne.3).and.(ipin.ne.4)) then
          do insp=1,nspec
            nsp=insp
            if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
          end do
        else
          nsp=1
        endif
        if (density(nsp).gt.0.) then
          call get_settling(xts,yts,zts,nsp,part(ipart)%settling)  !bugfix
          if ((ldirect.eq.1).and.(part(ipart)%z+part(ipart)%settling*dt.lt.0.)) then
            part(ipart)%settling=-part(ipart)%z/dt
          endif
          w=w+part(ipart)%settling
        end if
      end if
    endif

  ! Horizontal displacements during time step dt are small real values compared
  ! to the position; adding the two, would result in large numerical errors.
  ! Thus, displacements are accumulated during lsynctime and are added to the
  ! position at the end
  !****************************************************************************

    dxsave=dxsave+u*dt
    dysave=dysave+v*dt
    dawsave=dawsave+part(ipart)%turbvel%u*dt
    dcwsave=dcwsave+part(ipart)%turbvel%v*dt
    ! How can I change the w to w(eta) efficiently?

#ifdef ETA
    call update_z(ipart,w*dt*real(ldirect))
    zts=real(part(ipart)%z)
    ! HSO/AL: Particle managed to go over highest level -> interpolation
    ! error in goto 700
    !          alias interpol_wind (division by zero)
    if (zts.ge.height(nz)) call set_z(ipart,height(nz)-100.*eps) 
     ! Manually for z instead
#else
    call update_z(ipart,w*dt*real(ldirect))
    call pushpartdown(ipart)
#endif
    ! end select
    zts=real(part(ipart)%z)
    
    if (zts.gt.h) then
#ifdef ETA
      call update_z_to_zeta(itime,ipart)
#endif
      if (itimec.ne.itime+lsynctime) abovePBL=.true. 
        ! complete the current interval above PBL
      return 
    endif
    
  ! Determine probability of deposition
  !************************************
    if (DRYDEP) then
      call drydepo_probability(ipart,dt,zts,vdepo,ithread+1)
    endif

    if (zts.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))    
      ! if particle below ground -> reflection

    if (itimec.eq.(itime+lsynctime)) then
      ! Convert z position that changed by turbulent motions to eta coords
#ifdef ETA
      call update_z_to_zeta(itime,ipart)
#endif
      return  ! finished
    endif

  end do pbl_loop

#ifdef ETA
  call update_z_to_zeta(itime,ipart)
#endif

end subroutine adv_in_pbl

subroutine petterssen_corr(itime,ipart)

  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  integer ::                      &
    nsp,insp                        ! loop variables for number of species
  real ::                         &
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    uold,vold
#ifdef ETA
  real :: woldeta,weta_settling ! eta equivalents
#else
  real :: wold
#endif

  xts=real(part(ipart)%xlon)
  yts=real(part(ipart)%ylat)
  zts=real(part(ipart)%z)
#ifdef ETA
  ztseta=real(part(ipart)%zeta)
#else
  ztseta=0.
#endif
  if (lsettling) part(ipart)%settling=0.

  ! Memorize the old wind
  !**********************

  uold=u
  vold=v

#ifdef ETA
    woldeta=weta
#else
    wold=w
#endif

  ! Interpolate wind at new position and time
  !******************************************
  call interpol_wind_short(itime+part(ipart)%idt*ldirect,xts,yts,zts,ztseta)

  if (mdomainfill.eq.0) then
    if (lsettling) then
      if ((ipin.ne.3).and.(ipin.ne.4)) then
        do insp=1,nspec
          nsp=insp
          if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
        end do
      else
        nsp=1
      endif
      if (density(nsp).gt.0.) then
#ifdef ETA
        call update_zeta_to_z(itime+part(ipart)%idt,ipart)
        call update_z_to_zeta(itime+part(ipart)%idt,ipart)
        zts=real(part(ipart)%z)
        call get_settling(xts,yts,zts,nsp,part(ipart)%settling) !bugfix
        if ((ldirect.eq.1).and. &
          (part(ipart)%z+part(ipart)%settling*part(ipart)%idt.lt.0)) then
          part(ipart)%settling=-part(ipart)%z/part(ipart)%idt
        endif
        call w_to_weta( &
          itime+part(ipart)%idt, real(part(ipart)%idt), part(ipart)%xlon, &
          part(ipart)%ylat, part(ipart)%z, part(ipart)%zeta, &
          part(ipart)%settling, weta_settling)
        weta=weta+weta_settling
         !woldeta=
   !real(part(ipart)%zeta-part(ipart)%zeta_prev)/real(part(ipart)%idt*ldirect)
#else
        call get_settling(xts,yts,zts,nsp,part(ipart)%settling)
        if ((ldirect.eq.1).and. &
          (part(ipart)%z+part(ipart)%settling*part(ipart)%idt.lt.0)) then
          part(ipart)%settling=-part(ipart)%z/part(ipart)%idt
        endif
        w=w+part(ipart)%settling
#endif
      end if
    endif
  end if

  ! Determine the difference vector between new and old wind
  ! (use half of it to correct position according to Petterssen)
  !*************************************************************

  u=(u-uold)*0.5
  v=(v-vold)*0.5

#ifdef ETA
    weta=(weta-woldeta)*0.5
    call update_zeta(ipart,weta*real(part(ipart)%idt*ldirect))
    if (part(ipart)%zeta.ge.1.) call set_zeta(ipart,1.-(part(ipart)%zeta-1.))
    if (part(ipart)%zeta.eq.1.) call update_zeta(ipart,-eps_eta)
#else
    w=(w-wold)*0.5
    call update_z(ipart,w*real(part(ipart)%idt*ldirect))
    if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))          ! if particle below ground -> reflection
#endif

  ! Finally, correct the old position
  !**********************************

  call update_xy(u*part(ipart)%idt,v*part(ipart)%idt,ipart)

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************
  call pushpartdown(ipart)
  
end subroutine petterssen_corr

subroutine update_xy(xchange,ychange,ipart)

  implicit none

  integer, intent(in) ::          &
    ipart                           ! particle number
  real, intent(in) ::             &
    xchange,ychange                 ! change in position
  real ::                         &
    xlon,ylat,xpol,ypol,          & ! temporarily storing new particle positions
    gridsize,cosfact                ! used to compute new positions of particles

  eps=nxmax/3.e5

  if (ngrid.ge.0) then

    cosfact=dxconst/cos((real(part(ipart)%ylat)*dy+ylat0)*pi180)
    call update_xlon(ipart,real(xchange*cosfact*real(ldirect),kind=dp))
    call update_ylat(ipart,real(ychange*dyconst*real(ldirect),kind=dp))

  else if (ngrid.eq.-1) then      ! around north pole

    xlon=xlon0+real(part(ipart)%xlon)*dx !comment by MC: compute old part pos.
    ylat=ylat0+real(part(ipart)%ylat)*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)   
      !convert old particle position in polar stereographic
    gridsize=1000.*cgszll(northpolemap,ylat)   
      !calculate size in m of grid element in polar stereographic coordinate
    xpol=xpol+xchange/gridsize*real(ldirect)        
      !position in grid unit polar stereographic
    ypol=ypol+ychange/gridsize*real(ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)   
      !convert to lat long coordinate
    call set_xlon(ipart,real((xlon-xlon0)/dx,kind=dp))
      !convert to grid units in lat long coordinate, comment by mc
    call set_ylat(ipart,real((ylat-ylat0)/dy,kind=dp))

  else if (ngrid.eq.-2) then    ! around south pole

    xlon=xlon0+real(part(ipart)%xlon)*dx
    ylat=ylat0+real(part(ipart)%ylat)*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat)
    xpol=xpol+xchange/gridsize*real(ldirect)
    ypol=ypol+ychange/gridsize*real(ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    call set_xlon(ipart,real((xlon-xlon0)/dx,kind=dp))
    call set_ylat(ipart,real((ylat-ylat0)/dy,kind=dp))

  endif

  ! If global data are available, use cyclic boundary condition
  !************************************************************

  if (xglobal) then
    if (part(ipart)%xlon .ge. real(nxmin1, kind=dp)) &
      call update_xlon(ipart,-real(nxmin1, kind=dp))
    if (part(ipart)%xlon .lt. 0.) call update_xlon(ipart,real(nxmin1, kind=dp))
    if (part(ipart)%xlon .le. real(eps, kind=dp)) &
      call set_xlon(ipart,real(eps, kind=dp))
    if (abs( part(ipart)%xlon - real(nxmin1, kind=dp)) .le. eps) &
      call set_xlon(ipart,real(nxmin1-eps,kind=dp))
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************
  if (sglobal .and. part(ipart)%ylat.lt.0. ) then
    call set_xlon(ipart, &
      mod( part(ipart)%xlon + real(nxmin1*0.5, kind=dp), real(nxmin1, kind=dp)))
    call set_ylat(ipart,-part(ipart)%ylat)
   ! In extremely rare cases, the ylat exceeds the bounds, 
   ! so we set it back into the domain here
    if ( part(ipart)%ylat.gt.real(nymin1,kind=dp) ) &
      call set_ylat(ipart, &
        real(nymin1, kind=dp) - mod( part(ipart)%ylat, real(nymin1, kind=dp)))
  else if (nglobal .and. part(ipart)%ylat .gt. real(nymin1, kind=dp) ) then
    call set_xlon(ipart, &
      mod( part(ipart)%xlon + real(nxmin1*0.5, kind=dp), real(nxmin1, kind=dp)))
    call set_ylat(ipart,2.*real(nymin1,kind=dp) - part(ipart)%ylat)
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************
  ! Not necessary to check when using global domain, but some problems in the
  ! meteo data could cause particles to go crazy.
  ! if (gdomainfill) return 

  if (part(ipart)%xlon.lt.0. .or. part(ipart)%xlon.ge.real(nxmin1,kind=dp) &
 .or. part(ipart)%ylat.lt.0. .or. part(ipart)%ylat.gt.real(nymin1,kind=dp)) then
    part(ipart)%nstop=.true.
    return
  endif
  
end subroutine update_xy

subroutine pushpartdown(ipart)

  implicit none 

  integer, intent(in) ::          &
    ipart                           ! particle index

  eps=nxmax/3.e5

#ifdef ETA
  if (part(ipart)%zeta.le.real(uvheight(nz),kind=dp)) then
    if ((ldirect.eq.-1) .and. (lsettling)) then
      part(ipart)%nstop=.true.
    else
      call set_zeta(ipart,uvheight(nz)+eps_eta)
    endif
  endif
#else
  if (part(ipart)%z.ge.real(height(nz),kind=dp)) then
    if ((ldirect.eq.-1) .and. (lsettling)) then
      part(ipart)%nstop=.true.
    else
      call set_z(ipart,height(nz)-100.*eps)
    endif
  endif
#endif
  
end subroutine pushpartdown

end module advance_mod
