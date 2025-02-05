  !*****************************************************************************
  !                                                                            *
  ! L. Bakels 2021: This module contains dry deposition related subroutines    *
  !                                                                            *
  !*****************************************************************************

module wetdepo_mod
  use point_mod
  use par_mod
  use com_mod
  use particle_mod

  implicit none

contains

subroutine wetdepo(itime,ltsample,loutnext)
  !                  i      i        i
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  ! 2021 Andreas Plach: - moved backward wet depo. calc. here from timemanager *
  !                     - bugfix in-cloud scavenging                           *
  !                                                                            *
  ! PS, AP 2021: followed up on some variable renaming and                     *
  !                corrected get_wetscav subroutine parameter list             *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! ldeltat [s]        interval since radioactive decay was computed           *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! wetdeposit         mass that is wet deposited                              *
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
#ifdef _OPENMP
  use omp_lib
#endif
  use unc_mod

  implicit none

  integer :: i,jpart,itime,ltsample,loutnext,ldeltat
  integer :: itage,nage,inage,ithread,thread
  integer :: ks, kp, stat
  real :: gridfract,wetscav,restmass
  real,allocatable,dimension(:) :: wettmp
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

  ! Loop over all particles
  !************************

#ifdef _OPENMP
  call omp_set_num_threads(numthreads_grid)
#endif

!$OMP PARALLEL PRIVATE(jpart,itage,nage,inage,ks,kp,thread,wetscav,wettmp, &
!$OMP restmass, gridfract)

#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM() ! Starts with 0
#else
    thread = 0
#endif

  allocate( wettmp(nspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate wettmp inside of OMP loop'

!$OMP DO 
  do i=1,count%alive

    jpart=count%ialive(i)
    
  ! Determine age class of the particle - nage is used for the kernel
  !******************************************************************
    itage=abs(itime-part(jpart)%tstart)
    nage=1
    if (lagespectra.eq.1) then
      do inage=1,nageclass
        nage=inage
        if (itage.lt.lage(nage)) exit
      end do
    endif

    do ks=1,nspec      ! loop over species

      if (.not. WETDEPSPEC(ks)) cycle 

  !**************************************************
  ! CALCULATE DEPOSITION 
  !**************************************************

      call get_wetscav(itime,jpart,ks,gridfract,wetscav)

      if (WETBKDEP) then
        if ((xscav_frac1(jpart,ks).lt.-0.1)) then   ! particle marked as starting particle
          if (wetscav.gt.0.) then
             xscav_frac1(jpart,ks)=wetscav*(zpoint2(part(jpart)%npoint)-&
               zpoint1(part(jpart)%npoint))*gridfract
          else
            mass(jpart,ks)=0.
            xscav_frac1(jpart,ks)=0.
          endif
        endif
      endif

      if (wetscav.gt.0.) then
        wettmp(ks)=mass(jpart,ks)* &
             (1.-exp(-wetscav*abs(ltsample)))*gridfract  ! wet deposition

      else ! if no scavenging
        wettmp(ks)=0.
      endif
      wetdeposit(jpart,ks)=wetdeposit(jpart,ks)+wettmp(ks)
      restmass = mass(jpart,ks)-wettmp(ks)
      if (ioutputforeachrelease.eq.1) then
        kp=part(jpart)%npoint
      else
        kp=1
      endif
      if (restmass .gt. smallnum) then
        mass(jpart,ks)=restmass
      else
        mass(jpart,ks)=0.
      endif
  !   Correct deposited mass to the last time step when radioactive decay of
  !   gridded deposited mass was calculated
      if (decay(ks).gt.0.) then
        wettmp(ks)=wettmp(ks)*exp(abs(ldeltat)*decay(ks))
      endif

    end do ! loop over species

    !************************************************************************
    ! Sabine Eckhardt, June 2008 create deposition only for forward runs
    ! Add the wet deposition to accumulated amount on output grid 
    !                                      and nested output grid

    if ((ldirect.eq.1).and.(iout.ne.0)) then !OMP reduction necessary for wetgridunc
      call wetdepokernel(part(jpart)%nclass,wettmp,real(part(jpart)%xlon), &
           real(part(jpart)%ylat),nage,kp,thread+1)
      if (nested_output.eq.1) call wetdepokernel_nest(part(jpart)%nclass, &
           wettmp,real(part(jpart)%xlon),real(part(jpart)%ylat),nage,kp,thread+1)
    endif

  end do ! all particles

!$OMP END DO
  deallocate(wettmp)
!$OMP END PARALLEL

#ifdef _OPENMP
  call omp_set_num_threads(numthreads)
#endif

#ifdef _OPENMP
    if ((ldirect.eq.1).and.(iout.ne.0)) then
      do ithread=1,numthreads_grid
        wetgridunc(:,:,:,:,:,:)=wetgridunc(:,:,:,:,:,:)+gridunc_omp(:,:,1,:,:,:,:,ithread)
        gridunc_omp(:,:,1,:,:,:,:,ithread)=0.
      end do
      if (nested_output.eq.1) then
        do ithread=1,numthreads_grid
          wetgriduncn(:,:,:,:,:,:)=wetgriduncn(:,:,:,:,:,:)+griduncn_omp(:,:,1,:,:,:,:,ithread)
          griduncn_omp(:,:,1,:,:,:,:,ithread)=0.
        end do
      endif
    endif
#endif

end subroutine wetdepo

subroutine get_wetscav(itime,jpart,ks,gridfract,wetscav)
  !                      i      i        i       i    i    o         o
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  ! ZHG, for v10: use below-cloud scavenging according to Laakso et al (2003)  *
  !   and Kyro et al (2009) as described in Grytte et al (2017, GMD)           *
  !                                                                            *
  ! PS, AP 04/2019: - put back temporal interpolation of rain, from v10.01     *
  !                 - tansferred BCSCHEME parameters to par_mod.f90            *
  !                 - added call to rain interpolation subroutine with new     *
  !                    interpolation for rain and all cloud params             *
  !                 - cleaned up scavenging determination algorithm            *
  !                 - added new below-cloud scavenging scheme                  *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! cc [0-1]           total cloud cover                                       *
  ! convp [mm/h]       convective precipitation rate                           *
  ! grfraction [0-1]   fraction of grid, for which precipitation occurs        *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! lfr, cfr           area fraction covered by precipitation for large scale  *
  !                    and convective precipitation (dependent on prec. rate)  *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! lsp [mm/h]         large-scale precipitation rate                          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! prec [mm/h]        precipitation rate in subgrid, where precipitation occurs*
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use interpol_mod
  use windfields_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif

  implicit none

  integer,intent(in) :: jpart,itime,ks
  real,intent(out) :: gridfract,wetscav

  integer :: hz,interp_time, n,i,j,kz
  integer(kind=1) :: clouds_v
  integer(selected_int_kind(16)), dimension(nspec) :: blc_count, inc_count

  integer :: indcloud
  integer :: icbot,ictop
  real :: t_particle, si, cl, cle ! in cloud scavenging
  real :: lsp,convp,cc,prec
  !save lfr,cfr
  real :: xts,yts

  real, parameter :: precsub = 0.01 ! minimum precip rate (mm/h)

  real :: f

  logical :: readclouds_this_nest

  wetscav=0.

  ! Interpolate large-scale precipitation, convective precipitation,
  ! total cloud cover, particle temperature, cloud water content, 
  ! cloud bottom and top 
  !********************************************************************
  interp_time=itime

  xts=real(part(jpart)%xlon)
  yts=real(part(jpart)%ylat)

  ! Determine which nesting level to be used
  !*****************************************
  readclouds_this_nest=.false.

  call find_ngrid(xts,yts)
  if ( (ngrid.gt.0) ) then
    if (lcw_nest(ngrid)) readclouds_this_nest=.true.
  endif

  ! If point at border of grid -> small displacement into grid
  !***********************************************************
  if (ngrid.le.0) then
    if (xts.ge.real(nx-1)) xts=real(nx-1)-0.00001
    if (yts.ge.real(ny-1)) yts=real(ny-1)-0.00001
  else
    if (xts.ge.real(nx-1)) xts=real(nx-1)-0.00001
    if (yts.ge.real(ny-1)) yts=real(ny-1)-0.00001
  endif

  call find_grid_indices(xts,yts)
  call find_grid_distances(xts,yts)
  
#ifdef ETA
  call update_zeta_to_z(itime,jpart)
  call find_z_level_eta_uv(real(part(jpart)%zeta))
  kz=induv
#else
  call find_z_level_meters(real(part(jpart)%z))
  kz=indz
#endif
  
  ! Interpolate cloud information
  call interpol_rain(itime,kz,lsp,convp,cc,t_particle,cl,icbot,ictop,icmv)
  ! cc = total cloud cover
  ! cl = ctwc

! If total precipitation is less than precsub=0.01 mm/h - no scavenging
! Note: original PS version (in order avoid step at 0.01)
!-----------------------------------------------------------------------
  prec = lsp+convp
  if (prec .le. precsub) then
    return
  endif

  ! Remove the minimum 0.01 from the large scale precipitation (lsp) and 
  ! convective precipitation (convp). Why 0.01???
  f = (prec-precsub)/prec
  lsp   = f*lsp
  convp = f*convp

  if (abs(memtime(1)-interp_time) .lt. abs(memtime(2)-interp_time)) then
    n=memind(1)
  else
    n=memind(2)
  endif

  ! if particle is above the clouds no scavenging is done
  !------------------------------------------------------
  ! PS: part of 2011/2012 fix 
  ! NOTE this is just for z coordinate
  ! Reverse sign for eta
#ifdef ETA
  if   (part(jpart)%zeta .gt. real(ictop)/eta_convert) then
    if (part(jpart)%zeta .le. real(icbot)/eta_convert) then
#else
  if   (part(jpart)%z .le. real(ictop)) then
    if (part(jpart)%z .gt. real(icbot)) then
#endif
      indcloud = 2 ! in-cloud
    else
      indcloud = 1 ! below-cloud
    endif
  elseif (ictop .eq. icmv) then
    indcloud = 0 ! no cloud found, use old scheme
  else
    return ! above cloud
  endif

  ! 1) Parameterization of the the area fraction of the grid cell where the
  !    precipitation occurs: the absolute limit is the total cloud cover, but
  !    for low precipitation rates, an even smaller fraction of the grid cell
  !    is used. Large scale precipitation occurs over larger areas than
  !    convective precipitation.
  !**************************************************************************
  call get_gridfract(lsp,convp,cc,gridfract)

  ! 2) Computation of precipitation rate in sub-grid cell
  !******************************************************
  prec=(lsp+convp)/gridfract

  ! 3) Computation of scavenging coefficients for all species
  !**********************************************************
  !-------------------------------------------------------
  if (indcloud .eq. 0) then ! NO CLOUD FOUND
  !-------------------------------------------------------
  ! Note: more complex formulation using H or particle diametre
  !       may be introduced later
    wetscav=wet_a*prec**wet_b

  !-------------------------------------------------------
  else if (indcloud .eq. 1) then ! BELOW CLOUD SCAVENGING
  !-------------------------------------------------------
    if (dquer(ks).le.0. .and. &
      weta_gas(ks).gt.0. .and. wetb_gas(ks).gt.0.) then
      ! gas: if positive below-cloud parameters (A or B), and dquer<=0
      call get_wetscav_belowcld_gas(ks,prec,wetscav)

    else if (dquer(ks).gt.0. .and. &
      (crain_aero(ks).gt.0. .or. csnow_aero(ks).gt.0.)) then
      ! aerosols: if positive below-cloud parameters (Crain/Csnow or B), and dquer>0
      
      if (t_particle .ge. 273. .and. crain_aero(ks).gt.0.) then ! Rain:
        call get_wetscav_belowcld_aerosol_rain(ks,prec,wetscav)
      
      else if (t_particle .lt. 273. .and. csnow_aero(ks).gt.0.) then ! Snow:
        call get_wetscav_belowcld_aerosol_snow(ks,prec,wetscav)
      !else ????????
      endif      
    endif ! gas or particle
    ! positive below-cloud scavenging parameters given in Species file
    ! end below cloud scavening
  !---------------------------------------------------------
  elseif (indcloud .eq. 2) then !  IN CLOUD SCAVENGING
  !---------------------------------------------------------

    ! if negative coefficients (turned off) set to zero for use in equation
    if (ccn_aero(ks).lt.0.) ccn_aero(ks)=0.
    if (in_aero(ks).lt.0.) in_aero(ks)=0.

    if (dquer(ks).gt.0) then !aerosol
      call get_wetscav_incld_aerosol(ks,gridfract,prec,cl,cc,t_particle,wetscav)
    else !gas
      call get_wetscav_incld_gas(ks,gridfract,prec,cl,cc,t_particle,wetscav)
    endif
!---------------------------------------------------------
  endif ! incloud
!---------------------------------------------------------

end subroutine get_wetscav

subroutine get_wetscav_belowcld_gas(ks,prec,wetscav)
  implicit none

  integer,intent(in) :: ks ! Species index
  real, intent(in) :: prec ! precipitation in sub-grid cell
  real, intent(out) :: wetscav ! scavenging coefficient

  ! Sum number of particles below the cloud
  icnt_belowcld(ks)=icnt_belowcld(ks)+1
  !weta_gas and wetb_gas are set in the SPECIES file
  wetscav=weta_gas(ks)*prec**wetb_gas(ks)
end subroutine get_wetscav_belowcld_gas

subroutine get_wetscav_belowcld_aerosol_rain(ks,prec,wetscav)
  implicit none

  integer,intent(in) :: ks ! Species index
  real, intent(in) :: prec ! precipitation in sub-grid cell
  real, intent(out) :: wetscav ! scavenging coefficient


  real :: dquer_m, ldquer
  real :: wetscavlim, logAd, B

  ! Sum number of particles below the cloud
  icnt_belowcld(ks)=icnt_belowcld(ks)+1

  ! NIK 17.02.2015: local conversion particle diameter from um 
  ! (SPECIES file, readspecies) to meter
  dquer_m = dquer(ks)*1.e-6

  ! The parameterizations used by HG scheme are valid only for d < 10 um
  ! Therefore, locally the diametre is clipped. However, settling and dry dep
  ! are still calculated with the original diametre
  ! TODO check whether warning is written by readrelease for d > 10 um
  dquer_m = min( 10.e-6, dquer_m )

  ! PS note that solid precip is usually expected for T<Tf, not T<273
  ! also, if freezing/melting point is desired, why not 273.2?

  ! AT parameterization after WANG ET AL 2014   
  !    unit of dquer is in um
  !    unit of precip is in mm/h
  ! Wang et al. 2014: eq 6+7
  ldquer = log10(dquer(ks))
  if (dquer(ks) .le. 2.) then
   logAd = bclr_a(1)              + &
           bclr_a(2) * ldquer     + &
           bclr_a(3) * ldquer**2. + & 
           bclr_a(4) * ldquer**3.
     
   B     = bclr_c(1) + bclr_c(2)*ldquer
    
  else ! dquer .gt. 2.
   logAd = bclr_b(1)              + &
           bclr_b(2) * ldquer     + &
           bclr_b(3) * ldquer**2. + &
           bclr_b(4) * ldquer**3. + &
           bclr_b(5) * ldquer**4. + &
           bclr_b(6) * ldquer**5. + &
           bclr_b(7) * ldquer**6.
            
   B    = bclr_e(1)              + &
          bclr_e(2) * ldquer     + &
          bclr_e(3) * ldquer**2. + &
          bclr_e(4) * ldquer**3. + &
          bclr_e(5) * ldquer**4. + &
          bclr_e(6) * ldquer**5. + &
          bclr_e(7) * ldquer**6.

  endif ! dquer

  ! Wang et al. 2014: eq. 4
  wetscav = 10**(logAd+B*log10(prec))           

end subroutine get_wetscav_belowcld_aerosol_rain

subroutine get_wetscav_belowcld_aerosol_snow(ks,prec,wetscav)
  implicit none

  integer,intent(in) :: ks ! Species index
  real, intent(in) :: prec ! precipitation in sub-grid cell
  real, intent(out) :: wetscav ! scavenging coefficient


  real :: dquer_m, ldquer
  real :: wetscavlim, logAd, B

  ! Sum number of particles below the cloud
  icnt_belowcld(ks)=icnt_belowcld(ks)+1

  ! NIK 17.02.2015: local conversion particle diameter from um 
  ! (SPECIES file, readspecies) to meter
  dquer_m = dquer(ks)*1.e-6

  ! The parameterizations used by HG scheme are valid only for d < 10 um
  ! Therefore, locally the diametre is clipped. However, settling and dry dep
  ! are still calculated with the original diametre
  ! TODO check whether warning is written by readrelease for d > 10 um
  dquer_m = min( 10.e-6, dquer_m )

  ! AT parameterization after WANG ET AL 2014   
  !    unit of dquer is in um
  !    unit of precip is in mm/h
  ldquer = log10(dquer(ks))
  ! Wang et al. 2014: eq. 8+9
  if (dquer(ks) .le. 1.44) then   
   logAd = bcls_a(1)               + &
           bcls_a(2)  * ldquer     + &
           bcls_a(3)  * ldquer**2. + &
           bcls_a(4)  * ldquer**3. + &
           bcls_a(5)  * ldquer**4. + &
           bcls_a(6)  * ldquer**5. + &
           bcls_a(7)  * ldquer**6.
     
   B     = bcls_c(1)              + &
           bcls_c(2) * ldquer     + &
           bcls_c(3) * ldquer**2. + &
           bcls_c(4) * ldquer**3. + &
           bcls_c(5) * ldquer**4. + &
           bcls_c(6) * ldquer**5. + &
           bcls_c(7) * ldquer**6.
     
  else ! dquer .gt. 1.44
   logAd = bcls_b(1)              + &
           bcls_b(2) * ldquer     + &
           bcls_b(3) * ldquer**2. + &
           bcls_b(4) * ldquer**3. + &
           bcls_b(5) * ldquer**4. + &
           bcls_b(6) * ldquer**5. + &
           bcls_b(7) * ldquer**6.
             
   B    = bcls_e(1)              + &
          bcls_e(2) * ldquer     + &
          bcls_e(3) * ldquer**2. + &
          bcls_e(4) * ldquer**3. + &
          bcls_e(5) * ldquer**4. + &
          bcls_e(6) * ldquer**5. + &
          bcls_e(7) * ldquer**6.
   
  endif ! dquer

  ! Wang et al. 2014: eq. 4
  wetscav = 10**(logAd+B*log10(prec))

end subroutine get_wetscav_belowcld_aerosol_snow

subroutine get_wetscav_incld_aerosol(ks,gridfract,prec,cl,cc,t_particle,wetscav)
  implicit none

  integer,intent(in) :: ks
  real,intent(in) :: gridfract
  real,intent(in) :: t_particle ! temperature
  real, intent(in) :: prec,cc ! precipitation in sub-grid cell
  real, intent(inout) :: cl ! scavenging coefficient
  real,intent(out) :: wetscav

  real :: frac_act,si

  ! NIK 13 may 2015: do only if in-cloud aerosol scavenging parameters > 0
  ! (all defined in SPECIES)
  if (ccn_aero(ks)+in_aero(ks) .le. 0.) return

  icnt_incld(ks)=icnt_incld(ks)+1

  ! Compute cloud liquid and ice water
  call get_cloud_liquid(gridfract,prec,cc,cl)
  ! Compute actived fraction based on the temperature (rain vs. snow )
  call get_activated_frac(ks,t_particle, frac_act)
  !ZHG Use the activated fraction and the liqid water to calculate the washout
  si= frac_act/cl
  ! scavenging coefficient based on Hertel et al 1995 - 
  ! using si (S_i in paper) for both gas and aerosol

  ! wetscav = ratio_incloud*si*prec/3.6E6/cloud_thickness
  ! cloud_thickness cancels out since cl is computed without
  wetscav=ratio_incloud*si* prec/3.6E6 ! mm/h -> m/s
end subroutine get_wetscav_incld_aerosol

subroutine get_wetscav_incld_gas(ks,gridfract,prec,cl,cc,t_particle,wetscav)
  implicit none

  integer,intent(in) :: ks
  real,intent(in) :: gridfract
  real,intent(in) :: t_particle ! temperature
  real, intent(in) :: prec,cc ! precipitation in sub-grid cell
  real, intent(inout) :: cl ! scavenging coefficient
  real,intent(out) :: wetscav

  real :: frac_act,si,cle
  ! NIK 13 may 2015: do only if in-cloud aerosol scavenging parameters > 0
  ! and Henry's constant > 0 (all defined in SPECIES)
  if (ccn_aero(ks)+in_aero(ks)+henry(ks) .le. 0.) return

  icnt_incld(ks)=icnt_incld(ks)+1

  ! Compute cloud liquid and ice water
  call get_cloud_liquid(gridfract,prec,cc,cl)
  ! Compute actived fraction based on the temperature (rain vs. snow )
  call get_activated_frac(ks,t_particle, frac_act)

  !ZHG Use the activated fraction and the liqid water to calculate the washout
  cle=(1.-cl)/(henry(ks)*(r_air/3500.)*t_particle) + cl
  si=1./cle
  ! scavenging coefficient based on Hertel et al 1995 - 
  ! using si (S_i in paper) for both gas and aerosol

  ! wetscav = ratio_incloud*si*prec/3.6E6/cloud_thickness
  ! cloud_thickness cancels out since cl is computed without
  wetscav=ratio_incloud*si* prec/3.6E6 ! mm/h -> m/s
end subroutine get_wetscav_incld_gas

subroutine get_activated_frac(ks,t_particle,frac_act)
  implicit none
  integer, intent(in) :: ks
  real, intent(in) :: t_particle
  real, intent(out) :: frac_act
  real :: frac_liq, frac_ice
  ! AT use of correct Kelvin temperature for T_ice; after ECMWF
  if (t_particle .le. 250.16) then ! ice 
    frac_liq=0.
    frac_ice=1.
  ! AT use of correct Kelvin temperature for T_liquid; after ECMWF        
  elseif (t_particle .ge. 273.16) then ! liquid
    frac_liq=1.
    frac_ice=0.
  else ! mixed cloud
    ! Use exact value for the melting point and ice threshold as in ECMWF/IFS            
    frac_liq= ((t_particle-250.16)/(273.16-250.16))**2.
    frac_ice = max(0., 1. - frac_liq)
  endif
  frac_act = frac_liq*ccn_aero(ks) + frac_ice*in_aero(ks)
end subroutine get_activated_frac

subroutine get_gridfract(lsp,convp,cc,gridfract)
  ! 1) Parameterization of the the area fraction of the grid cell where the
  !    precipitation occurs: the absolute limit is the total cloud cover, but
  !    for low precipitation rates, an even smaller fraction of the grid cell
  !    is used. Large scale precipitation occurs over larger areas than
  !    convective precipitation.
  !**************************************************************************
  implicit none

  real,intent(in) :: lsp,cc,convp
  real, intent(out) :: gridfract
  ! Where do these numbers come from?
  real, parameter :: lfr(5) = (/ 0.5,0.65,0.8,0.9,0.95/)
  real, parameter :: cfr(5) = (/ 0.4,0.55,0.7,0.8,0.9 /)
  integer :: i, j

  if (.not. lgridfraction) then
    gridfract=1.0
    return
  endif

  if (lsp.gt.20.) then
    i=5
  else if (lsp.gt.8.) then
    i=4
  else if (lsp.gt.3.) then
    i=3
  else if (lsp.gt.1.) then
    i=2
  else
    i=1
  endif

  if (convp.gt.20.) then
    j=5
  else if (convp.gt.8.) then
    j=4
  else if (convp.gt.3.) then
    j=3
  else if (convp.gt.1.) then
    j=2
  else
    j=1
  endif

  ! In the future, we may differentiate the gridfract for lsp and convp
  ! for now they are still mixed
  gridfract=max( 0.05, cc* (lsp*lfr(i) + convp*cfr(j)) / (lsp+convp))
end subroutine get_gridfract

subroutine get_cloud_liquid(gridfract,prec,cc,cl)
  use interpol_mod

  implicit none

  real, intent(in) :: prec,cc,gridfract ! precipitation in sub-grid cell
  real, intent(inout) :: cl ! scavenging coefficient
  !ZHG 2015 use cloud liquid & ice water (CLWC+CIWC) from ECMWF

  ! MC -- Integrated water content:
  ! CTWC = SUM(CLWC * rho_water * cloud_thickness) [kg/kg * kg/m3 * m]
  ! -> Average water content: cl = CTWC/rho_water/cloud_thickness [m3(water)/m3(cloud)]

  ! Note that cloud_thickness is not included, since this will cancel out when
  ! computing the wetscavenging: Wetscav=ratio_incloud*S_i*(prec/3.6E6)/cloud_thickness
  !                              S_i=1/cl
  ! Mother grid
  if (ngrid.le.0) then
    if (lcw) then
      if (lgridfraction) then
        cl=cl/rho_water*(gridfract/cc)
      else
        cl=cl/rho_water ! Grythe et al. eq (1), no cloud_thickness since this will cancel out later
      endif
      ! cl = cl*(gridfract/cc)
      ! A.Plach 2021 cl should not become too small
      ! cl=max(0.2*prec**0.36, cl*(gridfract/cc))
    else ! no cloud water available, use parameterisation for cloud water [m2/m3]
      cl=0.2*prec**0.36 !max(0.2*prec**0.36, cl*(gridfract/cc))
      ! ZHG updated parameterization to better reproduce the values from ECMWF
      ! cl = 1.E6*2E-7*prec**0.36 ! SEC ECMWF new, is also suitable for GFS
      ! SEC test:
      ! cl=1.E6*1.E-7*prec**0.3  ! SEC GFS new
      ! cl=     2.E-7*prec**0.36 ! Andreas
      ! cl=    1.6E-6*prec**0.36 ! Henrik
    endif
  else
    if (lcw_nest(ngrid)) then
      if (lgridfraction) then
        cl=cl/rho_water*(gridfract/cc)
      else
        cl=cl/rho_water ! Grythe et al. eq (1), no cloud_thickness since this will cancel out later
      endif
      !cl = cl*(gridfract/cc)
    else
      cl=0.2*prec**0.36
    endif
  endif
end subroutine get_cloud_liquid

subroutine wetdepokernel(nunc,deposit,x,y,nage,kp,thread)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     deposition fields using a uniform kernel with bandwidths dxout and dyout.*
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

  use unc_mod

  implicit none

  integer,intent(in) :: thread
  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,nunc,nage,ks,kp

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
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
           (jy.le.numygrid-1)) then
#ifdef _OPENMP
        gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
             gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)
#else
        wetgridunc(ix,jy,ks,kp,nunc,nage)= &
             wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
#endif
      end if
    end do
  else ! use kernel 
    
  ! Determine mass fractions for four grid points
  !**********************************************

  do ks=1,nspec

    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=wx*wy
#ifdef _OPENMP
      gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
           gridunc_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ix,jy,ks,kp,nunc,nage)= &
           wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
      gridunc_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)= &
           gridunc_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ixp,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=(1.-wx)*wy
#ifdef _OPENMP
      gridunc_omp(ixp,jy,1,ks,kp,nunc,nage,thread)= &
           gridunc_omp(ixp,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ixp,jy,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=wx*(1.-wy)
#ifdef _OPENMP
      gridunc_omp(ix,jyp,1,ks,kp,nunc,nage,thread)= &
           gridunc_omp(ix,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ix,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

  end do
  end if
end subroutine wetdepokernel

subroutine wetdepokernel_nest(nunc,deposit,x,y,nage,kp,thread)
  !                           i    i       i i i    i
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
  !      2 September 2004: Adaptation from wetdepokernel.                      *
  !                                                                            *
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

  use unc_mod

  implicit none

  integer,intent(in) :: thread
  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,kp,nunc,nage

  xl=(x*dx+xoutshiftn)/dxoutn
  yl=(y*dy+youtshiftn)/dyoutn

  ! old:
  ! ix=int(xl) 
  ! jy=int(yl)

  ! ESO: for xl,yl in range <-.5,-1> we get ix,jy=0 and thus negative
  ! wx,wy as the int function rounds upwards for negative numbers.
  ! Either use the floor function, or (perhaps more correctly?) use "(xl.gt.-0.5)" 
  ! in place of "(ix.ge.0)" and similar for the upper boundary.

  ! new:
  ix=floor(xl)
  jy=floor(yl)

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
    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=wx*wy
#ifdef _OPENMP
      griduncn_omp(ix,jy,1,ks,kp,nunc,nage,thread)= &
           griduncn_omp(ix,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ix,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
      griduncn_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)= &
           griduncn_omp(ixp,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ixp,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=(1.-wx)*wy
#ifdef _OPENMP
      griduncn_omp(ixp,jy,1,ks,kp,nunc,nage,thread)= &
           griduncn_omp(ixp,jy,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ixp,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=wx*(1.-wy)
#ifdef _OPENMP
      griduncn_omp(ix,jyp,1,ks,kp,nunc,nage,thread)= &
           griduncn_omp(ix,jyp,1,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ix,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

  end do
end subroutine wetdepokernel_nest

subroutine writeprecip(itime,imem)

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file containing total precipitation for each      *
  !  releases point.                                                           *
  !                                                                            *
  !     Author: S. Eckhardt                                                    * 
  !     7 Mai 2017                                                             *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use date_mod
  use windfields_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,itime,i
  real(kind=dp) :: jul
  character :: adate*8,atime*6

  integer :: ix,jy,imem
  real :: xp1,yp1

  
  if (itime.eq.0) then
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       form='formatted',err=998)
  else
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       ACCESS='APPEND',form='formatted',err=998)
  endif

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0 !lat, long (real) coord
    yp1=ypoint1(i)*dy+ylat0 !lat, long (real) coord
    ix=int((xpoint1(i)+xpoint2(i))*0.5)
    jy=int((ypoint1(i)+ypoint2(i))*0.5)
    write(unitprecip,*)  jjjjmmdd, ihmmss,xp1,yp1, & 
      sum(lsprec(ix,jy,1,:,imem)),sum(convprec(ix,jy,1,:,imem))
 !time is the same as in the ECMWF windfield
! units mm/h, valid for the time given in the windfield
  enddo

  close(unitprecip)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header_txt'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  error stop
end subroutine writeprecip

end module wetdepo_mod
