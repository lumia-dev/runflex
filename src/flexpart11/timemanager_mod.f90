! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !   L. Bakels 2022: This module contains the timemanager                     *
  !                                                                            *
  !*****************************************************************************

module timemanager_mod

implicit none

contains

subroutine timemanager

  !*****************************************************************************
  !                                                                            *
  ! Handles the computation of trajectories, i.e. determines which             *
  ! trajectories have to be computed at what time.                             *
  ! Manages dry+wet deposition routines, radioactive decay and the computation *
  ! of concentrations.                                                         *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 May 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !        Call of convmix when new windfield is read                          *
  !------------------------------------                                        *
  !  Changes Petra Seibert, Sept 2002                                          *
  !     fix wet scavenging problem                                             *
  !     Code may not be correct for decay of deposition!                       *
  !  Changes Petra Seibert, Nov 2002                                           *
  !     call convection BEFORE new fields are read in BWD mode                 *
  !  Changes Caroline Forster, Feb 2005                                        *
  !   new interface between flexpart and convection scheme                     *
  !   Emanuel's latest subroutine convect43c.f is used                         *
  !  Changes Stefan Henne, Harald Sodemann, 2013-2014                          *
  !   added netcdf output code                                                 *
  !  Changes Espen Sollum 2014                                                 *
  !   For compatibility with MPI version,                                      *
  !   variables uap,ucp,uzp,us,vs,ws,cbt now in module com_mod                 *
  !  Unified ECMWF and GFS builds                                              *
  !   Marian Harustak, 12.5.2017                                               *
  !  Changes L Bakels 2022: - OpenMP parallelisation                           *
  !                         - converting input to ETA coordinates              *
  !                         - spawning particles from part_ic.nc               *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! DEP                .true. if either wet or dry deposition is switched on   *
  ! decay(maxspec) [1/s] decay constant for radioactive decay                  *
  ! DRYDEP             .true. if dry deposition is switched on                 *
  ! ideltas [s]        modelling period                                        *
  ! itime [s]          actual temporal position of calculation                 *
  ! ldeltat [s]        time since computation of radioact. decay of depositions*
  ! loutaver [s]       averaging period for concentration calculations         *
  ! loutend [s]        end of averaging for concentration calculations         *
  ! loutnext [s]       next time at which output fields shall be centered      *
  ! loutsample [s]     sampling interval for averaging of concentrations       *
  ! loutstart [s]      start of averaging for concentration calculations       *
  ! loutstep [s]       time interval for which concentrations shall be         *
  !                    calculated                                              *
  ! loutrestart [s]    time interval for which restart files will be produced  *
  ! npoint             index, which starting point the trajectory has          *
  !                    starting positions of trajectories                      *
  ! nstop              serves as indicator for fate of particles               *
  !                    in the particle loop                                    *
  ! nstop1             serves as indicator for wind fields (see getfields)     *
  ! outnum             number of samples for each concentration calculation    *
  ! recoutnum          number of samples for each receptor calculation         *
  ! recoutnumsat       number of samples for each satellite calculation        *
  ! prob               probability of absorption at ground due to dry          *
  !                    deposition                                              *
  ! WETDEP             .true. if wet deposition is switched on                 *
  ! weight             weight for each concentration sample (1/2 or 1)         *
  !                                                                            *
  !*****************************************************************************
  ! openmp change
  use omp_lib
  ! openmp change end
  use unc_mod
  use point_mod
  use xmass_mod
  use flux_mod
  use outgrid_mod
!  use ohr_mod
  use par_mod
  use com_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use particle_mod
  use conv_mod
  use windfields_mod
  use advance_mod, only: advance
  use drydepo_mod
  use wetdepo_mod
  use plume_mod
  use initialise_mod
  use getfields_mod
  use output_mod
  use restart_mod
  use interpol_mod, only: alloc_interpol,dealloc_interpol
#ifdef USE_NCF
  use chemistry_mod
  use initdomain_mod
  use receptor_netcdf_mod, only: readreceptors_satellite, verttransform_satellite
  use emissions_mod
  use totals_mod
#endif
  use receptor_mod, only: receptoroutput

  implicit none
  real, parameter ::        &
    e_inv = 1.0/exp(1.0)  
  integer ::                &
    j,i,                    & ! loop variables
    iterminate,             & ! Keep track of terminated particles per timestep
    ks,                     & ! loop variable species
    kp,                     & ! loop variable for maxpointspec_act
    itime=0,                & ! time index
    nstop1,                 & ! windfield existence flag
    loutnext,               & ! following timestep
    loutstart,loutend,      & ! concentration calculation starting and ending time
    lrecoutnext,            & ! following timestep for receptor output
    lrecoutstart,lrecoutend,& ! receptor calculation interval start and end time
    ldeltat,                & ! radioactive decay time
    itage,nage,inage,       & ! related to age classes
    i_nan=0,ii_nan,total_nan_intl=0, &  !added by mc to check instability in CBL scheme 
    stat,                   & ! Check if allocation was successful
    thread                    ! openmp thread number
  ! logical ::                &
  !   active_per_rel(maxpoint)  ! are there particles active in each release
  real ::                   &
    filesize!(maxpoint)        ! Keeping track of the size of the particledump output, so it can be splitted
  ! real(kind=dp) ::          &
  !   jul
  ! integer ::                &
  !   jjjjmmdd,ihmmss
  real ::                   &
    outnum,                 & ! number of samples for grid concentration calculation
    prob_rec(maxspec),      & ! dry deposition related
    xmassfract                ! dry deposition related
  real(dep_prec),allocatable,dimension(:) ::         &
    drytmp       ! dry deposition related
  logical :: itsopen
  real, dimension(maxrecsample) :: recoutnum ! number of samples for receptor calculation
  real, allocatable, dimension(:,:) :: recoutnumsat ! number of samples for satellite receptor calculation

  ! First output for time 0
  !************************
  if (itime_init.ne.0) then
    loutnext=loutnext_init
    lrecoutnext=lrecoutnext_init
    outnum=outnum_init
  else
    loutnext=loutstep/2
    lrecoutnext=lrecoutstep/2
    outnum=0.
    recoutnum(:)=0.
!    recoutnumsat(:,:)=0.
  endif
  loutstart=loutnext-loutaver/2
  loutend=loutnext+loutaver/2
  lrecoutstart=lrecoutnext-lrecoutaver/2
  lrecoutend=lrecoutnext+lrecoutaver/2

  ! Initialise the nan count for CBL option
  !****************************************
  sum_nan_count(:) = 0
  nan_count(:) = 0

  !**********************************************************************
  ! Loop over the whole modelling period in time steps of mintime seconds
  !**********************************************************************

  write(*,46) float(itime)/3600,itime,numpart
46      format(' Simulated ',f7.1,' hours (',i13,' s), ',i13, ' particles')

  filesize=0.
  ! active_per_rel=.false.
  
  do itime=itime_init,ideltas,lsynctime

  ! Computation of wet deposition, OH reaction and mass transfer
  ! between two species every lsynctime seconds
  ! maybe wet depo frequency can be relaxed later but better be on safe side
  ! wetdepo must be called BEFORE new fields are read in but should not
  ! be called in the very beginning before any fields are loaded, or
  ! before particles are in the system
  ! Code may not be correct for decay of deposition
  ! changed by Petra Seibert 9/02
  !********************************************************************

  ! Write basic information on the simulation to a file "header" for the
  ! first time step and open files that are to be kept open throughout 
  ! the simulation.
  ! In addition, open new particle dump files if required and keep track
  ! of the size of these files.
  !*********************************************************************
    
    write(*,*) 'Time: ', itime, 'seconds. Total spawned:',count%spawned, 'alive:',count%alive, 'terminated:',count%terminated

    if (itime.eq.itime_init) then
      call SYSTEM_CLOCK(count_clock, count_rate, count_max)
      s_firstt = real(count_clock)/real(count_rate)
    endif

  ! Writing restart file
  !*********************
    if ((itime.ne.itime_init).and.(loutrestart.ne.-1).and.(mod(itime,loutrestart).eq.0)) then
      call output_restart(itime,loutnext,lrecoutnext,outnum)
    endif

    call init_output(itime,filesize)

  ! Get necessary wind fields if not available
  !*******************************************
    call getfields(itime,nstop1) !OMP on verttransform_ecmwf and readwind_ecmwf, getfields_mod.f90
    if (nstop1.gt.1) error stop 'NO METEO FIELDS AVAILABLE'

  ! In case of ETA coordinates being read from file, convert the z positions to zeta
  !*********************************************************************************
#ifdef ETA
    if ((itime.eq.itime_init).and.((ipin.eq.1).or.(ipin.eq.3).or.(ipin.eq.4))) then 
      
      if (count%allocated.le.0) error stop 'Something is going wrong reading the old particle file! &
        &No particles found.'
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      do i=1,count%allocated ! Also includes particles that are not spawned yet
        ! If kindz>1, vertical positions computation
        if (ipin.eq.3 .or. ipin.eq.4) call kindz_to_z(i) 
#ifdef ETA
        call z_to_zeta(itime,part(i)%xlon,part(i)%ylat,part(i)%z,part(i)%zeta)
        part(i)%etaupdate = .true.
        part(i)%meterupdate = .true.
#endif
      end do
!$OMP END DO
!$OMP END PARALLEL
    endif
#endif

    if (WETDEP .and. (itime.ne.0) .and. (numpart.gt.0)) then
      call wetdepo(itime,lsynctime,loutnext)
    endif

  ! compute chemical losses 
  !************************
#ifdef USE_NCF
    if (CLREA .and. (itime.ne.0) .and. (numpart.gt.0)) then
      call chemreaction(itime)
    endif
#endif

  ! compute convection for backward runs
  !*************************************

    if ((ldirect.eq.-1).and.(lconvection.eq.1).and.(itime.lt.0)) then
      call convmix(itime)
    endif

  ! Get chemical fields if not available 
  !*************************************
#ifdef USE_NCF
    if (CLREA) then
      call getchemfield(itime)
      call getchemhourly(itime)
    endif
#endif
        
  ! Get emission fields if not available
  !*************************************

#ifdef USE_NCF
    if (LEMIS.and.mdomainfill.eq.1) then
      call getemissions(itime)
    endif
#endif

  ! Read satellite receptors
  !*************************

#ifdef USE_NCF
  call readreceptors_satellite(itime)
#endif
  if (.not.allocated(recoutnumsat)) then
    allocate(recoutnumsat(nlayermax,maxrecsample))
    recoutnumsat(:,:)=0.
  endif

  ! Release particles
  !******************
    if (mdomainfill.ge.1) then
      if (itime.eq.itime_init) then 
        if (llcmoutput) then  
#ifdef USE_NCF
          call init_domainfill_ncf
#else
          call init_domainfill
#endif
        else
          call init_domainfill
        endif
        if (ipin.eq.2) then
          ! Particles initialized from partoutput
#ifdef ETA
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
          do i=1,count%alive
            j=count%ialive(i)
            call update_z_to_zeta(itime,j)
          end do
!$OMP END DO
!$OMP END PARALLEL
#endif
        endif
      else 
        call boundcond_domainfill(itime,loutend)
      endif
    else if ((ipin.eq.3).or.(ipin.eq.4)) then
      ! If reading from user defined initial conditions, check which particles are 
      ! to be activated
      if (count%allocated.le.0) error stop 'Something is going wrong reading the part_ic.nc file!'

      do i=1,count%allocated
        if (.not. part(i)%alive) then
          if (ldirect.lt.0) then
            if ((part(i)%tstart.le.itime).and.(part(i)%tstart.gt.itime+lsynctime)) then
              call spawn_particle(itime,i)
              call init_mass_conversion(i,part(i)%npoint)
            endif
          else if ((part(i)%tstart.ge.itime).and.(part(i)%tstart.lt.itime+lsynctime)) then
            call spawn_particle(itime,i)
            call init_mass_conversion(i,part(i)%npoint)
          endif
        endif
      end do

#ifdef ETA
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
      do i=1,count%alive
        j=count%ialive(i)
        if (part(j)%tstart.eq.itime) then
          call update_z_to_zeta(itime,j)
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
#endif
      call get_totalpart_num(numpart)
    else
      call releaseparticles(itime)
    endif

  ! Inject emissions
  !*****************

#ifdef USE_NCF
    if (LEMIS.and.mdomainfill.eq.1) then
      call emissions(itime)
    endif
#endif

  ! Compute convective mixing for forward runs
  ! for backward runs it is done before next windfield is read in
  !**************************************************************
    if ((ldirect.eq.1).and.(lconvection.eq.1)) then
      call convmix(itime) !OMP (not the nested part yet), conv_mod.f90
    endif

  ! If middle of averaging period of output fields is reached, accumulated
  ! deposited mass radioactively decays
  !***********************************************************************
    if (DEP.and.(itime.eq.loutnext).and.(ldirect.gt.0)) call deposit_decay() !OMP, unc_mod.f90 (needs test)


  ! Is the time within the computation interval, if not, skip
  !************************************************************
    if ((ldirect*itime.ge.ldirect*loutstart).and.(ldirect*itime.le.ldirect*loutend)) then
      call SYSTEM_CLOCK(count_clock, count_rate, count_max)
      s_temp = (count_clock - count_clock0)/real(count_rate)
      ! If it is not time yet to write outputs, skip
      !***********************************************
      if ((itime.eq.loutend).and.(outnum.gt.0).and.(itime.ne.0)) then

        if ((iout.eq.4).or.(iout.eq.5)) call plumetraj(itime)
        if (iflux.eq.1) call fluxoutput(itime)
        if (ipout.eq.1) then
          if (mod(itime,ipoutfac*loutstep).eq.0) then

            call output_particles(itime)!,active_per_rel) ! dump particle positions
          endif
        endif
      endif
      ! Check whether concentrations are to be calculated and outputted
      !****************************************************************
      call output_conc(itime,loutstart,loutend,loutnext,outnum)
      call SYSTEM_CLOCK(count_clock, count_rate, count_max)
      s_writepart = s_writepart + ((count_clock - count_clock0)/real(count_rate)-s_temp)
    endif

  ! Check whether receptor concentrations are to be calculated
  !***********************************************************

    if ((ldirect*itime.ge.ldirect*lrecoutstart).and. &
          ((numreceptor.gt.0.).or.(numsatreceptor.gt.0)).and. &
          (ldirect*itime.le.ldirect*lrecoutend)) then
      call receptoroutput(itime,lrecoutstart,lrecoutend,lrecoutnext,recoutnum,recoutnumsat)
    endif

    if (itime.eq.ideltas) exit         ! almost finished

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

    if (itime.lt.loutnext) then
      ldeltat=itime-(loutnext-loutstep)
    else                                  ! first half of next interval
      ldeltat=itime-loutnext
    endif


  ! Loop over all particles
  !************************
  ! Various variables for testing reason of CBL scheme, by mc
    well_mixed_vector=0. !erase vector to test well mixed condition: modified by mc
    well_mixed_norm=0.   !erase normalization to test well mixed condition: modified by mc
    avg_ol=0.
    avg_wst=0.
    avg_h=0.
    avg_air_dens=0.  !erase vector to obtain air density at particle positions: modified by mc
  !-----------------------------------------------------------------------------

  ! openmp change
  ! LB, openmp following CTM version
!$OMP PARALLEL PRIVATE(prob_rec,inage,nage,itage,ks,kp,thread,j,i,xmassfract)

#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM() ! Starts with 0
#else
    thread = 0
#endif

!$OMP DO SCHEDULE(dynamic,max(1,numpart/1000))
! SCHEDULE(dynamic, max(1,numpart/1000))
!max(1,int(real(numpart)/numthreads/20.)))
    do i=1,count%alive

      j=count%ialive(i)

  ! If integration step is due, do it
  !**********************************
  ! Determine age class of the particle
  !************************************
      itage=abs(itime-part(j)%tstart)
      nage=1
      if (lagespectra.eq.1) then
        do inage=1,nageclass
          nage=inage
          if (itage.lt.lage(nage)) exit
        end do
      endif

  ! Initialize newly released particle
  !***********************************
      if ((part(j)%tstart.eq.itime).or.(itime.eq.0)) then
#ifdef ETA
        call update_zeta_to_z(itime, j)
#endif
        call init_particle(itime,j,thread)
      endif

  ! Memorize particle positions
  !****************************
      part(j)%xlon_prev=part(j)%xlon
      part(j)%ylat_prev=part(j)%ylat
      part(j)%z_prev=part(j)%z
#ifdef ETA
      part(j)%zeta_prev=part(j)%zeta
#endif

  ! RECEPTOR: dry/wet depovel
  !****************************
  ! Before the particle is moved 
  ! the calculation of the scavenged mass shall only be done once after release
  ! xscav_frac1 was initialised with a negative value

      if  (DRYBKDEP) then
        do ks=1,nspec
          if  ((xscav_frac1(j,ks).lt.0)) then
#ifdef ETA
            call update_zeta_to_z(itime,j)
#endif
            call get_vdep_prob(itime,real(part(j)%xlon),real(part(j)%ylat), &
              real(part(j)%z),prob_rec,thread+1)
            if (DRYDEPSPEC(ks)) then        ! dry deposition
              xscav_frac1(j,ks)=prob_rec(ks)
            else
              mass(j,ks)=0.
              xscav_frac1(j,ks)=0.
            endif
          endif
        enddo
      endif

  ! Integrate Langevin equation for lsynctime seconds
  !*************************************************

      call advance(itime,j,thread)

      if (part(j)%nstop.eqv..true.) cycle
      if (n_average.gt.0) call partpos_avg(itime,j)

  ! Calculate the gross fluxes across layer interfaces
  !***************************************************
      if (iflux.eq.1) call calcfluxes(itime,nage,j,real(part(j)%xlon_prev), &
        real(part(j)%ylat_prev),real(part(j)%z_prev),thread+1)
    end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef _OPENMP
  call omp_set_num_threads(numthreads_grid)
#endif

  ! Terminating particles flagged in advance call
  iterminate=0
  do i=1,count%allocated
    if ((part(i)%nstop).and.(part(i)%alive)) then
      call terminate_particle(i,itime)
      iterminate=iterminate+1
    endif
  end do
  if (iterminate.gt.0) call rewrite_ialive()

  if (DRYDEP.or.WETDEP.or.LDECAY.or.(lagespectra.eq.1)) then
!$OMP PARALLEL PRIVATE(prob_rec,nage,inage,itage,ks,kp,thread,i,j,xmassfract,drytmp)

!num_threads(numthreads_grid)

#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM() ! Starts with 0
#else
    thread = 0
#endif
  if (DRYDEP) then
    allocate( drytmp(maxspec),stat=stat )
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drytmp inside of OMP loop'
  endif

!$OMP DO SCHEDULE(static)
!, max(1,numpart/1000))
! SCHEDULE(dynamic, max(1,numpart/1000))
!max(1,int(real(numpart)/numthreads/20.)))
    do i=1,count%alive

      j=count%ialive(i)

  ! If integration step is due, do it
  !**********************************
  ! Determine age class of the particle
  !************************************
      itage=abs(itime-part(j)%tstart)
      nage=1
      if (lagespectra.eq.1) then
        do inage=1,nageclass
          nage=inage
          if (itage.lt.lage(nage)) exit
        end do
      endif

! Dry deposition and radioactive decay for each species
! Also check maximum (of all species) of initial mass remaining on the particle;
! if it is below a threshold value, terminate particle
!*****************************************************************************
      if (DRYDEP.or.WETDEP.or.LDECAY) then
        xmassfract=0.
        do ks=1,nspec

          if (DRYDEPSPEC(ks)) then     ! dry deposition (and radioactive decay)

            call drydepo_massloss(j,ks,ldeltat,drytmp(ks))

          else if (decay(ks).gt.0.) then  ! no dry depo, but radioactive decay

            mass(j,ks) = mass(j,ks) * &
              exp( -real(abs(lsynctime)) * decay(ks) )

          endif

    ! Skip check on mass fraction when npoint represents particle number
          if (mdomainfill.eq.0.and.mquasilag.eq.0) then
            if (ipin.eq.3 .or. ipin.eq.4) then 
              if (mass_init(j,ks).gt.0) &
                xmassfract = max( xmassfract, &
                                  mass(j,ks) / mass_init(j,ks) )
            else if (xmass(part(j)%npoint,ks).gt.0.) then
              xmassfract = max( xmassfract, real( npart(part(j)%npoint) ) * &
                mass(j,ks) /  xmass(part(j)%npoint,ks) )
            endif
          else
            xmassfract=1.0
          end if

        end do
        
        if (xmassfract.le.minmassfrac) then 
          ! flag all particles carrying less mass for termination after parallel region
          part(j)%nstop=.true.
        endif
      endif
!        Sabine Eckhardt, June 2008
!        don't create depofield for backward runs
      if (DRYDEP.AND.(ldirect.eq.1).and.(iout.ne.0)) then

        if (ioutputforeachrelease.eq.1) then
            kp=part(j)%npoint
        else
            kp=1
        endif

        call drydepokernel(part(j)%nclass,drytmp,real(part(j)%xlon), &
             real(part(j)%ylat),nage,kp,thread+1)
        if (nested_output.eq.1) call drydepokernel_nest( &
             part(j)%nclass,drytmp,real(part(j)%xlon),real(part(j)%ylat), &
             nage,kp,thread+1)
      endif

  ! Terminate trajectories that are older than maximum allowed age
  !***************************************************************
      if (lagespectra.eq.1) then
        if ((part(j)%alive).and.(abs(itime-part(j)%tstart).ge.lage(nageclass))) then
          if (linit_cond.ge.1) call initcond_calc(itime+lsynctime,j,thread+1)
          ! flag all particles for termination after parallel region
          part(j)%nstop=.true.
        endif
      endif
    end do !loop over particles

!$OMP END DO
    if (DRYDEP) deallocate(drytmp)
!$OMP END PARALLEL
  
    ! Terminating particles flagged due to insufficient mass or exceeded max age
    iterminate=0
    do i=1,count%allocated
      if ((part(i)%nstop).and.(part(i)%alive)) then
        call terminate_particle(i,itime)
        iterminate=iterminate+1
      endif
    end do
    if (iterminate.gt.0) call rewrite_ialive()

  endif !(DRYDEP.or.WETDEP.or.LDECAY.or.(lagespectra.eq.1))

#ifdef _OPENMP
  call omp_set_num_threads(numthreads)
#endif
  ! OpenMP Reduction for dynamically allocated arrays. This is done manually since this
  ! is not yet supported in most OpenMP versions
  !************************************************************************************
#ifdef _OPENMP
    if (iflux.eq.1) then
      do i=1,numthreads
        flux(:,:,:,:,:,:,:)=flux(:,:,:,:,:,:,:)+flux_omp(:,:,:,:,:,:,:,i)
        flux_omp(:,:,:,:,:,:,:,i)=0.
      end do
    endif
    if (linit_cond.ge.1) then
      do i=1,numthreads_grid
        init_cond(:,:,:,:,:)=init_cond(:,:,:,:,:)+init_cond_omp(:,:,:,:,:,i)
        init_cond_omp(:,:,:,:,:,i)=0.
      end do
    endif
    if (DRYDEP.AND.(ldirect.eq.1).and.(iout.ne.0)) then
      do i=1,numthreads_grid
        drygridunc(:,:,:,:,:,:)=drygridunc(:,:,:,:,:,:)+gridunc_omp(:,:,1,:,:,:,:,i)
        gridunc_omp(:,:,1,:,:,:,:,i)=0.
      end do
      if (nested_output.eq.1) then
        do i=1,numthreads_grid
          drygriduncn(:,:,:,:,:,:)=drygriduncn(:,:,:,:,:,:)+griduncn_omp(:,:,1,:,:,:,:,i)
          griduncn_omp(:,:,1,:,:,:,:,i)=0.
        end do
      endif
    endif
#endif
  ! write(*,*) 'DRYGRIDUNC:',sum(drygridunc),drygridunc(20,270,1,1,1,1),drygridunc(19,269,1,1,1,1)
  ! Counter of "unstable" particle velocity during a time scale of
  ! maximumtl=20 minutes (defined in com_mod)
  !***************************************************************
    
    total_nan_intl=0
    i_nan=i_nan+1 ! added by mc to count nan during a time of maxtl (i.e. maximum tl fixed here to 20 minutes, see com_mod)
    do i=1,numthreads
      sum_nan_count(i_nan)=sum_nan_count(i_nan)+nan_count(i)
    end do
    if (i_nan > maxtl/lsynctime) i_nan=1 !lsynctime must be <= maxtl
    do ii_nan=1, (maxtl/lsynctime) 
      total_nan_intl=total_nan_intl+sum_nan_count(ii_nan)
    end do
  ! Output to keep track of the numerical instabilities in CBL simulation and if
  ! they are compromising the final result (or not)
    if (cblflag.eq.1) print *,j,itime,'nan_synctime',sum_nan_count(i_nan),'nan_tl',total_nan_intl  

    if (itime.eq.itime_init) then
      call SYSTEM_CLOCK(count_clock, count_rate, count_max)
      s_firstt = real(count_clock)/real(count_rate) - s_firstt
    endif

    ! Output totals
    !**************

#ifdef USE_NCF
    if ((mdomainfill.eq.1).and.(llcmoutput)) then
      do ks=1,nspec
        tot_mass(ks)=sum(real(mass(1:count%alive,ks),kind=dp))
      end do
      call totals_write(itime)
    endif
#endif

  end do

  ! Complete the calculation of initial conditions for particles not yet terminated
  !*****************************************************************************
  call finalise_output(itime)

  ! Output residual emissions
  !**************************

#ifdef USE_NCF
  if (LEMIS.and.(ipout.eq.2)) then
    call em_res_write
  endif
#endif

  ! De-allocate memory and end
  !***************************
  call dealloc_all_particles
  call dealloc_windfields
  call dealloc_domainfill
  call dealloc_drydepo
  call dealloc_convect
  call dealloc_getfields
  call dealloc_interpol
  call dealloc_random
  if (numbnests.ge.1) call dealloc_windfields_nest
  if (iflux.eq.1) deallocate(flux)
  if (ipin.ne.3 .and. ipin.ne.4) deallocate(xmasssave)
  if (CLREA) then
    deallocate(CL_field,lonCL,latCL,altCL)
  endif
  deallocate(reaccconst,reacdconst,reacnconst)
  deallocate(emis_path,emis_file,emis_name,emis_unit,emis_coeff)
  if (lnetcdfout.eq.1) then
#ifdef USE_NCF
    call dealloc_netcdf
    if (LEMIS) then
      deallocate(em_field,em_res,em_area,mass_field)
    endif
#endif 
  else
    inquire(unit=unitoutrecept, opened=itsopen)
    if (itsopen) close(unitoutrecept)
    inquire(unit=unitoutreceptppt, opened=itsopen)
    if (itsopen) close(unitoutreceptppt)
    inquire(unit=unitoutsatellite, opened=itsopen)
    if (itsopen) close(unitoutsatellite)
  endif
  deallocate(xpoint1,xpoint2,ypoint1,ypoint2,zpoint1,zpoint2)
  deallocate(xmass)
  deallocate(ireleasestart,ireleaseend,npart,kindz)
  deallocate(nan_count)
  if (ipout.ne.0) deallocate( partopt )
  if (iout.ne.0) then
    deallocate(outheight,outheighthalf)
    deallocate(oroout, area, volume)
    deallocate(gridunc)
    deallocate(gridcnt)
#ifdef _OPENMP
    deallocate(gridunc_omp)
    deallocate(gridcnt_omp)
#endif
    if (ldirect.gt.0) then
      deallocate(drygridunc,wetgridunc)
#ifdef _OPENMP
      deallocate(drygridunc_omp,wetgridunc_omp)
#endif
    endif
    if (nested_output.eq.1) then
      deallocate(orooutn, arean, volumen)
      if (ldirect.gt.0) then
        deallocate(griduncn,drygriduncn,wetgriduncn)
#ifdef _OPENMP
        deallocate(griduncn_omp,drygriduncn_omp,wetgriduncn_omp)
#endif
      endif
    endif
  endif
end subroutine timemanager

end module timemanager_mod
