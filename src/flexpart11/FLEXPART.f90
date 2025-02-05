! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

program flexpart

  !*****************************************************************************
  !                                                                            *
  !     This is the Lagrangian Particle Dispersion Model FLEXPART.             *
  !     The main program manages the reading of model run specifications, etc. *
  !     All actual computing is done within subroutine timemanager.            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Changes:                                                                   *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     (moved to read_options_and_initialise_flexpart by LB)                  *
  !     - Added detection of metdata format using gributils routines           *
  !     - Distinguished calls to ecmwf/gfs gridcheck versions based on         *
  !       detected metdata format                                              *
  !     - Passed metdata format down to timemanager                            *
  !   L. Bakels 2022                                                           *
  !     - OpenMP parallelisation                                               *
  !     - Added input options                                                  *
  !     - Restructuring into subroutines (below)                               *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use par_mod
  use com_mod
  use timemanager_mod
  use output_mod

  implicit none

  integer :: i
  real :: s_timemanager
  character(len=256) ::   &
    inline_options          ! pathfile, flexversion, arg2
  character(len=256) :: gitversion_tmp="undefined"

  ! Keeping track of the total running time of FLEXPART, printed out at the end.
  !*****************************************************************************
  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_total = (count_clock - count_clock0)/real(count_rate)

  ! FLEXPART version string
  flexversion_major = '11' ! Major version number, also used for species file names

  flexversion='Version '//trim(flexversion_major)//'.0 (2023-07-11)'
  verbosity=0

  call update_gitversion(gitversion_tmp)
  ! Read the pathnames where input/output files are stored
  !*******************************************************

  inline_options='none'
  select case (iargc())
  case (2)
    call getarg(1,arg1)
    pathfile=arg1
    call getarg(2,arg2)
    inline_options=arg2
  case (1)
    call getarg(1,arg1)
    pathfile=arg1
    if (arg1(1:1).eq.'-') then
      write(pathfile,'(a11)') './pathnames'
      inline_options=arg1 
    endif
  case (0)
    write(pathfile,'(a11)') './pathnames'
  end select

  ! Print the GPL License statement
  !******************************************************
  print*,'Welcome to FLEXPART ', trim(flexversion)
  print*,"Git: ", trim(gitversion)
  print*,'FLEXPART is free software released under the GNU General Public License.'
#ifdef ETA
  write(*,*) 'FLEXPART is running with ETA coordinates.'
#else
  write(*,*) 'FLEXPART is running with METER coordinates.'
#endif

  ! Reading user specified options, allocating fields and checking bounds
  !**********************************************************************
  call read_options_and_initialise_flexpart

  ! Inform whether output kernel is used or not
  !*********************************************
  if (lroot) then
    if (.not.lusekerneloutput) then
      write(*,*) "Concentrations are calculated without using kernel"
    else
      write(*,*) "Concentrations are calculated using kernel"
    end if
  end if

  if (lturbulence.eq.0) write(*,*) 'WARNING: turbulence switched off.'
  if (lmesoscale_turb) write(*,*) 'WARNING: mesoscale turbulence switched on.'
  if (log_interpol) write(*,*) 'WARNING: using logarithmical vertical interpolation.'
  ! Calculate particle trajectories
  !********************************
  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_timemanager = (count_clock - count_clock0)/real(count_rate)

  call timemanager

  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_timemanager = (count_clock - count_clock0)/real(count_rate) - s_timemanager

  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_total = (count_clock - count_clock0)/real(count_rate) - s_total
  
  if (verbosity.gt.0) then
! NIK 16.02.2005 
    do i=1,nspec
      if (icnt_incld(i).gt.0) then
         write(*,*) '**********************************************'
         write(*,*) 'Scavenging statistics for species ', species(i), ':'
         write(*,*) 'Total number of occurences of below-cloud scavenging', &
           & icnt_belowcld(i)
         write(*,*) 'Total number of occurences of in-cloud    scavenging', &
           & icnt_incld(i)
         write(*,*) '**********************************************'
      endif
    end do
  endif
  
  write(*,*) 'Read wind fields: ', s_readwind, ' seconds'
  write(*,*) 'Timemanager: ', s_timemanager, ' seconds,', 'first timestep: ',s_firstt, 'seconds'
  write(*,*) 'Write particle files: ', s_writepart, ' seconds'
  write(*,*) 'Total running time: ', s_total, ' seconds'
  write(*,*) 'tps,io,tot: ', (s_timemanager-s_firstt)/4.,(s_readwind+s_writepart)/5.,s_total
  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLE&
       &XPART MODEL RUN!'

end program flexpart


subroutine read_options_and_initialise_flexpart

  !*****************************************************************************
  !                                                                            *
  !   Moved from main flexpart program:                                        *
  !   Reading all option files, initialisation of random numbers, and          *
  !   allocating memory for windfields, grids, etc.                            *
  !                                                                            *
  !   L. Bakels 2022                                                           *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use random_mod
  use par_mod
  use com_mod
  use conv_mod
  use class_gribfile_mod
  use readoptions_mod
  use windfields_mod
  use plume_mod
  use initialise_mod
  use drydepo_mod
  use getfields_mod
  use interpol_mod,         only: alloc_interpol
  use outgrid_mod
  use binary_output_mod
  use omp_lib,              only: OMP_GET_MAX_THREADS
#ifdef USE_NCF
  use chemistry_mod,        only: readreagents
  use totals_mod
  use receptor_netcdf_mod,  only: read_satellite_info, receptorout_init
#endif
  use receptor_mod,         only: alloc_receptor

  implicit none

  integer ::              &
    inest                   ! loop variable for nested gridcells
  integer ::              &
    j,                    & ! loop variable for random numbers
    stat,                 & ! Check if allocation was successful
    idummy=-320             ! dummy value used by the random routine


  ! Read pathnames from file in working director that specify I/O directories
  !**************************************************************************
  call readpaths
  
  ! Read the user specifications for the current model run
  !*******************************************************
  call readcommand

  ! Reading the number of threads available and print them for user
  !****************************************************************
#ifdef _OPENMP
    numthreads = OMP_GET_MAX_THREADS()
    numthreads_grid = min(numthreads,maxthreadgrid)
    !numthreads = min(40,numthreads)
#else
    numthreads = 1
    numthreads_grid = 1
#endif

  if (numthreads.gt.1) then
    write(*,*)
    write(*,*) "*********** WARNING  **********************************"
    write(*,*) "* FLEXPART running in parallel mode                   *"
    write(*,*) "* Number of uncertainty classes in                    *"
    write(*,901) " * set to number of threads:            ", &
      numthreads_grid, "          *" 
    write(*,901) " * All other computations are done with ",&
      numthreads, " threads. *"
    write(*,*) "*******************************************************"
    write(*,*)
901 format (a,i5,a)
  endif

  call alloc_random(numthreads)

  ! Generate a large number of random numbers
  !******************************************
  do j=1,maxrand-1,2
    call gasdev1(idummy,rannumb(j),rannumb(j+1))
  end do
  call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))

  ! Read the age classes to be used
  !********************************
  call readageclasses

  ! ! Allocate memory for windfields
  ! !*******************************
  ! call alloc_windfields
  ! if (numbnests.ge.1) then
  !   ! If nested wind fields are used, allocate arrays
  !   !************************************************
  !   call alloc_windfields_nest
  ! endif

  ! Read, which wind fields are available within the modelling period
  !******************************************************************
  call readavailable

  ! Detect metdata format
  !**********************
  call detectformat

  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    print *,'ECMWF metdata detected'
    if (nxshift.eq.-9999) nxshift=359
  elseif (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
    print *,'NCEP metdata detected'
    if (nxshift.eq.-9999) nxshift=0
  else
    error stop 'Unknown metdata format'
  endif
  write(*,*) 'NXSHIFT is set to', nxshift

  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    call gridcheck_ecmwf
  else 
    call gridcheck_gfs
  endif

  ! Set the upper level for where the convection will be working
  !*************************************************************
  call set_conv_top

  if (numbnests.ge.1) then
  ! If nested wind fields are used, allocate arrays
  !************************************************
    call alloc_nest_properties
    call gridcheck_nest
  endif

  ! Read the output grid specifications if requested by user
  !*********************************************************
  if (iout.ne.0) then
    call readoutgrid

    if (nested_output.eq.1) then
      call readoutgrid_nest
    endif
  endif

  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************
  numreceptor=0
  numsatreceptor=0
  nlayermax=1
#ifdef USE_NCF
  call read_satellite_info
#endif
  call readreceptors

  ! Read the landuse inventory
  !***************************
  call readlanduse ! CHECK ETA

  ! Read chemical reagent information
  !**********************************
  ! default settings
  nreagent=0
  reagents(:)=""
#ifdef USE_NCF
  call readreagents
#endif 

  ! For continuation of previous run or from user defined initial 
  ! conditions, read in particle positions
  !*************************************************************************
  call initialise_particles

  ! Initialize variables for totals calculations
  !*********************************************
#ifdef USE_NCF
  call alloc_totals
  call totals_init
#endif 

  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************
  call coordtrafo(nxmin1,nymin1)

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************
  call readdepo

  ! Allocate dry deposition fields if necessary
  !*********************************************
  call alloc_drydepo
  call alloc_convect
  call alloc_getfields
  call alloc_interpol
#ifdef USE_NCF
  if (lnetcdfout.eq.1) call alloc_netcdf
#endif 

  ! Assign fractional cover of landuse classes to each ECMWF grid point
  !********************************************************************
  call assignland

  ! Calculate volume, surface area, etc., of all output grid cells
  !***************************************************************
  if (iout.ne.0) then
    call outgrid_init
    if (nested_output.eq.1) call outgrid_init_nest ! CHECK ETA
  endif

  ! Initialize receptor output
  !***************************

  call alloc_receptor
  if (lnetcdfout.eq.1) then
#ifdef USE_NCF
    call receptorout_init
#endif  
  else
    call receptorout_init_binary
  endif

  if ((iout.eq.4).or.(iout.eq.5)) call openouttraj ! CHECK ETA


  ! Initialize cloud-base mass fluxes for the convection scheme
  !************************************************************
  if (lconvection.eq.1) then
    cbaseflux(0:nxmin1,0:nymin1,:)=0.
    do inest=1,numbnests
      cbasefluxn(0:nxn(inest)-1,0:nyn(inest)-1,inest,:)=0.
    end do
  endif

  ! Allocating nan_count for CBL option
  !************************************
  allocate(nan_count(numthreads), stat=stat)
  if (stat.ne.0) error stop "Could not allocate nan_count"

end subroutine read_options_and_initialise_flexpart

subroutine initialise_particles

  !*****************************************************************************
  !                                                                            *
  !   This subroutine handles the different forms of starting FLEXPART         *
  !   depending on IPIN (set in COMMAND)                                       *
  !                                                                            *
  !   IPIN=0: this routine is not called and particles are read from the       *
  !           RELEASE option file                                              *
  !   IPIN=1: restarting from a restart.bin file, written by a previous run    *
  !   IPIN=2: restarting from a partoutput_xxx.nc file written by a previous   *
  !           run, depending on what PARTOPTIONS the user has chosen, this     *
  !           option might not be possible to use                              *
  !   IPIN=3: starting a run from a user defined initial particle conditions,  *
  !           more on how to create such a file can be found in the manual     *
  !   IPIN=4: restarting a run, while also reading in the initial particle     *
  !           conditions                                                       *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use com_mod
  use initialise_mod
#ifdef USE_NCF  
  use netcdf_output_mod
#endif
  use readoptions_mod
  use restart_mod
  use settling_mod

  implicit none

  integer :: i
  ! logical :: l_lookup=.false.

  ! Read the coordinates of the release locations
  !**********************************************
  itime_init=0

  if (ipin.le.2) then 
    call readreleases
    ! needs to be called after maxspec is defined in readreleases or readinitconditions
    if (ipout.ne.0) call readpartoptions 
  else
#ifdef USE_NCF
    call readinitconditions_netcdf
#else
    error stop 'Compile with netCDF if you would like to use the ipin=3,4 options.'
#endif
  endif

  if (iout.ne.0) then
    call alloc_grid
    call alloc_grid_unc
    if (nested_output.eq.1) call alloc_grid_unc_nest
  endif
  
  if ((ipin.eq.1).or.(ipin.eq.4)) then ! Restarting from restart.bin file
    call readrestart
  else if (ipin.eq.2) then ! Restarting from netcdf partoutput file
#ifdef USE_NCF
    call readpartpositions
#else
    error stop 'Compile with netCDF if you want to use the ipin=2 option.'
#endif
  else if (ipin.eq.0) then
    ! Releases can only start and end at discrete times (multiples of lsynctime)
    !***************************************************************************
    do i=1,numpoint
      ireleasestart(i)=nint(real(ireleasestart(i))/real(lsynctime))*lsynctime
      ireleaseend(i)=nint(real(ireleaseend(i))/real(lsynctime))*lsynctime
    end do
    numpart=0
    numparticlecount=0
  endif

  ! ! Initialise look-up table for drag coefficients when necessary
  ! !**************************************************************
  ! if (lsettling) then
  !   do i=1,nspec
  !     if (ishape(i).eq.0) l_lookup=.true.
  !   end do
  !   if (l_lookup) call init_dragcoeff_lookup()
  ! endif

end subroutine initialise_particles
