!*******************************************************************************
!        Include file for particle diffusion model FLEXPART                    *
!        This file contains a global common block used by FLEXPART             *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        June 1996                                                             *
!                                                                              *
!        Update: 15 August 2013 IP                                             *
!        PS 19 Nov 2020: correct comment about lcw                             *
!        Anne Tipka, Petra Seibert 2021-02: implement new interpolation        *
!           for precipitation according to #295 using 2 additional fields      *
!                                                                              *
!*******************************************************************************

module com_mod

  use par_mod, only: dp, numpath, maxnests, maxndia, &
       numclass, maxcolumn, maxrand, numwfmem, numpf, &
       maxreagent, maxrecsample

  implicit none

  ! Partoptions derived type. This decides which fields are being computed and output
  !**********************************************************************************
  type :: particleoptions
    character(2) :: name
    character(28) :: long_name
    character(7) :: short_name
    logical :: print
    logical :: average=.false.
    integer :: i_average=0
    integer :: ncid
  end type particleoptions

  integer :: num_partopt=34
  integer :: n_average
  type(particleoptions),allocatable :: partopt(:)

  !****************************************************************
  ! Variables defining where FLEXPART input/output files are stored
  !****************************************************************

  character :: path(numpath+2*maxnests)*120
  integer :: length(numpath+2*maxnests)
  character(len=256) :: pathfile, flexversion, flexversion_major, arg1, arg2
  character(len=256) :: ohfields_path
  character(len=256) :: gitversion
  
  ! path                    path names needed for trajectory model
  ! length                  length of path names needed for trajectory model
  ! pathfile                file where pathnames are stored
  ! flexversion             version of flexpart (descriptive long string)
  ! flexversion_major       version of flexpart (major version number)
  ! arg                     input arguments from launch at command line
  ! ohfields_path           path to binary files for OH fields

  !********************************************************
  ! Variables defining the general model run specifications
  !********************************************************

  integer :: ibdate,ibtime,iedate,ietime,itime_init,loutnext_init,lrecoutnext_init
  real :: outnum_init
  real(kind=dp) :: bdate,edate


  ! ibdate                  beginning date (YYYYMMDD)
  ! ibtime                  beginning time (HHMISS)
  ! iedate                  ending date (YYYYMMDD)
  ! ietime                  ending time (HHMISS)
  ! itime_init              starting time in [s] in case of a restart
  ! bdate                   beginning date of simulation (julian date)
  ! edate                   ending date of simulation (julian date)
  ! outnum_init             concentration calculation sample number after restart
  ! loutnext_init           first writing time after restart


  integer :: ldirect,ideltas

  ! ldirect                 1 for forward, -1 for backward simulation
  ! ideltas                 length of trajectory loop from beginning to
  !                    ending date (s)

  integer :: loutstep,loutaver,loutsample,loutrestart,method,lsynctime
  integer :: lrecoutstep,lrecoutaver,lrecoutsample
  real :: outstep

  ! loutstep [s]            gridded concentration output every loutstep seconds
  ! loutaver [s]            concentration output is an average over [s] seconds
  ! loutsample [s]          sampling interval of gridded concentration output
  ! lrecoutstep [s]         receptor concentration output every loutstep seconds
  ! lrecoutaver [s]         receptor concentration output is an average over [s] seconds
  ! lrecoutsample [s]       sampling interval of receptor concentration output
  ! loutrestart [s]         time interval for writing restart files
  ! lsynctime [s]           synchronisation time of all particles
  ! method                  indicator which dispersion method is to be used
  ! outstep = real(abs(loutstep))

  real :: ctl,fine
  integer :: ifine,iout,ipout,ipin,iflux,mdomainfill,ipoutfac
  integer :: mquasilag,nested_output,ind_source,ind_receptor,nxshift
  integer :: ind_rel,ind_samp,ioutputforeachrelease,linit_cond,sfc_only
  integer :: surf_only ! deprecated
  logical :: turbswitch
  integer :: cblflag !added by mc for cbl
  logical :: llcmoutput

  ! ctl      factor, by which time step must be smaller than Lagrangian time scale
  ! ifine    reduction factor for time step used for vertical wind
  !     Langevin equation for the vertical wind component
  ! ioutputforeachrelease Should each release be a seperate output field?
  ! iflux    flux calculation options: 1 calculation of fluxes, 2 no fluxes
  ! iout     output options: 1 conc. output (ng/m3), 2 mixing ratio (pptv), 3 both
  ! ipout    particle dump options: 0 no, 1 every output interval, 2 only at end
  ! ipoutfac increase particle dump interval by factor (default 1)
  ! ipin     read in particle positions from dumped file from a previous run
  ! fine     real(ifine)
  ! mdomainfill 0: normal run
  !        1: particles are initialized according to atmospheric mass distribution
  ! ind_source switches between different units for concentrations at the source
  !  NOTE that in backward simulations the release of computational particles
  !  takes place at the "receptor" and the sampling of particles at the "source".
  !     1= mass units
  !     2= mass mixing ratio units
  ! ind_receptor switches between different units for FLEXPART concentration at the receptor
  !     1= mass units
  !     2= mass mixing ratio units
  ! linit_cond  switch on the output of sensitivity to initial conditions for backward runs
  !     0=no, 1=mass unit, 2=mass mixing ratio unit
  ! mquasilag 0: normal run
  !      1: Particle position output is produced in a condensed format and particles are numbered
  ! sfc_only   switch output in grid_time files for surface only or full vertical resolution
  !      0=no (full vertical resolution), 1=yes (surface only)
  ! nested_output: 0 no, 1 yes
  ! turbswitch              determines how the Markov chain is formulated
  ! nxshift            for global grids (in x), the grid can be shifted by
  !                    nxshift grid points, in order to accomodate nested
  !                    grids, and output grids overlapping the domain "boundary"
  !                    nxshift must not be negative; "normal" setting would be 0

  ! ind_rel and ind_samp  are used within the code to change between mass and mass-mix (see readcommand.f)
  ! cblflag !: 1 activate cbl skewed pdf routines with bi-gaussina pdf whan OL<0 added by mc
  ! llcmoutput  switch for LCM output (uses mass ratio of species to air tracer) or normal output


  integer :: mintime,itsplit

  ! mintime                 minimum time step to be used by FLEXPART
  ! itsplit                 time constant for splitting particles

  integer :: lsubgrid,lconvection,lturbulence,lagespectra

  ! lsubgrid     1 if subgrid topography parameterization switched on, 2 if not
  ! lconvection  1 if convection parameterization switched on, 0 if not
  ! lturbulence  1 if turbulence parameterization switched on, 0 if not
  ! lagespectra  1 if age spectra calculation switched on, 2 if not

  ! mesoscale turbulence is found to give issues, so turned off by default
  !***********************************************************************
  logical :: lmesoscale_turb=.false.

  integer :: lnetcdfout
  ! lnetcdfout   1 for netcdf grid output, 0 if not. Set in COMMAND (namelist input)

  integer :: linversionout
  ! linversionout 1 for one grid_time output file for each release containing all timesteps

  integer :: nageclass
  integer,allocatable,dimension(:) :: lage

  ! nageclass               number of ageclasses for the age spectra calculation
  ! lage [s]                ageclasses for the age spectra calculation

 !ESO: Disable settling if more than 1 species per release point
  logical :: lsettling=.true.

  logical :: gdomainfill
  ! gdomainfill             .T., if domain-filling is global, .F. if not

  logical :: lcw=.false. ! ZHG Sep 2015 ! AT renamed
  ! whether or not cloud water data found in GRIB, overwritten if CW is found

  logical :: lcwsum=.false. ! ESO Dec 2015 ! AT renamed
  ! whether or not both clwc and ciwc are present (if so they are summed)

  logical :: lprecint ! AT, PS 2021
  ! true if new interpolation using additional precip fields is used
    
  logical,dimension(maxnests) :: lcw_nest=.false.
  logical,dimension(maxnests) :: lcwsum_nest=.false.
  logical,dimension(maxnests) :: lprecintn
 
  !NIK 16.02.2015
  integer(selected_int_kind(16)),allocatable,dimension(:) :: icnt_belowcld, &
       &icnt_incld

  !*********************************************************************
  ! Variables defining the release locations, released species and their
  ! properties, etc.
  !*********************************************************************

  !change Sabine Eckhardt, only save the first 1000 identifier for releasepoints
  character :: compoint(1001)*45
  integer :: numpoint
  !sec, now dynamically allocated:
  ! ireleasestart(maxpoint),ireleaseend(maxpoint)
  !      real xpoint1(maxpoint),ypoint1(maxpoint)
  !real xpoint2(maxpoint),ypoint2(maxpoint)
  !real zpoint1(maxpoint),zpoint2(maxpoint)
  !integer*2 kindz(maxpoint)
  integer,allocatable,dimension(:) :: specnum
  !real xmass(maxpoint,maxspec)
  real,allocatable,dimension(:) :: decay
  real,allocatable,dimension(:) :: weta_gas,wetb_gas
  real,allocatable,dimension(:) :: crain_aero,csnow_aero
! NIK: 31.01.2013- parameters for in-cloud scavening
  real,allocatable,dimension(:) :: ccn_aero,in_aero
  real,allocatable,dimension(:) :: reldiff,henry,f0
  real,allocatable,dimension(:) :: density,dquer,dsigma
  integer,allocatable,dimension(:) :: ndia
  real,allocatable,dimension(:) :: vsetaver,cunningham,weightmolar
  real,allocatable,dimension(:,:) :: vset,schmi,fract
  real,allocatable,dimension(:,:) :: ri,rac
  real,allocatable,dimension(:,:,:) :: rcl,rgs,rlu
  real,allocatable,dimension(:) :: rm,dryvel
  ! Daria Tatsii: species shape properties
  real,allocatable,dimension(:) :: Fn,Fs ! Newton and Stokes' regime
  real,allocatable,dimension(:) :: ks1,ks2,kn2
  integer,allocatable,dimension(:) :: ishape,orient
  ! chemical reagent variables
  character(len=256) :: reag_path(maxreagent)
  character(len=16)  :: reagents(maxreagent), reag_unit(maxreagent)
  integer :: reag_hourly(maxreagent), nreagent  
  ! reaction rates
  real,allocatable,dimension(:,:) :: reaccconst,reacdconst,reacnconst
  ! emissions variables for LCM
  character(len=256),allocatable,dimension(:) :: emis_path,emis_file,emis_name
  integer,allocatable,dimension(:) :: emis_unit
  real,allocatable,dimension(:) :: emis_coeff

  real,allocatable,dimension(:,:) :: area_hour,point_hour
  real,allocatable,dimension(:,:) :: area_dow,point_dow

  !integer npart(maxpoint)
  integer :: nspec,maxpointspec_act
  character(len=10),allocatable,dimension(:) :: species


  ! compoint                comment, also "name" of each starting point
  ! numpoint                actual number of trajectory starting/ending points
  ! ireleasestart,ireleaseend [s] starting and ending time of each release
  ! xmass                   total mass emitted
  ! xpoint1,ypoint1         lower left coordinates of release area
  ! xpoint2,ypoint2         upper right coordinates of release area
  ! zpoint1,zpoint2         min./max. z-coordinates of release points
  ! kindz                   1: zpoint is in m agl, 2: zpoint is in m asl
  ! npart                   number of particles per release point
  ! nspec                   number of different species allowed for one release
  ! maxpointspec_act        number of releaspoints for which a different output shall be created
  ! species                 name of species
  ! decay                   decay constant of radionuclide

  ! WET DEPOSITION
  ! weta_gas, wetb_gas     parameters for below-cloud wet scavenging coefficients (gasses)
  ! crain_aero, csnow_aero parameters for below-cloud wet scavenging coefficients (aerosols)
  ! ccn_aero, cin_aero     parameters for in-cloud wet scavenging coefficients (aerosols)

  ! GAS DEPOSITION
  ! reldiff                 diffusivitiy of species relative to diff. of H2O
  ! henry [M/atm]           Henry constant
  ! f0                      reactivity relative to that of O3
  ! ri [s/m]                stomatal resistance
  ! rcl [s/m]               lower canopy resistance
  ! rgs [s/m]               ground resistance
  ! rlu [s/m]               leaf cuticular resistance
  ! rm [s/m]                mesophyll resistance
  ! dryvel [m/s]            constant dry deposition velocity

  ! PARTICLE DEPOSITION
  ! density [kg/m3]         density of particles
  ! dquer [m]               mean diameter of particles
  ! dsigma                  dsigma=10 or dsigma=0.1 means that 68% of the
  !                         mass are between 0.1*dquer and 10*dquer

  ! fract                   mass fraction of each diameter interval
  ! vset [m/s]              gravitational settling velocity in ni intervals
  ! cunningham              Cunningham slip correction (strictly valid only near surface)
  ! vsetaver [m/s]          average gravitational settling velocity
  ! schmi                   Schmidt number**2/3 of each diameter interval
  ! weightmolar [g/mol]     molecular weight

  ! TIME VARIATION OF EMISSION
  ! area_hour, point_hour   daily variation of emission strengths for area and point sources
  ! area_dow, point_dow     day-of-week variation of emission strengths for area and point sources

  !******************************************************************************
  ! Variables associated with the ECMWF meteorological input data ("wind fields")
  !******************************************************************************

  integer :: memtime(numwfmem),memind(3),lwindinterv ! eso: or memind(numwfmem) 

  ! memtime [s]             validation times of wind fields in memory
  ! memind                  pointer to wind field, in order to avoid shuffling
  !                         of wind fields
  ! lwindinterv [s]         Interval between wind fields currently in memory

  !********************************************************************
  ! Variables associated with the ECMWF input data (nested wind fields)
  !********************************************************************

  ! NOTE: all nested variables have the same name as the variables used
  ! for the mother domain, except with a 'n' appended at the end
  !********************************************************************

  integer :: numbnests, nxmaxn, nymaxn

  ! nxmax,nymax        maximum dimension of wind fields in x and y
  !                    direction, respectively
  ! numbnests    number of nested grids

  !******************************************************
  ! Variables defining the polar stereographic projection
  !******************************************************

  logical :: xglobal,sglobal,nglobal
  real :: switchnorthg,switchsouthg

  !xglobal             T for global fields, F for limited area fields
  !sglobal             T if domain extends towards south pole
  !nglobal             T if domain extends towards north pole
  !switchnorthg,switchsouthg   same as parameters switchnorth,
  !                    switchsouth, but in grid units

  real :: southpolemap(9),northpolemap(9)

  !southpolemap,northpolemap   define stereographic projections
  !                    at the two poles


  !******************
  ! Landuse inventory
  ! Sabine Eckhardt Dec 06: change to new landuse inventary - 11 classes, 1200 x 600 global
  !******************

  integer(kind=1) :: landinvent(1200,600,6)
  real :: z0(numclass)

! !$OMP THREADPRIVATE (z0)

  ! landinvent         landuse inventory (numclass=11 classes)
  ! z0                  roughness length for the landuse classes



  !**************************************************************************
  ! Variables characterizing the output grid and containing the model results
  !**************************************************************************

  integer :: numxgrid,numygrid,numzgrid
  real :: dxout,dyout,outlon0,outlat0,xoutshift,youtshift
  integer :: numxgridn,numygridn
  real :: dxoutn,dyoutn,outlon0n,outlat0n,xoutshiftn,youtshiftn
  !real outheight(maxzgrid),outheighthalf(maxzgrid)

  logical :: DEP,DRYDEP,WETDEP,CLREA,ASSSPEC,LDECAY,LEMIS
  logical,allocatable,dimension(:) :: DRYDEPSPEC,WETDEPSPEC
  logical :: DRYBKDEP,WETBKDEP

  ! numxgrid,numygrid       number of grid points in x,y-direction
  ! numxgridn,numygridn     number of grid points in x,y-direction for nested output grid
  ! numzgrid                number of vertical levels of output grid
  ! dxout,dyout             grid distance of output grid
  ! dxoutn,dyoutn           grid distance of nested output grid
  ! outlon0,outlat0         lower left corner of output grid
  ! outlon0n,outlat0n       lower left corner of nested output grid
  ! xoutshift,youtshift     xlon0-outlon0, ylat0-outlat0
  ! xoutshiftn,youtshiftn   xlon0-outlon0n, ylat0-outlat0n
  ! outheight [m]           upper levels of the output grid
  ! outheighthalf [m]       half (middle) levels of the output grid cells
  ! DEP                     .true., if either dry or wet depos. is switched on
  ! DRYDEP                  .true., if dry deposition is switched on
  ! DRYDEPSPEC              .true., if dry deposition is switched on for that species
  ! WETDEP                  .true., if wet deposition is switched on
  ! WETDEPSPEC              .true., if wet deposition is switched on for that species
  ! CLREA                   .true., if chemical reactions is switched on
  ! ASSSPEC                 .true., if there are two species asscoiated
  ! DRYBKDEP,WETBKDEP        .true., for bkwd runs, where mass deposited and source regions is calculated - either for dry or for wet deposition
  !                    (i.e. transfer of mass between these two occurs
  ! LEMIS                   .true., if particle mass should change due to surface fluxes

  !  if output for each releasepoint shall be created maxpointspec=number of releasepoints
  !  else maxpointspec is 1 -> moved to unc_mod
  !  the OUTGRID is moved to the module outgrid_mod
  !******************************************************************************

  ! gridunc,griduncn        uncertainty of outputted concentrations
  ! wetgridunc,wetgriduncn  uncertainty of accumulated wet deposited mass on output grid
  ! drygridunc,drygriduncn  uncertainty of accumulated dry deposited mass on output grid
  ! oroout,orooutn [m]      height of model topography at output grid
  ! area,arean [m2]         area of each grid cell
  ! volume,volumen [m3]     volume of each grid cell
  ! ... field names with n at the end indicate a nested output grid


  !***********************************
  ! Variables defining receptor points
  !***********************************

  ! general receptors
  real, allocatable, dimension(:) :: xreceptor,yreceptor,zreceptor
  integer, allocatable, dimension(:) :: treceptor
  real, allocatable, dimension(:) :: receptorarea
  real, allocatable, dimension(:,:) :: creceptor,crecuncert
  real, allocatable, dimension(:) :: xkreceptor,nnreceptor
  character(len=16), allocatable, dimension(:) :: receptorname
  integer :: cpointer(maxrecsample)
  integer :: numreceptor, numcurrec
  logical :: lrecregular

  ! satellite receptors
  real, allocatable, dimension(:) :: xsatellite,ysatellite
  integer, allocatable, dimension(:) :: tsatellite
  real, allocatable, dimension(:) :: satellitearea
  real, allocatable, dimension(:,:) :: zsatellite
  real, allocatable, dimension(:,:,:) :: csatellite, csatuncert
  real, allocatable, dimension(:,:)   :: xksatellite, nnsatellite
  character(len=24), allocatable, dimension(:) :: satellitename
  integer :: numsatreceptor, nlayermax, numsatellite, numcursat
  integer, allocatable, dimension(:) :: nnsatlayer
  integer :: csatpointer(maxrecsample)

  ! xreceptor,yreceptor,zreceptor     receptor position
  ! creceptor                         concentrations at receptor points
  ! receptorarea                      area of 1*1 grid cell at receptor point
  ! numreceptor                       number of receptors (non-satellite)
  ! numcurrec                         number of receptors in current time interval (updated each time interval)
  ! lrecregular                       logical if receptor output should be at regular intervals (and not according to RECEPTORS namelist)
  ! numsatreceptor                    number of satellite receptors (aka. retrievals)
  ! numcursat                         number of satellite receptors in current time interval (updated each time interval)
  ! numsatellite                      number of satellite instruments
  ! nlayermax                         max number of vertical layers in satellite retrievals
  ! nnsatlayer                        actual number of vertical layers for each satellite

  !***************************************
  ! Variables characterizing each particle
  !***************************************

  integer :: numpart=0
  integer :: numparticlecount
  integer :: maxspec ! Number of chemical species per release
  !integer :: maxndia ! Number of diameter bins (now set in par_mod.f90)
  !real, allocatable, dimension(:,:) :: xscav_frac1

  !****************************************************************
  ! Variables used for writing out interval averages for partoutput
  !****************************************************************

  integer, allocatable, dimension(:) :: npart_av
  real, allocatable, dimension(:) :: part_av_cartx,part_av_carty,part_av_cartz,part_av_z,part_av_topo
  real, allocatable, dimension(:) :: part_av_pv,part_av_qv,part_av_tt,part_av_rho,part_av_tro,part_av_hmix
  real, allocatable, dimension(:) :: part_av_uu,part_av_vv,part_av_energy

  !CGZ-lifetime
  real, allocatable, dimension(:,:) ::checklifetime, species_lifetime
  !CGZ-lifetime

  ! numpart                 actual number of particles in memory
  ! itra1 (maxpart) [s]     temporal positions of the particles
  ! npoint(maxpart)         indicates the release point of each particle
  ! nclass (maxpart)        one of nclassunc classes to which the particle is attributed
  ! itramem (maxpart) [s]   memorized release times of the particles
  ! itrasplit (maxpart) [s] next time when particle is to be split into two
  ! idt(maxpart) [s]        time step to be used for next integration
  ! numparticlecount        counts the total number of particles that have been released
  ! xtra1,ytra1,ztra1       spatial positions of the particles
  ! xmass1 [kg]             particle masses
  ! xscav_frac1             fraction of particle masse which has been scavenged at receptor
 


  !*******************************************************
  ! Info table on available chemical species/radionuclides
  !*******************************************************

  !character*10 specname(maxtable)
  !real decaytime(maxtable),wetscava(maxtable),wetscavb(maxtable)
  !real drydiff(maxtable),dryhenry(maxtable),dryactiv(maxtable)
  !real partrho(maxtable),partmean(maxtable),partsig(maxtable)
  !real dryvelo(maxtable),weightmol(maxtable),ohreact(maxtable)

  ! specname            Name of chemical species/radionuclide
  ! decaytime           Half time of radionuclides
  ! wetscava, wetscavb  Parameters for calculating scavenging coefficients
  ! drydiff             diffusivitiy of species relative to diff. of H2O
  ! dryhenry [M/atm]    Henry constant
  ! dryactiv            reactivity relative to that of O3
  ! partrho [kg/m3]     density of particles
  ! partmean [m]        mean diameter of particles
  ! partsig [m]         mean stand. deviation of particle diameter
  ! dryvelo [cm/s]      constant dry deposition velocity
  ! weightmol [g/mol]   molecular weight
  ! ohreact             OH reaction rate


  !********************
  ! Random number field
  !********************

  real :: rannumb(maxrand+2)

  ! rannumb                 field of normally distributed random numbers
  
  !********************************************************************
  ! variables to control stability of CBL scheme under variation 
  ! of statistics in time and space 
  !********************************************************************
  integer :: sum_nan_count(3600),maxtl=1200
  integer,allocatable,dimension(:) :: nan_count
  !added by mc , note that for safety sum_nan_count(N) with N>maxtl

  !********************************************************************
  ! variables to test well-mixed state of CBL scheme not to be included in final release
  !********************************************************************
  real :: well_mixed_vector(50),h_well,well_mixed_norm,avg_air_dens(50),avg_ol,avg_wst,avg_h 
  ! modified by mc to test well-mixed for cbl

  !********************
  ! Verbosity, testing flags, namelist I/O
  !********************   
  logical :: debug_mode=.false.
  integer :: verbosity=0
  integer :: info_flag=0
  integer :: count_clock, count_clock0,  count_rate, count_max
  real    :: tins
  logical, parameter :: nmlout=.true.

  !**************************************************************
  ! These variables are used to avoid having separate versions of
  ! files in cases where differences with MPI version are minor (eso)
  !*****************************************************************
  integer :: mpi_mode=0 ! .gt. 0 if running MPI version
  logical :: lroot=.true. ! true if serial version, or if MPI .and. root process
  
  logical, parameter :: interpolhmix=.false. ! true if the hmix shall be interpolated

  integer :: numthreads,numthreads_grid  ! number of available threads in parallel sections
  !integer :: nclassunc2, nrecclunc, ngriclunc

  ! Set maximum number of threads for doing grid computations in COMMAND
  ! Recommended to set this to max 16
  ! High numbers create more overhead and a larger memory footprint
  !***********************************************************************
  integer :: maxthreadgrid

  !*******************************************************************************
  ! Maximum output of each partoutput NetCDF-4 file in Mb 
  ! before a new one is created
  !*******************************************************************************

  integer :: maxfilesize

  ! This flag sets all vertical interpolation to logarithmic instead of linear
  !***************************************************************************
  integer :: logvertinterp
  logical :: log_interpol=.false.

  !*********************************************************
  !LB 04.05.2021, simple timing of IO and total running time
  !*********************************************************
  real :: s_readwind=0, s_writepartav=0, s_writepart=0, s_temp=0, s_total=0, s_firstt=0
  real, parameter :: eta_convert=1000000., zfac=100.


contains

  subroutine alloc_com
    implicit none
    integer :: stat

    allocate( icnt_belowcld(maxspec),icnt_incld(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate cnt_belowcld or icnt_incld"
    allocate( specnum(maxspec),decay(maxspec),weta_gas(maxspec), &
      wetb_gas(maxspec),crain_aero(maxspec),csnow_aero(maxspec), &
      ccn_aero(maxspec),in_aero(maxspec),ndia(maxspec), &
      reldiff(maxspec),henry(maxspec),f0(maxspec),density(maxspec), &
      dquer(maxspec),dsigma(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate particle property arrays 1"
    allocate( vsetaver(maxspec),cunningham(maxspec), &
      weightmolar(maxspec),ri(5,numclass),rac(5,numclass), &
      rcl(maxspec,5,numclass),rgs(maxspec,5,numclass), &
      rlu(maxspec,5,numclass),rm(maxspec),dryvel(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate particle property arrays 2"
    allocate( Fn(maxspec),Fs(maxspec),ks1(maxspec),ks2(maxspec), &
      kn2(maxspec),ishape(maxspec),orient(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate particle shape arrays"
    allocate( area_hour(maxspec,24),point_hour(maxspec,24), &
      area_dow(maxspec,7),point_dow(maxspec,7), &
      species(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate species arrays"
    allocate( DRYDEPSPEC(maxspec),WETDEPSPEC(maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate DRYDEPSPEC or WETDEPSPEC"
!    allocate( creceptor(numreceptor,maxspec),stat=stat)
!    if (stat.ne.0) error stop "Could not allocate creceptor"

    icnt_belowcld=0
    icnt_incld=0
  end subroutine alloc_com

  subroutine alloc_com_ndia
    implicit none
    integer :: stat

    allocate(vset(maxspec,maxndia),schmi(maxspec,maxndia),fract(maxspec,maxndia),stat=stat)
    if (stat.ne.0) error stop "Could not allocate vset,schmi or fract"
  end subroutine alloc_com_ndia

  subroutine dealloc_com
    deallocate(icnt_belowcld,icnt_incld,specnum,decay,weta_gas,wetb_gas, &
      crain_aero,csnow_aero,ccn_aero,in_aero,reldiff,henry,f0,density,dquer, &
      dsigma,ndia,vsetaver,cunningham,weightmolar,vset,schmi,fract,ri,rac,rcl, &
      rgs,rlu,rm,dryvel,Fn,Fs,ks1,ks2,kn2,ishape, &
      orient,area_hour,point_hour,area_dow,point_dow,species)
    deallocate(DRYDEPSPEC,WETDEPSPEC)
    deallocate(creceptor,xreceptor,yreceptor,receptorarea,receptorname)
    deallocate(lage)
  end subroutine dealloc_com

  subroutine mpi_alloc_part(nmpart)
  !*******************************************************************************    
  ! Dynamic allocation of arrays
  !
  ! For FLEXPART version 9.2 and earlier these arrays were statically declared
  ! with size maxpart. This function is introduced so that the MPI version
  ! can declare these arrays with smaller size ("maxpart per process"), while
  ! the serial version allocate at run-time with size maxpart 
  !
  !*******************************************************************************
    implicit none 

    integer, intent(in) :: nmpart ! maximum number of particles (per process)

    if (ipout.eq.3) then
      allocate(npart_av(nmpart),part_av_cartx(nmpart),part_av_carty(nmpart),&
           & part_av_cartz(nmpart),part_av_z(nmpart),part_av_topo(nmpart))
      allocate(part_av_pv(nmpart),part_av_qv(nmpart),part_av_tt(nmpart),&
           & part_av_rho(nmpart),part_av_tro(nmpart),part_av_hmix(nmpart))
      allocate(part_av_uu(nmpart),part_av_vv(nmpart),part_av_energy(nmpart))
    end if

  end subroutine mpi_alloc_part

  subroutine update_gitversion(gitversion_tmp)
    character(len=256) :: gitversion_tmp
    gitversion=gitversion_tmp
  end subroutine
end module com_mod
