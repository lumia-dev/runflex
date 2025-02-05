! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !  This module handles all gridded netcdf output for concentration or        *
  !  residence time and wet and dry deposition output.                         *
  !                                                                            *
  !  - writeheader_netcdf generates file including all information previously  *
  !    stored in separate header files                                         *
  !  - concoutput_netcdf write concentration output and wet and dry deposition *
  !                                                                            *
  !     Author: D. Brunner                                                     *
  !                                                                            *
  !     12 April 2013                                                          *
  !                                                                            *
  ! HSO: 21 Oct 2014                                                           *
  !  - added option to not writeout releases information by changing           *
  !    switch write_releases                                                   *
  !  - additional updates for FLEXPART 9.x                                     *
  !                                                                            *
  ! ESO 2016                                                                   *
  !  - Deposition fields can be calculated in double precision, see variable   *
  !    'dep_prec' in par_mod                                                   *
  !  - Hardcoded options 'write_vol' and 'write_area' for grid cell            *
  !    volume and area                                                         *
  !                                                                            *
  ! LB: 2021                                                                   *
  !  - Particle dump and initial particle positions in NetCDF                  *
  !  - Receptor files in NetCDF format                                         *
  !                                                                            *
  ! RLT: 2024                                                                  *
  !  - Moved receptor output new module                                        *
  !*****************************************************************************


module netcdf_output_mod

  use netcdf

  use point_mod, only: ireleasestart,ireleaseend,kindz,dx,xlon0,dy,ylat0,&
                       xpoint1,ypoint1,xpoint2,ypoint2,zpoint1,zpoint2,npart,xmass
  use outgrid_mod,  only: outheight,oroout,densityoutgrid,factor3d,volume,&
                       wetgrid,wetgridsigma,drygrid,drygridsigma,grid,gridsigma,&
                       area,arean,volumen,orooutn
  use par_mod,   only: dep_prec, sp, dp, nclassunc,&
                       unitoutrecept,unitoutreceptppt,unittmp,lpartoutputperfield
  use com_mod

  use windfields_mod, only: oro,rho,nxmax,height,nxmin1,nymin1,nz,nx,ny,hmix, &
                       ! for concoutput_netcdf and concoutput_nest_netcdf 
                       tropopause,oron,rhon,xresoln,yresoln,xrn,xln,yrn,yln,nxn,nyn
  use mean_mod

  implicit none

  !  include 'netcdf.inc'

  ! parameter for data compression (1-9, 9 = most aggressive)
  integer, parameter :: deflate_level = 5
  logical, parameter :: min_size = .false.   ! if set true, redundant fields (topography) are not written to minimize file size

  integer            :: tpointer=0,tpointer_part=0,ppointer_part=0,partinitpointer=0,partinitpointer1=0
  character(len=255) :: ncfname, ncfnamen, ncfname_part, ncfname_partinit, ncfname_part_end

  ! netcdf dimension and variable IDs for main and nested output grid
  integer,allocatable,dimension(:) :: specID,specIDppt,wdspecID,ddspecID
  integer,allocatable,dimension(:) :: specIDn,specIDnppt,wdspecIDn,ddspecIDn, &
    recconcID,recpptvID
  integer                     :: timeID, timeIDn, timeIDpart
  integer, dimension(6)       :: dimids, dimidsn
  integer, dimension(5)       :: depdimids, depdimidsn

  ! For initial particle outputs
  integer  :: partIDi,tIDi,lonIDi,latIDi,levIDi,relIDi,pvIDi,prIDi,qvIDi, &
    rhoIDi,ttIDi,uIDi,vIDi,wIDi,topoIDi,trIDi,hmixIDi
  integer,allocatable,dimension(:) :: massIDi

  real :: eps
  !  private:: writemetadata, output_units, nf90_err

  ! switch output of release point information on/off
  logical, parameter :: write_releases = .true.

  ! switch output of grid cell volume and area on/off
  logical, parameter :: write_vol = .false.
  logical, parameter :: write_area = .false.

  ! switch for first time topo output in case of domainfill
  logical :: topo_written=.false.
  logical :: mass_written=.false.
  logical :: massav_written=.false.

  ! coordinate transformation from internal to world coord
  real :: xp1,yp1,xp2,yp2


  private

  public :: writeheader_netcdf,concoutput_netcdf,&
       concoutput_nest_netcdf,writeheader_partoutput,partoutput_netcdf,&
       open_partoutput_file,close_partoutput_file,create_particles_initialoutput,&
       topo_written,mass_written,wrt_part_initialpos,partinit_netcdf,open_partinit_file, &
       readpartpositions_netcdf,readinitconditions_netcdf,partinitpointer1,tpointer, &
       alloc_netcdf,dealloc_netcdf,nf90_err,update_partoutput_pointers,ppointer_part

contains

subroutine alloc_netcdf
  implicit none
  integer :: stat

  allocate( specID(maxspec),specIDppt(maxspec), wdspecID(maxspec),ddspecID(maxspec), &
    specIDn(maxspec),specIDnppt(maxspec), wdspecIDn(maxspec),ddspecIDn(maxspec), &
    recconcID(maxspec),recpptvID(maxspec), stat=stat)
  if (stat.ne.0) error stop "Could not allocate netcdf fields"

  allocate( massIDi(maxspec), stat=stat)
  if (stat.ne.0) error stop "Could not allocate netcdf fields 2"
end subroutine alloc_netcdf

subroutine dealloc_netcdf
  deallocate( specID,specIDppt,wdspecID,ddspecID,specIDn,specIDnppt,wdspecIDn,ddspecIDn, &
    recconcID,recpptvID )
  deallocate( massIDi )
end subroutine dealloc_netcdf

!****************************************************************
! determine output units (see table 1 in Stohl et al., ACP 2005
!****************************************************************
subroutine output_units(units)
  implicit none
  character(len=15), intent(out) :: units
  if (ldirect.eq.1) then
     ! forward simulation
     if (ind_source.eq.1) then
        if (ind_receptor.eq.1) then
           units = 'ng m-3'   ! hes the kg in Tab1 is only indicating the units of the relase not the output
        else
           units = 'ng kg-1'
        endif
     else
        if (ind_receptor.eq.1) then
           units = 'ng m-3'
        else
           units = 'ng kg-1'
        endif
     endif
  else
     ! backward simulation
     if (ind_source.eq.1) then
        if (ind_receptor.eq.1) then
           units = 's'
        else 
           units = 's m3 kg-1'
        endif
     else
        if (ind_receptor.eq.1) then
           units = 's kg m-3'
        else
           units = 's'
        endif
     endif
  endif
end subroutine output_units


!****************************************************************
! write metadata to netCDF file 
!****************************************************************
subroutine writemetadata(ncid,lnest)
  
  implicit none 

  integer, intent(in) :: ncid
  logical, intent(in) :: lnest
  character           :: time*10,date*8,adate*8,atime*6
  character(5)        :: zone
  character(255)      :: login_name, host_name

  ! gather system information 
  call date_and_time(date,time,zone)
  call getlog(login_name)
  call hostnm(host_name)
  
  ! hes CF convention requires these attributes
  call nf90_err(nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.6'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'title', 'FLEXPART model output'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'git', trim(gitversion)))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'source', trim(flexversion)//' model output'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'history', date(1:4)//'-'//date(5:6)// &
       '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//' '//zone//'  created by '//  &
       trim(login_name)//' on '//trim(host_name)))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'references', &
       'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200'))

  ! attributes describing model run
  !************************************************************************************

  if (lnest) then
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlon0', outlon0n))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlat0', outlat0n))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dxout', dxoutn))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dyout', dyoutn))
  else
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlon0', outlon0))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlat0', outlat0))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dxout', dxout))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dyout', dyout))
  endif
  !   vertical levels stored in grid structure

  ! COMMAND file settings
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ldirect', ldirect))
  write(adate,'(i8.8)') ibdate
  write(atime,'(i6.6)') ibtime
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ibdate', adate))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ibtime', atime))
  write(adate,'(i8.8)') iedate
  write(atime,'(i6.6)') ietime
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iedate', adate))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ietime', atime))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutstep', loutstep))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutaver', loutaver))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutsample', loutsample))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutrestart', loutrestart))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lsynctime', lsynctime))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ctl', ctl))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ifine', ifine))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iout', iout))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ipout', ipout))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lsubgrid', lsubgrid))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lconvection', lconvection))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lagespectra', lagespectra))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ipin', ipin))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ioutputforeachrelease', ioutputforeachrelease))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iflux', iflux))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'mdomainfill', mdomainfill))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ind_source', ind_source))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ind_receptor', ind_receptor))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'mquasilag', mquasilag))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'nested_output', nested_output))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'sfc_only', sfc_only))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'linit_cond', linit_cond))
end subroutine writemetadata


subroutine nf90_err(status)
  !****************************************************************
  ! netcdf error message handling
  !****************************************************************
  implicit none
  integer, intent (in) :: status
   if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'Stopped'
    end if
end subroutine nf90_err

subroutine writeheader_netcdf(lnest)

  !****************************************************************
  ! Create netcdf file and write header/metadata information
  ! lnest = .false. : Create main output file
  ! lnest = .true.  : Create nested output file
  !****************************************************************
  implicit none

  logical, intent(in) :: lnest

  integer :: ncid, sID, wdsID, ddsID
  integer :: timeDimID, latDimID, lonDimID, levDimID
  integer :: nspecDimID, npointDimID, nageclassDimID, ncharDimID, pointspecDimID
  integer :: tID, lonID, latID, levID, lageID, oroID, ncharrecDimID
  integer :: volID, areaID
  integer :: rellng1ID, rellng2ID, rellat1ID, rellat2ID, relzz1ID, relzz2ID
  integer :: relcomID, relkindzID, relstartID, relendID, relpartID, relxmassID
  integer :: nnx, nny 
  integer, dimension(6)       :: dIDs
  integer, dimension(5)       :: depdIDs
  character(len=255)          :: fname
  character(len=15)           :: units
  character(len=20)           :: fprefix
  character(len=3)            :: anspec
  CHARACTER                   :: adate*8,atime*6,timeunit*32
  !REAL, DIMENSION(1000)       :: coord
  real, allocatable, dimension(:) :: coord

  integer                     :: cache_size
  integer, dimension(6)       :: chunksizes
  integer, dimension(5)       :: dep_chunksizes

  integer                     :: i
  integer                     :: numzwrite


  ! Check if output directory exists (the netcdf library will
  ! otherwise give an error which can look confusing). 
  ! *********************************************************************
  open(unit=unittmp,file=trim(path(2)(1:length(2)))//'test_dir.txt',status='replace',&
       &err=100)
  close (unittmp, status='delete')
  goto 101
100 write(*,FMT='(80("#"))') 
  write(*,*) 'ERROR: output directory ', trim(path(2)(1:length(2))), ' does not exist&
       & (or failed to write there).' 
  write(*,*) 'EXITING' 
  write(*,FMT='(80("#"))')
  error stop
101 continue

  !************************
  ! Create netcdf file
  !************************

  numzwrite=numzgrid
  if (sfc_only.eq.1) numzwrite=1

  if (ldirect.eq.1) then
     write(adate,'(i8.8)') ibdate
     write(atime,'(i6.6)') ibtime
     fprefix = 'grid_conc_'
  else
     write(adate,'(i8.8)') iedate
     write(atime,'(i6.6)') ietime
     fprefix = 'grid_time_'
  endif
  if (DRYBKDEP) fprefix='grid_drydep_'
  if (WETBKDEP) fprefix='grid_wetdep_'

  if (lnest) then
     fname = path(2)(1:length(2))//trim(fprefix)//adate//atime//'_nest.nc'
     ncfnamen = fname
     nnx = numxgridn
     nny = numygridn
  else
     fname = path(2)(1:length(2))//trim(fprefix)//adate//atime//'.nc'
     ncfname = fname
     nnx = numxgrid
     nny = numygrid
  endif

  cache_size = 16 * nnx * nny * numzgrid

  ! If starting from a restart file, new data will be added to the existing grid file
  if ((ipin.eq.1).or.(ipin.eq.4)) then
    call read_grid_id(lnest)
    return
  endif

  ! setting cache size in bytes. It is set to 4 times the largest data block that is written
  !   size_type x nx x ny x nz
  ! create file
  call nf90_err(nf90_create(trim(fname), cmode = nf90_hdf5, ncid = ncid, &
    cache_size = cache_size))  

  ! create dimensions:
  !*************************
  ! time
  call nf90_err(nf90_def_dim(ncid, 'time', nf90_unlimited, timeDimID))
  timeunit = 'seconds since '//adate(1:4)//'-'//adate(5:6)// &
     '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)

  ! lon
  call nf90_err(nf90_def_dim(ncid, 'longitude', nnx, lonDimID))
  ! lat
  call nf90_err(nf90_def_dim(ncid, 'latitude', nny, latDimID))
  ! level
!  call nf90_err(nf90_def_dim(ncid, 'height', numzgrid, levDimID))
  call nf90_err(nf90_def_dim(ncid, 'height', numzwrite, levDimID))
  ! number of species
  call nf90_err(nf90_def_dim(ncid, 'numspec', nspec, nspecDimID))
  ! number of release points
  call nf90_err(nf90_def_dim(ncid, 'pointspec', maxpointspec_act, pointspecDimID))
  ! number of age classes
  call nf90_err(nf90_def_dim(ncid, 'nageclass', nageclass, nageclassDimID))
  ! dimension for release point characters
  call nf90_err(nf90_def_dim(ncid, 'nchar', 45, ncharDimID))
  ! number of actual release points
  call nf90_err(nf90_def_dim(ncid, 'numpoint', numpoint, npointDimID))


  ! create variables
  !*************************

  ! time
  call nf90_err(nf90_def_var(ncid, 'time', nf90_int, (/ timeDimID /), tID))
  call nf90_err(nf90_put_att(ncid, tID, 'units', timeunit))
  call nf90_err(nf90_put_att(ncid, tID, 'calendar', 'proleptic_gregorian'))
  if (lnest) then
     timeIDn = tID
  else
     timeID = tID
  endif

  ! lon
  call nf90_err(nf90_def_var(ncid, 'longitude', nf90_float, (/ lonDimID /), lonID))
  call nf90_err(nf90_put_att(ncid, lonID, 'long_name', 'longitude in degree east'))
  call nf90_err(nf90_put_att(ncid, lonID, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonID, 'units', 'degrees_east'))
  call nf90_err(nf90_put_att(ncid, lonID, 'standard_name', 'grid_longitude'))
  call nf90_err(nf90_put_att(ncid, lonID, 'description', 'grid cell centers'))

  ! lat
  call nf90_err(nf90_def_var(ncid, 'latitude', nf90_float, (/ latDimID /), latID))
  call nf90_err(nf90_put_att(ncid, latID, 'long_name', 'latitude in degree north'))
  call nf90_err(nf90_put_att(ncid, latID, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latID, 'units', 'degrees_north'))
  call nf90_err(nf90_put_att(ncid, latID, 'standard_name', 'grid_latitude'))
  call nf90_err(nf90_put_att(ncid, latID, 'description', 'grid cell centers'))


  ! height
  call nf90_err(nf90_def_var(ncid, 'height', nf90_float, (/ levDimID /), levID))
  ! call nf90_err(nf90_put_att(ncid, levID, 'axis', 'Z'))
  call nf90_err(nf90_put_att(ncid, levID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, levID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, levID, 'standard_name', 'height'))
  call nf90_err(nf90_put_att(ncid, levID, 'long_name', 'height above ground'))

  ! volume
  if (write_vol) call nf90_err(nf90_def_var(ncid, 'volume', nf90_float, &
       &(/ lonDimID, latDimID, levDimID /), volID))
  ! area 
  if (write_area) call nf90_err(nf90_def_var(ncid, 'area', nf90_float, &
       &(/ lonDimID, latDimID /), areaID))


  if (write_releases.eqv..true.) then
    ! release comment
    call nf90_err(nf90_def_var(ncid, 'RELCOM', nf90_char, (/ ncharDimID,npointDimID /), &
         relcomID))
    call nf90_err(nf90_put_att(ncid, relcomID, 'long_name', 'release point name'))
    ! release longitude 1
    call nf90_err(nf90_def_var(ncid, 'RELLNG1', nf90_float, (/ npointDimID /), rellng1ID))
    call nf90_err(nf90_put_att(ncid, rellng1ID, 'units', 'degrees_east'))
    call nf90_err(nf90_put_att(ncid, rellng1ID, 'long_name', &
         'release longitude lower left corner'))
    ! release longitude 2
    call nf90_err(nf90_def_var(ncid, 'RELLNG2', nf90_float, (/ npointDimID /), rellng2ID))
    call nf90_err(nf90_put_att(ncid, rellng2ID, 'units', 'degrees_east'))
    call nf90_err(nf90_put_att(ncid, rellng2ID, 'long_name', &
         'release longitude upper right corner'))
    ! release latitude 1
    call nf90_err(nf90_def_var(ncid, 'RELLAT1', nf90_float, (/ npointDimID /), rellat1ID))
    call nf90_err(nf90_put_att(ncid, rellat1ID, 'units', 'degrees_north'))
    call nf90_err(nf90_put_att(ncid, rellat1ID, 'long_name', &
         'release latitude lower left corner'))
    ! release latitude 2
    call nf90_err(nf90_def_var(ncid, 'RELLAT2', nf90_float, (/ npointDimID /), rellat2ID))
    call nf90_err(nf90_put_att(ncid, rellat2ID, 'units', 'degrees_north'))
    call nf90_err(nf90_put_att(ncid, rellat2ID, 'long_name', &
         'release latitude upper right corner'))

    ! hes: if rotated_ll it would be convenient also to write the the release points in rotated_coordinates

    ! release height bottom
    call nf90_err(nf90_def_var(ncid, 'RELZZ1', nf90_float, (/ npointDimID /), relzz1ID))
    call nf90_err(nf90_put_att(ncid, relzz1ID, 'units', 'meters'))
    call nf90_err(nf90_put_att(ncid, relzz1ID, 'long_name', 'release height bottom'))
    ! release height top
    call nf90_err(nf90_def_var(ncid, 'RELZZ2', nf90_float, (/ npointDimID /), relzz2ID))
    call nf90_err(nf90_put_att(ncid, relzz2ID, 'units', 'meters'))
    call nf90_err(nf90_put_att(ncid, relzz2ID, 'long_name', 'release height top'))
    ! release kind
    call nf90_err(nf90_def_var(ncid, 'RELKINDZ', nf90_int, (/ npointDimID /), relkindzID))
    call nf90_err(nf90_put_att(ncid, relkindzID, 'long_name', 'release kind'))
    ! release start
    call nf90_err(nf90_def_var(ncid, 'RELSTART', nf90_int, (/ npointDimID /), relstartID))
    call nf90_err(nf90_put_att(ncid, relstartID, 'units', 'seconds'))
    call nf90_err(nf90_put_att(ncid, relstartID, 'long_name', &
         'release start relative to simulation start'))
    ! release end
    call nf90_err(nf90_def_var(ncid, 'RELEND', nf90_int, (/ npointDimID /), relendID))
    call nf90_err(nf90_put_att(ncid, relendID, 'units', 'seconds'))
    call nf90_err(nf90_put_att(ncid, relendID, 'long_name', &
         'release end relative to simulation start'))
    ! release particles
    call nf90_err(nf90_def_var(ncid, 'RELPART', nf90_int, (/ npointDimID /), relpartID))
    call nf90_err(nf90_put_att(ncid, relpartID, 'long_name', 'number of release particles'))
    ! release particle masses
    call nf90_err(nf90_def_var(ncid, 'RELXMASS', nf90_float, (/ npointDimID, nspecDimID /), &
         relxmassID))
    call nf90_err(nf90_put_att(ncid, relxmassID, 'long_name', 'total release particle mass'))
  end if
 
  ! age classes
  call nf90_err(nf90_def_var(ncid, 'LAGE', nf90_int, (/ nageclassDimID /), lageID))
  call nf90_err(nf90_put_att(ncid, lageID, 'units', 'seconds'))
  call nf90_err(nf90_put_att(ncid, lageID, 'long_name', 'age class'))

  ! output orography
  if (.not. min_size) then
    call nf90_err(nf90_def_var(ncid, 'ORO', nf90_int, (/ lonDimID, latDimID /), oroID,  &
      deflate_level=deflate_level, chunksizes= (/ nnx, nny /)))
    call nf90_err(nf90_put_att(ncid, oroID, 'standard_name', 'surface altitude'))
    call nf90_err(nf90_put_att(ncid, oroID, 'long_name', 'outgrid surface altitude'))
    call nf90_err(nf90_put_att(ncid, oroID, 'units', 'm'))
  end if

  ! concentration output, wet and dry deposition variables (one per species)
  call output_units(units)

  dIDs = (/ londimid, latdimid, levdimid, timedimid, pointspecdimid, nageclassdimid /)
  depdIDs = (/ londimid, latdimid, timedimid, pointspecdimid, nageclassdimid /)
  if (lnest) then
     dimidsn    = dIDs
     depdimidsn = depdIDs
  else
     dimids    = dIDs
     depdimids = depdIDs
  endif

  ! set chunksizes according to largest written portion of data in an individual call to 
  ! nf90_put_var
  if (int(nnx,kind=8)*int(nny,kind=8)*int(numzgrid,kind=8).gt.2147483647) then ! Larger than an 
    chunksizes = (/ nnx, nny, 1, 1, 1, 1 /)
  else
!    chunksizes = (/ nnx, nny, numzgrid, 1, 1, 1 /)
    chunksizes = (/ nnx, nny, numzwrite, 1, 1, 1 /)
  endif
  dep_chunksizes = (/ nnx, nny, 1, 1, 1 /)

  do i = 1,nspec
    write(anspec,'(i3.3)') i

    ! concentration output
    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      call nf90_err(nf90_def_var(ncid,'spec'//anspec//'_mr', nf90_float, dIDs, sID , &
           deflate_level = deflate_level,  &
           chunksizes = chunksizes ))
      call nf90_err(nf90_put_att(ncid, sID, 'units', units))
      call nf90_err(nf90_put_att(ncid, sID, 'long_name', species(i)))
      call nf90_err(nf90_put_att(ncid, sID, 'decay', decay(i)))
      call nf90_err(nf90_put_att(ncid, sID, 'weightmolar', weightmolar(i)))
    !        call nf90_err(nf90_put_att(ncid, sID, 'ohreact', ohreact(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'ohcconst', ohcconst(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'ohdconst', ohdconst(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'vsetaver', vsetaver(i)))

      if (lnest) then
         specIDn(i) = sID
      else
         specID(i) = sID
      endif
    endif

    ! mixing ratio output
    if ((iout.eq.2).or.(iout.eq.3)) then
      call nf90_err(nf90_def_var(ncid,'spec'//anspec//'_pptv', nf90_float, dIDs, sID , &
           deflate_level = deflate_level,  &
           chunksizes = chunksizes ))
      call nf90_err(nf90_put_att(ncid, sID, 'units', 'pptv'))
      call nf90_err(nf90_put_att(ncid, sID, 'long_name', species(i)))
      call nf90_err(nf90_put_att(ncid, sID, 'decay', decay(i)))
      call nf90_err(nf90_put_att(ncid, sID, 'weightmolar', weightmolar(i)))
    !        call nf90_err(nf90_put_att(ncid, sID, 'ohreact', ohreact(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'ohcconst', ohcconst(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'ohdconst', ohdconst(i)))
!      call nf90_err(nf90_put_att(ncid, sID, 'vsetaver', vsetaver(i)))

      if (lnest) then
         specIDnppt(i) = sID
      else
         specIDppt(i) = sID
      endif
    endif

    ! wet and dry deposition fields for forward runs
    if ((ldirect.eq.1).and.(wetdep)) then
      call nf90_err(nf90_def_var(ncid,'WD_spec'//anspec, nf90_float, depdIDs, &
           wdsID, deflate_level = deflate_level, &
           chunksizes = dep_chunksizes))
      call nf90_err(nf90_put_att(ncid, wdsID, 'units', '1e-12 kg m-2'))
      call nf90_err(nf90_put_att(ncid, wdsID, 'weta_gas', weta_gas(i)))
      call nf90_err(nf90_put_att(ncid, wdsID, 'wetb_gas', wetb_gas(i)))
      call nf90_err(nf90_put_att(ncid, wdsID, 'ccn_aero', ccn_aero(i)))
      call nf90_err(nf90_put_att(ncid, wdsID, 'in_aero', in_aero(i)))
      ! call nf90_err(nf90_put_att(ncid, wdsID, 'wetc_in', wetc_in(i)))
      ! call nf90_err(nf90_put_att(ncid, wdsID, 'wetd_in', wetd_in(i)))
      call nf90_err(nf90_put_att(ncid, wdsID, 'dquer', dquer(i)))
      call nf90_err(nf90_put_att(ncid, wdsID, 'henry', henry(i)))
      if (lnest) then
         wdspecIDn(i) = wdsID
      else
         wdspecID(i) = wdsID
      endif
    endif
    if ((ldirect.eq.1).and.(drydep)) then
      call nf90_err(nf90_def_var(ncid,'DD_spec'//anspec, nf90_float, depdIDs, &
           ddsID, deflate_level = deflate_level, &
           chunksizes = dep_chunksizes))
      call nf90_err(nf90_put_att(ncid, ddsID, 'units', '1e-12 kg m-2'))
      call nf90_err(nf90_put_att(ncid, ddsID, 'dryvel', dryvel(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'reldiff', reldiff(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'henry', henry(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'f0', f0(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'dquer', dquer(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'density', density(i)))
      call nf90_err(nf90_put_att(ncid, ddsID, 'dsigma', dsigma(i)))
      if (lnest) then
         ddspecIDn(i) = ddsID
      else
         ddspecID(i) = ddsID
      endif
    endif
  end do

  ! global (metadata) attributes
  !*******************************
  call writemetadata(ncid,lnest)


  ! moves the file from define to data mode
  call nf90_err(nf90_enddef(ncid))

  ! fill with data
  !******************************
  ! longitudes (grid cell centers)
  if (lnest) then
    if (.not.allocated(coord)) allocate(coord(numxgridn))
     do i = 1,numxgridn
        coord(i) = outlon0n + (i-0.5)*dxoutn
     enddo
     call nf90_err(nf90_put_var(ncid, lonID, coord(1:numxgridn)))
     deallocate(coord)
  else
    if (.not.allocated(coord)) allocate(coord(numxgrid))
     do i = 1,numxgrid
        coord(i) = outlon0 + (i-0.5)*dxout
     enddo
     call nf90_err(nf90_put_var(ncid, lonID, coord(1:numxgrid)))
     deallocate(coord)
  endif
  ! latitudes (grid cell centers)
  if (lnest) then
    if (.not.allocated(coord)) allocate(coord(numygridn))
     do i = 1,numygridn
        coord(i) = outlat0n + (i-0.5)*dyoutn
     enddo
     call nf90_err(nf90_put_var(ncid, latID, coord(1:numygridn)))
     deallocate(coord)
  else
    if (.not.allocated(coord)) allocate(coord(numygrid))
     do i = 1,numygrid
        coord(i) = outlat0 + (i-0.5)*dyout
     enddo
     call nf90_err(nf90_put_var(ncid, latID, coord(1:numygrid)))
     deallocate(coord)
  endif
  ! levels
!  call nf90_err(nf90_put_var(ncid, levID, outheight(1:numzgrid)))
  call nf90_err(nf90_put_var(ncid, levID, outheight(1:numzwrite)))

  ! volume
  if (write_vol) then
    if (lnest) then
!      call nf90_err(nf90_put_var(ncid, volID, volumen(:,:,:)))
      call nf90_err(nf90_put_var(ncid, volID, volumen(:,:,1:numzwrite)))
    else
!      call nf90_err(nf90_put_var(ncid, volID, volume(:,:,:)))
      call nf90_err(nf90_put_var(ncid, volID, volume(:,:,1:numzwrite)))
    end if
  end if

  ! area
  if (write_area) then
    if (lnest) then
      call nf90_err(nf90_put_var(ncid, areaID, arean(:,:)))
    else
      call nf90_err(nf90_put_var(ncid, areaID, area(:,:)))
    end if
  end if

  if ((write_releases.eqv..true.).and.(ipin.ne.3).and.(ipin.ne.4)) then
    ! release point information
    do i = 1,numpoint
       call nf90_err(nf90_put_var(ncid, relstartID, ireleasestart(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relendID, ireleaseend(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relkindzID, kindz(i), (/i/)))
       xp1=xpoint1(i)*dx+xlon0
       yp1=ypoint1(i)*dy+ylat0
       xp2=xpoint2(i)*dx+xlon0
       yp2=ypoint2(i)*dy+ylat0
       call nf90_err(nf90_put_var(ncid, rellng1ID, xp1, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellng2ID, xp2, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellat1ID, yp1, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellat2ID, yp2, (/i/)))
       call nf90_err(nf90_put_var(ncid, relzz1ID, zpoint1(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relzz2ID, zpoint2(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relpartID, npart(i), (/i/)))
       if ((i .le. 1000).and.(ipin.ne.3).and.(ipin.ne.4)) then
         call nf90_err(nf90_put_var(ncid, relcomID, compoint(i), (/1,i/), (/45,1/)))
       else
         call nf90_err(nf90_put_var(ncid, relcomID, 'NA', (/1,i/), (/45,1/)))
       endif 
       call nf90_err(nf90_put_var(ncid, relxmassID, xmass(i,1:nspec), (/i,1/), (/1,nspec/)))
    end do
  end if

  ! age classes
  call nf90_err(nf90_put_var(ncid, lageID, lage(1:nageclass)))

  ! orography 
  if (.not. min_size) then
    if (lnest) then
      call nf90_err(nf90_put_var(ncid, oroID, orooutn(0:(nnx-1), 0:(nny-1))))
    else
      call nf90_err(nf90_put_var(ncid, oroID, oroout(0:(nnx-1), 0:(nny-1))))
    endif
  end if

  call nf90_err(nf90_close(ncid))

  return
end subroutine writeheader_netcdf

subroutine read_grid_id(lnest)
  
  implicit none
  logical, intent(in) :: lnest

  integer :: ncid,i
  character(len=3)            :: anspec

  if (.not. lnest) then
    ! open output file
    call nf90_err(nf90_open(trim(ncfname), nf90_write, ncid))

    call nf90_err(nf90_inq_varid(ncid=ncid,name='time',varid=timeID))

    do i = 1,nspec
      write(anspec,'(i3.3)') i

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='spec'//anspec//'_mr',varid=specID(i)))
      endif
      if ((iout.eq.2).or.(iout.eq.3)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='spec'//anspec//'_pptv',varid=specIDppt(i)))
      endif
      if ((ldirect.eq.1).and.(wetdep)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='WD_spec'//anspec,varid=wdspecID(i)))
      endif
      if ((ldirect.eq.1).and.(drydep)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='DD_spec'//anspec,varid=ddspecID(i)))
      endif
    end do

  else

    ! open output file
    call nf90_err(nf90_open(trim(ncfnamen), nf90_write, ncid))

    call nf90_err(nf90_inq_varid(ncid=ncid,name='time',varid=timeIDn))

    do i = 1,nspec
      write(anspec,'(i3.3)') i

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='spec'//anspec//'_mr',varid=specIDn(i)))
      endif
      if ((iout.eq.2).or.(iout.eq.3)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='spec'//anspec//'_pptv',varid=specIDnppt(i)))
      endif
      if ((ldirect.eq.1).and.(wetdep)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='WD_spec'//anspec,varid=wdspecIDn(i)))
      endif
      if ((ldirect.eq.1).and.(drydep)) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name='DD_spec'//anspec,varid=ddspecIDn(i)))
      endif
    end do 
  endif

  call nf90_err(nf90_close(ncid))

end subroutine read_grid_id

subroutine concoutput_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
     
  !                          i     i          o             o
  !       o
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the concentrations.               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output in *
  !                    sparse matrix format; additional output of uncertainty  *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in module par_mod                  *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !     February 2010, Dominik Brunner, Empa                                   *
  !                    Adapted for COSMO                                       *
  !                    Remark: calculation of density could be improved.       *
  !                    Currently, it is calculated for the lower left corner   *
  !                    of each output grid cell rather than for its center.    *
  !                    Furthermore, the average density could be calculated    *
  !                    from the difference in pressure at the top and bottom   *
  !                    of each cell rather than by interpolation.              *
  !                                                                            *
  !     April 2013, Dominik Brunner, Empa                                      *
  !                    Adapted for netcdf output                               *
  !                                                                            *
  !     2022, Lucie Bakels:                                                    *
  !           - OpenMP parallelisation                                         *
  !           - Receptor output to NetCDF instead of binary format             *
  !                                                                            *
  !     January, 2024, Rona Thompson                                           *
  !           - removed output of receptors (new module)                       *
  !           - introduced option for LCM output                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum
  real(dep_prec),intent(out):: wetgridtotalunc,drygridtotalunc
  real, intent(out)   :: gridtotalunc
  integer             :: ncid,kp,ks,kz,ix,jy,iix,jjy,kzz,ngrid
  integer             :: ks_start
  integer             :: nage,i,l,jj
  real                :: tot_mu(maxspec,maxpointspec_act)
  real                :: halfheight,dz,dz1,dz2
  real                :: xl,yl
  real(dep_prec)      :: auxgrid(nclassunc)
  real(dep_prec)      :: gridtotal,gridsigmatotal
  real(dep_prec)      :: wetgridtotal,wetgridsigmatotal
  real(dep_prec)      :: drygridtotal,drygridsigmatotal
  ! real(sp)            :: gridtotal,gridsigmatotal
  ! real(sp)            :: wetgridtotal,wetgridsigmatotal
  ! real(sp)            :: drygridtotal,drygridsigmatotal
  integer             :: numzwrite

  real, parameter     :: weightair=28.97

  eps=nxmax/3.e5

  numzwrite=numzgrid
  if (sfc_only.eq.1 ) numzwrite=1

  ! open output file
  call nf90_err(nf90_open(trim(ncfname), nf90_write, ncid))

  ! write time
  tpointer = tpointer + 1
  call nf90_err(nf90_put_var( ncid, timeID, itime, (/ tpointer /)))
  
  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************

  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1.0
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0._dep_prec
  wetgridsigmatotal=0._dep_prec
  wetgridtotalunc=0._dep_prec
  drygridtotal=0._dep_prec
  drygridsigmatotal=0._dep_prec
  drygridtotalunc=0._dep_prec

  !*******************************************************************
  ! Compute air density:
  ! brd134: we now take into account whether we are in the mother or in
  !    a nested domain (before only from mother domain)
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !
  ! Note:
  !  llcmoutput = true: grid is mass_spec/mass_air  
  !                     for iout 1,3, or 5 multiply by rho
  !                     for iout 2 multiply by 1
  !  llcmoutput = false: grid is mass_spec/V
  !                     for iout 1,3, or 5 multiply by 1
  !                     for iout 2 multiply by 1/rho
  !*******************************************************************

!$OMP PARALLEL PRIVATE(halfheight,kzz,dz1,dz2,dz,xl,yl,ngrid,iix,jjy, &
!$OMP kz,ix,jy,l,ks,kp,nage,auxgrid) REDUCTION(+:wetgridtotal,wetgridsigmatotal, &
!$OMP drygridtotal,drygridsigmatotal,gridtotal,gridsigmatotal)

  if (((.not.llcmoutput).and.(iout.eq.2)).or.&
      (llcmoutput.and.((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)))) then
    ! compute density
!$OMP DO
    do kz=1,numzgrid
      if (kz.eq.1) then
        halfheight=outheight(1)*0.5
      else
        halfheight=(outheight(kz)+outheight(kz-1))*0.5
      endif
      do kzz=2,nz
        if ((height(kzz-1).lt.halfheight).and. &
             (height(kzz).gt.halfheight)) exit
      end do
      kzz=max(min(kzz,nz),2)
      dz1=halfheight-height(kzz-1)
      dz2=height(kzz)-halfheight
      dz=dz1+dz2

      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          xl=outlon0+real(ix)*dxout
          yl=outlat0+real(jy)*dyout
          ! grid index in mother domain
          xl=(xl-xlon0)/dx
          yl=(yl-ylat0)/dx

          ngrid=0
          do jj=numbnests,1,-1
            if ( xl.gt.xln(jj)+eps .and. xl.lt.xrn(jj)-eps .and. &
                   yl.gt.yln(jj)+eps .and. yl.lt.yrn(jj)-eps ) then
              ngrid=jj
              exit 
            end if
          end do

          if (ngrid.eq.0) then
            iix=max(min(nint(xl),nxmin1),0) ! if output grid cell is outside mother domain
            jjy=max(min(nint(yl),nymin1),0)

            densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,memind(2))*dz1+ &
              rho(iix,jjy,kzz-1,memind(2))*dz2)/dz
          else
            xl=(xl-xln(ngrid))*xresoln(ngrid)
            yl=(yl-yln(ngrid))*yresoln(ngrid)
            iix=max(min(nint(xl),nxn(ngrid)-1),0)
            jjy=max(min(nint(yl),nyn(ngrid)-1),0)

            densityoutgrid(ix,jy,kz)=(rhon(iix,jjy,kzz,memind(2), ngrid)*dz1+ &
              rhon(iix,jjy,kzz-1,memind(2), ngrid)*dz2)/dz
          endif
        end do
      end do
    end do
!$OMP END DO NOWAIT
    if (llcmoutput) then
      ! because divide grid by densityoutgrid
      densityoutgrid=1./densityoutgrid
    endif
  else
    ! no division by density
    densityoutgrid(:,:,:)=1.
  endif ! llcmoutput

  ! Output is different for forward and backward simulations
  if (ldirect.eq.1) then
!$OMP DO
    do kz=1,numzgrid
      do jy=0,numygrid-1
         do ix=0,numxgrid-1
            if (llcmoutput) then
              factor3d(ix,jy,kz)=1.e12/gridcnt(ix,jy,kz)
            else
              factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
            endif
         end do
      end do
    end do
!$OMP END DO
  else
!$OMP DO
    do kz=1,numzgrid
      do jy=0,numygrid-1
         do ix=0,numxgrid-1
            factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
         end do
      end do
    end do
!$OMP END DO
  endif

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  if (llcmoutput) then
    ks_start=2
  else
    ks_start=1
  endif

  do ks=ks_start,nspec

    do kp=1,maxpointspec_act
      do nage=1,nageclass
!$OMP DO
        do jy=0,numygrid-1
          do ix=0,numxgrid-1

            ! WET DEPOSITION
            if ((wetdep).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=wetgridunc0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=wetgridunc(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy)*real(nclassunc,kind=sp)
              wetgridtotal=wetgridtotal+wetgrid(ix,jy)
              ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
              wetgridsigmatotal=wetgridsigmatotal+ &
                   wetgridsigma(ix,jy)
            endif

            ! DRY DEPOSITION
            if ((drydep).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=drygridunc0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=drygridunc(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)*real(nclassunc,kind=sp)
              drygridtotal=drygridtotal+drygrid(ix,jy)
              ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc, kind=dep_prec))
              drygridsigmatotal=drygridsigmatotal+ &
                   drygridsigma(ix,jy)
            endif

            ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
              ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*real(nclassunc)
              gridtotal=gridtotal+grid(ix,jy,kz)
              ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
              gridsigmatotal=gridsigmatotal+ &
                   gridsigma(ix,jy,kz)
            end do
          end do
        end do
!$OMP END DO
  !       print*,gridtotal,maxpointspec_act

        !*******************************************************************
        ! Generate output: may be in concentration (ng/m3) or in mixing
        ! ratio (ppt) or both
        ! Output the position and the values alternated multiplied by
        ! 1 or -1, first line is number of values, number of positions
        ! For backward simulations, the unit is seconds, stored in grid_time
        !*******************************************************************

        ! Concentration output
        !*********************
!$OMP BARRIER
!$OMP SINGLE
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
            call nf90_err(nf90_put_var(ncid,wdspecID(ks),1.e12*&
                 wetgrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          end if

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
            call nf90_err(nf90_put_var(ncid,ddspecID(ks),1.e12*&
                 drygrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          endif

          ! Concentrations
!          call nf90_err(nf90_put_var(ncid,specID(ks),grid(0:numxgrid-1,0:numygrid-1,&
!             1:numzgrid)*factor3d(0:numxgrid-1,0:numygrid-1,1:numzgrid)/tot_mu(ks,kp),&
!               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzgrid,1,1,1 /) ))
          call nf90_err(nf90_put_var(ncid,specID(ks),grid(0:numxgrid-1,0:numygrid-1,&
             1:numzwrite)*factor3d(0:numxgrid-1,0:numygrid-1,1:numzwrite)/tot_mu(ks,kp),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzwrite,1,1,1 /) ))
 
        endif !  concentration output

        ! Mixing ratio output
        !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
            call nf90_err(nf90_put_var(ncid,wdspecID(ks),1.e12*&
                 wetgrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))

          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
            call nf90_err(nf90_put_var(ncid,ddspecID(ks),1.e12*&
                 drygrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          endif

          ! Mixing ratios
!          call nf90_err(nf90_put_var(ncid,specIDppt(ks),weightair/weightmolar(ks)*&
!               grid(0:numxgrid-1,0:numygrid-1,1:numzgrid)*&
!               factor3d(0:numxgrid-1,0:numygrid-1,1:numzgrid)/&
!               densityoutgrid(0:numxgrid-1,0:numygrid-1,1:numzgrid),&
!               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzgrid,1,1,1 /)))
          call nf90_err(nf90_put_var(ncid,specIDppt(ks),weightair/weightmolar(ks)*&
               grid(0:numxgrid-1,0:numygrid-1,1:numzwrite)*&
               factor3d(0:numxgrid-1,0:numygrid-1,1:numzwrite)/&
               densityoutgrid(0:numxgrid-1,0:numygrid-1,1:numzwrite),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzwrite,1,1,1 /)))

        endif ! output for ppt
!$OMP END SINGLE
!$OMP BARRIER
      end do ! nageclass
    end do ! maxpointspec_act

  end do ! nspec
!$OMP END PARALLEL

  if (gridtotal.gt.0.) gridtotalunc=real(gridsigmatotal/gridtotal,kind=sp)
  if (wetgridtotal.gt.0.) wetgridtotalunc=wetgridsigmatotal/ &
       wetgridtotal
  if (drygridtotal.gt.0.) drygridtotalunc=real(drygridsigmatotal/ &
       drygridtotal, kind=dep_prec)

  ! Close netCDF file
  !**************************
  call nf90_err(nf90_close(ncid))

  ! Reinitialization of grid
  !*************************
  gridunc(:,:,:,1:nspec,:,:,1:nageclass) = 0.  
  gridcnt(:,:,:) = 0.
#ifdef _OPENMP
  gridunc_omp(:,:,:,:,:,:,:,:) = 0.  
  gridcnt_omp(:,:,:,:) = 0.
#endif

end subroutine concoutput_netcdf


subroutine concoutput_nest_netcdf(itime,outnum)
  !                               i     i 
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the concentrations.               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output in *
  !                    sparse matrix format; additional output of uncertainty  *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                    are now specified in module par_mod                     *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                    in order to save disk space                             *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !     19 February 2010, Dominik Brunner, Empa: Adapted for COSMO             *
  !                                                                            *
  !     April 2013, Dominik Brunner, Empa                                      *
  !                    Adapted for netcdf output                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime           current simulation time                                    *
  ! outnum          number of samples                                          *
  !                                                                            *
  !*****************************************************************************

  use unc_mod, only: griduncn,drygriduncn,wetgriduncn,drygriduncn0,wetgriduncn0
 
  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum
  integer             :: ncid,kp,ks,kz,ix,jy,iix,jjy,kzz,ngrid
  integer             :: nage,i,l,jj
  real                :: tot_mu(maxspec,maxpointspec_act)
  real                :: halfheight,dz,dz1,dz2
  real                :: xl,yl
  real(dep_prec)      :: auxgrid(nclassunc)
  real                :: gridtotal
  real, parameter     :: weightair=28.97
  integer             :: numzwrite

  eps=nxmax/3.e5

  numzwrite=numzgrid
  if (sfc_only.eq.1 ) numzwrite=1

  ! open output file
  call nf90_err(nf90_open(trim(ncfnamen), nf90_write, ncid))

  ! write time (do not increase time counter here, done in main output domain)
  call nf90_err(nf90_put_var( ncid, timeID, itime, (/ tpointer /)))
  
  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************

  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1.0
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif

  gridtotal=0.
  !*******************************************************************
  ! Compute air density:
  ! brd134: we now take into account whether we are in the mother or in
  !    a nested domain (before only from mother domain)
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************
!$OMP PARALLEL PRIVATE(halfheight,kzz,dz1,dz2,dz,xl,yl,ngrid,iix,jjy, &
!$OMP kz,ix,jy,l,ks,kp,nage,auxgrid) REDUCTION(+:gridtotal)
!$OMP DO
  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)*0.5
    else
      halfheight=(outheight(kz)+outheight(kz-1))*0.5
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) exit
    end do
    kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2

    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy

        ngrid=0
        do jj=numbnests,1,-1
          if ( xl.gt.xln(jj)+eps .and. xl.lt.xrn(jj)-eps .and. &
                 yl.gt.yln(jj)+eps .and. yl.lt.yrn(jj)-eps ) then
            ngrid=jj
            exit 
          end if
        end do

        if (ngrid.eq.0) then
          iix=max(min(nint(xl),nxmin1),0)
          jjy=max(min(nint(yl),nymin1),0)

          densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,memind(2))*dz1+ &
             rho(iix,jjy,kzz-1,memind(2))*dz2)/dz
        else
          xl=(xl-xln(ngrid))*xresoln(ngrid)
          yl=(yl-yln(ngrid))*yresoln(ngrid)
          iix=max(min(nint(xl),nxn(ngrid)-1),0)
          jjy=max(min(nint(yl),nyn(ngrid)-1),0)
          densityoutgrid(ix,jy,kz)=(rhon(iix,jjy,kzz,memind(2), ngrid)*dz1+ &
             rhon(iix,jjy,kzz-1,memind(2), ngrid)*dz2)/dz
        endif

      end do
    end do
  end do
!$OMP END DO NOWAIT

  ! Output is different for forward and backward simulations
  if (ldirect.eq.1) then
!$OMP DO
     do kz=1,numzgrid
        do jy=0,numygridn-1
           do ix=0,numxgridn-1
              factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
           end do
        end do
     end do
!$OMP END DO
  else
!$OMP DO
     do kz=1,numzgrid
        do jy=0,numygridn-1
           do ix=0,numxgridn-1
              factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
           end do
        end do
     end do
!$OMP END DO
  endif

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

    do kp=1,maxpointspec_act
      do nage=1,nageclass
!$OMP DO
        do jy=0,numygridn-1
          do ix=0,numxgridn-1
            ! WET DEPOSITION
            if ((WETDEP).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=wetgriduncn0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=wetgriduncn(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy)*real(nclassunc)
              ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
            endif

            ! DRY DEPOSITION
            if ((DRYDEP).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=drygriduncn0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=drygriduncn(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)*real(nclassunc)
              ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
            endif

            ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
              ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*real(nclassunc)
              gridtotal=gridtotal+grid(ix,jy,kz)
              ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
            end do
          end do
        end do
!$OMP END DO
  !       print*,gridtotal,maxpointspec_act

        !*******************************************************************
        ! Generate output: may be in concentration (ng/m3) or in mixing
        ! ratio (ppt) or both
        ! Output the position and the values alternated multiplied by
        ! 1 or -1, first line is number of values, number of positions
        ! For backward simulations, the unit is seconds, stored in grid_time
        !*******************************************************************

        ! Concentration output
        !*********************
!$OMP SINGLE
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
             call nf90_err(nf90_put_var(ncid,wdspecIDn(ks),1.e12*&
                  wetgrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
             call nf90_err(nf90_put_var(ncid,ddspecIDn(ks),1.e12*&
                  drygrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Concentrations
!          call nf90_err(nf90_put_var(ncid,specIDn(ks),grid(0:numxgridn-1,0:numygridn-1,&
!             1:numzgrid)*factor3d(0:numxgridn-1,0:numygridn-1,1:numzgrid)/tot_mu(ks,kp),&
!               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzgrid,1,1,1 /)))
          call nf90_err(nf90_put_var(ncid,specIDn(ks),grid(0:numxgridn-1,0:numygridn-1,&
             1:numzwrite)*factor3d(0:numxgridn-1,0:numygridn-1,1:numzwrite)/tot_mu(ks,kp),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzwrite,1,1,1 /)))
 
        endif !  concentration output

        ! Mixing ratio output
        !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
             call nf90_err(nf90_put_var(ncid,wdspecIDn(ks),1.e12*&
                  wetgrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
             call nf90_err(nf90_put_var(ncid,ddspecIDn(ks),1.e12*&
                  drygrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Mixing ratios
!          call nf90_err(nf90_put_var(ncid,specIDnppt(ks),weightair/weightmolar(ks)*&
!               grid(0:numxgridn-1,0:numygridn-1,1:numzgrid)*&
!               factor3d(0:numxgridn-1,0:numygridn-1,1:numzgrid)/&
!               densityoutgrid(0:numxgridn-1,0:numygridn-1,1:numzgrid),&
!               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzgrid,1,1,1 /)))
          call nf90_err(nf90_put_var(ncid,specIDnppt(ks),weightair/weightmolar(ks)*&
               grid(0:numxgridn-1,0:numygridn-1,1:numzwrite)*&
               factor3d(0:numxgridn-1,0:numygridn-1,1:numzwrite)/&
               densityoutgrid(0:numxgridn-1,0:numygridn-1,1:numzwrite),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzwrite,1,1,1 /)))

        endif ! output for ppt
!$OMP END SINGLE
!$OMP BARRIER
      end do
    end do

  end do
!$OMP END PARALLEL
  ! Close netCDF file
  !**************************
  call nf90_err(nf90_close(ncid))

  ! Reinitialization of grid
  !*************************

  griduncn(:,:,:,1:nspec,:,:,1:nageclass) = 0.  

end subroutine concoutput_nest_netcdf

! subroutine concoutput_sfc_nest_netcdf(itime,outnum)

!   implicit none

!   integer, intent(in) :: itime
!   real, intent(in)    :: outnum

!   print*,'Netcdf output for surface only not yet implemented'
! end subroutine concoutput_sfc_nest_netcdf

subroutine create_particles_initialoutput(itime,idate,itime_start,idate_start)

  !*****************************************************************************
  !                                                                            *
  !   This subroutine creates an initial particle positions and properties     *
  !   NetCDF file: partinit_xxx.nc                                             *
  !   The release time, release number and positions, together with all fields *
  !   specified in the PARTOPTIONS option file will saved.                     *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in) :: itime,idate,itime_start,idate_start
  ! integer, intent(in) :: irelease
  integer             :: ncid,j,np
  integer             :: partDimID
  character(len=11)   :: fprefix
  character(len=3)    :: anspec
  character           :: adate*8,atime*6,adate_start*8,atime_start*6,timeunit*32
  character(len=255)  :: fname_partoutput
  real                :: fillval

  write(adate,'(i8.8)') idate
  write(atime,'(i6.6)') itime
  write(adate_start,'(i8.8)') idate_start
  write(atime_start,'(i6.6)') itime_start
  ! write(arelease, '(i3.3)') irelease
  fprefix = 'partinit_'!rel'//arelease//'_'

  fname_partoutput = path(2)(1:length(2))//trim(fprefix)//adate//atime//'.nc'
  !ncfname_part(irelease) = fname_partoutput
  ncfname_partinit = fname_partoutput

  call nf90_err(nf90_create(trim(fname_partoutput), cmode = nf90_hdf5, ncid = ncid))!, &
    ! cache_size = cache_size))

  ! create dimensions:
  !*************************
  
  ! particle
  partinitpointer=0
  call nf90_err(nf90_def_dim(ncid, 'particle', nf90_unlimited, partDimID))

  ! create variables
  !*************************

  ! particles
  call nf90_err(nf90_def_var(ncid, 'particle', nf90_int, (/ partDimID/), partIDi))
  call nf90_err(nf90_put_att(ncid, partIDi, 'long_name', 'particle index'))

  fillval = -1.
  ! time
  timeunit = 'seconds since '//adate_start(1:4)//'-'//adate_start(5:6)// &
     '-'//adate_start(7:8)//' '//atime_start(1:2)//':'//atime_start(3:4)

  call write_to_file(ncid,'time',nf90_int,(/ partDimID /),tIDi,(/ 1 /), &
    timeunit,.false.,'time','time of release')
  call nf90_err(nf90_put_att(ncid, tIDi, 'axis', 't'))
  call nf90_err(nf90_put_att(ncid, tIDi, 'calendar', 'proleptic_gregorian'))
  call nf90_err(nf90_put_att(ncid, tIDi, 'description', 'time of release'))

  ! lon  
  call write_to_file(ncid,'lon',nf90_float,(/ partDimID /),lonIDi,(/ 1 /), &
    'degrees_east',.false.,'longitude','longitude in degree east')
  call nf90_err(nf90_put_att(ncid, lonIDi, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'description', 'longitude of particles'))

  ! lat
  call write_to_file(ncid,'lat',nf90_float,(/ partDimID /),latIDi,(/ 1 /), &
    'degrees_north',.false.,'latitude','latitude in degree north')
  call nf90_err(nf90_put_att(ncid, latIDi, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latIDi, 'description', 'latitude of particles'))

  ! height
  call write_to_file(ncid,'z',nf90_float,(/ partDimID /),levIDi,(/ 1 /), &
   'meters',.true.,'height','height above ground')

  ! release
  call write_to_file(ncid,'release',nf90_int,(/ partDimID /),relIDi,(/ 1 /), &
   '',.true.,'release','particle release')

  do np=1,num_partopt
    if (.not. partopt(np)%print) cycle
    select case(partopt(np)%name)
      case ('PV') ! Potential vorticity
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),pvIDi,(/ 1 /), &
          'pvu',.false.,'potential_vorticity','potential vorticity')
      case ('PR') ! Pressure
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),prIDi,(/ 1 /), &
          'Pa',.false.,'pressure','pressure')
      case ('QV') ! Specific humidity
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),qvIDi,(/ 1 /), &
          '',.false.,'specific_humidity','specific humidity')
      case ('RH') ! Density
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),rhoIDi,(/ 1 /), &
          'kg/m3',.true.,'density','density')
      case ('TT') ! Temperature
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),ttIDi,(/ 1 /), &
          'K',.true.,'temperature','temperature')
      case ('UU')
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),uIDi,(/ 1 /), &
          'm/s',.false.,'u','longitudinal velocity')
      case ('VV')
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),vIDi,(/ 1 /), &
          'm/s',.false.,'v','latitudinal velocity')
      case ('WW')
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),wIDi,(/ 1 /), &
          'm/s',.false.,'w','vertical velocity')
      case ('MA')
        do j=1,nspec
          ! Masses
          write(anspec, '(i3.3)') j
          call write_to_file(ncid,trim(partopt(np)%short_name)//anspec,nf90_float,(/ partDimID /),massIDi(j), &
            (/ 1 /),'kg',.true.,'mass'//anspec,'mass for nspec'//anspec) 
        end do        
      case ('TO')
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),topoIDi,(/ 1 /), &
          'meters',.false.,'topography','topography above sealevel')
      case ('TR')
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),trIDi,(/ 1 /), &
          'meters',.true.,'htropo','height above ground of tropopause')
      case ('HM') ! Mixing layer height
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ partDimID /),hmixIDi,(/ 1 /), &
          'meters',.true.,'hmix','height above ground of mixing layer')        
      case default
        cycle
    end select
  end do

  ! moves the file from define to data mode
  call nf90_err(nf90_enddef(ncid))

  call nf90_err(nf90_close(ncid))
end subroutine create_particles_initialoutput

subroutine wrt_part_initialpos(itime,istart,iend)

  !*****************************************************************************
  !                                                                            *
  !   This subroutine saves initial particle positions, release time and       *
  !   releasenumber to a NetCDF file created in create_particles_initialoutput *
  !   evertime a new particle is spawned.                                      *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  use particle_mod

  implicit none

  integer, intent(in) ::  &
    itime,               & ! time of particle release
    istart,              & ! index of first newly released particle
    iend                   ! index of last newly released partile
  integer, allocatable    :: partindices(:),releasetimes(:)
  integer :: newpart,ncid,j

  newpart = iend-istart
  if (newpart.eq.0) return
  write(*,*) newpart, ' particles are being added to partinit.'
  call nf90_err(nf90_open(trim(ncfname_partinit), nf90_write, ncid))

  allocate ( partindices(newpart) )

  do j=1,newpart
    partindices(j)=j+partinitpointer
  end do

  partinitpointer1= partinitpointer+1 ! this is also used in partinit_netcdf
  call nf90_err(nf90_put_var(ncid,partIDi,partindices,(/ partinitpointer1 /),(/ newpart /)))
  deallocate (partindices)

  allocate ( releasetimes(newpart) )
  releasetimes=itime
  call nf90_err(nf90_put_var(ncid,tIDi,releasetimes,(/ partinitpointer1 /),(/ newpart /)))
  deallocate (releasetimes)
  call nf90_err(nf90_put_var(ncid,lonIDi,xlon0+part(partinitpointer1:iend)%xlon*dx, (/ partinitpointer1 /),(/ newpart /)))
  call nf90_err(nf90_put_var(ncid,latIDi,ylat0+part(partinitpointer1:iend)%ylat*dy, (/ partinitpointer1 /),(/ newpart /)))
  call nf90_err(nf90_put_var(ncid,levIDi,part(partinitpointer1:iend)%z, (/ partinitpointer1 /),(/ newpart /)))
  call nf90_err(nf90_put_var(ncid,relIDi,part(partinitpointer1:iend)%npoint, (/ partinitpointer1 /),(/ newpart /)))

  call nf90_err(nf90_close(ncid))

  partinitpointer = partinitpointer+newpart
end subroutine wrt_part_initialpos

subroutine partinit_netcdf(field,fieldname,imass,ncid)

  !*****************************************************************************
  !                                                                            *
  !   This subroutine saves properties chosen by the user in PARTOPTIONS       *
  !   to a NetCDF file created in create_particles_initialoutput.              *
  !   This happens whenever a new particle is spawned.                         *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in)            :: imass
  real, intent(in)               :: field(:)
  character(2), intent(in)       :: fieldname  ! input field to interpolate over
  integer                        :: ncid,newpart

  newpart = partinitpointer - (partinitpointer1-1)

  select case(fieldname)
    case('TO') ! Topography
      call nf90_err(nf90_put_var(ncid,topoIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('PV') ! Potential vorticity
      call nf90_err(nf90_put_var(ncid,pvIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('PR') ! Pressure
      call nf90_err(nf90_put_var(ncid,prIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('QV') ! Specific humidity
      call nf90_err(nf90_put_var(ncid,qvIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('RH') ! Air density
      call nf90_err(nf90_put_var(ncid,rhoIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('UU') ! Longitudinal velocity
      call nf90_err(nf90_put_var(ncid,uIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('VV') ! Latitudinal velocity
      call nf90_err(nf90_put_var(ncid,vIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('WW') ! Vertical velocity
      call nf90_err(nf90_put_var(ncid,wIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('TT') ! Temperature
      call nf90_err(nf90_put_var(ncid,ttIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('MA') ! Mass
      call nf90_err(nf90_put_var(ncid,massIDi(imass),field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('TR') ! Tropopause
      call nf90_err(nf90_put_var(ncid,trIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case('HM') ! Mixing height
      call nf90_err(nf90_put_var(ncid,hmixIDi,field(partinitpointer1:partinitpointer), &
        (/ partinitpointer1 /),(/ newpart /)))
    case default
      return
  end select
end subroutine partinit_netcdf

subroutine writeheader_partoutput(itime,idate,itime_start,idate_start)!,irelease)

  !*****************************************************************************
  !                                                                            *
  !   This subroutine creates a file (partoutput_xxx.nc), where every time     *
  !   interval particle properties specified in the PARTOPTIONS option file    *
  !   are saved to. Running options are saved as header informtion to this     *
  !   file as well.                                                            *
  !                                                                            *
  !   Author: L. Bakels 2021                                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in) :: itime,idate,itime_start,idate_start
  ! integer, intent(in) :: irelease
  integer             :: ncid,j,i,totpart,np
  integer             :: timeDimID,partDimID,tID
  integer             :: latDimID, lonDimID, lonID, latID
  character(len=255)  :: fprefix_part
  character(len=3)    :: anspec
  character           :: adate*8,atime*6,adate_start*8,atime_start*6,timeunit*32
  real, allocatable, dimension(:) :: coord

  logical,save :: first_time=.true.

  open(unit=unittmp,file=trim(path(2)(1:length(2)))//'test_dir.txt',status='replace',&
       &err=110)
  close (unittmp, status='delete')

  write(adate,'(i8.8)') idate
  write(atime,'(i6.6)') itime
  write(adate_start,'(i8.8)') idate_start
  write(atime_start,'(i6.6)') itime_start
  
  timeunit = 'seconds since '//adate_start(1:4)//'-'//adate_start(5:6)// &
       '-'//adate_start(7:8)//' '//atime_start(1:2)//':'//atime_start(3:4)

  ! write(arelease, '(i3.3)') irelease
  fprefix_part = 'partoutput_'//adate//atime !rel'//arelease//'_'

  ! Reset logicals that ensure ony 1 write out in case of domainfill
  topo_written=.false.
  mass_written=.false.
  massav_written=.false.

  totpart=0
  if (ipin.gt.1) then ! Not reading from a release has no npart
    totpart=numpart
  else
    do j=1,numpoint
      totpart = totpart+npart(j)
    end do
  endif
  !totpart = maxpart!max(numpart,totpart)
  !cache_size = 4 * 1 * (12+nspec)
  ncfname_part = path(2)(1:length(2))//trim(fprefix_part)
  if (lpartoutputperfield) then
    do np=1,num_partopt
      if (.not. partopt(np)%print ) cycle
      if (first_time) then
        call nf90_err(nf90_create(trim(ncfname_part)//'_'//trim(partopt(np)%long_name)//'_init.nc', &
          cmode = nf90_hdf5, ncid = partopt(np)%ncid))
        ncfname_part_end = '_init.nc'
      else
        call nf90_err(nf90_create(trim(ncfname_part)//'_'//trim(partopt(np)%long_name)//'.nc', &
          cmode = nf90_hdf5, ncid = partopt(np)%ncid))
        ncfname_part_end = '.nc'
      endif
    end do
    first_time=.false.
  else 
    if (first_time) then
      ncfname_part = path(2)(1:length(2))//trim(fprefix_part)//'_init.nc'
      first_time=.false.
    else
      ncfname_part = path(2)(1:length(2))//trim(fprefix_part)//'.nc'
    endif
    call nf90_err(nf90_create(trim(ncfname_part), cmode = nf90_hdf5, ncid = ncid))!, &
      ! cache_size = cache_size))    
  endif

  write(*,*) 'Write header, nspec,numpart,totpart: ', nspec,numpart,totpart

  if (lpartoutputperfield) then
    do np=1,num_partopt
      if (.not. partopt(np)%print) cycle
      call writeheader_partoutput_dims(np,partopt(np)%ncid,timeunit,timeDimID,partDimID,latDimID,lonDimID)
      call writeheader_partoutput_vars(np,partopt(np)%ncid,totpart,timeDimID,partDimID,latDimID,lonDimID)

      ! moves the file from define to data mode
      call nf90_err(nf90_enddef(partopt(np)%ncid))
      call nf90_err(nf90_close(partopt(np)%ncid))
    end do
  else
    call writeheader_partoutput_dims(1,ncid,timeunit,timeDimID,partDimID,latDimID,lonDimID)
    do np=1,num_partopt
      if (.not. partopt(np)%print) cycle
      call writeheader_partoutput_vars(np,ncid,totpart,timeDimID,partDimID,latDimID,lonDimID)
    end do

    ! moves the file from define to data mode
    call nf90_err(nf90_enddef(ncid))
    call nf90_err(nf90_close(ncid))
  endif

  return
110 write(*,FMT='(80("#"))') 
  write(*,*) 'ERROR: output directory ', trim(path(2)(1:length(2))), ' does not exist&
       & (or failed to write there).' 
  write(*,*) 'EXITING' 
  write(*,FMT='(80("#"))')
  error stop
end subroutine writeheader_partoutput

subroutine writeheader_partoutput_dims(np,ncid,timeunit,timeDimID,partDimID,latDimID,lonDimID)

  implicit none
  integer,intent(in)  :: ncid,np
  character,intent(in) :: timeunit*32
  integer,intent(out) :: timeDimID,partDimID
  integer,intent(out) :: latDimID, lonDimID
  integer             :: tID,partID

  logical,save :: first_time=.true.

  ! create dimensions:
  !*************************
  ! time
  call nf90_err(nf90_def_dim(ncid, 'time', nf90_unlimited, timeDimID))

  ! particle

  ! If domainfill, save topo, hmix, and htropo to grid to save space
  !*****************************************************************
  if (lpartoutputperfield.and.(mdomainfill.eq.1).and. &
    ((partopt(np)%name.eq.'TO') .or. &
    (partopt(np)%name.eq.'HM') .or. &
    (partopt(np)%name.eq.'TR'))) then
    call writeheader_partoutput_grid(ncid,lonDimID,latDimID)
  else 
    if (.not. lpartoutputperfield .and. (mdomainfill.eq.1)) then
      call writeheader_partoutput_grid(ncid,lonDimID,latDimID)
    endif
    call nf90_err(nf90_def_dim(ncid, 'particle', nf90_unlimited, partDimID))
    ! particles variables
    call nf90_err(nf90_def_var(ncid, 'particle', nf90_int, (/ partDimID/), partID))
    call nf90_err(nf90_put_att(ncid, partID, 'long_name', 'particle index'))
  endif
  ! create variables
  !*************************

  ! time
  tpointer_part=0
  call nf90_err(nf90_def_var(ncid, 'time', nf90_int, (/ timeDimID /), tID))
  call nf90_err(nf90_put_att(ncid, tID, 'units', timeunit))
  call nf90_err(nf90_put_att(ncid, tID, 'calendar', 'proleptic_gregorian'))

  timeIDpart=tID

  ! global (metadata) attributes
  !*******************************
  call writemetadata(ncid,lnest=.false.)

end subroutine writeheader_partoutput_dims

subroutine writeheader_partoutput_vars(np,ncid,totpart,timeDimID,partDimID,latDimID,lonDimID)

  implicit none
  integer,intent(in)  :: ncid,totpart,np
  integer,intent(in)  :: timeDimID,partDimID
  integer,intent(in)  :: latDimID, lonDimID
  integer             :: j,i,varid
  character(len=3)    :: anspec
  real                :: fillval

  fillval = -1.
  select case(partopt(np)%name)
    case ('LO') ! Longitude
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'degrees_east',.false.,'longitude','longitude of particles')
      call nf90_err(nf90_put_att(ncid, varid, 'axis', 'Lon'))
      call nf90_err(nf90_put_att(ncid, varid, 'description', 'longitude of particles'))
    case ('lo') ! Longitude averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'degrees_east',.false.,'longitude_average','averaged longitude of particles')
      call nf90_err(nf90_put_att(ncid, varid, 'axis', 'Lon'))
      call nf90_err(nf90_put_att(ncid, varid, 'description', 'averaged longitude of particles'))
    case ('LA') ! Latitude
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'degrees_north',.false.,'latitude','latitude in degree north')
      call nf90_err(nf90_put_att(ncid, varid, 'axis', 'Lat'))
      call nf90_err(nf90_put_att(ncid, varid, 'description', 'latitude of particles'))
    case ('la') ! Latitude averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'degrees_north',.false.,'latitude_average','averaged latitude in degree north')
      call nf90_err(nf90_put_att(ncid, varid, 'axis', 'Lat'))
      call nf90_err(nf90_put_att(ncid, varid, 'description', 'averaged latitude of particles'))
    case ('ZZ') ! Height
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'meters',.false.,'height','height above ground')
    case ('zz') ! Heights averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'meters',.false.,'height_average','averaged height above ground')
    case ('PV') ! Potential vorticity
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'pvu',.false.,'potential_vorticity','potential vorticity')
    case ('pv') ! Potential vorticity averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'pvu',.false.,'potential_vorticity_average','averaged potential vorticity')
    case ('PR') ! Pressure
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'Pa',.false.,'pressure','pressure')
    case ('pr') ! Pressure averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'Pa',.false.,'pressure_average','averaged pressure')
    case ('QV') ! Specific humidity
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'kg/kg',.false.,'specific_humidity','specific humidity')
    case ('qv') ! Specific humidity averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'kg/kg',.false.,'specific_humidity_average','averaged specific humidity')
    case ('RH') ! Density
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'kg/m3',.true.,'density','density')
    case ('rh') ! Density averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'kg/m3',.true.,'density_average','averaged density')
    case ('TT') ! Temperature
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'K',.true.,'temperature','temperature') 
    case ('tt') ! Temperature averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'K',.true.,'temperature_average','averaged temperature') 
    case ('UU')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'u','longitudinal velocity')    
    case ('uu')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'u_av','averaged longitudinal velocity')
    case ('VV')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'v','latitudinal velocity')
    case ('vv')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'v_average','latitudinal velocity averaged')
    case ('WW')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'w','vertical velocity')
    case ('ww')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'w_average','vertical velocity averaged')
    case ('VS')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'settling_velocity','settling velocity')
    case ('vs')
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /), &
        varid,(/ 1,totpart /),'m/s',.false.,'settling_velocity_average','settling velocity averaged')
    case ('MA') ! Mass
      if ((mdomainfill.ge.1).and.(nspec.eq.1)) then
        call nf90_err(nf90_def_var(ncid=ncid, name=trim(partopt(np)%short_name), xtype=nf90_float, &
          dimids=1, varid=varid))
        call nf90_err(nf90_put_att(ncid, varid, 'units', 'kg'))
        call nf90_err(nf90_put_att(ncid, varid, '_FillValue', fillval))
        call nf90_err(nf90_put_att(ncid, varid, 'positive', 'up'))
        call nf90_err(nf90_put_att(ncid, varid, 'standard_name', 'mass'))
        call nf90_err(nf90_put_att(ncid, varid, 'long_name', 'mass of each particle'))
      else
        do j=1,nspec
          ! Masses
          write(anspec, '(i3.3)') j
          call write_to_file(ncid,trim(partopt(np)%short_name)//anspec,nf90_float, &
            (/ timeDimID,partDimID /),varid, &
            (/ 1,totpart /),'kg',.true.,'mass'//anspec,'mass for nspec'//anspec)
        end do
      endif
    case ('ma') ! Mass averaged
      if ((mdomainfill.ge.1).and.(nspec.eq.1)) then
        call nf90_err(nf90_def_var(ncid=ncid, name=trim(partopt(np)%short_name), xtype=nf90_float, dimids=1, varid=varid))
        call nf90_err(nf90_put_att(ncid, varid, 'units', 'kg'))
        call nf90_err(nf90_put_att(ncid, varid, '_FillValue', fillval))
        call nf90_err(nf90_put_att(ncid, varid, 'positive', 'up'))
        call nf90_err(nf90_put_att(ncid, varid, 'standard_name', 'mass'))
        call nf90_err(nf90_put_att(ncid, varid, 'long_name', 'averaged mass of each particle'))
      else
        do j=1,nspec
          ! Masses averaged
          write(anspec, '(i3.3)') j
          call write_to_file(ncid,trim(partopt(np)%short_name)//anspec,nf90_float,(/ timeDimID,partDimID /),varid, &
            (/ 1,totpart /),'kg',.true.,'mass'//anspec,'averaged mass for nspec'//anspec) 
        end do
      endif
    case ('WD') ! Cumulative mass of wet deposition
      do j=1,nspec
        ! Masses
        write(anspec, '(i3.3)') j
        call write_to_file(ncid,trim(partopt(np)%short_name)//anspec,nf90_float,(/ timeDimID,partDimID /),varid, &
          (/ 1,totpart /),'kg',.true.,'mass'//anspec,'cumulative wet deposition for nspec'//anspec) 
      end do
    case ('DD') ! Cumulative mass of dry deposition
      do j=1,nspec
        ! Masses
        write(anspec, '(i3.3)') j
        call write_to_file(ncid,trim(partopt(np)%short_name)//anspec,nf90_float,(/ timeDimID,partDimID /),varid, &
          (/ 1,totpart /),'kg',.true.,'mass'//anspec,'cumulative dry deposition for nspec'//anspec) 
      end do
    case ('TO')  ! Topography, written to grid if domainfill
      if (mdomainfill.lt.1) then
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
          'meters',.false.,'topography','topography above sealevel')
      else
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ lonDimID,latDimID /),varid,(/ nx,ny /), &
          'meters',.false.,'topography','topography above sealevel')
      endif
    case ('to') ! Topography averaged, no grid when domainfill
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
        'meters',.false.,'topography','averaged topography above sealevel')
    case ('HM') ! Mixing layer height
      if (mdomainfill.lt.1) then
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
          'meters',.true.,'hmix','height above ground of mixing layer')
      else
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,lonDimID,latDimID /),varid,(/ 1,nx,ny /), &
          'meters',.true.,'hmix','height above ground of mixing layer')  
      endif
    case ('hm') ! Mixing layer height averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
        'meters',.true.,'hmix_average','averaged height above ground of mixing layer')
    case ('TR') ! Tropopause
      if (mdomainfill.lt.1) then
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
          'meters',.true.,'htropo','height above ground of tropopause')
      else
        call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,lonDimID,latDimID /),varid,(/ 1,nx,ny /), &
          'meters',.true.,'htropo','height above ground of tropopause')
      endif
    case ('tr') ! Tropopause averaged
      call write_to_file(ncid,trim(partopt(np)%short_name),nf90_float,(/ timeDimID,partDimID /),varid,(/ 1,totpart /), &
        'meters',.true.,'htropo_average','averaged height above ground of tropopause')
    case default
      write(*,*) 'The field you are trying to write to file is not coded in yet: ', partopt(np)%name,partopt(np)%long_name
      error stop
  end select

end subroutine writeheader_partoutput_vars

subroutine writeheader_partoutput_grid(ncid,lonDimID,latDimID)

  implicit none

  integer,intent(in)  :: ncid
  integer,intent(out) :: lonDimID,latDimID
  real, allocatable, dimension(:) :: coord
  integer :: lonID, latID, i

  call nf90_err(nf90_def_dim(ncid, 'longitude', nx, lonDimID))
  call nf90_err(nf90_def_dim(ncid, 'latitude', ny, latDimID))

  ! lon
  call write_to_file(ncid,'longitude',nf90_float,(/ lonDimID /),lonID,(/ 1 /), &
    'degrees_east',.false.,'grid_longitude','longitude in degree east')
  call nf90_err(nf90_put_att(ncid, lonID, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonID, 'description', 'grid cell centers'))

  ! lat
  call write_to_file(ncid,'latitude',nf90_float,(/ latDimID /),latID,(/ 1 /), &
    'degrees_north',.false.,'grid_latitude','latitude in degree north')
  call nf90_err(nf90_put_att(ncid, latID, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latID, 'description', 'grid cell centers'))

  if (.not.allocated(coord)) allocate(coord(nx))
  do i = 1,nx
    coord(i) = xlon0 + (i-1)*dx
  enddo
  call nf90_err(nf90_put_var(ncid, lonID, coord(1:nx)))
  deallocate(coord)

  if (.not.allocated(coord)) allocate(coord(ny))
  do i = 1,ny
    coord(i) = ylat0 + (i-1)*dy
  enddo
  call nf90_err(nf90_put_var(ncid, latID, coord(1:ny)))
  deallocate(coord)

end subroutine writeheader_partoutput_grid

subroutine write_to_file(ncid,short_name,xtype,dimids,varid,chunksizes,units,l_positive, &
  standard_name,long_name)

  !*****************************************************************************
  !                                                                            *
  !   Generalised writing data to netcdf file                                  *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in) :: ncid, xtype
  integer, intent(out) :: varid
  character(len = *), intent(in) :: short_name,standard_name,long_name,units
  integer, dimension(:), intent(in) :: dimids,chunksizes
  logical, intent(in) :: l_positive

  call nf90_err(nf90_def_var(ncid, short_name, xtype, dimids, varid))
  call nf90_err(nf90_def_var_chunking(ncid,varid,NF90_CHUNKED,chunksizes=chunksizes))
  call nf90_err(nf90_def_var_deflate(ncid,varid,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, varid, 'units', units))
  if(xtype.eq.nf90_float) then
    call nf90_err(nf90_put_att(ncid, varid, '_FillValue', -1.))
  else 
    call nf90_err(nf90_put_att(ncid, varid, '_FillValue', -1))
  endif
  if(l_positive) call nf90_err(nf90_put_att(ncid, varid, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, varid, 'standard_name', standard_name))
  call nf90_err(nf90_put_att(ncid, varid, 'long_name', long_name))
end subroutine write_to_file

subroutine open_partoutput_file(ncid,np)
  
  implicit none 

  integer, intent(out)               :: ncid
  integer, intent(in),optional         :: np

  if (lpartoutputperfield) then
    call nf90_err(nf90_open(trim(ncfname_part)//'_'//trim(partopt(np)%long_name)//trim(ncfname_part_end), &
      nf90_write, ncid))
  else
    call nf90_err(nf90_open(trim(ncfname_part), nf90_write, ncid))
  endif
end subroutine open_partoutput_file

subroutine close_partoutput_file(ncid)
  
  implicit none 

  integer                        :: ncid

  call nf90_err(nf90_close(ncid))
end subroutine close_partoutput_file

subroutine open_partinit_file(ncid)
  
  implicit none 

  integer, intent(inout)         :: ncid

  call nf90_err(nf90_open(trim(ncfname_partinit), nf90_write, ncid))
end subroutine open_partinit_file

subroutine update_partoutput_pointers(itime,ncid)

  use particle_mod

  implicit none

  integer, intent(in)         :: itime
  integer, intent(in),optional :: ncid
  integer, allocatable           :: partindices(:)
  integer                     :: j,tempIDend,newpart,np

  ! Time
  tpointer_part = tpointer_part + 1

  if (lpartoutputperfield) then
    do np=1,num_partopt
      if (.not. partopt(np)%print) cycle
      call nf90_err(nf90_inq_varid(ncid=partopt(np)%ncid,name='time',varid=tempIDend))
      call nf90_err(nf90_put_var(partopt(np)%ncid, tempIDend, itime, (/ tpointer_part /)))
    end do
  else
    call nf90_err(nf90_inq_varid(ncid=ncid,name='time',varid=tempIDend))
    call nf90_err(nf90_put_var(ncid, tempIDend, itime, (/ tpointer_part /)))
  endif

  ! Particles
  newpart = count%allocated - ppointer_part

  if (tpointer_part.eq.1) then 
    allocate ( partindices(count%allocated) )
    do j=1,count%allocated 
      partindices(j)=j
    end do 
    if (lpartoutputperfield) then
      do np=1,num_partopt
        if (.not. partopt(np)%print) cycle
        if ((mdomainfill.eq.1).and. &
          ((partopt(np)%name.eq.'TO') .or. &
          (partopt(np)%name.eq.'HM') .or. &
          (partopt(np)%name.eq.'TR'))) cycle
        call nf90_err(nf90_inq_varid(ncid=partopt(np)%ncid,name='particle',varid=tempIDend))
        call nf90_err(nf90_put_var(partopt(np)%ncid, tempIDend,partindices, (/ 1 /),(/ count%allocated /)))
      end do
    else
      call nf90_err(nf90_inq_varid(ncid=ncid,name='particle',varid=tempIDend))
      call nf90_err(nf90_put_var(ncid, tempIDend,partindices, (/ 1 /),(/ count%allocated /)))
    endif
    deallocate (partindices)

    ppointer_part = count%allocated

  else if (newpart.ge.0) then

    allocate ( partindices(newpart) )
    do j=1,newpart
      partindices(j)=j+ppointer_part
    end do
    if (lpartoutputperfield) then
      do np=1,num_partopt
        if (.not. partopt(np)%print) cycle
        if ((mdomainfill.eq.1).and. &
          ((partopt(np)%name.eq.'TO') .or. &
          (partopt(np)%name.eq.'HM') .or. &
          (partopt(np)%name.eq.'TR'))) cycle
        call nf90_err(nf90_inq_varid(ncid=partopt(np)%ncid,name='particle',varid=tempIDend))
        call nf90_err(nf90_put_var(partopt(np)%ncid, tempIDend,partindices, (/ ppointer_part+1 /),(/ newpart /)))
      end do
    else
      call nf90_err(nf90_inq_varid(ncid=ncid,name='particle',varid=tempIDend))
      call nf90_err(nf90_put_var(ncid, tempIDend,partindices, (/ ppointer_part+1 /),(/ newpart /)))
    endif
    deallocate (partindices)

    ppointer_part = count%allocated
  endif 

end subroutine update_partoutput_pointers

subroutine partoutput_netcdf(itime,field,np,imass,ncid)

  
  use particle_mod
  !*****************************************************************************
  !                                                                            *
  !   Writing a field from PARTOPTIONS to partoutput_xxx.nc created in         *
  !   writeheader_partoutput                                                   *
  !                                                                            *  
  !   Author: L. Bakels 2021                                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer, intent(in)            :: itime,imass,ncid
  real, intent(in)               :: field(:)
  integer, intent(in)            :: np  ! input field to interpolate over
  integer                        :: tempIDend
  character(len=3)               :: anspec

  ! ! open output file
  ! call nf90_err(nf90_open(trim(ncfname_part), nf90_write, ncid))
  if ((mdomainfill.ge.1).and. ((partopt(np)%name.eq.'TO').or. &
    (partopt(np)%name.eq.'HM').or.(partopt(np)%name.eq.'TR'))) then
    if (partopt(np)%name.eq.'TO')  then
      if (topo_written.eqv..false.) then
        call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
        call nf90_err(nf90_put_var(ncid,tempIDend,oro(0:nx-1,0:ny-1), (/ 1,1 /),(/ nx,ny /)))
        topo_written=.true.
      endif
    else if (partopt(np)%name.eq.'HM') then !HM
      call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
      call nf90_err(nf90_put_var(ncid,tempIDend,hmix(0:nx-1,0:ny-1,1,memind(1)), &
        (/ tpointer_part,1,1 /),(/ 1,nx,ny /)))
    else !TR
      call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
      call nf90_err(nf90_put_var(ncid,tempIDend,tropopause(0:nx-1,0:ny-1,1,memind(1)), &
        (/ tpointer_part,1,1 /),(/ 1,nx,ny /)))      
    endif

  else if (partopt(np)%name.eq.'MA') then
    if ((mdomainfill.ge.1).and.(imass.eq.1).and.(nspec.eq.1)) then
      if (mass_written.eqv..false.) then 
        call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
        call nf90_err(nf90_put_var(ncid=ncid,varid=tempIDend,values=field(1)))
      endif
      mass_written=.true.
    else
      write(anspec, '(i3.3)') imass
      call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name)//anspec,varid=tempIDend))
      call nf90_err(nf90_put_var(ncid,tempIDend,field, (/ tpointer_part,1 /),(/ 1,count%allocated /)))
    endif
  else if (partopt(np)%name.eq.'ma') then
    if ((mdomainfill.ge.1).and.(imass.eq.1).and.(nspec.eq.1)) then
      if (mass_written.eqv..false.) then 
        call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
        call nf90_err(nf90_put_var(ncid,tempIDend,field, (/ tpointer_part,1 /),(/ 1,count%allocated /)))
      endif
      massav_written=.true.
    else
      write(anspec, '(i3.3)') imass
      call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name)//anspec,varid=tempIDend))
      call nf90_err(nf90_put_var(ncid,tempIDend,field, (/ tpointer_part,1 /),(/ 1,count%allocated /)))
    endif
  else if ((partopt(np)%name.eq.'WD').or.(partopt(np)%name.eq.'DD')) then
    write(anspec, '(i3.3)') imass
    call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name)//anspec,varid=tempIDend))
    call nf90_err(nf90_put_var(ncid,tempIDend,field, (/ tpointer_part,1 /),(/ 1,count%allocated /)))
  else 
    call nf90_err(nf90_inq_varid(ncid=ncid,name=trim(partopt(np)%short_name),varid=tempIDend))
    call nf90_err(nf90_put_var(ncid,tempIDend,field, (/ tpointer_part,1 /),(/ 1,count%allocated /)))
  endif

  ! call nf90_err(nf90_close(ncid))
end subroutine partoutput_netcdf

subroutine readpartpositions_netcdf(ibtime,ibdate)

  !*****************************************************************************
  !                                                                            *
  !   IPIN=2: restarting from a partoutput_xxx.nc file written by a previous   *
  !           run, depending on what PARTOPTIONS the user has chosen, this     *
  !           option might not be possible to use                              *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  use random_mod
  use particle_mod
  use date_mod

  implicit none 

  integer, intent(in) :: ibtime,ibdate
  integer             :: ncidend,tIDend,pIDend,tempIDend
  integer             :: tlen,plen,tend,i,j,stat,iterminate
  integer             :: idate_start,itime_start
  character           :: adate*8,atime*6,timeunit*32,adate_start*8,atime_start*6
  character(len=3)    :: anspec
  real(kind=dp)       :: julin,julcommand
  ! real,allocatable,dimension(:) :: mass_temp
  integer :: idummy = -8

  write(adate,'(i8.8)') ibdate
  write(atime,'(i6.6)') ibtime
  
  if (mquasilag.ne.0) then 
    write(*,*) 'Combination of ipin, netcdf partoutput, and mquasilag!=0 does not work yet'
    error stop 
  endif

  ! Open partoutput_end.nc file
  call nf90_err(nf90_open(path(2)(1:length(2))//trim('partoutput_end.nc'), mode=NF90_NOWRITE,ncid=ncidend))

  ! Take the positions of the particles at the last timestep in the file
  ! It needs to be the same as given in the COMMAND file, this is arbitrary
  ! and should be removed in the future for easier use

  ! First get the time dimension
  call nf90_err(nf90_inq_dimid(ncid=ncidend,name='time',dimid=tIDend))
  call nf90_err(nf90_inquire_dimension(ncid=ncidend,dimid=tIDend,len=tlen))

  ! Check if the time corresponds to the one given in the COMMAND file
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='time',varid=tIDend))
  call nf90_err(nf90_get_att(ncid=ncidend,varid=tIDend,name='units',values=timeunit))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tIDend,values=tend,start=(/ tlen /)))!,count=(/ 1 /)))
  adate_start(1:4) = timeunit(15:18)
  adate_start(5:6) = timeunit(20:21)
  adate_start(7:8) = timeunit(23:24)
  atime_start = '000000'
  atime_start(1:2) = timeunit(26:27)
  atime_start(3:4) = timeunit(29:30)
  read(adate_start,*) idate_start
  read(atime_start,*) itime_start
  julin = juldate(idate_start,itime_start)+real(tend,kind=dp)/86400._dp
  julcommand = juldate(ibdate,ibtime)
  if (abs(julin-julcommand).gt.1.e-5) then 
    write(*,*) 'ERROR: The given starting time and date do not correspond to'
    write(*,*) 'the last timestep of partoutput_end.nc:'
    write(*,*) julin,julcommand,tend
    error stop 
  endif

  !! testing
!  print*, 'readpartpositions_netcdf: julin, julcommand = ',julin, julcommand

  ! Then the particle dimension
  call nf90_err(nf90_inq_dimid(ncid=ncidend,name='particle',dimid=pIDend))
  call nf90_err(nf90_inquire_dimension(ncid=ncidend,dimid=pIDend,len=plen))

  ! Now spawn the correct number of particles
  write(*,*) 'Npart:',plen
  call spawn_particles(0,plen)

  ! And give them the correct positions
  ! Longitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='lon',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%xlon, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))
  part(:)%xlon=(part(:)%xlon-xlon0)/dx
  ! Latitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='lat',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%ylat, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))
  part(:)%ylat=(part(:)%ylat-ylat0)/dx
  ! Height
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='z',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%z, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))
  ! Mass
  ! allocate(mass_temp(count%allocated), stat=stat)
  ! if (stat.ne.0) error stop "Could not allocate mass_temp"
  if ((mdomainfill.eq.0).or.(nspec.gt.1)) then
    do j=1,nspec
      write(anspec, '(i3.3)') j
      call nf90_err(nf90_inq_varid(ncid=ncidend,name='m'//anspec,varid=tempIDend))
      call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=mass(:,j), & 
        start=(/ tlen, 1 /),count=(/ 1, plen /)))
      ! do i=1,count%allocated
      !   part(i)%mass(j)=mass_temp(i)
      ! end do
    end do
  endif
  ! deallocate( mass_temp )

  iterminate=0
  do i=1,plen
    if (part(i)%z.lt.0) then 
      call terminate_particle(i,0)
      if (mdomainfill.eq.0) then
        write(*,*) 'Particle ',i,'is not alive in the restart file.'
      endif
      iterminate=iterminate+1
    endif
    part(i)%nclass=min(int(ran1(idummy,0)*real(nclassunc))+1, &
         nclassunc)    
    part(i)%idt=mintime
    part(i)%npoint=1
  end do

  if (iterminate.gt.0) call rewrite_ialive()
  
  call nf90_err(nf90_close(ncidend))

  !! testing
!  print*, 'readpartpositions_netcdf: number alive = ',count%alive
!  print*, 'readpartpositions_netcdf: range(part%z) = ',minval(part(1:count%alive)%z),maxval(part(1:count%alive)%z)
!  print*, 'readpartpositions_netcdf: part(1)%tstart = ',part(1)%tstart

end subroutine readpartpositions_netcdf

subroutine readinitconditions_netcdf()

  !*****************************************************************************
  !                                                                            *
  !   IPIN=3: starting a run from a user defined initial particle conditions,  *
  !           more on how to create such a file can be found in the manual     *
  !   IPIN=4: restarting a run, while also reading in the initial particle     *
  !           conditions                                                       *
  !                                                                            *
  !   Author: L. Bakels 2022                                                   *
  !                                                                            *
  !*****************************************************************************

  use random_mod
  use particle_mod
  use date_mod
  use readoptions_mod
  use drydepo_mod

  implicit none 

  integer             :: ncidend,pIDend,tempIDend,stat
  integer             :: plen,i,j,release_max,nsp
  integer(kind=2)     :: zkind
  real                :: cun
  integer,allocatable, dimension (:) :: specnum_rel,numpoint_max
  ! real,allocatable,dimension(:,:) :: mass_temp
  real,allocatable,dimension(:) :: vsh,fracth,schmih
  logical :: lstart=.true.

  integer :: idummy = -8
  
  if (mquasilag.ne.0) then 
    write(*,*) 'Combination of ipin, netcdf partoutput, and mquasilag!=0 does not work yet'
    error stop 
  endif

  ! Open part_ic.nc file
  call nf90_err(nf90_open(trim(path(1)(1:length(1))//'part_ic.nc'), mode=NF90_NOWRITE,ncid=ncidend))

  ! How many species are contained in each particle?
  call nf90_err(nf90_inquire_attribute(ncid=ncidend,name='nspecies',varid=NF90_GLOBAL))
  call nf90_err(nf90_get_att(ncid=ncidend,varid=NF90_GLOBAL,name='nspecies',values=nspec))

  maxspec=nspec
  call alloc_com()

  ! Read number of fields that need to be output. This needs to happen after maxspec is defined
  ! but before particles are allocated (n_average is necessary).
  if (ipout.ne.0) call readpartoptions 

  ! allocate with maxspec for first input loop
  allocate(specnum_rel(maxspec),stat=stat)
  if (stat.ne.0) error stop 'ERROR: could not allocate specnum_rel'

  if (nspec.gt.maxspec) then
    error stop 'number of species in part_ic.nc is larger than the allowed maxspec set in the par_mod.f90'
  endif
  ! Which species?
  call nf90_err(nf90_inquire_attribute(ncid=ncidend,name='species',varid=NF90_GLOBAL))
  call nf90_err(nf90_get_att(ncid=ncidend,varid=NF90_GLOBAL,name='species',values=specnum_rel(1:nspec)))

  ! Read species and derive initial conditions
  !****************************************************
  DEP=.false.
  DRYDEP=.false.
  WETDEP=.false.
  CLREA=.false.
  do nsp=1,maxspec
    DRYDEPSPEC(nsp)=.false.
    WETDEPSPEC(nsp)=.false.
  end do

  do nsp=1,nspec
    call readspecies(specnum_rel(nsp),nsp)
  end do

  ! Allocate fields that depend on ndia
  call alloc_com_ndia

  do nsp=1,nspec
    ! Allocate temporary memory necessary for the different diameter bins
    !********************************************************************
    allocate(vsh(ndia(nsp)),fracth(ndia(nsp)),schmih(ndia(nsp)), stat=stat)
    if (stat.ne.0) error stop "Could not allocate vsh,fracth,schmih"

    ! Molecular weight
    !*****************
    if (((iout.eq.2).or.(iout.eq.3)).and.(weightmolar(nsp).lt.0.)) then
      write(*,*) 'For mixing ratio output, valid molar weight'
      write(*,*) 'must be specified for all simulated species.'
      write(*,*) 'Check table SPECIES or choose concentration'
      write(*,*) 'output instead if molar weight is not known.'
      error stop
    endif

    ! Radioactive decay
    !******************
    decay(nsp)=0.693147/decay(nsp) !conversion half life to decay constant

  ! Dry deposition of gases
  !************************

    if (reldiff(nsp).gt.0.) rm(nsp)=1./(henry(nsp)/3000.+100.*f0(nsp))    ! mesophyll resistance

  ! Dry deposition of particles
  !****************************

    vsetaver(nsp)=0.
    cunningham(nsp)=0.
    dquer(nsp)=dquer(nsp)*1000000.         ! Conversion m to um
    if (density(nsp).gt.0.) then         ! Additional parameters
      call part0(dquer(nsp),dsigma(nsp),density(nsp),ndia(nsp),fracth,schmih,cun,vsh)
      do j=1,ndia(nsp)
        fract(nsp,j)=fracth(j)
        schmi(nsp,j)=schmih(j)
        vset(nsp,j)=vsh(j)
        cunningham(nsp)=cunningham(nsp)+cun*fract(nsp,j)
        vsetaver(nsp)=vsetaver(nsp)-vset(nsp,j)*fract(nsp,j)
      end do
      if (lroot) write(*,*) 'Average settling velocity: ',i,vsetaver(nsp)
    endif

    ! Dry deposition for constant deposition velocity
    !************************************************

    dryvel(nsp)=dryvel(nsp)*0.01         ! conversion to m/s

    ! Check if wet deposition or OH reaction shall be calculated
    !***********************************************************

    ! ESO 04.2016 check for below-cloud scavenging (gas or aerosol)
    if ((dquer(nsp).le.0..and.(weta_gas(nsp).gt.0. .or. wetb_gas(nsp).gt.0.)) .or. &
         &(dquer(nsp).gt.0. .and. (crain_aero(nsp) .gt. 0. .or. csnow_aero(nsp).gt.0.)))  then
      WETDEP=.true.
      WETDEPSPEC(nsp)=.true.
      if (lroot) then
        write (*,*) '  Below-cloud scavenging: ON'
      end if
    else
      if (lroot) write (*,*) '  Below-cloud scavenging: OFF'
    endif

    ! NIK 31.01.2013 + 10.12.2013 + 15.02.2015
    if (dquer(nsp).gt.0..and.(ccn_aero(nsp).gt.0. .or. in_aero(nsp).gt.0.))  then
      WETDEP=.true.
      WETDEPSPEC(nsp)=.true.
      if (lroot) then
        write (*,*) '  In-cloud scavenging: ON'
      end if
    else
      if (lroot) write (*,*) '  In-cloud scavenging: OFF' 
    endif

    if (any(reaccconst(:,:).gt.0.)) then
      CLREA=.true.
      if (lroot) write (*,*) '  Chemical reactions switched on'
    endif

    if ((reldiff(nsp).gt.0.).or.(density(nsp).gt.0.).or.(dryvel(nsp).gt.0.)) then
      DRYDEP=.true.
      DRYDEPSPEC(nsp)=.true.
    endif

    deallocate(vsh,fracth,schmih)
  end do ! end loop over species

  if (WETDEP.or.DRYDEP) then 
   DEP=.true.
  endif

  deallocate(specnum_rel)
  !********************************* END READING SPECIES

  ! Get the particle dimension
  call nf90_err(nf90_inq_dimid(ncid=ncidend,name='particle',dimid=pIDend))
  call nf90_err(nf90_inquire_dimension(ncid=ncidend,dimid=pIDend,len=plen))

  ! Now spawn the correct number of particles
  write(*,*) 'Npart:',plen
  call alloc_particles( plen )
  ! allocate temporary mass array
  ! allocate(mass_temp(plen,nspec))

  ! And give them the correct positions
  ! Longitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='longitude',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%xlon, & 
    start=(/ 1 /),count=(/ plen /)))
  part(:)%xlon=(part(:)%xlon-xlon0)/dx
  ! Latitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='latitude',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%ylat, & 
    start=(/ 1 /),count=(/ plen /)))
  part(:)%ylat=(part(:)%ylat-ylat0)/dx
  ! Height
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='height',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%z, & 
    start=(/ 1 /),count=(/ plen /)))
  ! Spawning time
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='time',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%tstart, & 
    start=(/ 1 /),count=(/ plen /)))
  ! Mass
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='mass',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=mass, & 
    start=(/ 1,1 /),count=(/ plen,nspec /)))
  ! do i=1,plen
  !   do nsp=1,nspec
  !     part(i)%mass(nsp)=mass_temp(i,nsp)
  !   end do
  ! end do
  ! deallocate(mass_temp)

  ! Check if they are within the bounds
  do i=1,plen
    if ((part(i)%xlon.lt.0).or.(part(i)%xlon.gt.nx)) then
      write(*,*) 'Dimensions (nx,ny): ',nx,ny
      write(*,*) "Particle", i, "with xlon", part(i)%xlon
      error stop "Initial latitude particle outside of domain."
    endif
    if ((part(i)%ylat.lt.0).or.(part(i)%ylat.gt.ny)) then
        write(*,*) 'Dimensions (nx,ny): ',nx,ny
      write(*,*) "Particle", i, "with ylat", part(i)%ylat
      error stop "Initial latitude particle outside of domain."
    endif
    if (part(i)%z.lt.0) then
      write(*,*) "Particle", i, "with height", part(i)%z
      error stop "Initial height particle below surface/sea level."
    endif
    do nsp=1,nspec
      if (mass(i,nsp).lt.0) then
        write(*,*) "Particle", i, "of species", nsp, "with mass", mass(i,nsp)
        error stop "Negative initial mass."
      endif
    end do
    if (part(i)%tstart*ldirect.lt.0) then
      if (lstart) then
        write(*,*) "WARNING: (some) particles have a starting time of", part(i)%tstart, "in"
        if (ldirect.le.0) then 
          write(*,*) "backward mode, please only use negative values."
          write(*,*) "time array in part_ic.nc should be given in seconds after the"
          write(*,*) "start of your simulation. Positive values will be converted to"
          write(*,*) "negative starting times..."
        else 
          write(*,*) "forward mode, please only use positive values."
          write(*,*) "time array in part_ic.nc should be given in seconds after the"
          write(*,*) "start of your simulation. Negative values will be converted to"
          write(*,*) "positive starting times..."
        endif
        lstart=.false.
      endif
      part(i)%tstart = part(i)%tstart*(-1)
    endif
  end do
  ! Release
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='release',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%npoint, & 
    start=(/ 1 /),count=(/ plen /)))
  ! ! Species
  ! call nf90_err(nf90_inq_varid(ncid=ncidend,name='species',varid=tempIDend))
  ! call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%species, & 
  !   start=(/ 1 /),count=(/ plen /)))

  ! Count number of releases
  numpoint=0
  ! Allocate plen of numpoint_max, since each particle could in principle have 
  ! a unique release number
  allocate(numpoint_max(plen), stat=stat)
  if (stat.ne.0) error stop "Could not allocate numpoint_max"
  numpoint_max=0
  release_max=0

  ! Count number of releases
  l1: do i=1,plen
    ! See if the release number already exists
    l2: do j=1,numpoint
      if (part(i)%npoint.eq.numpoint_max(numpoint)) then
        cycle l1
      endif
    end do l2
    numpoint = numpoint+1 ! Counting the number of releases
    numpoint_max(numpoint)=part(i)%npoint ! Save the release numbers
    ! Save maximum release number
    if (part(i)%npoint.gt.release_max) release_max=part(i)%npoint 
  end do l1

  if (numpoint.eq.0) numpoint=1

  allocate(kindz(numpoint),stat=stat)
  kindz=-1
  if (stat.ne.0) write(*,*)'ERROR: could not allocate kindz'
  ! Above sea-level or ground?
  call nf90_err(nf90_inquire_attribute(ncid=ncidend,name='kindz',varid=NF90_GLOBAL))
  call nf90_err(nf90_get_att(ncid=ncidend,varid=NF90_GLOBAL,name='kindz',values=zkind))
  
  kindz=zkind
  do j=1,numpoint
    if ((kindz(j).le.0).or.(kindz(j).ge.4)) then
      write(*,*) 'ERROR: kindz should be an integer between 1 and 3, not', kindz(nsp)
      error stop
    endif
    if (kindz(j).eq.3) then
      do i=1,plen
        if (part(i)%z.gt.1500.) then
          error stop 'Pressure heights should be given in hPa units. Input value exceeds surface pressure!'
        endif
      end do
    endif
  end do

  if (ioutputforeachrelease.eq.1) then
    maxpointspec_act=numpoint
  else
    maxpointspec_act=1
  endif

  if (release_max.gt.numpoint) then
    write(*,*) "WARNING: release numbers in part_ic.nc are not consecutive:", &
      release_max, "is larger than the total number of releases:", numpoint, &
      " Releases will be renumbered starting from 1."

    do j=1,numpoint
      do i=1,plen
        if (part(i)%npoint.eq.numpoint_max(j)) then
          part(i)%npoint=j
        endif
      end do
    end do
  endif
  deallocate(numpoint_max)

  ! Setting zpoint1 and zpoint2 necessary for wet backward deposition and plumes
  allocate(zpoint1(numpoint),zpoint2(numpoint), stat=stat)
  if (stat.ne.0) error stop "Could not allocate zpoint"
  zpoint2(:)=0.
  zpoint1(:)=1.e8
  do i=1,plen
    if (part(i)%npoint.ne.1) cycle ! This will be computed after information about
                                   ! topography (2) or pressure (3) is known (kindz_to_z)
    if (part(i)%z.gt.zpoint2(part(i)%npoint)) zpoint2(part(i)%npoint)=real(part(i)%z)
    if (part(i)%z.lt.zpoint1(part(i)%npoint)) zpoint1(part(i)%npoint)=real(part(i)%z)
  end do

  ! Setting xpoint1, ypoint1, xpoint2, ypoint2
  allocate(ypoint1(numpoint),ypoint2(numpoint), stat=stat)
  if (stat.ne.0) error stop "Could not allocate ypoint"
  allocate(xpoint1(numpoint),xpoint2(numpoint), stat=stat)
  if (stat.ne.0) error stop "Could not allocate xpoint"
  xpoint2(:)=0.
  xpoint1(:)=1.e8
  ypoint2(:)=0.
  ypoint1(:)=1.e8
  do i=1,plen
    if (part(i)%xlon.gt.xpoint2(part(i)%npoint)) xpoint2(part(i)%npoint)=real(part(i)%xlon)
    if (part(i)%xlon.lt.xpoint1(part(i)%npoint)) xpoint1(part(i)%npoint)=real(part(i)%xlon)
    if (part(i)%ylat.gt.ypoint2(part(i)%npoint)) ypoint2(part(i)%npoint)=real(part(i)%ylat)
    if (part(i)%ylat.lt.ypoint1(part(i)%npoint)) ypoint1(part(i)%npoint)=real(part(i)%ylat)
  end do
  do j=1,numpoint
    xpoint1(j)=(xpoint1(j)-xlon0)/dx
    xpoint2(j)=(xpoint2(j)-xlon0)/dx
    ypoint1(j)=(ypoint1(j)-ylat0)/dy
    ypoint2(j)=(ypoint2(j)-ylat0)/dy
  end do

  allocate(xmass(numpoint,nspec), npart(numpoint),ireleasestart(numpoint), &
    ireleaseend(numpoint), stat=stat)
  if (stat.ne.0) error stop "Could not allocate xmass,npart,ireleasestart,ireleaseend"
  xmass=0
  npart=0
  ireleasestart=-1
  ireleaseend=-1
  do i=1,plen
    do j=1,numpoint
      if (part(i)%npoint.eq.j) then
         do nsp=1,nspec
           xmass(j,nsp) = xmass(j,nsp)+mass(i,nsp)
         end do 
      endif
      if (part(i)%npoint.eq.j) then 
        npart(j)=npart(j)+1
        if ((ireleasestart(j).gt.part(i)%tstart*ldirect).or.(ireleasestart(j).eq.-1)) &
          ireleasestart(j)=part(i)%tstart
        if ((ireleaseend(j).le.part(i)%tstart*ldirect).or.(ireleaseend(j).eq.-1)) &
          ireleaseend(j)=part(i)%tstart
      endif
    end do
  end do

  part(:)%idt=mintime
  mass_init(:,:)=mass(:,:)
  do i=1,plen
    part(i)%nclass=min(int(ran1(idummy,0)*real(nclassunc))+1, &
         nclassunc)
  end do

  allocate(rho_rel(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate rho_rel in readinitconditions_netcdf'

  write(*,FMT='(A,ES14.7)') ' Total mass to be released:', sum(xmass(1:numpoint,1:nspec))
  call get_totalpart_num(numpart)
  numparticlecount=numpart
  call nf90_err(nf90_close(ncidend))

end subroutine readinitconditions_netcdf

end module netcdf_output_mod
