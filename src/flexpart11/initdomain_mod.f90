! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module initdomain_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for initializing         *
  !    particles and their mass for global domain-filling runs                 *
  !                                                                            *
  !*****************************************************************************

  use netcdf
  use par_mod
  use com_mod

  implicit none

  integer, allocatable, dimension(:)    :: specnini ! spec number in initial fields info
  integer, allocatable, dimension(:)    :: nxini    ! number grid cells longitude
  integer, allocatable, dimension(:)    :: nyini    ! number grid cells latitude
  integer, allocatable, dimension(:)    :: nzini    ! number grid cell vertical
  real, allocatable, dimension(:,:)     :: lonini   ! longitudes of initial fields 
  real, allocatable, dimension(:,:)     :: latini   ! latitudes of initial fields
  real, allocatable, dimension(:)       :: dxini, dyini ! longitude, latitude resolution of initial fields
  real, allocatable, dimension(:,:)     :: altini   ! altitudes of initial fields
  real, allocatable, dimension(:,:,:,:) :: prsini   ! pressure of initial fields
  real, allocatable, dimension(:,:,:,:) :: gridini  ! initial mixing ratios (dry air, ppbv)

  contains

  subroutine readgridini(ks, lexist)

  !******************************************************************************
  !                                                                             *
  !      Reads mixing ratios from a netcdf file with general format             *
  !                                                                             *
  !      Author: Rona Thompson, Sep-2023                                        *
  !                                                                             *
  !******************************************************************************
  !                                                                             *
  ! Variables:                                                                  *
  ! ks                 relative number of species                               *
  ! lexist             logical to indicate if INITCONC file specified           *
  !                                                                             *
  !******************************************************************************

    use date_mod
    use windfields_mod, only: nxmax,nymax,nzmax
    use netcdf_output_mod, only: nf90_err

    implicit none

    integer :: ks
    logical :: lexist
    integer :: ninit
    integer, dimension(:), allocatable :: specnum_rel
    character(len=256) :: path_name, file_name, unitinfo, strtmp1, strtmp2, nameout
    character(len=10)  :: var_name, hya_name, hyb_name, ps_name, prs_name, q_name, alt_name
    real :: coeff
    integer :: readerror
    integer :: ncid, dimid, varid, ndim, nvar, natt, unlimid, ndims, xtype, natts, len
    integer :: yyyymmdd, hhmiss, yyyy, mm
    integer, dimension(:), allocatable :: dimids
    character(len=256) :: dimname, attname
    character(len=2)  :: amonth
    character(len=4)  :: ayear
    integer :: nn, indxn, ntini, indxt, ix, jy, kz
    real, dimension(:), allocatable :: time, hya, hyb
    real, dimension(:,:), allocatable :: psurf
    real, dimension(:,:,:), allocatable :: shumid, pw
    real(kind=dp) :: julstart
    real(kind=dp), dimension(:), allocatable :: jdate
    real :: sclfact, offset
    logical :: lfexist
    integer,parameter :: unitinitconc=103, unitinitconcout=104

    ! declare namelists
    namelist /initconc_ctrl/ &
         ninit, &
         specnum_rel

    namelist /initconc/ &
         path_name, file_name, var_name, &
         hya_name, hyb_name, &
         ps_name, q_name, &
         prs_name, alt_name, &
         coeff

    ! Read initconc input info
    ! ************************

    allocate(specnum_rel(maxspec))

    ! presetting namelist initconc_ctrl
    ninit = -1  ! use negative value to determine failed namelist input
    specnum_rel(:) = 0

    ! read initconc_ctrl to find how many fields are specified 
    open(unitinitconc,file=path(1)(1:length(1))//'INITCONC',status='old',form='formatted',iostat=readerror)
    if (readerror.ne.0) then
      write(*,*) 'WARNING: cannot open file INITCONC'
      write(*,*) 'Trying to initialize using latitude profiles'
      lexist=.false.
      return
    endif

    ! check if namelist input provided
    read(unitinitconc,initconc_ctrl,iostat=readerror)
    if (readerror.ne.0 ) then
      write(*,*) 'ERROR in readgridini: cannot read INITCONC namelist'
      error stop
    endif
    write(*,*) 'ninit, specnum = ',ninit, specnum_rel

    ! namelist output
    if (nmlout.and.lroot) then
      inquire(file=path(2)(1:length(2))//'INITCONC.namelist',exist=lfexist)
      if (lfexist) then
        open(unitinitconcout,file=path(2)(1:length(2))//'INITCONC.namelist',status='old',&
              access='append',iostat=readerror)
      else
        open(unitinitconcout,file=path(2)(1:length(2))//'INITCONC.namelist',status='new',&
              iostat=readerror)
      endif
      if (readerror.ne.0) then
        write(*,*) 'ERROR in readgridini: file INITCONC cannot be opened'
        write(*,*) 'in the directory',path(2)(1:length(2))
        error stop
      endif
      if (.not.lfexist) then
        ! write this only once
        write(unitinitconcout,nml=initconc_ctrl)
      endif
    endif

    ! allocate variables used for initial mixing ratios
    if (.not.allocated(gridini)) then
      allocate( nxini(ninit), nyini(ninit), nzini(ninit) )
      allocate( dxini(ninit), dyini(ninit) )
      allocate( lonini(nxmax,ninit), latini(nymax,ninit), altini(nzmax,ninit) )
      allocate( prsini(nxmax,nymax,nzmax,ninit), gridini(nxmax,nymax,nzmax,ninit) )
      allocate( specnini(ninit) )
      nxini(:)=0
      nyini(:)=0
      nzini(:)=0
      lonini(:,:)=0.
      latini(:,:)=0.
      altini(:,:)=0.
      prsini(:,:,:,:)=0.
      gridini(:,:,:,:)=0.
      specnini(:)=specnum_rel(1:ninit)
    endif

    ! presetting namelist initconc
    path_name=""
    file_name=""
    var_name=""
    hya_name=""
    hyb_name=""
    ps_name=""
    q_name=""
    prs_name=""
    alt_name=""
    coeff=-1.

    ! read initconc file info
    do nn=1,ninit
      read(unitinitconc,initconc,iostat=readerror)
      if (readerror.ne.0) then
        write(*,*) 'ERROR in readgridini: cannot read file info for ',specnum_rel(nn)
        error stop
      endif
      if (specnum_rel(nn).eq.specnum(ks)) then
        indxn=nn ! index to gridini for this species
        exit
      endif
    end do

    ! namelist output
    if (nmlout.and.lroot) then
      write(unitinitconcout,nml=initconc)
    endif

    write(*,*) 'readgridini: path, filename, varname = ',path_name,file_name,var_name
    write(*,*) 'readgridini: specnum_rel(indxn), specnum(ks), indxn = ',specnum_rel(indxn), specnum(ks), indxn

    close(unitinitconc)
    close(unitinitconcout)

    ! get file name for current year and month
    call caldate(bdate,yyyymmdd,hhmiss)
    yyyy=yyyymmdd/10000
    mm=(yyyymmdd-(yyyymmdd/10000)*10000)/100
    write(ayear,'(I4)') yyyy
    write(amonth,'(I2.2)') mm
    nn=index(file_name,'YYYY',back=.false.)
    if (nn.ne.0) then
      strtmp1=file_name(1:nn-1)
      nn=index(file_name,'YYYY',back=.true.)
      strtmp2=file_name(nn+4:len_trim(file_name))
      file_name=trim(strtmp1)//ayear//trim(strtmp2)
      julstart=juldate((yyyymmdd/10000)*10000+101,0)
    endif
    nn=index(file_name,'MM',back=.false.)
    if (nn.ne.0) then
      strtmp1=file_name(1:nn-1)
      nn=index(file_name,'MM',back=.true.)
      strtmp2=file_name(nn+2:len_trim(file_name))
      file_name=trim(strtmp1)//amonth//trim(strtmp2)
      julstart=juldate((yyyymmdd/100)*100+1,0)
    endif
    write(*,*) 'readgridini: julstart = ',julstart
    write(*,*) 'readgridini: initial mixing ratio file to read: '//trim(path_name)//trim(file_name)

    ! check file exists
    inquire(file=trim(path_name)//trim(file_name), exist=lexist)
    if (.not.lexist) then
      write(*,*) 'ERROR readgridini: file not found: '//trim(path_name)//trim(file_name)
      error stop
    endif

    ! Read netcdf file
    !******************

    ! open file
    call nf90_err( nf90_open( trim(path_name)//trim(file_name), nf90_nowrite, ncid ) )

    ! inquire about dims and vars
    call nf90_err( nf90_inquire( ncid, ndim, nvar, natt, unlimid ) )
    write(*,*) 'ndim, nvar:', ndim, nvar
    allocate(dimids(ndim))

    ! read dimension info
    !********************
    do nn=1,ndim
      call nf90_err( nf90_inquire_dimension( ncid, nn, dimname, len ) )
      if ((index(dimname,'lon').ne.0).or.(index(dimname,'LON').ne.0) &
            .or.(index(dimname,'Lon').ne.0)) then
        if (len.gt.nxmax) then
          write(*,*) 'ERROR in readgridini: length longitude exceeds nxmax'
          error stop
        endif
        nxini(indxn)=len
        call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
        call nf90_err( nf90_get_var(ncid,varid,lonini(1:nxini(indxn),indxn)) )
        dxini(indxn)=abs(lonini(2,indxn)-lonini(1,indxn))
      else if ((index(dimname,'lat').ne.0).or.(index(dimname,'LAT').ne.0) &
            .or.(index(dimname,'Lat').ne.0)) then
        if (len.gt.nymax) then
          write(*,*) 'ERROR in readgridini: length latitude exceeds nymax'
          error stop
        endif
        nyini(indxn)=len
        call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
        call nf90_err( nf90_get_var(ncid,varid,latini(1:nyini(indxn),indxn)) )
        dyini(indxn)=abs(latini(2,indxn)-latini(1,indxn))
      else if (((index(dimname,'lev').ne.0).or.(index(dimname,'LEV').ne.0) &
            .or.(index(dimname,'Lev').ne.0)).and.(index(dimname,'hlevel').eq.0)) then
        if (len.gt.nzmax) then
          write(*,*) 'ERROR in readgridini: length vertical coord exceeds nzmax'
          error stop
        endif
        nzini(indxn)=len
      else if ((index(dimname,'time').ne.0).or.(index(dimname,'Time').ne.0) &
            .or.(index(dimname,'Date').ne.0).or.(index(dimname,'date').ne.0)) then
        ntini=len
        allocate( time(ntini), jdate(ntini) )
        call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
        call nf90_err( nf90_get_var(ncid,varid,time) )
        call nf90_err( nf90_get_att(ncid,varid,'units',unitinfo) )
        write(*,*) 'Time units: ',trim(unitinfo)
        if (index(unitinfo,'sec').ne.0) then
          jdate=real(time-time(1),kind=dp)/3600._dp/24._dp+julstart
          indxt=minloc(abs(jdate-bdate),dim=1)
        else if (index(unitinfo,'hour').ne.0) then
          jdate=real(time-time(1),kind=dp)/24._dp+julstart
          indxt=minloc(abs(jdate-bdate),dim=1)
        else if (index(unitinfo,'day').ne.0) then
          jdate=real(time-time(1),kind=dp)+julstart
          indxt=minloc(abs(jdate-bdate),dim=1)
        else if (index(unitinfo,'month').ne.0) then
          indxt=minloc(abs(time-mm),dim=1)
        else
          write(*,*) 'ERROR in readgridini: unknown time units in file: '//trim(path_name)//trim(file_name)
          error stop
        endif
        write(*,*) 'readgridini: time index in file = ',indxt
        deallocate( time, jdate )
      endif
    enddo

    if ((nxini(indxn).eq.0).or.(nyini(indxn).eq.0)) then
      write(*,*) 'ERROR in reagridini: unable to find lat and lon dimensions in file: '//trim(path_name)//trim(file_name)
      error stop
    endif
    write(*,*) 'readgridini: nxini(indxn), nyini(indxn), nzini(indxn) = ',nxini(indxn), nyini(indxn), nzini(indxn)
    write(*,*) 'readgridini: lonini(1:nxini(indxn),indxn) = ',lonini(1:nxini(indxn),indxn)
    write(*,*) 'readgridini: latini(1:nyini(indxn),indxn) = ',latini(1:nyini(indxn),indxn)

    ! read vertical coordinates
    !**************************
    if ((alt_name.eq."").and.(prs_name.eq."")) then
      ! hybrid pressure coordinates
      write(*,*) 'readgridini: reading hybrid pressure coordinates'
      if ((hya_name.eq."").or.(hyb_name.eq."").or.(ps_name.eq."")) then
        write(*,*) 'ERROR in readgridini: hybrid pressure coordinates and/or surface pressure missing'
        error stop
      endif
      ! read hybrid pressure coordinates
      allocate( hya(nzini(indxn)), hyb(nzini(indxn)), psurf(nxini(indxn),nyini(indxn)) )
      call nf90_err( nf90_inq_varid(ncid,trim(hya_name),varid) )
      call nf90_err( nf90_get_var(ncid,varid,hya) )
      call nf90_err( nf90_inq_varid(ncid,trim(hyb_name),varid) )
      call nf90_err( nf90_get_var(ncid,varid,hyb) )
      write(*,*) 'read hybrid pressure coords'
      ! read surface pressure
      call nf90_err( nf90_inq_varid(ncid,trim(ps_name),varid) )
      call nf90_err( nf90_inquire_variable(ncid,varid,nameout,xtype,ndims,dimids,natts) )
      write(*,*) 'psurf: ndims, natts = ',ndims, natts
      sclfact=0.
      offset=0.
      if (natts.gt.0 ) then
        do nn=1,natts
          call nf90_err( nf90_inq_attname(ncid,varid,nn,attname) )
          if (index(attname,'add_offset').ne.0) then
            call nf90_err( nf90_get_att(ncid,varid,'add_offset',offset) )
          else if (index(attname,'scale_factor').ne.0) then
            call nf90_err( nf90_get_att(ncid,varid,'scale_factor',sclfact) )
          endif
        end do
      endif
      write(*,*) 'psurf: sclfact, offset = ',sclfact, offset
      if (ndims.eq.2) then
        call nf90_err( nf90_get_var(ncid,varid,psurf,start=(/1,1/),count=(/nxini(indxn),nyini(indxn)/)) )
      else if (ndims.eq.3) then
        ! read from first time step (assume this is good enough)
        call nf90_err( nf90_get_var(ncid,varid,psurf,start=(/1,1,1/),count=(/nxini(indxn),nyini(indxn),1/)) )
      endif
      if ((sclfact.ne.0).and.(offset.ne.0)) then
        psurf=psurf*sclfact+offset
      endif
      ! calculate pressure
      do kz=1,nzini(indxn)
        prsini(1:nxini(indxn),1:nyini(indxn),kz,indxn)=hya(kz)+hyb(kz)*psurf(:,:)
      end do
      deallocate( psurf, hya, hyb )
    else if (alt_name.ne."") then
      ! height coordinates (assume metres above ground)
      call nf90_err( nf90_inq_varid(ncid,trim(alt_name),varid) )
      call nf90_err( nf90_get_var(ncid,varid,altini(1:nzini(indxn),indxn)) )
    else if (prs_name.ne."") then
      ! pressure coordinates
      call nf90_err( nf90_inq_varid(ncid,trim(prs_name),varid) )
      call nf90_err( nf90_get_var(ncid,varid,prsini(1,1,1:nzini(indxn),indxn)) )
      do jy=1,nyini(indxn)
        do ix=1,nxini(indxn)
          prsini(ix,jy,:,indxn)=prsini(1,1,:,indxn)
        end do
      end do
    endif

    ! read mixing ratio variables
    !****************************
    if (var_name.eq."") then
      write(*,*) 'ERROR in readgridini: mixing ratios missing in file'//trim(path_name)//trim(file_name)
      error stop
    endif
    call nf90_err( nf90_inq_varid(ncid,trim(var_name),varid) )
    call nf90_err( nf90_inquire_variable(ncid,varid,nameout,xtype,ndims,dimids,natts) )
    write(*,*) 'conc: ndims, natts, dimids = ',ndims, natts, dimids
    sclfact=0.
    offset=0.
    if (natts.gt.0 ) then
      do nn=1,natts
        call nf90_err( nf90_inq_attname(ncid,varid,nn,attname) )
        if (index(attname,'add_offset').ne.0) then
          call nf90_err( nf90_get_att(ncid,varid,'add_offset',offset) )
        else if (index(attname,'scale_factor').ne.0) then
          call nf90_err( nf90_get_att(ncid,varid,'scale_factor',sclfact) )
        endif
      end do
    endif
    write(*,*) 'conc: sclfact, offset = ',sclfact, offset
    if (ndims.eq.3) then
      call nf90_err( nf90_get_var(ncid,varid,gridini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn), &
                         start=(/1,1,1/),count=(/nxini(indxn),nyini(indxn),nzini(indxn)/)) )
    else if (ndims.eq.4) then
      call nf90_err( nf90_get_var(ncid,varid,gridini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn), &
                         start=(/1,1,1,indxt/),count=(/nxini(indxn),nyini(indxn),nzini(indxn),1/)) )
    endif
    if ((sclfact.ne.0).and.(offset.ne.0)) then
      gridini(:,:,:,indxn)=gridini(:,:,:,indxn)*sclfact+offset
    endif

    ! other variables
    !****************
    if (q_name.ne."") then
      ! read specific humidity
      allocate( shumid(nxini(indxn),nyini(indxn),nzini(indxn)), pw(nxini(indxn),nyini(indxn),nzini(indxn)) )
      call nf90_err( nf90_inq_varid(ncid,trim(q_name),varid) )
      call nf90_err( nf90_inquire_variable(ncid,varid,nameout,xtype,ndims,dimids,natts) )
      sclfact=0.
      offset=0.
      if (natts.gt.0 ) then
        do nn=1,natts
          call nf90_err( nf90_inq_attname(ncid,varid,nn,attname) )
          if (index(attname,'add_offset').ne.0) then
            call nf90_err( nf90_get_att(ncid,varid,'add_offset',offset) )
          else if (index(attname,'scale_factor').ne.0) then
            call nf90_err( nf90_get_att(ncid,varid,'scale_factor',sclfact) )
          endif
        end do
      endif
      if (ndims.eq.3) then
        call nf90_err( nf90_get_var(ncid,varid,shumid(:,:,:), &
                           start=(/1,1,1/),count=(/nxini(indxn),nyini(indxn),nzini(indxn)/)) )
      else if (ndims.eq.4) then
        call nf90_err( nf90_get_var(ncid,varid,shumid(:,:,:), &
                           start=(/1,1,1,indxt/),count=(/nxini(indxn),nyini(indxn),nzini(indxn),1/)) )
      endif
      write(*,*) 'conc: sclfact, offset = ',sclfact, offset
      if ((sclfact.ne.0).and.(offset.ne.0)) then
        shumid=shumid*sclfact+offset
      endif
      ! partial pressure of water from specific humidity
      pw = shumid*prsini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn)/(0.622 + 0.378*shumid)
      ! correct mixing ratio to dry air (needed in init_domainfill)
      gridini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn)= &
                   gridini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn)* &
                   prsini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn)/ &
                   (prsini(1:nxini(indxn),1:nyini(indxn),1:nzini(indxn),indxn) - pw)
      deallocate( shumid, pw )
    endif

    ! convert units to ppbv
    !**********************
    if (coeff.gt.0) gridini(:,:,:,indxn)=gridini(:,:,:,indxn)*coeff
    write(*,*) 'range(gridini) = ',minval(gridini(1:nxini(indxn),1:nyini(indxn),1,indxn)),&
                    maxval(gridini(1:nxini(indxn),1:nyini(indxn),1,indxn))

  end subroutine readgridini


  subroutine init_domainfill_ncf

  !******************************************************************************
  !                                                                             *
  !  Initializes particles equally distributed over the first release location  *
  !  specified in file RELEASES. This box is assumed to be the domain for doing *
  !  domain-filling trajectory calculations.                                    *
  !  All particles carry the same amount of mass which alltogether comprises the*
  !  mass of air within the box.                                                *
  !                                                                             *
  !  Author: A. Stohl, Oct 2002                                                 *
  !  Modifications:                                                             * 
  !    R. Thompson, Sep 2023: added initialization of mass from grid based      *
  !                           on code of S. Henne for Flexpart-CTM              *
  !                                                                             *
  !******************************************************************************
  !                                                                             *
  ! Variables:                                                                  *
  !                                                                             *
  ! numparticlecount    consecutively counts the number of particles released   *
  ! nx_we(2)       grid indices for western and eastern boundary of domain-     *
  !                filling trajectory calculations                              *
  ! ny_sn(2)       grid indices for southern and northern boundary of domain-   *
  !                filling trajectory calculations                              *
  !                                                                             *
  !******************************************************************************

    use point_mod
    use random_mod
    use outgrid_mod
    use particle_mod
#ifdef ETA
    use coord_ecmwf_mod
#endif
    use initialise_mod, only: nx_we,ny_sn,numcolumn,numcolumn_we,numcolumn_sn, &
                              zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn, &
                              xmassperparticle,alloc_domainfill
    use totals_mod,     only: tot_mass

    implicit none

    integer :: j,ix,jy,kz,ncolumn,numparttot,iterminate,stat
    real :: ylat,ylatp,ylatm,hzone
    real :: cosfactm,cosfactp,deltacol,dz1,dz2,dz,pnew,pnew_temp,fractus
    real,parameter :: pih=pi/180.
    real :: colmasstotal,zposition
    real :: hgt_tmp
    real, allocatable,dimension(:) :: pp

    integer :: ixm,ixp,jym,jyp,indzm,indzp,nn,indzh,i,ii,jj,ks,indxn
    real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2)
    integer :: idummy = -11
    logical :: deall
    real,parameter :: eps=1.e-6

    ! dry and moist air density at particle position
    real :: rho_d_i, rho_m_i

    ! variables for column mass calculation 
    integer,allocatable,dimension(:,:) :: nncolumn
    real,allocatable,dimension(:)   :: gridarea 
    real,allocatable,dimension(:,:) :: colmass 
    real,parameter :: weightair=28.97

    ! variables to store source profiles
    real :: bg_lat(maxspec,ny), dummy1, dummy2
    real :: ppbvpart, presspart
    character(256) :: filename
    logical :: lexist, lgridini
    real :: xl, yl

    ! io variables
    character(30) :: frmt
    integer, parameter :: unitcolmass=98, unitinitconc=99

    ! Determine the release region (only full grid cells), over which particles
    ! shall be initialized
    ! Use 2 fields for west/east and south/north boundary
    !**************************************************************************
    call alloc_domainfill
    nx_we(1)=max(int(xpoint1(1)),0)
    nx_we(2)=min((int(xpoint2(1))+1),nxmin1)
    ny_sn(1)=max(int(ypoint1(1)),0)
    ny_sn(2)=min((int(ypoint2(1))+1),nymin1)
 
    ! For global simulations (both global wind data and global domain-filling),
    ! set a switch, such that no boundary conditions are used
    !**************************************************************************
    if (xglobal.and.sglobal.and.nglobal) then
      if ((nx_we(1).eq.0).and.(nx_we(2).eq.nxmin1).and. &
           (ny_sn(1).eq.0).and.(ny_sn(2).eq.nymin1)) then
        gdomainfill=.true.
      else
        gdomainfill=.false.
      endif
    endif

    ! If resuming a run from particle dump
    ! calculate total mass each species then exit
    !********************************************
    if (gdomainfill.and.ipin.ne.0) then
      write(*,*) 'Initialising particles from partoutput'
      tot_mass(:)=0.
      do ks=2,nspec
        tot_mass(ks)=sum(mass(1:count%alive,ks))
        write(*,'(A,E12.4,A)') 'Species '//species(ks)//': ',tot_mass(ks),' (kg)'
      end do
      return
    endif

    ! Allocate fields used within this subroutine
    !*********************************************
    allocate( nncolumn(0:nxmax-1, 0:nymax-1),stat=stat )
    if (stat.ne.0) write(*,*)'ERROR: could not allocate nncolumn'
    allocate(gridarea(0:nymax-1),colmass(0:nxmax-2,0:nymax-1),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridarea or colmass'

    ! Do not release particles twice (i.e., not at both in the leftmost and rightmost
    ! grid cell) for a global domain
    !*****************************************************************************
    if (xglobal) nx_we(2)=min(nx_we(2),nx-2)
    write(*,*) 'init_domainfill: nx_we, ny_sn = ',nx_we, ny_sn
    
    ! Try reading initial mixing ratio fields into gridini
    !*****************************************************
    ! this is to restart from 3D fields
    ! grid_pptv files have to be present for all species except 1 (airtracer)
    lgridini = .true.
    do ks=2,nspec ! species 1 is assumed to be air tracer
      call readgridini(ks,lexist)
      if (.not.lexist) then
        lgridini = .false.
        exit
      end if
    end do

    if (lgridini) then
      write(*,*) "Initialising all species with gridded input"
      !! test
      open(unitoutgrid,file=path(2)(1:length(2))//'gridini.txt',status='replace',action='write')
      write(frmt,fmt='(A,I4,A)') '(',nxini(1),'(E14.6))'
      do kz=1,nzini(1)
        do jj=1,nyini(1)
          write(unitoutgrid,frmt) gridini(1:nxini(1),jj,kz,1)
        end do
      end do
      close(unitoutgrid)
      !!
    else
      ! read latitude profiles for initialization
      do ks=2,nspec  ! species 1 is assumed to be air tracer
        filename=trim(path(2)(1:length(2)))//'latitude_profile_'//trim(species(ks))//'.txt'
        print*, 'init_domainfill: filename = ',filename
        inquire(file=filename,exist=lexist)
        if (lexist) then
          open(unitinitconc,file=filename,action='read')
          do jj=1,ny
            read(unitinitconc,'(3F13.6)') dummy1, bg_lat(ks,jj), dummy2
            write(*,*) 'read values; ', dummy1, bg_lat(ks,jj), dummy2
          enddo
          close(unitinitconc)
          write(*,'(A)') "Initialising species "//species(ks)//" with latitudinal profile"
        else
          bg_lat(ks,:) = 1.
          write(*,'(A)') "Initialising species "//species(ks)//" with value in RELEASES"
        endif
      end do ! nspec
    endif ! lgridini

    ! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
    ! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
    ! Note gridarea is for meteo grid
    !************************************************************

    write(*,*) 'init_domainfill: nxmin1, nxmax, nymin1, nymax = ',nxmin1, nxmax, nymin1, nymax

    ! for the south pole
    if (sglobal) then
      ylat=ylat0
      ylatp=ylat+0.5*dy
      ylatm=ylat
      cosfactm=0.
      cosfactp=cos(ylatp*pih)*r_earth
      hzone=sqrt(r_earth**2-cosfactm**2)- &
           sqrt(r_earth**2-cosfactp**2)
      gridarea(0)=2.*pi*r_earth*hzone*dx/360.
    endif

    ! do the same for the north pole
    if (nglobal) then
      ylat=ylat0+real(nymin1)*dy
      ylatp=ylat
      ylatm=ylat-0.5*dy
      cosfactp=0.
      cosfactm=cos(ylatm*pih)*r_earth
      hzone=sqrt(r_earth**2-cosfactp**2)- &
           sqrt(r_earth**2-cosfactm**2)
      gridarea(nymin1)=2.*pi*r_earth*hzone*dx/360.
    endif

    ! Initialise the sum over the total mass of the atmosphere
    colmasstotal=0.
    colmass(:,:)=0.

    write(*,*) 'init_domainfill: ny_sn, nx_we = ',ny_sn, nx_we
    write(*,*) 'init_domainfill: dxini, dyini = ',dxini, dyini

    allocate( pp(nzmax),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate pp'

!$OMP PARALLEL PRIVATE(jy,ix,ylat,ylatp,ylatm,hzone,cosfactp,cosfactm,pp) 

!$OMP DO
    ! loop over latitudes
    do jy=ny_sn(1),ny_sn(2)
      if (sglobal.and.(jy.eq.ny_sn(1))) cycle
      if (nglobal.and.(jy.eq.ny_sn(2))) cycle
      ylat=ylat0+real(jy)*dy
      ylatp=ylat+0.5*dy
      ylatm=ylat-0.5*dy
      if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
        hzone=1./dyconst
      else
        cosfactp=cos(ylatp*pih)*r_earth
        cosfactm=cos(ylatm*pih)*r_earth
        if (cosfactp.lt.cosfactm) then
          hzone=sqrt(r_earth**2-cosfactp**2)- &
               sqrt(r_earth**2-cosfactm**2)
        else
          hzone=sqrt(r_earth**2-cosfactm**2)- &
               sqrt(r_earth**2-cosfactp**2)
        endif
      endif
      gridarea(jy)=2.*pi*r_earth*hzone*dx/360.
    end do
!$OMP END DO
!$OMP BARRIER

    ! Calculate total mass of each grid column and of the whole atmosphere
    !*********************************************************************

!$OMP DO
    do jy=ny_sn(1),ny_sn(2)        
      do ix=nx_we(1),nx_we(2)      
        pp(1)=prs(ix,jy,1,1)
        pp(nz)=prs(ix,jy,nz,1)
        colmass(ix,jy)=(pp(1)-pp(nz))/ga*gridarea(jy)
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    deallocate(pp)

    colmasstotal=sum(colmass)
    write(*,*) 'Atmospheric mass air = ',colmasstotal

    ! Output of colmass distribution
    !********************************

    open(unitcolmass,file=path(2)(1:length(2))//'colmass.dat',action='write')
    write(frmt, '(A, I4, A)') '(', ny_sn(2)-ny_sn(1)+1, 'E12.3)'
    do ix=nx_we(1),nx_we(2)
      write(unitcolmass, frmt) (colmass(ix,i),i=0,(nymax-1))
    end do
    close(unitcolmass)

    ! If not continuing from particle dump
    !*************************************

    if (ipin.eq.0) numpart=0

    ! Determine the particle positions
    !*********************************
    allocate( pp(nzmax),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate pp'
    numparttot=0
    numcolumn=0
    iterminate=0

    ! allocate all particles before loop
    call spawn_particles(0, npart(1))
    write(*,*)  'init_domainfill: count%alive = ',count%alive

    do jy=ny_sn(1),ny_sn(2)        ! loop about latitudes
      ylat=ylat0+real(jy)*dy
      do ix=nx_we(1),nx_we(2)      ! loop about longitudes
        ncolumn=nint(0.999*real(npart(1))*colmass(ix,jy)/colmasstotal)
        ! this condition means with 0.5 degrees need around 200 million particles 
        ! to avoid any grid cells having zero particles
        ncolumn=max(ncolumn,1)
        nncolumn(ix,jy) = ncolumn
        if (ncolumn.gt.numcolumn) numcolumn=ncolumn

        ! Calculate pressure at the altitudes of model surfaces, using the air density
        ! information, which is stored as a 3-d field
        pp(:)=prs(ix,jy,:,1)

        ! Loop over number of particles in grid column
        deltacol=(pp(1)-pp(nz))/real(ncolumn)
        pnew=pp(1)+deltacol/2.
        jj=0
        do j=1,ncolumn      
          jj=jj+1

          ! For columns with many particles (i.e. around the equator), distribute
          ! the particles equally, for columns with few particles (i.e. around the
          ! poles), distribute the particles randomly. When only few particles are
          ! left distribute them randomly

          if ((ncolumn.gt.20).and.(ncolumn-j.gt.20)) then
            pnew_temp=pnew-ran1(idummy,0)*deltacol
            pnew=pnew-deltacol
          else if (ncolumn.gt.20) then
            pnew_temp=pnew-ran1(idummy,0)*(pnew-pp(nz))
          else
            pnew_temp=pp(1)-ran1(idummy,0)*(pp(1)-pp(nz))
          endif
          pnew_temp=min(pp(1),pnew_temp) 
          pnew_temp=max(pp(nz)+eps,pnew_temp)

          ! find vertical layer
          do kz=1,nz-1
            if ((pp(kz).ge.pnew_temp).and.(pp(kz+1).lt.pnew_temp)) then
              dz1=log(pp(kz))-log(pnew_temp)
              dz2=log(pnew_temp)-log(pp(kz+1))
              dz=1./(dz1+dz2)

              ! Assign particle position
              !*************************

              ! Do the following steps only if particles are not read in from 
              ! previous model run
              if (ipin.eq.0) then
                
                ! Horizontal position
                call set_xlon(numpart+jj,real(real(ix)-0.5+ran1(idummy,0),kind=dp))
                if (ix.eq.0) call set_xlon(numpart+jj,real(ran1(idummy,0),kind=dp))
                if (ix.eq.nxmin1) &
                     call set_xlon(numpart+jj,real(real(nxmin1)-ran1(idummy,0),kind=dp))
                call set_ylat(numpart+jj,real(real(jy)-0.5+ran1(idummy,0),kind=dp))
                if (jy.eq.0) call set_ylat(numpart+jj,real(ran1(idummy,0),kind=dp))
                if (jy.eq.nymin1) &
                     call set_ylat(numpart+jj,real(real(nymin1)-ran1(idummy,0),kind=dp))  
                ! Vertical position
                hgt_tmp=(height(kz)*dz2+height(kz+1)*dz1)*dz
                call set_z(numpart+jj,hgt_tmp)
                if (real(part(numpart+jj)%z).gt.(height(nz)-0.5)) &
                     call set_z(numpart+jj,height(nz)-0.5)
#ifdef ETA
                call update_z_to_zeta(0, numpart+jj)
#endif
                ! Interpolate PV to the particle position
                !****************************************
                ixm=int(part(numpart+jj)%xlon)
                jym=int(part(numpart+jj)%ylat)
                ixp=ixm+1
                jyp=jym+1
                ddx=real(part(numpart+jj)%xlon)-real(ixm)
                ddy=real(part(numpart+jj)%ylat)-real(jym)
                rddx=1.-ddx
                rddy=1.-ddy
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy
                indzm=nz-1
                indzp=nz
                do ii=2,nz
                  if (real(height(ii),kind=dp).gt.part(numpart+jj)%z) then
                    indzm=ii-1
                    indzp=ii
                    exit
                  endif
                end do
                dz1=real(part(numpart+jj)%z)-height(indzm)
                dz2=height(indzp)-real(part(numpart+jj)%z)
                dz=1./(dz1+dz2)
                do ii=1,2
                  indzh=indzm+ii-1
                  y1(ii)=p1*pv(ixm,jym,indzh,1) &
                       + p2*pv(ixp,jym,indzh,1) &
                       + p3*pv(ixm,jyp,indzh,1) &
                       + p4*pv(ixp,jyp,indzh,1)
                end do
                pvpart=(dz2*y1(1)+dz1*y1(2))*dz
                if (ylat.lt.0.) pvpart=-1.*pvpart

                ! Interpolate moist air density to the particle position
                !*******************************************************
                do ii=1,2
                  indzh=indzm+ii-1
                  y1(ii)=p1*rho(ixm,jym,indzh,1) &
                       +p2*rho(ixp,jym,indzh,1) &
                       +p3*rho(ixm,jyp,indzh,1) &
                       +p4*rho(ixp,jyp,indzh,1)
                end do
                rho_m_i=(dz2*y1(1)+dz1*y1(2))*dz

                ! Interpolate dry air density to the particle position
                !*****************************************************
                do ii=1,2
                  indzh=indzm+ii-1
                  y1(ii)=p1*rho(ixm,jym,indzh,1)*(1-qv(ixm,jym,indzh,1)) &
                       +p2*rho(ixp,jym,indzh,1)*(1-qv(ixp,jym,indzh,1)) &
                       +p3*rho(ixm,jyp,indzh,1)*(1-qv(ixm,jyp,indzh,1)) &
                       +p4*rho(ixp,jyp,indzh,1)*(1-qv(ixp,jyp,indzh,1))
                end do
                rho_d_i=(dz2*y1(1)+dz1*y1(2))*dz

                ! For domain-filling option 2 (stratospheric O3), 
                ! do the rest only in the stratosphere
                !************************************************

                if (((part(numpart+jj)%z.gt.3000.).and. &
                     (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then

                  ! Assign certain properties to the particle
                  !******************************************
                  part(numpart+jj)%nclass=min(int(ran1(idummy,0)*real(nclassunc))+1,nclassunc)
                  numparticlecount=numparticlecount+1
                  part(numpart+jj)%npoint=numparticlecount
                  part(numpart+jj)%idt=mintime
                  mass(numpart+jj,1)=colmass(ix,jy)/real(ncolumn)

                  if (lgridini) then

                    ! Initialize particle mass using gridded input
                    !*********************************************
                    ! Assume input is in units of ppv and species 1 carries airmass tracer so
                    ! mass of other species is easily determined
                    ! loop over all species, assuming species 1 is airtracer
                    do ks=2, nspec
                      indxn=minloc(abs(specnini-specnum(ks)),dim=1)
                      ! lon and lat of particle
                      xl=real(part(numpart+jj)%xlon)*dx+xlon0
                      yl=real(part(numpart+jj)%ylat)*dy+ylat0
                      ! get coordinates in gridini
                      ! Assumes lon and lat dimensions are midpoints
                      ixm=int((xl-(lonini(1,indxn)-0.5*dxini(indxn)))/dxini(indxn))+1
                      jym=int((yl-(latini(1,indxn)-0.5*dyini(indxn)))/dyini(indxn))+1
                      ixm=min(ixm,nxini(indxn))
                      jym=min(jym,nyini(indxn))
                      !! testing
!                      if (jj.eq.1.and.jy.lt.5.and.ix.lt.5) then 
!                        print*, 'init_domainfill: lonini, xl, latini, yl = ',lonini(ixm,indxn),xl,latini(jym,indxn),yl
!                      endif !!
                      ! Get vertical position in gridini
                      if (any(altini(:,indxn).gt.0)) then
                        ! vertical coordinate in metres above ground
                        indzm=nzini(indxn)-1
                        indzp=nzini(indxn)
                        do ii=2,nzini(indxn)
                          if (altini(ii,indxn).gt.real(part(numpart+jj)%z)) then
                            indzm=ii-1
                            indzp=ii
                            exit
                          endif
                        enddo
                        dz1=real(part(numpart+jj)%z)-altini(indzm,indxn)
                        dz2=altini(indzp,indxn)-real(part(numpart+jj)%z)
                        dz=1./(dz1+dz2)
                        ppbvpart=(dz2*gridini(ixm,jym,indzm,indxn)+ &
                                    dz1*gridini(ixm,jym,indzp,indxn))*dz
                      else if (any(prsini(:,:,:,indxn).gt.0)) then
                        ! vertical coordinate in pressure (Pa)
                        presspart=pnew_temp
                        indzm=nzini(indxn)-1
                        indzp=nzini(indxn)
                        do ii=2,nzini(indxn)
                          if (presspart.gt.prsini(ixm,jym,ii,indxn)) then
                            indzm=ii-1
                            indzp=ii
                            exit
                          endif
                        end do
                        dz1=presspart-prsini(ixm,jym,indzm,indxn)
                        dz2=prsini(ixm,jym,indzp,indxn)-presspart
                        dz=1./(dz1+dz2)
                        ppbvpart=(dz2*gridini(ixm,jym,indzm,indxn)+ &
                                  dz1*gridini(ixm,jym,indzp,indxn))*dz
                      endif 
                      !! test
!                      if (numpart.lt.100) write(*,*) 'init_domainfill: ratio dry/moist density =',rho_d_i/rho_m_i
                      mass(numpart+jj,ks)=mass(numpart+jj,1) * &
                           weightmolar(ks)/weightair * &
                           rho_d_i/rho_m_i*ppbvpart/1.E9
                    end do ! nspec

                  else

                    ! Initialize with latitude profile
                    !*********************************
                    ! loop over all species, assuming species 1 is airtracer
                    do ks=2, nspec
                      mass(numpart+jj,ks)=mass(numpart+jj,1)* &
                           weightmolar(ks)/weightair * &
                           rho_d_i/rho_m_i*bg_lat(ks,jy)/1.E9
                    end do

                  endif ! lgridini

                  ! Assign ozone mass if domain-filling option 2
                  !*********************************************
                  if (mdomainfill.eq.2) mass(numpart+jj,1)= &
                       mass(numpart+jj,1)*pvpart*48./29.*ozonescale/10.**9

                  mass_init(numpart+jj,1)=mass(numpart+jj,1)

                else

                  ! Particle in stratosphere and not domain-filling option 1
                  !*********************************************************
                  call terminate_particle(numpart+jj, 0)
                  jj=jj-1
                  iterminate=iterminate+1

                endif ! domainfill option
              endif ! if initialization
            endif ! if in layer
          end do ! loop over layers
        end do ! loop over column

        numparttot=numparttot+ncolumn
        if (ipin.eq.0) numpart=numpart+jj

      end do ! loop over longitude

    end do ! loop over latitude

    write(*,*) 'init_domainfill: numpart, numparttot = ',numpart, numparttot

    ! Terminate unused particles
    !***************************
    do j=(numpart+1),count%alive
      call terminate_particle(j,0) ! Cannot be within an OMP region
      iterminate=iterminate+1
    end do  
    write(*,*) 'init_domainfill: after terminating extra particles count%alive = ',count%alive

    ! Total mass each species
    !************************
    tot_mass(:)=0.
    do ks=2,nspec
      tot_mass(ks)=sum(mass(1:count%alive,ks))
      write(*,'(A,E12.4,A)') 'Species '//species(ks)//': ',tot_mass(ks),' (kg)'
    end do

    xmassperparticle=colmasstotal/real(numparttot)

    ! Output colmass distribution
    !****************************

    open(unitcolmass,file=path(2)(1:length(2))//'ncolumn.dat',action='write')
    write(frmt, '(A, I4, A)') '(', ny_sn(2)-ny_sn(1)+1, 'I5)'
    do ix=nx_we(1), nx_we(2)
      write(unitcolmass,frmt) (nncolumn(ix, jj), jj=ny_sn(1),ny_sn(2))
    end do
    close(unitcolmass)

    ! For boundary conditions, we need fewer particle release heights per column,
    ! because otherwise it takes too long until enough mass has accumulated to
    ! release a particle at the boundary (would take dx/u seconds), leading to
    ! relatively large position errors of the order of one grid distance.
    ! It's better to release fewer particles per column, but to do so more often.
    ! Thus, use on the order of nz starting heights per column.
    ! We thus repeat the above to determine fewer starting heights, that are
    ! used furtheron in subroutine boundcond_domainfill.f90
    !****************************************************************************

    fractus=real(numcolumn)/real(nz)
    write(*,*) 'Total number of particles at model start: ',numpart
    write(*,*) 'Maximum number of particles per column: ',numcolumn
    write(*,*) 'If ',fractus,' <1, better use more particles'
    fractus=sqrt(max(fractus,1.))/2.

    do jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
      do ix=nx_we(1),nx_we(2)      ! loop about longitudes
        ncolumn=nint(0.999/fractus*real(npart(1))*colmass(ix,jy) &
           /colmasstotal)
        if (ncolumn.gt.maxcolumn) error stop 'maxcolumn too small'
        if (ncolumn.eq.0) cycle

        ! Memorize how many particles per column shall be used for all boundaries
        ! This is further used in subroutine boundcond_domainfill.f
        ! Use 2 fields for west/east and south/north boundary
        !************************************************************************

        if (ix.eq.nx_we(1)) numcolumn_we(1,jy)=ncolumn
        if (ix.eq.nx_we(2)) numcolumn_we(2,jy)=ncolumn
        if (jy.eq.ny_sn(1)) numcolumn_sn(1,ix)=ncolumn
        if (jy.eq.ny_sn(2)) numcolumn_sn(2,ix)=ncolumn

        ! Calculate pressure at the altitudes of model surfaces, using the air density
        ! information, which is stored as a 3-d field
        !*****************************************************************************

        do kz=1,nz
          pp(kz)=prs(ix,jy,kz,1) 
        end do

        ! Determine the reference starting altitudes
        !*******************************************

        deltacol=(pp(1)-pp(nz))/real(ncolumn)
        pnew=pp(1)+deltacol/2.
        do j=1,ncolumn
          pnew=pnew-deltacol
          do kz=1,nz-1
            if ((pp(kz).ge.pnew).and.(pp(kz+1).lt.pnew)) then
              dz1=pp(kz)-pnew
              dz2=pnew-pp(kz+1)
              dz=1./(dz1+dz2)
              zposition=(height(kz)*dz2+height(kz+1)*dz1)*dz
              if (zposition.gt.height(nz)-0.5) zposition=height(nz)-0.5

        ! Memorize vertical positions where particles are introduced
        ! This is further used in subroutine boundcond_domainfill.f
        !***********************************************************

              if (ix.eq.nx_we(1)) zcolumn_we(1,jy,j)=zposition
              if (ix.eq.nx_we(2)) zcolumn_we(2,jy,j)=zposition
              if (jy.eq.ny_sn(1)) zcolumn_sn(1,ix,j)=zposition
              if (jy.eq.ny_sn(2)) zcolumn_sn(2,ix,j)=zposition

        ! Initialize mass that has accumulated at boundary to zero
        !*********************************************************

              acc_mass_we(1,jy,j)=0.
              acc_mass_we(2,jy,j)=0.
              acc_mass_sn(1,jy,j)=0.
              acc_mass_sn(2,jy,j)=0.

            endif
          end do
        end do
      end do
    end do

    ! If there were more particles allocated than used,
    ! Deallocate unused memory and update numpart
    !**************************************************

    deall=.false.
    do i=numpart, 1, -1
      if (.not. part(i)%alive) then
        deall=.true.
        numpart = numpart - 1
      else
        exit
      endif
    end do

    if (deall) call dealloc_particle(numpart) ! deallocates everything above numpart
    write(*,*) 'init_domainfill: after dealloc count%alive = ',count%alive
    write(*,*) 'init_domainfill: count%allocated = ',count%allocated


    ! If particles shall be read in to continue an existing run,
    ! then the accumulated masses at the domain boundaries must be read in, too.
    ! This overrides any previous calculations.
    !***************************************************************************

    if ((ipin.eq.1).and.(.not.gdomainfill)) then
      open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
           form='unformatted')
      read(unitboundcond) numcolumn_we,numcolumn_sn, &
           zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
      close(unitboundcond)
    endif

    deallocate(pp,nncolumn,gridarea)

  end subroutine init_domainfill_ncf

  !*****************************************************************************
  !                                                                            *
  !    netcdf error message handling                                           *
  !                                                                            *
  !*****************************************************************************

!  subroutine nf90_err(status)
!
!    integer, intent (in) :: status
!
!    if(status /= nf90_noerr) then
!      print*, trim(nf90_strerror(status))
!      error stop 'Stopped'
!    end if
!
!  end subroutine nf90_err

end module initdomain_mod
