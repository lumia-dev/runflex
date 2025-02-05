! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for  calculating         *
  !    chemical loss of species                                                *
  !                                                                            *
  !*****************************************************************************

module chemistry_mod

  use netcdf
  use par_mod
  use com_mod
  use date_mod
  use windfields_mod,    only: rho,nxmax,nymax,nzmax
  use totals_mod,        only: chem_loss
  use netcdf_output_mod, only: nf90_err

  ! reagent field variables 

  implicit none

  integer, allocatable, dimension(:)      :: nxCL, nyCL, nzCL
  integer                                 :: nxjr, nyjr
  real, allocatable, dimension(:,:)       :: lonCL, latCL, altCL
  real, allocatable, dimension(:)         :: dxCL, dyCL
  real                                    :: dxjr, dyjr
  real, allocatable, dimension(:,:,:,:,:) :: CL_field           ! chemical loss fields at input resolution
  real, allocatable, dimension(:,:,:,:)   :: clfield_cur        ! chemical loss fields current hour
  integer, dimension(2)                   :: memCLtime          ! time of fields in memory (sec)  
  integer                                 :: curCLhour          ! current hour since start of simulation
  real(kind=dp)                           :: memCLday           ! current day of simulation
  real(kind=dp), dimension(2)             :: CL_time            ! julian date of fields in memory
  real, allocatable, dimension(:,:)       :: jrate_average      ! daily average actinic flux 
  real, allocatable, dimension(:)         :: lonjr, latjr

  private :: zenithangle, photo_O1D

  contains

  subroutine readreagents()

  !*****************************************************************************
  !                                                                            *
  !   This routine reads names and input paths of chemical reagents            *
  !                                                                            *
  !   Author: R. Thompson, Sep 2023                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! preagent              Reagent name                                         *
  ! preag_path            Path to reagent fields                               *
  ! phourly               Logical for hourly interpolate (0 = no, 1 = yes)     *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: preag_path
    character(len=16)  :: preagent, punit
    integer :: phourly
    integer,parameter :: unitreagents=2, unitreagentsout=3
    integer :: readerror
    integer :: j

    ! declare namelist
    namelist /reagent_params/ &
         preagent, preag_path, phourly, punit

    preagent="" ! read failure indicator value
    preag_path=""
    punit=""
    phourly=0
    reag_hourly(:)=0

    open(unitreagents,file=path(1)(1:length(1))//'REAGENTS',status='old',form='formatted',iostat=readerror)
    if (readerror.ne.0) then
      ! no REAGENTS file
      nreagent=0
      go to 999
    endif

    ! prepare namelist output if requested
    if (nmlout.and.lroot) then
      open(unitreagentsout,file=path(2)(1:length(2))//'REAGENTS.namelist',access='append',status='replace',iostat=readerror)
      if (readerror.ne.0) then
        write(*,*) '#### FLEXPART MODEL ERROR CANNOT CREATE FILE:   ####'
        write(*,*) '#### ',path(2)(1:length(2))//'REAGENTS.namelist ####'
      endif
    endif

    ! try namelist input
    read(unitreagents,reagent_params,iostat=readerror)
    close(unitreagents)
    if (readerror.ne.0) then
      ! not namelist format
      nreagent=0
      go to 999
    endif

    ! namelist input
    open(unitreagents,file=path(1)(1:length(1))//'REAGENTS',status='old',form='formatted',iostat=readerror)
    j=0
    do while (readerror.eq.0)
      j=j+1
      if (j.gt.maxreagent) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY REAGENTS #### '
        write(*,*) ' #### MAXIMUM NUMBER IS ',maxreagent,'        #### '
        write(*,*) ' #### PLEASE MAKE CHANGES IN FILE REAGENTS    #### '
        stop
      endif
      read(unitreagents,reagent_params,iostat=readerror)
      reagents(j)=preagent
      reag_path(j)=preag_path
      reag_hourly(j)=phourly
      reag_unit(j)=punit
      ! namelist output
      if (nmlout.and.lroot) then
        write(unitreagentsout,nml=reagent_params)
      endif
    end do
    nreagent=j-1
    if (lroot) then
      write(*,*) 'Number of reagents: ',nreagent
      write(*,*) 'Reagent names: ',reagents(1:nreagent)
      write(*,*) 'Reagent units: ',reag_unit(1:nreagent)
    endif
    close(unitreagents)
    close(unitreagentsout)

999 continue ! no reagents file

  end subroutine readreagents

  subroutine getchemfield(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine reads the chemical reagent fields into memory              *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! itime            time since start of simulation in sec                     *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer          :: itime
    integer          :: jjjjmmdd, hhmmss, mm, eomday
    integer          :: nr, memid
    real(kind=dp)    :: jd, jdmid
    character(len=2) :: amonth
    character(len=256):: CL_name, jr_name
    logical          :: lexist

    print*, 'getchemfield: ldirect*memCLtime(1) = ',ldirect*memCLtime(1)
    print*, 'getchemfield: ldirect*memCLtime(2) = ',ldirect*memCLtime(2)

    ! Check fields are available for the current time step
    !*****************************************************

    if ((ldirect*memCLtime(1).le.ldirect*itime).and. &
         (ldirect*memCLtime(2).gt.ldirect*itime)) then

      ! The rightfields are already in memory -> don't do anything
      return 

    else if ((ldirect*memCLtime(2).le.ldirect*itime).and. &
         (memCLtime(2).ne.0.)) then

      ! Current time is after 2nd chem field
      !*************************************

      write(*,*) 'Updating CL fields... '

      memCLtime(1)=memCLtime(2)  ! time in sec
      CL_time(1)=CL_time(2)      ! julian date
      memid=2

      ! Get date/time of next chem field 
      !*********************************
      ! assumes fields are monthly
     
      call caldate(CL_time(1), jjjjmmdd, hhmmss)
      eomday=calceomday(jjjjmmdd/100)
      memCLtime(2)=memCLtime(1)+ldirect*eomday*24*3600   ! time in sec
      CL_time(2)=CL_time(1)+real(ldirect*eomday,kind=dp) ! julian date 
      write(*,*) 'getchemfield: memCLtime = ',memCLtime(1), memCLtime(2)
      write(*,*) 'getchemfield: CL_time = ',CL_time(1), CL_time(2)
      call caldate(CL_time(2), jjjjmmdd,hhmmss)
      mm=int((jjjjmmdd-(jjjjmmdd/10000)*10000)/100)
      write(amonth,'(I2.2)') mm

      do nr=1,nreagent

        write(*,*) 'Updating CL field for: ',trim(reagents(nr))

        CL_field(:,:,:,nr,1)=CL_field(:,:,:,nr,2)

        ! Read new chem field and store in 2nd position
        !**********************************************

        write(CL_name,'(A)') trim(reag_path(nr))//trim(reagents(nr))//'_'//amonth//'.nc'
        inquire(file=CL_name,exist=lexist)
        if (lexist) then
          write(*,*) 'Reading CL field: ',CL_name
          call readchemfield(CL_name, memid, nr)
        else
          write(*,*) '#### FLEXPART ERROR                ####'
          write(*,*) '#### CHEMISTRY FIELD NOT FOUND     ####'
          write(*,*) '#### '//trim(CL_name)//' ####'
          error stop
        endif

      end do ! nreagent

    else

      ! No chem field in memory that can be used
      !******************************************  

      write(*,*) 'Reading two new CL fields...'

      ! Get date/time of both chem fields 
      !**********************************
      ! assumes fields are monthly

      do memid=1,2
        if (memid.eq.1) then
          jd=bdate+real(ldirect*itime,kind=dp)/86400._dp
          call caldate(jd, jjjjmmdd, hhmmss)
          ! middle of month day
          jdmid=juldate(int(jjjjmmdd/100)*100+15,0) 
          !! testing
          print*, 'getchemfield: jjjjmmdd, jjjjmmdd_mid = ',jjjjmmdd, int(jjjjmmdd/100)*100+15
          ! julian date of fields
          if (ldirect.gt.0) then
            if (jd.ge.jdmid) then
              ! use current month
              CL_time(memid)=jdmid 
            else
              ! use last month
              eomday=calceomday(jjjjmmdd/100)
              CL_time(memid)=jdmid-real(eomday,kind=dp)
            endif
          else
            if (jd.ge.jdmid) then
              ! use next month
              CL_time(memid)=jdmid+real(eomday,kind=dp)
            else
              ! use current month
              CL_time(memid)=jdmid
            endif
          endif
        else
          call caldate(jd, jjjjmmdd, hhmmss)
          eomday=calceomday(jjjjmmdd/100)
          CL_time(memid)=CL_time(memid-1)+real(ldirect*eomday,kind=dp)
        endif
        ! time of field in seconds from start
        memCLtime(memid)=int((CL_time(memid)-bdate)*86400._dp)

        write(*,*) 'getchemfield: memid, memCLtime = ',memCLtime(memid)
        write(*,*) 'getchemfield: memid, CL_time = ',CL_time(memid)

        call caldate(CL_time(memid), jjjjmmdd, hhmmss)
        mm=int((jjjjmmdd-(jjjjmmdd/10000)*10000)/100)
        write(amonth,'(I2.2)') mm

        do nr=1,nreagent

          write(*,*) 'Reading two new CL fields for: ',trim(reagents(nr))

          ! Read new chem field
          !********************

          write(CL_name,'(A)') trim(reag_path(nr))//trim(reagents(nr))//'_'//amonth//'.nc'
          inquire(file=CL_name,exist=lexist)
          if (lexist) then
            write(*,*) 'Reading CL field: ',CL_name
            call readchemfield(CL_name, memid, nr)
          else
            write(*,*) '#### FLEXPART ERROR                ####'
            write(*,*) '#### CHEMISTRY FIELD NOT FOUND     ####'
            write(*,*) '#### '//trim(CL_name)//' ####'
            error stop
          endif
 
        end do ! nreagent

      end do ! memid
      
    endif ! update hourly fields


  end subroutine getchemfield


  subroutine readchemfield(CL_name,memid,nr)

  !*****************************************************************************
  !                                                                            *
  !    Reads the chemical reagent fields into memory                           *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! CL_name           name of chemical reagent file                            *
  ! memid             time index to chemical field variable                    *
  ! nr                reagent index to chemical field variable                 *
  !                                                                            * 
  !*****************************************************************************

    implicit none

    character(len=256) :: CL_name
    integer            :: memid,nr,len
    integer            :: ncid,dimid,varid

    ! Read netcdf file
    !******************

    ! open file
    call nf90_err( nf90_open(trim(CL_name),nf90_NOWRITE,ncid) )

    ! longitude
    call nf90_err( nf90_inq_dimid(ncid,'lon',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=len) )
    if (.not.allocated(lonCL)) then
      allocate(nxCL(nreagent))
      allocate(dxCL(nreagent))
      allocate(lonCL(nxmax,nreagent))
    endif
    nxCL(nr)=len
    call nf90_err( nf90_inq_varid(ncid,'lon',varid) )
    call nf90_err( nf90_get_var(ncid,varid,lonCL(1:nxCL(nr),nr)) )
    dxCL(nr)=abs(lonCL(2,nr)-lonCL(1,nr))

    ! latitude
    call nf90_err( nf90_inq_dimid(ncid,'lat',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=len) )
    if (.not.allocated(latCL)) then
      allocate(nyCL(nreagent))
      allocate(dyCL(nreagent))
      allocate(latCL(nymax,nreagent))
    endif
    nyCL(nr)=len
    call nf90_err( nf90_inq_varid(ncid,'lat',varid) )
    call nf90_err( nf90_get_var(ncid,varid,latCL(1:nyCL(nr),nr)) )
    dyCL(nr)=abs(latCL(2,nr)-latCL(1,nr))

    ! altitude
    call nf90_err( nf90_inq_dimid(ncid,'lev',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=len) )
    if (.not.allocated(altCL)) then
      allocate(nzCL(nreagent))
      allocate(altCL(nzmax,nreagent))
    endif
    nzCL(nr)=len
    call nf90_err( nf90_inq_varid(ncid,'lev',varid) )
    call nf90_err( nf90_get_var(ncid,varid,altCL(1:nzCL(nr),nr)) )

    ! chemical field
    call nf90_err( nf90_inq_varid(ncid,trim(reagents(nr)),varid) )
    if (.not.allocated(CL_field)) then
      allocate(CL_field(nxmax,nymax,nzmax,nreagent,2))
      allocate(clfield_cur(nxmax,nymax,nzmax,nreagent))
    endif
    call nf90_err( nf90_get_var(ncid,varid,CL_field(1:nxCL(nr),1:nyCL(nr),1:nzCL(nr),nr,memid)) )

    ! close file
    call nf90_err( nf90_close(ncid) )

    !! testing
!    print*, 'readchemfield: nxCL, nyCL, nzCL = ',nxCL(nr), nyCL(nr), nzCL(nr)
!    print*, 'readchemfield: lonCL = ',lonCL(1:nxCL(nr),nr)
!    print*, 'readchemfield: latCL = ',latCL(1:nyCL(nr),nr)

    return

  end subroutine readchemfield

  subroutine getchemhourly(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine interpolates the chemistry fields to current hour and      *
  !    if required using information on solar zenith angle                     *
  !                                                                            *
  !    Author: Rona Thompson, Mar 2024                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            actual simulation time [s]                            *
  !                                                                            *
  !*****************************************************************************

    use point_mod,      only: xlon0, ylat0, dx, dy
    use windfields_mod, only : height, nz

    implicit none

    integer            :: itime, curhour, interp_time
    integer            :: nr, kz, ix, jy, ixm, jym, ixp, jyp, indz, indzm, ii
    real               :: dt1, dt2, dtt, sza, jrate, r
    real, dimension(2) :: r1
    real               :: dz1, dz2, dz, ddx, ddy, rddx, rddy, p1, p2, p3, p4
    integer            :: jjjjmmdd, hhmmss, hh
    real(kind=dp)      :: jul, curday
    real, parameter    :: avog = 6.02214e23 ! Avogadro constant (1/mol)
    !! testing
    character(len=4) :: atime
    real, dimension(nxCL(1),nyCL(1)) :: jscalar
    character(len=20) :: frmt

    !! testing
    jscalar(:,:)=0.

    ! current hour of simulation
    curhour=itime/3600
    write(*,*) 'getchemhourly: curhour, curCLhour = ',curhour, curCLhour

    if ((ldirect*curCLhour.eq.ldirect*curhour).and.(ldirect*itime.gt.0)) then
      ! chemical field is for current hour -> don't do anything
      return
    else

      ! interpolate chemical field
      !*****************************

      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,hhmmss)
      ! interpolate to middle of hour
      curCLhour=curhour
      interp_time=curhour*3600+1800
      dt1=float(interp_time-memCLtime(1))
      dt2=float(memCLtime(2)-interp_time)
      dtt=1./(dt1+dt2)

      if (any(reag_hourly.gt.0)) then
        ! interpolate using rate of O3 + hv -> O1D 
        ! check if daily mean rate exists for current day
        ! note: assume that all chemical fields using hourly interpolation
        ! are on the same grid
        do nr=1,nreagent
          if (reag_hourly(nr).gt.0) exit
        end do
        curday=bdate+real(itime,kind=dp)/86400._dp
        curday=floor(curday)
        write(*,*) 'getchemhourly: curday, memCLday = ',curday, memCLday
        if (memCLday.ne.curday) then
          if (.not.allocated(jrate_average)) then
            allocate(jrate_average(nxCL(1),nyCL(1)))
          endif
          memCLday=curday
          jrate_average(:,:)=0.
!$OMP PARALLEL &
!$OMP PRIVATE(jy,ix,hh,sza,jrate) &
!$OMP REDUCTION(+:jrate_average)
!$OMP DO
          do jy=1,nyCL(nr)
            do ix = 1,nxCL(nr)
              do hh=0,23
                sza=zenithangle(latCL(jy,nr),lonCL(ix,nr),jjjjmmdd,hh*10000)
                jrate=photo_O1D(sza)
                jrate_average(ix,jy)=jrate_average(ix,jy)+jrate
              end do
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
          jrate_average=jrate_average/24.
        end if
        !! testing
!        print*, 'getchemhourly: range(jrate_average) = ',minval(jrate_average),maxval(jrate_average)
      endif

      ! loop over reagents
      do nr=1,nreagent

        write(*,*) 'Interpolating fields for reagent: ',reagents(nr)
        clfield_cur(:,:,:,nr)=(dt2*CL_field(:,:,:,nr,1) + dt1*CL_field(:,:,:,nr,2))*dtt

        if (reag_unit(nr).eq.'mol/mol') then
          ! convert to molecule/cm3
!$OMP PARALLEL &
!$OMP PRIVATE(kz,ix,jy,indz,indzm,ixm,jym,ixp,jyp,ddx,ddy,rddx,rddy,&
!$OMP p1,p2,p3,p4,r1,ii,dz1,dz2,dz,r)
!$OMP DO
          do kz=1,nzCL(nr)
            ! assume chem fields vertical coordinate is in meters
            indzm=nz-1
            do indz=1,nz
              if (height(indz).gt.altCL(kz,nr)) then
                indzm=indz-1
                exit
              endif
            end do
            dz1=altCL(kz,nr)-height(indzm)
            dz2=height(indz)-altCL(kz,nr)
            dz=1./(dz1+dz2)
            do jy=1,nyCL(nr)
              jym=int((latCL(jy,nr)-ylat0)/dy)
              jyp=jym+1
              ddy=(latCL(jy,nr)-ylat0)/dy-real(jym)
              rddy=1.-ddy
              do ix=1,nxCL(nr)
                ixm=int((lonCL(ix,nr)-xlon0)/dx)
                ixp=ixm+1
                ddx=(lonCL(ix,nr)-xlon0)/dx-real(ixm)
                rddx=1.-ddx
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy
                ! take density from first field (accurate enough)
                do ii=1,2
                  indz=indzm+ii-1
                  r1(ii)=p1*rho(ixm,jym,indz,1)+&
                         p2*rho(ixp,jym,indz,1)+&
                         p3*rho(ixm,jyp,indz,1)+&
                         p4*rho(ixp,jyp,indz,1)
                end do
                r=(dz2*r1(1)+dz1*r1(2))*dz
                ! vmr*Avog*P/(RT)/1e6
                ! using P/T = rho*rair
                clfield_cur(ix,jy,kz,nr)=clfield_cur(ix,jy,kz,nr)*avog*r*r_air/rgas/1.e6 
              end do
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
        endif

        if (reag_hourly(nr).gt.0) then
          ! interpolate using rate of O3 + hv -> O1D 
!$OMP PARALLEL &
!$OMP PRIVATE(kz,jy,ix,sza,jrate) 
!$OMP DO
          do kz=1,nzCL(nr)
            do jy=1,nyCL(nr)
              do ix=1,nxCL(nr)
                ! solar zenith angle
                sza=zenithangle(latCL(jy,nr),lonCL(ix,nr),jjjjmmdd,hhmmss)
                ! calculate J(O1D) (jrate)
                jrate=photo_O1D(sza)
                ! hourly correction
                if (jrate_average(ix,jy).ne.0) then
                  clfield_cur(ix,jy,kz,nr)=clfield_cur(ix,jy,kz,nr)*jrate/jrate_average(ix,jy)
                  !! testing
!                  jscalar(ix,jy)=jrate/jrate_average(ix,jy)
                endif
              end do
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
        endif ! reag_hourly
      end do ! nreagent
    endif ! curhour

    !! testing
!    write(atime,fmt='(I4.4)') curhour
!    write(frmt,fmt='(A,I4,A)') '(',nxCL(1),'(E14.6))'
!    open(600,file=path(2)(1:length(2))//'jscalar_'//atime//'.txt',action='write',status='replace')
!    do jy=1,nyCL(1)
!      write(600,fmt=frmt) jscalar(:,jy)
!    end do
!    close(600)
    !!

  end subroutine getchemhourly


  subroutine chemreaction(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine computes the chemical loss of each species and updates     *
  !    the particle mass                                                       *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !    updated for v11, Jan 2024                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            actual simulation time [s]                            *
  !                                                                            *
  !*****************************************************************************

    use particle_mod,   only: count, part, mass
    use point_mod,      only: xlon0, ylat0, dx, dy
    use windfields_mod, only: xresoln, yresoln, xln, yln, xrn, yrn, height, tt, nz
    use omp_lib

    implicit none

    integer               :: jpart,itime
    integer               :: ii,i,j,ks,ix,jy,kz,jrx,jry,nr
    integer               :: ngrid,interp_time,n,indz
    integer               :: jjjjmmdd,hhmmss,mm,hh,m1,m2
    integer               :: clx,cly,clz,clzm,ithread
    real, dimension(nzmax) :: altCLtop
    real, dimension(2)    :: cl_tmp
    real                  :: xlon,ylat
    real                  :: xtn,ytn
    real                  :: dt1,dt2,dtt,dtt1,dtt2,dttt,dz1,dz2,dzz
    real                  :: restmass,clreacted,cl_cur
    real                  :: jrate,jrate_cur,sza
    real                  :: clrate,temp
    real, parameter       :: smallnum = tiny(0.0) ! smallest number that can be handled
    real(kind=dp)         :: jul
    real                  :: lonjrx,latjry
#ifdef _OPENMP
    real(kind=dp), allocatable, dimension(:,:,:) :: chem_loss_tmp
#endif

    ! use middle of synchronisation time step
    interp_time=nint(itime+0.5*lsynctime)
    dtt1=float(interp_time-memtime(1))
    dtt2=float(memtime(2)-interp_time)
    dttt=1./(dtt1+dtt2)

    ! initialization
    chem_loss(:,:)=0d0

#ifdef _OPENMP
    allocate( chem_loss_tmp(nreagent,nspec,numthreads) )
    chem_loss_tmp(:,:,:) = 0d0
#endif

    ! Loop over particles
    !*****************************************

!$OMP PARALLEL &
!$OMP PRIVATE(ii,jpart,ngrid,j,xtn,ytn,ix,jy, &
!$OMP   xlon,ylat,clx,cly,clz,clzm,kz,altCLtop,dz1,dz2,dzz,nr,i, &
!$OMP   cl_cur,indz,temp,ks,clrate,restmass,clreacted,ithread) 

#ifdef _OPENMP
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    ithread = 1
#endif

!$OMP DO
    do ii=1,count%alive

      jpart=count%ialive(ii)
      ! Determine which nesting level to be used
      ngrid=0
      do j=numbnests,1,-1 ! check if need +/- dxn below
        if ((part(jpart)%xlon.gt.xln(j)).and.(part(jpart)%xlon.lt.xrn(j)).and. &
             (part(jpart)%ylat.gt.yln(j)).and.(part(jpart)%ylat.lt.yrn(j))) then
          ngrid=j
          exit
        endif
      end do

      ! Determine nested grid coordinates
      if (ngrid.gt.0) then
        xtn=(real(part(jpart)%xlon)-xln(ngrid))*xresoln(ngrid)
        ytn=(real(part(jpart)%ylat)-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
      else
        ix=int(part(jpart)%xlon)
        jy=int(part(jpart)%ylat)
      endif

      ! Get CL from nearest grid-cell for current hour
      !***********************************************

      ! world coordinates
      xlon=real(part(jpart)%xlon)*dx+xlon0
      if (xlon.gt.180.) then
        xlon=xlon-360.
      endif
      ylat=real(part(jpart)%ylat)*dy+ylat0
      if (ylat.gt.90.) then
        ylat=90.
      endif

      do nr=1,nreagent

        ! get position in the chem field
        ! assumes chem field dimensions given as grid midpoints
        clx=min(nxCL(nr),int((xlon-(lonCL(1,nr)-0.5*dxCL(nr)))/dxCL(nr))+1)
        cly=min(nyCL(nr),int((ylat-(latCL(1,nr)-0.5*dyCL(nr)))/dyCL(nr))+1)   

        ! get the level of the chem field for the particle
        ! z is the z-coord of the trajectory above model orography in metres
        ! altCL is the height of the centre of the level in the chem field above orography
        do kz=2,nzCL(nr)
          altCLtop(kz-1)=altCL(kz-1,nr)+0.5*(altCL(kz,nr)-altCL(kz-1,nr))
        end do
        altCLtop(nzCL(nr))=altCL(nzCL(nr),nr)+0.5*(altCL(nzCL(nr),nr)-altCL(nzCL(nr)-1,nr))
        clzm=nzCL(nr)-1
        do clz=1,nzCL(nr)
          if (real(part(jpart)%z).lt.altCLtop(clz)) then
            clzm=clz-1
            exit
          endif
        end do
        clz=min(nzCL(nr),clz)
        if (clzm.eq.0 ) then
          dz1=1.
          dz2=1.
          clzm=clz
        else
          dz1=real(part(jpart)%z)-altCL(clzm,nr)
          dz2=altCL(clz,nr)-real(part(jpart)%z)
        endif
        if (dz1.eq.(-1.*dz2)) then
          dz1=1.
          dz2=1.
        endif
        dzz=1./(dz1+dz2)

        !! testing
!        if(ii.lt.100) print*, 'chemreaction: nr, clx, cly, clz, clzm = ',nr, clx, cly, clz, clzm

        ! chem reagent for particle time and location 
        cl_cur=(dz2*clfield_cur(clx,cly,clzm,nr) + &
                dz1*clfield_cur(clx,cly,clz,nr))*dzz

        ! Compute chemical loss
        !**********************

        if (cl_cur.gt.smallnum) then
          indz=nz
          do kz=2,nz
            if (height(kz).gt.part(jpart)%z) then
              indz=kz-1
              exit
            endif
          end do
          temp=(tt(ix,jy,indz,1)*dtt2 &
                 + tt(ix,jy,indz,2)*dtt1)*dttt
          do ks=1,nspec
            if (reaccconst(nr,ks).gt.0.) then
              ! k = CT^Nexp(-D/T)[reagent]
              clrate=reaccconst(nr,ks)*(temp**reacnconst(nr,ks))*exp(-1.*reacdconst(nr,ks)/temp)*cl_cur
              ! new particle mass
              restmass=mass(jpart,ks)*exp(-1.*clrate*abs(lsynctime))
              if (restmass.gt.smallnum) then
                clreacted=mass(jpart,ks)-restmass
                mass(jpart,ks)=restmass
              else
                clreacted=mass(jpart,ks)
                mass(jpart,ks)=0.
              endif
#ifdef _OPENMP
              chem_loss_tmp(nr,ks,ithread)=chem_loss_tmp(nr,ks,ithread)+real(clreacted,kind=dp)
#else
              chem_loss(nr,ks)=chem_loss(nr,ks)+real(clreacted,kind=dp)
#endif
            endif
          end do  ! nspec
        endif   

      end do  ! nreagent

    end do  ! loop over all particles
!$OMP END DO

!$OMP END PARALLEL

#ifdef _OPENMP
  do ithread=1,numthreads
    chem_loss(:,:) = chem_loss(:,:)+chem_loss_tmp(:,:,ithread)
  end do
  deallocate( chem_loss_tmp )
#endif

  end subroutine chemreaction


  function photo_O1D(sza) result(jrate)

  !*****************************************************************************
  !                                                                            *
  !    Calculates J(O1D) photolysis rate based on solar zenith angle           *
  !                                                                            *
  !    Ref: Hamer, P. D., et al. Atmos. Chem. Phys. 2015
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !    INPUT:                                                                  *
  !    sza        solar zenith angle (degrees)                                 *
  !                                                                            *
  !    OUTPUT:                                                                 *
  !    photo_O1D  J(O1D) photoylsis rate                                       *
  !                                                                            *
  !*****************************************************************************

    implicit none

    real, intent(in)  :: sza
    real :: sun, jrate
    real, parameter :: pi=3.1415927

    sun=sza*pi/180.
    if (sun.le.(pi/2.)) then
      sun=cos(sun)
    else
      sun=0.
    endif

    jrate=(8.978e-5*(sun**1.436)*exp(-0.936/sun))

  end function photo_O1D


  function zenithangle(ylat,xlon,jjjjmmdd,hhmmss) result(sza)

  !*****************************************************************************
  !                                                                            *
  !    This function returns the sinus of solar elevation as a function        *
  !    of geographic longitde, latitude and GMT-Time.                          *
  !                                                                            *
  !    Author: G. Wotawa, Nov 1993                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !     INPUT:                                                                 *
  !                                                                            *
  !            ylat          geographical latitude  [DEG]                      *
  !            xlon          geographical longitude [DEG]                      *
  !            jjjj          Year                                              *
  !            mm            Month                                             *
  !            dd            Day                                               *
  !            hh            Hour                                              *
  !            minute        Minute                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer,       intent(in)  :: jjjjmmdd,hhmmss
    real,          intent(in)  :: ylat,xlon
    real    :: sza
    integer :: jjjj,mm,id,iu,minute 
    integer :: ndaynum
    real :: sinsol,solelev
    real :: rnum,rylat,ttime,dekl,rdekl,eq
    real,parameter :: pi=3.1415927

    jjjj=jjjjmmdd/10000
    mm=jjjjmmdd/100-jjjj*100
    id=jjjjmmdd-jjjj*10000-mm*100
    iu=hhmmss/10000
    minute=hhmmss/100-100*iu

    ndaynum=31*(mm-1)+id
    if(mm.gt.2) ndaynum=ndaynum-int(0.4*mm+2.3)
    if((mm.gt.2).and.(jjjj/4*4.eq.jjjj)) ndaynum=ndaynum+1

    rnum=2.*pi*ndaynum/365.
    rylat=pi*ylat/180.
    ttime=real(iu)+real(minute)/60.

    dekl=0.396+3.631*sin(rnum)+0.038*sin(2.*rnum)+0.077*sin(3.*rnum)- &
         22.97*cos(rnum)-0.389*cos(2.*rnum)-0.158*cos(3.*rnum)
    rdekl=pi*dekl/180.

    eq=(0.003-7.343*sin(rnum)-9.47*sin(2.*rnum)- &
         0.329*sin(3.*rnum)-0.196*sin(4.*rnum)+ &
         0.552*cos(rnum)-3.020*cos(2.*rnum)- &
         0.076*cos(3.*rnum)-0.125*cos(4.*rnum))/60.
    sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)* &
         cos((ttime-12.+xlon/15.+eq)*pi/12.)
    ! Calculate the maximum solar elevation on that day
    !sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)*
    !    &       cos((eq)*pi/12.)
    solelev=asin(sinsol)*180./pi
    sza=90.-solelev

  end function zenithangle


end module chemistry_mod

