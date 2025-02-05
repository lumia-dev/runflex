! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module emissions_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for injecting mass       *
  !    into particles based on gridded emissions estimates                     *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use point_mod,         only: xlon0, ylat0, dx, dy, npart
  use particle_mod
  use date_mod
  use totals_mod,        only: tot_em_up, tot_em_field, tot_em_res
  use netcdf
  use netcdf_output_mod, only: nf90_err

  implicit none

  real, parameter    :: tau_em_r=4000. ! time scale of residual emission release (s)
  real, parameter    :: tau_ipm=900.   ! time scale of micro mixing (s)
  logical, parameter :: ABLMIX=.true.  ! mass exchange for particles in same PBL grid cell

  integer                               :: nxem, nyem
  real                                  :: dxem, dyem
  integer, dimension(2)                 :: em_memtime 
  real, allocatable, dimension(:,:,:,:) :: em_field                ! emission fields all species
  real, allocatable, dimension(:)       :: lonem, latem            ! emission field lon and lat
  real, allocatable, dimension(:,:,:)   :: em_res                  ! residual emissions
  real, allocatable, dimension(:,:,:)   :: mass_field
  real, allocatable, dimension(:,:)     :: em_area                 ! area for emissions grid
  integer, allocatable, dimension(:,:)  :: nn_field
  real(kind=dp), dimension(2)           :: em_time
  character(len=32), allocatable, dimension(:) :: emf_name

  contains

  subroutine emissions(itime)

  !*****************************************************************************
  !                                                                            *
  !    This subroutine calculates the fraction of emission to be released      *
  !    in each timestep and adds this mass to the particles in the PBL         *
  !                                                                            *
  !    Author: S. Henne, Mar-2009                                              *
  !    Adapted by R. Thompson for v10.4, Oct-2023                              *
  !                                                                            *
  !*****************************************************************************

    use windfields_mod, only: hmix
    use omp_lib

    implicit none

    integer       :: itime
    real          :: xlon,ylat
    integer       :: ix, jy, ixp, jyp, ii, ks, em_ix, em_jy, ithread
    real          :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
    real          :: em_dt1, em_dt2, em_dtt
    real, dimension(2) :: hm 
    real          :: hmixi
    real          :: tmp, max_val, em_cur
    logical, dimension(npart(1)) :: em_cond
    integer       :: mm 
    character(len=20) :: frmt
    real          :: em_frac
    real, allocatable, dimension(:,:,:) :: em_neg
    real, allocatable, dimension(:,:) :: tot_em_up_tmp
    real, allocatable, dimension(:,:,:,:) :: mass_field_tmp, em_neg_tmp
    real          :: f_m
    logical       :: lexist
!    integer, parameter :: unittest=120


    ! fraction of stored emissions that is released in a time step
    em_frac = 1. - exp(-1.*real(lsynctime)/tau_em_r)

    ! distance of emission fields in memory from current time
    em_dt1 = float(itime-em_memtime(1))
    em_dt2 = float(em_memtime(2)-itime)
    em_dtt = 1./(em_dt1+em_dt2)

    ! determine temporal distances for interpolation of meteo fields
    dt1=float(itime-memtime(1))
    dt2=float(memtime(2)-itime)
    dtt=1./(dt1+dt2)

    tot_em_up(:) = 0.

    ! estimate mass in PBL from particle positions
    !**********************************************

    mass_field(:,:,:) = 0.
    allocate( em_neg(nspec-1,nxem,nyem) )
#if _OPENMP
    allocate( mass_field_tmp(nspec,nxem,nyem,numthreads))
    allocate( em_neg_tmp(nspec-1,nxem,nyem,numthreads))
    allocate( tot_em_up_tmp(nspec,numthreads) )
    mass_field_tmp(:,:,:,:) = 0.
    em_neg_tmp(:,:,:,:) = 0.
    tot_em_up_tmp(:,:) = 0.
#endif
    tot_em_up(:) = 0.
    mass_field(:,:,:) = 0.
    em_neg(:,:,:) = 0.

!$OMP PARALLEL &
!$OMP PRIVATE(ii,xlon,ylat,em_ix,em_jy,ix,jy,ixp,jyp,ddx,ddy, &
!$OMP rddx,rddy,p1,p2,p3,p4,mm,hm,hmixi,ks,ithread) 

#ifdef _OPENMP
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    ithread = 1
#endif

!$OMP DO 
    do ii=1,count%alive  ! loop over all particles

      xlon=xlon0+part(ii)%xlon*dx
      ylat=ylat0+part(ii)%ylat*dy

      ! assume emission dimensions given as grid midpoints
      em_ix=min(nxem, int((xlon-(lonem(1)-0.5*dxem))/dxem)+1)
      em_jy=min(nyem, int((ylat-(latem(1)-0.5*dyem))/dyem)+1)
      !! testing
!      if (ii.lt.20) print*, 'lonem, lon, latem, lat = ',lonem(em_ix),xlon,latem(em_jy),ylat

      ! interpolate to particle position
      ix=int(part(ii)%xlon)
      jy=int(part(ii)%ylat)
      ixp=ix+1
      jyp=jy+1
      ddx=part(ii)%xlon-float(ix)
      ddy=part(ii)%ylat-float(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

      do mm=1,2
        ! PBL height at particle position
        hm(mm)=p1*hmix(ix,jy ,1,memind(mm)) + &
               p2*hmix(ixp,jy ,1,memind(mm)) + &
               p3*hmix(ix ,jyp,1,memind(mm)) + &
               p4*hmix(ixp,jyp,1,memind(mm))
      end do
      hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
      ! set minimum PBL height to dampen day/night amplitude of emission
      hmixi=max(hmixi,hmixmin)

      ! determine if particle is in PBL
      em_cond(ii) = part(ii)%z.le.hmixi

      if (em_cond(ii)) then
#ifdef _OPENMP
        mass_field_tmp(1:nspec,em_ix,em_jy,ithread)= &
              mass_field_tmp(1:nspec,em_ix,em_jy,ithread) + &
              mass(ii,1:nspec)
#else
        mass_field(1:nspec,em_ix,em_jy)=mass_field(1:nspec,em_ix,em_jy) + &
                                        mass(ii,1:nspec)
#endif
      endif

    end do   ! end of particle loop
!$OMP END DO
!$OMP END PARALLEL

#ifdef _OPENMP
    ! manual reduction of mass_field
    do ithread=1,numthreads
      mass_field(:,:,:) = mass_field(:,:,:)+mass_field_tmp(:,:,:,ithread)
    end do
#endif

    f_m = exp(-1.*real(lsynctime)/tau_ipm)

    ! estimate emissions for each particle
    !**************************************

!$OMP PARALLEL &
!$OMP PRIVATE(ii,xlon,ylat,em_ix,em_jy, &
!$OMP ks,em_cur,tmp,ithread) 

#ifdef _OPENMP
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    ithread = 1
#endif

!$OMP DO
    do ii=1,count%alive ! loop over particles

      if (.not.em_cond(ii)) cycle ! skip particles not in PBL

      xlon=xlon0+part(ii)%xlon*dx
      ylat=ylat0+part(ii)%ylat*dy
      ! assume emission dimensions given as grid midpoints
      em_ix=min(nxem, int((xlon-(lonem(1)-0.5*dxem))/dxem)+1)
      em_jy=min(nyem, int((ylat-(latem(1)-0.5*dyem))/dyem)+1)

      ! loop over species
      ! skip species 1 as it is always air tracer with no emission
      do ks=2,nspec
        em_cur=(em_field(ks-1,em_ix,em_jy,1)*em_dt2 + &
                em_field(ks-1,em_ix,em_jy,2)*em_dt1)*em_dtt
        tmp=(em_cur*real(lsynctime) + em_res(ks-1,em_ix,em_jy)*em_frac ) * &
               mass(ii,1)/mass_field(1,em_ix,em_jy)

        ! micro mixing scheme executed before new emissions are taken up
        !*****************************************************************
        if (ABLMIX) then
          mass(ii,ks)=mass(ii,1) * &
                         (f_m * mass(ii,ks)/mass(ii,1) + &
                         (1.-f_m) * mass_field(ks,em_ix,em_jy)/mass_field(1,em_ix,em_jy))
        endif

        ! deal with negative emissions
        if (tmp.lt.0.) then
          if (-1.*tmp.gt.mass(ii,ks)) then
            tmp = tmp + mass(ii,ks)
            ! subtract mass from atmosphere by setting it to zero below
#ifdef _OPENMP
            tot_em_up_tmp(ks,ithread) = tot_em_up_tmp(ks,ithread) - real(mass(ii,ks),kind=dp)
#else
            tot_em_up(ks) = tot_em_up(ks) - real(mass(ii,ks),kind=dp)
#endif
            mass(ii,ks) = 0.
            ! add remaining uptake to em_neg 
#ifdef _OPENMP
            em_neg_tmp(ks-1,em_ix,em_jy,ithread) = em_neg_tmp(ks-1,em_ix,em_jy,ithread) + &
              tmp/mass(ii,1)*mass_field(1,em_ix,em_jy)
#else
            em_neg(ks-1, em_ix, em_jy) = em_neg(ks-1, em_ix, em_jy) + &
              tmp/mass(ii,1)*mass_field(1,em_ix,em_jy)
#endif
          else
            mass(ii,ks)=mass(ii,ks)+tmp
#ifdef _OPENMP
            tot_em_up_tmp(ks,ithread) = tot_em_up_tmp(ks,ithread) + real(tmp,kind=dp)
#else
            tot_em_up(ks) = tot_em_up(ks) + real(tmp,kind=dp)
#endif
          endif
        else
          mass(ii,ks)=mass(ii,ks)+tmp
#ifdef _OPENMP
          tot_em_up_tmp(ks,ithread) = tot_em_up_tmp(ks,ithread) + real(tmp,kind=dp)
#else
          tot_em_up(ks) = tot_em_up(ks) + real(tmp,kind=dp)
#endif
        endif ! negative emissions

      end do ! nspec

    end do ! loop over particles
!$OMP END DO

    ! loop over grid points to update residual emissions 
    ! skip species 1 as it is always air tracer
    !***************************************************

!$OMP DO 
    do jy=1,nyem
      do ix=1,nxem
        do ks=2,nspec
          if (mass_field(1,ix,jy).eq.0.) then
            em_cur=(em_field(ks-1,ix,jy,1)*em_dt2 + &
                    em_field(ks-1,ix,jy,2)*em_dt1)*em_dtt
            em_res(ks-1,ix,jy)=em_res(ks-1,ix,jy) + &
                             em_cur*real(lsynctime)
          else
            em_res(ks-1,ix,jy)=(1.-em_frac) * em_res(ks-1,ix,jy)
          endif
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! manual reduction of em_neg and tot_em_up
#ifdef _OPENMP
    do ithread=1,numthreads
      em_neg(:,:,:) = em_neg(:,:,:)+em_neg_tmp(:,:,:,ithread)
      tot_em_up(:) = tot_em_up(:) + tot_em_up_tmp(:,ithread)
    end do
    deallocate( mass_field_tmp, em_neg_tmp, tot_em_up_tmp )
#endif

    ! update for negative emissions
    em_res(:,:,:) = em_res(:,:,:)+em_neg(:,:,:)
    deallocate(em_neg)

    ! calculate total emission flux and field
!$OMP PARALLEL IF(nspec>99) PRIVATE(ks)
!$OMP DO
    do ks=2,nspec
      tot_em_field(ks)=sum((em_field(ks-1,:,:,1)*em_dt2 + &
                         em_field(ks-1,:,:,2)*em_dt1)*em_dtt)*real(lsynctime)
      tot_em_res(ks)=sum(em_res(ks-1,:,:))
    end do
!$OMP END DO
!$OMP END PARALLEL

    !! test
!    inquire(file=path(2)(1:length(2))//'mass_field.txt',exist=lexist)
!    write(*,*) 'emissions: lexist = ',lexist
!    if (lexist) then
!      open(unittest,file=path(2)(1:length(2))//'mass_field.txt',ACCESS='APPEND',STATUS='OLD')
!    else
!      open(unittest,file=path(2)(1:length(2))//'mass_field.txt',STATUS='NEW')
!    endif
!    write(frmt, '(A, I4, A)') '(', nxem, 'E12.3)'
!    do jy=1,nyem-1
!      write(unittest,frmt) (mass_field(2,ix,jy),ix=1,nxem)
!    end do
!    close(unittest)
    !!

  end subroutine emissions

  subroutine getemissions(itime)

  !*****************************************************************************
  !                                                                            *
  !    This subroutine checks which emission fields need to be read and        *
  !    interpolates these to the current time step                             *
  !                                                                            *
  !    Author: Rona Thompson, Oct-2023                                         *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer       :: itime
    real(kind=dp) :: julstart, jd, jdmid
    integer       :: jjjjmmdd, hhmmss, dd, mm, yyyy
    integer       :: nn, ks, eomday, memid
    character(len=4)   :: ayear
    character(len=2)   :: amonth, aday
    character(len=256) :: em_name, file_name, strtmp1, strtmp2
    logical       :: lexist


    ! Check fields are available for the current time step
    !*****************************************************

    if ((ldirect*em_memtime(1).le.ldirect*itime).and. &
           (ldirect*em_memtime(2).gt.ldirect*itime)) then

      ! The rightfields are already in memory -> don't do anything
      return 

    else if ((ldirect*em_memtime(2).le.ldirect*itime).and. &
         (em_memtime(2).ne.0.)) then

      ! Current time is after 2nd field
      !*********************************

      ! update dates of emission fields
      em_memtime(1)=em_memtime(2)  ! time in sec
      em_time(1)=em_time(2)      ! julian date
      memid=2

      ! julian date of next emission field assuming monthly fields
      call caldate(em_time(1), jjjjmmdd, hhmmss)
      eomday=calceomday(jjjjmmdd/100)
      em_memtime(2)=em_memtime(1)+ldirect*eomday*24*3600   ! time in sec
      em_time(2)=em_time(1)+real(ldirect*eomday,kind=dp) ! julian date 
      call caldate(em_time(2), jjjjmmdd,hhmmss)
      yyyy=jjjjmmdd/10000
      mm=(jjjjmmdd-yyyy*10000)/100
      dd=jjjjmmdd-yyyy*10000-mm*100
      write(amonth,'(I2.2)') mm
      write(aday,'(I2.2)') dd
      write(ayear,'(I4)') yyyy

      ! skip species 1 as is always air tracer with no emissions
!$OMP PARALLEL IF(nspec>99) &
!$OMP PRIVATE(ks,file_name,nn,strtmp1,strtmp2,julstart,em_name,lexist)
!$OMP DO
      do ks=2,nspec

        write(*,*) 'Updating emission fields for: ',trim(species(ks)) 
    
        em_field(ks-1,:,:,1)=em_field(ks-1,:,:,2)

        ! Read new emission field and store in 2nd position
        !***************************************************

        ! TO DO: make more general to fit any field frequency

        file_name=emis_file(ks)
        nn=index(file_name,'YYYY',back=.false.)
        if (nn.ne.0) then
          strtmp1=file_name(1:nn-1)
          nn=index(file_name,'YYYY',back=.true.)
          strtmp2=file_name(nn+4:len_trim(file_name))
          file_name=trim(strtmp1)//ayear//trim(strtmp2)
          julstart=juldate((jjjjmmdd/10000)*10000+101,0)
        endif
        nn=index(file_name,'MM',back=.false.)
        if (nn.ne.0) then
          strtmp1=file_name(1:nn-1)
          nn=index(file_name,'MM',back=.true.)
          strtmp2=file_name(nn+2:len_trim(file_name))
          file_name=trim(strtmp1)//amonth//trim(strtmp2)
          julstart=juldate((jjjjmmdd/100)*100+1,0)
        endif
        nn=index(file_name,'DD',back=.false.)
        if (nn.ne.0) then
          strtmp1=file_name(1:nn-1)
          nn=index(file_name,'DD',back=.true.)
          strtmp2=file_name(nn+2:len_trim(file_name))
          file_name=trim(strtmp1)//aday//trim(strtmp2)
          julstart=juldate(jjjjmmdd,0)
        endif

        em_name=trim(emis_path(ks))//trim(file_name)
        inquire(file=em_name,exist=lexist)
        if (lexist) then
          write(*,*) 'Reading emissions field: ',trim(em_name)
          call reademissions(em_name, julstart, itime, memid, ks-1)
        else
          write(*,*) '#### FLEXPART ERROR                ####'
          write(*,*) '#### EMISSION FIELD NOT FOUND     ####'
          write(*,*) '#### '//trim(em_name)//' ####'
          error stop
        endif

      end do ! nspec
!$OMP END DO
!$OMP END PARALLEL

    else

      ! No emission field in memory that can be used
      !**********************************************

      ! read both fields into memory
      do memid=1,2

        if (memid.eq.1) then
          jd=bdate+real(ldirect*itime,kind=dp)/86400._dp 
          call caldate(jd, jjjjmmdd, hhmmss)
          ! middle of month day
          jdmid=juldate(int(jjjjmmdd/100)*100+15,0)
          if (jd.ge.jdmid) then
            ! use current month
            em_time(memid)=jdmid
          else
            ! use last month
            eomday=calceomday(jjjjmmdd/100)
            em_time(memid)=jdmid-real(eomday,kind=dp)          
          endif
        else
          call caldate(jd, jjjjmmdd, hhmmss)
          eomday=calceomday(jjjjmmdd/100)
          em_time(memid)=em_time(memid-1)+real(ldirect*eomday,kind=dp)
        endif
        em_memtime(memid)=int((em_time(memid)-bdate)*86400._dp)

        write(*,*) 'getemissions: memid, em_time = ',memid, em_time(memid)
        write(*,*) 'getemissions: memid, em_memtime = ',memid, em_memtime(memid)

        call caldate(em_time(memid), jjjjmmdd, hhmmss)
        yyyy=jjjjmmdd/10000
        mm=(jjjjmmdd-yyyy*10000)/100
        dd=jjjjmmdd-yyyy*10000-mm*100
        write(amonth,'(I2.2)') mm
        write(aday,'(I2.2)') dd
        write(ayear,'(I4)') yyyy

!$OMP PARALLEL IF(nspec>99) &
!$OMP PRIVATE(ks,file_name,nn,strtmp1,strtmp2,julstart,em_name,lexist)
!$OMP DO
        do ks=2,nspec

          write(*,*) 'Reading two new emission fields for: ',trim(species(ks))

          ! TO DO: make more general to fit any field frequency

          file_name=emis_file(ks)
          nn=index(file_name,'YYYY',back=.false.)
          if (nn.ne.0) then
            strtmp1=file_name(1:nn-1)
            nn=index(file_name,'YYYY',back=.true.)
            strtmp2=file_name(nn+4:len_trim(file_name))
            file_name=trim(strtmp1)//ayear//trim(strtmp2)
            julstart=juldate((jjjjmmdd/10000)*10000+101,0)
          endif
          nn=index(file_name,'MM',back=.false.)
          if (nn.ne.0) then
            strtmp1=file_name(1:nn-1)
            nn=index(file_name,'MM',back=.true.)
            strtmp2=file_name(nn+2:len_trim(file_name))
            file_name=trim(strtmp1)//amonth//trim(strtmp2)
            julstart=juldate((jjjjmmdd/100)*100+1,0)
          endif
          nn=index(file_name,'DD',back=.false.)
          if (nn.ne.0) then
            strtmp1=file_name(1:nn-1)
            nn=index(file_name,'DD',back=.true.)
            strtmp2=file_name(nn+2:len_trim(file_name))
            file_name=trim(strtmp1)//aday//trim(strtmp2)
            julstart=juldate(jjjjmmdd,0)
          endif

          em_name=trim(emis_path(ks))//trim(file_name)
          inquire(file=em_name,exist=lexist)
          if (lexist) then
            write(*,*) 'Reading emissions field: ',trim(em_name)
            call reademissions(em_name, julstart, itime, memid, ks-1)
          else
            write(*,*) '#### FLEXPART ERROR                ####'
            write(*,*) '#### EMISSION FIELD NOT FOUND     ####'
            write(*,*) '#### '//trim(em_name)//' ####'
            error stop
          endif

        end do ! nspec
!$OMP END DO
!$OMP END PARALLEL

      end do ! memid

    endif ! update fields


  end subroutine getemissions


  subroutine reademissions(em_name, julstart, itime, memid, kk)

  !*****************************************************************************
  !                                                                            *
  !    This subroutine reads the emission fields                               *
  !                                                                            *
  !    Author: Rona Thompson, Oct-2023                                         *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: em_name, unitinfo
    integer            :: memid, kk, itime
    real(kind=dp)      :: julstart, jtime
    integer            :: jjjjmmdd, ihmmss, mm
    integer            :: nn, ntem, len, ix, jy
    integer            :: ncid, dimid, varid, ndim, nvar, natt, xtype, indxt, unlimid
    real               :: ylatp, ylatm, hzone, cosfactm, cosfactp
    real               :: sclfact, offset
    integer, dimension(:), allocatable :: dimids
    real(kind=dp), dimension(:), allocatable :: jdate
    real, dimension(:), allocatable :: time
    real, dimension(:,:), allocatable :: emis
    character(len=32)  :: dimname, nameout, attname

    ! current time in julian days
    jtime=real(itime,kind=dp)/86400._dp+bdate

    ! Read netcdf file
    !******************

    ! open file
    call nf90_err( nf90_open(trim(em_name),nf90_NOWRITE,ncid) )

    ! inquire about dims and vars
    call nf90_err( nf90_inquire( ncid, ndim, nvar, natt, unlimid ) )
    allocate(dimids(ndim))

    if (.not.allocated(lonem)) then
      ! read dimension info
      do nn=1,ndim
        call nf90_err( nf90_inquire_dimension( ncid, nn, dimname, len ) )
        if ((index(dimname,'lon').ne.0).or.(index(dimname,'LON').ne.0) &
              .or.(index(dimname,'Lon').ne.0)) then
          ! longitude
          nxem=len
          allocate( lonem(nxem) )
          call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
          call nf90_err( nf90_get_var(ncid,varid,lonem) )
          dxem=abs(lonem(2)-lonem(1))
        endif
        if ((index(dimname,'lat').ne.0).or.(index(dimname,'LAT').ne.0) &
              .or.(index(dimname,'Lat').ne.0)) then
          ! latitude
          nyem=len
          allocate( latem(nyem) )
          call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
          call nf90_err( nf90_get_var(ncid,varid,latem) )
          dyem=abs(latem(2)-latem(1))
        endif
      end do ! ndim
      ! check dimensions read correctly
      if (.not.allocated(lonem)) then
        write(*,*) 'ERROR in reademissions: longitude dimension not found in file: '//trim(em_name)
        error stop
      endif
      if (.not.allocated(latem)) then
        write(*,*) 'ERROR in reademissions: latitude dimension not found in file: '//trim(em_name)
        error stop
      endif
      ! allocate emission variables
      allocate( em_field(nspec-1,nxem,nyem,2) )
      em_field(:,:,:,:)=0.
      allocate( em_res(nspec-1,nxem,nyem) )
      if (ipin.eq.2) then
        ! read residual emissions from end of previous run
        write(*,*) 'Reading residual emissions from previous run'
        call em_res_read
      else
        em_res(:,:,:)=0.
      endif
      allocate( mass_field(nspec,nxem,nyem) )
      mass_field(:,:,:)=0.
      ! calculate area for emissions grid
      allocate( em_area(nxem,nyem) )
      if (emis_unit(kk+1).eq.1) then
        ! emissions given per sqm
        do jy=1,nyem
          ylatp=latem(jy)+0.5*abs(latem(2)-latem(1))
          ylatm=latem(jy)-0.5*abs(latem(2)-latem(1))
          if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
            hzone=abs(latem(2)-latem(1))*r_earth*pi180
          else
            cosfactp=cos(ylatp*pi180)
            cosfactm=cos(ylatm*pi180)
            if (cosfactp.lt.cosfactm) then
              hzone=sqrt(1-cosfactp**2)-sqrt(1-cosfactm**2)
              hzone=hzone*r_earth
            else
              hzone=sqrt(1-cosfactm**2)-sqrt(1-cosfactp**2)
              hzone=hzone*r_earth
            endif
          endif
          em_area(:,jy)=pi180*r_earth*hzone*abs(lonem(2)-lonem(1))
        end do ! nyem
      else
        ! emissions given per gridcell
        em_area(:,:)=1.
      endif
    endif 
    
    ! read time dimension
    do nn=1,ndim
      call nf90_err( nf90_inquire_dimension( ncid, nn, dimname, len ) )
      if ((index(dimname,'time').ne.0).or.(index(dimname,'Time').ne.0).or. &
            (index(dimname,'TIME').ne.0).or.(index(dimname,'Date').ne.0).or. &
            (index(dimname,'date').ne.0).or.(index(dimname,'DATE').ne.0)) then
        ntem=len
        allocate( time(ntem), jdate(ntem) )
        call nf90_err( nf90_inq_varid(ncid,trim(dimname),varid) )
        call nf90_err( nf90_get_var(ncid,varid,time) )
        call nf90_err( nf90_get_att(ncid,varid,'units',unitinfo) )
        write(*,*) 'Time units: ',trim(unitinfo)
        if (index(unitinfo,'sec').ne.0) then
          ! seconds
          jdate=real(time-time(1),kind=dp)/3600._dp/24._dp+julstart
          indxt=minloc(abs(jdate-jtime),dim=1)
        else if (index(unitinfo,'hour').ne.0) then
          ! hours
          jdate=real(time-time(1),kind=dp)/24._dp+julstart
          indxt=minloc(abs(jdate-jtime),dim=1)
        else if (index(unitinfo,'day').ne.0) then
          ! days
          jdate=real(time-time(1),kind=dp)+julstart
          indxt=minloc(abs(jdate-jtime),dim=1)
        else if (index(unitinfo,'month').ne.0) then
          ! months
          call caldate(jtime,jjjjmmdd,ihmmss)
          mm=jjjjmmdd/10000
          indxt=minloc(abs(time-mm),dim=1)
        else
          write(*,*) 'ERROR in reademissions: unknown time units in file: '//trim(em_name)
          error stop
        endif
        deallocate( time, jdate )
        exit
      endif 
    end do ! ndim

    ! emission field
    allocate( emis(nxem,nyem) )
    write(*,*) 'Reading emissions for species '//trim(species(kk+1))
    call nf90_err( nf90_inq_varid(ncid,trim(emis_name(kk+1)),varid) )
    call nf90_err( nf90_inquire_variable(ncid,varid,nameout,xtype,ndim,dimids,natt) )
    sclfact=1.
    offset=0.
    if (natt.gt.0 ) then
      do nn=1,natt
        call nf90_err( nf90_inq_attname(ncid,varid,nn,attname) )
        if (index(attname,'add_offset').ne.0) then
          call nf90_err( nf90_get_att(ncid,varid,'add_offset',offset) )
        else if (index(attname,'scale_factor').ne.0) then
          call nf90_err( nf90_get_att(ncid,varid,'scale_factor',sclfact) )
        endif
      end do
    endif
    if (ndim.eq.2) then
      call nf90_err( nf90_get_var(ncid,varid,emis) )
    else if (ndim.eq.3) then
      write(*,*) 'reademissions: time index in file = ',indxt
      call nf90_err( nf90_get_var(ncid,varid,emis,start=(/1,1,indxt/),count=(/nxem,nyem,1/)) )
    endif
    emis=emis*sclfact+offset
    !! test
    print*, 'reademissions: sclfact, offset = ',sclfact, offset
    print*, 'reademissions: range(emis) = ',minval(emis), maxval(emis)
    print*, 'reademissions: range(em_area) = ',minval(em_area), maxval(em_area)
    print*, 'reademissions: emis_coeff, emis_unit = ',emis_coeff(kk+1), emis_unit(kk+1)
    print*, 'reademissions: emis total = ',sum(emis*em_area*emis_coeff(kk+1))*3600.*24.*365./1.e9
    em_field(kk,:,:,memid)=emis*em_area*emis_coeff(kk+1)

    ! close file
    call nf90_err( nf90_close(ncid) )

    deallocate(emis)

    return

  end subroutine reademissions

  subroutine em_res_write()

  !*****************************************************************************
  !                                                                            *
  !    This subroutine outputs the residual emissions to be used if            *
  !    picking up a run from where this one left off (ipin = 1)                *
  !                                                                            *
  !    Author: R. Thompson, Oct-2023                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: file_name
    integer :: nc_id, londim_id, latdim_id, specdim_id, nchardim_id
    integer :: spec_id, lon_id, lat_id, emres_id
    
    file_name=trim(path(2)(1:length(2)))//'emis_residual.nc' 
    
    call nf90_err( nf90_create(trim(file_name), nf90_clobber, ncid=nc_id) )

    ! define dimensions
    call nf90_err( nf90_def_dim(nc_id, 'species', nspec-1, specdim_id) )
    call nf90_err( nf90_def_dim(nc_id, 'longitude', nxem, londim_id) )
    call nf90_err( nf90_def_dim(nc_id, 'latitude', nyem, latdim_id) )
    call nf90_err( nf90_def_dim(nc_id, 'nchar', 18, nchardim_id) )

    ! define variables
    call nf90_err( nf90_def_var(nc_id, 'species', nf90_char, (/ nchardim_id, specdim_id /), spec_id) )
    call nf90_err( nf90_put_att(nc_id, spec_id, 'long_name', 'Species names') )
    call nf90_err( nf90_def_var(nc_id, 'longitude', nf90_float, (/ londim_id /), lon_id) )
    call nf90_err( nf90_put_att(nc_id, lon_id, 'units', 'degrees') )
    call nf90_err( nf90_def_var(nc_id, 'latitude', nf90_float, (/ latdim_id /), lat_id) )
    call nf90_err( nf90_put_att(nc_id, lat_id, 'units', 'degrees') )
    call nf90_err( nf90_def_var(nc_id, 'em_res', nf90_double, (/ specdim_id, londim_id, latdim_id /), emres_id) )
    call nf90_err( nf90_put_att(nc_id, emres_id, 'long_name', 'Emission residuals') )
    call nf90_err( nf90_put_att(nc_id, emres_id, 'units', 'kg') )

    call nf90_err( nf90_enddef(nc_id) )

    ! write variables
    call nf90_err( nf90_put_var(nc_id, spec_id, species(2:nspec)) )
    call nf90_err( nf90_put_var(nc_id, lon_id, lonem) )
    call nf90_err( nf90_put_var(nc_id, lat_id, latem) )
    call nf90_err( nf90_put_var(nc_id, emres_id, em_res) )

    call nf90_err( nf90_close(nc_id) )

  end subroutine em_res_write

  subroutine em_res_read()

  !*****************************************************************************
  !                                                                            *
  !    This subroutine reads the residual emissions from a previous run        *
  !    which have to be copied into the output directory of this run           *
  !                                                                            *
  !    Author: R. Thompson, Oct-2023                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: file_name
    integer :: nc_id, specdim_id, londim_id, latdim_id, emres_id
    integer :: xlen, ylen, splen
    logical :: lexist

    file_name=trim(path(2)(1:length(2)))//'emis_residual.nc'

    inquire(file=trim(file_name),exist=lexist)
    if (.not.lexist) then
      write(*,*) 'FLEXPART ERROR: cannot find file '//trim(file_name)
      error stop
    endif

    call nf90_err( nf90_open(file_name, nf90_nowrite, nc_id) )

    ! read dimensions
    call nf90_err( nf90_inq_dimid(nc_id, 'longitude', londim_id) )
    call nf90_err( nf90_inquire_dimension(nc_id, londim_id, len=xlen) )
    call nf90_err( nf90_inq_dimid(nc_id, 'latitude', latdim_id) )
    call nf90_err( nf90_inquire_dimension(nc_id, latdim_id, len=ylen) )
    call nf90_err( nf90_inq_dimid(nc_id, 'species', specdim_id) )
    call nf90_err( nf90_inquire_dimension(nc_id, specdim_id, len=splen) )

    ! check dimensions match input emissions
    if ((xlen.ne.nxem).or.(ylen.ne.nyem).or.(splen.ne.(nspec-1))) then
      write(*,*) 'FLEXPART ERROR: emis_residual dimensions do not match input emissions'
      error stop
    endif

    ! read em_res
    call nf90_err( nf90_inq_varid(nc_id, 'em_res', emres_id) )
    call nf90_err( nf90_get_var(nc_id, emres_id, em_res) )

    call nf90_err( nf90_close(nc_id) )

  end subroutine em_res_read


end module emissions_mod
