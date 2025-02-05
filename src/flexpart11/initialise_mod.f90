! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !   L. Bakels 2021: This module contains subroutines related to the          *
  !                   initialisation of the particles                          *
  !                                                                            *
  !*****************************************************************************

module initialise_mod

  use com_mod
  use par_mod
  use date_mod
  use particle_mod
  use windfields_mod
  use random_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
#ifdef USE_NCF
  use netcdf_output_mod
#endif

  implicit none

  !**********************************************************
  ! Variables used for domain-filling trajectory calculations
  !**********************************************************

  integer ::    &
    nx_we(2),   & ! x indices of western and eastern boundary of domain-filling.
    ny_sn(2),   & ! y indices of southern and northern boundary of domain-filling.
    numcolumn     ! Maximum number of particles to be released within a single column.
  integer,allocatable,dimension(:,:) :: &
    numcolumn_we,   & ! Number of particles to be released within one column
                      ! at the western and eastern boundary surfaces.
    numcolumn_sn      ! Same as numcolumn_we, but for southern and northern domain boundary.
  real,allocatable,dimension(:,:,:) ::  &
    zcolumn_we,     & ! Altitudes where particles are to be released
                      ! at the western and eastern boundary surfaces.
    zcolumn_sn,     & ! Same as zcolumn_we, but for southern and northern domain boundary.
    acc_mass_we,    & ! Mass that has accumulated at the western and eastern boundary;
                      ! if it exceeds xmassperparticle, a particle is released and
                      ! acc_mass_we is reduced accordingly.
    acc_mass_sn       ! Same as acc_mass_we, but for southern and northern domain boundary
  real ::           &
    xmassperparticle  ! Air mass per particle in the domain-filling traj. option.

contains

subroutine alloc_domainfill
  implicit none
  allocate(numcolumn_we(2,0:nymax-1),numcolumn_sn(2,0:nxmax-1))
  allocate(zcolumn_we(2,0:nymax-1,maxcolumn),zcolumn_sn(2,0:nxmax-1,maxcolumn), &
    acc_mass_we(2,0:nymax-1,maxcolumn),acc_mass_sn(2,0:nxmax-1,maxcolumn))              
end subroutine alloc_domainfill

subroutine dealloc_domainfill
  if (mdomainfill.lt.1) return
  deallocate(numcolumn_we,numcolumn_sn,zcolumn_sn,zcolumn_we,acc_mass_sn, &
    acc_mass_we)
end subroutine dealloc_domainfill

subroutine releaseparticles(itime)
  !                              o
  !*****************************************************************************
  !                                                                            *
  !     This subroutine releases particles from the release locations.         *
  !                                                                            *
  !     It searches for a "vacant" storage space and assigns all particle      *
  !     information to that space. A space is vacant either when no particle   *
  !     is yet assigned to it, or when it's particle is expired and, thus,     *
  !     the storage space is made available to a new particle.                 *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     29 June 2002                                                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            current time                                          *
  ! ireleasestart, ireleaseend          start and end times of all releases    *
  ! npart(maxpoint)      number of particles to be released in total           *
  ! numrel               number of particles to be released during this time   *
  !                      step                                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use xmass_mod
  use output_mod
  use interpol_mod

  implicit none

  !real xaux,yaux,zaux,ran1,rfraction,xmasssave(maxpoint)
  real :: xaux,yaux,zaux,rfraction
  real :: xlonav,timecorrect(maxspec)
  real :: average_timecorrect
  integer :: itime,numrel,i,j,k,ipart,minpart
  integer :: istart,iend,totpart,iterm_index
  integer :: nweeks,ndayofweek,nhour,jjjjmmdd,ihmmss,mm
  real(kind=dp) :: julmonday,jul,jullocal,juldiff

  integer :: idummy = -7
  !save idummy,xmasssave
  !data idummy/-7/,xmasssave/maxpoint*0./

  real :: eps
  eps=nxmax/3.e5


  ! Determine the actual date and time in Greenwich 
  ! (i.e., UTC + correction for daylight savings time)
  !***************************************************

  julmonday=juldate(19000101,0)          ! this is a Monday
  jul=bdate+real(itime,kind=dp)/86400._dp    ! this is the current day
  call caldate(jul,jjjjmmdd,ihmmss)
  mm=(jjjjmmdd-10000*(jjjjmmdd/10000))/100
  if ((mm.ge.4).and.(mm.le.9)) jul=jul+1._dp/24._dp   ! daylight savings time


  ! For every release point, check whether we are in the release time interval
  !***************************************************************************
  ! First allocate all particles that are going to be in the simulation
  ! If ipin==0,1, and ipout==0, then dead particles can be overwritten to save memory
  if (ipin.gt.1 .or. ipout.ne.0 .or. ispeed.eq.1) then
    if (itime.eq.0) then
      totpart=0
      do i=1,numpoint
        totpart = totpart+npart(i)
      end do
      call alloc_particles(totpart)
    else if (itime.eq.itime_init) then !From restart point only allocate particles that are yet to be born
      totpart=0
      do i=1,numpoint
        totpart = totpart+npart(i)
      end do
      if (totpart.gt.count%allocated) call alloc_particles(totpart-count%allocated)
    end if 
  endif

  call get_totalpart_num(istart)
  if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0) call rewrite_iterm()
  minpart=1
  do i=1,numpoint
    if ((itime.ge.ireleasestart(i)).and. &! are we within release interval?
         (itime.le.ireleaseend(i))) then

  ! Determine the local day and time
  !*********************************

      xlonav=xlon0+(xpoint2(i)+xpoint1(i))*0.5*dx  ! longitude needed to determine local time
      if (xlonav.lt.-180.) xlonav=xlonav+360.
      if (xlonav.gt.180.) xlonav=xlonav-360.
      jullocal=jul+real(xlonav,kind=dp)/360._dp   ! correct approximately for time zone to obtain local time

      juldiff=jullocal-julmonday
      nweeks=int(juldiff/7._dp)
      juldiff=juldiff-real(nweeks,kind=dp)*7._dp
      ndayofweek=int(juldiff)+1              ! this is the current day of week, starting with Monday
      nhour=nint((juldiff-real(ndayofweek-1,kind=dp))*24._dp)    ! this is the current hour
      if (nhour.eq.0) then
        nhour=24
        ndayofweek=ndayofweek-1
        if (ndayofweek.eq.0) ndayofweek=7
      endif

  ! Calculate a species- and time-dependent correction factor, distinguishing between
  ! area (those with release starting at surface) and point (release starting above surface) sources
  ! Also, calculate an average time correction factor (species independent)
  !*****************************************************************************
      average_timecorrect=0.
      do k=1,nspec
        if(abs(xpoint2(i)-xpoint1(i)).lt.1.E-4.and.abs(ypoint2(i)-ypoint1(i)).lt.1.E-4) then
  !        if (zpoint1(i).gt.0.5) then      ! point source
          timecorrect(k)=point_hour(k,nhour)*point_dow(k,ndayofweek)
        else                             ! area source
          timecorrect(k)=area_hour(k,nhour)*area_dow(k,ndayofweek)
        endif
        average_timecorrect=average_timecorrect+timecorrect(k)
      end do
      average_timecorrect=average_timecorrect/real(nspec)

  ! Determine number of particles to be released this time; at start and at end of release,
  ! only half the particles are released
  !*****************************************************************************

      if (ireleasestart(i).ne.ireleaseend(i)) then
        rfraction=abs(real(npart(i))*real(lsynctime)/ &
             real(ireleaseend(i)-ireleasestart(i)))
        if ((itime.eq.ireleasestart(i)).or. &
             (itime.eq.ireleaseend(i))) rfraction=rfraction*0.5

  ! Take the species-average time correction factor in order to scale the
  ! number of particles released this time
  !**********************************************************************
        rfraction=rfraction*average_timecorrect

        rfraction=rfraction+xmasssave(i)  ! number to be released at this time
        numrel=int(rfraction)
        xmasssave(i)=rfraction-real(numrel)
      else
        numrel=npart(i)
      endif

      xaux=xpoint2(i)-xpoint1(i)
      yaux=ypoint2(i)-ypoint1(i)
      zaux=zpoint2(i)-zpoint1(i)

      if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0 ) then
        call rewrite_iterm()
        totpart = numrel-count%iterm_max
        if (totpart.gt.0) call alloc_particles(totpart)
        call rewrite_iterm()
        iterm_index=1
      endif
      do j=1,numrel             ! loop over particles to be released this time
        call get_newpart_index(ipart,iterm_index)
        call spawn_particle(itime, ipart)

  ! Particle coordinates are determined by using a random position within the release volume
  !*****************************************************************************

  ! Determine horizontal particle position
  !***************************************
        call set_xlon(ipart,real(xpoint1(i)+ran1(idummy,0)*xaux,kind=dp))
        if (xglobal) then
          if (part(ipart)%xlon.gt.real(nxmin1,kind=dp)) &
            call set_xlon(ipart,-real(nxmin1,kind=dp))
          if (part(ipart)%xlon.lt.0.) &
            call set_xlon(ipart,real(nxmin1,kind=dp))
        endif
        call set_ylat(ipart,real(ypoint1(i)+ran1(idummy,0)*yaux,kind=dp))

  ! Assign mass to particle: Total mass divided by total number of particles.
  ! Time variation has partly been taken into account already by a species-average
  ! correction factor, by which the number of particles released this time has been
  ! scaled. Adjust the mass per particle by the species-dependent time correction factor
  ! divided by the species-average one
  ! for the scavenging calculation the mass needs to be multiplied with rho of the particle layer and
  ! divided by the sum of rho of all particles.
  !*****************************************************************************
        do k=1,nspec
          mass(ipart,k)=xmass(i,k)/real(npart(i)) &
                *timecorrect(k)/average_timecorrect
          mass_init(ipart,k)=mass(ipart,k)
        end do
  ! Assign certain properties to particle
  !**************************************
        part(ipart)%nclass=min(int(ran1(idummy,0)*real(nclassunc))+1, &
             nclassunc)
        numparticlecount=numparticlecount+1
        if (mquasilag.eq.0) then
          part(ipart)%npoint=i
        else
          part(ipart)%npoint=numparticlecount
        endif
        part(ipart)%idt=mintime               ! first time step

        ! Determine vertical particle position
        !*************************************
        call set_z(ipart,zpoint1(i)+ran1(idummy,0)*zaux)
        ! Interpolation of topography and density
        !****************************************

        ! Transform the verticle particle position from pressure or sea level to above ground
        ! if necessary
        !************************************************************************************
        call kindz_to_z(ipart)
#ifdef ETA
        call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat,part(ipart)%z,part(ipart)%zeta)
        part(ipart)%etaupdate = .true. ! The z(meter) coordinate is up to date
        part(ipart)%meterupdate = .true.
#endif
        
        call init_mass_conversion(ipart,i)

        call get_totalpart_num(numpart)

      end do  ! numrel 
    endif ! releasepoint
  end do ! numpoint

  if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0 ) call rewrite_iterm()

  call get_totalpart_num(iend)

  ! NetCDF only: write initial positions of new particles
! #ifdef USE_NCF
!   if ((iend-istart.gt.0).and.(ipout.ge.1)) then 
!     call wrt_part_initialpos(itime,istart,iend)
!     call output_particles(itime,.true.)
!   endif
! #endif
  return

! 996   continue
!   write(*,*) '#####################################################'
!   write(*,*) '#### FLEXPART MODEL SUBROUTINE RELEASEPARTICLES: ####'
!   write(*,*) '####                                             ####'
!   write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES REQUIRED  ####'
!   write(*,*) '#### EXCEEDS THE MAXIMUM ALLOWED NUMBER. REDUCE  ####'
!   write(*,*) '#### EITHER NUMBER OF PARTICLES PER RELEASE POINT####'
!   write(*,*) '#### OR REDUCE NUMBER OF RELEASE POINTS.         ####'
!   write(*,*) '#####################################################'
!   stop

end subroutine releaseparticles

subroutine kindz_to_z(ipart)
  use point_mod
  use xmass_mod
  use output_mod
  use interpol_mod

  implicit none

  integer,intent(in) :: ipart
  integer :: kz
  real :: dp1,dp2,press,presspart,pressold,topo,r,t
  real,parameter :: eps2=1.e-6


  ! Determine the nest we are in
  !*****************************
  call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)

  ! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
  !*****************************************************************************
  call find_grid_indices(real(part(ipart)%xlon),real(part(ipart)%ylat))
  call find_grid_distances(real(part(ipart)%xlon),real(part(ipart)%ylat))


  if (kindz(part(ipart)%npoint).eq.1) return ! Nothing needs to happen

  ! If starting height is in pressure coordinates, retrieve pressure profile and 
  ! convert zpart1 to meters
  !*****************************************************************************

  if (kindz(part(ipart)%npoint).eq.3) then
    presspart=real(part(ipart)%z)
    do kz=1,nz
      if (ngrid.gt.0) then
        r=p1*rhon(ix ,jy ,kz,2,ngrid) &
             +p2*rhon(ixp,jy ,kz,2,ngrid) &
             +p3*rhon(ix ,jyp,kz,2,ngrid) &
             +p4*rhon(ixp,jyp,kz,2,ngrid)
        t=p1*ttn(ix ,jy ,kz,2,ngrid) &
             +p2*ttn(ixp,jy ,kz,2,ngrid) &
             +p3*ttn(ix ,jyp,kz,2,ngrid) &
             +p4*ttn(ixp,jyp,kz,2,ngrid)
      else
        r=p1*rho(ix ,jy ,kz,2) &
             +p2*rho(ixp,jy ,kz,2) &
             +p3*rho(ix ,jyp,kz,2) &
             +p4*rho(ixp,jyp,kz,2)
        t=p1*tt(ix ,jy ,kz,2) &
             +p2*tt(ixp,jy ,kz,2) &
             +p3*tt(ix ,jyp,kz,2) &
             +p4*tt(ixp,jyp,kz,2)
      endif
      press=r*r_air*t/100.
      if (kz.eq.1) pressold=press

      if (press.lt.presspart) then
        if (kz.eq.1) then
          call set_z(ipart,height(1)*0.5)
        else
          dp1=pressold-presspart
          dp2=presspart-press
          call set_z(ipart,(height(kz-1)*dp2+height(kz)*dp1) &
               /(dp1+dp2))
        endif
        exit
      endif
      pressold=press
    end do

  ! If release positions are given in meters above sea level, subtract the
  ! topography from the starting height
  !***********************************************************************

  else if (kindz(part(ipart)%npoint).eq.2) then
    if (ngrid.gt.0) then
      topo=p1*oron(ix ,jy ,ngrid) &
           + p2*oron(ixp,jy ,ngrid) &
           + p3*oron(ix ,jyp,ngrid) &
           + p4*oron(ixp,jyp,ngrid)
    else
      topo=p1*oro(ix ,jy) &
           + p2*oro(ixp,jy) &
           + p3*oro(ix ,jyp) &
           + p4*oro(ixp,jyp)
    endif
    call update_z(ipart,-topo)
  endif
  if (part(ipart)%z.lt.eps2) call set_z(ipart,eps2)   ! Minimum starting height is eps2
  if (part(ipart)%z.gt.height(nz)-0.5) &
    call set_z(ipart,height(nz)-0.5) ! Maximum starting height is uppermost level - 0.5 meters

  if (ipin.eq.3 .or. ipin.eq.4) then
    if (part(ipart)%z.gt.zpoint2(part(ipart)%npoint)) &
      zpoint2(part(ipart)%npoint)=real(part(ipart)%z)
    if (part(ipart)%z.lt.zpoint1(part(ipart)%npoint)) &
      zpoint1(part(ipart)%npoint)=real(part(ipart)%z)
  endif    

end subroutine kindz_to_z

subroutine init_mass_conversion(ipart,ipoint)
  ! For special simulations, multiply particle concentration air density;
  ! Simply take the 2nd field in memory to do this (accurate enough)
  !***********************************************************************
  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af    NOTE that in backward simulations the release of particles takes place at the
  !Af         receptor and the sampling at the source.
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !            3 = wet deposition in outputfield
  !            4 = dry deposition in outputfield

  !Af switches for the releasefile:
  !Af IND_REL =  1 : xmass * rho
  !Af IND_REL =  0 : xmass * 1

  !Af ind_rel is defined in readcommand.f
  use point_mod
  use xmass_mod
  use interpol_mod

  implicit none

  integer,intent(in) :: ipart,ipoint
  integer :: n
  real :: rhoaux(2),rhoout
  real :: dz1,dz2


  if ((ind_rel .eq. 1).or.(ind_rel .eq. 3).or.(ind_rel .eq. 4)) then

    ! Interpolate the air density, horizontal values are computed in kindz_to_z
    !**************************************************************************
    call find_z_level_meters(real(part(ipart)%z))
    call find_vert_vars_lin(height,real(part(ipart)%z),indz,dz1,dz2,lbounds)

    if (ngrid.gt.0) then
      do n=1,2
        rhoaux(n)=p1*rhon(ix ,jy ,indz+n-1,2,ngrid) &
             +p2*rhon(ixp,jy ,indz+n-1,2,ngrid) &
             +p3*rhon(ix ,jyp,indz+n-1,2,ngrid) &
             +p4*rhon(ixp,jyp,indz+n-1,2,ngrid)
      end do
    else
      do n=1,2
        rhoaux(n)=p1*rho(ix ,jy ,indz+n-1,2) &
             +p2*rho(ixp,jy ,indz+n-1,2) &
             +p3*rho(ix ,jyp,indz+n-1,2) &
             +p4*rho(ixp,jyp,indz+n-1,2)
      end do
    endif
    rhoout=dz2*rhoaux(1)+dz1*rhoaux(2)
    rho_rel(ipoint)=rhoout


    ! Multiply "mass" (i.e., mass mixing ratio in forward runs) with density
    !***********************************************************************
    mass(ipart,:)=mass(ipart,:)*rhoout
    mass_init(ipart,:)=mass(ipart,:)
  endif
end subroutine init_mass_conversion

subroutine readpartpositions

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the particle dump file and reads all the particle     *
  !   positions from a previous run to initialize the current run.             *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 March 2000                                                          *
  !                                                                            *
  !  Changes                                                                   *
  !     2022, L. Bakels: NetCDF option for reading particle information        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************
  implicit none

  integer :: ibdatein,ibtimein,nspecin,itimein,numpointin,i,j,lix,ios
  integer :: id1,id2,it1,it2
  real :: xlonin,ylatin,topo,hmixi,pvi,qvi,rhoi,tri,tti
  character :: specin*7
  real(kind=dp) :: julin,julpartin

  integer :: idummy = -8

  numparticlecount=0

  ! Open header file of dumped particle data
  !*****************************************
#ifdef USE_NCF
    call readpartpositions_netcdf(ibtime,ibdate)
    call get_totalpart_num(numpart)
    numparticlecount=numpart
    return
#endif

  open(unitpartin,file=path(2)(1:length(2))//'header', &
       form='unformatted',err=998)

  read(unitpartin) ibdatein,ibtimein
  read(unitpartin)
  read(unitpartin)

  read(unitpartin)
  read(unitpartin)
  read(unitpartin) nspecin
  nspecin=nspecin/3
  if ((ldirect.eq.1).and.(nspec.ne.nspecin)) then
    write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
    write(*,*) ' #### THE NUMBER OF SPECIES TO BE READ IN DOES #### '
    write(*,*) ' #### NOT AGREE WITH CURRENT SETTINGS!         #### '
    stop
  end if

  do i=1,nspecin
    read(unitpartin)
    read(unitpartin)
    read(unitpartin) j,specin
    if ((ldirect.eq.1).and.(species(i)(1:7).ne.specin)) then
      write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
      write(*,*) ' #### SPECIES NAMES TO BE READ IN DO NOT       #### '
      write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
      stop
    end if
  end do

  read(unitpartin) numpointin
  if (numpointin.ne.numpoint) then
    write(*,*) ' #### FLEXPART MODEL WARNING IN READPARTPOSITIONS#### '
    write(*,*) ' #### NUMBER OF RELEASE LOCATIONS DOES NOT     #### '
    write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
  end if 
  do i=1,numpointin
    read(unitpartin)
    read(unitpartin)
    read(unitpartin)
    read(unitpartin)
    do j=1,nspec
      read(unitpartin)
      read(unitpartin)
      read(unitpartin)
    end do
  end do
  read(unitpartin)
  read(unitpartin)

  do lix=0,numxgrid-1
    read(unitpartin)
  end do


  ! Open data file of dumped particle data
  !***************************************

  close(unitpartin)
  open(unitpartin,file=path(2)(1:length(2))//'partposit_end', &
       form='unformatted',err=998)
  

  do 
    read(unitpartin,iostat=ios) itimein
    if (ios.lt.0) exit
    i=0
    do
      i=i+1
      read(unitpartin) part(i)%npoint,xlonin,ylatin,part(i)%z,part(i)%tstart, &
           topo,pvi,qvi,rhoi,hmixi,tri,tti,(mass(i,j),j=1,nspec)
      ! For switching coordinates: this happens in timemanager.f90 after the first fields are read
      if (xlonin.eq.-9999.9) exit
      call set_xlon(i,real((xlonin-xlon0)/dx,kind=dp))
      call set_ylat(i,real((ylatin-ylat0)/dy,kind=dp))
      numparticlecount=max(numparticlecount,part(i)%npoint)
    end do
  end do

  numpart=i-1

  close(unitpartin)

  julin=juldate(ibdatein,ibtimein)+real(itimein,kind=dp)/86400._dp
  if (abs(julin-bdate).gt.1.e-5) then
    write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
    write(*,*) ' #### ENDING TIME OF PREVIOUS MODEL RUN DOES   #### '
    write(*,*) ' #### NOT AGREE WITH STARTING TIME OF THIS RUN.#### '
    call caldate(julin,id1,it1)
    call caldate(bdate,id2,it2)
    write(*,*) 'julin: ',julin,id1,it1
    write(*,*) 'bdate: ',bdate,id2,it2
    stop
  end if
  do i=1,numpart
    julpartin=juldate(ibdatein,ibtimein)+ &
         real(part(i)%tstart,kind=dp)/86400._dp
    part(i)%nclass=min(int(ran1(idummy,0)*real(nclassunc))+1, &
         nclassunc)
    part(i)%idt=mintime
    part(i)%tstart=nint((julpartin-bdate)*86400.)
  end do

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'partposit'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine readpartpositions

subroutine init_particle(itime,ipart,ithread)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of trajectories utilizing a zero-acceleration scheme. The time*
  !  step is determined by the Courant-Friedrichs-Lewy (CFL) criterion. This   *
  !  means that the time step must be so small that the displacement within    *
  !  this time step is smaller than 1 grid distance. Additionally, a temporal  *
  !  CFL criterion is introduced: the time step must be smaller than the time  *
  !  interval of the wind fields used for interpolation.                       *
  !  For random walk simulations, these are the only time step criteria.       *
  !  For the other options, the time step is also limited by the Lagrangian    *
  !  time scale.                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Literature:                                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          current temporal position                               *
  ! ldt [s]            Suggested time step for next integration                *
  ! ladvance [s]       Total integration time period                           *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! usig,vsig,wsig     uncertainties of wind velocities due to interpolation   *
  ! xt,yt,zt           Next time step's spatial position of trajectory         *
  !                                                                            *
  !                                                                            *
  ! Constants:                                                                 *
  ! cfl                factor, by which the time step has to be smaller than   *
  !                    the spatial CFL-criterion                               *
  ! cflt               factor, by which the time step has to be smaller than   *
  !                    the temporal CFL-criterion                              *
  !                                                                            *
  !*****************************************************************************

  use turbulence_mod
  use random_mod, only: ran3
  use omp_lib
  use interpol_mod
  use cbl_mod

  implicit none

  integer,intent(in) ::  &
    ithread,             & ! OMP thread starting at 0
    itime,               &
    ipart
  integer :: nrand
  real :: dummy1,dummy2
  real :: xt,yt,zt,zteta

  ! ! Initialise scavenging for backward runs
  ! !****************************************
  ! if (DRYBKDEP.or.WETBKDEP) then ! if there is no scavenging in wetdepo it will be set to 0
  !   do k=1,nspec
  !     xscav_frac1(ipart,k)=-1.
  !   end do
  ! endif


  part(ipart)%icbt=1           ! initialize particle to no "reflection"

  nrand=int(ran3(iseed1(ithread),ithread)*real(maxrand-1))+1

  xt = real(part(ipart)%xlon)
  yt = real(part(ipart)%ylat)
  zt = real(part(ipart)%z)
#ifdef ETA
  zteta = real(part(ipart)%zeta)
#else
  zteta = 0.
#endif

  !******************************
  ! 2. Interpolate necessary data
  !******************************

  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(xt,yt)
  ! Compute maximum mixing height around particle position
  !*******************************************************
  call find_grid_indices(xt,yt)

  ! If part_ic.nc is used, convert particle mass to mixing mass in forward,
  ! or internal mass conversion for wet and dry backward depostion.
  !************************************************************************
  if ((iout.eq.4).or.(iout.eq.5)) call init_mass_conversion(ipart,part(ipart)%npoint)
  
  h=max(hmix(ix ,jy,1,memind(1)), &
       hmix(ixp,jy ,1,memind(1)), &
       hmix(ix ,jyp,1,memind(1)), &
       hmix(ixp,jyp,1,memind(1)), &
       hmix(ix ,jy ,1,memind(2)), &
       hmix(ixp,jy ,1,memind(2)), &
       hmix(ix ,jyp,1,memind(2)), &
       hmix(ixp,jyp,1,memind(2)))

  zeta=zt/h


  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************

  if (zeta.le.1.) then

    call interpol_pbl(itime,xt,yt,zt,zteta,ithread+1)

  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************
    call interpol_pbl_short(zt,dummy1,dummy2,ithread+1)

  ! Compute the turbulent disturbances

  ! Determine the sigmas and the timescales
  !****************************************

    if (turbswitch) then
      call hanna(zt)
    else
      call hanna1(zt)
    endif


  ! Determine the new diffusivity velocities
  !*****************************************

    if (nrand+2.gt.maxrand) nrand=1
    part(ipart)%turbvel%u=rannumb(nrand)*sigu
    part(ipart)%turbvel%v=rannumb(nrand+1)*sigv
    part(ipart)%turbvel%w=rannumb(nrand+2)
    if (.not.turbswitch) then     ! modified by mc
      part(ipart)%turbvel%w=part(ipart)%turbvel%w*sigw
    else if (cblflag.eq.1) then   ! modified by mc
      if(-h/ol.gt.5) then
  !if (ol.lt.0.) then
  !if (ol.gt.0.) then !by mc : only for test correct is lt.0
        call init_cbl_vel( &
          iseed1(ithread),zt,wst,h,sigw,part(ipart)%turbvel%w,ol,ithread)
      else
        part(ipart)%turbvel%w=part(ipart)%turbvel%w*sigw
      end if
    end if


  ! Determine time step for next integration
  !*****************************************

    if (turbswitch) then
      part(ipart)%idt = int( &
        min(tlw, &
          h/max(2.*abs(part(ipart)%turbvel%w*sigw),1.e-5), &
          0.5/abs(dsigwdz), &
          600.) &
        *ctl)
    else
      part(ipart)%idt = int( &
        min(tlw, &
          h/max(2.*abs(part(ipart)%turbvel%w),1.e-5), &
          600.) &
        *ctl)
    endif
    part(ipart)%idt=max(part(ipart)%idt,mintime)

    ! call interpol_average()
    ! usig=(usigprof(indzp)+usigprof(indz))/2.
    ! vsig=(vsigprof(indzp)+vsigprof(indz))/2.
    ! wsig=(wsigprof(indzp)+wsigprof(indz))/2.

    ! wsigeta=(wsigprofeta(indzpeta)+wsigprofeta(indzeta))/2.

  else



  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************


  ! Interpolate the wind
  !*********************

    call interpol_wind(itime,xt,yt,zt,zteta)


  ! Compute everything for above the PBL

  ! Assume constant turbulent perturbations
  !****************************************

    part(ipart)%idt=abs(lsynctime)

    if (nrand+1.gt.maxrand) nrand=1
    part(ipart)%turbvel%u=rannumb(nrand)*0.3
    part(ipart)%turbvel%v=rannumb(nrand+1)*0.3
    nrand=nrand+2
    part(ipart)%turbvel%w=0.
    sigw=0.

  endif

  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************


  ! It is assumed that the average interpolation error is 1/2 sigma
  ! of the surrounding points, autocorrelation time constant is
  ! 1/2 of time interval between wind fields
  !****************************************************************
  if (lmesoscale_turb) then
    call interpol_mesoscale(xt,yt,zt,zteta)
    if (nrand+2.gt.maxrand) nrand=1
    part(ipart)%mesovel%u=rannumb(nrand)*usig
    part(ipart)%mesovel%v=rannumb(nrand+1)*vsig
#ifdef ETA
    part(ipart)%mesovel%w=rannumb(nrand+2)*wsigeta
#else
    part(ipart)%mesovel%w=rannumb(nrand+2)*wsig
#endif
  endif
end subroutine init_particle

subroutine init_domainfill
  !
  !*****************************************************************************
  !                                                                            *
  ! Initializes particles equally distributed over the first release location  *
  ! specified in file RELEASES. This box is assumed to be the domain for doing *
  ! domain-filling trajectory calculations.                                    *
  ! All particles carry the same amount of mass which alltogether comprises the*
  ! mass of air within the box.                                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     15 October 2002                                                        *
  !                                                                            *
  !  Changes                                                                   *
  !     2022, L. Bakels: OpenMP parallelisation                                *
  !     2023, L. Bakels: smooth vertical particle distribution instead of      *
  !                      distributing particles on fixed vertical layers       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! numparticlecount    consecutively counts the number of particles released  *
  ! nx_we(2)       grid indices for western and eastern boundary of domain-    *
  !                filling trajectory calculations                             *
  ! ny_sn(2)       grid indices for southern and northern boundary of domain-  *
  !                filling trajectory calculations                             *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use particle_mod

  implicit none

  integer :: j,kz,lix,ljy,ncolumn,numparttot,stat,iterminate
  real :: ylat,ylatp,ylatm,hzone
  real :: cosfactm,cosfactp,deltacol,dz1,dz2,dz,pnew,pnew_temp,fractus
  real,parameter :: pih=pi/180.
  real :: colmasstotal,zposition
  real,allocatable,dimension(:) :: pp

  integer :: ixm,ixp,jym,jyp,indzm,indzh,indzp,i,jj,ii
  ! integer :: alive_tmp,allocated_tmp,spawned_tmp,terminated_tmp
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2)
  integer :: idummy = -11

  real :: height_tmp

  logical :: deall

  real,allocatable,dimension(:)   :: gridarea !
  real,allocatable,dimension(:,:) :: colmass !

  ! Determine the release region (only full grid cells), over which particles
  ! shall be initialized
  ! Use 2 fields for west/east and south/north boundary
  !**************************************************************************
  call alloc_domainfill

  nx_we(1)=max(int(xpoint1(1)),0)
  nx_we(2)=min((ceiling(xpoint2(1))),nxmin1)
  ny_sn(1)=max(int(ypoint1(1)),0)
  ny_sn(2)=min((ceiling(ypoint2(1))),nymin1)

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
  write(*,*) 'Global domain: ', gdomainfill

  ! Exit here if resuming a run from particle dump
  !***********************************************
  if (gdomainfill.and.ipin.ne.0) return

  ! Allocate grid and column mass
  !*******************************
  allocate(gridarea(0:nymax-1),colmass(0:nxmax-1,0:nymax-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate gridarea or colmass'

  ! Do not release particles twice (i.e., not at both in the leftmost and rightmost
  ! grid cell) for a global domain
  !*****************************************************************************
  if (xglobal) nx_we(2)=min(nx_we(2),nx-2)


  ! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
  ! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
  !************************************************************
  ! First for the south pole
  
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

  ! Do the same for the north pole

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



  ! Allocate memory for storing the particles
  !******************************************
  call alloc_particles(int(npart(1)*1.1)) ! A bit more to avoid single part alloc

  ! Initialise total particle number
  numparttot=0
  ! Initialise max column number
  numcolumn=0

  ! Initialise the sum over the total mass of the atmosphere
  colmasstotal=0.

!$OMP PARALLEL PRIVATE(ljy,ylat,ylatp,ylatm,hzone,cosfactp,cosfactm,pp,lix) &
!$OMP REDUCTION(+:colmasstotal)
  
  allocate( pp(nzmax),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate pp inside of OMP loop'

!$OMP DO
  do ljy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    ylat=ylat0+real(ljy)*dy
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
    gridarea(ljy)=2.*pi*r_earth*hzone*dx/360.
  end do
!$OMP END DO
!$OMP BARRIER

  ! Calculate total mass of each grid column and of the whole atmosphere
  !*********************************************************************
!$OMP DO 
  do ljy=ny_sn(1),ny_sn(2)          ! loop about latitudes
    do lix=nx_we(1),nx_we(2)      ! loop about longitudes
      pp(1)=prs(lix,ljy,1,1) !rho(lix,ljy,1,1)*r_air*tt(lix,ljy,1,1)
      pp(nz)=prs(lix,ljy,nz,1) !rho(lix,ljy,nz,1)*r_air*tt(lix,ljy,nz,1)
      colmass(lix,ljy)=(pp(1)-pp(nz))/ga*gridarea(ljy)
      colmasstotal=colmasstotal+colmass(lix,ljy)
    end do
  end do
!$OMP END DO
  deallocate(pp)
!$OMP END PARALLEL

  write(*,*) 'Atm. mass: ',colmasstotal

  allocate( pp(nzmax),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate pp'

  if (ipin.eq.0) numpart=0

  ! Determine the particle positions
  !*********************************
  iterminate=0
  do ljy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    ylat=ylat0+real(ljy)*dy
    do lix=nx_we(1),nx_we(2)      ! loop about longitudes
      ncolumn=nint(0.999*real(npart(1))*colmass(lix,ljy)/colmasstotal)
      if (ncolumn.eq.0) cycle
      if (ncolumn.gt.numcolumn) numcolumn=ncolumn

  ! Calculate pressure at the altitudes of model surfaces, using the air density
  ! information, which is stored as a 3-d field
  !*****************************************************************************

      pp(:)=prs(lix,ljy,:,1)!rho(lix,ljy,kz,1)*r_air*tt(lix,ljy,kz,1)


      deltacol=(pp(1)-pp(nz))/real(ncolumn)
      pnew=pp(1)+deltacol*0.5
      jj=0
      do j=1,ncolumn ! looping over the number of particles within the column

        ! For columns with many particles (i.e. around the equator), distribute
        ! the particles equally (1 on a random position within the deltacol range), 
        ! for columns with few particles (i.e. around the poles), 
        ! distribute the particles randomly
        !***********************************************************************

        if ((ncolumn.gt.20).and.(ncolumn-j.gt.20)) then
          pnew_temp=pnew-ran1(idummy,0)*deltacol
          pnew=pnew-deltacol
        else if ((ncolumn.gt.20).and.(ncolumn-j.le.20)) then
          ! When only few particles are left, distribute them randomly above pnew
          pnew_temp=pnew-ran1(idummy,0)*(pnew-pp(nz))
        else
          pnew_temp=pp(1)-ran1(idummy,0)*(pp(1)-pp(nz))
        endif

        do kz=1,nz-1
          if ((pp(kz).ge.pnew_temp).and.(pp(kz+1).lt.pnew_temp)) then
            dz1=log(pnew_temp)-log(pp(kz))
            dz=1./log(pp(kz+1)/pp(kz))

            ! Assign particle position
            !*************************
            ! Do the following steps only if particles are not read in 
            ! from previous model run
            !**********************************************************
            if (ipin.eq.0) then
              ! First spawn the particle into existence
              !****************************************
              jj=jj+1
              !THIS WILL CAUSE PROBLEMS WITH OMP! because of dynamical allocation
              call spawn_particle(0,numpart+jj)
              call set_xlon(numpart+jj,real(real(lix)-0.5+ran1(idummy,0),kind=dp))
              if (lix.eq.0) call set_xlon(numpart+jj,real(ran1(idummy,0),kind=dp))
              if (lix.eq.nxmin1) &
                call set_xlon(numpart+jj,real(real(nxmin1)-ran1(idummy,0),kind=dp))
              call set_ylat(numpart+jj,real(real(ljy)-0.5+ran1(idummy,0),kind=dp))
              ! Logarithmic distribution of particles along pressure levels:
              ! hx=h1+(h2-h1)/log(p2/p1)*log(px/p1)
              height_tmp=height(kz)+(height(kz+1)-height(kz))*dz*dz1
              call set_z(numpart+jj,height_tmp)
              if (real(part(numpart+jj)%z).gt.height(nz)-0.5) &
                call set_z(numpart+jj, height(nz)-0.5)
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

              !***********************************************
              indzm=nz-1
              indzp=nz
              do i=2,nz
                if (real(height(i),kind=dp).gt.part(numpart+jj)%z) then
                  indzm=i-1
                  indzp=i
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


              ! For domain-filling option 2 (stratospheric O3), 
              ! do the rest only in the stratosphere
              !************************************************

              if (((part(numpart+jj)%z.gt.3000.).and. &
                   (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then

                ! Assign certain properties to the particle
                !******************************************
                part(numpart+jj)%nclass=min( &
                  int(ran1(idummy,0)*real(nclassunc))+1,nclassunc)
                numparticlecount=numparticlecount+1
                part(numpart+jj)%npoint=numparticlecount
                part(numpart+jj)%idt=mintime
                mass(numpart+jj,1)=colmass(lix,ljy)/real(ncolumn)
                if (mdomainfill.eq.2) mass(numpart+jj,1)= &
                     mass(numpart+jj,1)*pvpart*48./29.*ozonescale/10.**9
                mass_init(numpart+jj,1)=mass(numpart+jj,1)
              else
                call terminate_particle(numpart+jj, 0)
                jj=jj-1
                iterminate=iterminate+1
              endif
            endif
          endif
        end do
      end do
      numparttot=numparttot+ncolumn
      if (ipin.eq.0) numpart=numpart+jj
    end do
  end do


  ! alive_tmp=count%alive
  ! spawned_tmp=count%spawned
  ! terminated_tmp=count%terminated

! !$OMP PARALLEL PRIVATE(j) REDUCTION(+:alive_tmp,spawned_tmp,allocated_tmp,terminated_tmp)
 
  ! Make sure that all particles are within domain
  !***********************************************
! !$OMP DO
  do j=1,numpart
    if ((part(j)%xlon.lt.0.).or.(part(j)%xlon.ge.real(nxmin1,kind=dp)).or. &
         (part(j)%ylat.lt.0.).or.(part(j)%ylat.ge.real(nymin1,kind=dp))) then
      call terminate_particle(j,0) ! Cannot be within an OMP region
      iterminate=iterminate+1
      ! alive_tmp=alive_tmp-1
      ! terminated_tmp=terminated_tmp+1
    endif
  end do

  if (iterminate.gt.0) call rewrite_ialive()
! !$OMP END DO
! !$OMP END PARALLEL

  ! count%alive=alive_tmp
  ! count%spawned=spawned_tmp
  ! count%terminated=terminated_tmp
  ! Check whether numpart is really smaller than maxpart
  !*****************************************************

  ! ! ESO :TODO: this warning need to be moved further up, else out-of-bounds error earlier
  !   if (numpart.gt.maxpart) then
  !     write(*,*) 'numpart too large: change source in init_atm_mass.f'
  !     write(*,*) 'numpart: ',numpart,' maxpart: ',maxpart
  !   endif


  xmassperparticle=colmasstotal/real(numparttot)


  ! For boundary conditions, we need fewer particle release heights per column,
  ! because otherwise it takes too long until enough mass has accumulated to
  ! release a particle at the boundary (would take dx/u seconds), leading to
  ! relatively large position errors of the order of one grid distance.
  ! It's better to release fewer particles per column, but to do so more often.
  ! Thus, use on the order of nz starting heights per column.
  ! We thus repeat the above to determine fewer starting heights, that are
  ! used furtheron in subroutine boundcond_domainfill.f.
  !****************************************************************************

  fractus=real(numcolumn)/real(nz)
  write(*,*) 'Total number of particles at model start: ',numpart
  write(*,*) 'Maximum number of particles per column: ',numcolumn
  write(*,*) 'If ',fractus,' <1, better use more particles'
  fractus=sqrt(max(fractus,1.))*0.5

  do ljy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    do lix=nx_we(1),nx_we(2)      ! loop about longitudes
      ncolumn=nint(0.999/fractus*real(npart(1))*colmass(lix,ljy) &
           /colmasstotal)
      if (ncolumn.gt.maxcolumn) stop 'maxcolumn too small'
      if (ncolumn.eq.0) cycle


  ! Memorize how many particles per column shall be used for all boundaries
  ! This is further used in subroutine boundcond_domainfill.f
  ! Use 2 fields for west/east and south/north boundary
  !************************************************************************

      if (lix.eq.nx_we(1)) numcolumn_we(1,ljy)=ncolumn
      if (lix.eq.nx_we(2)) numcolumn_we(2,ljy)=ncolumn
      if (ljy.eq.ny_sn(1)) numcolumn_sn(1,lix)=ncolumn
      if (ljy.eq.ny_sn(2)) numcolumn_sn(2,lix)=ncolumn

  ! Calculate pressure at the altitudes of model surfaces, using the air density
  ! information, which is stored as a 3-d field
  !*****************************************************************************

      do kz=1,nz
        pp(kz)=prs(lix,ljy,kz,1) !rho(lix,ljy,kz,1)*r_air*tt(lix,ljy,kz,1)
      end do

  ! Determine the reference starting altitudes
  !*******************************************

      deltacol=(pp(1)-pp(nz))/real(ncolumn)
      pnew=pp(1)+deltacol*0.5
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

            if (lix.eq.nx_we(1)) zcolumn_we(1,ljy,j)=zposition
            if (lix.eq.nx_we(2)) zcolumn_we(2,ljy,j)=zposition
            if (ljy.eq.ny_sn(1)) zcolumn_sn(1,lix,j)=zposition
            if (ljy.eq.ny_sn(2)) zcolumn_sn(2,lix,j)=zposition

  ! Initialize mass that has accumulated at boundary to zero
  !*********************************************************

            acc_mass_we(1,ljy,j)=0.
            acc_mass_we(2,ljy,j)=0.
            acc_mass_sn(1,ljy,j)=0.
            acc_mass_sn(2,ljy,j)=0.
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

  if (deall) call dealloc_particle(numpart) !Deallocates everything above numpart (F2008)


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

  deallocate(gridarea,colmass)
end subroutine init_domainfill

subroutine boundcond_domainfill(itime,loutend)
  !                                  i      i
  !*****************************************************************************
  !                                                                            *
  ! Particles are created by this subroutine continuously throughout the       *
  ! simulation at the boundaries of the domain-filling box.                    *
  ! All particles carry the same amount of mass which alltogether comprises the*
  ! mass of air within the box, which remains (more or less) constant.         *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 October 2002                                                        *
  !                                                                            *
  !  Changes                                                                   *
  !     2022, L. Bakels: OpenMP parallelisation                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nx_we(2)       grid indices for western and eastern boundary of domain-    *
  !                filling trajectory calculations                             *
  ! ny_sn(2)       grid indices for southern and northern boundary of domain-  *
  !                filling trajectory calculations                             *
  !                                                                            *
  !*****************************************************************************

  use point_mod
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  real :: dz,dz1,dz2,dt1,dt2,dtt,ylat,cosfact,accmasst
  integer :: itime,ii,indz,indzp,i,loutend,numparticlecount_tmp
  integer :: j,k,ix,jy,m,indzh,indexh,ipart,mmass,ithread
  integer :: numactiveparticles,iterminate

  real :: windl(2),rhol(2)
  real :: windhl(2),rhohl(2)
  real :: windx,rhox
  real :: deltaz,boundarea,fluxofmass

  integer :: ixm,ixp,jym,jyp,indzm,mm,iterm_index
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2),yh1(2)

  integer :: idummy = -11


  ! If domain-filling is global, no boundary conditions are needed
  !***************************************************************

  if (gdomainfill) return

  ! Determine auxiliary variables for time interpolation
  !*****************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  numactiveparticles=0
  numparticlecount_tmp=numparticlecount
  accmasst=0.
  ! Terminate trajectories that have left the domain, if domain-filling
  ! trajectory calculation domain is not global
  !********************************************************************
  if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0) call rewrite_iterm()

  iterminate=0
  do i=1,count%allocated
    if (.not. part(i)%alive) cycle

    if ((part(i)%ylat.gt.real(ny_sn(2))).or.(part(i)%ylat.lt.real(ny_sn(1)))) then
      call terminate_particle(i,itime)
      iterminate=iterminate+1
    endif
    if (((.not.xglobal).or.(nx_we(2).ne.(nx-2))).and. &
         ((part(i)%xlon.lt.real(nx_we(1))).or. &
         (part(i)%xlon.gt.real(nx_we(2))))) then
      call terminate_particle(i,itime)
      iterminate=iterminate+1
    endif
    if (part(i)%alive) numactiveparticles = numactiveparticles+1
  end do
  if (iterminate.gt.0) call rewrite_ialive()
  !***************************************
  ! Western and eastern boundary condition
  !***************************************

  ! Loop from south to north
  !*************************
! OMP doesn't work here yet because particles are being spawned inside of the loop
! !$OMP PARALLEL PRIVATE(i,jy,k,j,ii,deltaz,boundarea,indz,indzp,indexh,windl,rhol, &
! !$OMP windhl,rhohl,windx,rhox,fluxofmass,mmass,ixm,jym,ixp,jyp,ddx,ddy,rddx, &
! !$OMP rddy,p1,p2,p3,p4,indzm,mm,indzh,pvpart,ylat,ix,cosfact,ipart) &
! !$OMP REDUCTION(+:numactiveparticles,numparticlecount_tmp,accmasst)

! #ifdef _OPENMP
!   ithread = OMP_GET_THREAD_NUM()
! #else
!   ithread = 0
! #endif
  ithread=0
  iterm_index=1
! !$OMP DO
  do jy=ny_sn(1),ny_sn(2)

  ! Loop over western (index 1) and eastern (index 2) boundary
  !***********************************************************

    do k=1,2

  ! Loop over all release locations in a column
  !********************************************

      do j=1,numcolumn_we(k,jy)

  ! Determine, for each release location, the area of the corresponding boundary
  !*****************************************************************************

        if (j.eq.1) then
          deltaz=(zcolumn_we(k,jy,2)+zcolumn_we(k,jy,1))*0.5
        else if (j.eq.numcolumn_we(k,jy)) then
  ! In order to avoid taking a very high column for very many particles,
  ! use the deltaz from one particle below instead
          deltaz=(zcolumn_we(k,jy,j)-zcolumn_we(k,jy,j-2))*0.5
        else
          deltaz=(zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))*0.5
        endif
        if ((jy.eq.ny_sn(1)).or.(jy.eq.ny_sn(2))) then
          boundarea=deltaz*111198.5/2.*dy
        else
          boundarea=deltaz*111198.5*dy
        endif


  ! Interpolate the wind velocity and density to the release location
  !******************************************************************

  ! Determine the model level below the release position
  !*****************************************************
        indz=nz-1
        indzp=nz
        do i=2,nz
          if (height(i).gt.zcolumn_we(k,jy,j)) then
            indz=i-1
            indzp=i
            exit
          endif
        end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

        dz1=zcolumn_we(k,jy,j)-height(indz)
        dz2=height(indzp)-zcolumn_we(k,jy,j)
        dz=1./(dz1+dz2)

  ! Vertical and temporal interpolation
  !************************************

        do m=1,2
          indexh=memind(m)
          do ii=1,2
            indzh=indz+ii-1
            windl(ii)=uu(nx_we(k),jy,indzh,indexh)
            rhol(ii)=rho(nx_we(k),jy,indzh,indexh)
          end do

          windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
          rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
        end do

        windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
        rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

  ! Calculate mass flux
  !********************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)


  ! If the mass flux is directed into the domain, add it to previous mass fluxes;
  ! if it is out of the domain, set accumulated mass flux to zero
  !******************************************************************************

        if (k.eq.1) then
          if (fluxofmass.ge.0.) then
            acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+fluxofmass
          else
            acc_mass_we(k,jy,j)=0.
          endif
        else
          if (fluxofmass.le.0.) then
            acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+abs(fluxofmass)
          else
            acc_mass_we(k,jy,j)=0.
          endif
        endif
        accmasst=accmasst+acc_mass_we(k,jy,j)

  ! If the accumulated mass exceeds half the mass that each particle shall carry,
  ! one (or more) particle(s) is (are) released and the accumulated mass is
  ! reduced by the mass of this (these) particle(s)
  !******************************************************************************

        if (acc_mass_we(k,jy,j).ge.xmassperparticle*0.5) then
          mmass=int((acc_mass_we(k,jy,j)+xmassperparticle*0.5)/ &
               xmassperparticle)
          acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)- &
               real(mmass)*xmassperparticle
        else
          mmass=0
        endif

        do m=1,mmass
          !THIS WILL CAUSE PROBLEMS WITH OMP! because of dynamical allocation
          call get_newpart_index(ipart,iterm_index)
          call spawn_particle(itime, ipart)

  ! Assign particle positions
  !**************************

          call set_xlon(ipart,real(nx_we(k),kind=dp))
          if (jy.eq.ny_sn(1)) then
            call set_ylat(ipart,real(real(jy)+0.5*ran1(idummy,ithread),kind=dp))
          else if (jy.eq.ny_sn(2)) then
            call set_ylat(ipart,real(real(jy)-0.5*ran1(idummy,ithread),kind=dp))
          else
            call set_ylat(ipart,real(real(jy)+(ran1(idummy,ithread)-.5),kind=dp))
          endif
          if (j.eq.1) then
            call set_z(ipart,zcolumn_we(k,jy,1)+(zcolumn_we(k,jy,2)- &
                 zcolumn_we(k,jy,1))/4.)
          else if (j.eq.numcolumn_we(k,jy)) then
            call set_z(ipart,(2.*zcolumn_we(k,jy,j)+ &
                 zcolumn_we(k,jy,j-1)+height(nz))/4.)
          else
            call set_z(ipart,zcolumn_we(k,jy,j-1)+ran1(idummy,ithread)* &
                 (zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1)))
          endif
#ifdef ETA
          call update_z_to_zeta(itime, ipart)
#endif
  ! Interpolate PV to the particle position
  !****************************************
          ixm=int(part(ipart)%xlon)
          jym=int(part(ipart)%ylat)
          ixp=ixm+1
          jyp=jym+1
          ddx=real(part(ipart)%xlon)-real(ixm)
          ddy=real(part(ipart)%ylat)-real(jym)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy
          indzm=nz-1
          indzp=nz
          do i=2,nz
            if (real(height(i),kind=dp).gt.part(ipart)%z) then
              indzm=i-1
              indzp=i
              exit
            endif
          end do
          dz1=real(part(ipart)%z)-height(indzm)
          dz2=height(indzp)-real(part(ipart)%z)
          dz=1./(dz1+dz2)
          do mm=1,2
            indexh=memind(mm)
            do ii=1,2
              indzh=indzm+ii-1
              y1(ii)=p1*pv(ixm,jym,indzh,indexh) &
                   +p2*pv(ixp,jym,indzh,indexh) &
                   +p3*pv(ixm,jyp,indzh,indexh) &
                   +p4*pv(ixp,jyp,indzh,indexh)
            end do
            yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
          end do
          pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
          ylat=ylat0+real(part(ipart)%ylat)*dy
          if (ylat.lt.0.) pvpart=-1.*pvpart


  ! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
  !*****************************************************************************

          if (((part(ipart)%z.gt.3000.).and. &
               (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
            part(ipart)%nclass=min(int(ran1(idummy,ithread)* &
                 real(nclassunc))+1,nclassunc)
            numactiveparticles=numactiveparticles+1
            numparticlecount_tmp=numparticlecount_tmp+1
            part(ipart)%npoint=numparticlecount_tmp
            part(ipart)%idt=mintime
            part(ipart)%tstart=itime
            mass(ipart,1)=xmassperparticle
            if (mdomainfill.eq.2) mass(ipart,1)= &
                 mass(ipart,1)*pvpart*48./29.*ozonescale/10.**9
            mass_init(ipart,1)=mass(ipart,1)
          else
            stop 'boundcond_domainfill error: look into original to understand what should happen here'
          endif
        end do ! particles
      end do ! release locations in column
    end do ! western and eastern boundary
  end do ! south to north
! !$OMP END DO

  !*****************************************
  ! Southern and northern boundary condition
  !*****************************************

  ! Loop from west to east
  !***********************
! !$OMP DO
  do ix=nx_we(1),nx_we(2)

  ! Loop over southern (index 1) and northern (index 2) boundary
  !*************************************************************

    do k=1,2
      ylat=ylat0+real(ny_sn(k))*dy
      cosfact=cos(ylat*pi180)

  ! Loop over all release locations in a column
  !********************************************

      do j=1,numcolumn_sn(k,ix)

  ! Determine, for each release location, the area of the corresponding boundary
  !*****************************************************************************

        if (j.eq.1) then
          deltaz=(zcolumn_sn(k,ix,2)+zcolumn_sn(k,ix,1))*0.5
        else if (j.eq.numcolumn_sn(k,ix)) then
  !        deltaz=height(nz)-(zcolumn_sn(k,ix,j-1)+
  !    +        zcolumn_sn(k,ix,j))/2.
  ! In order to avoid taking a very high column for very many particles,
  ! use the deltaz from one particle below instead
          deltaz=(zcolumn_sn(k,ix,j)-zcolumn_sn(k,ix,j-2))*0.5
        else
          deltaz=(zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))*0.5
        endif
        if ((ix.eq.nx_we(1)).or.(ix.eq.nx_we(2))) then
          boundarea=deltaz*111198.5*0.5*cosfact*dx
        else
          boundarea=deltaz*111198.5*cosfact*dx
        endif


  ! Interpolate the wind velocity and density to the release location
  !******************************************************************

  ! Determine the model level below the release position
  !*****************************************************
        indz=nz-1
        indzp=nz
        do i=2,nz
          if (height(i).gt.zcolumn_sn(k,ix,j)) then
            indz=i-1
            indzp=i
            exit
          endif
        end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

        dz1=zcolumn_sn(k,ix,j)-height(indz)
        dz2=height(indzp)-zcolumn_sn(k,ix,j)
        dz=1./(dz1+dz2)

  ! Vertical and temporal interpolation
  !************************************

        do m=1,2
          indexh=memind(m)
          do ii=1,2
            indzh=indz+ii-1
            windl(ii)=vv(ix,ny_sn(k),indzh,indexh)
            rhol(ii)=rho(ix,ny_sn(k),indzh,indexh)
          end do

          windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
          rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
        end do

        windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
        rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

  ! Calculate mass flux
  !********************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)

  ! If the mass flux is directed into the domain, add it to previous mass fluxes;
  ! if it is out of the domain, set accumulated mass flux to zero
  !******************************************************************************

        if (k.eq.1) then
          if (fluxofmass.ge.0.) then
            acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+fluxofmass
          else
            acc_mass_sn(k,ix,j)=0.
          endif
        else
          if (fluxofmass.le.0.) then
            acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+abs(fluxofmass)
          else
            acc_mass_sn(k,ix,j)=0.
          endif
        endif
        accmasst=accmasst+acc_mass_sn(k,ix,j)

  ! If the accumulated mass exceeds half the mass that each particle shall carry,
  ! one (or more) particle(s) is (are) released and the accumulated mass is
  ! reduced by the mass of this (these) particle(s)
  !******************************************************************************

        if (acc_mass_sn(k,ix,j).ge.xmassperparticle*0.5) then
          mmass=int((acc_mass_sn(k,ix,j)+xmassperparticle*0.5)/ &
               xmassperparticle)
          acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)- &
               real(mmass)*xmassperparticle
        else
          mmass=0
        endif

        do m=1,mmass
          call get_newpart_index(ipart,iterm_index)
          call spawn_particle(itime, ipart)
  
  ! Assign particle positions
  !**************************
          call set_ylat(ipart,real(ny_sn(k),kind=dp))
          if (ix.eq.nx_we(1)) then
            call set_xlon(ipart,real(real(ix)+0.5*ran1(idummy,ithread),kind=dp))
          else if (ix.eq.nx_we(2)) then
            call set_xlon(ipart,real(real(ix)-0.5*ran1(idummy,ithread),kind=dp))
          else
            call set_xlon(ipart,real(real(ix)+(ran1(idummy,ithread)-.5),kind=dp))
          endif
          if (j.eq.1) then
            call set_z(ipart,zcolumn_sn(k,ix,1)+(zcolumn_sn(k,ix,2)- &
                 zcolumn_sn(k,ix,1))/4.)
          else if (j.eq.numcolumn_sn(k,ix)) then
            call set_z(ipart,(2.*zcolumn_sn(k,ix,j)+ &
                 zcolumn_sn(k,ix,j-1)+height(nz))/4.)
          else
            call set_z(ipart,zcolumn_sn(k,ix,j-1)+ran1(idummy,ithread)* &
                 (zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1)))
          endif
#ifdef ETA
          call update_z_to_zeta(itime, ipart)
#endif
  ! Interpolate PV to the particle position
  !****************************************
          ixm=int(part(ipart)%xlon)
          jym=int(part(ipart)%ylat)
          ixp=ixm+1
          jyp=jym+1
          ddx=real(part(ipart)%xlon)-real(ixm)
          ddy=real(part(ipart)%ylat)-real(jym)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy
          indzm=nz-1
          indzp=nz
          do i=2,nz
            if (real(height(i),kind=dp).gt.part(ipart)%z) then
              indzm=i-1
              indzp=i
              exit
            endif
          end do
          dz1=real(part(ipart)%z)-height(indzm)
          dz2=height(indzp)-real(part(ipart)%z)
          dz=1./(dz1+dz2)
          do mm=1,2
            indexh=memind(mm)
            do ii=1,2
              indzh=indzm+ii-1
              y1(ii)=p1*pv(ixm,jym,indzh,indexh) &
                   +p2*pv(ixp,jym,indzh,indexh) &
                   +p3*pv(ixm,jyp,indzh,indexh) &
                   +p4*pv(ixp,jyp,indzh,indexh)
            end do
            yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
          end do
          pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
          if (ylat.lt.0.) pvpart=-1.*pvpart


  ! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
  !*****************************************************************************

          if (((part(ipart)%z.gt.3000.).and. &
               (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
            part(ipart)%nclass=min(int(ran1(idummy,ithread)* &
                 real(nclassunc))+1,nclassunc)
            numactiveparticles=numactiveparticles+1
            numparticlecount_tmp=numparticlecount_tmp+1
            part(ipart)%npoint=numparticlecount_tmp
            part(ipart)%idt=mintime
            mass(ipart,1)=xmassperparticle
            if (mdomainfill.eq.2) mass(ipart,1)= &
                 mass(ipart,1)*pvpart*48./29.*ozonescale/10.**9
            mass_init(ipart,1)=mass(ipart,1)
          else
            stop 'boundcond_domainfill error: look into original to understand what should happen here'
          endif
        end do ! particles
      end do ! releases per column
    end do ! east west
  end do ! north south
! !$OMP END DO
! !$OMP END PARALLEL
  if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0) call rewrite_iterm()
  numparticlecount = numparticlecount_tmp
  ! If particles shall be dumped, then accumulated masses at the domain boundaries
  ! must be dumped, too, to be used for later runs
  !*****************************************************************************

  if ((ipout.gt.0).and.(itime.eq.loutend)) then
    open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
         form='unformatted')
    write(unitboundcond) numcolumn_we,numcolumn_sn, &
         zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
    close(unitboundcond)
  endif
end subroutine boundcond_domainfill

end module initialise_mod
