
  !*****************************************************************************
  !                                                                            *
  ! Module that organises particle information in derived types                *
  ! Particles are terminated and spawned using routines from this module so    *
  ! that global particle counts (spawned, allocated, terminated) are kept      *
  ! up to date                                                                 *
  !                                                                            *
  ! Author: L. Bakels 2022                                                     *
  !                                                                            *
  !*****************************************************************************

module particle_mod
  use com_mod, only: maxspec,DRYDEP,WETDEP,DRYBKDEP,WETBKDEP,iout, &
    n_average,nspec,ipout,ipin
  use par_mod, only: dp,ispeed

  implicit none
  
  type :: coordinates
    real(kind=dp) ::              &
      xlon,                       & ! longitude in grid coordinates
      ylat                          ! latitude in grid coordinates
    real          ::              &
      z                             ! height in meters
#ifdef ETA
    real          ::              &
      zeta                          ! height in eta (ECMWF) coordinates
#endif
  end type coordinates

  type :: velocities
    real               ::         &
      u,                          & ! x velocity
      v,                          & ! y velocity
      w                             ! z velocity
#ifdef ETA
    real               ::         &
      weta                          ! z velocity in eta (ECMWF) coordinates
#endif
  end type velocities

  type :: particle
    real(kind=dp)      ::         &
      xlon,                       & ! Longitude in grid coordinates
      ylat,                       & ! Latitude in grid coordinates
      xlon_prev, ylat_prev,       & ! Keeping the previous positions in memory
      z,                          & ! height in meters
      z_prev                        ! Previous position
#ifdef ETA
    real(kind=dp)      ::         &
      zeta,                       & ! Height in eta (ECMWF) coordinates
      zeta_prev                     ! Previous position
#endif
    type(velocities)   ::         &
      !vel,                        & ! Velocities from interpolated windfields
      turbvel,                    & ! Random turbulent velocities
      mesovel                       ! Mesoscale turbulent velocities
    real               ::         &
      settling                      ! Settling velocity for dry and wet(?) deposit
    logical            ::         &
      alive=.false.,              & ! Flag to show if the particle is still in the running
      spawned=.false.,            & ! Flag to show that particle has spawned (includes terminated parts)
      nstop=.false.                 ! Flag to stop particle (used in advance, stopped in timemanager)
#ifdef ETA
    logical            ::         &
      etaupdate=.false.,          & ! If false, z(meter) is more up-to-date than z(eta)
      meterupdate=.false.           ! If false, z(eta) is more up-to-date than z(meter)
#endif
    integer(kind=2)    ::         &
      icbt                          ! Forbidden state flag   
    integer            ::         &
      tstart,                     & ! spawning time in seconds after start
      tend,                       & ! termination time in seconds after start
      npoint,                     & ! release point
      nclass,                     &
      !species(maxspec),           & ! the number of the corresponding species file of the particle
      idt                           ! internal time of the particle
    ! real,allocatable,dimension(:) ::  &
    !   mass,                       & ! Particle mass for each particle species
    !   mass_init,                  & ! Initial mass of each particle
    !   wetdepo,                    & ! Wet deposition (cumulative)
    !   drydepo,                    & ! Dry deposition (cumulative)
    !   prob                          ! Probability of absorption at ground due to dry deposition
    
    ! real,allocatable   ::         &
    !   val_av(:)                     ! Averaged values; only used when average_output=.true.
    real               ::         &
      ntime=0.,                   & ! Number of timesteps to average over
      cartx_av=0.,                & ! Averaged x pos;
      carty_av=0.,                & ! Averaged y pos;
      cartz_av=0.                   ! Averaged z pos;

  end type particle

  type :: particlecount          
    integer              ::       &
      alive=0,                    & ! Number of particles that are alive
      spawned=0,                  & ! Total number of spawned particles
      terminated=0,               & ! Total number of particles that have been terminated
      allocated=0,                & ! Number of total allocated particle spaces
      iterm_max=0,                & ! Number of empty spaces for overwriting particles
      ninmem=0                      ! Number of particles currently in memory
    logical,allocatable  ::       &
      inmem(:)                      ! Logical to keep track which particle numbers are allocated
    integer,allocatable  ::       &                  
      ialive(:),                  & ! Array that stores alive particle numbers up to count%alive for OMP loops 
      iterm(:)                      ! Array that stores terminated particle numbers up to count%allocated
  end type

  type(particle), allocatable ::  &
    part(:)                         ! This is where all particles are being stored
  type(particlecount)         ::  &
    count                           ! Keeping track of global particle number within the simulation
  real,allocatable            ::  &
    val_av(:,:),                  &  ! Averaged values; only used when average_output=.true.
    xscav_frac1(:,:),             &  ! Only allocated when wet or dry deposit backward mode is switched on
    mass(:,:),                    &  ! mass 
    mass_init(:,:),               &
    wetdeposit(:,:),                 &
    drydeposit(:,:),                 &
    prob(:,:)

  real,allocatable            ::  &
    xplum(:),yplum(:),zplum(:)      ! Only allocated for iout=4 or 5 (plumetraj)
  integer,allocatable         ::  &
    nclust(:)                       ! Only allocated for iout=4 or 5 (plumetraj)
  ! private ::                      &
  !   count             
  public ::                       &
    particle,                     &
    part,                         &
    alloc_particles,              &
    dealloc_particle_range,       &
    dealloc_particle,             &
    dealloc_all_particles,        &
    terminate_particle,           &
    rewrite_ialive,               &
    rewrite_iterm,                &
    spawn_particle,               &
    spawn_particles,              &
    get_totalpart_num,            &
    get_alivepart_num,            &
    get_newpart_index,            &
    particle_allocated,           &
    update_xlon,                  &
    update_ylat,                  &
    update_z,                     &
    count,                        &
    val_av,                       &
    xscav_frac1,                  &
    mass,                         &
    mass_init,                    &
    wetdeposit,                   &
    drydeposit,                   &
    prob,                         &
    xplum,                        &
    yplum,                        &
    zplum,                        &
    nclust                       

  interface update_xlon
    procedure update_xlon_dp, update_xlon_sp, update_xlon_int
  end interface update_xlon

  interface set_xlon
    procedure set_xlon_dp, set_xlon_sp, set_xlon_int
  end interface set_xlon

  interface update_ylat
    procedure update_ylat_dp, update_ylat_sp, update_ylat_int
  end interface update_ylat

  interface set_ylat
    procedure set_ylat_dp, set_ylat_sp, set_ylat_int
  end interface set_ylat

  interface update_z
    procedure update_z_dp,update_z_sp
  end interface update_z

  interface set_z
    procedure set_z_dp,set_z_sp
  end interface set_z

#ifdef ETA
  interface update_zeta
    procedure update_zeta_dp,update_zeta_sp
  end interface update_zeta

  interface set_zeta
    procedure set_zeta_dp,set_zeta_sp
  end interface set_zeta
#endif
contains

  logical function particle_allocated(ipart)
    !******************************************
    ! Checks if the memory of the particle is *
    ! still allocated                         *
    !******************************************

    implicit none 

    integer, intent(in)    :: ipart   ! Particle index
    !logical :: is_particle_allocated
    
    if (ipart.gt.count%allocated) then
      particle_allocated = .false.
    else
      particle_allocated = count%inmem(ipart)
    endif
  end function particle_allocated

  subroutine get_newpart_index(ipart,iterm_index)
    !**************************************************
    ! Returns the first free spot to put a new particle
    !**************************************************
    implicit none

    integer, intent(inout) :: ipart   ! First free index
    integer, intent(inout) :: iterm_index
    integer :: i

    if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0 ) then
      if ((ipin.eq.0 .and. count%terminated.eq.0) .or. &
        (count%allocated.gt.count%spawned)) then
        ipart = count%spawned + 1
      else if (iterm_index.le.count%iterm_max) then
        ! Find dead particles to replace
        if (count%iterm(iterm_index).eq.-1) then
          error stop 'BUG: Attempting to overwrite particle: get_newpart_index.'
        endif
        ipart=count%iterm(iterm_index)
        count%iterm(iterm_index) = -1
        iterm_index = iterm_index+1
      else
        write(*,*) ipart
        ipart=count%allocated + 1
      endif
    else
      ipart = count%spawned + 1
    endif
  end subroutine get_newpart_index

  subroutine get_totalpart_num(npart)
    !********************************************
    ! Returns total number of particles spawned *
    !********************************************
    implicit none 

    integer, intent(inout) :: npart   ! Number of particles

    npart = count%spawned
  end subroutine get_totalpart_num

  subroutine get_alivepart_num(npart)
    !**********************************************
    ! Returns number of particles currently alive *
    !**********************************************
    implicit none 

    integer, intent(inout) :: npart   ! Number of particles

    npart = count%alive
  end subroutine get_alivepart_num

  subroutine spawn_particles(itime, nmpart)
    !******************************************************
    ! Spawning particles
    !
    ! This routine spawns new particles and allocates the 
    ! memory if necessary.
    !******************************************************
    implicit none 

    integer, intent(in) :: &
      itime,               &  ! spawning time
      nmpart                  ! number of particles that are being spawned
    integer ::             &
      i ,j                    ! loop variable

    ! Check if new memory needs to be allocated 
    !*******************************************
    if (nmpart+count%spawned.gt.count%allocated) then
      call alloc_particles( (nmpart+count%spawned) - count%allocated )
    endif

    ! Set the spawning time for each new particle and mark it as alive
    !*****************************************************************
    part(count%spawned+1:count%spawned+nmpart)%tstart = itime
    part(count%spawned+1:count%spawned+nmpart)%alive = .true.
    part(count%spawned+1:count%spawned+nmpart)%spawned = .true.

    ! Updating the list with alive particle numbers that is used to
    ! loop over when doing particle computations
    !*************************************************************
    j=count%spawned+1
    do i=count%alive+1,count%alive+nmpart
      count%ialive(i)=j
      j = j+1
    end do

    ! Update the number of particles that are currently alive
    !********************************************************
    count%alive = count%alive + nmpart

    ! Update the total number of spawned particles
    !*********************************************
    count%spawned = count%spawned + nmpart

  end subroutine spawn_particles

  subroutine spawn_particle(itime, ipart)
    !******************************************************
    ! Spawning particles
    !
    ! This routine spawns new particles and allocates the 
    ! memory if necessary.
    !******************************************************
    implicit none 

    integer, intent(in) :: &
      itime,               & ! spawning time
      ipart                  ! number of particles that are being spawned

    ! Check if new memory needs to be allocated 
    !*******************************************
    if (.not. particle_allocated(ipart)) call alloc_particle(ipart)

    if (part(ipart)%alive) then
      write(*,*) ipart, count%alive, count%terminated, count%allocated
      error stop 'Attempting to overwrite existing particle'
    endif

    ! Update the number of particles that are currently alive
    !********************************************************
    count%alive = count%alive + 1

    ! Set the spawning time for each new particle and mark it as alive
    !*****************************************************************
    part(ipart)%tstart = itime
    part(ipart)%alive = .true.
    part(ipart)%spawned = .true.

    ! Updating the list with alive particle numbers that is used to
    ! loop over when doing particle computations
    !*************************************************************
    count%ialive(count%alive) = ipart

    ! Update the total number of spawned particles
    !*********************************************
    count%spawned = count%spawned + 1

  end subroutine spawn_particle

  subroutine terminate_particle(ipart,itime)
    !*****************************************************
    ! Terminating specified particle
    !
    ! This routine terminates a selected particle
    !***************************************************** 
    implicit none

    integer, intent(in) :: &
      ipart,               & ! to be terminated particle index
      itime                  ! Time at which particle is terminated

    ! Flagging the particle as having been terminated
    !************************************************
    part(ipart)%alive=.false.  
    part(ipart)%tend=itime

    ! Update the number of current particles that are alive
    !******************************************************
    count%alive = count%alive - 1

    ! Update the total number of terminated particles during the whole run
    !**********************************************************************
    count%terminated = count%terminated + 1
  end subroutine terminate_particle

  subroutine rewrite_ialive()
    implicit none

    integer :: i,j

    j=1
    do i=1,count%allocated
      if (part(i)%alive) then
        count%ialive(j)=i
        j=j+1
      endif
    end do

    count%alive=j-1

    if (ipin.le.1 .and. ipout.eq.0 .and. ispeed.eq.0) call rewrite_iterm
  end subroutine rewrite_ialive

  subroutine rewrite_iterm()
    implicit none

    integer :: i,j

    j=1
    do i=1,count%allocated
      if (.not. part(i)%alive) then
        count%iterm(j)=i
        j=j+1
      endif
    end do

    count%iterm_max=j-1

  end subroutine rewrite_iterm

  subroutine rewrite_ialive_single(ipart)
    implicit none

    integer, intent(in) :: &
      ipart                  ! to be terminated particle index    
    integer ::             &
      i,                   & ! loop variable
      iloc                   ! location of ipart in count%ialive

    ! And remove from the ialive array
    !*********************************
    ! iloc=findloc(count%ialive,ipart,1) ! findloc not supported in gcc<v9
    iloc=count%allocated
    do i=1,count%alive+2
      if (count%ialive(i).eq.ipart) then
        iloc=i
        exit
      endif
    end do
    if (iloc.ne.count%allocated) then
!$OMP PARALLEL
!$OMP WORKSHARE
      count%ialive(iloc:count%allocated-1)=count%ialive(iloc+1:count%allocated)
!$OMP END WORKSHARE
!$OMP END PARALLEL
    endif
  end subroutine rewrite_ialive_single
 
  subroutine alloc_particles(nmpart)

    implicit none 

    integer, intent(in)        :: nmpart
    type(particle),allocatable :: tmppart(:)
    logical, allocatable       :: tmpcount(:)
    real, allocatable          :: tmpxscav(:,:),tmpval_av(:,:),tmpmass(:,:),tmpdepo(:,:)
    real, allocatable          :: tmpxl(:),tmpyl(:),tmpzl(:)
    integer, allocatable       :: tmpnclust(:)
    integer                    :: i,stat

    if (nmpart.gt.100) &
      write(*,*) 'Allocating ',nmpart,' particles', count%allocated, count%terminated, count%spawned

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    allocate( tmpcount(count%allocated+nmpart),stat=stat)
    if (stat.ne.0) error stop "Could not allocate tmpcount"
    if (count%allocated.gt.0) tmpcount(1:count%allocated) = count%inmem
    call move_alloc(tmpcount,count%inmem)
    allocate( tmpnclust(count%allocated+nmpart),stat=stat)
    if (stat.ne.0) error stop "Could not allocate tmpnclust"
    if (count%allocated.gt.0) tmpnclust(1:count%allocated) = count%ialive
    call move_alloc(tmpnclust,count%ialive)

    if (ipin.le.1 .and. ipout.eq.0) then
      allocate( tmpnclust(count%allocated+nmpart),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpnclust"
      if (count%allocated.gt.0) tmpnclust(1:count%allocated) = count%iterm
      call move_alloc(tmpnclust,count%iterm)
    endif 

    count%inmem(count%allocated+1:count%allocated+nmpart) = .true.

    ! Allocating new particle spaces
    !*******************************
    allocate( tmppart(count%allocated+nmpart),stat=stat)
    if (stat.ne.0) error stop "Could not allocate tmppart"

    do i=1,count%allocated+nmpart
      ! allocate( tmppart(i)%mass(maxspec),tmppart(i)%mass_init(maxspec),stat=stat)
      ! if (stat.ne.0) error stop "Could not allocate tmppart"
      ! if (DRYDEP) then
      !   allocate( tmppart(i)%drydepo(maxspec),tmppart(i)%prob(maxspec),stat=stat)
      !   if (stat.ne.0) error stop "Could not allocate tmppart"
      !   tmppart(i)%drydepo(maxspec)=0.
      ! endif
      ! if (WETDEP) then 
      !   allocate( tmppart(i)%wetdepo(maxspec),stat=stat)
      !   if (stat.ne.0) error stop "Could not allocate tmppart"
      !   tmppart(i)%wetdepo(maxspec)=0.
      ! endif
      tmppart(i)%ntime=0 ! Preventing particles to be written to partoutput when they just spawned
    end do
    if (count%allocated.gt.0) tmppart(1:count%allocated) = part
    call move_alloc(tmppart,part)

    allocate( tmpmass(count%allocated+nmpart,maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate tmpmass"
    if (count%allocated.gt.0) tmpmass(1:count%allocated,:) = mass
    call move_alloc(tmpmass,mass)

    allocate( tmpmass(count%allocated+nmpart,maxspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate tmpmass_init"
    if (count%allocated.gt.0) tmpmass(1:count%allocated,:) = mass_init
    call move_alloc(tmpmass,mass_init)

    if (DRYDEP) then
      allocate( tmpdepo(count%allocated+nmpart,maxspec),stat=stat)
      if (stat.ne.0) error stop "Could not allocate drydeposit"
      if (count%allocated.gt.0) tmpdepo(1:count%allocated,:) = drydeposit
      call move_alloc(tmpdepo,drydeposit)
      drydeposit(count%allocated+1:count%allocated+nmpart,:)=0.

      allocate( tmpdepo(count%allocated+nmpart,maxspec),stat=stat)
      if (stat.ne.0) error stop "Could not allocate prob"
      if (count%allocated.gt.0) tmpdepo(1:count%allocated,:) = prob
      call move_alloc(tmpdepo,prob)
      prob(count%allocated+1:count%allocated+nmpart,:)=0.
    endif
    if (WETDEP) then 
      allocate( tmpdepo(count%allocated+nmpart,maxspec),stat=stat)
      if (stat.ne.0) error stop "Could not allocate wetdeposit"
      if (count%allocated.gt.0) tmpdepo(1:count%allocated,:) = wetdeposit
      call move_alloc(tmpdepo,wetdeposit)
      wetdeposit(count%allocated+1:count%allocated+nmpart,:)=0.
    endif

    ! If wet or dry deposition backward mode is switched on, xscav_frac1
    ! needs to be allocated
    !*******************************************************************
    if (WETBKDEP.or.DRYBKDEP) then
      allocate( tmpxscav(count%allocated+nmpart,maxspec),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpxscav"
      if (count%allocated.gt.0) tmpxscav(1:count%allocated,:) = xscav_frac1
      call move_alloc(tmpxscav,xscav_frac1)
      ! Initialise it here
      xscav_frac1(count%allocated+1:count%allocated+nmpart,:) = -1.
    endif

    if (n_average.gt.0) then 
      allocate( tmpval_av(count%allocated+nmpart,n_average),stat=stat)
      if (stat.ne.0) error stop "Could not allocate val_av"
      if (count%allocated.gt.0) tmpval_av(1:count%allocated,:) = val_av
      tmpval_av(count%allocated+1:count%allocated+nmpart,:) = 0
      call move_alloc(tmpval_av, val_av)
    endif

    if ((iout.eq.4).or.(iout.eq.5)) then
      allocate( tmpxl(count%allocated+nmpart),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpxl"
      if (count%allocated.gt.0) tmpxl(1:count%allocated) = xplum
      call move_alloc(tmpxl,xplum)
      
      allocate( tmpyl(count%allocated+nmpart),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpyl"
      if (count%allocated.gt.0) tmpyl(1:count%allocated) = yplum
      call move_alloc(tmpyl,yplum)

      allocate( tmpzl(count%allocated+nmpart),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpzl"
      if (count%allocated.gt.0) tmpzl(1:count%allocated) = zplum
      call move_alloc(tmpzl,zplum)

      allocate( tmpnclust(count%allocated+nmpart),stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmpnclust"
      if (count%allocated.gt.0) tmpnclust(1:count%allocated) = nclust
      call move_alloc(tmpnclust,nclust)
    endif

    count%allocated = count%allocated+nmpart
    if (nmpart.gt.100) write(*,*) 'Finished allocation'
  end subroutine alloc_particles

  subroutine alloc_particle(ipart)

    implicit none 

    integer, intent(in) :: ipart

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    if (ipart.gt.count%allocated) then 
      call alloc_particles(ipart-count%allocated)
    else
      error stop 'Error: You are trying to allocate an already existing particle'
    endif

  end subroutine alloc_particle

  subroutine dealloc_particle_range(istart,iend)

    implicit none

    integer, intent(in) :: istart,iend

    !deallocate( part(istart:iend) )
    count%inmem(istart:iend) = .false.
  end subroutine dealloc_particle_range

  subroutine dealloc_particle(ipart)

    implicit none

    integer, intent(in) :: ipart ! particle index

    !deallocate( part(ipart) )
    part = part(1:ipart) ! FORTRAN 2008 only
    count%inmem(ipart+1:) = .false.
  end subroutine dealloc_particle

  subroutine dealloc_all_particles()

    implicit none

    integer :: i

    if (n_average.gt.0) deallocate( val_av )
    deallocate( part )
    deallocate( count%inmem )
    deallocate( count%ialive )
    if (ipin.le.1 .and. ipout.eq.0) deallocate( count%iterm )
    deallocate( mass, mass_init )

    if (WETBKDEP.or.DRYBKDEP) then
      deallocate( xscav_frac1 )
    endif

    if (WETDEP) deallocate( wetdeposit )
    if (DRYDEP) deallocate( drydeposit,prob )

    if ((iout.eq.4).or.(iout.eq.5)) then
      deallocate( xplum )
      deallocate( yplum )
      deallocate( zplum )
      deallocate( nclust )
    endif
  end subroutine dealloc_all_particles

! Update_xlon
  subroutine update_xlon_dp(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart   ! particle index
    real(kind=dp), intent(in) :: xchange

    part(ipart)%xlon = part(ipart)%xlon + xchange
  end subroutine update_xlon_dp

  subroutine update_xlon_sp(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in) :: ipart    ! particle index
    real, intent(in)    :: xchange

    part(ipart)%xlon = part(ipart)%xlon + real(xchange,kind=dp)
  end subroutine update_xlon_sp

  subroutine update_xlon_int(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in) :: ipart  ! particle index
    integer, intent(in) :: xchange

    part(ipart)%xlon = part(ipart)%xlon + real(xchange,kind=dp)
  end subroutine update_xlon_int
! End Update_xlon

! Set_xlon
  subroutine set_xlon_dp(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: xvalue

    part(ipart)%xlon = xvalue
  end subroutine set_xlon_dp

  subroutine set_xlon_sp(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart   ! particle index
    real, intent(in)       :: xvalue

    part(ipart)%xlon = real(xvalue,kind=dp)
  end subroutine set_xlon_sp

  subroutine set_xlon_int(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    integer, intent(in)    :: xvalue

    part(ipart)%xlon = real(xvalue,kind=dp)
  end subroutine set_xlon_int
! End Set_xlon 

! Update_ylat
  subroutine update_ylat_dp(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: ychange

    part(ipart)%ylat = part(ipart)%ylat + ychange
  end subroutine update_ylat_dp

  subroutine update_ylat_sp(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: ychange

    part(ipart)%ylat = part(ipart)%ylat + real(ychange,kind=dp)
  end subroutine update_ylat_sp

  subroutine update_ylat_int(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart ! particle index
    integer, intent(in)    :: ychange

    part(ipart)%ylat = part(ipart)%ylat + real(ychange,kind=dp)
  end subroutine update_ylat_int
! End Update_ylat

! Set_ylat
  subroutine set_ylat_dp(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: yvalue

    part(ipart)%ylat = yvalue
  end subroutine set_ylat_dp

  subroutine set_ylat_sp(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: yvalue

    part(ipart)%ylat = real(yvalue,kind=dp)
  end subroutine set_ylat_sp

  subroutine set_ylat_int(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    integer, intent(in) :: yvalue

    part(ipart)%ylat = real(yvalue,kind=dp)
  end subroutine set_ylat_int
! End Set_ylat

! Update z positions
  subroutine update_z_dp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zchange

    part(ipart)%z = part(ipart)%z + zchange
#ifdef ETA
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
#endif
  end subroutine update_z_dp

  subroutine update_z_sp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zchange

    part(ipart)%z = part(ipart)%z + real(zchange,kind=dp)
#ifdef ETA
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
#endif
  end subroutine update_z_sp  

#ifdef ETA
  subroutine update_zeta_dp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zchange

    part(ipart)%zeta = part(ipart)%zeta + zchange
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine update_zeta_dp

  subroutine update_zeta_sp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zchange

    part(ipart)%zeta = part(ipart)%zeta + real(zchange,kind=dp)
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine update_zeta_sp
#endif
! End update z positions

! Update z positions
  subroutine set_z_dp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zvalue

    part(ipart)%z = zvalue
#ifdef ETA
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
#endif
  end subroutine set_z_dp  

  subroutine set_z_sp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zvalue

    part(ipart)%z = real(zvalue,kind=dp)
#ifdef ETA
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
#endif
  end subroutine set_z_sp

#ifdef ETA
  subroutine set_zeta_dp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zvalue

    part(ipart)%zeta = zvalue
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine set_zeta_dp

  subroutine set_zeta_sp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zvalue

    part(ipart)%zeta = real(zvalue,kind=dp)
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine set_zeta_sp
#endif
! End update z positions

end module particle_mod
