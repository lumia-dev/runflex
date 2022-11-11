module particles_mod

    use par_mod, only : dp

    type type_particles
        integer          :: number
        real(kind=dp)    :: x, y    ! horizontal position of the particle (xtra1, ytra1)
        real             :: z       ! vertical position of the particle (ztra1)
        integer          :: t       ! temporal position of the particle (itra1)
        integer          :: release ! release point of the particle (npoint)
        logical          :: active
        logical          :: free
        integer          :: release_time
        real, dimension(:), allocatable :: mass
    end type type_particles

    type(type_particles), dimension(:), allocatable, target :: particles

    type(type_particles), pointer :: pp

    contains

    subroutine init_particles(npart)

        integer, intent(in) :: npart

        allocate(particles(npart))

        particles(:)%t = -999999999
        particles(:)%free = .true.
        particles(:)%active = .false.

    end subroutine init_particles

end module particles_mod
