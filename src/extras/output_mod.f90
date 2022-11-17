module output_mod

    use netcdf
    use netcdf_tools,  only : nf90_err
    use particles_mod, only : particles, pp

    implicit none

    contains

    subroutine output_particles_final_position

        use com_mod, only : numpart, path

        ! NetCDF variables:
        integer             :: ncstat
        integer             :: fid, dimid, varid
        character(len=200)  :: filename

        ! Write the position of the particles at the end of their life (or when the simulation ends).
        write(filename, *) trim(path(2))//'particles_final.nc'

        ! Open file
        ncstat = nf90_create(trim(filename), NF90_NETCDF4, fid)

        ! Dimensions
        print*, numpart
        call nf90_err(nf90_def_dim(fid, 'particles', numpart, dimid))

        ! Variables:
        print*, 'lon'
        print*, allocated(particles)
        print*, shape(particles)
        print*, numpart

        call nf90_err(nf90_def_var(fid, 'lon', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%lon, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'lat', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%lat, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'height', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%z, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'surface_height', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%oro, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'potential_vorticity', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%potential_vorticity, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'specific_humidity', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%specific_humidity, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'air_density', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%density, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'pbl_height', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%pbl_height, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'tropopause_height', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%tropopause_height, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'temperature', nf90_float, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%temperature, .not. particles%free)))

        !call nf90_err(nf90_def_var(fid, 'mass', nf90_float, (/dimid/), varid))
        !call nf90_err(nf90_put_var(fid, varid, particles%mass))

        call nf90_err(nf90_def_var(fid, 'release_number', nf90_int, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%release, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'release_time', nf90_int, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%release_time, .not. particles%free)))

        call nf90_err(nf90_def_var(fid, 'itime', nf90_int, (/dimid/), varid))
        call nf90_err(nf90_put_var(fid, varid, pack(particles%t, .not. particles%free)))

        !call nf90_err(nf90_def_var(fid, 'active', nf90_byte, (/dimid/), varid))
        !call nf90_err(nf90_put_var(fid, varid, particles%active))

        !call nf90_err(nf90_def_var(fid, 'free', nf90_byte, (/dimid/), varid))
        !call nf90_err(nf90_put_var(fid, varid, particles%free))

        call nf90_err(nf90_close(fid))

    end subroutine output_particles_final_position

end module output_mod