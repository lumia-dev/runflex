! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module totals_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for injecting mass       *
  !    into particles based on gridded emissions estimates                     *
  !                                                                            *
  !*****************************************************************************

  use netcdf
  use par_mod,           only: dp
  use com_mod
  use netcdf_output_mod, only: nf90_err

  implicit none

    character(len=256) :: fn_totals
    integer            :: nc_id, specdim_id, reagdim_id, timedim_id, nchardim_id
    integer            :: time_id, spec_id, cl_id, emis_id, efld_id, eres_id, mass_id
    real(kind=dp), dimension(:,:), allocatable :: chem_loss
    real(kind=dp), dimension(:), allocatable   :: tot_mass
    real(kind=dp), dimension(:), allocatable   :: tot_em_up
    real(kind=dp), dimension(:), allocatable   :: tot_em_res
    real(kind=dp), dimension(:), allocatable   :: tot_em_field

  contains

  subroutine alloc_totals

  !*****************************************************************************
  !                                                                            *
  !    Allocate variables for totals                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: stat

    allocate(chem_loss(nreagent,nspec),stat=stat)
    if (stat.ne.0) error stop "Could not allocate totals arrays"
    chem_loss(:,:)=0.

    allocate( tot_mass(nspec) )
    allocate( tot_em_up(nspec) )
    allocate( tot_em_field(nspec) )
    allocate( tot_em_res(nspec) )

  end subroutine alloc_totals

  subroutine totals_init()

  !*****************************************************************************
  !                                                                            *
  !    This subroutine initializes the totals output                           *
  !                                                                            *
  !    Author: S. Henne                                                        *
  !    Adapted by R. Thompson for v11, Feb-2024                                *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=10)  :: time
    character(len=8)   :: date
    character(len=5)   :: zone
    character(len=8)   :: adate
    character(len=6)   :: atime
    character(len=32)  :: timeunit
    character(len=256) :: host_name, login_name
    integer :: ks

    ! get string of start time
    write(adate,'(i8.8)') ibdate
    write(atime,'(i6.6)') ibtime
    timeunit = 'seconds since '//adate(1:4)//'-'//adate(5:6)// &
                 '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)

    ! get run info 
    call date_and_time(date,time,zone)
    call getlog(login_name)
    call hostnm(host_name)

    ! file name
    fn_totals=trim(path(2)(1:length(2)))//'totals.nc'

    ! open new file handle
    call nf90_err( nf90_create(trim(fn_totals), cmode=nf90_hdf5, ncid=nc_id) )

    ! define dimensions
    !******************

    call nf90_err( nf90_def_dim(nc_id, "species", nspec, specdim_id) )
    call nf90_err( nf90_def_dim(nc_id, "reagents", nreagent, reagdim_id) )
    call nf90_err( nf90_def_dim(nc_id, "time", nf90_unlimited, timedim_id) )
    call nf90_err( nf90_def_dim(nc_id, "nchar", 10, nchardim_id) )

    ! define variables
    !*****************

    ! time
    call nf90_err( nf90_def_var(nc_id, 'time', nf90_int, (/ timedim_id /), time_id) )
    call nf90_err( nf90_put_att(nc_id, time_id, 'units', timeunit) )
    ! species names
    call nf90_err( nf90_def_var(nc_id, 'species', nf90_char, (/ nchardim_id, specdim_id /), spec_id) )
    call nf90_err( nf90_put_att(nc_id, spec_id, 'long_name', 'Species names') )
    ! total masses
    call nf90_err( nf90_def_var(nc_id, 'mass', nf90_double, (/ specdim_id, timedim_id /), mass_id) )
    call nf90_err( nf90_put_att(nc_id, mass_id, 'units', 'kg') )
    call nf90_err( nf90_put_att(nc_id, mass_id, 'long_name', 'Total global mass') )
    ! emission uptake
    call nf90_err( nf90_def_var(nc_id, 'emissions', nf90_double, (/ specdim_id, timedim_id /), emis_id) )
    call nf90_err( nf90_put_att(nc_id, emis_id, 'units', 'kg s-1') )
    call nf90_err( nf90_put_att(nc_id, emis_id, 'long_name', 'Actual emission flux') )
    ! emission field
    call nf90_err( nf90_def_var(nc_id, 'em_field', nf90_double, (/ specdim_id, timedim_id /), efld_id) )
    call nf90_err( nf90_put_att(nc_id, efld_id, 'units', 'kg s-1') )
    call nf90_err(  nf90_put_att(nc_id, efld_id, 'long_name', 'Emission flux in fields') )
    ! emission residual
    call nf90_err( nf90_def_var(nc_id, 'em_res', nf90_double, (/ specdim_id, timedim_id /), eres_id) )
    call nf90_err( nf90_put_att(nc_id, eres_id, 'units', 'kg') )
    call nf90_err( nf90_put_att(nc_id, eres_id, 'long_name', 'Emission residuals') )
    ! chemical loss
    call nf90_err( nf90_def_var(nc_id, 'chem_loss', nf90_double, (/ reagdim_id, specdim_id, timedim_id /), cl_id) )
    call nf90_err( nf90_put_att(nc_id, cl_id, 'units', 'kg s-1') )
    call nf90_err( nf90_put_att(nc_id, cl_id, 'long_name', 'Mass loss through chemical reactions') )

    ! write global attributes
    !************************

    call nf90_err( nf90_put_att(nc_id, nf90_global, 'title', 'FLEXPART total mass and flux output') )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'source', trim(flexversion)//' model output') )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'history', date(1:4)//'-'//date(5:6)// &
       '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//' '//zone//'  created by '//  &
       trim(login_name)//' on '//trim(host_name)) )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'references', &
       'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200;'//&
       'Henne et al., in Lagrangian Modeling of the Atmosphere, 2012, doi:10.1029/2012GM001247') )

    ! end definition of file 
    call nf90_err( nf90_enddef(nc_id) )

    ! write species info
    do ks=1,nspec
      print*, 'totals: species = ',trim(species(ks))
      call nf90_err( nf90_put_var(nc_id, spec_id, species(ks), (/1,ks/), (/10,1/)) )
    end do

    ! close file
    call nf90_err( nf90_close(nc_id) )

  end subroutine totals_init


  subroutine totals_write(itime)

  !*****************************************************************************
  !                                                                            *
  !    This subroutine writes the totals to file                               *
  !                                                                            *
  !    Author: S. Henne                                                        *
  !    Adapted by R. Thompson for v11, Feb-2024                                *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: itime, tidx
    integer :: var_id

    ! open file
    call nf90_err( nf90_open(trim(fn_totals), nf90_write, nc_id) )

    ! get length of time dimension -> increase index by one to write new data
    call nf90_err( nf90_inq_dimid(nc_id, "time", timedim_id) )
    call nf90_err( nf90_inquire_dimension(nc_id, timedim_id, len=tidx) )
    tidx = tidx + 1

    ! add to time variable
    call nf90_err( nf90_inq_varid(nc_id, "time", time_id) )
    call nf90_err( nf90_put_var(nc_id, time_id, itime, (/ tidx /)) )

    ! write variables
    !****************

    call nf90_err( nf90_inq_varid(nc_id, "mass", var_id) )
    call nf90_err( nf90_put_var(nc_id, var_id, tot_mass(1:nspec), (/ 1, tidx /), (/ nspec, 1/)) )

    call nf90_err( nf90_inq_varid(nc_id, "emissions", var_id) )
    call nf90_err( nf90_put_var(nc_id, var_id, tot_em_up(1:nspec)/real(lsynctime), &
                    (/ 1, tidx /), (/ nspec, 1/)) )

    call nf90_err( nf90_inq_varid(nc_id, "em_field", var_id) )
    call nf90_err( nf90_put_var(nc_id, var_id, tot_em_field(1:nspec)/real(lsynctime), &
                    (/ 1, tidx /), (/ nspec, 1/)) )

    ! em_res accumulated over all time steps (units kg) -> no division by lsynctime
    call nf90_err( nf90_inq_varid(nc_id, "em_res", var_id) )
    call nf90_err( nf90_put_var(nc_id, var_id, tot_em_res(1:nspec),  &
                    (/ 1, tidx /), (/ nspec, 1/)) )

    if (nreagent.gt.0) then
      call nf90_err( nf90_inq_varid(nc_id, "chem_loss", var_id) )
      call nf90_err( nf90_put_var(nc_id, var_id, chem_loss(1:nreagent,1:nspec)/real(lsynctime),   &
                      (/ 1, 1, tidx /), (/ nreagent, nspec, 1/)) )
    endif

    call nf90_err( nf90_close(nc_id) )

  end subroutine totals_write

end module totals_mod
