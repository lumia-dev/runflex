! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

! DJM - 2017-05-09 - added #ifdef USE_MPIINPLACE cpp directive to     *
! enable declaration of a gridunc0 array if required by MPI code in   *
! mpi_mod.f90                                                         *
!                                                                     *
!**********************************************************************

module unc_mod

  use par_mod, only:dp,dep_prec,nclassunc

  implicit none

  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: gridunc
  real(dep_prec),allocatable, dimension (:,:,:)         :: gridcnt
#ifdef USE_MPIINPLACE
#else
  ! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
  ! then an aux array is needed for parallel grid reduction
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: gridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: griduncn0
#endif
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: griduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn
#ifdef _OPENMP
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:,:) :: gridunc_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:,:) :: griduncn_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: drygridunc_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: drygriduncn_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: wetgridunc_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: wetgriduncn_omp
  real(dep_prec),allocatable, dimension (:,:,:,:)       :: gridcnt_omp
#endif
! For sum of individual contributions, used for the MPI version
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn0

contains

subroutine alloc_grid_unc()
  use com_mod

  implicit none 

  integer :: stat
  
  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(gridunc(0:numxgrid-1,0:numygrid-1,numzgrid,nspec, &
       maxpointspec_act,nclassunc,nageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(gridcnt(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate gridcnt'
#ifdef _OPENMP
  allocate(gridunc_omp(0:numxgrid-1,0:numygrid-1,numzgrid,nspec, &
       maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
  if (stat.ne.0) then
    write(*,*)'ERROR: could not allocate gridunc_omp'
    write(*,*)'increase the memory or reduce MAXTHREADGRID in COMMAND.'
    error stop
  endif
  allocate(gridcnt_omp(0:numxgrid-1,0:numygrid-1,numzgrid,numthreads_grid),stat=stat)
  if (stat.ne.0) then
    write(*,*)'ERROR: could not allocate gridcnt_omp'
    write(*,*)'increase the memory or reduce MAXTHREADGRID in COMMAND.'
    error stop
  endif
#endif
  if (ldirect.gt.0) then
    allocate(wetgridunc(0:numxgrid-1,0:numygrid-1,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc'
    allocate(drygridunc(0:numxgrid-1,0:numygrid-1,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc'
#ifdef _OPENMP
    allocate(wetgridunc_omp(0:numxgrid-1,0:numygrid-1,nspec, &
         maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc_omp'
    allocate(drygridunc_omp(0:numxgrid-1,0:numygrid-1,nspec, &
         maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc_omp'
#endif
  endif

#ifdef USE_MPIINPLACE
#else
! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
! then an aux array is needed for parallel grid reduction
    allocate(gridunc0(0:numxgrid-1,0:numygrid-1,numzgrid,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(gridunc0(1,1,1,1,1,1,1),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  end if
#endif
  if (ldirect.gt.0) then
    if (lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(0:numxgrid-1,0:numygrid-1,nspec, &
           maxpointspec_act,nclassunc,nageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc0'
      allocate(drygridunc0(0:numxgrid-1,0:numygrid-1,nspec, &
           maxpointspec_act,nclassunc,nageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc0'

  ! allocate a dummy to avoid compilator complaints
    else if (.not.lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(1,1,1,1,1,1),stat=stat)
      allocate(drygridunc0(1,1,1,1,1,1),stat=stat)
    end if
  end if

  if (lroot) then
    write (*,*) 'Allocating fields for global output (x,y): ', &
      numxgrid,numygrid
  end if

end subroutine alloc_grid_unc

subroutine alloc_grid_unc_nest()
  use com_mod

  implicit none

  integer :: stat

  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(griduncn(0:numxgridn-1,0:numygridn-1,numzgrid,nspec, &
       maxpointspec_act,nclassunc,nageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
#ifdef _OPENMP
  allocate(griduncn_omp(0:numxgridn-1,0:numygridn-1,numzgrid,nspec, &
       maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc_omp'
#endif

  if (ldirect.gt.0) then
    allocate(wetgriduncn(0:numxgridn-1,0:numygridn-1,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested wetgridunc'
    allocate(drygriduncn(0:numxgridn-1,0:numygridn-1,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested drygridunc'
#ifdef _OPENMP
    allocate(wetgriduncn_omp(0:numxgridn-1,0:numygridn-1,nspec, &
         maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested wetgridunc_omp'
    allocate(drygriduncn_omp(0:numxgridn-1,0:numygridn-1,nspec, &
         maxpointspec_act,nclassunc,nageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested drygriduncn_omp'
#endif
  endif

#ifdef USE_MPIINPLACE
#else
  ! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
  ! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid_nest(),
  ! then an aux array is needed for parallel grid reduction
    allocate(griduncn0(0:numxgridn-1,0:numygridn-1,numzgrid,nspec, &
         maxpointspec_act,nclassunc,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  ! allocate a dummy to avoid compilator complaints
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(griduncn0(1,1,1,1,1,1,1),stat=stat)
  end if
#endif
  if (ldirect.gt.0) then
    if (lroot.and.mpi_mode.gt.0) then
      allocate(wetgriduncn0(0:numxgridn-1,0:numygridn-1,nspec, &
           maxpointspec_act,nclassunc,nageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
      allocate(drygriduncn0(0:numxgridn-1,0:numygridn-1,nspec, &
           maxpointspec_act,nclassunc,nageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  !  endif
  ! allocate a dummy to avoid compilator complaints
    else if (.not.lroot.and.mpi_mode.gt.0) then
      allocate(wetgriduncn0(1,1,1,1,1,1),stat=stat)
      allocate(drygriduncn0(1,1,1,1,1,1),stat=stat)
    end if
  end if

  if (lroot) then
    write (*,*) 'Allocating fields for nested output (x,y): ', &
      numxgridn,numygridn
  end if

end subroutine alloc_grid_unc_nest

subroutine deposit_decay()
  ! Accumulated deposited mass radioactively decays
  use com_mod

  implicit none

  integer ::                &
    j,i,                    & ! loop variable over grid
    ks,                     & ! loop variable species
    kp,                     & ! loop variable for maxpointspec_act
    l,                      & ! loop variable over nclassunc
    nage                      ! loop variable over age classes

!$OMP PARALLEL PRIVATE(ks,kp,nage,l,j,i)
!$OMP DO COLLAPSE(2)
  do ks=1,nspec
  do kp=1,maxpointspec_act
    if (decay(ks).gt.0.) then
      do nage=1,nageclass
        do l=1,nclassunc
  ! Mother output grid
          do j=0,numygrid-1
            do i=0,numxgrid-1
              wetgridunc(i,j,ks,kp,l,nage)= &
                   wetgridunc(i,j,ks,kp,l,nage)* &
                   exp(-1.*outstep*decay(ks))
              drygridunc(i,j,ks,kp,l,nage)= &
                   drygridunc(i,j,ks,kp,l,nage)* &
                   exp(-1.*outstep*decay(ks))
            end do
          end do
  ! Nested output grid
          if (nested_output.eq.1) then
            do j=0,numygridn-1
              do i=0,numxgridn-1
                wetgriduncn(i,j,ks,kp,l,nage)= &
                     wetgriduncn(i,j,ks,kp,l,nage)* &
                     exp(-1.*outstep*decay(ks))
                drygriduncn(i,j,ks,kp,l,nage)= &
                     drygriduncn(i,j,ks,kp,l,nage)* &
                     exp(-1.*outstep*decay(ks))
              end do
            end do
          endif
        end do
      end do
    endif
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine deposit_decay

end module unc_mod
