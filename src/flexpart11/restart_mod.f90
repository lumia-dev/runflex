! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module restart_mod

  !*****************************************************************************
  !                                                                            *
  !    This module write variables to file for eventual restart, and reads     *
  !    these variables from file in case of restart                            *
  !                                                                            *
  !*****************************************************************************

  use particle_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use outgrid_mod
  use unc_mod
  use date_mod
#ifdef USE_NCF
  use netcdf_output_mod
#endif

  character(len=256) :: restart_filename1,restart_filename2,restart_filename3

contains

subroutine output_restart(itime,loutnext,lrecoutnext,outnum)

  implicit none

  integer, intent(in) :: itime,loutnext,lrecoutnext
  real, intent(in) :: outnum
  integer :: imax,i,j,jjjjmmdd,ihmmss,ipart,iwritten
  integer :: ks,kp,kz,nage,jy,ix,l,n
  real(kind=dp) :: jul
  character :: adate*8,atime*6


  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  restart_filename3 = restart_filename2
  restart_filename2 = restart_filename1
  restart_filename1 = path(2)(1:length(2))//'restart_'//adate//atime

  write(*,*) 'Writing Restart file:', trim(restart_filename1)

#ifdef ETA
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
  do j=1,count%alive
    i=count%ialive(j)
    if (part(i)%alive) then
      call update_zeta_to_z(itime,i)
      call update_z_to_zeta(itime,i)
    endif
  end do
!$OMP END DO
!$OMP END PARALLEL
#endif

  open(unitrestart,file=restart_filename1,form='unformatted')

  !Get largest live particle number for allocation new start
  
  imax=0
  iwritten=0
  do ipart=1,count%allocated
    if ((ipout.gt.0).and.(n_average.gt.0)) then
      if((.not. part(ipart)%alive).and.(abs(part(ipart)%tend-itime).ge.ipoutfac*loutstep)) &
        cycle
    else
      if (.not. part(ipart)%alive) cycle
    endif
    if (ipart.gt.imax) imax=ipart
    iwritten=iwritten+1
  end do
  ! Write current time to file
  !***************************
  write(unitrestart) itime
  write(unitrestart) imax
  write(unitrestart) iwritten
  write(unitrestart) loutnext
  write(unitrestart) lrecoutnext
  write(unitrestart) outnum

  do ipart=1,imax
    if ((ipout.gt.0).and.(n_average.gt.0)) then
      if((.not. part(ipart)%alive).and.(abs(part(ipart)%tend-itime).ge.ipoutfac*loutstep)) &
        cycle
    else
      if (.not. part(ipart)%alive) cycle
    endif
    write(unitrestart) ipart
    write(unitrestart) part(ipart)%xlon,part(ipart)%ylat,part(ipart)%z, &
#ifdef ETA
      part(ipart)%zeta, &
#endif
      part(ipart)%npoint,part(ipart)%nclass,part(ipart)%idt,part(ipart)%tend, &
      part(ipart)%tstart,part(ipart)%alive,part(ipart)%turbvel%u, &
      part(ipart)%turbvel%v,part(ipart)%turbvel%w,part(ipart)%mesovel%u, &
      part(ipart)%mesovel%v,part(ipart)%mesovel%w,(mass(ipart,j),j=1,nspec), &
      (mass_init(ipart,j),j=1,nspec)
    if (wetdep) write(unitrestart) (wetdeposit(ipart,j),j=1,nspec)
    if (drydep) write(unitrestart) (drydeposit(ipart,j),j=1,nspec)
    if ((drybkdep).or.(wetbkdep)) write(unitrestart) (xscav_frac1(ipart,j),j=1,nspec)
  end do
  if (iout.gt.0) then
#ifdef USE_NCF
    if (lnetcdfout.eq.1) write(unitrestart) tpointer
#endif

    do ks=1,nspec
      do kp=1,maxpointspec_act
        do nage=1,nageclass
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              do l=1,nclassunc
                do kz=1,numzgrid
                  write(unitrestart) gridunc(ix,jy,kz,ks,kp,l,nage)
                end do
                if ((wetdep).and.(ldirect.gt.0)) then
                  write(unitrestart) wetgridunc(ix,jy,ks,kp,l,nage)
                endif
                if ((drydep).and.(ldirect.gt.0)) then
                  write(unitrestart) drygridunc(ix,jy,ks,kp,l,nage)
                endif
              end do
            end do
          end do
          if (nested_output.eq.1) then
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                do l=1,nclassunc
                  do kz=1,numzgrid
                    write(unitrestart) griduncn(ix,jy,kz,ks,kp,l,nage)
                  end do
                  if ((wetdep).and.(ldirect.gt.0)) then
                    write(unitrestart) wetgriduncn(ix,jy,ks,kp,l,nage)
                  endif
                  if ((drydep).and.(ldirect.gt.0)) then
                    write(unitrestart) drygriduncn(ix,jy,ks,kp,l,nage)
                  endif
                end do
              end do
            end do
          endif
          if (iflux.eq.1) then
            do kz=1,numzgrid
              do jy=0,numygridn-1
                do ix=0,numxgridn-1
                  do i=1,5
                    write(unitrestart) flux(i,ix,jy,kz,ks,kp,nage)
                  end do
                end do
              end do
            end do        
          endif
        end do
        if (linit_cond.gt.0) then
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                write(unitrestart) init_cond(ix,jy,kz,ks,kp)
              end do 
            end do
          end do 
        endif
      end do
      if (numreceptor.gt.0) then 
        do n=1,numreceptor
          read(unitpartin) creceptor(n,ks)
        end do
      endif
    end do
  endif
  close(unitrestart)

  ! open(unit=1234, iostat=stat, file=restart_filename3, status='old')
  ! if(stat == 0) close(1234, status='delete')
end subroutine output_restart

subroutine readrestart

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the particle dump file and reads all the particle     *
  !   positions and gridded information from a previous run to initialize      *
  !   the current run.                                                         *
  !                                                                            *
  !     Author: L. Bakels 2022                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,j,ipart,ios,iterminate
  integer :: id1,id2,it1,it2,imax
  integer :: ks,kp,kz,nage,jy,ix,l,n
  real(kind=dp) :: julin
  real :: a

  numparticlecount=0


  open(unitpartin,file=path(2)(1:length(2))//'restart.bin', &
       form='unformatted',err=9989)

  write(*,*) 'Reading Restart file:', path(2)(1:length(2))//'restart.bin'
  
  read(unitpartin,iostat=ios) itime_init
  read(unitpartin) imax ! count%alive
  read(unitpartin) numpart
  read(unitpartin) loutnext_init
  read(unitpartin) lrecoutnext_init
  read(unitpartin) outnum_init

  count%alive=numpart
  if (ipin.eq.1) then
    count%spawned=imax
  else
    count%spawned=numpart
  endif
  if (count%allocated.lt.imax) call alloc_particles(imax-count%allocated)
  do i=1,numpart
    read(unitpartin) ipart
    if (ipout.gt.0) ipart=i ! No need to keep dead particle spots when no part dump
    read(unitpartin) part(ipart)%xlon,part(ipart)%ylat,part(ipart)%z, &
#ifdef ETA
      part(ipart)%zeta, &
#endif
      part(ipart)%npoint,part(ipart)%nclass,part(ipart)%idt,part(ipart)%tend, &
      part(ipart)%tstart,part(ipart)%alive,part(ipart)%turbvel%u, &
      part(ipart)%turbvel%v,part(ipart)%turbvel%w,part(ipart)%mesovel%u, &
      part(ipart)%mesovel%v,part(ipart)%mesovel%w,(mass(ipart,j),j=1,nspec), &
      (mass_init(ipart,j),j=1,nspec)
    part(ipart)%spawned = .true.
    if (wetdep) read(unitpartin) (wetdeposit(ipart,j),j=1,nspec)
    if (drydep) read(unitpartin) (drydeposit(ipart,j),j=1,nspec)
    if ((drybkdep).or.(wetbkdep)) read(unitpartin) (xscav_frac1(ipart,j),j=1,nspec)

#ifdef ETA
    part(ipart)%etaupdate=.true.
    part(ipart)%meterupdate=.true.
#endif
  end do

  if (iout.gt.0) then 
#ifdef USE_NCF
    if (lnetcdfout.eq.1) read(unitpartin) tpointer
#endif
    do ks=1,nspec
      do kp=1,maxpointspec_act
        do nage=1,nageclass
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              do l=1,nclassunc
                do kz=1,numzgrid
                  read(unitpartin) gridunc(ix,jy,kz,ks,kp,l,nage)
                end do
                if ((wetdep).and.(ldirect.gt.0)) then
                  read(unitpartin) wetgridunc(ix,jy,ks,kp,l,nage)
                endif
                if ((drydep).and.(ldirect.gt.0)) then
                  read(unitpartin) drygridunc(ix,jy,ks,kp,l,nage)
                endif
              end do
            end do
          end do
          if (nested_output.eq.1) then
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                do l=1,nclassunc
                  do kz=1,numzgrid
                    read(unitpartin) griduncn(ix,jy,kz,ks,kp,l,nage)
                  end do
                  if ((wetdep).and.(ldirect.gt.0)) then
                    read(unitpartin) wetgriduncn(ix,jy,ks,kp,l,nage)
                  endif
                  if ((drydep).and.(ldirect.gt.0)) then
                    read(unitpartin) drygriduncn(ix,jy,ks,kp,l,nage)
                  endif
                end do
              end do
            end do
          endif
          if (iflux.eq.1) then
            do kz=1,numzgrid
              do jy=0,numygridn-1
                do ix=0,numxgridn-1
                  do i=1,5
                    read(unitpartin) flux(i,ix,jy,kz,ks,kp,nage)
                  end do
                end do
              end do
            end do        
          endif
        end do
        if (linit_cond.gt.0) then
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                read(unitpartin) init_cond(ix,jy,kz,ks,kp)
              end do 
            end do
          end do 
        endif
      end do
      if (numreceptor.gt.0) then 
        do n=1,numreceptor
          read(unitpartin) creceptor(n,ks)
        end do
      endif
    end do
  endif
  close(unitpartin)

  iterminate=0
  do i=1,imax
    if ((part(i)%spawned .eqv. .true.) .and. (.not. part(i)%alive)) then
      if (part(i)%tstart.le.itime_init) then
        call terminate_particle(i,part(i)%tend)
        iterminate=iterminate+1
      endif
    endif
  end do

  call rewrite_ialive()
  !count%spawned=count%spawned-iterminate
  numpart=count%spawned
  
  julin=juldate(ibdate,ibtime)+real(itime_init,kind=dp)/86400._dp
  if (abs(julin-bdate).le.1.e-5) then
    write(*,*) ' #### FLEXPART ERROR: PLEASE KEEP IBDATE     #### '
    write(*,*) ' #### AND IBTIME INTACT FROM THE INITIAL RUN!#### '
    error stop
  endif
  call caldate(julin,id1,it1)
  call caldate(bdate,id2,it2)
  write(*,*) ' #### Restarting Flexpart from restart.bin.    #### '
  write(*,*) ' #### Original run started on                  #### '
  write(*,*) 'bdate: ',bdate,id2,it2
  write(*,*) ' #### Restarting run starts on                 #### '
  write(*,*) 'julin: ',julin,id1,it1

  return

9989   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE             #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'restart.bin'//'    #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS        #### '
  write(*,*) ' #### NAME DOES NOT EXISTS, RENAME THE APPROPRIATE #### '
  write(*,*) ' #### RESTART FILE TO restart.bin.                 #### '
end subroutine readrestart

end module restart_mod
