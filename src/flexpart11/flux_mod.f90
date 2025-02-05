! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module flux_mod

  ! flux eastward, westward, northward, southward, upward and downward
  ! fluxes of all species and all ageclasses
  ! areaeast,areanorth [m2] side areas of each grid cell
  use outgrid_mod
  use par_mod
  use com_mod
  use windfields_mod

  implicit none

  !Moved to outgrid_mod, because of dependencies
  ! real,allocatable, dimension (:,:,:,:,:,:,:) :: flux 

  !1 fluxw west - east
  !2 fluxe east - west
  !3 fluxs south - north
  !4 fluxn north - south
  !5 fluxu upward
  !6 fluxd downward
  !real,allocatable, dimension (:,:,:) :: areanorth
  !real,allocatable, dimension (:,:,:) :: areaeast

contains

subroutine calcfluxes(itime,nage,jpart,xold,yold,zold,thread)
  !                       i     i    i    i    i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the gross fluxes across horizontal, eastward and        *
  !     northward facing surfaces. The routine calculates the mass flux        *
  !     due to the motion of only one particle. The fluxes of subsequent calls *
  !     to this subroutine are accumulated until the next output is due.       *
  !     Upon output, flux fields are re-set to zero in subroutine fluxoutput.f.*
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     04 April 2000                                                          *
  !                                                                            *
  !     Changes                                                                *
  !        2021 L. Bakels: OpenMP parallelisation                              *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nage                  Age class of the particle considered                 *
  ! jpart                 Index of the particle considered                     *
  ! xold,yold,zold        "Memorized" old positions of the particle            *
  !                                                                            *
  !*****************************************************************************
  
  use particle_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif

  implicit none
  integer, intent(in) :: thread ! for OMP, number of thread
  integer :: itime,jpart,nage,ixave,jyave,kz,kzave,kp
  integer :: k,k1,k2,ix,ix1,ix2,ixs,jy,jy1,jy2
  real :: xold,yold,zold,xmean,ymean


  ! Determine average positions
  !****************************

  if ((ioutputforeachrelease.eq.1).and.(mdomainfill.eq.0)) then
     kp=part(jpart)%npoint
  else
     kp=1
  endif

#ifdef ETA
  call update_zeta_to_z(itime,jpart)
#endif
  xmean=(xold+real(part(jpart)%xlon))*0.5
  ymean=(yold+real(part(jpart)%ylat))*0.5

  ixave=int((xmean*dx+xoutshift)/dxout)
  jyave=int((ymean*dy+youtshift)/dyout)
  do kz=1,numzgrid                ! determine height of cell
    if (outheight(kz).gt.part(jpart)%z) exit
  end do
  kzave=kz


  ! Determine vertical fluxes
  !**************************

  if ((ixave.ge.0).and.(jyave.ge.0).and.(ixave.le.numxgrid-1).and. &
       (jyave.le.numygrid-1)) then
    do kz=1,numzgrid                ! determine height of cell
      if (outheighthalf(kz).gt.zold) exit
    end do
    k1=min(numzgrid,kz)
    do kz=1,numzgrid                ! determine height of cell
      if (outheighthalf(kz).gt.part(jpart)%z) exit
    end do
    k2=min(numzgrid,kz)

    do k=1,nspec
      do kz=k1,k2-1
#ifdef _OPENMP
        flux_omp(5,ixave,jyave,kz,k,kp,nage,thread)= &
             flux_omp(5,ixave,jyave,kz,k,kp,nage,thread)+ &
             mass(jpart,k)
#else
        flux(5,ixave,jyave,kz,k,kp,nage)= &
             flux(5,ixave,jyave,kz,k,kp,nage)+ &
             mass(jpart,k)
#endif
      end do
      do kz=k2,k1-1
#ifdef _OPENMP
        flux_omp(6,ixave,jyave,kz,k,kp,nage,thread)= &
             flux_omp(6,ixave,jyave,kz,k,kp,nage,thread)+ &
             mass(jpart,k)
#else
        flux(6,ixave,jyave,kz,k,kp,nage)= &
             flux(6,ixave,jyave,kz,k,kp,nage)+ &
             mass(jpart,k)
#endif
      end do
    end do
  endif


  ! Determine west-east fluxes (fluxw) and east-west fluxes (fluxe)
  !****************************************************************

  if ((kzave.le.numzgrid).and.(jyave.ge.0).and. &
       (jyave.le.numygrid-1)) then

  ! 1) Particle does not cross domain boundary

    if (abs(xold-part(jpart)%xlon).lt.real(nx)*0.5) then
      ix1=int((xold*dx+xoutshift)/dxout+0.5)
      ix2=int((part(jpart)%xlon*dx+xoutshift)/dxout+0.5)
      do k=1,nspec
        do ix=ix1,ix2-1
          if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
#ifdef _OPENMP
            flux_omp(1,ix,jyave,kzave,k,kp,nage,thread)= &
                 flux_omp(1,ix,jyave,kzave,k,kp,nage,thread) &
                 +mass(jpart,k)
#else
            flux(1,ix,jyave,kzave,k,kp,nage)= &
                 flux(1,ix,jyave,kzave,k,kp,nage) &
                 +mass(jpart,k)
#endif
          endif
        end do
        do ix=ix2,ix1-1
          if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
#ifdef _OPENMP
            flux_omp(2,ix,jyave,kzave,k,kp,nage,thread)= &
                 flux_omp(2,ix,jyave,kzave,k,kp,nage,thread) &
                 +mass(jpart,k)
#else
            flux(2,ix,jyave,kzave,k,kp,nage)= &
                 flux(2,ix,jyave,kzave,k,kp,nage) &
                 +mass(jpart,k)
#endif
          endif
        end do
      end do

  ! 2) Particle crosses domain boundary: use cyclic boundary condition
  !    and attribute flux to easternmost grid row only (approximation valid
  !    for relatively slow motions compared to output grid cell size)

    else
      ixs=int(((real(nxmin1)-1.e5)*dx+xoutshift)/dxout)
      if ((ixs.ge.0).and.(ixs.le.numxgrid-1)) then
        if (xold.gt.part(jpart)%xlon) then       ! west-east flux
          do k=1,nspec
#ifdef _OPENMP
            flux_omp(1,ixs,jyave,kzave,k,kp,nage,thread)= &
                 flux_omp(1,ixs,jyave,kzave,k,kp,nage,thread) &
                 +mass(jpart,k)
#else
            flux(1,ixs,jyave,kzave,k,kp,nage)= &
                 flux(1,ixs,jyave,kzave,k,kp,nage) &
                 +mass(jpart,k)
#endif
          end do
        else                                 ! east-west flux
          do k=1,nspec
#ifdef _OPENMP
            flux_omp(2,ixs,jyave,kzave,k,kp,nage,thread)= &
                 flux_omp(2,ixs,jyave,kzave,k,kp,nage,thread) &
                 +mass(jpart,k)
#else
            flux(2,ixs,jyave,kzave,k,kp,nage)= &
                 flux(2,ixs,jyave,kzave,k,kp,nage) &
                 +mass(jpart,k)
#endif
          end do
        endif
      endif
    endif
  endif


  ! Determine south-north fluxes (fluxs) and north-south fluxes (fluxn)
  !********************************************************************

  if ((kzave.le.numzgrid).and.(ixave.ge.0).and. &
       (ixave.le.numxgrid-1)) then
    jy1=int((yold*dy+youtshift)/dyout+0.5)
    jy2=int((part(jpart)%ylat*dy+youtshift)/dyout+0.5)

    do k=1,nspec
      do jy=jy1,jy2-1
        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
#ifdef _OPENMP
          flux_omp(3,ixave,jy,kzave,k,kp,nage,thread)= &
               flux_omp(3,ixave,jy,kzave,k,kp,nage,thread) &
               +mass(jpart,k)
#else
          flux(3,ixave,jy,kzave,k,kp,nage)= &
               flux(3,ixave,jy,kzave,k,kp,nage) &
               +mass(jpart,k)
#endif
        endif
      end do
      do jy=jy2,jy1-1
        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
#ifdef _OPENMP
          flux_omp(4,ixave,jy,kzave,k,kp,nage,thread)= &
               flux_omp(4,ixave,jy,kzave,k,kp,nage,thread) &
               +mass(jpart,k)
#else
          flux(4,ixave,jy,kzave,k,kp,nage)= &
               flux(4,ixave,jy,kzave,k,kp,nage) &
               +mass(jpart,k)
#endif
        endif
      end do
    end do
  endif
end subroutine calcfluxes

subroutine fluxoutput(itime)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Output of the gridded fluxes.                                          *
  !     Eastward, westward, northward, southward, upward and downward gross    *
  !     fluxes are written to output file in either sparse matrix or grid dump *
  !     format, whichever is more efficient.                                   *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     04 April 2000                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncellse         number of cells with non-zero values for eastward fluxes   *
  ! sparsee         .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************
  use date_mod
  
  implicit none

  real(kind=dp) :: jul
  integer :: itime,ix,jy,kz,k,nage,jjjjmmdd,ihmmss,kp,i
  integer :: ncellse(maxspec,nageclass),ncellsw(maxspec,nageclass)
  integer :: ncellss(maxspec,nageclass),ncellsn(maxspec,nageclass)
  integer :: ncellsu(maxspec,nageclass),ncellsd(maxspec,nageclass)
  logical :: sparsee(maxspec,nageclass),sparsew(maxspec,nageclass)
  logical :: sparses(maxspec,nageclass),sparsen(maxspec,nageclass)
  logical :: sparseu(maxspec,nageclass),sparsed(maxspec,nageclass)
  character :: adate*8,atime*6


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  open(unitflux,file=path(2)(1:length(2))//'grid_flux_'//adate// &
       atime,form='unformatted')

  !**************************************************************
  ! Check, whether output of full grid or sparse matrix format is
  ! more efficient in terms of storage space. This is checked for
  ! every species and for every age class
  !**************************************************************

  do k=1,nspec
    do nage=1,nageclass
      ncellse(k,nage)=0
      ncellsw(k,nage)=0
      ncellsn(k,nage)=0
      ncellss(k,nage)=0
      ncellsu(k,nage)=0
      ncellsd(k,nage)=0
    end do
  end do

  do k=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do kz=1,numzgrid
            if (flux(2,ix,jy,kz,k,kp,nage).gt.0) ncellse(k,nage)= &
                 ncellse(k,nage)+1
            if (flux(1,ix,jy,kz,k,kp,nage).gt.0) ncellsw(k,nage)= &
                 ncellsw(k,nage)+1
            if (flux(4,ix,jy,kz,k,kp,nage).gt.0) ncellsn(k,nage)= &
                 ncellsn(k,nage)+1
            if (flux(3,ix,jy,kz,k,kp,nage).gt.0) ncellss(k,nage)= &
                 ncellss(k,nage)+1
            if (flux(5,ix,jy,kz,k,kp,nage).gt.0) ncellsu(k,nage)= &
                 ncellsu(k,nage)+1
            if (flux(6,ix,jy,kz,k,kp,nage).gt.0) ncellsd(k,nage)= &
                 ncellsd(k,nage)+1
          end do
        end do
      end do
    end do
  end do
  end do

  ! Output in sparse matrix format more efficient, if less than
  ! 2/5 of all cells contains concentrations>0
  !************************************************************

  do k=1,nspec
    do nage=1,nageclass
      if (4*ncellse(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsee(k,nage)=.true.
      else
        sparsee(k,nage)=.false.
      endif
      if (4*ncellsw(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsew(k,nage)=.true.
      else
        sparsew(k,nage)=.false.
      endif
      if (4*ncellsn(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsen(k,nage)=.true.
      else
        sparsen(k,nage)=.false.
      endif
      if (4*ncellss(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparses(k,nage)=.true.
      else
        sparses(k,nage)=.false.
      endif
      if (4*ncellsu(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparseu(k,nage)=.true.
      else
        sparseu(k,nage)=.false.
      endif
      if (4*ncellsd(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsed(k,nage)=.true.
      else
        sparsed(k,nage)=.false.
      endif
    end do
  end do



  ! Flux output: divide by area and time to get flux in ng/m2/s
  !************************************************************

  write(unitflux) itime
  do k=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass

      if (sparsee(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(2,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(2,ix,jy,kz,k,kp,nage)/areaeast(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(2,ix,jy,kz,k,kp,nage)/ &
                 areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsew(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(1,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(1,ix,jy,kz,k,kp,nage)/areaeast(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(1,ix,jy,kz,k,kp,nage)/ &
                 areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparses(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(3,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(3,ix,jy,kz,k,kp,nage)/areanorth(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(3,ix,jy,kz,k,kp,nage)/ &
                 areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsen(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1 ! north
              if (flux(4,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(4,ix,jy,kz,k,kp,nage)/areanorth(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(4,ix,jy,kz,k,kp,nage)/ &
                 areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparseu(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(5,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(5,ix,jy,kz,k,kp,nage)/area(ix,jy)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(5,ix,jy,kz,k,kp,nage)/ &
                 area(ix,jy)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsed(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(6,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(6,ix,jy,kz,k,kp,nage)/area(ix,jy)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(6,ix,jy,kz,k,kp,nage)/ &
                 area(ix,jy)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

    end do
  end do
  end do


  close(unitflux)

  write(*,*) 'Flux:', itime, flux(:,1,1,1,1,1,1)
  ! Reinitialization of grid
  !*************************

  do k=1,nspec
  do kp=1,maxpointspec_act
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
          do kz=1,numzgrid
            do nage=1,nageclass
              do i=1,6
                flux(i,ix,jy,kz,k,kp,nage)=0.
#ifdef _OPENMP
                flux_omp(i,ix,jy,kz,k,kp,nage,:)=0.
#endif
              end do
            end do
          end do
      end do
    end do
  end do
  end do
end subroutine fluxoutput

end module flux_mod
