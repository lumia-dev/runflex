! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module receptor_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for  calculating         *
  !    receptor concentrations and writing these to file                       *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use point_mod
  use particle_mod
  use date_mod
  use windfields_mod,      only: rho, prs, height, nzmax, nz !, qv
#if USE_NCF
  use receptor_netcdf_mod, only: nc_id, concvar_id, uncvar_id, &
                                recnamevar_id, timevar_id, &
                                reclonvar_id, reclatvar_id, recaltvar_id, &
                                nnvar_id, xkvar_id, rpointer, &
                                ncsat_id, satvar_id, satuncvar_id, &
                                satnamevar_id, sattimevar_id, &
                                satlonvar_id, satlatvar_id, sataltvar_id, &
                                satnnvar_id, satxkvar_id, spointer, sat_name, &
                                receptor_output_netcdf, write_receptor_netcdf, &
                                satellite_output_netcdf, write_satellite_netcdf, &
                                close_receptor_netcdf, close_satellite_netcdf
#endif
  use binary_output_mod,  only: receptorout_init_binary, write_receptor_binary, &
                                satelliteout_init_binary, write_satellite_binary

  implicit none

  contains

 
  subroutine alloc_receptor

  !*****************************************************************************
  !                                                                            *
  !    This routine allocates variables for receptor concentrations            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    if (numreceptor.gt.0) then
      allocate( creceptor(nspec,maxrecsample) )
      allocate( crecuncert(nspec,maxrecsample) )
      allocate( nnreceptor(maxrecsample) )
      allocate( xkreceptor(maxrecsample) )
      creceptor(:,:)=0.
      crecuncert(:,:)=0.
      nnreceptor(:)=0.
      xkreceptor(:)=0.
    endif

  end subroutine alloc_receptor


  subroutine receptoroutput(itime,lrecoutstart,lrecoutend,lrecoutnext,recoutnum,recoutnumsat)

  !*****************************************************************************
  !                                                                            *
  !    Output concentrations at receptors                                      *
  !                                                                            *
  !    Author: R. Thompson, Sep-2023                                           *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

    use omp_lib, only: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM

    implicit none

    integer   :: itime, lrecoutstart, lrecoutend, lrecoutnext
    real, dimension(maxrecsample) :: recoutnum
    real, dimension(nlayermax,maxrecsample) :: recoutnumsat
    real      :: weight
    real(kind=dp) :: jul
    character(len=256) :: fn, fnsat
    character :: adate*8, atime*6
    integer   :: jjjjmmdd, ihmmss
    integer   :: numsatlayer, nchar, ks_start
    integer   :: n, nn, nr, ix, jy, ixp, jyp, indz, indzp, il, ind, ks, k, kz
    real      :: ddx, ddy, rddx, rddy, p1, p2, p3, p4, dz1, dz2, dz
    real      :: rhoi, zmid
    real, dimension(2) :: rho_p
    real, dimension(:), allocatable        :: densityoutrecept
    real, dimension(:,:), allocatable      :: nnrec, xkrec, altrec
    real, dimension(:), allocatable        :: lonrec, latrec
    integer, dimension(:), allocatable     :: timerec
    real, dimension(:,:,:), allocatable    :: crec, cunc
    character(len=16), dimension(:), allocatable :: namerec
    character(len=24), dimension(:), allocatable :: namesatrec
    real, dimension(:,:,:), allocatable    :: nnrec_omp, xkrec_omp, altrec_omp
    real, dimension(:,:), allocatable      :: lonrec_omp, latrec_omp
    integer, dimension(:,:), allocatable   :: timerec_omp
    real, dimension(:,:,:,:), allocatable  :: crec_omp, cunc_omp
    character(len=16), dimension(:,:), allocatable :: namerec_omp
    character(len=24), dimension(:,:), allocatable :: namesatrec_omp
    integer, dimension(:), allocatable     :: nr_omp
    real, parameter :: weightair=28.97
    integer, parameter :: unitoutrecdates=109
    logical :: lexist
    integer :: nthreads, thread, ithread


    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    if (mod(itime-lrecoutstart,lrecoutsample).eq.0) then
      if ((itime.eq.lrecoutstart).or.(itime.eq.lrecoutend)) then
        weight=0.5
      else
        weight=1.0
      endif
      call receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)
    endif

    !! testing
    print*, 'receptor_mod: itime, lrecoutstart, lrecoutend = ',itime, lrecoutstart, lrecoutend

    if ((itime.eq.lrecoutend).and.(any(recoutnum.gt.0.).or.any(recoutnumsat.gt.0.))) then

      ! output receptor concentrations
      !*******************************

      if ((iout.le.3.).or.(iout.eq.5)) then

        ! Determine current date for output in dates_receptor file
        !**********************************************************

        jul=bdate+dble(float(itime))/86400.
        call caldate(jul,jjjjmmdd,ihmmss)
        write(adate,'(i8.8)') jjjjmmdd
        write(atime,'(i6.6)') ihmmss

        inquire(file=path(2)(1:length(2))//'dates_receptors',exist=lexist)
        if (.not.lexist) then
          ! initialize dates output file
          open(unitoutrecdates,file=path(2)(1:length(2))//'dates_receptors',STATUS='REPLACE')
          close(unitoutrecdates)
        endif
        open(unitoutrecdates,file=path(2)(1:length(2))//'dates_receptors', &
              ACCESS='APPEND', STATUS='OLD')
        write(unitoutrecdates,'(a)') adate//atime
        close(unitoutrecdates)

        ! For netcdf output open files and get variable info
        !***************************************************

        if (lnetcdfout.eq.1) then
#ifdef USE_NCF
          if (numreceptor.gt.0) then
            call receptor_output_netcdf()
          endif
          if (numsatreceptor.gt.0) then
            call satellite_output_netcdf()
          endif
#endif
        endif

        !! testing
!        print*, 'receptor_mod: concvar_id = ',concvar_id
!        print*, 'receptor_mod: satvar_id = ',satvar_id


        ! Initialize variables
        !*********************

#if (defined _OPENMP)
        nthreads=OMP_GET_MAX_THREADS()
#else
        nthreads=1
#endif

        !! testing
!        print*, 'receptor_mod: nthread = ',nthreads

        allocate(densityoutrecept(nlayermax))
        allocate(nnrec(maxrecsample,nlayermax),&
                 xkrec(maxrecsample,nlayermax),&
                 lonrec(maxrecsample),&
                 latrec(maxrecsample),&
                 altrec(maxrecsample,nlayermax),&
                 namerec(maxrecsample),&
                 namesatrec(maxrecsample),&
                 timerec(maxrecsample),&
                 crec(nspec,maxrecsample,nlayermax),&
                 cunc(nspec,maxrecsample,nlayermax))
        allocate(nnrec_omp(maxrecsample,nlayermax,nthreads),&
                 xkrec_omp(maxrecsample,nlayermax,nthreads),&
                 lonrec_omp(maxrecsample,nthreads),&
                 latrec_omp(maxrecsample,nthreads),&
                 altrec_omp(maxrecsample,nlayermax,nthreads),&
                 namerec_omp(maxrecsample,nthreads),&
                 namesatrec_omp(maxrecsample,nthreads),&
                 timerec_omp(maxrecsample,nthreads),&
                 crec_omp(nspec,maxrecsample,nlayermax,nthreads),&
                 cunc_omp(nspec,maxrecsample,nlayermax,nthreads))
        allocate(nr_omp(nthreads))
        nnrec(:,:)=0.
        xkrec(:,:)=0.
        crec(:,:,:)=0.
        cunc(:,:,:)=0.
        nnrec_omp(:,:,:)=0.
        xkrec_omp(:,:,:)=0.
        crec_omp(:,:,:,:)=0.
        cunc_omp(:,:,:,:)=0.
        nr_omp(:)=0

        ! Loop over general receptors
        !**************************** 
 
        write(*,fmt='(A,1X,I8,1X,A,1X,I3)') 'Number of receptors output at itime ',itime,'is',numcurrec

!$OMP PARALLEL &
!$OMP PRIVATE(n,k,nn,nr,ix,jy,ixp,jyp,ddx,ddy,rddx,rddy,p1,p2,p3,p4,indz,indzp,il,dz1,dz2,dz, &
!$OMP         ind,rho_p,rhoi,densityoutrecept,ks,thread) 

#if (defined _OPENMP)
        thread=OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
        thread=1
#endif

        nr=0
!$OMP DO
        do k=1,numcurrec

          if (cpointer(k).eq.0.) cycle

          ! number of the receptor
          n=cpointer(k)

          ! counter of receptor values this time interval
          nr=nr+1
          nr_omp(thread)=nr
          
          !! testing
!          print*, 'receptor_mod: n, thread, nr, cpointer(k), rpointer = ',n, thread, nr, cpointer(k), rpointer

          if (((.not.llcmoutput).and.(iout.eq.2)).or.&
              (llcmoutput.and.((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)))) then

            ! Compute air density
            !*********************

            ix=int(xreceptor(n))
            jy=int(yreceptor(n))
            ixp=ix+1
            jyp=jy+1
            ddx=xreceptor(n)-float(ix)
            ddy=yreceptor(n)-float(jy)
            rddx=1.-ddx
            rddy=1.-ddy
            p1=rddx*rddy
            p2=ddx*rddy
            p3=rddx*ddy
            p4=ddx*ddy

            indz=nzmax-1
            indzp=nzmax
            do il=2,nzmax
              if (height(il).gt.zreceptor(n)) then
                indz=il-1
                indzp=il
                exit
              endif
            end do

            dz1=zreceptor(n)-height(indz)
            dz2=height(indzp)-zreceptor(n)
            dz=1./(dz1+dz2)

            ! Take density from 2nd wind field in memory
            !********************************************

            do ind=indz,indzp
              ! assume moist air density
              rho_p(ind-indz+1)=p1*rho(ix ,jy ,ind,2) + &
                                p2*rho(ixp,jy ,ind,2) + &
                                p3*rho(ix ,jyp,ind,2) + &
                                p4*rho(ixp,jyp,ind,2)
              ! dry air density             
!              rho_p(ind-indz+1)= &
!                 p1*rho(ix ,jy ,ind,2)*(1. - qv(ix ,jy ,ind,2)) + &
!                 p2*rho(ixp,jy ,ind,2)*(1. - qv(ixp,jy ,ind,2)) + &
!                 p3*rho(ix ,jyp,ind,2)*(1. - qv(ix ,jyp,ind,2)) + &
!                 p4*rho(ixp,jyp,ind,2)*(1. - qv(ixp,jyp,ind,2))
            end do
            rhoi=(dz1*rho_p(2)+dz2*rho_p(1))*dz

            if (.not.llcmoutput) then
              densityoutrecept(1)=1./rhoi
            else
              densityoutrecept(1)=rhoi
            endif

          else 

            ! no multiplication or division by air density
            densityoutrecept(1)=1.

          endif ! llcmoutput

          ! Write receptor output
          !**********************
          ! llcmoutput = true: creceptor = m_spec/m_air
          ! llcmoutput = false: creceptor = m_spec/V                   

          if (recoutnum(k).gt.0.) then
            do ks = ks_start,nspec
              ! write mass concentration
              if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
                ! concentration (ng/m3)
                crec_omp(ks,nr,1,thread)=creceptor(ks,k)*1.E12*densityoutrecept(1)/recoutnum(k)
                cunc_omp(ks,nr,1,thread)=crecuncert(ks,k)*1.E12*densityoutrecept(1)/recoutnum(k)
              else if (iout.eq.2) then
                ! mixing ratio (ppt)
                ! note: for iout=3 (both conc and mixing ratio) only output conc at receptors
                crec_omp(ks,nr,1,thread)=creceptor(ks,k)*1.E12*densityoutrecept(1)* &
                              weightair/weightmolar(ks)/recoutnum(k)
                cunc_omp(ks,nr,1,thread)=crecuncert(ks,k)*1.E12*densityoutrecept(1)* &
                              weightair/weightmolar(ks)/recoutnum(k)
              endif
            end do
            nnrec_omp(nr,1,thread)=nnreceptor(k)/recoutnum(k)
            xkrec_omp(nr,1,thread)=xkreceptor(k)/recoutnum(k)
          else
            ! no particles found in kernel for this receptor
            crec_omp(:,nr,1,thread)=-999.
            cunc_omp(:,nr,1,thread)=-999.
            nnrec_omp(nr,1,thread)=-999.
            xkrec_omp(nr,1,thread)=-999.
          endif
          lonrec_omp(nr,thread)=xreceptor(n)*dx+xlon0
          latrec_omp(nr,thread)=yreceptor(n)*dy+ylat0
          altrec_omp(nr,1,thread)=zreceptor(n)
          if ( lrecregular ) then
            timerec_omp(nr,thread)=itime
          else
            timerec_omp(nr,thread)=treceptor(n)
          endif
          namerec_omp(nr,thread)=receptorname(n)

        end do ! numcurrec
!$OMP END DO
!$OMP END PARALLEL

        ! concatentate from all threads
        nr=1
        do ithread=1,nthreads
          !! testing
!          print*, 'receptor_mod: nr_omp(ithread) = ',nr_omp(ithread)
          if (nr_omp(ithread).eq.0) cycle
          crec(:,nr:(nr+nr_omp(ithread)-1),1)=crec_omp(:,1:nr_omp(ithread),1,ithread)
          cunc(:,nr:(nr+nr_omp(ithread)-1),1)=cunc_omp(:,1:nr_omp(ithread),1,ithread)
          nnrec(nr:(nr+nr_omp(ithread)-1),1)=nnrec_omp(1:nr_omp(ithread),1,ithread)
          xkrec(nr:(nr+nr_omp(ithread)-1),1)=xkrec_omp(1:nr_omp(ithread),1,ithread)
          lonrec(nr:(nr+nr_omp(ithread)-1))=lonrec_omp(1:nr_omp(ithread),ithread)
          latrec(nr:(nr+nr_omp(ithread)-1))=latrec_omp(1:nr_omp(ithread),ithread)
          altrec(nr:(nr+nr_omp(ithread)-1),1)=altrec_omp(1:nr_omp(ithread),1,ithread)
          timerec(nr:(nr+nr_omp(ithread)-1))=timerec_omp(1:nr_omp(ithread),ithread)
          namerec(nr:(nr+nr_omp(ithread)-1))=namerec_omp(1:nr_omp(ithread),ithread)
          nr=nr+nr_omp(ithread)
          !! testing
!          print*, 'receptor_mod: nr = ',nr
        end do

        ! total number of receptors to write this interval
        nr=nr-1

        ! write receptor output this time interval
        if (nr.gt.0) then
          if (lnetcdfout.eq.1) then
#if USE_NCF
            call write_receptor_netcdf(crec,cunc,nnrec,xkrec,lonrec,latrec,altrec,timerec,namerec,nr)
#endif
          else
            call write_receptor_binary(crec,cunc,nnrec,xkrec,lonrec,latrec,altrec,timerec,namerec,nr)
          endif
        endif
        ! advance output index
        rpointer=rpointer+nr
        
        !! testing
!        print*, 'receptor_mod: nr_omp = ',nr_omp
!        print*, 'receptor_mod: rpointer = ',rpointer

        ! Loop over satellites
        !*********************

        nnrec(:,:)=0.
        xkrec(:,:)=0.
        crec(:,:,:)=0.
        cunc(:,:,:)=0.
        nnrec_omp(:,:,:)=0.
        xkrec_omp(:,:,:)=0.
        crec_omp(:,:,:,:)=0.
        cunc_omp(:,:,:,:)=0.
        nr_omp(:)=0

        write(*,fmt='(A,1X,I8,1X,A,1X,I3)') 'Number of satellite receptors output at itime ',itime,'is',numcursat

!$OMP PARALLEL &
!$OMP PRIVATE(n,nr,nn,nchar,k,ix,jy,ixp,jyp,ddx,ddy,rddx,rddy,p1,p2,p3,p4,kz,numsatlayer,zmid,indz,indzp, &
!$OMP         il,dz1,dz2,dz,ind,rho_p,rhoi,densityoutrecept,ks,thread)

#if (defined _OPENMP)
        thread=OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
        thread=1
#endif

        nr=0
!$OMP DO
        do k=1,numcursat

          if (csatpointer(k).eq.0) cycle

          ! number of satellite receptor
          n=csatpointer(k)

          ! counter of receptor values this time interval
          nr=nr+1
          nr_omp(thread)=nr

          !! testing
!          print*, 'receptor_mod: n, thread, nr, csatpointer(k) = ',n, thread, nr, csatpointer(k)

          ! get actual number vertical layers for this retrieval
          do nn=1,numsatellite
            nchar=len_trim(sat_name(nn))
            if (satellitename(n)(1:nchar).eq.trim(sat_name(nn))) then
              numsatlayer=nnsatlayer(nn)
              exit
            endif
          end do

          if (.not.llcmoutput) then

            ! Compute air density
            !*********************

            ix=int(xsatellite(n))
            jy=int(ysatellite(n))
            ixp=ix+1
            jyp=jy+1
            ddx=xsatellite(n)-float(ix)
            ddy=ysatellite(n)-float(jy)
            rddx=1.-ddx
            rddy=1.-ddy
            p1=rddx*rddy
            p2=ddx*rddy
            p3=rddx*ddy
            p4=ddx*ddy

            do kz=1,numsatlayer
  
              zmid=0.5*(zsatellite(kz,n)+zsatellite(kz+1,n))
              indz=nzmax-1
              indzp=nzmax
              do il=2,nzmax
                if (height(il).gt.zmid) then
                  indz=il-1
                  indzp=il
                  exit
                endif
              end do

              dz1=zmid-height(indz)
              dz2=height(indzp)-zmid
              dz=1./(dz1+dz2)

              ! Take density from 2nd wind field in memory
              !********************************************

              do ind=indz,indzp
                ! assume moist air density
                rho_p(ind-indz+1)=p1*rho(ix ,jy ,ind,2) + &
                                  p2*rho(ixp,jy ,ind,2) + &
                                  p3*rho(ix ,jyp,ind,2) + &
                                  p4*rho(ixp,jyp,ind,2)
                ! dry air density             
!                rho_p(ind-indz+1)= &
!                   p1*rho(ix ,jy ,ind,2)*(1. - qv(ix ,jy ,ind,2)) + &
!                   p2*rho(ixp,jy ,ind,2)*(1. - qv(ixp,jy ,ind,2)) + &
!                   p3*rho(ix ,jyp,ind,2)*(1. - qv(ix ,jyp,ind,2)) + &
!                   p4*rho(ixp,jyp,ind,2)*(1. - qv(ixp,jyp,ind,2))
              end do
              rhoi=(dz1*rho_p(2)+dz2*rho_p(1))*dz
              densityoutrecept(kz)=1./rhoi
 
            end do ! numsatlayer  

          else

            ! no normalization by density
            densityoutrecept(:)=1.

          endif ! llcmoutput

          ! Write receptor output
          !**********************

          do kz=1,numsatlayer
            if (recoutnumsat(kz,k).gt.0.) then
              do ks = ks_start,nspec
                ! only mixing ratio (ppt) for satellites
                crec_omp(ks,nr,kz,thread)=csatellite(ks,kz,k)*1.E12*densityoutrecept(kz)* &
                              weightair/weightmolar(ks)/recoutnumsat(kz,k)
                cunc_omp(ks,nr,kz,thread)=csatuncert(ks,kz,k)*1.E12*densityoutrecept(kz)* &
                              weightair/weightmolar(ks)/recoutnumsat(kz,k)
              end do
              nnrec_omp(nr,kz,thread)=nnsatellite(kz,k)/recoutnumsat(kz,k)
              xkrec_omp(nr,kz,thread)=xksatellite(kz,k)/recoutnumsat(kz,k)
            else
              crec_omp(:,nr,kz,thread)=-999.
              cunc_omp(:,nr,kz,thread)=-999.
              nnrec_omp(nr,kz,thread)=-999.
              xkrec_omp(nr,kz,thread)=-999.
            endif
          end do
          lonrec_omp(nr,thread)=xsatellite(n)*dx+xlon0
          latrec_omp(nr,thread)=ysatellite(n)*dy+ylat0
          altrec_omp(nr,:,thread)=zsatellite(:,n)
          timerec_omp(nr,thread)=tsatellite(n)
          namesatrec_omp(nr,thread)=satellitename(n)

        end do ! numcursat
!$OMP END DO
!$OMP END PARALLEL

        ! concatentate from all threads
        nr=1
        do ithread=1,nthreads
          !! testing
!          print*, 'receptor_mod: satellites, nr_omp(ithread) = ',nr_omp(ithread)
          if (nr_omp(ithread).eq.0) cycle
          crec(:,nr:nr+nr_omp(ithread)-1,:)=crec_omp(:,1:nr_omp(ithread),:,ithread)
          cunc(:,nr:nr+nr_omp(ithread)-1,:)=cunc_omp(:,1:nr_omp(ithread),:,ithread)
          nnrec(nr:nr+nr_omp(ithread)-1,:)=nnrec_omp(1:nr_omp(ithread),:,ithread)
          xkrec(nr:nr+nr_omp(ithread)-1,:)=xkrec_omp(1:nr_omp(ithread),:,ithread)
          lonrec(nr:nr+nr_omp(ithread)-1)=lonrec_omp(1:nr_omp(ithread),ithread)
          latrec(nr:nr+nr_omp(ithread)-1)=latrec_omp(1:nr_omp(ithread),ithread)
          altrec(nr:nr+nr_omp(ithread)-1,:)=altrec_omp(1:nr_omp(ithread),:,ithread)
          timerec(nr:nr+nr_omp(ithread)-1)=timerec_omp(1:nr_omp(ithread),ithread)
          namesatrec(nr:nr+nr_omp(ithread)-1)=namesatrec_omp(1:nr_omp(ithread),ithread)
          nr=nr+nr_omp(ithread)
          !! testing
!          print*, 'receptor_mod: satellites, nr = ',nr
        end do

        ! total number of receptors to write this interval
        nr=nr-1

        ! write satellite output this time interval
        if (nr.gt.0) then
          if (lnetcdfout.eq.1) then
#if USE_NCF
            call write_satellite_netcdf(crec,cunc,nnrec,xkrec,lonrec,latrec,altrec,timerec,namesatrec,nr)
#endif  
          else
            call write_satellite_binary(crec,cunc,nnrec,xkrec,lonrec,latrec,altrec,timerec,namesatrec,nr)
          endif
        endif
        ! advance output index
        spointer=spointer+nr

        !! testing
        print*, 'receptor_mod: satellites, nr = ',nr
        print*, 'receptor_mod: spointer = ',spointer

        ! close files
        if (lnetcdfout.eq.1) then
#ifdef USE_NCF
          if (numreceptor.gt.0) then
            call close_receptor_netcdf
          endif
          if (numsatreceptor.gt.0) then
            call close_satellite_netcdf
          endif
#endif
        endif

        deallocate(densityoutrecept)
        deallocate(crec,cunc,nnrec,xkrec)
        deallocate(crec_omp,cunc_omp,nnrec_omp,xkrec_omp)
        deallocate(lonrec,latrec,altrec,timerec,namerec,namesatrec)
        deallocate(lonrec_omp,latrec_omp,altrec_omp,timerec_omp,namesatrec_omp)
        deallocate(nr_omp)

        ! End of receptor output
        !***********************

        recoutnum(:)=0.
        recoutnumsat(:,:)=0.

      endif ! output receptor conc

    endif ! receptor output

    if (itime.eq.lrecoutend) then

      ! Reinitialize
      !*************
 
      creceptor(:,:)=0.
      crecuncert(:,:)=0.
      nnreceptor(:)=0.
      xkreceptor(:)=0.
      if (numsatellite.gt.0) then
        csatellite(:,:,:)=0.
        csatuncert(:,:,:)=0.
        nnsatellite(:,:)=0.
        xksatellite(:,:)=0.
      endif

      recoutnum(:)=0.
      recoutnumsat(:,:)=0.

      ! Update output timesteps
      !************************
      lrecoutnext=lrecoutnext+lrecoutstep
      lrecoutstart=lrecoutnext-lrecoutaver/2
      lrecoutend=lrecoutnext+lrecoutaver/2

      if (itime.eq.lrecoutstart) then
        weight=0.5
        call receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)
      endif

    endif 

  end subroutine receptoroutput


  subroutine receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the concentrations at receptor points using the         *
  !     kernel method                                                          *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     Modifications:                                                         * 
  !       Sep-2023, R. Thompson: added option for domain filling mode          *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: itime, itage, lrecoutstart, lrecoutend
    real, dimension(maxrecsample) :: recoutnum
    real, dimension(nlayermax,maxrecsample) :: recoutnumsat
    real    :: weight
    real    :: xd, yd, zd, hx, hy, hz, h, r2, xkern
    real, dimension(nspec) :: conc 
    real, dimension(:,:), allocatable :: unc 
    integer :: n, j, k, ks, kz, i, nn, ks_start, nchar, numsatlayer
    real    :: hxmax, hymax, hzmax, hxsat, hysat, rec_ff, xksum, eta, zmid
    real, parameter :: factor=.596831  ! 15/(8*pi)
    real, parameter :: zref=2000.      ! normalizing height for calculating eta
    integer, parameter :: jmax=4000    ! max number of particles in kernel

    ! initialization
    !***************

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif
    allocate(unc(nspec,jmax))

    ! hxmax and hymax in degrees, hzmax in metres
    !********************************************

    if (mdomainfill.ne.1) then
      hxmax=6.0
      hymax=4.0
      hzmax=150.
    else
      hxmax=2.0
      hymax=1.5
      hzmax=300.
      hxsat=1.75
      hysat=1.25
    endif

    ! convert h-values to grid coordinates
    !*************************************

    hxmax=hxmax/dx
    hymax=hymax/dy
    hxsat=hxsat/dx
    hysat=hysat/dy

    ! Loop over receptors
    !********************

    ! pointer for creceptor, xkreceptor and nnreceptor
    numcurrec=0
    cpointer(:)=0
    k=0

    do n=1,numreceptor

      if ((.not.lrecregular).and.((treceptor(n).lt.lrecoutstart).or. &
           (treceptor(n).ge.lrecoutend))) cycle  ! skip if not in current sampling time interval

      ! update pointer
      k=k+1
      if (k.gt.maxrecsample) then
        write(*,*) 'FLEXPART ERROR in receptorcalc: maxrecsample too small'
        error stop
      endif
      cpointer(k)=n

      !! testing
!      print*, 'receptorcalc: n, receptorname(n) = ',n, receptorname(n)

      ! Reset concentrations for new receptor
      conc(:)=0.
      unc(:,:)=0.
      xksum=0.
      j=0

      if (zreceptor(n).lt.hzmax) then
        rec_ff=0.5 + 0.9375*(zreceptor(n)/hzmax) - &
                 0.625*(zreceptor(n)/hzmax)**3 + &
                 0.1875*(zreceptor(n)/hzmax)**5
        rec_ff=1./rec_ff
      else
        rec_ff=1.
      endif
      h=hxmax*hymax*hzmax

      !! testing
!      print*, 'receptorcalc: rec_ff = ',rec_ff

      ! Estimate concentration at receptor
      !***********************************

!$OMP PARALLEL &
!$OMP PRIVATE(i,itage,hz,zd,hx,xd,hy,yd,h,r2,xkern,ks) &
!$OMP SHARED(j,unc,conc) &
!$OMP REDUCTION(+:xksum,xkreceptor,nnreceptor)

!$OMP DO
      do i=1,count%alive

        if (mdomainfill.ne.1) then

          ! not domain filling run so consider age of particle
          itage=abs(itime-part(i)%tstart)

          hz=min(50.+0.3*sqrt(real(itage)),hzmax)
          zd=part(i)%z/hz
          if (zd.gt.1.) cycle           ! save computing time

          hx=min((0.29+2.222e-3*sqrt(real(itage)))*dx+ &
               real(itage)*1.2e-5,hxmax)                     ! 80 km/day
          xd=(part(i)%xlon-xreceptor(n))/hx
          if (xd*xd.gt.1.) cycle        ! save computing time

          hy=min((0.18+1.389e-3*sqrt(real(itage)))*dy+ &
               real(itage)*7.5e-6,hymax)                     ! 80 km/day
          yd=(part(i)%ylat-yreceptor(n))/hy
          if (yd*yd.gt.1.) cycle        ! save computing time
          h=hx*hy*hz

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time
          xkern=2.*factor*(1.-r2)       ! parabolic kernel

        else

          ! domain filling run 
          xd=(part(i)%xlon-xreceptor(n))/hxmax
          if (xd*xd.gt.1) cycle         ! save computing time

          yd=(part(i)%ylat-yreceptor(n))/hymax
          if (yd*yd.gt.1.) cycle        ! save computing time

          zd=(part(i)%z-zreceptor(n))/hzmax
          if (zd*zd.gt.1.) cycle        ! save computing time

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time
          xkern=rec_ff*factor*(1.-r2)   ! parabolic kernel

        endif ! mdomainfill

        ! counter of particles used in conc calculation
        ! note: use of atomic may not be efficient, alternative split unc, j over threads
        ! and combine all threads after end of parallel section, need to save j from each thread
!$OMP ATOMIC
        j=j+1
        if (j.gt.jmax) then
          write(*,*) 'FLEXPART ERROR in receptorcalc: size of jmax too small'
          error stop
        endif

        do ks=ks_start,nspec
          if (llcmoutput) then
            ! special case LCM output use mass ratio species to airtracer
            ! species 1 is always airtracer
            conc(ks)=conc(ks) + mass(i,ks)/mass(i,1) * &
                      weight * xkern
            unc(ks,j)=mass(i,ks)/mass(i,1)
          else
            ! normal case
            conc(ks)=conc(ks) + mass(i,ks) * &
                       weight * xkern/h/receptorarea(n)
            unc(ks,j)=mass(i,ks)/h/receptorarea(n)
          endif
        end do
        nnreceptor(k)=nnreceptor(k) + 1.
        xkreceptor(k)=xkreceptor(k) + xkern
        xksum=xksum + xkern

      end do ! count%alive
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL IF(nspec>99) PRIVATE(ks)
!$OMP DO
      do ks=ks_start,nspec
        if (conc(ks).gt.0.) then
          ! only do if conc could be calculated for this receptor at this time
          creceptor(ks,k)=creceptor(ks,k) + conc(ks)/xksum
          crecuncert(ks,k)=crecuncert(ks,k) + &
                               sqrt(sum((unc(ks,1:j)-conc(ks)/xksum/weight)**2)/real(j))
        endif
      end do
!$OMP END DO
!$OMP END PARALLEL
      if (any(conc(:).gt.0.)) then
        recoutnum(k)=recoutnum(k)+weight
      endif

      ! update number receptors this time interval
      numcurrec=k

      !! testing
!      print*, 'receptorcalc: j, conc, xksum = ',j, conc(2), xksum
!      print*, 'receptorcalc: n, k, cpointer(k) = ',n, k, cpointer(k)
!      print*, 'receptorcalc: nnreceptor(k) = ',nnreceptor(k)
!      print*, 'receptorcalc: xkreceptor(k) = ',xkreceptor(k)
!      print*, 'receptorcalc: creceptor(2,k) = ',creceptor(2,k)
!      print*, 'receptorcalc: crecuncert(2,k) = ',crecuncert(2,k)

    end do ! numreceptor

    ! Loop over satellites
    !*********************

    ! pointer for csatellite, xksatellite and nnsatellite
    numcursat=0
    csatpointer(:)=0
    k=0

    do n=1,numsatreceptor

      if ((tsatellite(n).lt.lrecoutstart).or. &
           (tsatellite(n).ge.lrecoutend)) cycle  ! skip if not in current sampling time interval
   
       !! testing
!       print*, 'receptorcalc: lrecoutstart, lrecoutend, tsatellite(n) = ',lrecoutstart, lrecoutend, tsatellite(n)

      ! update pointer
      k=k+1
      if (k.gt.maxrecsample) then
        write(*,*) 'FLEXPART ERROR in receptorcalc: maxrecsample too small'
        error stop
      endif
      csatpointer(k)=n

      ! get actual number vertical layers for this retrieval
      do nn=1,numsatellite
        nchar=len_trim(sat_name(nn))
        if (satellitename(n)(1:nchar).eq.trim(sat_name(nn))) then
          numsatlayer=nnsatlayer(nn)
          exit
        endif
      end do

      !! testing
!      print*, 'receptorcalc: n, satellitename(n) = ',n, satellitename(n)

      do kz=1,numsatlayer

        ! Reset concentrations for new receptor
        conc(:)=0.
        unc(:,:)=0.
        xksum=0.
        j=0

        ! height of layer
        hz=zsatellite(kz+1,n)-zsatellite(kz,n)

        ! rare cases when 2 satellite pressure layers fall in same meteo layer
        ! set minimum height between layers 
        hz=max(200.,hz)

        ! midpoint of layer
        zmid=zsatellite(kz,n)+0.5*hz

        ! factor by which to expand horizontal threshhold
        ! as expect particle density to decrease with altitude
        eta=max(1.,sqrt(0.5*zmid/zref))

        h=hxsat*hysat*(0.5*hz)

        !! testing
!        print*, 'zmid, hz, eta =',zmid, hz, eta

!$OMP PARALLEL &
!$OMP PRIVATE(i,zd,xd,yd,r2,xkern,ks) &
!$OMP SHARED(j,unc,conc) &
!$OMP REDUCTION(+:xksum,xksatellite,nnsatellite)

!$OMP DO
        do i=1,count%alive

          ! sample satellite retrievals for domain-filling and 
          ! non-domain-filling modes in the same way

          zd=(part(i)%z-zmid)/(0.5*hz)
          if (zd*zd.gt.1.) cycle        ! save computing time

          xd=(part(i)%xlon-xsatellite(n))/(eta*hxsat)
          if (xd*xd.gt.1) cycle         ! save computing time

          yd=(part(i)%ylat-ysatellite(n))/(eta*hysat)
          if (yd*yd.gt.1.) cycle        ! save computing time

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time

          xkern=factor*(1.-r2)          ! parabolic kernel

          ! counter of particles used in conc calculation
!$OMP ATOMIC
          j=j+1
          if (j.gt.jmax) then
            write(*,*) 'FLEXPART ERROR in receptorcalc: size of jmax too small'
            error stop
          endif

          do ks=ks_start,nspec
            if (llcmoutput) then
              ! special case LCM output use mass ratio species to airtracer
              ! species 1 is always airtracer
              conc(ks)=conc(ks) + mass(i,ks)/mass(i,1) * &
                        weight * xkern
              unc(ks,j)=mass(i,ks)/mass(i,1)
            else
              ! normal case
              conc(ks)=conc(ks) + mass(i,ks) * &
                         weight * xkern/h/satellitearea(n)
              unc(ks,j)=mass(i,ks)/h/satellitearea(n)
            endif
          end do
          nnsatellite(kz,k)=nnsatellite(kz,k) + 1.
          xksatellite(kz,k)=xksatellite(kz,k) + xkern
          xksum=xksum + xkern

        end do ! count%alive
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL IF(nspec>99) PRIVATE(ks)
!$OMP DO
        do ks=ks_start,nspec
          if (conc(ks).gt.0.) then
            ! only do if conc could be calculated for this receptor and time
            csatellite(ks,kz,k)=csatellite(ks,kz,k) + conc(ks)/xksum
            csatuncert(ks,kz,k)=csatuncert(ks,kz,k) + &
                                 sqrt(sum((unc(ks,1:j)-conc(ks)/xksum/weight)**2)/real(j))
          endif
        end do
!$OMP END DO
!$OMP END PARALLEL
        if (any(conc(:).gt.0.)) then
          recoutnumsat(kz,k)=recoutnumsat(kz,k)+weight
        endif

        !! testing
!        print*, 'receptorcalc: for satellite: j, conc, xksum = ',j, conc(2), xksum
!        print*, 'receptorcalc: n, k, csatpointer(k) = ',n, k, csatpointer(k)
!        print*, 'receptorcalc: nnsatellite(:,k) = ',nnsatellite(:,k)
!        print*, 'receptorcalc: xksatellite(:,k) = ',xksatellite(:,k)

      end do ! numsatlayer

      ! update current number of satellite receptors this time interval
      numcursat=k

    end do ! numsatreceptor

  end subroutine receptorcalc


end module receptor_mod

