!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine convmix(itime,metdata_format)
    !                     i
    !**************************************************************
    !handles all the calculations related to convective mixing
    !Petra Seibert, Bernd C. Krueger, Feb 2001
    !nested grids included, Bernd C. Krueger, May 2001
    !
    !Changes by Caroline Forster, April 2004 - February 2005:
    !  convmix called every lsynctime seconds
    !CHANGES by A. Stohl:
    !  various run-time optimizations - February 2005
    ! CHANGES by C. Forster, November 2005, NCEP GFS version
    !      in the ECMWF version convection is calculated on the
    !      original eta-levels
    !      in the GFS version convection is calculated on the
    !      FLEXPART levels
    !
    !   Unified ECMWF and GFS builds                                             
    !   Marian Harustak, 12.5.2017                                              
    !     - Merged convmix and convmix_gfs into one routine using if-then           
    !       for meteo-type dependent code                                        
    !**************************************************************

    use flux_mod
    use par_mod
    use com_mod
    use conv_mod
    use class_gribfile

    use particles_mod, only : particles, pp
    use settings,      only : config

    implicit none

    integer :: igr,igrold, ipart, itime, ix, j, inest
    integer :: ipconv
    integer :: jy, kpart, ktop, ngrid,kz
    integer :: igrid(config%maxpart), ipoint(config%maxpart), igridn(config%maxpart, maxnests)
    integer :: metdata_format

    ! itime [s]                 current time
    ! igrid(maxpart)            horizontal grid position of each particle
    ! igridn(maxpart,maxnests)  dto. for nested grids
    ! ipoint(maxpart)           pointer to access particles according to grid position

    logical :: lconv
    real :: x, y, xtn,ytn, ztold, delt
    real :: dt1,dt2,dtt
    integer :: mind1,mind2
    ! dt1,dt2,dtt,mind1,mind2       variables used for time interpolation
    integer :: itage,nage
    real,parameter :: eps=nxmax/3.e5


    !monitoring variables
    !real sumconv,sumall


    ! Calculate auxiliary variables for time interpolation
    !*****************************************************

    dt1=real(itime-memtime(1))
    dt2=real(memtime(2)-itime)
    dtt=1./(dt1+dt2)
    mind1=memind(1)
    mind2=memind(2)
    delt=real(abs(lsynctime))


    lconv = .false.

    ! if no particles are present return after initialization
    !********************************************************

    if (numpart<=0) return

    ! Assign igrid and igridn, which are pseudo grid numbers indicating particles
    ! that are outside the part of the grid under consideration
    ! (e.g. particles near the poles or particles in other nests).
    ! Do this for all nests but use only the innermost nest; for all others
    ! igrid shall be -1
    ! Also, initialize index vector ipoint
    !************************************************************************

    do ipart=1,numpart
        igrid(ipart)=-1
        do j=numbnests,1,-1
            igridn(ipart,j)=-1
        end do
        ipoint(ipart)=ipart
        ! do not consider particles that are (yet) not part of simulation

        pp => particles(ipart)
        if (pp%active) then
            x = xtra1(ipart)
            y = ytra1(ipart)

            ! Determine which nesting level to be used
            !**********************************************************

            ngrid=0
            if (metdata_format==GRIBFILE_CENTRE_ECMWF) then
                do j=numbnests,1,-1
                    if ( x>xln(j)+eps .and. x<xrn(j)-eps .and. &
                        y>yln(j)+eps .and. y<yrn(j)-eps ) then
                        ngrid=j
                        goto 23
                    endif
                end do
            else
                do j=numbnests,1,-1
                    if ( x>xln(j) .and. x<xrn(j) .and. &
                        y>yln(j) .and. y<yrn(j) ) then
                        ngrid=j
                        goto 23
                    endif
                end do
            endif
            23   continue

            ! Determine nested grid coordinates
            !**********************************

            if (ngrid>0) then
                ! nested grids
                xtn=(x-xln(ngrid))*xresoln(ngrid)
                ytn=(y-yln(ngrid))*yresoln(ngrid)
                ix=nint(xtn)
                jy=nint(ytn)
                igridn(ipart,ngrid) = 1 + jy*nxn(ngrid) + ix
            else if(ngrid==0) then
                ! mother grid
                ix=nint(x)
                jy=nint(y)
                igrid(ipart) = 1 + jy*nx + ix
            endif

        end if

        nullify(pp)
    end do

    ! sumall = 0.
    ! sumconv = 0.

    !*****************************************************************************
    ! 1. Now, do everything for the mother domain and, later, for all of the nested domains
    ! While all particles have to be considered for redistribution, the Emanuel convection
    ! scheme only needs to be called once for every grid column where particles are present.
    ! Therefore, particles are sorted according to their grid position. Whenever a new grid
    ! cell is encountered by looping through the sorted particles, the convection scheme is called.
    !*****************************************************************************

    ! sort particles according to horizontal position and calculate index vector IPOINT

    call sort2(numpart,igrid,ipoint)

    ! Now visit all grid columns where particles are present
    ! by going through the sorted particles

    igrold = -1
    do kpart=1,numpart
        igr = igrid(kpart)
        if (igr /= -1) then
            ipart = ipoint(kpart)

            pp => particles(ipart)

            ! sumall = sumall + 1
            if (igr /= igrold) then
                ! we are in a new grid column
                jy = (igr-1)/nx
                ix = igr - jy*nx - 1

                ! Interpolate all meteorological data needed for the convection scheme
                psconv=(ps(ix,jy,1,mind1)*dt2+ps(ix,jy,1,mind2)*dt1)*dtt
                tt2conv=(tt2(ix,jy,1,mind1)*dt2+tt2(ix,jy,1,mind2)*dt1)*dtt
                td2conv=(td2(ix,jy,1,mind1)*dt2+td2(ix,jy,1,mind2)*dt1)*dtt
                !!$      do kz=1,nconvlev+1      !old
                if (metdata_format==GRIBFILE_CENTRE_ECMWF) then
                    do kz=1,nuvz-1           !bugfix
                        tconv(kz)=(tth(ix,jy,kz+1,mind1)*dt2+ &
                            tth(ix,jy,kz+1,mind2)*dt1)*dtt
                        qconv(kz)=(qvh(ix,jy,kz+1,mind1)*dt2+ &
                            qvh(ix,jy,kz+1,mind2)*dt1)*dtt
                    end do
                else
                    do kz=1,nuvz-1           !bugfix
                        pconv(kz)=(pplev(ix,jy,kz,mind1)*dt2+ &
                            pplev(ix,jy,kz,mind2)*dt1)*dtt
                        tconv(kz)=(tt(ix,jy,kz,mind1)*dt2+ &
                            tt(ix,jy,kz,mind2)*dt1)*dtt
                        qconv(kz)=(qv(ix,jy,kz,mind1)*dt2+ &
                            qv(ix,jy,kz,mind2)*dt1)*dtt
                    end do
                end if

                ! Calculate translocation matrix
                call calcmatrix(lconv,delt,cbaseflux(ix,jy),metdata_format)
                igrold = igr
                ktop = 0
            endif

            ! treat particle only if column has convection
            if (lconv .eqv. .true.) then
                ! assign new vertical position to particle

                ztold=ztra1(ipart)
                call redist(ipart,ktop,ipconv)
                !    if (ipconv.le.0) sumconv = sumconv+1

                ! Calculate the gross fluxes across layer interfaces
                !***************************************************

                if (iflux==1) then
                    itage = abs(pp%t - itramem(ipart))
                    do nage=1,nageclass
                        if (itage<lage(nage)) goto 37
                    end do
                    37     continue

                    if (nage<=nageclass) &
                        call calcfluxes(nage,ipart,real(xtra1(ipart)), &
                        real(ytra1(ipart)),ztold)
                endif

            endif   !(lconv .eqv. .true)
            nullify(pp)
        end if
    end do


    !*****************************************************************************
    ! 2. Nested domains
    !*****************************************************************************

    ! sort particles according to horizontal position and calculate index vector IPOINT

    do inest=1,numbnests
        do ipart=1,numpart
            ipoint(ipart)=ipart
            igrid(ipart) = igridn(ipart,inest)
        enddo
        call sort2(numpart,igrid,ipoint)

        ! Now visit all grid columns where particles are present
        ! by going through the sorted particles

        igrold = -1
        do kpart=1,numpart
            igr = igrid(kpart)
            if (igr /= -1) then
                ipart = ipoint(kpart)

                pp => particles(ipart)

                ! sumall = sumall + 1
                if (igr /= igrold) then
                    ! we are in a new grid column
                    jy = (igr-1)/nxn(inest)
                    ix = igr - jy*nxn(inest) - 1

                    ! Interpolate all meteorological data needed for the convection scheme
                    psconv=(psn(ix,jy,1,mind1,inest)*dt2+ &
                        psn(ix,jy,1,mind2,inest)*dt1)*dtt
                    tt2conv=(tt2n(ix,jy,1,mind1,inest)*dt2+ &
                        tt2n(ix,jy,1,mind2,inest)*dt1)*dtt
                    td2conv=(td2n(ix,jy,1,mind1,inest)*dt2+ &
                        td2n(ix,jy,1,mind2,inest)*dt1)*dtt
                    !!$        do kz=1,nconvlev+1    !old
                    do kz=1,nuvz-1           !bugfix
                        tconv(kz)=(tthn(ix,jy,kz+1,mind1,inest)*dt2+ &
                            tthn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
                        qconv(kz)=(qvhn(ix,jy,kz+1,mind1,inest)*dt2+ &
                            qvhn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
                    end do

                    ! calculate translocation matrix
                    !*******************************
                    call calcmatrix(lconv,delt,cbasefluxn(ix,jy,inest),metdata_format)
                    igrold = igr
                    ktop = 0
                endif

                ! treat particle only if column has convection
                if (lconv .eqv. .true.) then
                    ! assign new vertical position to particle
                    ztold=ztra1(ipart)
                    call redist(ipart,ktop,ipconv)
                    !      if (ipconv.le.0) sumconv = sumconv+1

                    ! Calculate the gross fluxes across layer interfaces
                    !***************************************************

                    if (iflux==1) then
                        itage = abs(pp%t - itramem(ipart))
                        do nage=1,nageclass
                            if (itage<lage(nage)) goto 47
                        end do
                        47       continue

                        if (nage<=nageclass) &
                            call calcfluxes(nage,ipart,real(xtra1(ipart)), &
                            real(ytra1(ipart)),ztold)
                    endif

                endif !(lconv .eqv. .true.)


            end if
        end do
    end do
    !--------------------------------------------------------------------------
    ! write(*,*)'############################################'
    ! write(*,*)'TIME=',&
    !    &  itime
    ! write(*,*)'fraction of particles under convection',&
    !    &  sumconv/(sumall+0.001)
    ! write(*,*)'total number of particles',&
    !    &  sumall
    ! write(*,*)'number of particles under convection',&
    !    &  sumconv
    ! write(*,*)'############################################'

    return
end subroutine convmix
