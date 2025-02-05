! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  ! L. Bakels 2022, This module contains routines that are relevant for the    *
  !                 plume trajectory computations                              *
  !                                                                            *
  !*****************************************************************************
module plume_mod
  
  implicit none
  private :: centerofmass,clustering,distance,distance2

  public :: plumetraj,openouttraj
contains

subroutine plumetraj(itime)
  !                       i
  !*****************************************************************************
  !                                                                            *
  ! Determines a plume centroid trajectory for each release site, and manages  *
  ! clustering of particle locations. Certain parameters (average PV,          *
  ! tropopause height, etc., are provided along the plume trajectories.        *
  ! At the end, output is written to file 'trajectories.txt'.                  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 January 2002                                                        *
  !                                                                            *
  ! Variables:                                                                 *
  ! fclust          fraction of particles belonging to each cluster            *
  ! hmixcenter      mean mixing height for all particles                       *
  ! ncluster        number of clusters to be used                              *
  ! pvcenter        mean PV for all particles                                  *
  ! pvfract         fraction of particles with PV<2pvu                         *
  ! rms             total horizontal rms distance after clustering             *
  ! rmsdist         total horizontal rms distance before clustering            *
  ! rmsclust        horizontal rms distance for each individual cluster        *
  ! topocenter      mean topography underlying all particles                   *
  ! tropocenter     mean tropopause height at the positions of particles       *
  ! tropofract      fraction of particles within the troposphere               *
  ! zrms            total vertical rms distance after clustering               *
  ! zrmsdist        total vertical rms distance before clustering              *
  ! xclust,yclust,  Cluster centroid positions                                 *
  ! zclust                                                                     *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use mean_mod
  use particle_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use windfields_mod

  implicit none

  integer :: ii,itime,ix,jy,ixp,jyp,indexh,i,j,k,m,n,il,ind,indz,indzp
  ! real :: xl(maxpart),yl(maxpart),zl(maxpart) ! moved to particle_mod and now xplum,yplum,zplum
  real :: xcenter,ycenter,zcenter,dist,rmsdist,zrmsdist

  real :: xclust(ncluster),yclust(ncluster),zclust(ncluster)
  real :: fclust(ncluster),rms,rmsclust(ncluster),zrms

  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,topocenter,hm(2),hmixi,hmixfract,hmixcenter
  real :: pv1(2),pvprof(2),pvi,pvcenter,pvfract,tr(2),tri,tropofract
  real :: tropocenter


  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  ! Loop about all release points
  !******************************

  do j=1,numpoint
    if (abs(ireleasestart(j)-itime).gt.lage(nageclass)) cycle
    topocenter=0.
    hmixcenter=0.
    hmixfract=0.
    tropocenter=0.
    tropofract=0.
    pvfract=0.
    pvcenter=0.
    rmsdist=0.
    zrmsdist=0.

    n=0
    do ii=1,count%alive
      i=count%ialive(ii)
      if (part(i)%npoint.ne.j) cycle
      n=n+1
      xplum(n)=xlon0+real(part(i)%xlon)*dx
      yplum(n)=ylat0+real(part(i)%ylat)*dy
#ifdef ETA
      call update_zeta_to_z(itime,i)
#endif
      zplum(n)=real(part(i)%z)

  ! Interpolate PBL height, PV, and tropopause height to each
  ! particle position in order to determine fraction of particles
  ! within the PBL, above tropopause height, and average PV.
  ! Interpolate topography, too, and convert to altitude asl
  !**************************************************************

      ix=int(part(i)%xlon)
      jy=int(part(i)%ylat)
      ixp=ix+1
      jyp=jy+1

      ! eso: Temporary fix for particle exactly at north pole
      if (jyp >= nymax) then
        write(*,*) 'WARNING: plume_mod.f90 jyp >= nymax. xt,yt:',part(i)%xlon,part(i)%ylat
        jyp=jyp-1
      end if

      if (ixp >= nxmax) then
        write(*,*) 'WARNING: plume_mod.f90 ixp >= nxmax. xt,yt:',part(i)%xlon,part(i)%ylat
        ixp=ixp-nxmax
      end if

      ddx=real(part(i)%xlon)-real(ix)
      ddy=real(part(i)%ylat)-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

  ! Topography
  !***********

      topo=p1*oro(ix ,jy) &
           + p2*oro(ixp,jy) &
           + p3*oro(ix ,jyp) &
           + p4*oro(ixp,jyp)
      topocenter=topocenter+topo

  ! Potential vorticity
  !********************
      indz=nz-1
      indzp=nz
      do il=2,nz
        if (height(il).gt.zplum(n)) then
          indz=il-1
          indzp=il
          exit
        endif
      end do

      dz1=zplum(n)-height(indz)
      dz2=height(indzp)-zplum(n)
      dz=1./(dz1+dz2)


      do ind=indz,indzp
        do m=1,2
          indexh=memind(m)
          pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
               +p2*pv(ixp,jy ,ind,indexh) &
               +p3*pv(ix ,jyp,ind,indexh) &
               +p4*pv(ixp,jyp,ind,indexh)
        end do
        pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
      end do
      pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
      pvcenter=pvcenter+pvi
      if (yplum(n).gt.0.) then
        if (pvi.lt.2.) pvfract=pvfract+1.
      else
        if (pvi.gt.-2.) pvfract=pvfract+1.
      endif


  ! Tropopause and PBL height
  !**************************

      do m=1,2
        indexh=memind(m)

        tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
             + p2*tropopause(ixp,jy ,1,indexh) &
             + p3*tropopause(ix ,jyp,1,indexh) &
             + p4*tropopause(ixp,jyp,1,indexh)

        hm(m)=p1*hmix(ix ,jy ,1,indexh) &
             + p2*hmix(ixp,jy ,1,indexh) &
             + p3*hmix(ix ,jyp,1,indexh) &
             + p4*hmix(ixp,jyp,1,indexh)
      end do

      hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
      tri=(tr(1)*dt2+tr(2)*dt1)*dtt
      if (zplum(n).lt.tri) tropofract=tropofract+1.
      tropocenter=tropocenter+tri+topo
      if (zplum(n).lt.hmixi) hmixfract=hmixfract+1.
      zplum(n)=zplum(n)+topo        ! convert to height asl
      hmixcenter=hmixcenter+hmixi

    end do


  ! Make statistics for all plumes with n>0 particles
  !**************************************************

    if (n.gt.0) then
      topocenter=topocenter/real(n)
      hmixcenter=hmixcenter/real(n)
      pvcenter=pvcenter/real(n)
      tropocenter=tropocenter/real(n)
      hmixfract=100.*hmixfract/real(n)
      pvfract=100.*pvfract/real(n)
      tropofract=100.*tropofract/real(n)

  ! Cluster the particle positions
  !*******************************

      call clustering(n,xclust,yclust,zclust,fclust,rms, &
           rmsclust,zrms)


  ! Determine center of mass position on earth and average height
  !**************************************************************

      call centerofmass(xplum,yplum,n,xcenter,ycenter)
      call mean(zplum,zcenter,zrmsdist,n)

  ! Root mean square distance from center of mass
  !**********************************************

      do k=1,n
        dist=distance(yplum(k),xplum(k),ycenter,xcenter)
        rmsdist=rmsdist+dist*dist
      end do
      if (rmsdist.gt.0.) rmsdist=sqrt(rmsdist/real(n))
      rmsdist=max(rmsdist,0.)

  ! Write out results in trajectory data file
  !******************************************

      write(unitouttraj,'(i5,i8,2f9.4,4f8.1,f8.2,4f8.1,3f6.1,&
           &5(2f8.3,f7.0,f6.1,f8.1))')&
           &j,itime-(ireleasestart(j)+ireleaseend(j))/2, &
           xcenter,ycenter,zcenter,topocenter,hmixcenter,tropocenter, &
           pvcenter,rmsdist,rms,zrmsdist,zrms,hmixfract,pvfract, &
           tropofract, &
           (xclust(k),yclust(k),zclust(k),fclust(k),rmsclust(k), &
           k=1,ncluster)
    endif

  end do
end subroutine plumetraj

subroutine centerofmass(xl,yl,n,xcenter,ycenter)
  !                        i  i  i    o       o
  !*****************************************************************************
  !                                                                            *
  !   This routine calculates the center of mass of n points on the Earth.     *
  !   Input are the longitudes (xl) and latitudes (yl) of the individual       *
  !   points, output is the longitude and latitude of the centre of mass.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 January 2002                                                        *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  integer :: n,l
  real :: xl(n),yl(n),xll,yll,xav,yav,zav,x,y,z,xcenter,ycenter


  xav=0.
  yav=0.
  zav=0.

  do l=1,n

  ! Convert longitude and latitude from degrees to radians
  !*******************************************************

    xll=xl(l)*pi180
    yll=yl(l)*pi180

  ! Calculate 3D coordinates from longitude and latitude
  !*****************************************************

    x = cos(yll)*sin(xll)
    y = -1.*cos(yll)*cos(xll)
    z = sin(yll)


  ! Find the mean location in Cartesian coordinates
  !************************************************

    xav=xav+x
    yav=yav+y
    zav=zav+z
  end do

  xav=xav/real(n)
  yav=yav/real(n)
  zav=zav/real(n)


  ! Project the point back onto Earth's surface
  !********************************************

  xcenter=atan2(xav,-1.*yav)
  ycenter=atan2(zav,sqrt(xav*xav+yav*yav))

  ! Convert back to degrees
  !************************

  xcenter=xcenter/pi180
  ycenter=ycenter/pi180
end subroutine centerofmass

subroutine clustering(n,xclust,yclust,zclust,fclust,rms, &
       rmsclust,zrms)
  !                      i  i  i  i   o      o      o      o     o
  !   o      o
  !*****************************************************************************
  !                                                                            *
  !   This routine clusters the particle position into ncluster custers.       *
  !   Input are the longitudes (xl) and latitudes (yl) of the individual       *
  !   points, output are the cluster mean positions (xclust,yclust).           *
  !   Vertical positions are not directly used for the clustering.             *
  !                                                                            *
  !   For clustering, the procedure described in Dorling et al. (1992) is used.*
  !                                                                            *
  !   Dorling, S.R., Davies, T.D. and Pierce, C.E. (1992):                     *
  !   Cluster analysis: a technique for estimating the synoptic meteorological *
  !   controls on air and precipitation chemistry - method and applications.   *
  !   Atmospheric Environment 26A, 2575-2581.                                  *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 February 2002                                                        *
  !                                                                            *
  ! Variables:                                                                 *
  ! fclust          fraction of particles belonging to each cluster            *
  ! ncluster        number of clusters to be used                              *
  ! rms             total horizontal rms distance after clustering             *
  ! rmsclust        horizontal rms distance for each individual cluster        *
  ! zrms            total vertical rms distance after clustering               *
  ! xclust,yclust,  Cluster centroid positions                                 *
  ! zclust                                                                     *
  ! xl,yl,zl        particle positions                                         *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use particle_mod

  implicit none

  integer :: n,i,j,l,numb(ncluster),ncl
  real :: xclust(ncluster),yclust(ncluster),x,y,z
  real :: zclust(ncluster),distances,distancemin,rms,rmsold
  real :: xav(ncluster),yav(ncluster),zav(ncluster),fclust(ncluster)
  real :: rmsclust(ncluster)
  real :: zdist,zrms



  if (n.lt.ncluster) return
  rmsold=-5.

  ! Convert longitude and latitude from degrees to radians
  !*******************************************************

  do i=1,n
    nclust(i)=i
    xplum(i)=xplum(i)*pi180
    yplum(i)=yplum(i)*pi180
  end do


  ! Generate a seed for each cluster
  !*********************************

  do j=1,ncluster
    zclust(j)=0.
    xclust(j)=xplum(j*n/ncluster)
    yclust(j)=yplum(j*n/ncluster)
  end do


  ! Iterative loop to compute the cluster means
  !********************************************

  do l=1,100

  ! Assign each particle to a cluster: criterion minimum distance to the
  ! cluster mean position
  !*********************************************************************


    do i=1,n
      distancemin=10.**10.
      ncl=1
      do j=1,ncluster
        distances=distance2(yplum(i),xplum(i),yclust(j),xclust(j))
        if (distances.lt.distancemin) then
          distancemin=distances
          ncl=j
        endif
      end do
      nclust(i)=ncl
    end do


  ! Recalculate the cluster centroid position: convert to 3D Cartesian coordinates,
  ! calculate mean position, and re-project this point onto the Earth's surface
  !*****************************************************************************

    do j=1,ncluster
      xav(j)=0.
      yav(j)=0.
      zav(j)=0.
      rmsclust(j)=0.
      numb(j)=0
    end do
    rms=0.

    do i=1,n
      numb(nclust(i))=numb(nclust(i))+1
      distances=distance2(yplum(i),xplum(i), &
           yclust(nclust(i)),xclust(nclust(i)))

  ! rms is the total rms of all particles
  ! rmsclust is the rms for a particular cluster
  !*********************************************

      rms=rms+distances*distances
      rmsclust(nclust(i))=rmsclust(nclust(i))+distances*distances

  ! Calculate Cartesian 3D coordinates from longitude and latitude
  !***************************************************************

      x = cos(yplum(i))*sin(xplum(i))
      y = -1.*cos(yplum(i))*cos(xplum(i))
      z = sin(yplum(i))
      xav(nclust(i))=xav(nclust(i))+x
      yav(nclust(i))=yav(nclust(i))+y
      zav(nclust(i))=zav(nclust(i))+z
    end do

    rms=sqrt(rms/real(n))


  ! Find the mean location in Cartesian coordinates
  !************************************************

    do j=1,ncluster
      if (numb(j).gt.0) then
        rmsclust(j)=sqrt(rmsclust(j)/real(numb(j)))
        xav(j)=xav(j)/real(numb(j))
        yav(j)=yav(j)/real(numb(j))
        zav(j)=zav(j)/real(numb(j))

  ! Project the point back onto Earth's surface
  !********************************************

        xclust(j)=atan2(xav(j),-1.*yav(j))
        yclust(j)=atan2(zav(j),sqrt(xav(j)*xav(j)+yav(j)*yav(j)))
      endif
    end do


  ! Leave the loop if the RMS distance decreases only slightly between 2 iterations
  !*****************************************************************************

    if ((l.gt.1).and.(abs(rms-rmsold)/rmsold.lt.0.005)) exit
    rmsold=rms

  end do

  ! Convert longitude and latitude from radians to degrees
  !*******************************************************

  do i=1,n
    xplum(i)=xplum(i)/pi180
    yplum(i)=yplum(i)/pi180
    zclust(nclust(i))=zclust(nclust(i))+zplum(i)
  end do

  do j=1,ncluster
    xclust(j)=xclust(j)/pi180
    yclust(j)=yclust(j)/pi180
    if (numb(j).gt.0) zclust(j)=zclust(j)/real(numb(j))
    fclust(j)=100.*real(numb(j))/real(n)
  end do

  ! Determine total vertical RMS deviation
  !***************************************

  zrms=0.
  do i=1,n
    zdist=zplum(i)-zclust(nclust(i))
    zrms=zrms+zdist*zdist
  end do
  if (zrms.gt.0.) zrms=sqrt(zrms/real(n))

end subroutine clustering

real function distance(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH.
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN DEGREES
  !rlon1    - REAL LONGITUDE OF POINT 1 IN DEGREES
  !rlat2    - REAL LATITUDE OF POINT 2 IN DEGREES
  !rlon2    - REAL LONGITUDE OF POINT 2 IN DEGREES
  !
  !   OUTPUT ARGUMENT LIST:
  !distance - REAL GREAT CIRCLE DISTANCE IN KILOMETERS
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real          :: rlat1,rlon1,rlat2,rlon2
  real(kind=dp) :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp, dpr=180.0_dp/pi
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ((abs(rlat1-rlat2).lt.0.03).and. &
       (abs(rlon1-rlon2).lt.0.03)) then
    distance=0.
  else
    clat1=cos(real(rlat1,kind=dp)/dpr)
    slat1=sin(real(rlat1,kind=dp)/dpr)
    clat2=cos(real(rlat2,kind=dp)/dpr)
    slat2=sin(real(rlat2,kind=dp)/dpr)
    cdlon=cos(real((rlon1-rlon2),kind=dp)/dpr)
    crd=slat1*slat2+clat1*clat2*cdlon
    distance=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance

real function distance2(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH. COORDINATES ARE GIVEN IN RADIANS!
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN RADIANS
  !rlon1    - REAL LONGITUDE OF POINT 1 IN RADIANS
  !rlat2    - REAL LATITUDE OF POINT 2 IN RADIANS
  !rlon2    - REAL LONGITUDE OF POINT 2 IN RADIANS
  !
  !   OUTPUT ARGUMENT LIST:
  !distance2   - REAL GREAT CIRCLE DISTANCE IN KM
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real                    :: rlat1,rlon1,rlat2,rlon2
  real(kind=dp)           :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp

  if ((abs(rlat1-rlat2).lt.0.0003).and. &
       (abs(rlon1-rlon2).lt.0.0003)) then
    distance2=0.0_dp
  else

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    clat1=cos(real(rlat1,kind=dp))
    slat1=sin(real(rlat1,kind=dp))
    clat2=cos(real(rlat2,kind=dp))
    slat2=sin(real(rlat2,kind=dp))
    cdlon=cos(real(rlon1-rlon2,kind=dp))
    crd=slat1*slat2+clat1*clat2*cdlon
    distance2=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance2

subroutine openouttraj

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the output file for the plume trajectory output       *
  !   produced by the cluster analysis.                                        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     27 January 2001                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i
  real :: xp1,yp1,xp2,yp2


  ! Open output file for trajectory output
  !***************************************

  open(unitouttraj,file=path(2)(1:length(2))//'trajectories.txt', &
       form='formatted',err=998)

  if (ldirect.eq.1) then
  write(unitouttraj,'(i8,1x,i6,1x,a)') ibdate,ibtime, trim(flexversion)
  else
  write(unitouttraj,'(i8,1x,i6,1x,a)') iedate,ietime, trim(flexversion)
  endif
  write(unitouttraj,*) method,lsubgrid,lconvection
  write(unitouttraj,*) numpoint
  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitouttraj,*) ireleasestart(i),ireleaseend(i), &
         xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i),kindz(i),npart(i)
    if (numpoint.le.1000) then
      write(unitouttraj,'(a)') compoint(i)(1:40)
    else
      write(unitouttraj,'(a)') compoint(1001)(1:40)
    endif
  end do

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### trajectories.txt                         #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  error stop

end subroutine openouttraj

end module plume_mod