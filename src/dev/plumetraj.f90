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
    use particles_mod, only : particles, pp
    use settings,      only : config

    implicit none

    integer :: itime,ix,jy,ixp,jyp,indexh,i,j,k,m,n,il,ind,indz,indzp
    real :: xl(config%maxpart), yl(config%maxpart), zl(config%maxpart)
    real :: xcenter,ycenter,zcenter,dist,distance,rmsdist,zrmsdist

    real :: xclust(ncluster),yclust(ncluster),zclust(ncluster)
    real :: fclust(ncluster),rms,rmsclust(ncluster),zrms

    real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
    real :: topo,topocenter,hm(2),hmixi,hmixfract,hmixcenter
    real :: pv1(2),pvprof(2),pvi,pvcenter,pvfract,tr(2),tri,tropofract
    real :: tropocenter

    real :: xcentertra1, ycentertra1
    real :: tt2_interp_s(2), tt2_interp_s_t
    real :: td2_interp_s(2), td2_interp_s_t
    real :: ustar_interp_s(2), ustar_interp_s_t
    real :: wstar_interp_s(2), wstar_interp_s_t
    real :: hmix_interp_s(2), hmix_interp_s_t
    real :: sshf_interp_s(2), sshf_interp_s_t
    real :: lsprec_interp_s(2), lsprec_interp_s_t
    real :: convprec_interp_s(2), convprec_interp_s_t
    real :: ssr_interp_s(2), ssr_interp_s_t
    real :: u10_interp_s(2), u10_interp_s_t
    real :: v10_interp_s(2), v10_interp_s_t
    real :: ws10_interp_s_t
    real :: ps_interp_s(2), ps_interp_s_t
    real :: msl_interp_s(2), msl_interp_s_t
    real :: oli_interp_s(2), oli_interp_s_t
    real :: rh2_interp_s_t

    real :: tt_interp_s(2), tt_interp_s_t(nzmeteo)
    real :: prs_interp_s(2), prs_interp_s_t(nzmeteo)
    real :: qv_interp_s(2), qv_interp_s_t(nzmeteo)
    real :: uu_interp_s(2), uu_interp_s_t(nzmeteo)
    real :: vv_interp_s(2), vv_interp_s_t(nzmeteo)
    real :: clwc_interp_s(2), clwc_interp_s_t(nzmeteo)

    real :: e_s, w_s, w
    real :: rh_interp_s_t(nzmeteo)
    real :: height_meteo(nzmeteo)

    dt1=real(itime-memtime(1))
    dt2=real(memtime(2)-itime)
    dtt=1./(dt1+dt2)

    ! Loop about all release points
    !******************************

    do j=1,numpoint
        if (abs(ireleasestart(j) -itime ) <= lage(nageclass)) then
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
            do i=1,numpart

                pp => particles(i)

                if ((pp%active) .and. (npoint(i) == j)) then
                    n=n+1
                    xl(n)=xlon0+xtra1(i)*dx
                    yl(n)=ylat0+ytra1(i)*dy
                    zl(n)=ztra1(i)


                    ! Interpolate PBL height, PV, and tropopause height to each
                    ! particle position in order to determine fraction of particles
                    ! within the PBL, above tropopause height, and average PV.
                    ! Interpolate topography, too, and convert to altitude asl
                    !**************************************************************

                    ix=int(xtra1(i))
                    jy=int(ytra1(i))
                    ixp=ix+1
                    jyp=jy+1
                    ddx=xtra1(i)-real(ix)
                    ddy=ytra1(i)-real(jy)
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

                    do il=2,nz
                        if (height(il).gt.zl(n)) then
                            indz=il-1
                            indzp=il
                            exit
                        endif
                    end do

                    dz1=zl(n)-height(indz)
                    dz2=height(indzp)-zl(n)
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
                    if (yl(n).gt.0.) then
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
                    if (zl(n).lt.tri) tropofract=tropofract+1.
                    tropocenter=tropocenter+tri+topo
                    if (zl(n).lt.hmixi) hmixfract=hmixfract+1.
                    zl(n)=zl(n)+topo        ! convert to height asl
                    hmixcenter=hmixcenter+hmixi


                end if
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

                call clustering(xl,yl,zl,n,xclust,yclust,zclust,fclust,rms, &
                    rmsclust,zrms)


                ! Determine center of mass position on earth and average height
                !**************************************************************

                call centerofmass(xl,yl,n,xcenter,ycenter)
                call mean(zl,zcenter,zrmsdist,n)

                ! Root mean square distance from center of mass
                !**********************************************

                do k=1,n
                    dist=distance(yl(k),xl(k),ycenter,xcenter)
                    rmsdist=rmsdist+dist*dist
                end do
                if (rmsdist.gt.0.) rmsdist=sqrt(rmsdist/real(n))
                rmsdist=max(rmsdist,0.)

                ! Write out results in trajectory data file
                !******************************************
                
                write(unitouttraj,'(i5,i9,2f10.4,4f8.1,f8.2,4f8.1,3f6.1,&
                    &5(2f9.3,f7.0,f6.1,f8.1))')&
                    &j,itime-(ireleasestart(j)+ireleaseend(j))/2, &
                    xcenter,ycenter,zcenter,topocenter,hmixcenter,tropocenter, &
                    pvcenter,rmsdist,rms,zrmsdist,zrms,hmixfract,pvfract, &
                    tropofract,&
                    (xclust(k),yclust(k),zclust(k),fclust(k),rmsclust(k), &
                    k=1,ncluster)

                ! Get other meteorological variables from center plume
                ! Interpolate based on the centerofmass and mean height
                !****************************************************
                xcentertra1=(xcenter-xlon0)/dx
                ycentertra1=(ycenter-ylat0)/dy

                ix=int(xcentertra1)
                jy=int(ycentertra1)
                ixp=ix+1
                jyp=jy+1
                ddx=xcentertra1-real(ix)
                ddy=ycentertra1-real(jy)
                rddx=1.-ddx
                rddy=1.-ddy
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy

                topo=p1*oro(ix ,jy) &
                    + p2*oro(ixp,jy) &
                    + p3*oro(ix ,jyp) &
                    + p4*oro(ixp,jyp)

                ! variables from 2d field
                ! Want ustar, hmix, sshf, rain, windspeed, ssr, tt2, ps, rh (from td)
                !************************

                do m=1,2
                    indexh=memind(m)

                    tt2_interp_s(m)=p1*tt2(ix ,jy ,1,indexh) &
                                    + p2*tt2(ixp,jy ,1,indexh) &
                                    + p3*tt2(ix ,jyp,1,indexh) &
                                    + p4*tt2(ixp,jyp,1,indexh)

                    td2_interp_s(m)=p1*td2(ix ,jy ,1,indexh) &
                                    + p2*td2(ixp,jy ,1,indexh) &
                                    + p3*td2(ix ,jyp,1,indexh) &
                                    + p4*td2(ixp,jyp,1,indexh)
                    
                    ustar_interp_s(m)=p1*ustar(ix ,jy ,1,indexh) &
                                    + p2*ustar(ixp,jy ,1,indexh) &
                                    + p3*ustar(ix ,jyp,1,indexh) &
                                    + p4*ustar(ixp,jyp,1,indexh) 
                                    
                    wstar_interp_s(m)=p1*wstar(ix ,jy ,1,indexh) &
                                    + p2*wstar(ixp,jy ,1,indexh) &
                                    + p3*wstar(ix ,jyp,1,indexh) &
                                    + p4*wstar(ixp,jyp,1,indexh)

                    hmix_interp_s(m)=p1*hmix(ix ,jy ,1,indexh) &
                                    + p2*hmix(ixp,jy ,1,indexh) &
                                    + p3*hmix(ix ,jyp,1,indexh) &
                                    + p4*hmix(ixp,jyp,1,indexh)

                    sshf_interp_s(m)=p1*sshf(ix ,jy ,1,indexh) &
                                    + p2*sshf(ixp,jy ,1,indexh) &
                                    + p3*sshf(ix ,jyp,1,indexh) &
                                    + p4*sshf(ixp,jyp,1,indexh)

                    lsprec_interp_s(m)=p1*lsprec(ix ,jy ,1,indexh) &
                                    + p2*lsprec(ixp,jy ,1,indexh) &
                                    + p3*lsprec(ix ,jyp,1,indexh) &
                                    + p4*lsprec(ixp,jyp,1,indexh)

                    convprec_interp_s(m)=p1*convprec(ix ,jy ,1,indexh) &
                                    + p2*convprec(ixp,jy ,1,indexh) &
                                    + p3*convprec(ix ,jyp,1,indexh) &
                                    + p4*convprec(ixp,jyp,1,indexh)

                    ssr_interp_s(m)=p1*ssr(ix ,jy ,1,indexh) &
                                    + p2*ssr(ixp,jy ,1,indexh) &
                                    + p3*ssr(ix ,jyp,1,indexh) &
                                    + p4*ssr(ixp,jyp,1,indexh)

                    u10_interp_s(m)=p1*u10(ix ,jy ,1,indexh) &
                                    + p2*u10(ixp,jy ,1,indexh) &
                                    + p3*u10(ix ,jyp,1,indexh) &
                                    + p4*u10(ixp,jyp,1,indexh)

                    v10_interp_s(m)=p1*v10(ix ,jy ,1,indexh) &
                                    + p2*v10(ixp,jy ,1,indexh) &
                                    + p3*v10(ix ,jyp,1,indexh) &
                                    + p4*v10(ixp,jyp,1,indexh)
                    
                    ps_interp_s(m)=p1*ps(ix ,jy ,1,indexh) &
                                    + p2*ps(ixp,jy ,1,indexh) &
                                    + p3*ps(ix ,jyp,1,indexh) &
                                    + p4*ps(ixp,jyp,1,indexh)
                    
                    msl_interp_s(m)=p1*msl(ix ,jy ,1,indexh) &
                                    + p2*msl(ixp,jy ,1,indexh) &
                                    + p3*msl(ix ,jyp,1,indexh) &
                                    + p4*msl(ixp,jyp,1,indexh)

                    oli_interp_s(m)=p1*oli(ix ,jy ,1,indexh) &
                                    + p2*oli(ixp,jy ,1,indexh) &
                                    + p3*oli(ix ,jyp,1,indexh) &
                                    + p4*oli(ixp,jyp,1,indexh)


                end do

                tt2_interp_s_t=(tt2_interp_s(1)*dt2+tt2_interp_s(2)*dt1)*dtt
                
                td2_interp_s_t=(td2_interp_s(1)*dt2+td2_interp_s(2)*dt1)*dtt

                ustar_interp_s_t=(ustar_interp_s(1)*dt2+ustar_interp_s(2)*dt1)*dtt

                wstar_interp_s_t=(wstar_interp_s(1)*dt2+wstar_interp_s(2)*dt1)*dtt

                hmix_interp_s_t=(hmix_interp_s(1)*dt2+hmix_interp_s(2)*dt1)*dtt

                sshf_interp_s_t=(sshf_interp_s(1)*dt2+sshf_interp_s(2)*dt1)*dtt

                lsprec_interp_s_t=(lsprec_interp_s(1)*dt2+lsprec_interp_s(2)*dt1)*dtt

                convprec_interp_s_t=(convprec_interp_s(1)*dt2+convprec_interp_s(2)*dt1)*dtt

                ssr_interp_s_t=(ssr_interp_s(1)*dt2+ssr_interp_s(2)*dt1)*dtt
                
                u10_interp_s_t=(u10_interp_s(1)*dt2+u10_interp_s(2)*dt1)*dtt
                
                v10_interp_s_t=(v10_interp_s(1)*dt2+v10_interp_s(2)*dt1)*dtt

                ws10_interp_s_t=sqrt(u10_interp_s_t*u10_interp_s_t+v10_interp_s_t*v10_interp_s_T)
                
                ps_interp_s_t=(ps_interp_s(1)*dt2+ps_interp_s(2)*dt1)*dtt/100 !hPa
                
                msl_interp_s_t=(msl_interp_s(1)*dt2+ps_interp_s(2)*dt1)*dtt/100 !hPa
                
                oli_interp_s_t=(oli_interp_s(1)*dt2+oli_interp_s(2)*dt1)*dtt

                rh2_interp_s_t = 100*exp(17.625*(td2_interp_s_t-273.15)/(243.04+td2_interp_s_t-273.15))/exp(17.625*(tt2_interp_s_t-273.15)/(243.04+tt2_interp_s_t-273.15))

                ! variable from 3d fields
                ! want tt, prs, qv, rho, rho_dry, uu, vv, clwc, RHs, u_vel, v_vel, cloudLWC
                !*****************************************

                !Loop over height indicies, upp to maximum pblh  
                do il=1,nzmeteo

                    do m=1,2
                        indexh=memind(m)
    
                        tt_interp_s(m)=p1*tt(ix ,jy ,il,indexh) &
                                        + p2*tt(ixp,jy ,il,indexh) &
                                        + p3*tt(ix ,jyp,il,indexh) &
                                        + p4*tt(ixp,jyp,il,indexh)

                        prs_interp_s(m)=p1*prs(ix ,jy ,il,indexh) &
                                        + p2*prs(ixp,jy ,il,indexh) &
                                        + p3*prs(ix ,jyp,il,indexh) &
                                        + p4*prs(ixp,jyp,il,indexh)

                        qv_interp_s(m)=p1*qv(ix ,jy ,il,indexh) &
                                        + p2*qv(ixp,jy ,il,indexh) &
                                        + p3*qv(ix ,jyp,il,indexh) &
                                        + p4*qv(ixp,jyp,il,indexh)

                        uu_interp_s(m)=p1*uu(ix ,jy ,il,indexh) &
                                        + p2*uu(ixp,jy ,il,indexh) &
                                        + p3*uu(ix ,jyp,il,indexh) &
                                        + p4*uu(ixp,jyp,il,indexh)

                        vv_interp_s(m)=p1*vv(ix ,jy ,il,indexh) &
                                        + p2*vv(ixp,jy ,il,indexh) &
                                        + p3*vv(ix ,jyp,il,indexh) &
                                        + p4*vv(ixp,jyp,il,indexh)

                        clwc_interp_s(m)=p1*clwc(ix ,jy ,il,indexh) &
                                        + p2*clwc(ixp,jy ,il,indexh) &
                                        + p3*clwc(ix ,jyp,il,indexh) &
                                        + p4*clwc(ixp,jyp,il,indexh)
                        
                        
                    end do

                    height_meteo(il) = height(il)

                    tt_interp_s_t(il)=(tt_interp_s(1)*dt2+tt_interp_s(2)*dt1)*dtt
                    
                    prs_interp_s_t(il)=(prs_interp_s(1)*dt2+prs_interp_s(2)*dt1)*dtt/100 !hPa
                    
                    qv_interp_s_t(il)=(qv_interp_s(1)*dt2+qv_interp_s(2)*dt1)*dtt
                    
                    uu_interp_s_t(il)=(uu_interp_s(1)*dt2+uu_interp_s(2)*dt1)*dtt
                    
                    vv_interp_s_t(il)=(vv_interp_s(1)*dt2+vv_interp_s(2)*dt1)*dtt
                    
                    clwc_interp_s_t(il)=(clwc_interp_s(1)*dt2+clwc_interp_s(2)*dt1)*dtt

                    e_s =  611*exp(17.67*(tt_interp_s_t(il)-273.16)/(tt_interp_s_t(il)-29.65))

                    w_s = e_s*r_air/(r_water*(prs_interp_s_t(il)-e_s))
                    
                    w = qv_interp_s_t(il)/(1-qv_interp_s_t(il))

                    rh_interp_s_t(il) = 100*w/(r_air/r_water+w)*(r_air/r_water+w_s)/w_s


                end do



                ! Write out results in trajectory meteo data file
                !NOTE 
                !******************************************
                write(unitouttrajmeteo,&
                    '(i5,i9,&
                    f8.2,& 
                    f8.2,&
                    f8.2,&
                    f7.4,&
                    f8.4,&
                    f9.2,&
                    f9.3,&
                    f10.4,&
                    f10.4,&
                    f9.2,&
                    f7.2,&
                    f7.2,&
                    f7.2,&
                    f9.2,&
                    f9.2,&
                    f10.4,&
                    f9.2,&
                    500f10.3)')&
                    &j,itime-(ireleasestart(j)+ireleaseend(j))/2, &
                    tt2_interp_s_t,&
                    td2_interp_s_t,&
                    rh2_interp_s_t,&
                    ustar_interp_s_t,&
                    wstar_interp_s_t,&
                    hmix_interp_s_t,&
                    sshf_interp_s_t,&
                    lsprec_interp_s_t,&
                    convprec_interp_s_t,&
                    ssr_interp_s_t,&
                    u10_interp_s_t,&
                    v10_interp_s_t,&
                    ws10_interp_s_t,&
                    ps_interp_s_t,&
                    msl_interp_s_t,&
                    oli_interp_s_t,&
                    topo,&
                    ! variables with height levels:
                    height_meteo,&
                    tt_interp_s_t,&
                    prs_interp_s_t,&
                    qv_interp_s_t,&
                    uu_interp_s_t,&
                    vv_interp_s_t,&
                    clwc_interp_s_t,&
                    rh_interp_s_t
                    


            endif


        end if
    end do



end subroutine plumetraj
