! -------------------------------------------------------------------------------
!    Copyright (C) 2006. GPL - General Public Licence
!    Author: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)

!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License along
!    with this program; if not, write to the Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! -------------------------------------------------------------------------------

! Subroutina za formiranje matrica [Yu] i [Yp]. Ova matrice formiraju se na sljedeci nacin:
!   [Yu] - matrica uzduznih admitancija lokalnih cvorova segmenata uzemljivaca. Ona
!          se dobiva na sljedeci nacin:

!                           [Yu] = trans([Au]) * [Yuv] * [Au]

!          Pri cemu su:
!               [Au] - pravokutna matrica veze (matrica incidencije).
!               trans([Au]) - transponirana matrica veze [Au]
!               [Yuv] - matrica uzduznih admitancija segmenata uzemljivaca (grana).
!                       Ona se racuna kako slijedi:

!                                  [Yuv] = inv([Zuv])

!                       Dakle, inverzijom matrice uzduznih impedancija segmenata 
!                       uzemljivaca. Ova matrica se takodjer ovdje formira (matrica
!                       impedancija grana).
!   [Yp] - matrica poprecnih admitancija segmenata uzemljivaca. Ona se dobiva kako
!          slijedi:

!                           [Yp] = trans([Ap]) * [Ypv] * [Ap]

!          Pri cemu su:
!               [Ap] - pravokutna matrica veze (matrica incidencije).
!               trans([Ap]) - transponirana matrica veze [Ap]
!               [Ypv] - matrica poprecnih admitancija segmenata uzemljivaca (grana). 
!                       Ona se racuna kako slijedi:

!                                  [Ypv] =  inv([Zpv])

!                       Dakle, inverzijom matrice poprecnih impedancija segmenata 
!                       uzemljivaca. Ova matrica se takodjer forimra na ovom mjestu
!                       (u ovoj subroutini).

! Osnovni problem je, dakle, formirati matrice uzduznih i poprecnih impedancija
! segmenata uzemljivaca (matrice impedancija grana [Zuv] i [Zpv]). Matrice veze
! (incidencije) [Au] i [Ap] su poznate realne ili kompleksne matrice.


subroutine matrica(f,L,xp,yp,zp,xk,yk,zk,Ns,Nsi,r0,kapa,gama,sigmav,Zunjs,&
iso,N,Fr,HD,h,Yuv,Ypv,Yu,Ypp)
    use funkcije
    implicit none
    
    real(8) f
    real(8),dimension(:) :: L
    real(8),dimension(:) :: xp,yp,zp,xk,yk,zk
    integer Ns,Nsi,N
    real(8),dimension(:) :: r0,sigmav
    complex(8),dimension(:) :: kapa,gama,Zunjs
    integer,dimension(:) :: iso
    complex(8),dimension(:) :: Fr
    real(8),dimension(:) :: HD,h
    !Output
    complex(8),dimension(:,:) :: Yuv,Ypv,Yu,Ypp

    !Lokalne varijable
    integer Nss
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    real(8),parameter :: pi = 3.141592653589793238d0
    complex(8),dimension(:,:),allocatable :: Zu,Zpp
    complex(8),dimension(:,:),allocatable :: Au,Ap,CC
    real(8) x1,y1,z1,x2,y2,z2,pol
    real(8) x3,y3,z3,x4,y4,z4
    real(8) xsi,ysi,zsi,xsj,ysj,zsj
    complex(8) Ft
    complex(8) gamae
    complex(8) faktoru,faktorp
    real(8) temp,ttemp
    real(8) andi
    complex(8) nudi,integral
    real(8) cosfi,cosfiu
    complex(8),dimension(33) :: C
    real(8),dimension(33) :: zm
    real(8),dimension(3) :: zm1,zm2
    real(8) teta,hi
    real(8),dimension(15) :: eta
    real(8),dimension(25) :: w
    complex(8),dimension(15) :: alfa,beta
    real(8) xps,xks,yps,yks,zps,zks
    integer ngt1
    real(8),dimension(5) :: ug1,Hug1
    real(8) xs,ys,zs,xt,yt,zt
    complex(8) Fprig
    integer si,isl
    integer iseg,jseg
    real(8) zt1,zt2
    integer i,j,ii,ig
    integer,dimension(:),allocatable :: ipiv
    integer info,Lwork,NB
    complex(8),dimension(:),allocatable :: work
    complex(8) alpha,betha
    character(1) transa,transb

    
    ! Opis ulaznih varijabli:
    !   f - frekvencija za koju se racuna, [Hz]
    !   L - (realni vektor) duljine svih segmenata sustava uzemljenja, [m]
    !   xp - (realni vektor) pocetna x koordinata svih segmenata sustava uzemljenja
    !   yp - (realni vektor) pocetna y koordinata svih segmenata sustava uzemljenja
    !   zp - (realni vektor) pocetna z koordinata svih segmenata sustava uzemljenja
    !   xk - (realni vektor) krajnja x koordinata svih segmenata sustava uzemljenja
    !   yk - (realni vektor) krajnja y koordinata svih segmenata sustava uzemljenja
    !   zk - (realni vektor) krajnja z koordinata svih segmenata sustava uzemljenja
    !   r0 - (realni vektor) radijus svih segmenata sustava uzemljenja, [m]
    !   Ns - ukupni broj segmenata sustava uzemljenja
    !   Nsi - ukupni broj izoliranih segmenata
    !   N - ukupni broj slojeva modela (viseslojno tlo plus zrak)
    !   kapa - (kompleksni vektor) vrijednosti specificne kompleksne elektricne 
    !          vodljivosti svih slojeva modela
    !   gama - (kompleksni vektor) vrijednosti valne konstante svih slojeva modela
    !   Zunjs - (kompleksni vektor) vrijednosti jedinicnih unutarnjih impedancija
    !           svih segmenata sustava uzemljenja (izracunava se u fileu: zunjs.f90)
    !   iso - (cjelobrojni vektor) redni broj sloja u kojem se nalazi svaki od
    !         segmenata sustava uzemljenja
    !   Fr - (kompleksni vektor) faktori refleksije svih slojeva modela
    !   HD - (realni vektor) z koordinate donjih granicnih ploha svih slojeva modela
    !   h - (realni vektor) debljine svih slojeva modela, [m]
    !   Yuv (output) - kompleksna matrica koja predstavlja matricu [Yuv]
    !   Ypv (output) - kompleksna matrica koja predstavlja matricu [Ypv]
    !   Y (output) - (kompleksna matrica) matrica [Y].
    
        
    ! ================================================================================
    !                           POCETAK PRORACUNA
    ! ================================================================================
    temp = 2.d0*pi*f*1.0e-7
    faktoru = dcmplx(0.d0,temp)
    Nss = Ns - Nsi
    
    ! ********************************************************************************
    ! IMPEDANCIJE SEGMENATA (UZDUZNE I POPRECNE) U HOMOGENOM I NEOGRANICENOM SREDSTVU
    ! ********************************************************************************
    allocate(Zu(Ns,Ns))
    allocate(Zpp(Nss,Nss))
    do iseg = 1,Nss
        si = iso(iseg)
        x1 = xp(iseg); x2 = xk(iseg)
        y1 = yp(iseg); y2 = yk(iseg)
        z1 = zp(iseg); z2 = zk(iseg)
        do jseg = iseg,Nss
            isl = iso(jseg)
            x3 = xp(jseg); x4 = xk(jseg)
            y3 = yp(jseg); y4 = yk(jseg)
            z3 = zp(jseg); z4 = zk(jseg)
            pol = r0(jseg)
            faktorp = one/(4.d0*pi*kapa(isl))
            if (iseg==jseg) then
                ! -------------------------------------------
                !      Vlastita impedancija segmenta
                ! -------------------------------------------
                ! ANALITICKI DIO INTEGRALA
                call otpss(L(jseg),pol,andi) 
                ! NUMERICKI DIO INTEGRALA
                xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                call sdnivl(gamae,L(jseg),pol,nudi)
                ! Faktor prigusenja (gamae*r0)
                Fprig = (andi + nudi)/andi
                integral = andi * Fprig
                ! Uzduzna impedancija
                if (f<1.e-16) then
                    ttemp = one/(sigmav(jseg)*pol**2*pi)
                    Zu(iseg,jseg) = dcmplx(ttemp,0.d0)
                else
                    Zu(iseg,jseg) = Zunjs(jseg)*L(jseg) + faktoru*integral
                end if
                ! Poprecna impedancija
                Zpp(iseg,jseg) = (faktorp/L(jseg)**2)*integral
            else
                ! ----------------------------------------------------
                !      Medjusobna impedancija izmedju segmenata
                ! ----------------------------------------------------
                call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                if (cosfi<0.999995) then
                    ! ................................................
                    !         Segment je u kosom polozaju
                    ! ................................................
                    call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                else
                    ! ................................................
                    !      Segment je u horizontalnom polozaju
                    ! ................................................
                    call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                end if
                xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                integral = andi * Fprig
                call koskutsegmenata1(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfiu)
                ! Uzduzna impedancija
                if (f<1.e-16) then
                    Zu(iseg,jseg) = dcmplx(0.d0,0.d0)
                    Zu(jseg,iseg) = Zu(iseg,jseg)
                else
                    Zu(iseg,jseg) = faktoru*cosfiu*integral
                    Zu(jseg,iseg) = Zu(iseg,jseg)
                end if
                ! Poprecna impedancija
                Zpp(iseg,jseg) = (faktorp/(L(iseg)*L(jseg)))*integral
                Zpp(jseg,iseg) = Zpp(iseg,jseg)
            end if
        end do
    end do
    
    ! Postoje i izolirani segmenti u sustavu uzemljenja
    IF (Nsi /= 0) THEN
        do iseg = Nss+1,Ns
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            do jseg = iseg,Ns
                isl = iso(jseg)
                x3 = xp(jseg); x4 = xk(jseg)
                y3 = yp(jseg); y4 = yk(jseg)
                z3 = zp(jseg); z4 = zk(jseg)
                pol = r0(jseg)
                faktorp = one/(4.d0*pi*kapa(isl))
                if (iseg==jseg) then
                    ! -------------------------------------------
                    !      Vlastita impedancija segmenta
                    ! -------------------------------------------
                    ! ANALITICKI DIO INTEGRALA
                    call otpss(L(jseg),pol,andi) 
                    ! NUMERICKI DIO INTEGRALA
                    xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                    ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                    zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                    call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                    call sdnivl(gamae,L(jseg),pol,nudi)
                    ! Faktor prigusenja (gamae*r0)
                    Fprig = (andi + nudi)/andi
                    integral = andi * Fprig
                    ! Uzduzna impedancija
                    if (f<1.e-16) then
                        ttemp = one/(sigmav(jseg)*pol**2*pi)
                        Zu(iseg,jseg) = dcmplx(ttemp,0.d0)
                    else
                        Zu(iseg,jseg) = Zunjs(jseg)*L(jseg) + faktoru*integral
                    end if
                else
                    ! ----------------------------------------------------
                    !      Medjusobna impedancija izmedju segmenata
                    ! ----------------------------------------------------
                    call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                    if (cosfi<0.999995) then
                        ! ................................................
                        !         Segment je u kosom polozaju
                        ! ................................................
                        call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                    else
                        ! ................................................
                        !      Segment je u horizontalnom polozaju
                        ! ................................................
                        call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                    end if
                    xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                    ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                    zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                    call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                    Fprig = zexp(-gamae)
                    integral = andi * Fprig
                    call koskutsegmenata1(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfiu)
                    ! Uzduzna impedancija
                    if (f<1.e-16) then
                        Zu(iseg,jseg) = dcmplx(0.d0,0.d0)
                        Zu(jseg,iseg) = Zu(iseg,jseg)
                    else
                        Zu(iseg,jseg) = faktoru*cosfiu*integral
                        Zu(jseg,iseg) = Zu(iseg,jseg)
                    end if
                end if
            end do
        end do
    END IF

    ! ********************************************************************************
    !          IMPEDANCIJE SEGMENATA (POPRECNE) U VISESLOJNOM SREDSTVU
    ! ********************************************************************************
    ! Uzduzne impedancije vise nije potrebno racunati jer su one upravo jednake
    ! uzduznim impedancijama u homogenom i neogranicenom sredstvu (za sada).
    IF (N==2) THEN
        ! --------------------------------------------------------------------------------
        !               DVOSLOJNO SREDSTVO -> ZRAK + HOMOGENO TLO
        ! --------------------------------------------------------------------------------
        do iseg = 1,Nss
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            do jseg = iseg,Nss
                isl = iso(jseg)
                x3 = xp(jseg); x4 = xk(jseg)
                y3 = yp(jseg); y4 = yk(jseg)
                z3 = zp(jseg); z4 = zk(jseg)
                pol = r0(jseg)
                Ft = vektor_A(si,isl,N,Fr)
                faktorp = one/(4.d0*pi*kapa(isl))
                if (iseg==jseg) then
                    ! =================================================
                    !        Vlastita impedancija segmenta
                    ! =================================================
                    if (dabs(z4-z3)<=0.001) then
                        ! SEGMENT JE U HORIZONTALNOM POLOZAJU
                        C(1:3) = dcmplx(0.d0,0.d0)
                        zsj = (z4+z3)/2.d0
                        zm(1) = zsj
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm(2) = 2.d0*HD(isl-1)-zsj
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm(2) = 2.d0*HD(si-1)-zsj
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm(3) = 2.d0*HD(si)-zsj
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm(3) = 2.d0*HD(isl)-zsj
                            C(3) = Ft*Fr(isl)
                        end if
                        ! -------------------------------------------------------
                        call otpss(L(jseg),pol,andi) 
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae) 
                        call sdnivl(gamae,L(jseg),pol,nudi)
                        ! Faktor prigusenja (gamae*r0)
                        Fprig = (andi + nudi)/andi
                        ! Vlastita impedancija segmenta
                        integral = dcmplx(0.d0,0.d0)
                        do j = 2,3
                            z3 = zm(j); z4 = zm(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                    else
                        ! SEGMENT JE U KOSOM POLOZAJU
                        C(1:3) = dcmplx(0.d0,0.d0)
                        zm1(1) = z3
                        zm2(1) = z4
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm1(2) = 2.d0*HD(isl-1)-z3
                            zm2(2) = 2.d0*HD(isl-1)-z4
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm1(2) = 2.d0*HD(si-1)-z3
                            zm2(2) = 2.d0*HD(si-1)-z4
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm1(3) = 2.d0*HD(si)-z3
                            zm2(3) = 2.d0*HD(si)-z4
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm1(3) = 2.d0*HD(isl)-z3
                            zm2(3) = 2.d0*HD(isl)-z4
                            C(3) = Ft*Fr(isl)
                        end if
                        ! -------------------------------------------------------
                        call otpss(L(jseg),pol,andi) 
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae) 
                        call sdnivl(gamae,L(jseg),pol,nudi)
                        ! Faktor prigusenja (gamae*r0)
                        Fprig = (andi + nudi)/andi
                        ! Vlastita impedancija segmenta
                        integral = dcmplx(0.d0,0.d0)
                        do j = 2,3
                            z3 = zm1(j); z4 = zm2(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                    end if
                else
                    ! =======================================================
                    !       Medjusobna impedancija izmedju segmenata
                    ! =======================================================
                    if (dabs(z4-z3)<=0.001) then
                        ! SEGMENT JE U HORIZONTALNOM POLOZAJU
                        C(1:3) = dcmplx(0.d0,0.d0)
                        zsj = (z4+z3)/2.d0
                        zm(1) = zsj
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm(2) = 2.d0*HD(isl-1)-zsj
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm(2) = 2.d0*HD(si-1)-zsj
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm(3) = 2.d0*HD(si)-zsj
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm(3) = 2.d0*HD(isl)-zsj
                            C(3) = Ft*Fr(isl)
                        end if
                        ! -------------------------------------------------------
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                        Fprig = zexp(-gamae)
                        ! Medjusobna impedancija izmedju segmenata
                        integral = dcmplx(0.d0,0.d0)
                        do j = 1,3
                            z3 = zm(j); z4 = zm(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = integral * Fprig
                        Zpp(jseg,iseg) = Zpp(iseg,jseg)
                    else
                        ! SEGMENT JE U KOSOM POLOZAJU
                        C(1:3) = dcmplx(0.d0,0.d0)
                        zm1(1) = z3
                        zm2(1) = z4
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm1(2) = 2.d0*HD(isl-1)-z3
                            zm2(2) = 2.d0*HD(isl-1)-z4
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm1(2) = 2.d0*HD(si-1)-z3
                            zm2(2) = 2.d0*HD(si-1)-z4
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm1(3) = 2.d0*HD(si)-z3
                            zm2(3) = 2.d0*HD(si)-z4
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm1(3) = 2.d0*HD(isl)-z3
                            zm2(3) = 2.d0*HD(isl)-z4
                            C(3) = Ft*Fr(isl)
                        end if
                        ! -------------------------------------------------------
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                        Fprig = zexp(-gamae)
                        ! Medjusobna impedancija izmedju segmenata
                        integral = dcmplx(0.d0,0.d0)
                        do j = 1,3
                            z3 = zm1(j); z4 = zm2(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = integral * Fprig
                        Zpp(jseg,iseg) = Zpp(iseg,jseg)
                    end if
                end if
            end do
        end do
    ELSE
        ! --------------------------------------------------------------------------------
        !           OPCENITO VISESLOJNO SREDSTVO -> ZRAK + VISESLOJNO TLO
        ! --------------------------------------------------------------------------------
        ngt1 = 5
        call kGint(ngt1,ug1,Hug1)
        do iseg = 1,Nss
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            do jseg = iseg,Nss  
                isl = iso(jseg)
                x3 = xp(jseg); x4 = xk(jseg)
                y3 = yp(jseg); y4 = yk(jseg)
                z3 = zp(jseg); z4 = zk(jseg)
                pol = r0(jseg)
                Ft = vektor_A(si,isl,N,Fr)
                faktorp = one/(4.d0*pi*kapa(isl))
                if (iseg==jseg) then
                    ! =============================================
                    !      Vlastita impedancija segmenta
                    ! =============================================
                    if (dabs(z4-z3)<0.001) then
                        ! -----------------------------------------
                        ! SEGMENT JE U HORIZONTALNOM POLOZAJU
                        ! -----------------------------------------
                        ! Proracun koeficijenata odslikavanja - C
                        C(1:33) = dcmplx(0.d0,0.d0)
                        zsj = (z4+z3)/2.d0
                        zm(1) = zsj
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm(2) = 2.d0*HD(isl-1)-zsj
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm(2) = 2.d0*HD(si-1)-zsj
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm(3) = 2.d0*HD(si)-zsj
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm(3) = 2.d0*HD(isl)-zsj
                            C(3) = Ft*Fr(isl)
                        end if
                        if (isl/=1) then
                            call parametri_teta_hi(si,isl,zsj,HD,h,N,teta,hi)
                            call vektor_eta(N,eta)
                            call vektor_w(N,w)
                            call alfauz(teta,w,eta,h,N,zsj,Fr,HD,si,kapa,isl,alfa)
                            do j = 4,18
                                zm(j) = HD(isl-1)-teta*eta(j-3)
                                C(j) = alfa(j-3)
                            end do
                        end if
                        if (isl/=N) then
                            call parametri_teta_hi(si,isl,zsj,HD,h,N,teta,hi)
                            call vektor_eta(N,eta)
                            call vektor_w(N,w)
                            call betauz(hi,w,eta,h,N,zsj,Fr,HD,si,kapa,isl,beta)
                            do j = 19,33
                                zm(j) = HD(isl)+hi*eta(j-18)
                                C(j) = beta(j-18)
                            end do
                        end if
                        ! Proracun vlastite (poprecne) impedancije horizontalnog segmenta
                        ! ---------------------------------------------------------------
                        call otpss(L(jseg),pol,andi) 
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae) 
                        call sdnivl(gamae,L(jseg),pol,nudi)
                        ! Faktor prigusenja (gamae*r0)
                        Fprig = (andi + nudi)/andi
                        integral = dcmplx(0.d0,0.d0)
                        do j = 2,33
                            z3 = zm(j); z4 = zm(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                    else
                        ! --------------------------------------------
                        !    SEGMENT JE U KOSOM POLOZAJU - vlastita
                        ! --------------------------------------------
                        C(1:33) = dcmplx(0.d0,0.d0)
                        zm1(1) = z3
                        zm2(1) = z4
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm1(2) = 2.d0*HD(isl-1)-z3
                            zm2(2) = 2.d0*HD(isl-1)-z4
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm1(2) = 2.d0*HD(si-1)-z3
                            zm2(2) = 2.d0*HD(si-1)-z4
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm1(3) = 2.d0*HD(si)-z3
                            zm2(3) = 2.d0*HD(si)-z4
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm1(3) = 2.d0*HD(isl)-z3
                            zm2(3) = 2.d0*HD(isl)-z4
                            C(3) = Ft*Fr(isl)
                        end if
                        ! Proracun vlastite (poprecne) impedancije kosog segmenta
                        ! -------------------------------------------------------
                        call otpss(L(jseg),pol,andi) 
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae) 
                        call sdnivl(gamae,L(jseg),pol,nudi)
                        ! Faktor prigusenja (gamae*r0)
                        Fprig = (andi + nudi)/andi
                        integral = dcmplx(0.d0,0.d0)
                        zt1 = z3; zt2 = z4
                        do j = 2,3
                            z3 = zm1(j); z4 = zm2(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j)*andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                        
                        if (isl/=1) then
                            integral = dcmplx(0.d0,0.d0)
                            z3 = zt1; z4 = zt2
                            do ig = 1,5
                                xt = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                                yt = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                                zt = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                                call parametri_teta_hi(si,isl,zt,HD,h,N,teta,hi)
                                call vektor_eta(N,eta)
                                call vektor_w(N,w)
                                call alfauz(teta,w,eta,h,N,zt,Fr,HD,si,kapa,isl,alfa)
                                do ii = 1,15
                                    zt = HD(isl)-teta*eta(ii)
                                    ! Analiticki dio
                                    call pot(x4,y4,z4,xt,yt,zt,L(jseg),xsj,ysj,zsj,pol,andi)
                                    integral = integral + alfa(ii)*Hug1(ig)*andi*(faktorp/L(jseg))
                                end do
                            end do
                            Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                        end if

                        if (isl/=N) then
                            integral = dcmplx(0.d0,0.d0)
                            z3 = zt1; z4 = zt2
                            do ig = 1,5
                                xt = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                                yt = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                                zt = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                                call parametri_teta_hi(si,isl,zt,HD,h,N,teta,hi)
                                call vektor_eta(N,eta)
                                call vektor_w(N,w)
                                call betauz(hi,w,eta,h,N,zt,Fr,HD,si,kapa,isl,beta)
                                do ii = 1,15
                                    zt = HD(isl+1)+hi*eta(ii)
                                    ! Analiticki dio
                                    call pot(x4,y4,z4,xt,yt,zt,L(jseg),xsj,ysj,zsj,pol,andi)
                                    integral = integral + beta(ii)*Hug1(ig)*andi*(faktorp/L(jseg))
                                end do
                            end do
                            Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                        end if
                    end if
                else
                    ! =======================================================
                    !       Medjusobna impedancija izmedju segmenata
                    ! =======================================================
                    if (dabs(z2-z1)<0.001) then
                        ! ---------------------------------------------------
                        !      SEGMENT JE U HORIZONTALNOM POLOZAJU
                        ! ---------------------------------------------------
                        ! Proracun koeficijenata odslikavanja - C
                        C(1:33) = dcmplx(0.d0,0.d0)
                        zsj = (z4+z3)/2.d0
                        zm(1) = zsj
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm(2) = 2.d0*HD(isl-1)-zsj
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm(2) = 2.d0*HD(si-1)-zsj
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm(3) = 2.d0*HD(si)-zsj
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm(3) = 2.d0*HD(isl)-zsj
                            C(3) = Ft*Fr(isl)
                        end if
                        if (isl/=1) then
                            call parametri_teta_hi(si,isl,zsj,HD,h,N,teta,hi)
                            call vektor_eta(N,eta)
                            call vektor_w(N,w)
                            call alfauz(teta,w,eta,h,N,zsj,Fr,HD,si,kapa,isl,alfa)
                            do j = 4,18
                                zm(j) = HD(isl-1)-teta*eta(j-3)
                                C(j) = alfa(j-3)
                            end do
                        end if
                        if (isl/=N) then
                            call parametri_teta_hi(si,isl,zsj,HD,h,N,teta,hi)
                            call vektor_eta(N,eta)
                            call vektor_w(N,w)
                            call betauz(hi,w,eta,h,N,zsj,Fr,HD,si,kapa,isl,beta)
                            do j = 19,33
                                zm(j) = HD(isl)+hi*eta(j-18)
                                C(j) = beta(j-18)
                            end do
                        end if
                        ! Proracun medjusobne impedancije segmenta u horizontalnom polozaju
                        ! -----------------------------------------------------------------
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                        Fprig = zexp(-gamae)
                        integral = dcmplx(0.d0,0.d0)
                        do j = 1,33
                            z3 = zm(j); z4 = zm(j)
                            ! Analiticki dio
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg)))
                        end do
                        Zpp(iseg,jseg) = integral * Fprig
                        Zpp(jseg,iseg) = Zpp(iseg,jseg)
                    else
                        ! -------------------------------------------------
                        !      SEGMENT JE U KOSOM POLOZAJU - medjusobna
                        ! -------------------------------------------------
                        C(1:33) = dcmplx(0.d0,0.d0)
                        zm1(1) = z3
                        zm2(1) = z4
                        C(1) = Ft
                        if ((isl<=si).and.(isl/=1)) then
                            zm1(2) = 2.d0*HD(isl-1)-z3
                            zm2(2) = 2.d0*HD(isl-1)-z4
                            C(2) = -Ft*Fr(isl-1)
                        else if((isl>=si).and.(si/=1)) then
                            zm1(2) = 2.d0*HD(si-1)-z3
                            zm2(2) = 2.d0*HD(si-1)-z4
                            C(2) = -Ft*Fr(si-1)
                        end if
                        if ((isl<=si).and.(si/=N)) then
                            zm1(3) = 2.d0*HD(si)-z3
                            zm2(3) = 2.d0*HD(si)-z4
                            C(3) = Ft*Fr(si)
                        else if((isl>=si).and.(isl/=N)) then
                            zm1(3) = 2.d0*HD(isl)-z3
                            zm2(3) = 2.d0*HD(isl)-z4
                            C(3) = Ft*Fr(isl)
                        end if
                        ! Proracun medjusobne impedancije segmenta u kosom polozaju
                        ! ---------------------------------------------------------
                        xsi = (x2+x1)/2.d0; xsj = (x4+x3)/2.d0
                        ysi = (y2+y1)/2.d0; ysj = (y4+y3)/2.d0
                        zsi = (z2+z1)/2.d0; zsj = (z4+z3)/2.d0
                        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                        Fprig = zexp(-gamae)
                        integral = dcmplx(0.d0,0.d0)
                        zt1 = z3; zt2 = z4
                        do j = 1,3
                            z3 = zm1(j); z4 = zm2(j)
                            call koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
                            ! Analiticki dio
                            if (cosfi<0.999995) then
                                call otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,cosfi,pol,andi)
                            else
                                call otpps(L(jseg),x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,pol,andi)
                            end if
                            integral = integral + C(j) * andi * (faktorp/(L(iseg)*L(jseg))) 
                        end do
                        Zpp(iseg,jseg) = integral * Fprig
                        Zpp(jseg,iseg) = Zpp(iseg,jseg)
                        
                        if (isl/=1) then
                            integral = dcmplx(0.d0,0.d0)
                            z3 = zt1; z4 = zt2
                            do ig = 1,5
                                xt = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                                yt = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                                zt = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                                call parametri_teta_hi(si,isl,zt,HD,h,N,teta,hi)
                                call vektor_eta(N,eta)
                                call vektor_w(N,w)
                                call alfauz(teta,w,eta,h,N,zt,Fr,HD,si,kapa,isl,alfa)
                                do ii = 1,15
                                    zt = HD(isl)-teta*eta(ii)
                                    ! Analiticki dio
                                    call pot(x4,y4,z4,xt,yt,zt,L(jseg),xsj,ysj,zsj,pol,andi)
                                    integral = integral + alfa(ii)*Hug1(ig)*andi*(faktorp/L(jseg))
                                end do
                            end do
                            Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                            Zpp(jseg,iseg) = Zpp(iseg,jseg)
                        end if
                        
                        if (isl/=N) then
                            integral = dcmplx(0.d0,0.d0)
                            z3 = zt1; z4 = zt2
                            do ig = 1,5
                                xt = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                                yt = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                                zt = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                                call parametri_teta_hi(si,isl,zt,HD,h,N,teta,hi)
                                call vektor_eta(N,eta)
                                call vektor_w(N,w)
                                call betauz(hi,w,eta,h,N,zt,Fr,HD,si,kapa,isl,beta)
                                do ii = 1,15
                                    zt = HD(isl+1)+hi*eta(ii)
                                    ! Analiticki dio
                                    call pot(x4,y4,z4,xt,yt,zt,L(jseg),xsj,ysj,zsj,pol,andi)
                                    integral = integral + beta(ii)*Hug1(ig)*andi*(faktorp/L(jseg))
                                end do
                            end do
                            Zpp(iseg,jseg) = Zpp(iseg,jseg) + integral * Fprig
                            Zpp(jseg,iseg) = Zpp(iseg,jseg)
                        end if
                    end if
                end if
            end do
        end do
    END IF

    ! *********************************************************************************
    !                     ODREDJIVANJE MATRICA [Yu] i [Ypp] 
    ! *********************************************************************************
    ! Kratak opis tijeka proracuna:
    ! -----------------------------
    ! 1) Invertira se matrica [Zu] da se dobije matrica [Yuv] - to je matrica
    !    segmenata sustava uzemljenja - matrica grana. Za to se koriste dvije LAPACK
    !    subrutine naslova: ZGETRF - LU faktorizacija matrice (s pivotiranjem), ZGETRI -
    !    inverzija matrice (koristenjem prethodno izracunate faktorizacije).
    ! 2) Formira se matrica incidencije [Au], koja je pravokutna matrica (u nasem slucaju
    !    kompleksna - ali moze biti i realna).
    ! 3) Slijedi racunanje sljedeceg izraza: [Yu] = trans([Au])*[Yuv]*[Au]. Za ovaj
    !    proracun koristi se BLAS3 rutina: ZGEMM.
    !    Napomena: Primjenjena oznaka trans([Au]) oznacava transponiranu matricu [Au].
    !              Ovo transponiranje nije potrebno exsplicite provoditi, jer je ta 
    !              opcija sadrzana u subrutinama, cime se ubrzava proces proracuna.
        
    ! 4) Invertira se matrica [Zp] da se dobije matrica [Ypv] na isti nacin kao sto je
    !    opisano pod tockom (1). 
    ! 5) Formira se nova matrica incidencije [Ap], koja je takodjer pravokutna i
    !    kompleksna.
    ! 6) Slijedi ponovno na potpuno analogan nacin kao sto je opisano pod tockom (3) 
    !    racunanje izraza: [Yp] = trans([Ap])*[Ypv]*[Ap] 
    
    ! Napomena: Potrebno je bilo sacuvati matrice [Yuv] i [Ypv] koje se koriste u 
    !           kasnijim proracunima. 
    
    ! ==============================================================================
    !        ODREDJIVANJE MATRICE [Yu] = trans([Au])*inv([Zu])*[Au] - LAPACK
    ! ==============================================================================
    allocate(ipiv(Ns))
    call zgetrf(Ns,Ns,Zu,Ns,ipiv,info)
    NB = ILAENV( 1, 'ZGETRI', ' ', Ns, -1, -1, -1 )
    Lwork = Ns*NB
    allocate(work(Lwork))
    call zgetri(Ns,Zu,Ns,ipiv,work,Lwork,info)
    Yuv = Zu
    deallocate(Zu)
    deallocate(ipiv,work)

    ! FORMIRANJE MATRICE [Au]
    allocate(Au(Ns,2*Ns))
    Au(:,:) = dcmplx(0.d0,0.d0)
    do i = 1,Ns
        do j = 1,Ns
            if (i==j) then      
                Au(i,j) = dcmplx(1.d0,0.d0)
            end if
        end do
        do j = Ns+1,2*Ns
            if ((i+Ns)==j) then
                Au(i,j) = dcmplx(-1.d0,0.d0)
            end if
        end do
    end do
    
    ! RACUNANJE IZRAZA [Yu] = trans([Au])*[Yuv]*[Au] 
    ! ------------------------------------------------------
    transa = 'T'
    transb = 'N'
    alpha = dcmplx(1.d0,0.d0)
    betha = dcmplx(0.d0,0.d0)
    allocate(CC(2*Ns,Ns))
    call zgemm(transa,transb,2*Ns,Ns,Ns,alpha,Au,Ns,Yuv,Ns,betha,CC,2*Ns)
    transa = 'N'
    transb = 'N'
    call zgemm(transa,transb,2*Ns,2*Ns,Ns,alpha,CC,2*Ns,Au,Ns,betha,Yu,2*Ns)
    
    deallocate(CC)
    deallocate(Au)

    ! ==============================================================================
    !        ODREDJIVANJE MATRICE [Yp] = trans([Ap])*inv([Zps])*[Ap] - LAPACK
    ! ==============================================================================
    allocate(ipiv(Nss))
    call zgetrf(Nss,Nss,Zpp,Nss,ipiv,info)
    NB = ILAENV( 1, 'ZGETRI', ' ', Nss, -1, -1, -1 )
    Lwork = Nss*NB
    allocate(work(Lwork))
    call zgetri(Nss,Zpp,Nss,ipiv,work,Lwork,info)
    Ypv = Zpp
    deallocate(Zpp)
    deallocate(ipiv,work)

    ! FORMIRANJE MATRICE [Ap]
    allocate(Ap(Nss,2*Nss))
    Ap(:,:) = dcmplx(0.d0,0.d0)
    do i = 1,Nss
        do j = 1,Nss
            if (i==j) then
                Ap(i,j) = dcmplx(0.5d0,0.d0)
            end if
        end do
        do j = Nss+1,2*Nss
            if ((i+Nss)==j) then
                Ap(i,j) = dcmplx(0.5d0,0.d0)
            end if
        end do
    end do
    
    ! RACUNANJE IZRAZA [Yp] = trans([Ap])*[Yp]*[Ap] - BLAS3
    ! ------------------------------------------------------
    transa = 'T'
    transb = 'N'
    alpha = dcmplx(1.d0,0.d0)
    betha = dcmplx(0.d0,0.d0)
    allocate(CC(2*Nss,Nss))
    call zgemm(transa,transb,2*Nss,Nss,Nss,alpha,Ap,Nss,Ypv,Nss,betha,CC,2*Nss)
    transa = 'N'
    transb = 'N'
    call zgemm(transa,transb,2*Nss,2*Nss,Nss,alpha,CC,2*Nss,Ap,Nss,betha,Ypp,2*Nss)
    
    deallocate(CC)
    deallocate(Ap)

    return
end subroutine
