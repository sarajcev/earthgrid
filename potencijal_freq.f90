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

! ==========================================================================
!       PRORACUN POTENCIJALA U JEDNOJ TOCKI BILO GDJE U VISESLOJNOM
!       SREDSTVU, USLJED SUSTAVA ENERGIZIRANIH VODICA. NAJVAZNIJE SU
!       PRITOM TOCKE KOJE SE NALAZE NA POVRSINI TLA. UKOLIKO SE ZELI
!       RASPODJELA POTENCIJALA DUZ NEKOG PRAVCA, OVA RUTINA SE POZIVA
!       ZA SVAKU TOCKU NA ODABRANOM PRAVCU, UZ NEKI DOVOLJNO MALI KORAK.
! ==========================================================================

subroutine potencijal_tocke_freq(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,&
    L,Ip,HD,h,kapa,gama,potencijal)
    use funkcije
    implicit none
    
    integer N,Ns,Nsi
    integer isl
    integer,dimension(:) :: iso
    real(8) :: xt,yt,zt
    real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
    complex(8),dimension(:) :: Fr,Ip
    real(8),dimension(:) :: HD,h,L
    complex(8),dimension(:) :: kapa,gama
    !Output
    complex(8) potencijal
    
    integer Nss,ngt1
    real(8),dimension(5) :: ug1,Hug1
    real(8) x1,y1,z1,x2,y2,z2,pol
    integer si,iseg
    complex(8) Ft,faktorp
    real(8),parameter :: pi = 3.141592653589793238d0
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    real(8) xtg,ytg,ztg
    complex(8),dimension(33) :: C
    real(8),dimension(33) :: zm
    real(8),dimension(3) :: zm1,zm2
    real(8) teta,hi
    complex(8),dimension(15) :: alfa,beta
    real(8),dimension(15) :: eta
    real(8),dimension(25) :: w
    real(8) xsi,ysi,zsi,xsj,ysj,zsj,xs,ys,zs
    real(8) andi,korjen
    complex(8) nudi,integral,gamae,faktore
    complex(8) Fprig,Pk
    real(8) zt1,zt2
    integer ii,ig,j
    
    !Opis ulaznih varijabli:
    !  N - ukupni broj slojeva modela (veseslojno tlo + zrak)
    !  Ns - ukupni broj segmenata uzemljivackih sustava
    !  Nsi - ukupni broj izoliranih segmenata
    !  xt - x koordinata tocke promatranja (tocke u kojoj se racuna raspodjela potenc.)
    !  yt - y koordinata tocke promatranja
    !  zt - z koordinata tocke promatranja
    !  isl - redni broj sloja u kojem se nalazi tocka promatranja
    !  xp (realni vektor) - x koordinate pocetnih tocaka svih segmenata
    !  yp (realni vektor) - y koordinate pocetnih tocaka svih segmenata
    !  zp (realni vektor) - z koordinate pocetnih tocaka svih segmenata
    !  xk (realni vektor) - x koordinate krajnjih tocaka svih segmenata
    !  yk (realni vektor) - y koordinate krajnjih tocaka svih segmenata
    !  zk (realni vektor) - z koordinate krajnjih tocaka svih segmenata
    !  r0 (realni vektor) - radijus svakog od segmenata
    !  iso (realni vektor) - redni broj sloja u kojem se nalazi svaki od segmenata
    !  L (realni vektor) - duljine svih segmenata sustava uzemljivaca
    !  Fr (complex vektor) - faktori refleksije svih slojeva
    !  HD (realni vektor) - z koordinate donjih granicnih ploha svih slojeva
    !  h (realni vektor) - debljine svih slojeva modela
    !  Ip (complex vektor) - vrijednosti poprecnih struja svih segmenata
    !  kapa (complex vektor) - kompleksne specificne elektr. vodljivosti svih slojeva
    !  gama (complex vektor) - kompleksne valne konstante svih slojeva modela
    

    ! --------------------------------------------------------------------------------
    !                POCETAK PRORACUNA POTENCIJALA U JEDNOJ TOCKI
    ! --------------------------------------------------------------------------------
    Nss = Ns - Nsi

    IF (N==2) THEN
        ! ............................................................................
        !             DVOSLOJNO SREDSTVO -> ZRAK + HOMOGENO TLO
        ! ............................................................................
        potencijal = dcmplx(0.d0,0.d0)
        do iseg = 1,Nss
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            pol = r0(iseg)
            Ft = vektor_A(si,isl,N,Fr)
            faktorp = (one/(4.d0*pi*kapa(isl)))*Ip(iseg)
            if (dabs(z2-z1)<=0.001) then
                ! ------------------------------------
                ! SEGMENT JE U HORIZONTALNOM POLOZAJU
                ! -------------------------------------
                ! Koeficijenti odslikavanja homogenog tla
                C(1:3) = dcmplx(0.d0,0.d0)
                zsj = (z2+z1)/2.d0
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
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                do j = 1,3
                    zsi = zm(j); z2 = zm(j)
                    ! Analiticki dio
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                potencijal = potencijal + integral*(faktorp/L(iseg)) * Fprig
            else
                ! ------------------------------------
                !     SEGMENT JE U KOSOM POLOZAJU
                ! -------------------------------------
                ! Koeficijenti odslikavanja homogenog tla
                C(1:3) = dcmplx(0.d0,0.d0)
                zsi = (z2+z1)/2.d0
                zm1(1) = zsi
                zm2(1) = z2
                C(1) = Ft
                if ((isl<=si).and.(isl/=1)) then
                    zm1(2) = 2.d0*HD(isl-1)-zsi
                    zm2(2) = 2.d0*HD(isl-1)-z2
                    C(2) = -Ft*Fr(isl-1)
                else if((isl>=si).and.(si/=1)) then
                    zm1(2) = 2.d0*HD(si-1)-zsi
                    zm2(2) = 2.d0*HD(si-1)-z2
                    C(2) = -Ft*Fr(si-1)
                end if
                if ((isl<=si).and.(si/=N)) then
                    zm1(3) = 2.d0*HD(si)-zsi
                    zm2(3) = 2.d0*HD(si)-z2
                    C(3) = Ft*Fr(si)
                else if((isl>=si).and.(isl/=N)) then
                    zm1(3) = 2.d0*HD(isl)-zsi
                    zm2(3) = 2.d0*HD(isl)-z2
                    C(3) = Ft*Fr(isl)
                end if
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j)
                    ! Analiticki dio
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                potencijal = potencijal + integral*(faktorp/L(iseg)) * Fprig
            end if
        end do
    ELSE
        ! ............................................................................
        !          OPCENITO VISESLOJNO SREDSTVO -> ZRAK + VISESLOJNO TLO
        ! ............................................................................
        potencijal = dcmplx(0.d0,0.d0)
        ngt1 = 5
        call kGint(ngt1,ug1,Hug1)
        do iseg = 1,Nss
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            pol = r0(iseg)
            Ft = vektor_A(si,isl,N,Fr)
            faktorp = one/(4.d0*pi*kapa(isl))
            if (dabs(z2-z1)<0.001) then
                ! -----------------------------------------
                !   SEGMENT JE U HORIZONTALNOM POLOZAJU
                ! -----------------------------------------
                ! Proracun koeficijenata odslikavanja - C
                C(1:33) = dcmplx(0.d0,0.d0)
                zsj = (z2+z1)/2.d0
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
                ! Proracun impedancije izmedju horizontalnog segmenta i tocke
                ! ---------------------------------------------------------------
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                integral = dcmplx(0.d0,0.d0)
                do j = 1,33
                    zsi = zm(j); z2 = zm(j)
                    ! Analiticki dio
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                potencijal = potencijal + integral * Fprig * (faktorp/L(iseg)) * Ip(iseg)
            else
                ! -----------------------------------------
                !    SEGMENT JE U KOSOM POLOZAJU
                ! -----------------------------------------
                Pk = dcmplx(0.d0,0.d0)
                C(1:33) = dcmplx(0.d0,0.d0)
                zsi = (z2+z1)/2.d0
                zm1(1) = zsi
                zm2(1) = z2
                C(1) = Ft
                if ((isl<=si).and.(isl/=1)) then
                    zm1(2) = 2.d0*HD(isl-1)-zsi
                    zm2(2) = 2.d0*HD(isl-1)-z2
                    C(2) = -Ft*Fr(isl-1)
                else if((isl>=si).and.(si/=1)) then
                    zm1(2) = 2.d0*HD(si-1)-zsi
                    zm2(2) = 2.d0*HD(si-1)-z2
                    C(2) = -Ft*Fr(si-1)
                end if
                if ((isl<=si).and.(si/=N)) then
                    zm1(3) = 2.d0*HD(si)-zsi
                    zm2(3) = 2.d0*HD(si)-z2
                    C(3) = Ft*Fr(si)
                else if((isl>=si).and.(isl/=N)) then
                    zm1(3) = 2.d0*HD(isl)-zsi
                    zm2(3) = 2.d0*HD(isl)-z2
                    C(3) = Ft*Fr(isl)
                end if
                ! Proracun impedancije izmedju kosog segmenta i tocke
                ! -------------------------------------------------------
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                integral = dcmplx(0.d0,0.d0)
                zt1 = zsi; zt2 = z2
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j)
                    ! Analiticki dio
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                Pk = Pk + integral * (faktorp/L(iseg)) * Fprig

                if (isl/=1) then
                    integral = dcmplx(0.d0,0.d0)
                    z2 = zt2; zsi = zt1
                    do ig = 1,5
                        xtg = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                        ytg = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                        ztg = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                        call parametri_teta_hi(si,isl,ztg,HD,h,N,teta,hi)
                        call vektor_eta(N,eta)
                        call vektor_w(N,w)
                        call alfauz(teta,w,eta,h,N,ztg,Fr,HD,si,kapa,isl,alfa)
                        do ii = 1,15
                            korjen = dsqrt((xt-xtg)**2+(yt-ytg)**2+(teta*eta(ii)+zt-HD(isl-1))**2)
                            faktore = one/korjen
                            integral = integral + alfa(ii)*Hug1(ig)*faktore
                        end do
                    end do
                    Pk = Pk + integral * faktorp * Fprig
                end if

                if (isl/=N) then
                    integral = dcmplx(0.d0,0.d0)
                    z2 = zt2; zsi = zt1
                    do ig = 1,5
                        xtg = 0.5d0*(x1+x2) + ug1(ig)*(x2-x1)
                        ytg = 0.5d0*(y1+y2) + ug1(ig)*(y2-y1)
                        ztg = 0.5d0*(z1+z2) + ug1(ig)*(z2-z1)
                        call parametri_teta_hi(si,isl,ztg,HD,h,N,teta,hi)
                        call vektor_eta(N,eta)
                        call vektor_w(N,w)
                        call betauz(hi,w,eta,h,N,ztg,Fr,HD,si,kapa,isl,beta)
                        do ii = 1,15
                            korjen = dsqrt((xt-xtg)**2+(yt-ytg)**2+(hi*eta(ii)+HD(isl)-zt)**2)
                            faktore = one/korjen
                            integral = integral + beta(ii)*Hug1(ig)*faktore
                        end do
                    end do
                    Pk = Pk + integral * faktorp * Fprig
                end if
                potencijal = potencijal + Pk * Ip(iseg)
            end if
        end do
    END IF

    return
end subroutine
