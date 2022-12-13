! ...
! Proracun Fiz komponente potrebne za proracun Ez komponente elektricnog polja.
! Naime, vrijedi opcenito sljedeci izraz:

!             @Fi
!     Ez = - ----- - j*omega*Az = -Fiz -j*omega*Az
!              @z

! Ova subroutina racuna Fiz za potrebe gornjeg izraza, i to od ukupnog broja
! segmenata, dakle, ukupni Fiz (od svih segmenata). U ovoj subroutini impleme-
! ntiran je dvosloj i visesloj, kao i horizontalni i kosi polozaj segmenata.

subroutine Epolje_Fiz(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,&
    L,Ip,HD,h,kapa,gama,Fiz)
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
    complex(8) Fiz

    ! Local variables
    integer Nss
    complex(8) prvi_dio,drugi_dio,treci_dio,cetvrti_dio
    real(8) u,v,DuDz,DvDz,DFDu,DFDv,DFDz,temp1,temp2
    integer ngt1
    real(8),dimension(5) :: ug1,Hug1
    real(8) x1,y1,z1,x2,y2,z2,pol
    integer si,iseg
    complex(8) Ft,faktorp
    real(8),parameter :: pi = 3.141592653589793238d0
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    real(8) xtg,ytg,ztg
    complex(8),dimension(33) :: C
    real(8),dimension(33) :: zm
    real(8),dimension(3) :: zm1,zm2,zm3
    real(8) teta,hi
    complex(8),dimension(15) :: alfa,beta
    real(8),dimension(15) :: eta
    real(8),dimension(25) :: w
    real(8) xsi,ysi,zsi,xsj,ysj,zsj,xs,ys,zs
    real(8) andi,korjen,Rs
    complex(8) nudi,integral,gamae,faktore,ggama
    complex(8) Fprig
    real(8) zt1,zt2,zt3
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
    !                POCETAK PRORACUNA Fiz ZA PRORACUN Ez
    ! --------------------------------------------------------------------------------
    Nss = Ns - Nsi

    IF (N==2) THEN
        ! ............................................................................
        !             DVOSLOJNO SREDSTVO -> ZRAK + HOMOGENO TLO
        ! ............................................................................
        Fiz = dcmplx(0.d0,0.d0)
        prvi_dio = dcmplx(0.d0,0.d0)
        drugi_dio = dcmplx(0.d0,0.d0)
        DO iseg = 1,Nss
            si = iso(iseg)
            x1 = xp(iseg); x2 = xk(iseg)
            y1 = yp(iseg); y2 = yk(iseg)
            z1 = zp(iseg); z2 = zk(iseg)
            pol = r0(iseg)
            Ft = vektor_A(si,isl,N,Fr)
            faktorp = (one/(4.d0*pi*kapa(isl)))*(Ip(iseg)/L(iseg))
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
                Rs = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2)
                ggama = gamae/Rs
                do j = 1,3
                    zsi = zm(j); z2 = zm(j)
                    ! Analiticki dio - ANDI
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                zsi = (zk(iseg)+zp(iseg))/2.d0
                prvi_dio = faktorp*(-ggama)*Fprig*((zt-zsi)/Rs)*integral
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                do j = 1,3
                    zsi = zm(j); z2 = zm(j); z1 = zm(j)
                    ! Derivacija - DFDz
                    u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
                    v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
                    if (v < pol) v = pol
                    DuDz = (z2-z1)/L(iseg)
                    DvDz = (1.d0/v) * (zt-zsi-u*((z2-z1)/L(iseg)))
                    temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
                    temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
                    DFDu = 1.d0/temp1 - 1.d0/temp2
                    DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
                    DFDz = DFDu*DuDz + DFDv*DvDz
                    integral = integral + C(j)*DFDz
                end do
                drugi_dio = faktorp*Fprig*integral
            else
                ! ------------------------------------
                !     SEGMENT JE U KOSOM POLOZAJU
                ! -------------------------------------
                ! Koeficijenti odslikavanja homogenog tla
                C(1:3) = dcmplx(0.d0,0.d0)
                zsi = (z2+z1)/2.d0
                zm1(1) = zsi
                zm2(1) = z2
                zm3(1) = z1
                C(1) = Ft
                if ((isl<=si).and.(isl/=1)) then
                    zm1(2) = 2.d0*HD(isl-1)-zsi
                    zm2(2) = 2.d0*HD(isl-1)-z2
                    zm3(2) = 2.d0*HD(isl-1)-z1
                    C(2) = -Ft*Fr(isl-1)
                else if((isl>=si).and.(si/=1)) then
                    zm1(2) = 2.d0*HD(si-1)-zsi
                    zm2(2) = 2.d0*HD(si-1)-z2
                    zm3(2) = 2.d0*HD(si-1)-z1
                    C(2) = -Ft*Fr(si-1)
                end if
                if ((isl<=si).and.(si/=N)) then
                    zm1(3) = 2.d0*HD(si)-zsi
                    zm2(3) = 2.d0*HD(si)-z2
                    zm3(3) = 2.d0*HD(si)-z1
                    C(3) = Ft*Fr(si)
                else if((isl>=si).and.(isl/=N)) then
                    zm1(3) = 2.d0*HD(isl)-zsi
                    zm2(3) = 2.d0*HD(isl)-z2
                    zm3(3) = 2.d0*HD(isl)-z1
                    C(3) = Ft*Fr(isl)
                end if
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                Rs = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2)
                ggama = gamae/Rs
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j)
                    ! Analiticki dio - ANDI
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                zsi = (zk(iseg)+zp(iseg))/2.d0
                prvi_dio = faktorp*(-ggama)*Fprig*((zt-zsi)/Rs)*integral
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j); z1 = zm3(j)
                    ! Derivacija - DFDz
                    u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
                    v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
                    if (v < pol) v = pol
                    DuDz = (z2-z1)/L(iseg)
                    DvDz = (1.d0/v) * (zt-zsi-u*((z2-z1)/L(iseg)))
                    temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
                    temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
                    DFDu = 1.d0/temp1 - 1.d0/temp2
                    DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
                    DFDz = DFDu*DuDz + DFDv*DvDz
                    integral = integral + C(j)*DFDz
                end do
                drugi_dio = faktorp*Fprig*integral
            end if
        Fiz = Fiz + (prvi_dio + drugi_dio)
        END DO
    ELSE
        ! ............................................................................
        !          OPCENITO VISESLOJNO SREDSTVO -> ZRAK + VISESLOJNO TLO
        ! ............................................................................
        Fiz = dcmplx(0.d0,0.d0)
        ngt1 = 5
        call kGint(ngt1,ug1,Hug1)
        DO iseg = 1,Nss
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
                prvi_dio = dcmplx(0.d0,0.d0)
                drugi_dio = dcmplx(0.d0,0.d0)
                treci_dio = dcmplx(0.d0,0.d0)
                cetvrti_dio = dcmplx(0.d0,0.d0)
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
                ! ------------------------------------------------------
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                Rs = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2)
                ggama = gamae/Rs
                integral = dcmplx(0.d0,0.d0)
                do j = 1,33
                    zsi = zm(j); z2 = zm(j)
                    ! Analiticki dio
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                zsi = (zk(iseg)+zp(iseg))/2.d0
                prvi_dio = faktorp*(Ip(iseg)/L(iseg))*(-ggama)*Fprig*((zt-zsi)/Rs)*integral
                ! ------------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                do j = 1,33
                    zsi = zm(j); z2 = zm(j); z1 = zm(j)
                    ! Derivacija - DFDz
                    u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
                    v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
                    if (v < pol) v = pol
                    DuDz = (z2-z1)/L(iseg)
                    DvDz = (1.d0/v) * (zt-zsi-u*((z2-z1)/L(iseg)))
                    temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
                    temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
                    DFDu = 1.d0/temp1 - 1.d0/temp2
                    DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
                    DFDz = DFDu*DuDz + DFDv*DvDz
                    integral = integral + C(j)*DFDz
                end do
                drugi_dio = faktorp*(Ip(iseg)/L(iseg))*Fprig*integral
            else
                ! -----------------------------------------
                !    SEGMENT JE U KOSOM POLOZAJU
                ! -----------------------------------------
                prvi_dio = dcmplx(0.d0,0.d0)
                drugi_dio = dcmplx(0.d0,0.d0)
                ! Koeficijenti odslikavanja homogenog tla
                C(1:33) = dcmplx(0.d0,0.d0)
                zsi = (z2+z1)/2.d0
                zm1(1) = zsi
                zm2(1) = z2
                zm3(1) = z1
                C(1) = Ft
                if ((isl<=si).and.(isl/=1)) then
                    zm1(2) = 2.d0*HD(isl-1)-zsi
                    zm2(2) = 2.d0*HD(isl-1)-z2
                    zm3(2) = 2.d0*HD(isl-1)-z1
                    C(2) = -Ft*Fr(isl-1)
                else if((isl>=si).and.(si/=1)) then
                    zm1(2) = 2.d0*HD(si-1)-zsi
                    zm2(2) = 2.d0*HD(si-1)-z2
                    zm3(2) = 2.d0*HD(si-1)-z1
                    C(2) = -Ft*Fr(si-1)
                end if
                if ((isl<=si).and.(si/=N)) then
                    zm1(3) = 2.d0*HD(si)-zsi
                    zm2(3) = 2.d0*HD(si)-z2
                    zm3(3) = 2.d0*HD(si)-z1
                    C(3) = Ft*Fr(si)
                else if((isl>=si).and.(isl/=N)) then
                    zm1(3) = 2.d0*HD(isl)-zsi
                    zm2(3) = 2.d0*HD(isl)-z2
                    zm3(3) = 2.d0*HD(isl)-z1
                    C(3) = Ft*Fr(isl)
                end if
                ! ----------------------------------------------------
                xsi = (x2+x1)/2.d0; xsj = xt
                ysi = (y2+y1)/2.d0; ysj = yt
                zsi = (z2+z1)/2.d0; zsj = zt
                call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
                Fprig = zexp(-gamae)
                Rs = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2)
                ggama = gamae/Rs
                integral = dcmplx(0.d0,0.d0)
                zt1 = zsi; zt2 = z2; zt3 = z1
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j)
                    ! Analiticki dio - ANDI
                    call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)
                    integral = integral + C(j)*andi
                end do
                zsi = (zk(iseg)+zp(iseg))/2.d0
                prvi_dio = faktorp*(Ip(iseg)/L(iseg))*(-ggama)*Fprig*((zt-zsi)/Rs)*integral
                ! ----------------------------------------------------
                integral = dcmplx(0.d0,0.d0)
                do j = 1,3
                    zsi = zm1(j); z2 = zm2(j); z1 = zm3(j)
                    ! Derivacija - DFDz
                    u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
                    v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
                    if (v < pol) v = pol
                    DuDz = (z2-z1)/L(iseg)
                    DvDz = (1.d0/v) * (zt-zsi-u*((z2-z1)/L(iseg)))
                    temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
                    temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
                    DFDu = 1.d0/temp1 - 1.d0/temp2
                    DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
                    DFDz = DFDu*DuDz + DFDv*DvDz
                    integral = integral + C(j)*DFDz
                end do
                drugi_dio = faktorp*(Ip(iseg)/L(iseg))*Fprig*integral
                ! *********************************************************************
                treci_dio = dcmplx(0.d0,0.d0)
                if (isl/=1) then
                    integral = dcmplx(0.d0,0.d0)
                    z2 = zt2; zsi = zt1; z1 = zt3
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
                            faktore = (zt-HD(isl-1)+teta*eta(ii))/(korjen**3)
                            integral = integral + alfa(ii)*Hug1(ig)*faktore
                        end do
                    end do
                    treci_dio = - faktorp * Fprig * integral * Ip(iseg)
                end if
                ! *********************************************************************
                cetvrti_dio = dcmplx(0.d0,0.d0)
                if (isl/=N) then
                    integral = dcmplx(0.d0,0.d0)
                    z2 = zt2; zsi = zt1; z1 = zt3
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
                            faktore = (zt-HD(isl)-hi*eta(ii))/(korjen**3)
                            integral = integral + beta(ii)*Hug1(ig)*faktore
                        end do
                    end do
                    cetvrti_dio = - faktorp * Fprig * integral * Ip(iseg)
                end if
            end if
        Fiz = Fiz + (prvi_dio + drugi_dio + treci_dio + cetvrti_dio)
        END DO
    END IF
    
    return
end subroutine
