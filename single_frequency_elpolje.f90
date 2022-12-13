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

! Proracun tranzijentnog elektricnog polja u jednoj tocki viseslojnog sredstva 
! za jednu (bilo koju) frekvenciju. Ova subroutina sadrzi sve sto je potrebno 
! za proracunati vrijednost tranzijentnog elek. polja u bilo kojoj tocki vise-
! slojnog sredstva za bilo koju frekvenciju. 

! Subroutina je organizirana na sljedeci nacin:
!   1) Proracunavaju se svi frekvencijski ovisni parametri, kao sto su npr.
!      kapa, gama i faktori refleksije svih slojeva viseslojnog modela, te
!      jedinicne unutarnje impedancije svih segmenata uzemljivaca. Za doticne
!      proracune pozivaju se odgovarajuce subroutine.
!   2) Formira se ukupna matrica admitancija uzemljivackog sustava pri danoj
!      frekvenciji. Racunaju se matrice [Yuv], [Ypv] i ukupna matrica [Y]. Za
!      to se poziva subroutina: matrica.
!   3) Vrsi se asembliranje matrice [Y] u globalnu matricu sustava [Yg] uz
!      pomoc vektora veze (koji je prethodno formiran u glavnom programu).
!   4) Rjesava se globalni sustav jedandzbi [Yg]*{Fig} = {Ig}, pomocu LAPACK  
!      rutine ZGESV. Rjesenjem globalnog sustava postaju poznati potencijali 
!      svih globalnih cvorova sustava uzemljivaca. {Ig} je poznati vektor nari-
!      nutih struja globalnih cvorova.
!   5) Na osnovu vektora veze i poznatih potencijala globalnih cvorova, racunaju
!      se potencijali svih (lokalnih) cvorova sustava uzemljivaca.
!   6) Racunaju se uzduzne struje svih segmenata sustava uzemljivaca. Uzduzne 
!      struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
!      {Iu} = [Yuv]*[Au]*{Fi}, pri cemu vektor {Fi} predstavlja vektor poten-
!      cijala lokalnih cvorova segmenata uzemljivaca. Za ovaj proracun koriste
!      se BLAS 3 rutine: ZGEMM i ZGEMV.
!   7) Racunaju se poprecne struje svih segmenata sustava uzemljivaca. Poprecne
!      struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
!      {Ip} = [Ypv]*[Ap]*{Fi}, pri cemu vektor {Fi} predstavlja vektor poten-
!      cijala lokalnih cvorova segmenata uzemljivaca. Za ovaj proracun koriste
!      se BLAS 3 rutine: ZGEMM i ZGEMV.
!   8) Racuna se tranzijentno totalno (ukupno) elektricno polje u jednoj tocki
!      viseslojnog sredstva, na datoj frekvenciji, pozivom sljedecih subroutina:
!      Epolje_AxAyAz, Epolje_Fix, Epolje_Fiy i Epolje_Fiz.


subroutine single_frequency_elpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,&
    h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
    use funkcije
    implicit none

    ! Input variables:
    integer N,Ns,Nsi,Nc,Ni,sloj,izbor
    integer,dimension(:) :: iveze,icv,iso
    real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
    real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
    real(8) xt,yt,zt,f
    complex(8),dimension(:) :: Ii
    ! Output variables:
    complex(8) potenc_freq   !Recikliranje varijable
    ! Local variables:
    integer Nss
    complex(8),dimension(:),allocatable :: kapa,gama,Fr,zunjv,Fig,Fi,Iu,Ip
    complex(8),dimension(:,:),allocatable :: Yuv,Ypv,Yu,Ypp,Yg,B,Au,Ap,CC
    integer,dimension(:),allocatable :: ipiv
    complex(8) alpha,betha
    complex(8) Fix,Fiy,Fiz,Ax,Ay,Az,Ex,Ey,Ez,El_polje
    character(1) transa,transb
    integer nrhs,info,ig,jg,i,j

    ! ---------------------------------------------------------------------------
    !            PRORACUN FREKVENCIJSKI OVISNIH PARAMETARA PROGRAMA
    ! ---------------------------------------------------------------------------
    Nss = Ns - Nsi

    ! Proracun kompleksnih specificnih vodljivosti slojeva modela
    allocate(kapa(N))
    call proracun_kapa(N,f,ro,epsr,kapa)
    ! Proracun valnih konstanti slojeva modela
    allocate(gama(N))
    call proracun_gama(f,N,kapa,gama)
    ! Proracun faktora refleksije svih slojeva modela
    allocate(Fr(N))
    call vektor_F(N,kapa,Fr)
    ! Proracun jedinicnih unutarnjih impedancija svih segmenata
    allocate(zunjv(Ns))
    call proracun_zunjv(Ns,f,mirv,sigmav,r0,zunjv)

    ! ---------------------------------------------------------------------------
    !             FORMIRANJE UKUPNE LOKALNE MATRICE SUSTAVA [Y] 
    ! ---------------------------------------------------------------------------
    allocate(Yuv(Ns,Ns))
    allocate(Ypv(Nss,Nss))
    allocate(Yu(2*Ns,2*Ns))
    allocate(Ypp(2*Nss,2*Nss))
    ! ................................................................
    call matrica(f,L,xp,yp,zp,xk,yk,zk,Ns,Nsi,r0,kapa,gama,sigmav,zunjv,&
         iso,N,Fr,HD,h,Yuv,Ypv,Yu,Ypp)
    ! ................................................................
    deallocate(zunjv)

    ! ---------------------------------------------------------------------------
    !         ASSEMBLING - FORMIRANJE UKUPNE GLOBALNE MATRICE SUSTAVA [Yg] 
    ! ---------------------------------------------------------------------------
    allocate(Yg(Nc,Nc))
    Yg(:,:) = dcmplx(0.d0,0.d0)
    ! Asembliranje uzduzne matrice admitancija [Yu]
    do j = 1,2*Ns
        jg = iveze(j)
        do i = 1,2*Ns
            ig = iveze(i)
            Yg(ig,jg) = Yg(ig,jg) + Yu(i,j)
        end do
    end do
    deallocate(Yu)
    ! Asembliranje poprecne matrice admitancija [Yp]
    do j = 1,2*Nss
        jg = iveze(j)
        do i = 1,2*Nss
            ig = iveze(i)
            Yg(ig,jg) = Yg(ig,jg) + Ypp(i,j)
        end do
    end do
    deallocate(Ypp)

    ! ---------------------------------------------------------------------------
    !         RJESAVANJE GLOBALNOG SUSTAVA JEDNADZBI [Yg]*{Fig} = {Ig} 
    !         ZA ODREDJIVANJE NEPOZNATIH POTENCIJALA GLOBALNIH CVOROVA
    ! ---------------------------------------------------------------------------
    ! Globalni sustav jednadzbi, complex(8), rjesava se pomocu LAPACK subrutine:
    ! ZGESV. Ova subroutina koristi punu LU dekompoziciju matrice s parcijalnim
    ! pivotiranjem (row interchanges), te Gaussovu elimincaiju.
    nrhs = 1
    allocate(B(Nc,nrhs))
    B(:,nrhs) = dcmplx(0.d0,0.d0)
    do i = 1,Nc
        do j = 1,Ni
            if (icv(j)==i) then
                B(i,nrhs) = Ii(j)
            end if
        end do
    end do 
    allocate(ipiv(Nc))
    ! Poziv Lapack subroutine: ZGESV
    call zgesv(Nc,nrhs,Yg,Nc,ipiv,B,Nc,info)
    deallocate(ipiv)
    allocate(Fig(Nc))
    Fig(:) = B(:,nrhs)
    deallocate(B)
    deallocate(Yg)

    ! ---------------------------------------------------------------------------
    !               PRORACUN POTENCIJALA LOKALNIH CVOROVA
    ! ---------------------------------------------------------------------------
    allocate(Fi(2*Ns))
    do i = 1,2*Ns
        Fi(i) = Fig(iveze(i))
    end do
    deallocate(Fig)

    ! ---------------------------------------------------------------------------
    !               PRORACUN UZDUZNIH STRUJA SVIH SEGMENATA
    ! ---------------------------------------------------------------------------
    ! Uzduzne struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
    ! {Iu} = [Yuv]*[Au]*{Fi}, pri cemu vektor {Fi} predstavlja vektor potencijala
    ! lokalnih cvorova segmenata uzemljivaca. Ovaj proracun odvija se u dva koraka,
    ! kako slijedi:
    !   1) Formira se matrica [Au] koja je pravokutna Au(Ns,2*Ns) i kompleksna
    !      (moze biti i realna).
    !   2) Racuna se izraz {Iu} = [Yuv]*[Au]*{Fi} koristeci BLAS3 i BLAS2 rutine:
    !      ZGEMM i ZGEMV, respektivno. Subroutina ZGEMM mnozi dvije generalne kom-
    !      pleksne matrice (korsteci blokovski algoritam), dok subroutina ZGEMV
    !      mnozi kompleksnu matricu s vektorom.
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
    ! RACUNANJE IZRAZA {Iu} = [Yuv]*[Au]*{Fi} 
    ! ----------------------------------------
    alpha = dcmplx(1.d0,0.d0)
    betha = dcmplx(0.d0,0.d0)
    transa = 'N'
    transb = 'N'
    allocate(CC(Ns,2*Ns))
    call zgemm(transa,transb,Ns,2*Ns,Ns,alpha,Yuv,Ns,Au,Ns,betha,CC,Ns)
    deallocate(Au)
    deallocate(Yuv)
    allocate(Iu(Ns))
    call zgemv(transa,Ns,2*Ns,alpha,CC,Ns,Fi,1,betha,Iu,1)
    deallocate(CC)

    ! ---------------------------------------------------------------------------
    !               PRORACUN POPRECNIH STRUJA SVIH SEGMENATA
    ! ---------------------------------------------------------------------------
    ! Poprecne struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
    ! {Ip} = [Ypv]*[Ap]*{Fi}, pri cemu vektor {Fi} predstavlja vektor potencijala
    ! lokalnih cvorova segmenata uzemljivaca. Ovaj proracun odvija se u dva koraka,
    ! kako slijedi:
    !   1) Formira se matrica [Ap] koja je pravokutna Ap(Ns,2*Ns) i kompleksna (moze
    !      biti i realna).
    !   2) Racuna se izraz {Ip} = [Ypv]*[Ap]*{Fi} koristeci BLAS3 rutine: ZGEMM i
    !      ZGEMV
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
    ! RACUNANJE IZRAZA {Ip} = [Yp]*[Ap]*{Fi} 
    ! ------------------------------------------------------
    alpha = dcmplx(1.d0,0.d0)
    betha = dcmplx(0.d0,0.d0)
    transa = 'N'
    transb = 'N'
    allocate(CC(Nss,2*Nss))
    call zgemm(transa,transb,Nss,2*Nss,Nss,alpha,Ypv,Nss,Ap,Nss,betha,CC,Nss)
    deallocate(Ap)
    deallocate(Ypv)
    allocate(Ip(Nss))
    call zgemv(transa,Nss,2*Nss,alpha,CC,Nss,Fi(1:2*Nss),1,betha,Ip,1)
    deallocate(CC)
    deallocate(Fi)

    ! ---------------------------------------------------------------------------
    !           PRORACUN ELEKTRICNOG POLJA ZA JEDNU FREKVENCIJU 
    ! ---------------------------------------------------------------------------
    select case(izbor)
        case(1)
            ! Poziv rutine za proracun Ex komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fix(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fix)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)
            ! ..........................................................................
            Ex = -Fix - Ax
            !Ex komponenta elektricnog polja
            El_polje = Ex
        case(2)
            ! Poziv rutine za proracun Ey komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fiy(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fiy)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)
            ! ..........................................................................
            Ey = -Fiy - Ay
            !Ey komponenta elektricnog polja
            El_polje = Ey
        case(3)
            ! Poziv rutine za proracun Ez komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fiz(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fiz)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)            
            ! ..........................................................................
            Ez = -Fiz - Az
            !Ez komponenta elektricnog polja
            El_polje = Ez
        case(4)
            ! Poziv rutine za proracun Ex komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fix(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fix)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)
            ! ..........................................................................
            Ex = -Fix - Ax
            ! Poziv rutine za proracun Ey komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fiy(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fiy)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)
            ! ..........................................................................
            Ey = -Fiy - Ay
            ! Poziv rutine za proracun Ez komponente elektricnog polja u odabranoj tocki
            ! ..........................................................................
            call Epolje_Fiz(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
            kapa,gama,Fiz)
            call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
            gama,f,Ax,Ay,Az)            
            ! ..........................................................................
            Ez = -Fiz - Az
            ! Proracun ukupnog/totalnog elektricnog polja
            El_polje = zsqrt(Ex**2 + Ey**2 + Ez**2)
    end select

    potenc_freq = El_polje  !Recikliranje varijable

    deallocate(Iu,Ip)
    deallocate(kapa,gama)
    deallocate(Fr)

    return
end subroutine