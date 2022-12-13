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

! Module sa definicijom vanjskih funkcija i subroutine-a
module funkcije
    implicit none
    
    interface
        ! Proracun modula i kuta kompleksnog broja, vodeci pritom
        ! racuna o orjentaciji kuta prema kvadrantima koord. sustava
        subroutine zmk(zz,zmod,akut)
            complex(8) zz
            real(8) zmod,akut
        end subroutine

        ! Funkcija: dops - dops.f90
        real(8) function dops(u,v)
            real(8) u,v
        end function
        
        ! Funkcija: doks - doks.f90
        real(8) function doks(x,z,rcos,delta)
            real(8) x,z,rcos
            real(8) delta
        end function
        
        ! Subroutine: lkkp - lkkp.for
        SUBROUTINE LKKP(X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,XL1,XL2,ZL1,ZL2,DELTA)
            real(8) x1,x2,x3,x4
            real(8) y1,y2,y3,y4
            real(8) z1,z2,z3,z4
            real(8) xl1,xl2,zl1,zl2
            real(8) delta
        END SUBROUTINE
        
        ! Funkcija: vektor_A - faktori_transmisije.f90
        complex(8) function vektor_A(iso,isl,n,Fr)
            integer iso,isl,n
            complex(8),dimension(:) :: Fr
        end function
        
        ! Vektor refleksija slojeva - faktori_refleksije.f90
        subroutine vektor_F(n,kapa,F)
            integer n
            complex(8),dimension(:) :: kapa,F
        end subroutine
        
        ! Matrica [E] - Pseudoinverzna_matrica.for
        subroutine invmat(E,m,n)
            integer m,n
            real(8),dimension(:,:) :: E
        end subroutine
        
        ! Vektor kapa - proracun_kapa.f90
        subroutine proracun_kapa(n,f,ro,epsr,kapa)
            integer n
            real(8) f
            real(8),dimension(:) :: ro,epsr
            complex(8),dimension(:) :: kapa
        end subroutine

        ! Vektor gama - proracun_gama.f90
        subroutine proracun_gama(f,n,kapa,gama)
            integer n
            real(8) f
            complex(8),dimension(:) :: kapa,gama
        end subroutine
        
        ! Koordinate (donje granicene plohe) slojeva - debljina_slojeva.f90
        subroutine vektor_H(n,h_sloja,HD)
            integer n
            real(8),dimension(:) :: h_sloja
            real(8),dimension(:) :: HD
        end subroutine
        
        ! Parametri teta i hi (za odslikavanje)
        subroutine parametri_teta_hi(iso,isl,d,HD,h_sloj,n,teta,hi)
            integer iso,isl
            real(8) d
            real(8),dimension(:) :: HD,h_sloj
            integer n
            real(8) teta,hi
        end subroutine
        
        ! Vektor: vektor_fi - vektor_fi.f90
        subroutine vektor_fi(n,k,d,lamda,iso,isl,HD,Ft,Fr,H,fi)
            integer n,k
            real(8),dimension(:) :: lamda
            real(8) d
            integer iso,isl
            complex(8),dimension(:) :: Fr
            complex(8) Ft
            real(8),dimension(:) :: HD,H
            complex(8) fi
        end subroutine
        
        ! Vektor: vektor_Gi - vektor_gi.f90
        subroutine vektor_Gi(n,k,d,lamda,iso,isl,HD,Ft,Fr,H,Gi)
            integer n,k
            real(8),dimension(:) :: lamda
            real(8) d
            integer iso,isl
            complex(8),dimension(:) :: Fr
            complex(8) Ft
            real(8),dimension(:) :: HD,H
            complex(8) Gi
        end subroutine
        
        !Vektori: qi, ai i bi - vektori_qi_ai_bi.f90
        subroutine vektori_q_a_b(n,kapa,qi,ai,bi)
            integer n
            complex(8),dimension(:) :: kapa
            complex(8),dimension(:) :: qi,ai,bi
        end subroutine
        
        !Vektor: tetan - vektor_teta_n.f90
        subroutine vektor_teta_n(n,dd,H,HD,lamda,kapa,Fr,isl,iso,tetan)
            integer n
            real(8) dd
            real(8),dimension(:) :: H,HD,lamda
            complex(8),dimension(:) :: kapa,Fr
            integer isl,iso
            complex(8),dimension(:) :: tetan
        end subroutine
        
        !Vektori: teta_k i hi_k - vektor_teta.f90
        subroutine vektor_teta_hi(n,d,tetan,kapa,lamda,isl,iso,Fr,H,HD,teta,hi)
            integer n
            real(8) d
            complex(8),dimension(:) :: tetan
            complex(8),dimension(:) :: kapa
            real(8),dimension(:) :: lamda
            integer isl,iso
            complex(8),dimension(:) :: Fr
            real(8),dimension(:) :: H, HD
            complex(8),dimension(:) :: teta,hi
        end subroutine
        
        ! Proracun jedinicne unutarnje impedancije segmenata uzemljivaca
        subroutine proracun_zunjv(Ns,f,mirv,sigmav,r0,zunjv)
            integer Ns
            real(8) f
            real(8),dimension(:) :: mirv,sigmav,r0
            complex(8),dimension(:) :: zunjv
        end subroutine

        ! Proracun ekvivalentne valne konstante gamae
        subroutine gamaekv(n,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
            integer n
            complex(8),dimension(:) :: gama
            real(8),dimension(:) :: HD
            real(8) xsi,ysi,zsi,xsj,ysj,zsj
            complex(8) gamae
        end subroutine

        subroutine gaussk(ngt,ug,Hug)
            integer ngt
            real(8),dimension(:) :: ug,Hug
        end subroutine

        subroutine kGint(ngt1,ug1,Hug1)
            integer ngt1
            real(8),dimension(:) :: ug1,Hug1
        end subroutine

        ! Proracun koeficijenata alfa - alfa.f90
        subroutine alfauz(teta,w,eta,h,n,d,Fr,HD,s,kapa,sloj,alfa)
            real(8) teta
            real(8),dimension(:) :: eta,w
            real(8),dimension(:) :: h,HD
            integer n
            real(8) d
            complex(8),dimension(:) :: Fr
            integer s,sloj
            complex(8),dimension(:) :: kapa
            complex(8),dimension(:) :: alfa
        end subroutine

        ! Proracun koeficijenata beta - beta.f90
        subroutine betauz(hi,w,eta,h,n,d,Fr,HD,s,kapa,sloj,beta)
            real(8) hi
            real(8),dimension(:) :: eta,w
            real(8),dimension(:) :: h,HD
            integer n
            real(8) d
            complex(8),dimension(:) :: Fr
            integer s,sloj
            complex(8),dimension(:) :: kapa
            complex(8),dimension(:) :: beta
        end subroutine

        ! Proracun matrice admitancija sustava segmenata uzemljivaca
        subroutine matrica(f,L,xp,yp,zp,xk,yk,zk,Ns,Nsi,r0,kapa,gama,sigmav,Zunjs,&
        iso,N,Fr,HD,h,Yuv,Ypv,Yu,Ypp)
            real(8) f
            real(8),dimension(:) :: L
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk
            integer Ns,Nsi,N
            real(8),dimension(:) :: r0,sigmav
            complex(8),dimension(:) :: kapa,gama,Zunjs
            integer,dimension(:) :: iso
            complex(8),dimension(:) :: Fr
            real(8),dimension(:) :: HD,h
            complex(8),dimension(:,:) :: Yuv,Ypv,Yu,Ypp
        end subroutine

        ! Proracun potencijala u jednoj tocki uz retardaciju potencijala
        subroutine potencijal_tocke_freq(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,&
        Fr,iso,r0,L,Ip,HD,h,kapa,gama,potencijal)
            integer N,Ns,Nsi
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Fr,Ip
            real(8),dimension(:) :: HD,h,L
            complex(8),dimension(:) :: kapa,gama
            complex(8) potencijal
        end subroutine

        ! Proracun potencijala u jednoj tocki - pri jednoj frekvenciji,
        ! za potrebe frekvencijske analize: FFT - IFFT
        subroutine single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,&
            r0,f,L,HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,&
            potenc_freq)
            integer N,Ns,Nsi,Nc,Ni,sloj
            integer,dimension(:) :: iveze,icv,iso
            real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Ii
            real(8) xt,yt,zt,f
            complex(8) potenc_freq
        end subroutine

        ! Proracun Y(s) - pri jednoj frekvenciji,
        ! za potrebe vector fitting analize
        subroutine Y_single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,&
            r0,f,L,HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
            integer N,Ns,Nsi,Nc,Ni,sloj
            integer,dimension(:) :: iveze,icv,iso
            real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Ii
            real(8) xt,yt,zt,f
            complex(8) potenc_freq
        end subroutine

        ! Proracun Fix komponente za proracun Ex komponente elektricnog polja
        subroutine Epolje_Fix(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,&
            L,Ip,HD,h,kapa,gama,Fix)
            integer N,Ns,Nsi
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Fr,Ip
            real(8),dimension(:) :: HD,h,L
            complex(8),dimension(:) :: kapa,gama
            complex(8) Fix
        end subroutine

        ! Proracun Fiy komponente za proracun Ex komponente elektricnog polja
        subroutine Epolje_Fiy(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,&
            L,Ip,HD,h,kapa,gama,Fiy)
            integer N,Ns,Nsi
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Fr,Ip
            real(8),dimension(:) :: HD,h,L
            complex(8),dimension(:) :: kapa,gama
            complex(8) Fiy
        end subroutine

        ! Proracun Fiz komponente za proracun Ex komponente elektricnog polja
        subroutine Epolje_Fiz(N,Ns,Nsi,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,&
            L,Ip,HD,h,kapa,gama,Fiz)
            integer N,Ns,Nsi
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Fr,Ip
            real(8),dimension(:) :: HD,h,L
            complex(8),dimension(:) :: kapa,gama
            complex(8) Fiz
        end subroutine

        ! Proracun Ax, Ay i Az komponente za proracun elektricnog polja
        subroutine Epolje_Ax_Ay_Az(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,&
        r0,L,Iu,HD,gama,f,Ax,Ay,Az)
            integer N,Ns
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt,f
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Iu
            real(8),dimension(:) :: HD,L
            complex(8),dimension(:) :: gama
            complex(8) Ax,Ay,Az
        end subroutine
        
        ! Proracun Bx komponente magnetskog polja
        subroutine Hpolje_Bx(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bx)
            integer N,Ns
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt,f
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Iu
            real(8),dimension(:) :: HD,L
            complex(8),dimension(:) :: gama
            complex(8) Bx
        end subroutine

        ! Proracun By komponente magnetskog polja
        subroutine Hpolje_By(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,By)
            integer N,Ns
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt,f
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Iu
            real(8),dimension(:) :: HD,L
            complex(8),dimension(:) :: gama
            complex(8) By
        end subroutine

        ! Proracun Bz komponente magnetskog polja
        subroutine Hpolje_Bz(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bz)
            integer N,Ns
            integer isl
            integer,dimension(:) :: iso
            real(8) :: xt,yt,zt,f
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            complex(8),dimension(:) :: Iu
            real(8),dimension(:) :: HD,L
            complex(8),dimension(:) :: gama
            complex(8) Bz
        end subroutine
        
        ! Proracun ukupnog elektricnog polja u jednoj tocki za jednu frekvenciju
        ! za potrebe tranzijentne (FFT-MoM-IFFT) analize
        subroutine single_frequency_elpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,&
            h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
            integer N,Ns,Nsi,Nc,Ni,sloj,izbor
            integer,dimension(:) :: iveze,icv,iso
            real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            real(8) xt,yt,zt,f
            complex(8),dimension(:) :: Ii
            complex(8) potenc_freq  !Recikliranje varijable
        end subroutine

        ! Proracun ukupnog magnetskog polja u jednoj tocki za jednu frekvenciju
        ! za potrebe tranzijentne (FFT-MoM-IFFT) analize
        subroutine single_frequency_magpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,&
            h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
            integer N,Ns,Nsi,Nc,Ni,sloj,izbor
            integer,dimension(:) :: iveze,icv,iso
            real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
            real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
            real(8) xt,yt,zt,f
            complex(8),dimension(:) :: Ii
            complex(8) potenc_freq   !Recikliranje varijable
        end subroutine

        ! Lapack funkcija
        INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
          CHARACTER*( * )    NAME, OPTS
          INTEGER            ISPEC, N1, N2, N3, N4
        END FUNCTION
        
    end interface
end module
