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

!================================================================================
!   Subroutina koja odreduje vrijednost vektora teta   u 25 uzorkovanih tocaka
!                                                   n
!================================================================================

!   Naime, potrebno je odrediti sljedece vrijednosti:

!           teta  (u ) , i = 1,2,...,25.
!               n   i

!   Dakle, ovdje se odreduje n-ta vrijednost vektora teta za svih 25
!   uzorkovanih vrijednosti.

subroutine vektor_teta_n(n,dd,H,HD,lamda,kapa,Fr,isl,iso,tetan)
    use funkcije
    implicit none
!   Input
    integer n
    real(8) dd
    real(8),dimension(:) :: H,HD,lamda
    complex(8),dimension(:),allocatable :: qi,ai,bi
    complex(8),dimension(:) :: kapa
    complex(8),dimension(:) :: Fr
    complex(8) fi,Gi
    complex(8) Ft
    integer isl,iso
!   Output
    complex(8),dimension(:) :: tetan
!   Local variables
    complex(8) A,B,C,D
    integer k,i,ii
    complex(8) faktor
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)

!   Opis lokalnih varijabli:
!       n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!       dd - dubina ukopavanja izvora (negativna vrijednost ako je izvor u zraku).
!            Ova vrijednost jos je oznacivana i sa slovom d, [m].
!       H - vektor debljina svih slojeva, [m]. (h_sloj)
!       HD - vektor z koordinata donjih granicnih ploha svih slojeva, [m]. (H)
!       lamda - vektor sljedecih vrijednosti: lamda(i) = 1 / u(i), gdje su u(i)
!               uzorkovane tocke.
!       kapa - vektor kompleksnih specificnih elektricnih vodljivosti svih slojeva.
!       Fr - vektor faktora refleksije svih slojeva (F).
!       isl - sloj u kojem se racuna raspodjela potencijala
!       iso - sloj u kojem se nalazi tockasti izvor struje
!       tetan - vektor vrijednosti teta  u 25 uzorkovanih tocaka
!                                      n
!   -------------------------------------------------------------------------------

    allocate(qi(n-1))
    allocate(ai(n-1))
    allocate(bi(n-1))

    call vektori_q_a_b(n,kapa,qi,ai,bi)

    do k = 1,25
        A = one
        B = -qi(1)
        C = -qi(1) * dexp(-lamda(k)*H(2))
        ii = 1
        Ft = vektor_A(iso,ii,n,Fr)
        call vektor_fi(n,k,dd,lamda,iso,ii,HD,Ft,Fr,H,fi)
        D = fi
        do i = 1,n-1
            if (i<(n-1)) then
                A = ai(i) - B
                B = -bi(i) * dexp(-lamda(k)*H(i+1)) - C
                C = dcmplx(0.0,0.0)
                Ft = vektor_A(iso,i,n,Fr)
                call vektor_Gi(n,k,dd,lamda,iso,i,HD,Ft,Fr,H,Gi)
                D = Gi - D
                faktor = dexp(-lamda(k)*H(i+1)) / A
                B = faktor * B
                D = faktor * D
                A = one - B
                B = -qi(i+1)
                if(i<(n-2)) then
                    C = -qi(i+1) * dexp(-lamda(k)*H(i+2))
                end if
                ii = i+1
                Ft = vektor_A(iso,ii,n,Fr)
                call vektor_fi(n,k,dd,lamda,iso,ii,HD,Ft,Fr,H,fi)
                D = fi - D
                faktor = one / A
                B = faktor * B
                if(i<(n-2)) then
                    C = faktor * C
                end if
                D = faktor * D
            else if(i==(n-1)) then
                A = ai(i) - B
                Ft = vektor_A(iso,i,n,Fr)
                call vektor_Gi(n,k,dd,lamda,iso,i,HD,Ft,Fr,H,Gi)
                D = Gi - D
            end if
        end do
        tetan(k) = D / A
    end do

    deallocate(qi)
    deallocate(ai)
    deallocate(bi)

    return
end subroutine