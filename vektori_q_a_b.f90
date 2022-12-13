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
!   Subroutina koja odreduje vrijednosti vektora q, a i b potrebnih za
!   odredjivanje vektora Gi i fi. Vidjeti Nada Magistarski rad str. 20 za
!   odredjivanje vektora q, te str. 22 za vektore a i b                                         
!================================================================================

!   Vektor q je definiran na sljedeci nacin:

!                        kapa
!                            i
!               q   =  -----------
!                i       kapa
!                            i+1

!   pri cemu je kapa - kompleksna specificna elektricna vodljivost
!   sredstva, [S/m].

!   Vaktor a definiran je na sljedeci nacin:

!                        1 - q
!                             i
!               a   =  -----------
!                i         2

!   pri cemu je q(i) prethodno definiran.

!   Vektor b je definiran na sljedeci nacin:

!                        1 + q
!                             i
!               b   =  -----------
!                i         2

!   Svi vektori imaju ukupno n-1 elemenata!

subroutine vektori_q_a_b(n,kapa,qi,ai,bi)
    implicit none

!   Input
    integer n
    complex(8),dimension(:) :: kapa
!   Output variables
    complex(8),dimension(:) :: qi,ai,bi
!   Local variable
    complex(8),parameter :: one = dcmplx(1.0,0.0)
    integer i

    do i = 1,n-1
        qi(i) = (kapa(i))/(kapa(i+1))
        ai(i) = (one - qi(i)) / 2.d0
        bi(i) = (one + qi(i)) / 2.d0
    end do


    return
end subroutine