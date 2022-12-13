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
!           Subroutina koja racuna velicinu H koja predstavlja z koordinatu
!           donje granicne plohe za svaki sloj n-slojnog modela.
!================================================================================

!Velicina H skicirana je na donjoj slici:

!       kapa(1)       h(1)=0       H(1)=0               zrak
!       ----------------------------------------------------
!                   ^           |       |       |       tlo
!       kapa(2)     | h(2)      |  H(2) |       |
!                   v           v       |       |
!       --------------------------------|-------|-----------
!                   ^                   | H(3)  |       
!                   |                   |       |
!       kapa(3)     | h(3)              |       |
!                   v                   v       |
!       ----------------------------------------|-----------
!                                               |
!                               .               | H(n-1)
!                               .               |
!                               .               |
!                                               |
!       ----------------------------------------|-----------
!                   ^                           |
!                   |                           |
!       kapa(n-1)   | h(n-1)                    |
!                   v                           v
!       ----------------------------------------------------

!       kapa(n)       h(n)=0                      H(n)=0

!pri cemu su:
!      n - ukupni broj slojeva modela (zrak + viseslojno tlo).
!   h(i) - debljine pojedinih slojeva n-slojnog modela, [m]. Vektor h ima
!          ukupno n elemenata. Ove vrijednosti se zadaju s ulaznim podacima.
!   H(i) - z koordinata donje grani�e plohe pojedinog i-tog sloja, [m]. Ovaj
!          vektor ima ukupno n elemenata, pri cemu je H(1) = 0. H(2) je z
!          koordinata prvog sloja, H(3) drugog itd. H(n) je beskonacno ali
!          se postavlja na vrijednost nula.
!================================================================================

subroutine vektor_H(n,h_sloja,HD)
    implicit none

!   Input
    integer n
    real(8),dimension(:) :: h_sloja
!   Output
    real(8),dimension(:) :: HD
!   Local variable
    integer i


    !Tijelo funkcije
    HD(1) = 0.0
    do i = 2,n-1
        HD(i) = HD(i-1) + h_sloja(i)
    end do
    HD(n) = 0.0

    return
end subroutine