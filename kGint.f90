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

! Subroutinea koja sadrzi koeficijente 5 Gaussovih tocaka
! za odslikavanje kosog segmenta u viseslojnom sredstvu

subroutine kGint(ngt1,ug1,Hug1)
    implicit none
    
    integer ngt1
    real(8),dimension(:) :: ug1,Hug1
    
    
    ug1(1) = 0.4530899
    ug1(2) = -ug1(1)
    ug1(3) = 0.26923465
    ug1(4) = -ug1(3)
    ug1(5) = 0.d0
    
    Hug1(1) = 0.11846345
    Hug1(2) = Hug1(1)
    Hug1(3) = 0.23931435
    Hug1(4) = Hug1(3)
    Hug1(5) = 0.28444445
    
    return
end subroutine
