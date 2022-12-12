! -------------------------------------------------------------------------------
!	 Copyright (C) 2006. GPL - General Public Licence
!	 Author: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)

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

! Definicija koeficijenata Gauss-ove integracije potrebnih pri
! izracunavanju numerickog dijela integrala (za impedancije).

subroutine gaussk(ngt,ug,Hug)
	implicit none
	
	integer ngt
	real(8),dimension(:) :: ug,Hug
	
	
	ug(1) = 0.949107912
	ug(2) = 0.741531186
	ug(3) = 0.405845151
	ug(4) = 0.d0
	ug(5) = -ug(3)
	ug(6) = -ug(2)
	ug(7) = -ug(1)
	
	Hug(1) = 0.129484966
	Hug(2) = 0.279705391
	Hug(3) = 0.381830051
	Hug(4) = 0.417959184
	Hug(5) = Hug(3)
	Hug(6) = Hug(2)
	Hug(7) = Hug(1)
	
	return
end subroutine
