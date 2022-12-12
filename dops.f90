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

! Funkcija koja racuna osnovni izraz za medjusobnu impedanciju
! (otpor) izmedju medjusobno paralelnih segmenata. Poziva se u
! subroutini: otpps

real(8) function dops(u,v)
	implicit none
	
	real(8) u,v
	real(8) r
	
	
	r = dsqrt(u**2 + v**2)
	if (u >= 0.d0) then
		dops = u*dlog(r+u)-r
	else
		dops = u*dlog(v**2/(r-u))-r
	end if
	
	return
end function