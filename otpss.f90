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

! Proracun analitickog dijela vlastite impedancije segmenta
! VLASTITA IMPEDANCIJA SEGMENTA - ANALITICKI DIO

subroutine otpss(duls,v,andi)
	implicit none
	
	real(8) duls,v
	real(8) andi
	! Lokalne varijable
	real(8) r,pc
	
	! Opis varijabli:
	! duls - duljina segmenta koji se razmatra (li)
	! v - radijus razmatranog segmenta (r0)
	! andi (output) - vrijednost analitickog dijela integrala
	!                 za proracun vlastite impedancije segmenta
	! ---------------------------------------------------------
	
	r = dsqrt(duls**2 + v**2)
	pc = (r + duls)/v
	andi = 2.d0*(duls*dlog(pc)-r+v)
	
	return
end subroutine
