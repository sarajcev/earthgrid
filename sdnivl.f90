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

! Proracun dvostruke numericke Gauss-ove integracije po
! vlastitom segmentu, u postupku proracuna numericog
! dijela integrala
! VLASTITA IMPEDANCIJA - NUMERICKI DIO

subroutine sdnivl(gamae,duls,vt,cgdni)
	use funkcije
	implicit none
	
	complex(8) gamae
	real(8) duls,vt
	complex(8) cgdni
	! Lokalne varijable
	integer :: ngt = 7
	real(8),dimension(:),allocatable :: ug,Hug
	real(8) duls2,u,pds
	complex(8) cpiugt
	real(8) x,r
	complex(8) cfut
	complex(8),parameter :: one = dcmplx(1.d0,0.d0)
	integer ig,ii
	
	! Opis varijabli:
	!  gamae - ekvivalentna valna konstanta
	!  duls - duljina razmatranog segmenta
	!  vt - radijus razmatranog segmenta (r0)
	!  cgdni (outpu) - izracunata vrijednost integrala
	
	duls2 = duls/2.d0
	allocate(ug(ngt))
	allocate(Hug(ngt))
	call GAUSSK(ngt,ug,Hug)
	
	cgdni = dcmplx(0.d0,0.d0)
	do ig = 1,ngt
		u = ug(ig)*duls2
		! Proracun prve Gauss-ove numericke integracije po
		! segmentu (uz podjelu segmenta na dva dijela)
		pds = (duls2+u)/2.d0
		cpiugt = dcmplx(0.d0,0.d0)
		do ii = 1,ngt
			x = pds*(1.d0-ug(ii))
			r = dsqrt(x**2+vt**2)
			cfut = (cdexp(-gamae*r)-one)/r
			cpiugt = cpiugt + Hug(ii)*pds*cfut
		end do
		pds = (duls2-u)/2.d0
		do ii = 1,ngt
			x = pds*(1.d0-ug(ii))
			r = dsqrt(x**2+vt**2)
			cfut = (cdexp(-gamae*r)-one)/r
			cpiugt = cpiugt + Hug(ii)*pds*cfut
		end do
		cgdni = cgdni + cpiugt*Hug(ig)*duls2
	end do
	
	deallocate(ug)
	deallocate(Hug)

	return
end subroutine