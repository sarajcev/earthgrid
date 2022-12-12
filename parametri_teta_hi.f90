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

!================================================================================
!			Funkcija koja racuna parametre teta i hi potrebne za odslikavanje
!			tockastog izvora u viseslojnom tlu. Ovi su faktori odredjeni 
!			posebnim proracunom, koji nije prikazan u radu: Nada magistarski.
!================================================================================

subroutine parametri_teta_hi(s,sloj,d,H,h_sloj,n,teta,hi)
	implicit none

!	Input variables
	integer s,sloj
	real(8) d
	real(8),dimension(:) :: H,h_sloj
	integer n
!	Output variables
	real(8) teta,hi

!	Objasnjenje varijabli:
!		s      - sloj u kojem se nalzi tockasti izvor struje (iso)
!		sloj   - sloj za koji se racunaju trazeni parametri teta i hi (isl)
!		d      - dubina ukopavanja tockastog izvora struje, [m]
!		H      - vektor koji sadrzi koordinate donjih granicnih ploha svih slojeva
!			     viseslojnog modela tla u [m], (ukupno n elemenata) - HD
!		h_sloj - vektor koji sadrzi debljine svih slojeva u [m] viseslojnog 
!				 modela tla (unupno n elemenata) - H.
!		n      - broj slojeva viseslojnog modela (zrak + viseslojno tlo)

!	Local variables
	real(8),dimension(2) :: dvve
	real(8),dimension(:),allocatable :: hp
	integer t

!	Interface valjske funkcije
	interface
		real(8) function minimum(hp,t)
			integer t
			real(8),dimension(:) :: hp
		end function
	end interface


!	Odredjivanje parametara teta i hi
!	---------------------------------
	if((sloj>=3).and.(sloj<s).and.(s<=(n-2))) then
		!Faktor teta
		t = (s-1) - (sloj-1) + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj-1:s-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(s)-d)
		teta = d - H(sloj-1) + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = s - sloj + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj:s)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(s+1)-d)
		hi = d - H(sloj) + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((sloj==2).and.(sloj<s).and.(s<=(n-2))) then
		!Faktor teta
		t = s - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:s-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(s)-d)
		teta = d + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = s - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:s)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(s+1)-d)
		hi = d - H(2) + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((sloj==1).and.(sloj<s).and.(s<=(n-2))) then
		!Faktor teta
		teta = 0.d0
		!Faktor hi
		t = s - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:s)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(s+1)-d)
		hi = d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((sloj>=3).and.(sloj<s).and.(s==(n-1))) then
		!Faktor teta
		t = (n-2) - (sloj-1) + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj-1:n-2)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(n-1)-d)
		teta = d - H(sloj-1) + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = (n-1) - sloj + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj:n-1)
		hi = d - H(sloj) + minimum(hp,t)
		deallocate(hp)

	else if((sloj==2).and.(sloj<s).and.(s==(n-1))) then
		!Faktor teta
		t = (n-2) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-2)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (H(n-1)-d)
		teta = d + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		hi = d - H(2) + minimum(hp,t)
		deallocate(hp)

	else if((sloj==1).and.(sloj<s).and.(s==(n-1))) then
		!Faktor teta
		teta = 0.d0
		!Faktor hi
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		hi = d + minimum(hp,t)
		deallocate(hp)

	else if((sloj>=3).and.(sloj<s).and.(s==n)) then
		!Faktor teta
		t = (n-1) - (sloj-1) + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj-1:n-1)
		teta = d - H(sloj-1) + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (n-1) - sloj + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(sloj:n-1)
		hi = d - H(sloj) + minimum(hp,t)
		deallocate(hp)

	else if((sloj==2).and.(sloj<s).and.(s==n)) then
		!Faktor teta
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		teta = d + minimum(hp,t)
		!Faktor hi
		hi = d - H(2) + minimum(hp,t)
		deallocate(hp)

	else if((sloj==1).and.(sloj<s).and.(s==n)) then
		!Faktor teta
		teta = 0.d0
		!Faktor hi
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		hi = d + minimum(hp,t)
		deallocate(hp)

	else if((sloj>=4).and.(sloj==s).and.(s<=(n-2))) then
		!Faktor teta
		dvve(1) = H(s) - d + h_sloj(s)
		dvve(2) = d - H(s-2) + h_sloj(s-1)
		teta = dmin1(dvve(1),dvve(2))
		!Faktor hi
		dvve(1) = d - H(s-1) + h_sloj(s)
		dvve(2) = H(s+1) - d + h_sloj(s+1)
		hi = dmin1(dvve(1),dvve(2))

	else if((sloj==3).and.(sloj==s).and.(s<=(n-2))) then
		!Faktor teta
		dvve(1) = H(3) - d + h_sloj(3)
		dvve(2) = d + h_sloj(2)
		teta = dmin1(dvve(1),dvve(2))
		!Faktor hi
		dvve(1) = d - H(2) + h_sloj(3)
		dvve(2) = H(4) - d + h_sloj(4)
		hi = dmin1(dvve(1),dvve(2))

	else if((sloj==2).and.(sloj==s).and.(s<=(n-2))) then
		!Faktor teta
		teta = 2.d0 * h_sloj(2) - d
		!Faktor hi
		dvve(1) = d + h_sloj(2)
		dvve(2) = H(3) - d + h_sloj(3)
		hi = dmin1(dvve(1),dvve(2))

	else if((sloj==1).and.(sloj==s).and.(s<=(n-2))) then
		!Faktor teta
		teta = 0.d0
		!Faktor hi
		hi = 2.d0 * h_sloj(2) - d

	else if((sloj>=4).and.(sloj==s).and.(s==(n-1))) then
		!Faktor teta
		dvve(1) = H(n-1) - d + h_sloj(n-1)
		dvve(2) = d - H(n-3) + h_sloj(n-2)
		teta = dmin1(dvve(1),dvve(2))
		!Faktor hi
		hi = d - H(n-2) + h_sloj(n-1)

	else if((sloj==3).and.(sloj==s).and.(s==(n-1))) then
		!Faktor teta
		dvve(1) = H(3) - d + h_sloj(3)
		dvve(2) = d + h_sloj(2)
		teta = dmin1(dvve(1),dvve(2))
		!Faktor hi
		hi = d - H(2) + h_sloj(3)

	else if((sloj==2).and.(sloj==s).and.(s==(n-1))) then
		!Faktor teta
		teta = 2.d0 * h_sloj(2) - d
		!Faktor hi
		hi = d + h_sloj(2)

	else if((sloj>=4).and.(sloj==s).and.(s==n)) then
		!Faktor teta
		teta = d - H(n-2) + h_sloj(n-1)
		!Faktor hi
		hi = 0.d0

	else if((sloj==3).and.(sloj==s).and.(s==n)) then
		!Faktor teta
		teta = d + h_sloj(2)
		!Faktor hi
		hi = 0.d0

	else if((s>=4).and.(s<sloj).and.(sloj<=(n-2))) then
		!Faktor teta
		t = sloj - s + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(s:sloj)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(s-2))
		teta = H(sloj-1) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = (sloj+1) - (s+1) + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(s+1:sloj+1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(s-1))
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==3).and.(s<sloj).and.(sloj<=(n-2))) then
		!Faktor teta
		t = sloj - 3 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(3:sloj)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (sloj+1) - 4 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(4:sloj+1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(2))
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==2).and.(s<sloj).and.(sloj<=(n-2))) then
		!Faktor teta
		t = sloj - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:sloj)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (sloj+1) - 3 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(3:sloj+1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * d
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==1).and.(s<sloj).and.(sloj<=(n-2))) then
		!Faktor teta
		t = sloj - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:sloj)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (sloj+1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:sloj+1)
		hi = H(sloj) - d + minimum(hp,t)
		deallocate(hp)

	else if((s>=4).and.(s<sloj).and.(sloj==(n-1))) then
		!Faktor teta
		t = (n-1) - s + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(s:n-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(s-2))
		teta = H(sloj-1) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		t = (n-1) - (s+1) + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(s+1:n-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(s-1))
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==3).and.(s<sloj).and.(sloj==(n-1))) then
		!Faktor teta
		t = (n-1) - 3 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(3:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (n-1) - 4 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(4:n-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(2))
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==2).and.(s<sloj).and.(sloj==(n-1))) then
		!Faktor teta
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		t = (n-1) - 3 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(3:n-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * d
		hi = H(sloj) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)

	else if((s==1).and.(s<sloj).and.(sloj==(n-1))) then
		!Faktor teta
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		!Faktor hi
		hi = H(sloj) - d + minimum(hp,t)
		deallocate(hp)

	else if((s>=4).and.(s<sloj).and.(sloj==n)) then
		!Faktor teta
		t = (n-1) - s + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(s:n-1)
		dvve(1) = minimum(hp,t)
		dvve(2) = 2.d0 * (d-H(s-2))
		teta = H(sloj-1) - d + dmin1(dvve(1),dvve(2))
		deallocate(hp)
		!Faktor hi
		hi = 0.d0

	else if((s==3).and.(s<sloj).and.(sloj==n)) then
		!Faktor teta
		t = (n-1) - 3 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(3:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		hi = 0.d0

	else if((s==2).and.(s<sloj).and.(sloj==n)) then
		!Faktor teta
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		hi = 0.d0

	else if((s==1).and.(s<sloj).and.(sloj==n)) then
		!Faktor teta
		t = (n-1) - 2 + 1
		allocate(hp(t))
		hp(:) = 2.d0 * h_sloj(2:n-1)
		teta = H(sloj-1) - d + minimum(hp,t)
		deallocate(hp)
		!Faktor hi
		hi = 0.d0
	end if

	if((teta>0.0).and.(teta<1.0)) then
		teta = 1.d0
	end if

	if((hi>0.0).and.(hi<1.0)) then
		hi = 1.d0
	end if

	return
end subroutine


!Funkcija za odredjivanje minimuma vektora hp
real(8) function minimum(hp,t)
	implicit none
	integer t
	real(8),dimension(:) :: hp
	real(8) min
	integer i


	min = hp(1)
	do i = 2,t
		if (hp(i)<min) min = hp(i)
	end do

	minimum = min

end function