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
!			Funkcija koja racuna faktor transmisije A za onaj sloj
!			u kojem se trazi raspodjela potencijala.
!================================================================================

!Faktori transmisije definiraju su na sljedeci nacin:

!		A = (1 + F(s-1))*(1 + F(s-2))*...*(1 + F(i)) ;  za: i < s 

!		A = 1 :  za: i = s

!		A = (1 - F(s))*(1 - F(s+1))*...*(1-F(i-1)) ;  za: i > s 

!pri cemu su:
!	i - trenutna vrijednost sloja za koju se racuna faktor transmisije (isl)
!	s - sloj u kojem se nalazi tockasti izvor struje (iso)
!	F - faktori refleksije, koji su prethodno odredeni. Vektor F je kompleksan
!		i sadrzi n elemenata, pri cemu je n - unupni broj slojeva modela.
!================================================================================

complex(8) function vektor_A(iso,isl,n,F)
	implicit none

!	Input
	integer iso,isl,n
	complex(8),dimension(:) :: F
!	Local variables
!	complex(8) prod
	complex(8),parameter :: one = dcmplx(1.d0,0.d0)
	complex(8) A
	integer i
!	Opis varijabli:
!	iso - sloj u kojem se nalazi tockasti izvor struje (s)
!	isl - sloj u kojem se trazi raspodjela potencijala (sloj)

	A = one

	if (isl<iso) then
		do i = isl,iso-1
			A = A * (one + F(i))
		end do

	else if(isl==iso) then
		A = one

	else if (isl>iso) then
		do i = iso,isl-1
			A = A * (one - F(i))
		end do
	end if
	
	vektor_A = A

	return
end function