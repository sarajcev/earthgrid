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
!	Subroutina koja odreduje vrijednosti vektora teta i hi. Ovaj vektor se definira
!	rekurzivnim relacijama. Vektor teta je jezgra funkcija koja se ovdje definira
!	u 25 uzorkovanih tocaka, za sve slojeve tla. Rekurzivni algoritam dan je u
!	Nada Magistarski rad, str. 25. (korigiran je u ovom slucaju)
!================================================================================

!	Dakle, ovdje se odreduju sljedece velicine:

!			teta (u )   , i = 1,2,...,25.
!				   i

!			hi (u )     , i = 1,2,...,25.
!                i

!	uvazavajuci cinjenicu da je:

!			teta (u )   , i = 1,2,...,25.
!				n  i

!	vec prethodno odredjen.

subroutine vektor_teta_hi(n,d,tetan,kapa,lamda,isl,iso,Fr,H,HD,teta,hi)
	use funkcije
	implicit none

!	Input
	integer n
	integer isl,iso
	real(8) d
	complex(8),dimension(:) :: tetan
	complex(8),dimension(:) :: kapa,Fr
	complex(8),dimension(:),allocatable :: qi,ai,bi
	real(8),dimension(:) :: H,HD
!	Output
	complex(8),dimension(:) :: teta,hi
!	Local variables
	real(8) ni,nip
	real(8),dimension(:) :: lamda
	complex(8) TI,HII,TIP,HIP
	complex(8) AF
	complex(8) fi,Gi
	integer k,i

!	Objasnjenje ulaznih varijabli:
!	------------------------------
!		n - ukupni broj slojeva viseslojnog modela (zrak + viseslojno tlo)
!		tetan - kompleksni vektor prethodno izracunatih vrijednosti teta(n) za
!				25 uzorkovanih vrijednosti. On je poznat.
!		qi,ai,bi - vektori koji su takodjer prethodno odredjeni (vidjeti file:
!				   vektori_q_a_b.f90).
!		u - vektor uzorkovanih tocaka
!		iso - sloj u kojem se nalazi tockasti izvor izmjenicne struje
!		isl - sloj za koji se racunaju parametri teta
!		H  - debljine svih slojeva modela (h_sloj), [m]
!		HD - koordinate donjih ploha svih slojeva, [m] - file: debljina_slojeva.f90


!	-------------------------------
!		Pocetak proracuna
!	-------------------------------
	allocate(qi(n-1))
	allocate(ai(n-1))
	allocate(bi(n-1))
	call vektori_q_a_b(n,kapa,qi,ai,bi)

	do k = 1,25
		TI = tetan(k)
		HII = dcmplx(0.0,0.0)
		if (isl/=n) then
			do i = n-1,isl,-1
				if (i==(n-1)) then
					nip = 0.d0
				else
					nip = dexp(-lamda(k)*H(i+1))
				end if
				TIP = TI
				HIP = HII
				AF = vektor_A(iso,i,n,Fr)
		
				call vektor_Gi(n,k,d,lamda,iso,i,HD,AF,Fr,H,Gi)
				HII = Gi - ai(i)*TIP + bi(i)*nip*HIP
		
				ni = dexp(-lamda(k)*H(i))
				call vektor_fi(n,k,d,lamda,iso,i,HD,AF,Fr,H,fi)
				TI = (fi - HII + qi(i)*TIP + qi(i)*nip*HIP) / ni
			end do
		end if
		teta(k) = TI
		hi(k) = HII
	end do

	deallocate(qi)
	deallocate(ai)
	deallocate(bi)

	return
end subroutine