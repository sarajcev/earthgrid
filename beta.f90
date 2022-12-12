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

! Odredjivanje koeficijenata BETA za odslikavanje segmenta u visesloju. Ova
! subroutina sadrzava kompletan proracun.

subroutine betauz(hi,w,eta,h,n,d,Fr,HD,s,kapa,sloj,beta)
	use funkcije
	implicit none
	
	real(8) hi
	real(8),dimension(:) :: eta,w
	real(8),dimension(:) :: h,HD
	integer n
	real(8) d
	complex(8),dimension(:) :: Fr
	integer s,sloj
	complex(8),dimension(:) :: kapa
	complex(8),dimension(:) :: beta
	! Lokalne varijable
	real(8),dimension(25) :: u
	integer m_E,n_E
	real(8),dimension(15,25) :: E
	real(8),dimension(25) :: lamda
	complex(8),dimension(25) :: tetan
	complex(8),dimension(25) :: teta_k,hi_k
	integer i
	! Opis varijabli:
	!  teta - koeficijent odslikavanja
	!  eta - nepromjenjene vrijednosti uzorkovanih tocaka
	!  w - nepromjenjene vrijednosti uzorkovanih tocaka
	!  gamae - ekvivalentna valna kosntanta
	!  h - vektor debljina svih slojeva
	!  HD - vektor koordinata donjih granicnih ploha svih slojeva
	!  n - ukupni broj slojeva modela (viseslojno tlo + zrak)
	!  d - dubina na kojoj se nalazi tockasti izvor (kosi segment
	!      se, naime, zamjenjuje sa 5 tockastih izvora
	!  Fr - faktori refleksija svih slojeva
	!  Ft - ukupni faktor transmisije
	!  s - sloj u kojem se nalazi tockasti izvor (kosi segment se
	!      zamjenjuje sa 5 tockastih izvora) - iso
	!  sloj - sloj u kojem se nalazi tocka promatranja (isl)
	!  kapa - vrijednosti specif. elektr. vodljivosti svih slojeva
	!  alfa (output) - izracunati koefcijenti ALFA (ukupno 15)
	! ------------------------------------------------------------
	
	! Uzorkovane tocke
	u(:) = hi * w(:)
	
	m_E = 15
	n_E = 25
	call INVMAT(E,m_E,n_E)
	
	lamda(1:24) = 1.d0/u(1:24)
	lamda(25) = 0.d0
	
	! Odredjivanje vrijednosti teta_n u 25 uzorkovanih tocaka
	call VEKTOR_TETA_N(n,d,h,HD,lamda,kapa,Fr,sloj,s,tetan)
	
	! Odredjivanje vektora TETA u 25 uzorkovanih tocaka
	! (pomocu rekurzivnih relacija)
	call VEKTOR_TETA_HI(n,d,tetan,kapa,lamda,sloj,s,Fr,h,HD,teta_k,hi_k)
	
	! Odredjivanje nepoznatih koeficijenata BETA mnozenjem s 
	! pseudoinverznom matricom [E]
	beta = matmul(E,hi_k)
	
	return
end subroutine
