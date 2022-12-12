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
!	Subroutina koja odreduje vrijednost vektora Gi. Ovaj vektor je definiran
!	u Nada Magistarski rad, str. 22. Prepravljen je u mojim materijalima
!	na str. 19, jer je alfa zamjenjen s lamda (lamda = 1/u)
!================================================================================

!	Vektor Gi definiran je na sljedeci nacin:

!		za  i < s:	

!                          2
!					      F  * F      -lamda*(2H  - H  - d)
!                          i    s               s    i
!			G  =  A  *  ---------- * e                  
!			 i	   i           2
!					      1 - F
!							   i

!		za i >= s:


!                      -lamda*(H  - d)                               -2*lamda*(d - H   )
!                               i                  2                                s-1
!			G  = A  * e                * [ F   * ni    - F * F   * e                    ]
!			 i    i                         i+1    i+1    i   s-1
!
!================================================================================

subroutine vektor_Gi(n,k,d,lamda,iso,isl,HD,Ft,Fr,H,Gi)
	implicit none
!	Input
	integer n,k
	real(8),dimension(:) :: lamda
	real(8) d
	integer isl,iso
	complex(8),dimension(:) :: Fr
	complex(8) Ft
	real(8),dimension(:) :: HD,H
!	Output
	complex(8) Gi
!	Local variables
	real(8),parameter :: one = (1.d0,0.d0)


!	Objasnjenje ulaznih varijabli:
!		n - bkupni roj slojeva viseslojnog modela (zrak + viseslojno tlo)
!		u - vektor uzorkovanih tocaka. Odredjen je na sljedeci nacin:

!							u  =  teta * w  ili  u  = hi * w
!							 i            i       i         i

!		    pri cemu su vrijednosti w prethodno izracunate (vidjeti:
!			file: parametri_eta_w.f90).
!		d - dubina ukopavalja tockastog izvora struje, [m]
!		iso - sloj u kojem se nalazi tockasti izvor struje
!		isl - sloj u kojem se trazi vrijednost Gi
!		Ft - faktor transmisije. Definiran je na sljedeci nacin:

!			A = (1 + F(s-1))*(1 + F(s-2))*...*(1 + F(i)) ;  za: i <= s - 1

!			A = 1 :  za: i = s

!			A = (1 - F(s))*(1 - F(s+1))*...*(1-F(i-1)) ;  za: i >= s + 1

!			 i odredjeni su u file-u: faktori_transmisije.f90.
!		Fr - vektor faktora refleksije. Definiran je na sljedeci nacin:

!				F(0) = ne postoji
	
!				        kapa(i) - kapa(i+1)
!				F(i) = ---------------------  ;   i = 1,2,3,...,n-1
!						kapa(i) + kapa(i+1)

!				F(n) = 0

!			 i odredjeni su u file-u: faktori_refleksije.f90.
!		H  - vektor debljina slojeva modela (h_sloj), [m].
!		HD - vektor z koordinata donjih granicnih ploha svih slojeva
!			 viseslojnog modela 
!		lamda - uzorkovana vrijednost koja se definira na sljedeci nacin:

!								     1
!						lamda(i) = ------
!								    u(i)

!				dok je u prethodno definirana uzorkovana vrijednost.
!		Gi - (output) vrijednost Gi
!--------------------------------------------------------------------------------

	if ((iso==1).and.(isl==1)) then
		Gi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * Fr(isl+1)*(dexp(-lamda(k)*h(isl+1)))**2
	else if((iso==1).and.(isl/=1)) then
		if(isl==(n-1)) then
			Gi = dcmplx(0.d0,0.d0)
		else	
			Gi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * Fr(isl+1)*(dexp(-lamda(k)*h(isl+1)))**2
		end if
	else if((iso==n).and.(isl==1)) then
		Gi = dcmplx(0.d0,0.d0)
	else if((iso==n).and.(isl/=1)) then
		Gi = dcmplx(0.d0,0.d0)
	else if((iso/=1).and.(iso/=n).and.(isl==1)) then
		Gi = Ft * ((Fr(isl)**2*Fr(iso))/(one-Fr(isl)**2)) * dexp(-lamda(k)*(2.d0*HD(iso)-HD(isl)-d))
	else if((iso/=1).and.(iso/=n).and.(isl/=1)) then
		if (isl<iso) then
			Gi = Ft * ((Fr(isl)**2*Fr(iso))/(one-Fr(isl)**2)) * &
			dexp(-lamda(k)*(2.d0*HD(iso)-HD(isl)-d))
		else
			if(isl==(n-1)) then
				Gi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * (-Fr(isl)*Fr(iso-1)*&
				dexp(-2.d0*lamda(k)*(d-HD(iso-1))))
			else
				Gi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * ( Fr(isl+1)*(dexp(-lamda(k)*h(isl+1)))**2 - &
				Fr(isl)*Fr(iso-1)*dexp(-2.d0*lamda(k)*(d-HD(iso-1))) )
			end if
		end if
	end if

	return
end subroutine