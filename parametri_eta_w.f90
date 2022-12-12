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
!			Funkcije koje racunaju nepromjenjive parametre eta i w.
!			Ovi parametri su definirani u: Nada Magistarski rad, str. 33
!================================================================================

!Parametri eta definirani su na sljedeci nacin:

!				eta  = 1.0
!				   1

!				                 14 ______
!								 --|               k-8.5
!				eta  =  eta    *   | 2000  * (1.05)        ;  k = 2,3,...,15
!				   k	   k-1	   |			 

!		Napomena:
!		---------
!		Dakle, postoji ukupno 15 parametara eta. Eta je stoga realni vektor
!		duljine 15.
!================================================================================

subroutine vektor_eta(n,eta)
	implicit none

!	Input
	integer n
!	Output
	real(8),dimension(15) :: eta
!	Local variable
	integer k


	!Tijelo funkcije
	eta(1) = 1.d0
	do k = 2,15
		eta(k) = eta(k-1) * ((2000.d0)**(1.d0/14.d0)) * ((1.05d0)**(k-8.5d0))
	end do

	return
end subroutine


!================================================================================
!Parametri w definirani su na sljedeci nacin:

!				  w  = 0.2
!				   1

!				                 23 _______
!								 --| 2500    
!				  w  =    w    *   | -----      ;  k = 2,3,...,24
!				   k	   k-1	   |  0.2			 


!				  w   = oo
!				   25 

!		Napomena:
!		---------
!		Dakle, postoji ukupno 25 nepromjenjivih parametara  w. w je stoga realni
!		vektor duljine 25 (posljednji clan koji je beskonacan - stavljen je nula
!		jer se kasnije korigira).
!================================================================================

subroutine vektor_w(n,w)
	implicit none

!	Input
	integer n
!	Output
	real(8),dimension(25) :: w
!	Local variable
	integer k


	!Tijelo funkcije
	w(1) = 0.2d0
	do k = 2,24
		w(k) = w(k-1) * ((2500.d0/0.2d0)**(1.d0/23.d0))
	end do
	w(25) = 0.d0

	return
end subroutine