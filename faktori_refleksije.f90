! -------------------------------------------------------------------------------
!    Copyright (C) 2006. GPL - General Public Licence
!    Author: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)

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
!           Funkcija koja racuna faktore refleksije Fi za sve slojeve
!           n-slojnog modela (zrak + viseslojno tlo).
!================================================================================

!Faktori refleksije definiraju su na sljedeci nacin:

!       F(0) = 0

!               kapa(i) - kapa(i+1)
!       F(i) = ---------------------  ;   i = 1,2,...,n-1
!               kapa(i) + kapa(i+1)

!       F(n) = 0

!pri cemu je:
!   i - trenutna vrijednost sloja za koji se odreduje koeficijent refleksije
!   n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!   kapa(i) - kompleksna specificna elektricana vodljivost i-tog sloja
!             modela tla [S/m]. Definira se na sljedeci nacin:

!                       kapa(i) = sigma(i) + j*omega*eps(i)

!             pri cemu su:
!                   sigma(i) - specificna elektricna vodljivost i-tog sloja
!                              modela tla [S/m]. Vrijedi pritom sljedeca 
!                              jednakost:

!                                                        1
!                                           sigma(i) = ------
!                                                       ro(i)

!                              gdje ro(i) - predstavlja specificnu elektricnu
!                              otpornost i-tog sloja modela tla [Ohm/m].
!                   omega - kruzna frekvencija. Jednaka je sljedecem izrazu:

!                                           omega = 2*pi*f

!                           pri cemu f - predstavlja frekvenciju, [Hz].
!                   eps(i) - relativna permitivnost i-tog sloja tla. Racuna se kao:

!                                           eps(i) = eps0 * epsr(i)

!                            gdje su:
!                               eps0 - relativna permitivnost slobodnog prostora.
!                                      Ona iznosi:

!                                           eps0 = 8.854e-12 [As/Vm]

!                               epsr(i) - relativna permitivnost i-tog sloja
!                                         viseslojnog modela tla, [As/Vm].
!                   j - imaginarna jedinica.

!             Napomena:
!             ---------
!             Vektor kapa je kompleksan i sadrzi n vrijednosti. kapa(1) se odnosi
!             na zrak. Varijable n; f; ro(i),i=1,n; epsr(i),i=1,n zadaju se kao 
!             ulazni podaci.

!Napomena:
!---------
!Vektor faktora refleksije je kompleksan, duljine n. Vrijednosti vektora F(0) i
!F(-1) ne postoje!!! 
!===================================================================================

subroutine vektor_F(n,kapa,F)
    implicit none
!   Input
    integer n
    complex(8),dimension(:) :: kapa
!   Output
    complex(8),dimension(:) :: F
!   Local variable
    integer i

!   Tijelo funkcije
    do i = 1,n-1
        F(i) = (kapa(i) - kapa(i+1)) / (kapa(i) + kapa(i+1))
    end do
    F(n) = dcmplx(0.0,0.0)

    return
end subroutine