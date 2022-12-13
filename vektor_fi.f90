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
!   Subroutina koja odreduje vrijednost fi za sve uzorkovane tocke.
!   Ova vrijednost je definirana u Nada Magistarski rad, str. 21. 
!   Prepravljen je u mojim materijalima na str. 18, jer je alfa 
!   zamjenjen s lamda (lamda = 1/u)
!================================================================================

!   Vektor fi definiran je na sljedeci nacin:

!       za  i < s:  

!                       -lamda*(d - H )                     F  * F      -2*lamda*(H  - d)
!                                    i                2      i    s                s
!           f  =  A  * e               * [ F    * ni   +  ---------- * e                  ]
!            i     i                        i-1     i       1 - F
!                                                                i

!       za i >= s:


!                      -lamda*(H  - d)                                         -2*lamda*(d - H   )
!                               i                             2                               s-1
!           f  = A  * e                * [ (1 + F ) * F   * ni    - F * F   * e                    ]
!            i    i                              i     i+1    i+1    i   s-1
!
!================================================================================

subroutine vektor_fi(n,k,d,lamda,iso,isl,HD,Ft,Fr,H,fi)
    implicit none
!   Input
    integer n,k
    real(8) d
    real(8),dimension(:) :: lamda
    integer iso,isl
    complex(8) Ft
    complex(8),dimension(:) :: Fr
    real(8),dimension(:) :: HD,H
!   Output
    complex(8) fi
!   Local variables
    complex(8),parameter :: one = (1.d0,0.d0)


!   Objasnjenje ulaznih varijabli:
!       n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!       u - vektor uzorkovanih tocaka. Odredjen je na sljedeci nacin:

!                   u  =  teta * w   ili  u  = hi * w
!                    i            i        i         i

!           pri cemu su vrijednosti w prethodno izracunate (vidjeti:
!           file: parametri_eta_w.f90).
!       d - dubina ukopavalja tockastog izvora struje, [m]
!       iso - sloj u kojem se nalazi tockasti izvor struje
!       isl - sloj za koji se racuna vrijednost fi
!       Ft -  faktor transmisije. Definiran je na sljedeci nacin:

!           A = (1 + F(s-1))*(1 + F(s-2))*...*(1 + F(i)) ;  za: i <= s - 1

!           A = 1 :  za: i = s

!           A = (1 - F(s))*(1 - F(s+1))*...*(1-F(i-1)) ;  za: i >= s + 1

!            i odredjen u file-u: faktori_transmisije.f90.
!       Fr - vektor faktora refleksije. Definiran je na sljedeci nacin:

!               F(0) ne postoji

!                       kapa(i) - kapa(i+1)
!               F(i) = ---------------------  ;   i = 1,2,3,...,n-1
!                       kapa(i) + kapa(i+1)

!               F(n) = 0

!            i odredjeni su u file-u: faktori_refleksije.f90.
!       H  - vektor debljina slojeva modela (h_sloj), [m]
!       HD - vektor z koordinata donjih granicnih ploha svih slojeva
!            viseslojnog modela. Ovaj vektor odredjen je u file-u:
!            debljina_slojeva.f90 
!       lamda - uzorkovana vrijednost koja se definira na sljedeci nacin:

!                                    1
!                       lamda(i) = ------
!                                   u(i)

!               dok je u prethodno definirana uzorkovana vrijednost.
!       fi (output) - vrijednost fi
!--------------------------------------------------------------------------------

    if ((iso==1).and.(isl==1)) then
        fi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * (one+Fr(isl))*Fr(isl+1)* &
        (dexp(-lamda(k)*h(isl+1)))**2
    else if((iso==1).and.(isl/=1)) then
        if (isl==(n-1)) then
            fi = dcmplx(0.d0,0.d0)
        else
            fi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * (one+Fr(isl))*Fr(isl+1)* &
            (dexp(-lamda(k)*h(isl+1)))**2 
        end if
    else if ((iso==n).and.(isl==1)) then
        fi = dcmplx(0.d0,0.d0)
    else if((iso==n).and.(isl/=1)) then
        fi = Ft * dexp(-lamda(k)*(d-HD(isl))) * Fr(isl-1) * (dexp(-lamda(k)*h(isl)))**2
    else if((iso/=1).and.(iso/=n).and.(isl==1)) then
        fi = Ft * dexp(-lamda(k)*(d-HD(isl))) * ((Fr(isl)*Fr(iso))/(one-Fr(isl))) * &
        dexp(-2.d0*lamda(k)*(HD(iso)-d))
    else if((iso/=1).and.(iso/=n).and.(isl/=1)) then
        if (isl<iso) then
            fi = Ft * dexp(-lamda(k)*(d-HD(isl))) * ( Fr(isl-1)*(dexp(-lamda(k)*h(isl)))**2 + &
            ((Fr(isl)*Fr(iso))/(one-Fr(isl))) * dexp(-2.d0*lamda(k)*(HD(iso)-d)) )
        else
            if (isl==(n-1)) then
                fi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * (-Fr(isl)*Fr(iso-1)*dexp(-2.d0*lamda(k)*(d-HD(iso-1))))
            else
                fi = Ft * dexp(-lamda(k)*(HD(isl)-d)) * ( (one+Fr(isl))*Fr(isl+1)*(dexp(-lamda(k)*h(isl+1)))**2 - &
                Fr(isl)*Fr(iso-1)*dexp(-2.d0*lamda(k)*(d-HD(iso-1))) )
            end if
        end if
    end if

    return
end subroutine