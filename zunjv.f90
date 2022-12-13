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

! PRORACUN JEDINICNE UNUTARNJE IMPEDANCIJE VODICA 'Zunjv'
! Formula po kojoj se racuna nalazi se u Nada, Magistarski rad.

subroutine proracun_zunjv(Ns,f,mirv,sigmav,r0,zunjv)
    implicit none
    
    integer Ns
    real(8) f
    real(8),dimension(:) :: mirv,sigmav,r0
    complex(8),dimension(:) :: zunjv
    ! Lokalne varijable
    real(8) omega
    real(8),parameter :: pi = 3.141592653589793238d0
    real(8) mi0,mir,t
    complex(8) z,p,k,pp
    real(8) xr,xi
    integer NB
    real(8),dimension(2) :: BR,BI
    complex(8) Jo,J1
    integer i
    
    external DBESCJ

    ! Opis ulaznih varijabli:
    ! -----------------------
    !   Ns - ukupni broj segmenata sustava uzemljenja
    !   f - frekvencija za koju se racuna impedancija
    !   mirv - (realni vektor) relativna permeabilnost svakog od segmenata
    !   sigmav - (realni vektor) specificna vodljivost svakog od segmenata
    !   r0 - (realni vektor) radijus svakog od segmenata
    !   zunjv (output) - (kompleksni vektor) jedinicna unutarnja impedancija
    !                    svakog od segmenata sustav uzemljenja

    ! Napomena: Za proracun jedinicne unutarnje impedancije segmenata uzemljivaca
    !           potrebno je izracunati vrijednosti BESSEL-ovih funkcija I vrste
    !           nultog i prvog reda, kompleksnog argumenta (Bessel function of
    !           the First kind of integer order and complex argument). Ovaj 
    !           proracun omogucava vanjska rutina: DBESCJ.
    ! ----------------------------------------------------------------------
    

    omega = 2.d0*pi*f
    mi0 = 4.d0*pi*1.e-7
    t = (3.d0/4.d0)*pi
    z = dcmplx(0.d0,t)
    p = zexp(z)
    NB = 2
    
    do i = 1,Ns
        mir = mi0*mirv(i)
        k = dsqrt(omega*mir*sigmav(i))*p
        pp = k*r0(i)
        xr = dreal(pp)
        xi = dimag(pp)
                
        ! Poziv subroutine za proracun BESSEL-ovih funkcija
        call DBESCJ(xr,xi,NB,BR,BI) 
            
        Jo = dcmplx(BR(1),BI(1))
        J1 = dcmplx(BR(2),BI(2))

        zunjv(i) = (k/(2.d0*r0(i)*pi*sigmav(i)))*(Jo/J1)
    end do
    
    return
end subroutine