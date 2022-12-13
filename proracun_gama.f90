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
!       Subroutina koja racuna valne konstante (gama) svih slojeva modela
!================================================================================

! Valna konstanta i-tog sloja definirana je na sljedeci nacin:

!               -----------------
!     gama  = -|j*omega*mi *kapa
!         i    |          0     i

! pri cemu su:
!    j - imaginarna jedinica
!    omega - kruzna frekvencija, koja glasi:

!                 omega = 2*pi*f

!            gdje je f - frekvencija.
!    mi - permeabilnost slobodnog prostora. Ona iznosi:
!      0

!              mi = 4*pi*1.e-7
!                0

!    kapa - kompleksna specificna elektricna vodljivost i-tog sloja
!        i
! ------------------------------------------------------------------------------

subroutine proracun_gama(f,n,kapa,gama)
    use funkcije
    implicit none

    real(8) f
    integer n
    complex(8),dimension(:) :: kapa
    complex(8),dimension(:) :: gama
    !Lokalne varijable
    real(8) omega
    real(8),parameter :: pi = 3.141592653589793238d0
    real(8) mi0
    real(8) temp
    complex(8) zemp
    integer i


    mi0 = 4.d0*pi*1.e-7
    omega = 2.d0*pi*f
    temp = omega*mi0
    do i = 1,n
        zemp = dcmplx(0.d0,temp)*kapa(i)
        gama(i) = zsqrt(zemp)
    end do

    return
end subroutine