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

! Funkcija za proracun osnovnog dijela izraza za impedanciju
! (otpor) izmedju medjusobno kosih segmenata, pri analitickom
! rjesavanju dijela integrala za medjusobne impedancije kosih
! segmenata. Ova funkcija poziva se iz subroutine: otpks

real(8) function doks(x,z,rcos,delta)
    implicit none
    
    real(8) x,z,rcos,delta
    ! LOkalne varijable
    real(8) xzc,zx,xz,xpz,d2
    real(8) rkor,sinf,pl1,pl2,xzk,tgf2
    
    
    xzc = 2.d0*x*z*rcos
    zx = z-x*rcos
    xz = x-z*rcos
    xpz = x+z
    d2 = delta**2
    rkor = dsqrt(x**2 + z**2 + d2 - xzc)
    sinf = dsqrt(dabs(1.d0 - rcos**2))
    if (zx >= 0.d0) then
        pl1 = rkor+zx
    else
        pl1 = (d2+x**2*sinf**2)/(rkor-zx)
    end if
    if (xz >= 0.d0) then
        pl2 = rkor+xz
    else
        pl2 = (d2+z**2*sinf**2)/(rkor-xz)
    end if
    if (xpz >= 0.d0) then
        xzk = rkor+xpz
    else
        xzk = (d2 - 2.d0*x*z*(1.d0+rcos))/(rkor-xpz)
    end if
    tgf2 = dsqrt(dabs((1.d0-rcos)/(1.d0+rcos)))
    
    doks = x*dlog(pl1) + z*dlog(pl2) + (2.d0*delta*datan((xzk*tgf2)/delta))/sinf
    
    return
end function
