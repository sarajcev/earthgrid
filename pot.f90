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

! Proracun potencijala izmedju tocke promatranja i segmenta,
! analiticki dio ukupnog potencijala
! POTENCIJAL TOCKE - ANALITICKI DIO

subroutine pot(x2,y2,z2,x,y,z,duls,xs1,ys1,zs1,pol,andi)
    implicit none
    
    real(8) x2,y2,z2
    real(8) x,y,z
    real(8) xs1,ys1,zs1
    real(8) duls,pol
    real(8) andi
    ! Lokalne varijable
    real(8) x1,y1,z1,d2
    real(8) ds,a,ve,ad1,ad2,br,rn

    ! Opis ulaznih varijabli:
    !    x2,y2,z2 - x, y i z koordinata krajnje tocke segmenta (K), respektivno
    !    x,y,z - x, y i z koordinata tocke promatranja (T), respektivno
    !    xs1,ys1,zs1 - x, y, i z koordinata srednje tocke segmenta (S), respektivno
    !    pol - radijus doticnog segmenta
    
    
    d2 = duls/2.d0
    x1 = dabs(x-xs1)
    y1 = dabs(y-ys1)
    z1 = dabs(z-zs1)
    if (x1 < 1.e-15) x1 = 0.d0
    if (y1 < 1.e-15) y1 = 0.d0
    if (z1 < 1.e-15) z1 = 0.d0

    ds = dsqrt(x1**2 + y1**2 + z1**2)
    a = dabs((x-xs1)*(x2-xs1) + (y-ys1)*(y2-ys1) + (z-zs1)*(z2-zs1))/d2
    ve = dsqrt(dabs(ds**2 - a**2))
    ad1 = a + d2
    ad2 = a - d2
    if ((ve < pol).and.(ad2 < pol)) ve = pol
    br = dlog(dsqrt(ve**2+ad1**2)+ad1)
    rn = dlog(dsqrt(ve**2+ad2**2)+dabs(ad2))
    if (ad2 < 0.d0) rn = 2.d0*dlog(ve)-rn
    
    andi = br - rn
    
    return
end subroutine
