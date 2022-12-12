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

! Proracun analitickog dijela medjusobne impedancije izmedju
! medjusobno paralelnih segmenata
! MEDJUSOBNA IMPEDANCIJA PARALELNIH SEGMENATA - ANALITICKI DIO

subroutine otpps(duls2,xp1,xk1,yp1,yk1,zp1,zk1,xp2,xk2,yp2,yk2,zp2,zk2,pol,andi)
	use funkcije
	implicit none
	
	real(8) duls2
	real(8) xp1,yp1,zp1,xk1,yk1,zk1
	real(8) xp2,yp2,zp2,xk2,yk2,zk2
	real(8) pol
	real(8) andi
	
	real(8) xs2,ys2,zs2
	real(8) xp1s,yp1s,zp1s,xk1s,yk1s,zk1s,xk2s,yk2s,zk2s,d2
	real(8) a1p,a2p,xp12,yp12,zp12,xk12,yk12,zk12,ak1,ak2
	real(8) v1,v2,ve
	integer indt
	real(8) a1a,a2a,a1,a2
	real(8) u1,u2,u3,u4
	
	! Opis varijabli:
	!  duls2 - duljina segmenta
	!  xp1,yp1,zp1 - x, y i z koordinate pocetne tocke prvog segmenta
	!  xk1,yk1,zk1 - x, y i z koordinate krajnje tocke prvog segmenta
	!  xp2,yp2,zp2 - x, y i z koordinate pocetne tocke drugog segmenta
	!  xk2,yk2,zk2 - x, y i z koordinate krajnje tocke drugog segmenta
	!  pol - polumjer (radijus) segmenata
	!  andi (output) - vrijednost analitickog dijela integrala, pri
	!                  proracunu medjusobne impedancije paralelnih segmenata
	! ----------------------------------------------------------------
	
	xs2 = (xp2+xk2)/2.d0
	ys2 = (yp2+yk2)/2.d0
	zs2 = (zp2+zk2)/2.d0
	
	xp1s = xp1-xs2
	yp1s = yp1-ys2
	zp1s = zp1-zs2
	xk1s = xk1-xs2
	yk1s = yk1-ys2
	zk1s = zk1-zs2
	xk2s = xk2-xs2
	yk2s = yk2-ys2
	zk2s = zk2-zs2
	d2 = duls2/2.d0
	a1p = (xp1s*xk2s + yp1s*yk2s + zp1s*zk2s)/d2
	a2p = (xk1s*xk2s + yk1s*yk2s + zk1s*zk2s)/d2
	xp12 = xp1s**2
	yp12 = yp1s**2
	zp12 = zp1s**2
	xk12 = xk1s**2
	yk12 = yk1s**2
	zk12 = zk1s**2
	ak1 = a1p**2
	ak2 = a2p**2
	v1 = dsqrt(dabs(xp12 + yp12 + zp12 - ak1))
	v2 = dsqrt(dabs(xk12 + yk12 + zk12 - ak2))
	ve = (v1+v2)/2.d0
	if (ve < pol) ve = pol
	indt = 0
	a1a = dabs(a1p)
	a2a = dabs(a2p)
	if (((a1p<0.d0).and.(a2p>0.d0)).or.((a2p<0.d0).and.(a1p>0.d0))) then
		a1 = a1a
		a2 = a2a
	else
		indt = 1
		if (a1a>a2a) then
			a1 = a2a
			a2 = a1a
		else
			a1 = a1a
			a2 = a2a
		end if
	end if
	u1 = a2+d2
	u2 = a1+d2
	u3 = a1-d2
	u4 = a2-d2
	if (indt==1) then
		andi = DOPS(u1,ve) + DOPS(u3,ve) - DOPS(u2,ve) - DOPS(u4,ve)
	else if(indt==0) then
		andi = DOPS(u1,ve) + DOPS(u2,ve) - DOPS(u3,ve) - dops(u4,ve) - 4.d0*d2*dlog(ve)
	end if
	
	return
end subroutine
