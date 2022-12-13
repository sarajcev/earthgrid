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

! Proracun kosinusa kuta izmedju pravaca na kojima leze
! segmenti - da se ustanovi njihov medjusobni polozaj

subroutine koskutsegmenata(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfi)
    implicit none
    
    real(8) x1,y1,z1,x2,y2,z2
    real(8) x3,y3,z3,x4,y4,z4
    real(8) cosfi
    
    real(8) brojnik,naz1,naz2,nazivnik
    
    ! Opis varijabli:
    !  x1,y1,z1 - x, y i z koordinate pocetne tocke prvog segmenta
    !  x2,y2,z2 - x, y i z koordinate krajnje tocke prvog segmenta
    !  x3,y3,z3 - x, y i z koordinate pocetne tocke drugog segmenta
    !  x4,y4,z4 - x, y i z koordinate krajnje tocke drugog segmenta
    !  cosfi (output) - kosinus kuta izmedju segmenata
    
    
    brojnik = (x2-x1)*(x4-x3) + (y2-y1)*(y4-y3) + (z2-z1)*(z4-z3)
    naz1 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
    naz2 = (x4-x3)**2 + (y4-y3)**2 + (z4-z3)**2
    nazivnik = dsqrt(naz1*naz2)
    
    cosfi = dabs(brojnik/nazivnik)
    
    return
end subroutine


! Proracun kosinusa kuta izmedju pravaca na kojima leze
! segmenti - da se ustanovi njihov medjusobni polozaj -
! ZA MATRICU UZDUZNIH IMPEDANCIJA

subroutine koskutsegmenata1(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cosfiu)
    implicit none
    
    real(8) x1,y1,z1,x2,y2,z2
    real(8) x3,y3,z3,x4,y4,z4
    real(8) cosfiu
    
    real(8) brojnik,naz1,naz2,nazivnik
    
    ! Opis varijabli:
    !  x1,y1,z1 - x, y i z koordinate pocetne tocke prvog segmenta
    !  x2,y2,z2 - x, y i z koordinate krajnje tocke prvog segmenta
    !  x3,y3,z3 - x, y i z koordinate pocetne tocke drugog segmenta
    !  x4,y4,z4 - x, y i z koordinate krajnje tocke drugog segmenta
    !  cosfiu (output) - kosinus kuta izmedju segmenata
    
    
    brojnik = (x2-x1)*(x4-x3) + (y2-y1)*(y4-y3) + (z2-z1)*(z4-z3)
    naz1 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
    naz2 = (x4-x3)**2 + (y4-y3)**2 + (z4-z3)**2
    nazivnik = dsqrt(naz1*naz2)
    
    cosfiu = brojnik/nazivnik
    
    return
end subroutine
