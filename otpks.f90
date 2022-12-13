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

! Proracun analitickog dijela integrala za proracun medjusobne
! impedancije medjusobno kosih segmenata
! MEDJUSOBNA IMPEDANCIJA KOSIH SEGMENATA - ANALITICKI DIO

subroutine otpks(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,rcosf,pol,andi)
    use funkcije
    implicit none
    
    real(8) x1,x2,x3,x4
    real(8) y1,y2,y3,y4
    real(8) z1,z2,z3,z4
    real(8) rcosf,pol
    real(8) andi
    
    real(8) xl1,xl2,zl1,zl2,delta
    
    
    ! Opis varijabli:
    !  x1,y1,z1 - x, y i z koordinate pocetne tocke prvog segmenta
    !  x2,y2,z2 - x, y i z koordinate krajnje tocke prvog segmenta
    !  x3,y3,z3 - x, y i z koordinate pocetne tocke drugog segmenta
    !  x4,y4,z4 - x, y i z koordinate krajnje tocke drugog segmenta
    !  rcosf - kosinus kuta izmedju segmenata
    !  pol - polumjer (radijus) segmenata
    !  andi (output) - vrijednost analitickog dijela integrala, pri
    !                  proracunu medjusobne impedancije paralelnih segmenata
    ! ----------------------------------------------------------------
    
    ! Poziv vanjske subroutine koja racuna lokalne koordinate kosih pravaca
    call LKKP(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xl1,xl2,zl1,zl2,delta)
    if (delta < pol) delta = pol
    
    andi = DOKS(xl1,zl1,rcosf,delta) + DOKS(xl2,zl2,rcosf,delta) - &
    DOKS(xl1,zl2,rcosf,delta) - DOKS(xl2,zl1,rcosf,delta)
    
    return
end subroutine
