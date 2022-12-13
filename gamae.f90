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

! Proracun ekvivalentne valne konstante izmedju dvaju segmenata
! uzemljivaca (uzimaju se njihove sredisnje tocke kao referentne)
! Ekvivalentna valna konstanta za segment uzemljivaca ili izmedju
! segmenata uzemljivaca ovisi o njihovim medjusobnim polozajima
! (broju slojeva koji ih dijeli) i racuna se za svaki par
! segmenata posebno, ovom subroutineom

! Potrebno je prije poziva ove funkcije izracunati koordinate
! sredine i-tog i j-tog segmenta, na sljedeci nacin:
! do i = 1,br_seg
!   xs(i) = (xk(i)+xp(i))/2.d0
!   ys(i) = (yk(i)+yp(i))/2.d0
!   zs(i) = (zk(i)+zp(i))/2.d0
! end do
! do i = 1,br_seg
!   do j = i,br_seg
!       call gamaekv(n,gama,H,xs(i),ys(i),zs(i),xs(j),ys(j),zs(j),gamae)
!       Z(i,j) = ......
!       Z(j,i) = Z(i,j)
!   end do
! end do

! Napomena: Ukoliko je rijec o vlastitom segmentu, ili se segmenti nalaze u 
!           istom sloju, tada je gamae upravo jednaka valnoj konstanti sre-
!           dstva u kojem se nalazi segment. Ova subroutina u tom slucaju 
!           vraca vrijednost gama*d, pri cemu su: gama valna konstanta tog
!           sredstva, a d udaljenost medju sredistima segmenata (ovo moze 
!           biti i slucaj segmenta i tocke promatranja).
!           Ukoliko se trazi gamae izmedju dvaju segmenata koji se nalaze u
!           razlicitim slojevima, tada ova rutina vraca vrijednost gamae*d0,
!           pri cemu je d0 ekvivalentna udaljenost izmedju segmenata. 
! -----------------------------------------------------------------------------

subroutine gamaekv(n,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
    implicit none

    integer n
    complex(8),dimension(:) :: gama
    real(8),dimension(:) :: HD
    real(8) xsi,ysi,zsi,xsj,ysj,zsj
    complex(8) gamae
    !Lokalne varijable
    integer is,js,i
    real(8) d,d0
    real(8) x1,y1,z1,x2,y2,z2,x0,y0,z0

    !Opis varijabli:
    !  n - ukupni broj slojeva modela (viseslojno tlo + zrak)
    !  gama - valne konstante svih n slojeva (kompleksni vektor)
    !  HD - koordinate donjih granicnih ploha slojeva modela
    !  xsi,ysi,zsi - koordinate srednje tocke i-tog segmenta
    !  xsj,ysj,zsj - koordinate srednje tocke j-tog segmenta
    !  gamae - ekvivalentna valna konstanta izmedju segmenata
    
    ! -----------------------------------------------------------
    !                 POCETAK PRORACUNA
    ! -----------------------------------------------------------
    is = 1
    js = 1
    !Odredjivanje sloja u kojem se nalaze i-ti i j-ti segment
    do i = 1,n-1
        if (zsi>HD(i)) is = is + 1
        if (zsj>HD(i)) js = js + 1
    end do
    
    if (is==js) then
        !Segmenti se nalaze u istom sloju
        d = dsqrt((xsi-xsj)**2+(ysi-ysj)**2+(zsi-zsj)**2)
        gamae = gama(is)*d

    else
        !Segmenti su u razlicitim slojevima
        if (is>js) then
            gamae = dcmplx(0.d0,0.d0)
            d = dsqrt((xsi-xsj)**2+(ysi-ysj)**2+(zsi-zsj)**2)
            
            x1 = xsj
            y1 = ysj
            z1 = zsj

            x2 = xsi
            y2 = ysi
            z2 = zsi
            
            do i = js,is-1
                z0 = HD(i)
                x0 = x1 + ((x2-x1)/(z2-z1))*(z0-z1)
                y0 = y1 + ((y2-y1)/(z2-z1))*(z0-z1)

                d0 = dsqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)

                gamae = gamae + gama(i)*d0

                x1 = x0
                y1 = y0
                z1 = z0
            end do

            d0 = dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

            gamae = gamae + gama(is)*d0
        
        else
            gamae = dcmplx(0.d0,0.d0)
            d = dsqrt((xsj-xsi)**2+(ysj-ysi)**2+(zsj-zsi)**2)

            x1 = xsi
            y1 = ysi
            z1 = zsi

            x2 = xsj
            y2 = ysj
            z2 = zsj

            do i = is,js-1
                z0 = HD(i)
                x0 = x1 + ((x2-x1)/(z2-z1))*(z0-z1)
                y0 = y1 + ((y2-y1)/(z2-z1))*(z0-z1)

                d0 = dsqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)

                gamae = gamae + gama(i)*d0

                x1 = x0
                y1 = y0
                z1 = z0
            end do

            d0 = dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

            gamae = gamae + gama(js)*d0

        end if
    end if

    return
end subroutine