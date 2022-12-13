! ...
! Odredjivanje By komponente magnetskog polja, preko vektorskog magnetskog potencijala.
! Vrijedi, naime, sljedeci izraz:

!             @Ax     @Az
!       Bx = ----- - ----- 
!              @z      @x

subroutine Hpolje_By(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,By)
    use funkcije
    implicit none
    
    integer N,Ns
    integer isl
    integer,dimension(:) :: iso
    real(8) :: xt,yt,zt,f
    real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
    complex(8),dimension(:) :: Iu
    real(8),dimension(:) :: HD,L
    complex(8),dimension(:) :: gama
    !Output
    complex(8) By
    ! Local variables
    integer iseg,si
    real(8) x1,y1,z1,x2,y2,z2,pol
    real(8) xsi,ysi,zsi,xsj,ysj,zsj,andi
    real(8) u,v,DuDz,DvDz,DFDu,DFDv,DFDz,DuDx,DvDx,DFDx,temp1,temp2
    real(8),parameter :: pi = 3.141592653589793238d0
    real(8) mio,konst,Rs
    complex(8) gamae,Fprig,ggama
    complex(8) DAxDz,DAzDx

    !Opis ulaznih varijabli:
    !  N - ukupni broj slojeva modela (veseslojno tlo + zrak)
    !  Ns - ukupni broj segmenata uzemljivackih sustava
    !  xt - x koordinata tocke promatranja (tocke u kojoj se racuna raspodjela potenc.)
    !  yt - y koordinata tocke promatranja
    !  zt - z koordinata tocke promatranja
    !  isl - redni broj sloja u kojem se nalazi tocka promatranja
    !  xp (realni vektor) - x koordinate pocetnih tocaka svih segmenata
    !  yp (realni vektor) - y koordinate pocetnih tocaka svih segmenata
    !  zp (realni vektor) - z koordinate pocetnih tocaka svih segmenata
    !  xk (realni vektor) - x koordinate krajnjih tocaka svih segmenata
    !  yk (realni vektor) - y koordinate krajnjih tocaka svih segmenata
    !  zk (realni vektor) - z koordinate krajnjih tocaka svih segmenata
    !  r0 (realni vektor) - radijus svakog od segmenata
    !  iso (realni vektor) - redni broj sloja u kojem se nalazi svaki od segmenata
    !  HD (realni vektor) - z koord. donjih granicnih ploha svih slojeva
    !  L (realni vektor) - duljine svih segmenata sustava uzemljivaca
    !  Iu (complex vektor) - vrijednosti poprecnih struja svih segmenata
    !  gama (complex vektor) - kompleksne valne konstante svih slojeva modela

    ! ---------------------------------------------------------------------------------
    !            POCETAK PRORACUNA KOMPONENTE Bx ZA MAGNETSKO POLJE
    ! ---------------------------------------------------------------------------------
    mio = 4.d0*pi*1.e-7
    konst = mio/(4.d0*pi)
    By = dcmplx(0.d0,0.d0)
    DO iseg = 1,Ns
        si = iso(iseg)
        x1 = xp(iseg); x2 = xk(iseg)
        y1 = yp(iseg); y2 = yk(iseg)
        z1 = zp(iseg); z2 = zk(iseg)
        pol = r0(iseg)

        xsi = (x2+x1)/2.d0; xsj = xt
        ysi = (y2+y1)/2.d0; ysj = yt
        zsi = (z2+z1)/2.d0; zsj = zt
        call gamaekv(N,gama,HD,xsi,ysi,zsi,xsj,ysj,zsj,gamae)
        Fprig = zexp(-gamae)

        Rs = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2)
        ggama = gamae/Rs
        
        call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi) !F

        ! Derivacija - DFDz
        u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
        v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
        DuDz = (z2-z1)/L(iseg)
        DvDz = (1.d0/v) * (zt-zsi-u*((z2-z1)/L(iseg)))
        temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
        temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
        DFDu = 1.d0/temp1 - 1.d0/temp2
        DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
        DFDz = DFDu*DuDz + DFDv*DvDz
        
        ! Derivacija Ax po z
        DAxDz = konst * (((x2-x1)/L(iseg))*Iu(iseg) * (-ggama*Fprig*((zt-zsi)/Rs)*andi + Fprig*DFDz))

        ! Derivacija - DFDx
        u = (2.d0/L(iseg)) * ((xt-xsi)*(x2-xsi) + (yt-ysi)*(y2-ysi) + (zt-zsi)*(z2-zsi))
        v = dsqrt((xt-xsi)**2+(yt-ysi)**2+(zt-zsi)**2-u**2)
        DuDx = (x2-x1)/L(iseg)
        DvDx = (1.d0/v) * (xt-xsi-u*((x2-x1)/L(iseg)))
        temp1 = dsqrt(v**2+(u+L(iseg)/2.d0)**2)
        temp2 = dsqrt(v**2+(u-L(iseg)/2.d0)**2)
        DFDu = 1.d0/temp1 - 1.d0/temp2
        DFDv = -((u+L(iseg)/2.d0)/v)*(1.d0/temp1) + ((u-L(iseg)/2.d0)/v)*(1.d0/temp2)
        DFDx = DFDu*DuDx + DFDv*DvDx
        
        ! Derivacija Az po x
        DAzDx = konst * (((z2-z1)/L(iseg))*Iu(iseg) * (-ggama*Fprig*((xt-xsi)/Rs)*andi + Fprig*DFDx))

        By = By + (DAxDz - DAzDx)
    END DO

    return
end subroutine