! ...
! Proracun Ax, Ay i Az komponente potrebne za proracun Ex, Ey i Ez komponente
! elektricnog polja. Naime, vrijede opcenito sljedeci izrazi:

!             @Fi
!     Ex = - ----- - j*omega*Ax
!              @x

!             @Fi
!     Ey = - ----- - j*omega*Ay
!              @y

!             @Fi
!     Ez = - ----- - j*omega*Az
!              @z

! Ova subroutina racuna Ax, Ay i Az za potrebe gornjih izraza, i to od ukupnog broja
! segmenata, dakle, ukupne Ax, Ay i Az (od svih segemenata). Komponente Ax, Ay i Az
! iste su, bilo da je rijec o dvosloju ili visesloju, kao i o tome da li je segment
! u horizontalnom ili kosom polozaju.

subroutine Epolje_Ax_Ay_Az(N,Ns,isl,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,f,Ax,Ay,Az)
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
	complex(8) Ax,Ay,Az
	! Local variables
	integer iseg,si
	real(8) x1,y1,z1,x2,y2,z2,pol
	real(8) xsi,ysi,zsi,xsj,ysj,zsj,andi
	real(8),parameter :: pi = 3.141592653589793238d0
	real(8) mio,konst,omega
	complex(8) jomega,gamae,Fprig,sumax,sumay,sumaz

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
	!  f - frekvencija za koju se racuna vrijednost Ax

	! ---------------------------------------------------------------------------------
	!						POCETAK PRORACUNA 
	! ---------------------------------------------------------------------------------
	mio = 4.d0*pi*1.e-7
	konst = mio/(4.d0*pi)
	omega = 2.d0*pi*f
	jomega = dcmplx(0.d0,omega)
	sumax = dcmplx(0.d0,0.d0)
	sumay = dcmplx(0.d0,0.d0)
	sumaz = dcmplx(0.d0,0.d0)
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

		call pot(x2,y2,z2,xt,yt,zt,L(iseg),xsi,ysi,zsi,pol,andi)

		sumax = sumax + ((x2-x1)/L(iseg)) * Iu(iseg) * andi * Fprig
		sumay = sumay + ((y2-y1)/L(iseg)) * Iu(iseg) * andi * Fprig
		sumaz = sumaz + ((z2-z1)/L(iseg)) * Iu(iseg) * andi * Fprig
	END DO
	Ax = jomega * konst * sumax  
	Ay = jomega * konst * sumay  
	Az = jomega * konst * sumaz  

	return
end subroutine