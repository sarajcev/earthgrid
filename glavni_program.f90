! ...
! ******************************************************************************

!						EartHGriD  Ver. 1.0  Linux 64bit

! ******************************************************************************

!                              GLAVNI PROGRAM

! -------------------------------------------------------------------------------
! DOUBLE PRECISION VERSION 
! -------------------------------------------------------------------------------
! PROGRAM ZA NUMERICKI PRORACUN HARMONICKOG STRUJNOG POLJA SUSTAVA UZEMLJIVACA
! POLOZENIH U VISESLOJNO TLO. VISESLOJNI MODEL UKLJUCUJE VISESLOJNO TLO I ZRAK.
! SEGMENTI SUSTAVA UZEMLJENJA MOGU SE NALAZITI BILO GDJE U TLU ILI PAK U ZRAKU,
! PRI CEMU MOGU BITI U BILO KOJEM POLOZAJU RELATIVNO U ODNOSU NA POVRSINU TLA.
! PROGRAM JE BAZIRAN NA PROSIRENJU METODE PREZENTIRANE U [1], UZ SAZNANJA IZ
! [2] I DODATNU NADOGRADNJU. OSNOVA JE METODA KONACNIH ELEMENATA PRIMJENJENA
! NA TANKO-ZICANU APROKSIMACIJU SEGMENATA SUSTAVA UZEMLJENJA. RJESAVA SE PRITOM
! INTEGRALNA FORMULACIJA PROBLEMA TEORIJE POLJA (HELMHOLTZOVA DIFERENCIJALNA
! JEDNADZBA) UZ ODREDJENE APROKSIMACIJE (IZBJEGAVANJE RJESAVANJA SOMMERFELDOVIH
! INTEGRALA - SVODJENJE NA POISSONOVU DIFERENCIJALNU JEDNADZBU). PRIGUSENJE JE
! PRITOM TRETIRANO NA POSEBAN NACIN. PRI FORMULACIJI SE KORISTI METODA SREDNJEG 
! POTENCIJALA I NUMERICKO ODSLIKAVANJE SEGMENATA U VISESLOJNOM SREDSTVU (NUMERI-
! CKA APROKSIMACIJA KERNEL FUNKCIJA TOCKASTOG IZVORA HARMONICKE STRUJE, [1]).
! PROGRAM JE NADOGRADJEN SA TRANZIJENTNIM PRORACUNOM SUSTAVA UZEMLJIVACA U VISE-
! SLOJNOM TLU UZ PRIMJENU ALGORITMA: FFT - MoM - IFFT, NPR. [3]. KORISTI SE PUNA
! FFT - INVERSE FFT TRANSFORMACIJA (TRANSFER FUNKCIJA SE RACUNA ZA SVE FREKVEN-
! CIJE - NEMA INTERPOLACIJE). PROGRAM KORISTI MNOGE NUMERICKE PAKETE IZ BOGATE
! FORTRANSKE BIBLIOTEKE (NETLIB REPOSITORY I SL.), [4-10]. 
! PROGRAMSKI PAKET "EartHGriD" NAPISAN JE KOMPLETNO U DVOSTRUKOJ PRECIZNOSTI 
! (DOUBLE PRECISION) REAL(8) I COMPLEX(8) S DINAMICKOM ALOKACIJOM MEMORIJE.

! [1] Vujevic, S.: "Kombinirani postupak proracuna uzemljivaca u krsevitom tlu",
!     Doktorska disertacija, FESB, Split, 1994.
! [2] Magazin, N.: "Elektromagnetski model za proracun harmonickog polja sustava
!     uzemljivaca", Magistarski rad, FESB, Split, 2003.
! [3] Grcev, L.: "Proracun tranzijentne impedancije uzemljivackih sistema",
!     Doktorska disertacija, FER, Zagreb, 1986.
! [4] BLAS: Level 1, Level 2 and Level 3 set of FORTRAN Basic Linear Algebra
!     Subprograms, Netlib, http://www.netlib.org/blas.
! [5] LAPACK 3.0: Linear Algebra PACKage, Netlib, http://www.netlib.org/lapack.
! [6] Visual Numerics: "IMSL", Fortran Subroutines for Mathematical Applica-
!     tions, Math/Library Vol. 1 and 2, Visual Numerics Inc., 1997.
! [7] Press, H. W.; Teukolsky, S. A.; Vetterling, W. T.; Flannery, B.P.: 
!     "Numerical Recipes in FORTRAN 77", The Art of Scientific Computing, Second
!     Edition, ISBN 0-521-43064-X, Cambridge University Press, 1995.
! [8] Swarztrauber, P. N.: "FFTPACK", National Center for Atmospheric Research,
!     Boulder, CO, USA.
! [9] SLATEC: Common Mathematical Library, Version 4.1, July 1993, http://www.
!     netlib.org/slatec.
![10] Boisvert, R. F.; Howe, S. E.; Kahaner, D. K.: "CMLIB" - NIST Core Math 
!     Library Ver. 3.0, http://lib.stat.cmu.edu/cmlib/index/index.

! PROGRAMSKI PAKET NAPISAN JE U PROGRAMSKOM JEZIKU "VISUAL FORTRAN F90/95" I 
! KORISTI SVE NJEGOVE PREDNOSTI (NPR. COLLON NOTATION, COMPLEX ARITHMETIC,
! DYNAMIC ARRAY ALLOCATION, ETC.).

! Author: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)
! -------------------------------------------------------------------------------

PROGRAM EartHGriD
	use funkcije
	implicit none
	
	! Definicija varijabli
	integer N,Ns,Nvi,Nsi,Nss
	real(8),dimension(:),allocatable :: ro,epsr,h
	real(8),dimension(:),allocatable :: xp,yp,zp,xk,yk,zk,r0,mirv,sigmav
	integer Nc,Ni
	real(8) rei,imi,f
	integer,dimension(:),allocatable :: icv
	complex(8),dimension(:),allocatable :: Ii
	integer,dimension(:),allocatable :: iveze
	
	integer,dimension(:),allocatable :: iso
	real(8),dimension(:),allocatable :: HD,L
	complex(8),dimension(:),allocatable :: kapa,gama,Fr,zunjv
	complex(8),dimension(:,:),allocatable :: Yuv,Ypv,Yu,Ypp,Yg
	complex(8),dimension(:,:),allocatable :: Au,Ap,CC
	complex(8),dimension(:),allocatable :: Iu,Ip,Fi,Fig,Fi1,Fi2
	real(8) xt,yt,zt
	complex(8) potencijal,potenc_freq
	real(8) modul,kut
	real(8),parameter :: pi = 3.141592653589793238d0
	integer sloj,Ntoc
	real(8) xp1,yp1,zp1,xk1,yk1,zk1,dx,dy,dz,t

	! Ostale varijable
	complex(8),dimension(:,:),allocatable :: B
	real(8) zsi
	character(128) ime_file
	integer Nv
	integer NPVI,NBS,LL
	real(8) X1,Y1,Z1,X2,Y2,Z2,POL,XD,YD,ZD
	real(8) mirvv,sigmavv
	integer i,j,ig,jg
	complex(8) Zuzem
	integer izbor
	real(8) xo,yo,zo,M1,M2
	real(8),dimension(:),allocatable :: Step,xx
	integer ik,korak,Nm
	real(8),dimension(:,:),allocatable :: FF
	real(8),dimension(:),allocatable :: Xos,Yos
	character(1) paralelan
	real(8) xtp,ytp,razmak
	integer Nprav,status,hv
	real(8),dimension(:,:),allocatable :: Lk
	real(8),dimension(:),allocatable :: Gx,Gy,Gz,Gt,xii,yii,zii
	real(8) xL,yL,zL,eps,xi,yi,zi
	integer ngc,pom
	integer,dimension(:),allocatable :: ipiv
	integer nrhs,info
	character(1) transa,transb,izb
	complex(8) alpha,betha,Nap_dod,FFig
	real(8),dimension(:,:),allocatable :: napon_koraka
	logical ponovi
	integer res,nG,broj_cvora
	
	! Elektricno i magentsko polje
	complex(8) Fix,Fiy,Fiz,Ax,Ay,Az,Ex,Ey,Ez
	complex(8),dimension(:),allocatable :: Exx,Eyy,Ezz
	real(8) modulx,moduly,modulz,mio
	complex(8) Bx,By,Bz
	complex(8),dimension(:),allocatable :: Bxx,Byy,Bzz
	
	! Tranzijentna analiza
	integer oblik_struje
	integer Nsample,N2,kk
	real(8) Tmax,deltaT,df
	real(8),dimension(:),allocatable :: ti,func
	integer kr_1,kr_2,kr_3
	real(8) Uo_1,Uo_2,Uo_3,T1_1,T2_1,a_2,b_2,a_3,omeg,psi_3,f_3
	real(8),dimension(:),allocatable :: wsave,freq
	complex(8),dimension(:),allocatable :: ccc,trans_pot
	complex(8),dimension(:,:),allocatable :: TransPot,TransPotP
	integer Ntip,tipv,tip_pr,IN
	integer,dimension(:),allocatable :: tipvod
	real(8),dimension(:),allocatable :: r0vod,mirvvod,sigmavvod
	
	! Mjerenje vremena
	real(8) start1,end1,time1,time2,time3,time4,time5,time6,time7
	real(8) time8,time9,time10,time11,time12,time13,time14,time15
	real(8) time16,time22,time33,time44,time55,time66,time77,time_uk
	integer ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,ti9,ti10,ti11,ti12
	integer ti13,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8
	
	complex(8),dimension(:),allocatable :: jomega

	! Vanjske subroutine - FFTPACK
	external zffti,zfftf,zfftb

	! Opis koristenih varijabli:
	!  ......................................................................................
	!  N - ukupni broj slojeva modela (viseslojno tlo plus zrak)
	!  ro (realni vektor) - specificna elektricna otpornost svakog od slojeva modela, [Ohm*m]
	!  epsr (realni vektor) - specificna relativna permitivnost svakog od slojeva modela
	!  h (realan vektor) - debljina svakog od slojeva modela, [m]
	!
	!  Nv - ukupni broj vodica sustava uzemljenja
	!  Ns - ukupni broj segmenata sustava uzemljenja
	!  Nvi - broj izoliranih vodica u sustavu uzemljenja
	!  Nsi - ukupni broj izoliranih segmenata u sustavu uzemljenja
	!  Napomena: Broj neizoliranih vodica i segmenata sustava uzemljenja racuna
	!  se respektivno kao: Nv-Nvi, odnosno, Ns-Nsi!
	!  X1 - x koordinata pocetne tocke vodica, [m]
	!  Y1 - y koordinata pocetne tocke vodica, [m]
	!  Z1 - z koordinata pocetne tocke vodica, [m]
	!  X2 - x koordinata krajnje tocke vodica, [m]
	!  Y2 - y koordinata krajnje tocke vodica, [m]
	!  Z2 - z koordinata krajnje tocke vodica, [m]
	!  POL - radijus (polumjer) vodica, [m]
	!  xp (realni vektor) - x koordinata pocetne tocke svih segmenata sustava uzemljenja, [m]
	!  yp (realni vektor) - y koordinata pocetne tocke svih segmenata sustava uzemljenja, [m]
	!  zp (realni vektor) - z koordinata pocetne tocke svih segmenata sustava uzemljenja, [m]
	!  xk (realni vektor) - x koordinata krajnje tocke svih segmenata sustava uzemljenja, [m]
	!  yk (realni vektor) - y koordinata krajnje tocke svih segmenata sustava uzemljenja, [m]
	!  zk (realni vektor) - z koordinata krajnje tocke svih segmenata sustava uzemljenja, [m]
	!  r0 (realni vektor) - radijus svakog od segmenata sustava uzemljenja, [m].
	!  mirv (realni vektor) - relativna permeabilnost svakog od segmenata uzemljivaca
	!  sigmav (realni vektor) - specificna elektricna vodljivost svakog od segmenata, [S/m]
	!  Nc - ukupni (stvarni) broj globalnih cvorova uzemljivackog sustava
	!  Ni - ukupni broj cvorova (globalni cv.) u kojima su injektirane struje
	!  icv (realni vektor) - globalni cvor u kojem je injektirana pojedina struja
	!  Ii (kompleksni vektor) - vrijednosti narinutih struja (Re + jIm), [A]
	!  iveze (cjelobr. vektor) - vektor veze lokalnih i globalnih cvorova
	!  ......................................................................................
	!  iso (cjelobrojni vektor) - sloj u kojem se nalazi svaki od segmenata uzemljivaca
	!  HD (realan vektor) - z koordinate donjih granicnih ploha svih slojeva modela, [m]
	!  kapa (kompleksni vektor) - specificna kompleksna vodljivost svakog od slojeva modela
	!  gama (kompleksni vektor) - valna konstanta svakog od slojeva modela
	!  L (relani vektor) - duljina svakog od segmenata sustava uzemljivaca, [m]
		
	! Napomena: Korisnik mora osigurati da se svaki segment uzemljivaca nalazi tocno 
	!			unutar granica jednog sloja, jer se to naknadno u programu ne provjerava.
		
	! ...........................................................................
	!		             PRIMJER FORME ULAZNE DATOTEKE
	! ...........................................................................
	! N=4
	! 0.0 1.0 0.0
	! ro2 epsr2 h2
	! ro3 epsr3 h3
	! ro4 epsr4 0.0
	! Ntip=2
	! 1 r0 mirv sigmavv
	! 2 r0 mirv sigmavv
	! Nv=6 Ns=6 Nvi=2 Nsi=2
	! X1 Y1 Z1 X2 Y2 Z2 1 1
	! X1 Y1 Z1 X2 Y2 Z2 1 1
	! X1 Y1 Z1 X2 Y2 Z2 2 1
	! X1 Y1 Z1 X2 Y2 Z2 2 1
	! X1 Y1 Z1 X2 Y2 Z2 1 1
	! X1 Y1 Z1 X2 Y2 Z2 1 1
	! tip_pr
	! Ni=1 f
	! xi yi zi Re Im
	! ...........................................................................

	! Objasnjenje primjera ulazne datoteke:
	! -------------------------------------
	! Ulazna datoteka moze se podijeliti u cetiri (4) funkcionalno odvojene kategorije,
	! kako slijedi:
	
	! (1) - podaci o visesloju:
	! Ulazna datoteka sadrzi primjer cetverosloja (zrak + troslojno tlo) -> N=4.
	! U drugom retku "ro" zraka je beskonacan pa se uvijek stavlja 0.0, kao i debljina
	! sloja (zraka) te se takodjer uvijek stavlja 0.0! Relativna permitivnost zraka 
	! je 1 (epsr zraka = 1). Debljina posljednjeg sloja (u petom retku ulazne datoteke)
	! je uvijek beskonacna pa se stavlja uvijek 0.0. 

	! Napomena: Program vodi racuna da su potonje vrijednosti beskonacne. Ukratko, 
	!           drugi redak ulazne datoteke ce uvijek imati vrijednosti: 0.0 1.0 0.0!

	! (2) - podaci o vodicima uzemljivaca:
	! Zadaje se broj tipova vodica od kojih je sastavljen uzemljivacki sustav (Ntip). 
	! Za svaki tip vodica zadaje se njegov redni broj, polumjer (r0), relativna 
	! permeabilnost materijala vodica (mirvv) i vodljivost materijala od kojeg je
	! izradjen doticni vodic (sigmavv). U gornjem primjeru postoje 2 tipa vodica, 
	! svaki sa svojim karakteristikama.

	! (3) - geometrisjki podaci uzemljivaca:
	! Uzemljivacki sustav moze se sastojati od neizoliranih (golih) i izoliranih
	! vodica. Zadaje se prvo ukupni broj svih vodica uzemljivackog sustava, te
	! ukupni broj svih segmenata od kojih se sastoji (na koji se dijele vodici)
	! doticni uzemljivacki sustav. Potom se u istom redku zadaje broj izoliranih
	! vodica, te broj izoliranih segmenata (na koji se izolirani vodici dijele).

	! U ovom primjeru uzemljivacki sustav sastoji se od cetiri neizolirana vodica i
	! dva izolirana vodica, dakle, Nv=6; Ns=6. Dakle, u ovom primjeru vodici uzem-
	! ljivaca se dodatno ne dijele na segmente (jer je stavljeno da je broj segmenata
	! jednak broju vodica). U protivnom varijablu Ns treba postaviti na ukupni broj
	! segmenata sustava uzemljenja. (Naime, svaki od vodica uzemljivaca moze se 
	! dijeliti na segmente). Ukupni broj svih segmenata (kad se sumiraju sve podjele
	! na segmente svih vodica uzemljivaca) definira upravo varijabla Ns.

	! Za svaki vodic zadaju se njegove globalne koordinate pocetne tocke (X1,Y1,Z1)
	! i globalne koordinate krajnje tocke (X2,Y2,Z2). Potom se zadaje vrijednost 
	! koja definira tip vodica (prema prethodno definiranim tipovima vodica). U 
	! konkretnom slucaju postoje dva tipa vodica i prva dva vodica uzemljivaca
	! pripadaju prvom tipu (1) dok druga dva vodica pripadaju drugom (2) tipu.

	! Napomena: Prvo se zadaju geometrijski podaci neizoliranih (golih) vodica kojih
	!			je ujedno i najvise u nekom uzemljivackom sustavu. Potom se na potpuno
	!			isti nacin zadaju podaci izoliranih vodica koji se eventualno nalaze
	!			u uzemljivackom sustavu. Ukoliko ovakvi vodici ne postoje, odgovarajuce
	!			varijable Nvi, Nsi se postavljaju na vrijednost nula!
	 
	! Konacno se zadaje varijabla koja definira na koliko se segmenata dijeli
	! doticni vodic (u ovom slucaju je to 1, jer se vodic dodatno ne dijeli na
	! segmente). U slucaju da se vodic dijeli na segmente, ovaj broj bi predstavljao
	! broj segmenata na koji ce se doticni vodic podijeliti (varijabla Ns bi pritom
	! morala imati vrijednost ukupnog broja svih segmenata od svih podjela vodica).

	! Napomena: Moguce je tretirati istovremeno vise glavanski odvojenih uzemljivackih
	!           sustava, bez posebnog oznacavanja kojem uzemljivacu pripada koji
	!           pojedini vodic (to se automatski definira pri asembliranju). Vodici
	!           vise razlicitih uzemljivaca zadaju se proizvoljnim redosljedom.

	! (4) - vrsta proracuna i narinuta struja:
	! Zadaje se varijabla koja definira tip/vrstu proracuna (tip_pr), kako slijedi:
	!    tip_pr = 1 -> provodi se proracun za samo jednu frekvenciju
	!    tip_pr = 2 -> provodi se tranzijentni proracun koristenjem FFT - IFFT
	
	! Napomena: Ukoliko je tip_pr=2, tada se za frekvenciju proracuna u nastavku
	!           zadavanja ulaznih podataka obicno stavlja vrijednost: 0.0 ili bilo
	!           koja druga vrijednost frekvencije (doticni podatak se ignorira). 
	!			Frekvencije po kojima ce se racunati dobivaju se FFT analizom.

	! Napomena: Ukoliko je tip_pr=2, tada se za narinutu struju obavezno stavlja
	!           vrijednost 1.0 + j0.0! S ovom strujom odredjuje se transfer funkci-
	!			ja u frekventnoj domeni (unmodulated transfer function), koja se 
	!           potom mnozi s faktorima od FFT proracuna (modulated transfer functi-
	!			on), te se konacno IFFT proracunom vraca u vremensku domenu. U ovom 
	!           slucaju broj struja (Ni) jednak je jedinici (Ni=1)!!!

	! Zadaje se ukupni broj cvorova uzemljivackog sustava (Ni) u kojima su
	! injektirane struje, kao i frekvencija injektiranih struja (tip_pr=2).
	! Potom se za svaki cvor u kojem je injektirana struja zadaju globalne
	! koordinate doticnog cvora, te potom i realni i imaginarni dio
	! injektirane struje. U ovom slucaju broj injektiranih struja je 1 (Ni=1),
	! te postoji samo jedan redak (xi yi zi Re Im). U protivnom, postojalo bi 
	! onoliko redaka koliki je broj Ni (koliko se struja injektira). Odredjene
	! struje koje ulaze u uzemljivac imaju pozitivni predznak, dok negativni
	! predznak struje pretpostavlja da struje napustaju uzemljivacki sustav.
	
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

	! ===============================================================================
	!                        POCETAK GLAVNOG PROGRAMA
	! ===============================================================================
	write(*,'(/,tr15,"EartHGriD Linux Ver. 1.0",/)')
	write(*,'("PROGRAM ZA  NUMERICKI PRORACUN  HARMONICKOG POLJA SUSTAVA")')
	write(*,'("UZEMLJIVACA POLOZENIH U VISESLOJNO SREDSTVO - STACIONARNA")')
	write(*,'(7X,"I TRANZIJENTNA (FFT - MoM - IFFT) ANALIZA")')
	write(*,'(/,"Autor: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)")')
	
	open(2,file='output.txt')
	write(2,'(/,tr15,"EartHGriD Linux Ver. 1.0",/)')
	write(2,'("PROGRAM ZA  NUMERICKI PRORACUN  HARMONICKOG POLJA SUSTAVA")')
	write(2,'("UZEMLJIVACA POLOZENIH U VISESLOJNO SREDSTVO - STACIONARNA")')
	write(2,'(7X,"I TRANZIJENTNA (FFT - MoM - IFFT) ANALIZA")')
	write(2,'(/,"Autor: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)")')
	
	! ---------------------------------------------------------------------------
	!       CITANJE ULAZNIH PODATAKA IZ VANJSKOG FILE-A (INPUT)
	!       I TISKANJE ULAZNIH PODATAKA U DRUGI VANJSKI FILE (OUTPUT)
	! ---------------------------------------------------------------------------
	! Tiskanje naslova programskog paketa, GPL licence i forme ulazne datoteke
	call naslov(ime_file)
	
	open(unit=1,file=ime_file)
	! Podaci o viseslojnom sredstvu
	! ------------------------------
	read(1,*) N
	allocate(ro(N),epsr(N),h(N))
	do i = 1,N
		read(1,*) ro(i),epsr(i),h(i)
	end do

	write(2,'(//,T5,"U L A Z N I   P O D A C I")')
	write(2,'(/,"Podaci o sredstvu:")')
	write(2,'(18("-"))')
	write(2,'("Ukupni broj slojeva sredstva (zrak + viseslojno tlo):",i2)') N
	write(2,*)
	write(2,'("Br.",7x,"Epsr",5x,"Ro",9x,"h")')
	write(2,'("sloja",3x,4x,4x,"[Ohm*m]",6x,"[m]")')
	write(2,'("-----------------------------------")')
	do i = 1,n
		write(2,'(i5,1x,f7.1,3x,f8.2,2x,f8.2)') i,epsr(i),ro(i),h(i)
	end do
	write(2,'("-----------------------------------")')
	write(2,'("Napomena:",/,&
	"Varijable u gornjoj tablici imaju sljedeca znacenja:",/,&
	"epsr - specificna permeabilnost sloja",/,&
	"ro   - specificna elektricna otpornost slojeva",/,&
	"h    - debljina sloja, u metrima.")')
	write(2,*)

	! Podaci o tipovima vodica uzemljivaca
	! -------------------------------------
	read(1,*) Ntip
	allocate(tipvod(Ntip))
	allocate(r0vod(Ntip))
	allocate(mirvvod(Ntip))
	allocate(sigmavvod(Ntip))
	do i = 1,Ntip
		read(1,*) tipvod(i),r0vod(i),mirvvod(i),sigmavvod(i)
	end do

	write(2,'(/,"Podaci o tipovima vodica uzemljivaca:")')
	write(2,'(44("-"))')
	write(2,'("Br.",7x,"r0 (m)",6x,"mir",7x,"sigma (S/m)")')
	write(2,'(44("-"))')
	do i = 1,Ntip
		write(2,'(i2,5x,f10.6,2x,f6.1,6x,e12.6)') tipvod(i),&
		r0vod(i),mirvvod(i),sigmavvod(i)
	end do

	! Podaci o uzemljivackim sustavima
	! ---------------------------------
	read(1,*) Nv,Ns,Nvi,Nsi
	allocate(xp(Ns),yp(Ns),zp(Ns))
	allocate(xk(Ns),yk(Ns),zk(Ns))
	allocate(r0(Ns),mirv(Ns),sigmav(Ns))
	NBS = 0

!	==============
	Nss = Ns - Nsi	!BROJ NEIZOLIRANIH SEGMENATA
!	==============

	! Tiskanje podataka o uzemljivackom sustavu u vanjski file
	write(2,'(/,"Podaci o uzemljivackom sustavu:")')
	write(2,'(31("-"))')
	write(2,'("Ukupni broj vodica svih uzemljivaca: ",i2)') Nv
	write(2,'("Ukupni broj segmenata svih uzemljivaca: ",i3)') Ns
	write(2,'("Ukupni broj izoliranih vodica uzemljivaca: ",i3)') Nvi
	write(2,'("Ukupni broj izoliranih segmenata uzemljivaca: ",i3)') Nsi
	write(2,*)
	write(2,'("Elem.",T12,"Xp (m)",T21,"Yp (m)",T30,"Zp (m)",T40,&
	"Xk (m)",T50,"Yk (m)",T60,"Zk (m)",T70,"r0 (m)",T80,"mirv",T90,&
	"sigmav",T100,"Seg.")')
	write(2,'(104("-"))')

	open(4,file='tempp.dat') ! Privremena datoteka
	DO i = 1,Nv
		read(1,*) X1,Y1,Z1,X2,Y2,Z2,tipv,NPVI
		do j = 1,Ntip
			if (tipv==tipvod(j)) then
				POL = r0vod(j)
				mirvv = mirvvod(j)
				sigmavv = sigmavvod(j)
			end if
		end do
		write(4,'(6f12.3)') x1,y1,z1,x2,y2,z2
		write(2,'(I4,T10,F6.2,T19,F6.2,T28,F6.2,T38,F6.2,T48,F6.2,T58,&
		F6.2,T69,F6.4,T78,F6.1,T88,e10.3,T100,i3)') i,x1,y1,z1,x2,y2,z2,&
		pol,mirvv,sigmavv,npvi
         XD = (X2-X1)/NPVI                 
         YD = (Y2-Y1)/NPVI                 
         ZD = (Z2-Z1)/NPVI                 
         DO LL = 1,NPVI                           
            NBS = NBS+1
			r0(NBS) = pol
			mirv(NBS) = mirvv
			sigmav(NBS) = sigmavv
            XP(NBS) = X1 + XD*(LL-1)                     
            XK(NBS) = XP(NBS) + XD
            YP(NBS) = Y1 + YD*(LL-1)
            YK(NBS) = YP(NBS) + YD
            ZP(NBS) = Z1 + ZD*(LL-1)
            ZK(NBS) = ZP(NBS) + ZD
         END DO
	END DO
	close(4)
	write(2,'(104("-"))')
	deallocate(tipvod)
	deallocate(r0vod)
	deallocate(mirvvod)
	deallocate(sigmavvod)

	! Tip proracuna koji ce se izvrsavati
	! ------------------------------------
	read(1,*) tip_pr

	! Podaci o narinutim strujama
	! ----------------------------
	read(1,*) Ni,f
	allocate(Ii(Ni))
	allocate(icv(Ni))
	allocate(xii(Ni),yii(Ni),zii(Ni))
	do i = 1,Ni
		read(1,*) xi,yi,zi,rei,imi
		Ii(i) = dcmplx(rei,imi)
		xii(i) = xi; yii(i) = yi; zii(i) = zi
	end do
	close(1)

	IF (tip_pr==2) THEN
		! ***************************************************************************
		!          TRANZIJENTNA ANALIZA UZEMLJIVACKIH SUSTAVA U VISESLOJU
		!                        (FFT - MoM - IFFT)
		! ***************************************************************************
		write(*,'(/,"Tranzijentna analiza uzemljivaca...")')
		! Brisem privremenu datoteku: tempp.dat jer u ovom slucaju
		! (tranzijentna analiza) ista nije potrebna!
		!res = DELFILESQQ('tempp.dat')

		write(*,'(/,"Odaberite oblik tranzijentne struje:")')
		write(*,'(59("-"))')
		write(*,'(4x,"1 - Oblik struje tipa ramp (npr. povratni preskok)")')
		write(*,'(4x,"2 - Dvostruko exponencijalni oblik struje (Lightning)")')
		write(*,'(4x,"3 - Prigusena sinusoida oblik struje (Sklopni prenapon)")')
		write(*,'(4x,"4 - Proracun Y(s): odabire se u nastavku proracun (8)")')
		write(*,'(4x,"0 - Exit (izlaz iz programa)")')
		write(*,'(59("-"))')
		write(*,'("Izbor: ")',advance='no')
		read(*,*) oblik_struje

		select case(oblik_struje)
			case(1)
				! Struja ima oblik "ramp type"
				! -----------------------------
				write(*,'(/,"Parametri ramp type struje:")')
				write(*,'("Amplituda (A): ")',advance='no')
				read(*,*) Uo_1
				write(*,'("Trajanje cela: T1 (s): ")',advance='no')
				read(*,*) T1_1
				write(*,'("Trajanje hrbta: T2 (s): ")',advance='no')
				read(*,*) T2_1
				write(*,'("Parametar kr: ")',advance='no')
				read(*,*) kr_1
				write(*,'(/,"Broj samplova: N = ")',advance='no')
				read(*,*) Nsample
				write(*,'("Ukupno vrijeme promatranja: Tmax (s): ")',advance='no')
				read(*,*) Tmax
				!Ramp type function
				deltaT = Tmax/Nsample
				allocate(ti(Nsample))
				ti(1) = 0.d0
				do i = 2,Nsample
					ti(i) = ti(i-1) + deltaT
				end do
				allocate(func(Nsample))
				do i = 1,Nsample
					if(ti(i)>0.0 .and. ti(i)<T1_1)then
						func(i)=Uo_1/T1_1*ti(i)
					else if(ti(i)>T1_1 .and. ti(i)<T2_1) then
						if(kr_1==1) then
							func(i)=Uo_1
						else if(kr_1==0) then
							func(i)=Uo_1-Uo_1/(T2_1-T1_1)*(ti(i)-T1_1)
						end if
					end if
				end do
			case(2)
				! Struja ima oblik "double exponential"
				! --------------------------------------
				write(*,'(/,"Parametri dvostruko exponencijalne struje:")')
				write(*,'("Amplituda (A): ")',advance='no')
				read(*,*) Uo_2
				write(*,'("Parametar: a (1/s): ")',advance='no')
				read(*,*) a_2
				write(*,'("Parametar: b (1/s): ")',advance='no')
				read(*,*) b_2
				write(*,'("Koeficijent kr: ")',advance='no')
				read(*,*) kr_2
				write(*,'(/,"Broj samplova: N = ")',advance='no')
				read(*,*) Nsample
				write(*,'("Ukupno vrijeme promatranja: Tmax (s): ")',advance='no')
				read(*,*) Tmax
				! Double exponential function
				deltaT = Tmax/Nsample
				allocate(ti(Nsample))
				ti(1) = 0.d0
				do i = 2,Nsample
					ti(i) = ti(i-1) + deltaT
				end do
				allocate(func(Nsample))
				do i = 1,Nsample
					if(kr_2==0) then
						func(i)=Uo_2*dexp(-a_2*ti(i))
					else if(kr_2==1) then
						func(i)=Uo_2*(dexp(-a_2*ti(i))-dexp(-b_2*ti(i)))
					end if
				end do
			case(3)
				! Struja ima oblik "dumped sinusoid"
				! -----------------------------------
				write(*,'(/,"Parametri oblika struje sklopnog prenapona:")')
				write(*,'("Amplituda (A): ")',advance='no')
				read(*,*) Uo_3
				write(*,'("Parametar: a (1/s): ")',advance='no')
				read(*,*) a_3
				write(*,'("Parametar psi: ")',advance='no')
				read(*,*) psi_3
				write(*,'("Frekvencija (Hz): ")',advance='no')
				read(*,*) f_3
				write(*,'("Koeficijent kr: ")',advance='no')
				read(*,*) kr_3
				write(*,'(/,"Broj samplova: N = ")',advance='no')
				read(*,*) Nsample
				write(*,'("Ukupno vrijeme promatranja: Tmax (s): ")',advance='no')
				read(*,*) Tmax
				! Dumped sinusoid function
				deltaT = Tmax/Nsample
				allocate(ti(Nsample))
				ti(1) = 0.d0
				do i = 2,Nsample
					ti(i) = ti(i-1) + deltaT
				end do
				allocate(func(Nsample))
				omeg = 2.d0*pi*f_3
				do i = 1,Nsample
					func(i)=Uo_3*(kr_3*dexp(-a_3*ti(i))-dcos(omeg*ti(i)+psi_3*pi/180.))
				end do
			case(4)
				write(*,'(/,"No. of frequency samples: N = ")',advance='no')
				read(*,*) Nsample
				write(*,'(/,"Max. frequency of interest: fmax = ")',advance='no')
				read(*,*) Tmax !recikliranje
				df = Tmax/Nsample
			case(0)
				write(*,'("Exit!")')
				stop
		end select
		!Tiskanje ulaznih podataka
		write(2,'(/,"Forward FFT calculation:",/)')
		write(2,'("Number of samples: N = ",i10)') Nsample
		write(2,'("Time of observation: Tmax = ",e16.6," [s]")') Tmax
		write(2,'(/,"Surge option:",i3)') oblik_struje
		select case(oblik_struje)
			case(1)
				write(2,'("Lightning Surge (simple):")')
				write(2,'("Uo = ",f12.3)') Uo_1
				write(2,'("T1 = ",e16.6)') T1_1
				write(2,'("T2 = ",e16.6)') T2_1
				write(2,'("kr = ",i3)') kr_1
			case(2)
				write(2,'("Lightning Surge (double exponential):")')
				write(2,'("Uo = ",f12.3)') Uo_2
				write(2,'("a = ",e16.6)') a_2
				write(2,'("b = ",e16.6)') b_2
				write(2,'("kr = ",i3)') kr_2
			case(3)
				write(2,'("Switching Surge:")')
				write(2,'("Uo = ",f12.3)') Uo_3
				write(2,'("a = ",e16.6)') a_3
				write(2,'("f = ",f12.5)') f_3
				write(2,'("psi = ",f8.3)') psi_3
				write(2,'("kr = ",i3)') kr_3
			case(4)
				continue
		end select
		select case(oblik_struje)
			case(1:3)
				!Tiskanje u file za PlotXY
				open(3,file='signal_orig.adf')
				write(3,'(/,5x,"t(sec)",10x,"func(t)")')
				do i=1,Nsample
					write(3,'(f16.10,2x,f16.10)') ti(i),func(i)
				end do
				close(3)
		end select
		! Pocetak prvog mjerenja CPU vremena trajanja proracuna
		call cpu_time(start1)

		select case(oblik_struje)
		case(1:3)
			! --------------------------------------------------------------------------
			!                     Frekvencije za FFT proracun
			! --------------------------------------------------------------------------
			allocate(freq(Nsample/2+1))
			freq(1) = 1.d-6		!Istosmjerna (f=0)
			df = 1.d0/Tmax
			do i = 2,Nsample/2+1
				freq(i) = freq(i-1) + df
			end do

			allocate(ccc(Nsample))
			do i = 1,Nsample
				ccc(i) = dcmplx(func(i),0.d0)
			end do
			! ---------------------------------------------------------------------------
			!                    FORWARD FFT PRORACUN - FFTPACK
			! ---------------------------------------------------------------------------
			allocate(wsave(4*Nsample+15))
			! Poziv rutine za Forward FFT proracun
			! ------------------------------------
			call zffti(Nsample,wsave)
			call zfftf(Nsample,ccc,wsave)
			! ------------------------------------
			! Normiranje
			ccc = ccc/Nsample
		case(4)
			allocate(freq(Nsample+1))
			freq(1) = 1.d-6
			do i = 2,Nsample+1
				freq(i) = freq(i-1) + df
			end do
		end select

		! ---------------------------------------------------------------------------
		!               PRORACUN KONSTANTNIH PARAMETARA PROGRAMA
		! ---------------------------------------------------------------------------
		allocate(HD(N))
		! Proracun z koordinate donjih granicnih ploha slojeva modela
		call vektor_H(N,h,HD)
		allocate(L(Ns))
		! Proracun duljine svih segmenata
		do i = 1,Ns
			L(i) = dsqrt((xk(i)-xp(i))**2+(yk(i)-yp(i))**2+(zk(i)-zp(i))**2)
		end do
		allocate(iso(Ns))
		! Odredjivanje sloja u kojem se nalazi svaki od segmenata
		do i = 1,Ns
			zsi = (zp(i)+zk(i))/2.d0
			iso(i) = 1
			do j = 1,N-1
				if (zsi>HD(j)) iso(i) = iso(i) + 1
			end do
		end do

		! ---------------------------------------------------------------------------
		!               FORMIRANJE VEKTORA VEZE ZA ASEMBLIRANJE
		! ---------------------------------------------------------------------------
		allocate(iveze(2*Ns))
		allocate(Lk(3,2*Ns))
		Lk(1,1:Ns) = xp(:)
		Lk(1,(Ns+1):(2*Ns)) = xk(:)
		Lk(2,1:Ns) = yp(:)
		Lk(2,(Ns+1):(2*Ns)) = yk(:)
		Lk(3,1:Ns) = zp(:)
		Lk(3,(Ns+1):(2*Ns)) = zk(:)
		ngc = 1
		allocate(Gx(ngc))
		allocate(Gy(ngc))
		allocate(Gz(ngc))
		eps = minval(r0)	! Tolerancija
		Gx(1) = Lk(1,1)
		Gy(1) = Lk(2,1)
		Gz(1) = Lk(3,1)
		iveze(1) = 1
		do i = 2,2*Ns
			xL = Lk(1,i)
			yL = Lk(2,i)
			zL = Lk(3,i)
			pom = 0
			do j = 1,ngc
				if ((dabs(xL-Gx(j))<=eps).and.(dabs(yL-Gy(j))<=eps)&
					.and.(dabs(zL-Gz(j))<=eps)) then
					pom = j
				end if
			end do
			if (pom/=0) then
				iveze(i) = pom
			else
				ngc = ngc + 1
				! Varijabla x
				allocate(Gt(ngc-1))
				Gt = Gx
				deallocate(Gx)
				allocate(Gx(ngc))
				Gx(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gx(ngc) = Lk(1,i)
				! Varijabla y
				allocate(Gt(ngc-1))
				Gt = Gy
				deallocate(Gy)
				allocate(Gy(ngc))
				Gy(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gy(ngc) = Lk(2,i)
				! Varijabla z
				allocate(Gt(ngc-1))
				Gt = Gz
				deallocate(Gz)
				allocate(Gz(ngc))
				Gz(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gz(ngc) = Lk(3,i)
				iveze(i) = ngc
			end if
		end do
		deallocate(Lk)
		! ---------------------------------------------------------------------------
		!          Ukupni broj globalnih cvorova 
		! ---------------------------------------------------------------------------
		Nc = ngc   !; write(*,'("Br. globalnih cv.:",i3)') Nc
		! ---------------------------------------------------------------------------
		! Odredjivanje vrijednosti globalnog cvora u kojem je injektirana struja
		! ---------------------------------------------------------------------------
		do j = 1,Ni
			do i = 1,Nc
				if ((dabs(Gx(i)-xi)<=eps).and.(dabs(Gy(i)-yi)<=eps)&
					.and.(dabs(Gz(i)-zi)<=eps)) then
					icv(j) = i
				end if
			end do
		end do
		deallocate(Gx,Gy,Gz)
		deallocate(xii,yii,zii)
		! -----------------------------------------
		! Tiskanje ulaznih podataka u vanjski file
		! -----------------------------------------
		write(2,'(/,"Podaci o narinutim strujama:")')
		write(2,'(28("-"))')
		write(2,'("Broj globalnih cvorova u kojima su narinute struje: ",i2)') Ni
		write(2,'("Frekvencija narinutih struja:",f10.2," [Hz]")') f
		write(2,*)
		write(2,'("Cvor.",T10,"Re(I) [A]",T25,"Im(I) [A]")')
		write(2,'(35("-"))')
		do i = 1,Ni
			write(2,'(T3,i3,T10,f10.3,T20,f10.3)') icv(i),rei,imi
		end do
		write(2,'(35("-"))')
		write(2,'(//,T5,"R E Z U L T A T I   P R O R A C U N A")')
		
		! Kraj prvog mjerenja CPU vremena trajanja proracuna
		call cpu_time(end1)
		time1 = end1 - start1

		! ==============================================================================
		!		TRANZIJENTNI PRORACUNI OD INTERESA KORISNIKU PROGRAMSKOG PAKETA
		! ==============================================================================
		tt1 = 0; tt2 = 0; tt3 = 0; tt4 = 0
		tt5 = 0; tt6 = 0; tt7 = 0; tt8 = 0

		write(*,'(/,"Odaberite sto zelite racunati:")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'("Skalarni potencijal:")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'(4x,"1 - Tranzijentni potencijal u bilo kojoj tocki viseslojnog")')
		write(*,'("        sredstva (FFT-MoM-IFFT)")')
		write(*,'(4x,"2 - Raspodjelu tranzijentnog potencijala duz pravca u")')
		write(*,'("        viseslojnom sredstvu (3D perspektiva) FFT-MoM-IFFT")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'("Elektricno i magnetsko polje/indukcija:")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'(4x,"3 - Tranzijentno elektricno polje Ex, Ey i Ez, te ukupno")')
		write(*,'("        el. polje u samo jednoj tocki viseslojnog sredstva")')
		write(*,'("        (FFT-MoM-IFFT)")')
		write(*,'(4x,"4 - Tranzijentno magnetsko polje/indukcija u jednoj tocki")')
		write(*,'("        viseslojnog sredstva (FFT-MoM-IFFT)")')
		write(*,'(4x,"5 - Raspodjelu tranzijentnog magnetskog polja Hx, Hy i Hz,")')
		write(*,'("        te ukupnog mag. polja; raspodjela tranzijentne magn.")')
		write(*,'("        indukcije Bx, By, Bz i ukupne mag. indukcije duz pravca")')
		write(*,'("        u viseslojnom sredstvu (3D perspektiva) FFT-MoM-IFFT")')
		write(*,'(4x,"6 - Raspodjelu tranzijentnog elektricnog polja Ex, Ey, Ez,")')
		write(*,'("        te ukupnog el. polja duz pravca u viseslojnom sredstvu")')
		write(*,'("        (3D perspektiva) FFT-MoM-IFFT")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'("Ulazna impedancija:")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'(4x,"7 - Tranzijentnu ulaznu impedanciju uzemljivackog sustava")')
		write(*,'("        (FFT-MoM-IFFT)")') 
		write(*,'(4x,"8 - Frekvencijski ovisnu nadomjesnu admintanciju sustava")')
		write(*,'("        Y(s) za potrebe vector fitting analize")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'(4x,"0 - Zavrsetak proracuna i izlaz iz programa.")')
		write(*,'("--------------------------------------------------------------")')
		write(*,'("Izbor: ")',advance='no')
		read(*,*) izbor

		if (oblik_struje == 4) then
			izbor = 8
			!write(*,'(/,"Izbor mora biti (8)!")')
		end if

		IF (izbor==1) THEN
			! ---------------------------------------------------------------------------
			!            PRORACUN TRANZIJENTNOG POTENCIJALA U JEDNOJ TOCKI 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Unesite koordinate tocke za koju zelite izracunati potencijal:")')
			write(*,'("Xt = ")',advance='no')
			read(*,*) xt
			write(*,'("Yt = ")',advance='no')
			read(*,*) yt
			write(*,'("Zt = ")',advance='no')
			read(*,*) zt
			write(*,'("Sloj u kojem je tocka: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Racunam ...")')

			! Pocetak drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt1 = 1

			! ==========================================================================
			!    PRORACUN TRANZIJENTNOG POTENCIJALA U PRETHODNO ZADANOJ TOCKI
			!        U VISESLOJNOM SREDSTVU PO SVIM FREKVENCIJAMA OD FFT
			! ==========================================================================
			allocate(trans_pot(Nsample/2+1))
			DO kk = 1, Nsample/2+1
				! Frekvencija za koju se racuna
				! -----------------------------
				f = freq(kk)
				! Narinuta tranzijentna struja (od FFT proracuna)
				! -----------------------------------------------
				Ii(:) = ccc(kk)
				! Proracun tranzijentnog potencijala u jednoj tocki viseslojnog sredstva
				! za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati onoliko
				! puta koliko se razlicitih frekvencija zeli racunati.
				! .......................................................................
				call single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,h,&
				iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
				! .......................................................................
				trans_pot(kk) = potenc_freq		
				if (kk == nint((Nsample/2+1)/10.)) then
					write(*,'(/,"...10%")')
				else if (kk == nint(2*(Nsample/2+1)/10.)) then
					write(*,'("...20%")')
				else if (kk == nint(3*(Nsample/2+1)/10.)) then
					write(*,'("...30%")')
				else if (kk == nint(4*(Nsample/2+1)/10.)) then
					write(*,'("...40%")')
				else if (kk == nint(5*(Nsample/2+1)/10.)) then
					write(*,'("...50%")')
				else if (kk == nint(6*(Nsample/2+1)/10.)) then
					write(*,'("...60%")')
				else if (kk == nint(7*(Nsample+1)/10.)) then
					write(*,'("...70%")')
				else if (kk == nint(8*(Nsample/2+1)/10.)) then
					write(*,'("...80%")')
				else if (kk == nint(9*(Nsample/2+1)/10.)) then
					write(*,'("...90%")')
				else if (kk == Nsample/2) then
					write(*,'("..100%")')
				end if
			END DO

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Priprema vektora: ccc za proracun IFFT. Prva polovica + jedan clan vektora
			! ccc su vrijednosti tranzijentnog potencijala u odabranoj tocki za sve pozi-
			! tivne frekvencije + Nyquist-ovu frekvenciju. Druga polovica vektora: ccc je
			! konjugirano-kompleksna i simetricna s obzirom na Nyquist-ovu frekvenciju.
			! Formiran vektor ccc se transformira u vremensku domenu pomocu IFFT algoritma.
			N2 = Nsample/2+1
			! Istosmjerni clan
			ccc(1) = trans_pot(1)
			! Pozitivne frekvencije
			do i = 2,Nsample/2 
				ccc(i) = trans_pot(i)
			end do
			! Nyquist-ova frekvencija
			ccc(N2) = trans_pot(N2)
			! Negativne frekvencije
			do i = 2,Nsample/2
				IN = Nsample - i + 2
				ccc(IN) = dconjg(ccc(i))
			end do
			deallocate(trans_pot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			! Poziv subroutine za Inverse FFT
			! -------------------------------
			call zfftb(Nsample,ccc,wsave)
			! -------------------------------
			deallocate(wsave)

			!Tiskanje signala (nakon Inverse FFT)
			open(3,file="signal_potenc.adf")
			write(3,'(/,5x,"t",10x,"Re[f(t)]",10x,"Im[f(t)]")')
			do i=1,Nsample
				write(3,'(f16.10,2x,e16.10,2x,e16.10)') ti(i),dreal(ccc(i)),dimag(ccc(i))
			end do
			close(3)
			deallocate(ccc)

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

			! Kraj drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time2 = end1 - start1

			write(*,'(/,"Napomena:")')
			write(*,'("Tranzijentni potencijal u odabranoj tocki visesloja snimljen")')
			write(*,'("je u file: signal_potenc.adf Za prikaz koristiti: PlotXY.")')
			
		ELSE IF(izbor==2) THEN
			! ---------------------------------------------------------------------------
			!            PRORACUN TRANZIJENTNOG POTENCIJALA DUZ PRAVCA (3D) 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Koordinate pocetne tocke pravca:")')
			write(*,'("Xp = ")',advance='no')
			read(*,*) xp1
			write(*,'("Yp = ")',advance='no')
			read(*,*) yp1
			write(*,'("Zp = ")',advance='no')
			read(*,*) zp1
			write(*,'("Koordinate krajnje tocke pravca:")')
			write(*,'("Xk = ")',advance='no')
			read(*,*) xk1
			write(*,'("Yk = ")',advance='no')
			read(*,*) yk1
			write(*,'("Zk = ")',advance='no')
			read(*,*) zk1
			write(*,'("Korak duz osi x: ")',advance='no')
			read(*,*) dx
			write(*,'("Korak duz osi y: ")',advance='no')
			read(*,*) dy
			write(*,'("Korak duz osi z: ")',advance='no')
			read(*,*) dz
			write(*,'("Sloj u kojem je pravac: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Racunam ...")')
				
			! Pocetak mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt2 = 1
				
			if (dx/=0.0) then
				Ntoc = (xk1-xp1)/dx
			else if(dy/=0.0) then
				Ntoc = (yk1-yp1)/dy
			else if(dz/=0.0) then
				Ntoc = (zk1-zp1)/dz
			end if

			allocate(TransPot(Nsample/2+1,Ntoc)) 
			allocate(trans_pot(Nsample/2+1))
			xt = xp1
			yt = yp1
			zt = zp1
			do i = 1,Ntoc
				write(*,'(" -- procesiram vrijednost: ",i3," / ",i3)') i,Ntoc
				! ==========================================================================
				!    PRORACUN TRANSFER FUNKCIJE POTENCIJALA U PRETHODNO ZADANOJ TOCKI NA
				!    ZADANOJ POVRSINI (FREQUENCY RESPONSE) - PO SVIM FREKVENCIJAMA OD FFT
				! ==========================================================================
				DO kk = 1, Nsample/2+1
					! Frekvencija za koju se racuna
					! -----------------------------
					f = freq(kk)
					! Narinuta tranzijentna struja (od FFT proracuna)
					! -----------------------------------------------
					Ii(:) = ccc(kk)
					! Proracun tranzijentnog potencijala u jednoj tocki viseslojnog sredstva
					! za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati onoliko
					! puta koliko se razlicitih frekvencija zeli racunati.
					! .......................................................................
					call single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,h,&
					iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
					! .......................................................................
					trans_pot(kk) = potenc_freq		
				END DO
				TransPot(:,i) = trans_pot(:)
					
				xt = xt + dx
				yt = yt + dy
				zt = zt + dz
			end do
			deallocate(trans_pot)

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Prva polovica vektora + jedan clan vektora su vrijednosti tranzijentnog 
			! potencijala u odabranoj tocki za sve pozitivne frekvencije + Nyquist-ovu
			! frekvenciju. Druga polovica vektora je konjugirano-kompleksna i simetricna
			! s obzirom na Nyquist-ovu frekvenciju. Formiran vektor se transformira u 
			! vremensku domenu pomocu IFFT algoritma (FFTPACK).
			N2 = Nsample/2+1
			allocate(TransPotP(Nsample,Ntoc))
			DO kk = 1,Ntoc
				! Istosmjerni clan
				TransPotP(1,kk) = TransPot(1,kk)
				! Pozitivne frekvencije
				do i = 2,Nsample/2 
					TransPotP(i,kk) = TransPot(i,kk)
				end do
				! Nyquist-ova frekvencija
				TransPotP(N2,kk) = TransPot(N2,kk)
				! Negativne frekvencije
				do i = 2,Nsample/2
					IN = Nsample - i + 2
					TransPotP(IN,kk) = dconjg(TransPot(i,kk))
				end do
			END DO
			deallocate(TransPot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			allocate(FF(Ntoc,Nsample)) !Recikliranje varijable FF
			DO kk = 1,Ntoc
				ccc(:) = TransPotP(:,kk)
				! Poziv subroutine za Inverse FFT
				!-------------------------------
				call zfftb(Nsample,ccc,wsave)
				!-------------------------------
				FF(kk,:) = dreal(ccc(:))
			END DO
			deallocate(ccc,wsave)
			deallocate(TransPotP)

			! Kraj mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time22 = end1 - start1

			!-----------------------------------------------------------------------
			!       3D prikaz izracunatih potencijala - Gnuplot
			!-----------------------------------------------------------------------
			open(7,file='3Dtranspot.dat')
			do j = 1,Ntoc
				do i = 1,Nsample
					write(7,'(e12.4)') FF(j,i)
				end do
				write(7,*)
			end do
			close(7)
			!-----------------------------------------------------------------------
			!       3D prikaz izracunatih potencijala - LabPlot
			!-----------------------------------------------------------------------
			open(7,file='3Dtrans_pot.dat')
			do i = 1,Nsample
				write(7,'(80e12.4)') (FF(j,i),j=1,Ntoc)
			end do
			close(7)						
			deallocate(FF)
			
			write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele tranzijentnog potencijala snimljeni su u file:")')
			write(*,'("3Dtranspot.dat. Za prikaz koristiti: Gnuplot. Za LabPlot: 3Dtrans_pot.dat")')

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

		ELSE IF (izbor==3) THEN
			! ---------------------------------------------------------------------------
			!       PRORACUN TRANZIJENTNOG ELEKTRICNOG POLJA U JEDNOJ TOCKI 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Unesite koord. tocke za koju zelite izracunati elek. polje:")')
			write(*,'("Xt = ")',advance='no')
			read(*,*) xt
			write(*,'("Yt = ")',advance='no')
			read(*,*) yt
			write(*,'("Zt = ")',advance='no')
			read(*,*) zt
			write(*,'("Sloj u kojem je tocka: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
			write(*,'("----------------------------------------------")')
			write(*,'("(1) Raspodjelu Ex komponente elektricnog polja")')
			write(*,'("(2) Raspodjelu Ey komponente elektricnog polja")')
			write(*,'("(3) Raspodjelu Ez komponente elektricnog polja")')
			write(*,'("(4) Raspodjelu ukupnog/totalnog  elektr. polja")')
			write(*,'("----------------------------------------------")')
			write(*,'("Izbor: ")',advance='no')
			read(*,*) izbor

			write(*,'(/,"Racunam ...")')

			! Pocetak drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt3 = 1

			! ==========================================================================
			!    PRORACUN TRANZIJENTNOG ELEKTRICNOG POLJA U PRETHODNO ZADANOJ TOCKI
			!          U VISESLOJNOM SREDSTVU PO SVIM FREKVENCIJAMA OD FFT
			! ==========================================================================
			allocate(trans_pot(Nsample/2+1)) !Recikliranje varijable
			DO kk = 1, Nsample/2+1
				! Frekvencija za koju se racuna
				! -----------------------------
				f = freq(kk)
				! Narinuta tranzijentna struja (od FFT proracuna)
				! -----------------------------------------------
				Ii(:) = ccc(kk)
				! Proracun tranzijentnog elektricnog polja u jednoj tocki viseslojnog 
				! sredstva za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati
				! onoliko puta koliko se razlicitih frekvencija zeli racunati.
				! .......................................................................
				call single_frequency_elpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,&
				HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
				! .......................................................................
				trans_pot(kk) = potenc_freq		!Recikliranje varijable
			END DO

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Priprema vektora: ccc za proracun IFFT. Prva polovica + jedan clan vektora
			! ccc su vrijednosti tranzijentnog elektr. polja u odabranoj tocki za sve pozi-
			! tivne frekvencije + Nyquist-ovu frekvenciju. Druga polovica vektora: ccc je
			! konjugirano-kompleksna i simetricna s obzirom na Nyquist-ovu frekvenciju.
			! Formiran vektor ccc se transformira u vremensku domenu pomocu IFFT algoritma.
			N2 = Nsample/2+1
			! Istosmjerni clan
			ccc(1) = trans_pot(1)
			! Pozitivne frekvencije
			do i = 2,Nsample/2 
				ccc(i) = trans_pot(i)
			end do
			! Nyquist-ova frekvencija
			ccc(N2) = trans_pot(N2)
			! Negativne frekvencije
			do i = 2,Nsample/2
				IN = Nsample - i + 2
				ccc(IN) = dconjg(ccc(i))
			end do
			deallocate(trans_pot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			! Poziv subroutine za Inverse FFT
			! -------------------------------
			call zfftb(Nsample,ccc,wsave)
			! -------------------------------
			deallocate(wsave)

			!Tiskanje signala (nakon Inverse FFT)
			open(3,file="signal_elpolje.adf")
			write(3,'(/,5x,"t",10x,"Re[f(t)]",10x,"Im[f(t)]")')
			do i=1,Nsample
				write(3,'(f16.10,2x,e16.10,2x,e16.10)') ti(i),dreal(ccc(i)),dimag(ccc(i))
			end do
			close(3)
			deallocate(ccc)

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

			! Kraj drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time33 = end1 - start1

			write(*,'(/,"Napomena:")')
			write(*,'("Tranzijentno elektricno polje u odabranoj tocki visesloja snimljeno")')
			write(*,'("je u file: signal_elpolje.adf. Za prikaz koristiti: PlotXY.")')
			
		ELSE IF (izbor==4) THEN
			! ---------------------------------------------------------------------------
			!       PRORACUN TRANZIJENTNOG MAGNETSKOG POLJA U JEDNOJ TOCKI 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Unesite koord. tocke za koju zelite izracunati mag. polje:")')
			write(*,'("Xt = ")',advance='no')
			read(*,*) xt
			write(*,'("Yt = ")',advance='no')
			read(*,*) yt
			write(*,'("Zt = ")',advance='no')
			read(*,*) zt
			write(*,'("Sloj u kojem je tocka: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
			write(*,'("-----------------------------------------------")')
			write(*,'("(1) Raspodjelu Hx komponente  magnetskog  polja")')
			write(*,'("(2) Raspodjelu Hy komponente  magnetskog  polja")')
			write(*,'("(3) Raspodjelu Hz komponente  magnetskog  polja")')
			write(*,'("(4) Raspodjelu ukupnog/total.  magnetskog polja")')
			write(*,'("(5) Raspodjelu Bx komponente  magnet. indukcije")')
			write(*,'("(6) Raspodjelu By komponente  magnet. indukcije")')
			write(*,'("(7) Raspodjelu Bz komponente  magnet. indukcije")')
			write(*,'("(8) Raspodjelu ukupne/totalne magnet. indukcije")')
			write(*,'("-----------------------------------------------")')
			write(*,'("Izbor: ")',advance='no')
			read(*,*) izbor

			write(*,'(/,"Racunam ...")')

			! Pocetak drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt4 = 1

			! ==========================================================================
			!    PRORACUN TRANZIJENTNOG MAGNETSKOG POLJA U PRETHODNO ZADANOJ TOCKI
			!          U VISESLOJNOM SREDSTVU PO SVIM FREKVENCIJAMA OD FFT
			! ==========================================================================
			allocate(trans_pot(Nsample/2+1)) !Recikliranje varijable
			DO kk = 1, Nsample/2+1
				! Frekvencija za koju se racuna
				! -----------------------------
				f = freq(kk)
				! Narinuta tranzijentna struja (od FFT proracuna)
				! -----------------------------------------------
				Ii(:) = ccc(kk)
				! Proracun tranzijentnog magnetskog polja u jednoj tocki viseslojnog 
				! sredstva za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati
				! onoliko puta koliko se razlicitih frekvencija zeli racunati.
				! .......................................................................
				call single_frequency_magpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,&
				HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
				! .......................................................................
				trans_pot(kk) = potenc_freq		!Recikliranje varijable
			END DO

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Priprema vektora: ccc za proracun IFFT. Prva polovica + jedan clan vektora
			! ccc su vrijednosti tranzijentnog magnet. polja u odabranoj tocki za sve pozi-
			! tivne frekvencije + Nyquist-ovu frekvenciju. Druga polovica vektora: ccc je
			! konjugirano-kompleksna i simetricna s obzirom na Nyquist-ovu frekvenciju.
			! Formiran vektor ccc se transformira u vremensku domenu pomocu IFFT algoritma.
			N2 = Nsample/2+1
			! Istosmjerni clan
			ccc(1) = trans_pot(1)
			! Pozitivne frekvencije
			do i = 2,Nsample/2 
				ccc(i) = trans_pot(i)
			end do
			! Nyquist-ova frekvencija
			ccc(N2) = trans_pot(N2)
			! Negativne frekvencije
			do i = 2,Nsample/2
				IN = Nsample - i + 2
				ccc(IN) = dconjg(ccc(i))
			end do
			deallocate(trans_pot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			! Poziv subroutine za Inverse FFT
			! -------------------------------
			call zfftb(Nsample,ccc,wsave)
			! -------------------------------
			deallocate(wsave)

			!Tiskanje signala (nakon Inverse FFT)
			open(3,file="signal_magpolje.adf")
			write(3,'(/,5x,"t",10x,"Re[f(t)]",10x,"Im[f(t)]")')
			do i=1,Nsample
				write(3,'(f16.10,2x,e16.10,2x,e16.10)') ti(i),dreal(ccc(i)),dimag(ccc(i))
			end do
			close(3)
			deallocate(ccc)

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

			! Kraj drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time44 = end1 - start1

			write(*,'(/,"Napomena:")')
			write(*,'("Tranzijentno magnetsko polje u odabranoj tocki visesloja snimljeno")')
			write(*,'("je u file: signal_magpolje.adf. Za prikaz koristiti: PlotXY.")')
			
		ELSE IF(izbor==5) THEN
			! ---------------------------------------------------------------------------
			!  PRORACUN TRANZIJENTNOG MAGNETSKOG POLJA I MAG. INDUKCIJE DUZ PRAVCA (3D) 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Koordinate pocetne tocke pravca:")')
			write(*,'("Xp = ")',advance='no')
			read(*,*) xp1
			write(*,'("Yp = ")',advance='no')
			read(*,*) yp1
			write(*,'("Zp = ")',advance='no')
			read(*,*) zp1
			write(*,'("Koordinate krajnje tocke pravca:")')
			write(*,'("Xk = ")',advance='no')
			read(*,*) xk1
			write(*,'("Yk = ")',advance='no')
			read(*,*) yk1
			write(*,'("Zk = ")',advance='no')
			read(*,*) zk1
			write(*,'("Korak duz osi x: ")',advance='no')
			read(*,*) dx
			write(*,'("Korak duz osi y: ")',advance='no')
			read(*,*) dy
			write(*,'("Korak duz osi z: ")',advance='no')
			read(*,*) dz
			write(*,'("Sloj u kojem je pravac: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
			write(*,'("-----------------------------------------------")')
			write(*,'("(1) Raspodjelu Hx komponente  magnetskog  polja")')
			write(*,'("(2) Raspodjelu Hy komponente  magnetskog  polja")')
			write(*,'("(3) Raspodjelu Hz komponente  magnetskog  polja")')
			write(*,'("(4) Raspodjelu ukupnog/total.  magnetskog polja")')
			write(*,'("(5) Raspodjelu Bx komponente  magnet. indukcije")')
			write(*,'("(6) Raspodjelu By komponente  magnet. indukcije")')
			write(*,'("(7) Raspodjelu Bz komponente  magnet. indukcije")')
			write(*,'("(8) Raspodjelu ukupne/totalne magnet. indukcije")')
			write(*,'("-----------------------------------------------")')
			write(*,'("Izbor: ")',advance='no')
			read(*,*) izbor

			write(*,'(/,"Racunam ...")')
				
			! Pocetak mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt5 = 1
				
			if (dx/=0.0) then
				Ntoc = (xk1-xp1)/dx
			else if(dy/=0.0) then
				Ntoc = (yk1-yp1)/dy
			else if(dz/=0.0) then
				Ntoc = (zk1-zp1)/dz
			end if

			allocate(TransPot(Nsample/2+1,Ntoc)) !Recikliranje varijable
			allocate(trans_pot(Nsample/2+1))	 !Recikliranje varijable
			xt = xp1
			yt = yp1
			zt = zp1
			do i = 1,Ntoc
				write(*,'(" -- procesiram vrijednost: ",i3," / ",i3)') i,Ntoc
				! ==========================================================================
				!    PRORACUN TRANZIJENTNOG MAGNETSKOG POLJA U PRETHODNO ZADANOJ TOCKI
				!          U VISESLOJNOM SREDSTVU PO SVIM FREKVENCIJAMA OD FFT
				! ==========================================================================
				DO kk = 1, Nsample/2+1
					! Frekvencija za koju se racuna
					! -----------------------------
					f = freq(kk)
					! Narinuta tranzijentna struja (od FFT proracuna)
					! -----------------------------------------------
					Ii(:) = ccc(kk)
					! Proracun tranzijentnog magnetskog polja u jednoj tocki viseslojnog 
					! sredstva za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati
					! onoliko puta koliko se razlicitih frekvencija zeli racunati.
					! .......................................................................
					call single_frequency_magpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,&
					HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
					! .......................................................................
					trans_pot(kk) = potenc_freq		!Recikliranje varijable
				END DO
				TransPot(:,i) = trans_pot(:)
					
				xt = xt + dx
				yt = yt + dy
				zt = zt + dz
			end do
			deallocate(trans_pot)

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Prva polovica vektora + jedan clan vektora su vrijednosti tranzijentnog 
			! magnet. polja u odabranoj tocki za sve pozitivne frekvencije + Nyquist-ovu
			! frekvenciju. Druga polovica vektora je konjugirano-kompleksna i simetricna
			! s obzirom na Nyquist-ovu frekvenciju. Formiran vektor se transformira u 
			! vremensku domenu pomocu IFFT algoritma (FFTPACK).
			N2 = Nsample/2+1
			allocate(TransPotP(Nsample,Ntoc))	!Recikliranje varijable
			DO kk = 1,Ntoc
				! Istosmjerni clan
				TransPotP(1,kk) = TransPot(1,kk)
				! Pozitivne frekvencije
				do i = 2,Nsample/2 
					TransPotP(i,kk) = TransPot(i,kk)
				end do
				! Nyquist-ova frekvencija
				TransPotP(N2,kk) = TransPot(N2,kk)
				! Negativne frekvencije
				do i = 2,Nsample/2
					IN = Nsample - i + 2
					TransPotP(IN,kk) = dconjg(TransPot(i,kk))
				end do
			END DO
			deallocate(TransPot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			allocate(FF(Ntoc,Nsample)) !Recikliranje varijable FF
			DO kk = 1,Ntoc
				ccc(:) = TransPotP(:,kk)
				! Poziv subroutine za Inverse FFT
				!-------------------------------
				call zfftb(Nsample,ccc,wsave)
				!-------------------------------
				FF(kk,:) = dreal(ccc(:))
			END DO
			deallocate(ccc,wsave)
			deallocate(TransPotP)

			! Kraj mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time55 = end1 - start1

			!-----------------------------------------------------------------------
			!       3D prikaz izracunatog magnetskog polja - Gnuplot
			!-----------------------------------------------------------------------
			open(7,file='3Dmagpolje.dat')
			do j = 1,Ntoc
				do i = 1,Nsample
					write(7,'(e12.4)') FF(j,i)
				end do
				write(7,*)
			end do
			close(7)
			!-----------------------------------------------------------------------
			!       3D prikaz izracunatog magnetskog polja - LabPlot
			!-----------------------------------------------------------------------
			open(7,file='3Dtrans_magpolje.dat')
			do i = 1,Nsample
				write(7,'(80e12.4)') (FF(j,i),j=1,Ntoc)
			end do
			close(7)
			deallocate(FF)
					
			write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele tranzijentnog magnet. polja/indukcije snimljeni")')
			write(*,'("su u file: 3Dmagpolje.dat. Za prikaz koristiti: Gnuplot. Za LabPlot: 3Dtrans_magpolje.dat")')

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

		ELSE IF(izbor==6) THEN
			! ---------------------------------------------------------------------------
			!        PRORACUN TRANZIJENTNOG ELEKTRICNOG POLJA DUZ PRAVCA (3D) 
			! ---------------------------------------------------------------------------
			write(*,'(/,"Koordinate pocetne tocke pravca:")')
			write(*,'("Xp = ")',advance='no')
			read(*,*) xp1
			write(*,'("Yp = ")',advance='no')
			read(*,*) yp1
			write(*,'("Zp = ")',advance='no')
			read(*,*) zp1
			write(*,'("Koordinate krajnje tocke pravca:")')
			write(*,'("Xk = ")',advance='no')
			read(*,*) xk1
			write(*,'("Yk = ")',advance='no')
			read(*,*) yk1
			write(*,'("Zk = ")',advance='no')
			read(*,*) zk1
			write(*,'("Korak duz osi x: ")',advance='no')
			read(*,*) dx
			write(*,'("Korak duz osi y: ")',advance='no')
			read(*,*) dy
			write(*,'("Korak duz osi z: ")',advance='no')
			read(*,*) dz
			write(*,'("Sloj u kojem je pravac: ")',advance='no')
			read(*,*) sloj

			write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
			write(*,'("----------------------------------------------")')
			write(*,'("(1) Raspodjelu Ex komponente elektricnog polja")')
			write(*,'("(2) Raspodjelu Ey komponente elektricnog polja")')
			write(*,'("(3) Raspodjelu Ez komponente elektricnog polja")')
			write(*,'("(4) Raspodjelu ukupnog/totalnog  elektr. polja")')
			write(*,'("----------------------------------------------")')
			write(*,'("Izbor: ")',advance='no')
			read(*,*) izbor

			write(*,'(/,"Racunam ...")')
				
			! Pocetak mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt6 = 1
				
			if (dx/=0.0) then
				Ntoc = (xk1-xp1)/dx
			else if(dy/=0.0) then
				Ntoc = (yk1-yp1)/dy
			else if(dz/=0.0) then
				Ntoc = (zk1-zp1)/dz
			end if

			allocate(TransPot(Nsample/2+1,Ntoc)) !Recikliranje varijable
			allocate(trans_pot(Nsample/2+1))	 !Recikliranje varijable
			xt = xp1
			yt = yp1
			zt = zp1
			do i = 1,Ntoc
				write(*,'(" -- procesiram vrijednost: ",i3," / ",i3)') i,Ntoc
				! ==========================================================================
				!    PRORACUN TRANZIJENTNOG ELEKTRICNOG POLJA U PRETHODNO ZADANOJ TOCKI
				!          U VISESLOJNOM SREDSTVU PO SVIM FREKVENCIJAMA OD FFT
				! ==========================================================================
				DO kk = 1, Nsample/2+1
					! Frekvencija za koju se racuna
					! -----------------------------
					f = freq(kk)
					! Narinuta tranzijentna struja (od FFT proracuna)
					! -----------------------------------------------
					Ii(:) = ccc(kk)
					! Proracun tranzijentnog elektricnog polja u jednoj tocki viseslojnog 
					! sredstva za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati
					! onoliko puta koliko se razlicitih frekvencija zeli racunati.
					! .......................................................................
					call single_frequency_elpolje(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,&
					HD,h,iso,sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,izbor,potenc_freq)
					! .......................................................................
					trans_pot(kk) = potenc_freq		!Recikliranje varijable
				END DO
				TransPot(:,i) = trans_pot(:)
					
				xt = xt + dx
				yt = yt + dy
				zt = zt + dz
			end do
			deallocate(trans_pot)

			! PRIPREMA ZA POZIV IFFT PRORACUNA
			! --------------------------------
			! Prva polovica vektora + jedan clan vektora su vrijednosti tranzijentnog 
			! elektr. polja u odabranoj tocki za sve pozitivne frekvencije + Nyquist-ovu
			! frekvenciju. Druga polovica vektora je konjugirano-kompleksna i simetricna
			! s obzirom na Nyquist-ovu frekvenciju. Formiran vektor se transformira u 
			! vremensku domenu pomocu IFFT algoritma (FFTPACK).
			N2 = Nsample/2+1
			allocate(TransPotP(Nsample,Ntoc))	!Recikliranje varijable
			DO kk = 1,Ntoc
				! Istosmjerni clan
				TransPotP(1,kk) = TransPot(1,kk)
				! Pozitivne frekvencije
				do i = 2,Nsample/2 
					TransPotP(i,kk) = TransPot(i,kk)
				end do
				! Nyquist-ova frekvencija
				TransPotP(N2,kk) = TransPot(N2,kk)
				! Negativne frekvencije
				do i = 2,Nsample/2
					IN = Nsample - i + 2
					TransPotP(IN,kk) = dconjg(TransPot(i,kk))
				end do
			END DO
			deallocate(TransPot)

			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			allocate(FF(Ntoc,Nsample)) !Recikliranje varijable FF
			DO kk = 1,Ntoc
				ccc(:) = TransPotP(:,kk)
				! Poziv subroutine za Inverse FFT
				!-------------------------------
				call zfftb(Nsample,ccc,wsave)
				!-------------------------------
				FF(kk,:) = dreal(ccc(:))
			END DO
			deallocate(ccc,wsave)
			deallocate(TransPotP)

			! Kraj mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time66 = end1 - start1

			!-----------------------------------------------------------------------
			!       3D prikaz izracunatog elektricnog polja - Gnuplot
			!-----------------------------------------------------------------------
			open(7,file='3Delpolje.dat')
			do j = 1,Ntoc
				do i = 1,Nsample
					write(7,'(e12.4)') FF(j,i)
				end do
				write(7,*)
			end do
			close(7)	
			!-----------------------------------------------------------------------
			!       3D prikaz izracunatog elektricnog polja - LabPlot
			!-----------------------------------------------------------------------
			open(7,file='3Dtrans_elpolje.dat')
			do i = 1,Nsample
				write(7,'(80e12.4)') (FF(j,i),j=1,Ntoc)
			end do
			close(7)		
			deallocate(FF)
			
			write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele tranzijentnog elektr. polja snimljeni su u file:")')
			write(*,'("3Delpolje.dat. Za prikaz koristiti: Gnuplot. Za LabPlot: 3Dtrans_elpolje.dat")')

			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

		ELSE IF (izbor==7) THEN
			! ---------------------------------------------------------------------------
			!           PRORACUN TRANZIJENTNE ULAZNE IMPEDANCIJE UZEMLJIVACA 
			! ---------------------------------------------------------------------------
			! Tranzijentna ulazna impedancija uzemljivackog sustava definirana je kao
			! omjer izmedju tranzijentnog napona na metalu u tocki u kojoj je narinuta
			! tranzijentna struja i tranzijetne narinute struje - Z(t) = u(t)/i(t)
			
			!Koordinate tocke u kojoj je narinuta tranzijentna struja
			xt = xi; yt = yi; zt = zi
			!Sloj u kojem se nalazi tocka sa narinutom tranz. strujom
			sloj = 1
			do j = 1,N-1
				if (zt>HD(j)) sloj = sloj + 1
			end do
			
			write(*,'(/,"Racunam ...")')
			
			! Pocetak drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt7 = 1
			
			! ==========================================================================
			!    PRORACUN TRANZIJENTNOG POTENCIJALA U PRETHODNO ZADANOJ TOCKI (TOCKA
			!    U KOJOJ JE NARINUTA TRANZIJENTNA STRUJA) PO SVIM FREKVENCIJAMA OD FFT
			! ==========================================================================
			allocate(trans_pot(Nsample/2+1))
			DO kk = 1, Nsample/2+1
				! Frekvencija za koju se racuna
				! -----------------------------
				f = freq(kk)
				! Struja za koju se racuna - FFT
				! ------------------------------
				Ii(:) = ccc(kk)
				! Proracun tranzijentnog potencijala u jednoj tocki viseslojnog sredstva
				! za bilo koju frekvenciju - ovu subroutinu potrebno je pozvati onoliko
				! puta koliko se razlicitih frekvencija zeli racunati.
				! .......................................................................
				call single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,h,iso,&
				sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
				! .......................................................................
				trans_pot(kk) = potenc_freq		
			END DO
	
			! ---------------------------------------------------------------------------
			!                   PRIPREMA ZA POZIV INVERSE FFT
			! ---------------------------------------------------------------------------
			N2 = Nsample/2+1
			! Istosmjerni clan
			ccc(1) = trans_pot(1)
			! Pozitivne frekvencije
			do i = 2,Nsample/2 
				ccc(i) = trans_pot(i)
			end do
			! Nyquist-ova frekvencija
			ccc(N2) = trans_pot(N2)
			! Negativne frekvencije
			do i = 2,Nsample/2
				IN = Nsample - i + 2
				ccc(IN) = dconjg(ccc(i))
			end do
			deallocate(trans_pot)
	
			! ---------------------------------------------------------------------------
			!            PRORACUN INVERZNE FFT TRANSFORMACIJE - FFTPACK
			! ---------------------------------------------------------------------------
			! Poziv subroutine za Inverse FFT
			! -------------------------------
			call zfftb(Nsample,ccc,wsave)
			! -------------------------------
			deallocate(wsave)

			! ---------------------------------------------------------------------------
			!           TRANZIJENTNA ULAZNA IMPEDANCIJA UZEMLJIVACA 
			! ---------------------------------------------------------------------------
			allocate(wsave(Nsample-1))  !Recikliranje varijable: wsave
			do i = 2,Nsample
				wsave(i-1) = dreal(ccc(i))/func(i)
			end do
			deallocate(ccc)

			!Tiskanje ulazne impedancije
			open(3,file="trans_imped.adf")
			write(3,'(/,5x,"t",10x,"Imped(t)")')
			do i=1,Nsample-1
				write(3,'(f16.10,2x,e16.10)') ti(i),wsave(i)
			end do
			close(3)
			deallocate(wsave)
	
			! Oslobadjanje zauzete memorije
			deallocate(ti,func)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)
			
			! Kraj drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time77 = end1 - start1
	
			write(*,'(/,"Napomena:")')
			write(*,'("Tranzijentna ulazna impedanicja uzemljivackog sustava snimljena")')
			write(*,'("je u file: trans_imped.adf. Za prikaz koristiti: PlotXY.")')

		ELSE IF (izbor==8) THEN
			! ---------------------------------------------------------------------------
			!           PRORACUN Y(s) ZA POTREBE VECTOR FITTING ANALIZE
			! ---------------------------------------------------------------------------
			!Koordinate tocke u kojoj je narinuta tranzijentna struja
			xt = xi; yt = yi; zt = zi
			!Sloj u kojem se nalazi tocka sa narinutom tranz. strujom
			sloj = 1
			do j = 1,N-1
				if (zt>HD(j)) sloj = sloj + 1
			end do

			write(*,'(/,"Racunam ...")')

			! Pocetak drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(start1)
			tt8 = 1

			! ==========================================================================
			!    PRORACUN Y(s) U PRETHODNO ZADANOJ TOCKI (TOCKA
			!    U KOJOJ JE NARINUTA STRUJA) PO SVIM FREKVENCIJAMA
			! ==========================================================================
			allocate(trans_pot(Nsample+1))  !recikliranje varijable
			allocate(jomega(Nsample+1))
			DO kk = 1, Nsample+1
				! Frekvencija za koju se racuna
				! -----------------------------
				f = freq(kk)
				jomega(kk) = 2.d0*PI*freq(kk) * dcmplx(0.d0,1.d0)
				! Struja za koju se racuna - FFT
				! ------------------------------
				!Ii(:) = ccc(kk)
				Ii(:) = cmplx(1.d0,0.d0)
				! Proracun Y(s)za bilo koju frekvenciju - ovu subroutinu potrebno je
				! pozvati onoliko puta koliko se razlicitih frekvencija zeli racunati.
				! .......................................................................
				call Y_single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,h,iso,&
				sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
				! .......................................................................
				trans_pot(kk) = potenc_freq  !recikliranje varijable
				if (kk == nint((Nsample+1)/10.)) then
					write(*,'(/,"...10%")')
				else if (kk == nint(2*(Nsample+1)/10.)) then
					write(*,'("...20%")')
				else if (kk == nint(3*(Nsample+1)/10.)) then
					write(*,'("...30%")')
				else if (kk == nint(4*(Nsample+1)/10.)) then
					write(*,'("...40%")')
				else if (kk == nint(5*(Nsample+1)/10.)) then
					write(*,'("...50%")')
				else if (kk == nint(6*(Nsample+1)/10.)) then
					write(*,'("...60%")')
				else if (kk == nint(7*(Nsample+1)/10.)) then
					write(*,'("...70%")')
				else if (kk == nint(8*(Nsample+1)/10.)) then
					write(*,'("...80%")')
				else if (kk == nint(9*(Nsample+1)/10.)) then
					write(*,'("...90%")')
				else if (kk == Nsample) then
					write(*,'("..100%")')
				end if
			END DO

			!Tiskanje ulazne impedancije
			open(3,file="admintancija.txt")
			!write(3,'(10x,"jomega",10x,"ReY(s)",10x,"ImY(s)")')
			do i=1,Nsample+1
				write(3,'(e16.6,3x,e16.6,3x,e16.6)') aimag(jomega(i)),real(trans_pot(i)),aimag(trans_pot(i))
			end do
			close(3)

			deallocate(trans_pot)
			deallocate(jomega)

			! Oslobadjanje zauzete memorije
			!deallocate(ti,func)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)

			! Kraj drugog mjerenja CPU vremena trajanja proracuna
			call cpu_time(end1)
			time77 = end1 - start1

			write(*,'(/,"Napomena:")')
			write(*,'("Rezultati Y(s) u fileu: admintancija.txt")')

		ELSE IF(izbor==0) THEN
			! Oslobadjanje zauzete memorije
			deallocate(ti,func,freq)
			deallocate(ccc,wsave)
			deallocate(iveze)
			deallocate(L,r0)
			deallocate(xp,yp,zp)
			deallocate(xk,yk,zk)
			deallocate(HD,h)
			deallocate(iso)
			deallocate(ro,epsr)
			deallocate(mirv,sigmav)
		END IF

	ELSE IF(tip_pr==1) THEN
		! ***************************************************************************
		!        ANALIZA UZEMLJIVACKIH SUSTAVA U VISESLOJU ZA 1 FREKVENCIJU
		!              (STACIONARNA ANALIZA UZEMLJIVACKIH SUSTAVA)
		! ***************************************************************************
		write(*,'(/,"Stacionarna analiza uzemljivaca...")')
		! Pocetak prvog mjerenja CPU vremena trajanja proracuna
		call cpu_time(start1)

		! ---------------------------------------------------------------------------
		!               PRORACUN KONSTANTNIH PARAMETARA PROGRAMA
		! ---------------------------------------------------------------------------
		allocate(HD(N))
		! Proracun z koordinate donjih granicnih ploha slojeva modela
		call vektor_H(N,h,HD)
		allocate(kapa(N))
		! Proracun kompleksnih specificnih vodljivosti slojeva modela
		call proracun_kapa(N,f,ro,epsr,kapa)
		! Oslobadjanje dijela zauzete memorije
		deallocate(ro,epsr)
		allocate(gama(N))
		! Proracun valnih konstanti slojeva modela
		call proracun_gama(f,N,kapa,gama)
		allocate(Fr(N))
		! Proracun faktora refleksije svih slojeva modela
		call vektor_F(N,kapa,Fr)
		allocate(L(Ns))
		! Proracun duljine svih segmenata
		do i = 1,Ns
			L(i) = dsqrt((xk(i)-xp(i))**2+(yk(i)-yp(i))**2+(zk(i)-zp(i))**2)
		end do
		allocate(zunjv(Ns))
		! Proracun jedinicnih unutarnjih impedancija svih segmenata
		call proracun_zunjv(Ns,f,mirv,sigmav,r0,zunjv)
		deallocate(mirv)
		allocate(iso(Ns))
		! Odredjivanje sloja u kojem se nalazi svaki od segmenata
		do i = 1,Ns
			zsi = (zp(i)+zk(i))/2.d0
			iso(i) = 1
			do j = 1,N-1
				if (zsi>HD(j)) iso(i) = iso(i) + 1
			end do
		end do

		write(*,'(/,"Racunam ...")')

		! ---------------------------------------------------------------------------
		!             FORMIRANJE UKUPNE LOKALNE MATRICE SUSTAVA [Y] 
		! ---------------------------------------------------------------------------
		allocate(Yuv(Ns,Ns))
		allocate(Ypv(Nss,Nss))
		allocate(Yu(2*Ns,2*Ns))
		allocate(Ypp(2*Nss,2*Nss))
		! ................................................................
		call matrica(f,L,xp,yp,zp,xk,yk,zk,Ns,Nsi,r0,kapa,gama,sigmav,zunjv,&
			 iso,N,Fr,HD,h,Yuv,Ypv,Yu,Ypp)
		! ................................................................
		deallocate(zunjv,sigmav)

		! ---------------------------------------------------------------------------
		!               FORMIRANJE VEKTORA VEZE ZA ASEMBLIRANJE
		! ---------------------------------------------------------------------------
		allocate(iveze(2*Ns))
		allocate(Lk(3,2*Ns))
		Lk(1,1:Ns) = xp(:)
		Lk(1,(Ns+1):(2*Ns)) = xk(:)
		Lk(2,1:Ns) = yp(:)
		Lk(2,(Ns+1):(2*Ns)) = yk(:)
		Lk(3,1:Ns) = zp(:)
		Lk(3,(Ns+1):(2*Ns)) = zk(:)
		ngc = 1
		allocate(Gx(ngc))
		allocate(Gy(ngc))
		allocate(Gz(ngc))
		eps = minval(r0)	! Tolerancija
		Gx(1) = Lk(1,1)
		Gy(1) = Lk(2,1)
		Gz(1) = Lk(3,1)
		iveze(1) = 1
		do i = 2,2*Ns
			xL = Lk(1,i)
			yL = Lk(2,i)
			zL = Lk(3,i)
			pom = 0
			do j = 1,ngc
				if ((dabs(xL-Gx(j))<=eps).and.(dabs(yL-Gy(j))<=eps)&
					.and.(dabs(zL-Gz(j))<=eps)) then
					pom = j
				end if
			end do
			if (pom/=0) then
				iveze(i) = pom
			else
				ngc = ngc + 1
				! Varijabla x
				allocate(Gt(ngc-1))
				Gt = Gx
				deallocate(Gx)
				allocate(Gx(ngc))
				Gx(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gx(ngc) = Lk(1,i)
				! Varijabla y
				allocate(Gt(ngc-1))
				Gt = Gy
				deallocate(Gy)
				allocate(Gy(ngc))
				Gy(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gy(ngc) = Lk(2,i)
				! Varijabla z
				allocate(Gt(ngc-1))
				Gt = Gz
				deallocate(Gz)
				allocate(Gz(ngc))
				Gz(1:(ngc-1)) = Gt
				deallocate(Gt)
				Gz(ngc) = Lk(3,i)
				iveze(i) = ngc
			end if
		end do
		deallocate(Lk)
		! ---------------------------------------------------------------------------
		!          Ukupni broj globalnih cvorova 
		! ---------------------------------------------------------------------------
		Nc = ngc
		! ---------------------------------------------------------------------------
		! Odredjivanje vrijednosti globalnog cvora u kojem je injektirana struja
		! ---------------------------------------------------------------------------
		do j = 1,Ni
			do i = 1,Nc
				if ((dabs(Gx(i)-xii(j))<=eps).and.(dabs(Gy(i)-yii(j))<=eps)&
					.and.(dabs(Gz(i)-zii(j))<=eps)) then
					icv(j) = i
				end if
			end do
		end do
		deallocate(xii,yii,zii)

		! Tiskanje koordinata globalnih cvorova u privremeni file
		open(3,file='temp.dat')
		do i = 1,Nc
			write(3,'(3f12.3)') Gx(i),Gy(i),Gz(i)
		end do
		close(3)
		deallocate(Gx,Gy,Gz)
		! Tiskanje ulaznih podataka u vanjski file
		! -----------------------------------------
		write(2,'(/,"Podaci o narinutim strujama:")')
		write(2,'(28("-"))')
		write(2,'("Broj globalnih cvorova u kojima su narinute struje: ",i2)') Ni
		write(2,'("Frekvencija narinutih struja:",f10.2," [Hz]")') f
		write(2,*)
		write(2,'("Cvor.",T10,"Re(I) [A]",T25,"Im(I) [A]")')
		write(2,'(35("-"))')
		do i = 1,Ni
			write(2,'(T3,i3,T10,f10.3,T20,f10.3)') icv(i),rei,imi
		end do
		write(2,'(35("-"))')
		write(2,'(//,T5,"R E Z U L T A T I   P R O R A C U N A")')
		
		! ---------------------------------------------------------------------------
		!         ASSEMBLING - FORMIRANJE UKUPNE GLOBALNE MATRICE SUSTAVA [Yg] 
		! ---------------------------------------------------------------------------
		allocate(Yg(Nc,Nc))
		Yg(:,:) = dcmplx(0.d0,0.d0)
		! Asembliranje uzduzne matrice admitancija [Yu]
		do j = 1,2*Ns
			jg = iveze(j)
			do i = 1,2*Ns
				ig = iveze(i)
				Yg(ig,jg) = Yg(ig,jg) + Yu(i,j)
			end do
		end do
		deallocate(Yu)
		! Asembliranje poprecne matrice admitancija [Yp]
		do j = 1,2*Nss
			jg = iveze(j)
			do i = 1,2*Nss
				ig = iveze(i)
				Yg(ig,jg) = Yg(ig,jg) + Ypp(i,j)
			end do
		end do
		deallocate(Ypp)
		
		! ---------------------------------------------------------------------------
		!        RJESAVANJE GLOBALNOG SUSTAVA JEDNADZBI [Yg]*{Fig} = {Ig} 
		!	     ZA ODREDJIVANJE NEPOZNATIH POTENCIJALA GLOBALNIH CVOROVA
		! ---------------------------------------------------------------------------
		! Globalni sustav jednadzbi, complex(8), rjesava se pomocu LAPACK rutine:
		! ZGESV. Ova subroutina koristi punu LU dekompoziciju matrice s parcijalnim
		! pivotiranjem (row interchanges), te Gaussovu elimincaiju.
		nrhs = 1
		allocate(B(Nc,nrhs))
		B(:,nrhs) = dcmplx(0.d0,0.d0)
		do i = 1,Nc
    		do j = 1,Ni
    			if (icv(j)==i) then
    				B(i,nrhs) = Ii(j)
    			end if
    		end do
		end do 
		allocate(ipiv(Nc))
		! Poziv Lapack subroutine: ZGESV
		! ...........................................
		call zgesv(Nc,nrhs,Yg,Nc,ipiv,B,Nc,info)
		! ...........................................
		deallocate(ipiv)
		allocate(Fig(Nc))
		Fig(:) = B(:,nrhs)
		deallocate(B)
		deallocate(Yg)

		! Tiskanje potencijala globalnih cvorova u izlaznu datoteku
		!---------------------------------------------------------------------------
		write(2,'(/,t10,"Potencijal globalnih cvorova uzemljivaca:")')
		write(2,'(62("-"))')
		write(2,'("Red.",t10,"x",t20,"y",t30,"z",t45,"Potencijal")')
		write(2,'("broj",t9,"(m)",t19,"(m)",t29,"(m)",t40,"Iznos (V)",t55,"kut (o)")')
		write(2,'(62("-"))')
		open(3,file='temp.dat')
		open(5,file='glob_pot.dat') !Za proracun napona dodira
		do i = 1,Nc
			read(3,'(3f12.3)') xl,yl,zl
			open(4,file='tempp.dat')
			do j = 1,Nv
				read(4,'(6f12.3)') x1,y1,z1,x2,y2,z2
				if (((dabs(xl-x1)<=eps).and.(dabs(yl-y1)<=eps).and.(dabs(zl-z1)<=eps)).or.&
					((dabs(xl-x2)<=eps).and.(dabs(yl-y2)<=eps).and.(dabs(zl-z2)<=eps))) then
					! Tiskaj samo globalne cvorove vodica (ne i segmenata)
					modul = dsqrt(dreal(Fig(i))**2+dimag(Fig(i))**2)
					kut = datan2d(dreal(Fig(i)),dimag(Fig(i)))
					write(2,'(i3,t5,f8.2,t15,f8.2,t25,f8.2,t35,f16.4,t53,f8.2)') i,xl,yl,zl,modul,kut
					! Spremanje potencijala globalnih cvorova za proracun napona dodira
					write(5,'(i3,3f8.2,2e16.6)') i,xl,yl,zl,dreal(Fig(i)),dimag(Fig(i))
				end if
			end do
			close(4)
		end do
		write(2,'(62("-"))')
		!res = DELFILESQQ('tempp.dat')
		close(3)
		!res = DELFILESQQ('temp.dat')
		close(5)

		! ---------------------------------------------------------------------------
		!                  PRORACUN ULAZNE IMPEDANCIJE UZEMLJIVACA
		! ---------------------------------------------------------------------------
		! Ulazna impedancija racuna se samo u slucaju da je narinuta struja samo
		! u jednom cvoru (kao omjer napona u tom cvoru prema narinutoj struji)
		if (Ni==1) then
			j = icv(Ni)
			Zuzem = Fig(j)/Ii(Ni)
			deallocate(icv)
			deallocate(Ii)
			call zmk(Zuzem,modul,kut)
			! Tiskanje rezultata u file
			write(2,'(/,5x,"Ulazna impedancija uzemljivaca:")')
			write(2,'(45("-"))')
			write(2,'(" Zu =",f10.4," /",f7.2," (o)","   [Ohm]")') modul,kut
			write(2,'(45("-"))')
		end if

		! ---------------------------------------------------------------------------
		!               PRORACUN POTENCIJALA LOKALNIH CVOROVA
		! ---------------------------------------------------------------------------
		allocate(Fi(2*Ns))
		do i = 1,2*Ns
			Fi(i) = Fig(iveze(i))
		end do
		deallocate(Fig)
		deallocate(iveze)

		! ---------------------------------------------------------------------------
		!               PRORACUN UZDUZNIH STRUJA SVIH SEGMENATA
		! ---------------------------------------------------------------------------
		! Uzduzne struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
		! {Iu} = [Yuv]*[Au]*{Fi}, pri cemu vektor {Fi} predstavlja vektor potencijala
		! lokalnih cvorova segmenata uzemljivaca. Ovaj proracun odvija se u dva koraka,
		! kako slijedi:
		!   1) Formira se matrica [Au] koja je pravokutna Au(Ns,2*Ns) i kompleksna
		!      (moze biti i realna).
		!   2) Racuna se izraz {Iu} = [Yuv]*[Au]*{Fi} koristeci BLAS3 i BLAS2 rutine:
		!      ZGEMM i ZGEMV, respektivno. Subroutina ZGEMM mnozi dvije generalne kom-
		!      pleksne matrice (korsteci blokovski algoritam), dok subroutina ZGEMV
		!      mnozi kompleksnu matricu s vektorom.
		! FORMIRANJE MATRICE [Au]
		allocate(Au(Ns,2*Ns))
		Au(:,:) = dcmplx(0.d0,0.d0)
		do i = 1,Ns
			do j = 1,Ns
				if (i==j) then      
					Au(i,j) = dcmplx(1.d0,0.d0)
				end if
			end do
			do j = Ns+1,2*Ns
				if ((i+Ns)==j) then
					Au(i,j) = dcmplx(-1.d0,0.d0)
				end if
			end do
		end do
		! RACUNANJE IZRAZA {Iu} = [Yuv]*[Au]*{Fi} 
		! ----------------------------------------
		alpha = dcmplx(1.d0,0.d0)
		betha = dcmplx(0.d0,0.d0)
		transa = 'N'
		transb = 'N'
		allocate(CC(Ns,2*Ns))
		call zgemm(transa,transb,Ns,2*Ns,Ns,alpha,Yuv,Ns,Au,Ns,betha,CC,Ns)
		deallocate(Au)
		deallocate(Yuv)
		allocate(Iu(Ns))
		call zgemv(transa,Ns,2*Ns,alpha,CC,Ns,Fi,1,betha,Iu,1)
		deallocate(CC)

		! ---------------------------------------------------------------------------
		!               PRORACUN POPRECNIH STRUJA SVIH SEGMENATA
		! ---------------------------------------------------------------------------
		! Poprecne struje svih segmenata racunaju se iz sljedece matricne jednadzbe:
		! {Ip} = [Ypv]*[Ap]*{Fi}, pri cemu vektor {Fi} predstavlja vektor potencijala
		! lokalnih cvorova segmenata uzemljivaca. Ovaj proracun odvija se u dva koraka,
		! kako slijedi:
		!   1) Formira se matrica [Ap] koja je pravokutna Ap(Ns,2*Ns) i kompleksna (moze
		!      biti i realna).
		!   2) Racuna se izraz {Ip} = [Ypv]*[Ap]*{Fi} koristeci BLAS3 rutine: ZGEMM i
		!      ZGEMV
		! FORMIRANJE MATRICE [Ap]
		allocate(Ap(Nss,2*Nss))
		Ap(:,:) = dcmplx(0.d0,0.d0)
		do i = 1,Nss
			do j = 1,Nss
				if (i==j) then
					Ap(i,j) = dcmplx(0.5d0,0.d0)
				end if
			end do
			do j = Nss+1,2*Nss
				if ((i+Nss)==j) then
					Ap(i,j) = dcmplx(0.5d0,0.d0)
				end if
			end do
		end do
		! RACUNANJE IZRAZA {Ip} = [Yp]*[Ap]*{Fi} 
		! ------------------------------------------------------
		alpha = dcmplx(1.d0,0.d0)
		betha = dcmplx(0.d0,0.d0)
		transa = 'N'
		transb = 'N'
		allocate(CC(Nss,2*Nss))
		call zgemm(transa,transb,Nss,2*Nss,Nss,alpha,Ypv,Nss,Ap,Nss,betha,CC,Nss)
		deallocate(Ap)
		deallocate(Ypv)
		allocate(Ip(Nss))
		call zgemv(transa,Nss,2*Nss,alpha,CC,Nss,Fi(1:2*Nss),1,betha,Ip,1)
		deallocate(CC)
		deallocate(Fi)

		! Tiskanje poprecnih i uzduznih struja svih segmenata u izlaznu datoteku
		! ----------------------------------------------------------------------
		write(2,'(/,10x,"Poprecne i uzduzne struje svih segmenata uzemljivaca:")')
		write(2,'(80("-"))')
		write(2,'(20x,"Poprecne struje",20x,"Uzduzne struje")')
		write(2,'(80("-"))')
		write(2,'("Br. seg.",10x,"Re + j*Im (A)",25x,"Re + j*Im (A)")')
		write(2,'(80("-"))')
		do i = 1,Nss
			write(2,'(i3,5x,2f16.6,5x,2f16.6)') i,Ip(i),Iu(i)
		end do
		do i = Nss+1,Ns
			write(2,'(i3,42x,2f16.6)') i,Iu(i)
		end do
		write(2,'(80("-"))')

		! Kraj prvog mjerenja CPU vremena trajanja proracuna
		call cpu_time(end1)
		time3 = end1 - start1
	
		! *********************************************************************************
		!                 PRORACUNI OD INTERESA KORISNIKU PROGRAMA
		! *********************************************************************************
		ponovi = .true.
		ti1 = 0; ti2 = 0; ti3 = 0; ti4 = 0
		ti5 = 0; ti6 = 0; ti7 = 0; ti8 = 0
		ti9 = 0; ti10 = 0; ti11 = 0; ti12 = 0; ti13 = 0

		DO WHILE(ponovi)
			write(*,'(/,"Odaberite sto zelite racunati:")')
			write(*,'("---------------------------------------------------------------")')
			write(*,'("Skalarni potencijal, naponi koraka i naponi dodira:")')
			write(*,'("---------------------------------------------------------------")')
			write(*,'(4x,"1 - Potencijal u samo jednoj tocki viseslojnog sredstva")')
			write(*,'(4x,"2 - Raspodjelu potencijala duz pravca u viseslojnom sredstvu")')
			write(*,'(4x,"3 - Raspodjelu potencijala duz nekog pravca na povrsini tla")')
			write(*,'(4x,"4 - Napon koraka duz nekog pravca")')
			write(*,'(4x,"5 - Raspodjelu potencijala duz vise paralelnih pravaca")')
			write(*,'("        na povrsini tla istovremeno (3D raspodjela potencijala)")')
			write(*,'(4x,"6 - Napon dodira u nekoj tocki")')
			write(*,'(4x,"7 - Napon koraka duz vise paralelnih pravaca istovremeno")')
			write(*,'("        (3D raspodjela napona koraka u smjeru osi x/y)")')
			write(*,'(3x,"14 - Napon dodira duz jednog pravca")')
			!write(*,'(3x,"15 - Napon dodira duz vise paralelnih pravaca istovremeno")')
			!write(*,'("         (3D raspodjela napona dodira)")')
			write(*,'("--------------------------------------------------------------")')
			write(*,'("Elektricno i magnetsko polje/indukcija:")')
			write(*,'("--------------------------------------------------------------")')
			write(*,'(4x,"8 - Elektricno polje u jednoj tocki viseslojnog sredstva")')
			write(*,'(4x,"9 - Elektricno polje duz pravca u viseslojnom sredstvu")')
			write(*,'(3x,"10 - 3D raspodjela elektricnog polja Ex, Ey i Ez, kao i 3D")')
			write(*,'("        raspodjela ukupnog elektricnog polja po bilo kojoj")')
			write(*,'("        plohi, ukljucujuci i povrsinu tla")')
			write(*,'(3x,"11 - Magnetsko polje u jednoj tocki viseslojnog sredstva")')
			write(*,'(3x,"12 - Raspodjela magnetskog polja Hx, Hy i Hz po bilo kojoj")')
			write(*,'("        plohi, ukljucujuci i povrsinu tla (3D raspodjela kompo-")')
			write(*,'("        nenti mag. polja). Ukljucen je i proracun ukupnog mag-")')
			write(*,'("        netskog polja.")')
			write(*,'(3x,"13 - Raspodjela magnet. indukcije Bx, By i Bz po bilo kojoj")')
			write(*,'("        plohi, ukljucujuci i povrsinu tla (3D raspodjela mag.")')
			write(*,'("        indukcije). Ukljucen je i proracun ukupne mag. induk.")')
			write(*,'("--------------------------------------------------------------")')
			write(*,'(4x,"0 - Zavrsetak proracuna i izlaz iz programa.")')
			write(*,'("--------------------------------------------------------------")')
			write(*,'("Izbor: ")',advance='no')
			read(*,*) izbor

			IF (izbor==1) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN RASPODJELE POTENCIJALA U ODABRANOJ TOCKI
				! ---------------------------------------------------------------------------
				write(*,'(/,"Unesite koordinate tocke za koju zelite izracunati potencijal:")')
				write(*,'("Xt = ")',advance='no')
				read(*,*) xt
				write(*,'("Yt = ")',advance='no')
				read(*,*) yt
				write(*,'("Zt = ")',advance='no')
				read(*,*) zt
				write(*,'("Sloj u kojem je tocka: ")',advance='no')
				read(*,*) sloj

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti1 = 1

				! Poziv rutine za proracun potencijala odabrane tocke
				! ..........................................................................
				call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
				r0,L,Ip,HD,h,kapa,gama,potencijal)
				! ..........................................................................
				
				call zmk(potencijal,modul,kut)
				write(*,'("Potencijal tocke (modul/kut):",f16.6," [V]",f8.2," (o)")') modul,kut
				
				! Tiskanje rezultata proracuna potencijala u odabranoj tocki
				write(2,'(/"Potencijal u odabranoj tocki (modul/kut):")')
				write(2,'(41("-"))')
				write(2,'("Sloj: ",i2,"  ->  ","T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V]")') sloj,xt,yt,zt,modul,kut

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time9 = end1 - start1
			
				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==2) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN RASPODJELE POTENCIJALA DUZ ODABRANOG PRAVCA 
				! ---------------------------------------------------------------------------
				write(*,'(/,"Koordinate pocetne tocke pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Zp = ")',advance='no')
				read(*,*) zp1
				write(*,'("Koordinate krajnje tocke pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				write(*,'("Zk = ")',advance='no')
				read(*,*) zk1
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				write(*,'("Korak duz osi z: ")',advance='no')
				read(*,*) dz
				write(*,'("Sloj u kojem je pravac: ")',advance='no')
				read(*,*) sloj

				write(*,'(/,"Racunam ...")')
				
				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti2 = 1
				
				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Fi(Ntoc)) !recikliranje varijable Fi
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					! Poziv rutine za proracun potencijala u jednoj tocki
					! ..................................................................
					call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
					r0,L,Ip,HD,h,kapa,gama,potencijal)
					! ..................................................................
					Fi(i) = potencijal
						
					xt = xt + dx
					yt = yt + dy
					zt = zt + dz
				end do

				! ................................................
				! PlotXY - graficki prikaz raspodjele potencijala
				! ................................................
				open(1,file='Potencijal_pravac.adf')  
				write(1,*)
				write(1,'("x",15x,"Fi(x)")')
				if (dx/=0.0) then
						t = xp1
				else if(dy/=0.0) then
						t = yp1
				else if(dz/=0.0) then
						t = zp1
				end if
				do i = 1,Ntoc
					call zmk(Fi(i),modul,kut)
					write(1,'(f10.3,f16.6)') t,modul
					if (dx/=0.0) then
						t = t + dx
					else if(dy/=0.0) then
						t = t + dy
					else if(dz/=0.0) then
						t = t + dz
					end if
				end do
				close(1)
				deallocate(Fi)

				write(*,'(/,"Napomena:",/,9("-"),/,"Rezultati za graficki prikaz snimljeni su u fileove:")')
				write(*,'(4x,"Potencijal_pravac.adf - graf potencijala duz pravca u visesloju")')
				write(*,'(4x,"Potencijal_tlo.adf    - graf potencijala duz pravca na povrsini tla")')
				write(*,'(4x,"Napon_koraka.adf      - graf napona koraka duz odabranog pravca")')
				write(*,'(4x,"El_polje_pravac.adf   - graf raspodjele iznosa elektricnog polja duz pravca")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time4 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==3) THEN
				! ---------------------------------------------------------------------------
				!        PRORACUN RASPODJELE POTENCIJALA DUZ PRAVCA NA POVRSINI TLA
				! ---------------------------------------------------------------------------
				write(*,'(/,"Koordinate pocetne tocke pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				zp1 = 0.d0
				write(*,'("Koordinate krajnje tocke pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				zk1 = 0.d0
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				dz = 0.d0
				sloj = 2

				write(*,'(/,"Racunam ...")')

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti3 = 1
				
				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Fi(Ntoc)) !recikliranje varijable Fi
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					! Poziv rutine za proracun potencijala u jednoj tocki
					! ..................................................................
					call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
					r0,L,Ip,HD,h,kapa,gama,potencijal)
					! ..................................................................
					Fi(i) = potencijal
						
					xt = xt + dx
					yt = yt + dy
				end do

				! ................................................
				! PlotXY - graficki prikaz raspodjele potencijala
				! ................................................
				open(1,file='Potencijal_tlo.adf')  
				write(1,*)
				write(1,'("x",15x,"Fi(x)")')
				if (dx/=0.0) then
						t = xp1
				else if(dy/=0.0) then
						t = yp1
				else if(dz/=0.0) then
						t = zp1
				end if
				do i = 1,Ntoc
					call zmk(Fi(i),modul,kut)
					write(1,'(f10.3,f16.6)') t,modul
					if (dx/=0.0) then
						t = t + dx
					else if(dy/=0.0) then
						t = t + dy
					else if(dz/=0.0) then
						t = t + dz
					end if
				end do
				close(1)
				deallocate(Fi)

				write(*,'(/,"Napomena:",/,9("-"),/,"Rezultati za graficki prikaz snimljeni su u fileove:")')
				write(*,'(4x,"Potencijal_pravac.adf - graf potencijala duz pravca u visesloju")')
				write(*,'(4x,"Potencijal_tlo.adf    - graf potencijala duz pravca na povrsini tla")')
				write(*,'(4x,"Napon_koraka.adf      - graf napona koraka duz odabranog pravca")')
				write(*,'(4x,"El_polje_pravac.adf   - graf raspodjele iznosa elektricnog polja duz pravca")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time5 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==4) THEN
				! ---------------------------------------------------------------------------
				!           PRORACUN NAPONA KORAKA DUZ PRAVCA NA POVRSINI TLA
				! ---------------------------------------------------------------------------
				write(*,'(/,"Koordinate pocetne tocke pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				zp1 = 0.0
				write(*,'("Koordinate krajnje tocke pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				zk1 = 0.0
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				dz = 0.0
				sloj = 2

				write(*,'("Koliko je vrijednosti dx (odnosno dy) sadrzano u 1 m (korak)?")')
				write(*,'("Korak: ")',advance='no')
				read(*,*) korak
				!Korak - varijabla koja definira koliko ima dx ili dy u jednom koraku koji
				!iznosi tocno jedan metar (1 m)!!!

				write(*,'(/,"Racunam ...")')

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti4 = 1

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Fi(Ntoc)) !recikliranje varijable Fi
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					! Poziv rutine za proracun potencijala u jednoj tocki
					! ..................................................................
					call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
					r0,L,Ip,HD,h,kapa,gama,potencijal)
					! ..................................................................
					Fi(i) = potencijal
						
					xt = xt + dx
					yt = yt + dy
					zt = zt + dz
				end do

				!PRORACUN NAPONA KORAKA DUZ ODABRANOG PRAVCA
				if(dx/=0.d0) then
					Nm = xk1 - xp1
				else if(dy/=0.d0) then
					Nm = yk1 - yp1
				end if
				allocate(Step(Nm))
				ik = korak
				M1 = dsqrt((dreal(Fi(1))**2+dimag(Fi(1))**2))
				M2 = dsqrt((dreal(Fi(ik))**2+dimag(Fi(ik))**2))
				Step(1) = dabs(M1-M2)
				do i = 2,Nm
					M1 = dsqrt((dreal(Fi(ik))**2+dimag(Fi(ik))**2))
					M2 = dsqrt((dreal(Fi(ik+korak))**2+dimag(Fi(ik+korak))**2))
					Step(i) = dabs(M1-M2)
					ik = ik + korak
				end do
				
				deallocate(Fi)

				!Priprema podataka za crtanje napona koraka
				xo = xp1
				yo = yp1
				zo = zp1
				allocate(xx(Nm))   !predstavlja udaljenost u razmaku od po 1 m
				if(dx/=0.d0) then
					xx(1) = xo
					do i = 2,Nm
						xx(i) = xx(i-1) + korak*dx
					end do
				else if(dy/=0.d0) then
					xx(1) = yo
					do i = 2,Nm
						xx(i) = xx(i-1) + korak*dy
					end do
				end if

				!Tiskanje izracunatih vrijednosti u file za PlotXY
				open(3,file='Napon_koraka.adf')
				write(3,'(/,5x,"x(m)",10x,"Vk(x)")')
				do i = 1,Nm
					write(3,'(e14.8,2x,e16.6)') xx(i),Step(i)
				end do
				close(3)
				
				deallocate(xx)
				deallocate(Step)

				write(*,'(/,"Napomena:",/,9("-"),/,"Rezultati za graficki prikaz snimljeni su u fileove:")')
				write(*,'(4x,"Potencijal_pravac.adf - graf potencijala duz pravca u visesloju")')
				write(*,'(4x,"Potencijal_tlo.adf    - graf potencijala duz pravca na povrsini tla")')
				write(*,'(4x,"Napon_koraka.adf      - graf napona koraka duz odabranog pravca")')
				write(*,'(4x,"El_polje_pravac.adf   - graf raspodjele iznosa elektricnog polja duz pravca")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time6 = end1 - start1
				
				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==5) THEN
				! ---------------------------------------------------------------------------
				! PRORACUN RASPODJELE POTENCIJALA DUZ VISE PARALELNIH PRAVACA NA POVRSINI TLA
				!                3D RASPODJELA POTENCIJALA PO POVRSINI TLA
				! ---------------------------------------------------------------------------
				write(*,'(/,"Ukupni broj paralelnih pravaca: ")',advance='no')
				read(*,*) Nprav
				write(*,'("Pravci paralelni koordinatnoj osi x ili y (x/y): ")',advance='no')
				read(*,*) paralelan
				write(*,'("Razmak medju paralelnim pravcima: ")',advance='no')
				read(*,*) razmak
				write(*,'(/,"Koordinate pocetne tocke prvog pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				zp1 = 0.d0
				write(*,'("Koordinate krajnje tocke prvog pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				zk1 = 0.d0

				if ((paralelan=='x').or.(paralelan=='X')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dx
					dy = 0.d0
				else if ((paralelan=='y').or.(paralelan=='Y')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dy
					dx = 0.d0
				end if
				dz = 0.d0
				sloj = 2

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti5 = 1

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if

				! ................................................
				!      Proracun raspodjele potencijala
				! ................................................
				write(*,'(/,"Racunam ...")')
				allocate(FF(Nprav,Ntoc))
				xtp = xp1
				ytp = yp1
				zt = 0.d0
				do j = 1,Nprav
					allocate(Fi(Ntoc)) !recikliranje varijable Fi
					xt = xtp
					yt = ytp
					do i = 1,Ntoc
						! Poziv rutine za proracun potencijala u jednoj tocki
						! ..................................................................
						call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
						r0,L,Ip,HD,h,kapa,gama,potencijal)
						! ..................................................................
						Fi(i) = potencijal
							
						xt = xt + dx
						yt = yt + dy
					end do
					do i = 1,Ntoc
						modul = dsqrt((dreal(Fi(i))**2+dimag(Fi(i))**2))
						FF(j,i) = modul
					end do
					deallocate(Fi)
					
					if ((paralelan=='x').or.(paralelan=='X')) then
						ytp = ytp + razmak
					else if ((paralelan=='y').or.(paralelan=='Y')) then
						xtp = xtp + razmak
					end if
				end do

				!-----------------------------------------------------------------------
				!       3D prikaz izracunatih potencijala - Gnuplot
				!-----------------------------------------------------------------------
				open(7,file='3Dpotencijal.dat')
				do j = 1,Nprav
					do i = 1,Ntoc
						write(7,'(e12.4)') FF(j,i)
					end do
					write(7,*)
				end do
				close(7)
				!-----------------------------------------------------------------------
				!       3D prikaz izracunatih potencijala - LabPlot
				!-----------------------------------------------------------------------
				open(7,file='3Dpotenc.dat')
				do i = 1,Ntoc
					write(7,'(100e12.4)') (FF(j,i),j=1,Nprav)
				end do
				close(7)
				! ...............................................
				! Snimanje podataka za graficki prikaz : ParaView
				! ...............................................
				! Priprema koord. osi x i y za crtanje 3D raspodjele potencijala
				if ((paralelan=='y').or.(paralelan=='Y')) then
					allocate(Yos(Ntoc))
					allocate(Xos(Nprav))
					yt = yp1
					do i = 1,Ntoc
						Yos(i) = yt
						yt = yt + dy
					end do
					xt = xp1
					do i = 1,Nprav
						Xos(i) = xt
						xt = xt + razmak
					end do
				else if((paralelan=='x').or.(paralelan=='X')) then
					allocate(Yos(Nprav))
					allocate(Xos(Ntoc))
					xt = xp1
					do i = 1,Ntoc
						Xos(i) = xt
						xt = xt + dx
					end do
					yt = yp1
					do i = 1,Nprav
						Yos(i) = yt
						yt = yt + razmak
					end do
				end if
				!ParaView format
				open(99,file='3DPotenc_para.vtk')
				write(99,'("# vtk DataFile Version 2.0")')
				write(99,'("3D prikaz raspodjele potencijala po plohi")')
				write(99,'("ASCII")')
				write(99,*)
				write(99,'("DATASET STRUCTURED_GRID")')
				write(99,'("DIMENSIONS",i5,i5,i5)') Ntoc,Nprav,1
				write(99,'("POINTS",i5," float")') Nprav*Ntoc
				do j = 1,Nprav
					do i = 1,Ntoc
						write(99,'(es10.3,2x,es10.3,2x,es10.3)') Xos(i),Yos(j),FF(j,i)*1.d-3 !(kV)
					end do
				end do
				write(99,*)
				write(99,'("POINT_DATA",i5)') Nprav*Ntoc
				write(99,*)
				write(99,'("SCALARS",2x,"potenc"," float"," 1")')
				write(99,'("LOOKUP_TABLE default")')
				do i = 1,Nprav
					do j = 1,Ntoc
						write(99,'(es10.3)') FF(i,j)*1.d-3 !(kV)
					end do
				end do
				close(99)
				deallocate(Xos)
				deallocate(Yos)
				deallocate(FF)
				
				write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D potencijala snimljeni su u file: 3Dpotencijal.dat.")')
				write(*,'("Za graficki prikaz koristiti: Gnuplot. Za LabPlot: 3Dpotenc.dat")')
				write(*,'("Za ParaView koristiti file: 3DPotenc_para.vtk")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time7 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==6) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN NAPONA DODIRA U ODABRANOJ TOCKI 
				! ---------------------------------------------------------------------------
				write(*,'(/,"Koordinate cvorova uzemljivaca:")')
				write(*,'(50("-"))')
				write(*,'("Br.",t10,"x",t20,"y",t30,"z",t40,"Potencijal")')
				write(*,'("cvor",t9,"(m)",t19,"(m)",t29,"(m)",t44,"(V)")')
				write(*,'(50("-"))')
				open(5,file='glob_pot.dat')
				nG = 0
				do while(.not.eof(5))
					read(5,'(i3,3f8.2,2e16.6)') i,xl,yl,zl,modul,kut !recikliranje modul,kut
					write(*,'(i3,t5,f8.2,t15,f8.2,t25,f8.2,t35,f16.6)') i,xl,yl,zl,&
					dsqrt(modul**2+kut**2)
					nG = nG + 1
				end do
				close(5)
				write(*,'(50("-"))')

				allocate(ipiv(nG))
				allocate(Fig(nG))
				open(5,file='glob_pot.dat')
				do i = 1,nG
					read(5,'(i3,3f8.2,2e16.6)') ipiv(i),xl,yl,zl,&
					modul,kut !recikliranje modul,kut
					Fig(i) = dcmplx(modul,kut)
				end do
				close(5)

				write(*,'("Upisite broj cvora za proracun napona dodira: ")',advance='no')
				read(*,*) broj_cvora
				
				do i = 1,nG
					if (broj_cvora==ipiv(i)) then
						FFig = Fig(broj_cvora)
					end if
				end do
				
				deallocate(ipiv)
				deallocate(Fig)

				write(*,'(/,"Koordinate tocke na povrsini tla za napon dodira:")')
				write(*,'("Xt = ")',advance='no')
				read(*,*) xt
				write(*,'("Yt = ")',advance='no')
				read(*,*) yt
				zt = 0.d0
				sloj = 2

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti6 = 1

				! Poziv rutine za proracun potencijala odabrane tocke
				! ..................................................................
				call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
				r0,L,Ip,HD,h,kapa,gama,potencijal)
				! ..................................................................
				
				! Napon dodira
				Nap_dod = FFig - potencijal
				
				call zmk(Nap_dod,modul,kut)

				write(*,'(/"Napon dodira u odabranoj tocki (modul/kut):")')
				write(*,'("Cvor: ",i2,"  ->  ","T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V]")') broj_cvora,xt,yt,zt,modul,kut
				
				! Tiskanje napona dodira u odabranoj tocki
				write(2,'(/"Napon dodira u odabranoj tocki (modul/kut):")')
				write(2,'(41("-"))')
				write(2,'("Cvor: ",i2,"  ->  ","T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V]")') broj_cvora,xt,yt,zt,modul,kut

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time10 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==7) THEN
				! ---------------------------------------------------------------------------
				! PRORACUN NAPONA KORAKA DUZ VISE PARALELNIH PRAVACA ISTOVREMENO (3D)
				! ---------------------------------------------------------------------------
				write(*,'(/,"Ukupni broj paralelnih pravaca: ")',advance='no')
				read(*,*) Nprav
				write(*,'("Pravci paralelni koordinatnoj osi x ili y (x/y): ")',advance='no')
				read(*,*) paralelan
				write(*,'("Razmak medju paralelnim pravcima: ")',advance='no')
				read(*,*) razmak
				write(*,'(/,"Koordinate pocetne tocke prvog pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				zp1 = 0.d0
				write(*,'("Koordinate krajnje tocke prvog pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				zk1 = 0.d0

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti7 = 1

				if ((paralelan=='x').or.(paralelan=='X')) then
					dx = 1.d0; dy = 0.d0
				else if ((paralelan=='y').or.(paralelan=='Y')) then
					dy = 1.d0; dx = 0.d0
				end if
				dz = 0.d0
				sloj = 2
						
				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if

				! ................................................
				!      Proracun raspodjele potencijala
				! ................................................
				write(*,'(/,"Racunam ...")')
				allocate(FF(Nprav+1,Ntoc+1))
				xtp = xp1
				ytp = yp1
				zt = 0.d0
				do j = 1,Nprav+1
					allocate(Fi(Ntoc+1)) !recikliranje varijable Fi
					xt = xtp
					yt = ytp
					do i = 1,Ntoc+1
						! Poziv rutine za proracun potencijala u jednoj tocki
						! ..................................................................
						call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
						r0,L,Ip,HD,h,kapa,gama,potencijal)
						! ..................................................................
						Fi(i) = potencijal
							
						xt = xt + dx
						yt = yt + dy
					end do
					do i = 1,Ntoc+1
						modul = dsqrt((dreal(Fi(i))**2+dimag(Fi(i))**2))
						FF(j,i) = modul
					end do
					deallocate(Fi)
					
					if ((paralelan=='x').or.(paralelan=='X')) then
						ytp = ytp + razmak
					else if ((paralelan=='y').or.(paralelan=='Y')) then
						xtp = xtp + razmak
					end if
				end do

				! --------------------------------------------------------------------
				! Proracun napona koraka duz paralelnih pravaca - proracun se provodi 
				! posebno ili u smjeru osi x ili u smjeru osi y.
				! --------------------------------------------------------------------
				allocate(napon_koraka(Nprav,Ntoc))
				do j = 1,Nprav
					do i = 1,Ntoc
						napon_koraka(j,i) = dabs(FF(j,i)-FF(j,i+1)) 
					end do
				end do

				deallocate(FF)

				!-----------------------------------------------------------------------
				!       3D prikaz izracunatih potencijala - Gnuplot
				!-----------------------------------------------------------------------
				open(7,file='3Dnaponkoraka.dat')
				do j = 1,Nprav
					do i = 1,Ntoc
						write(7,'(e12.4)') FF(j,i)
					end do
					write(7,*)
				end do
				close(7)		
				!-----------------------------------------------------------------------
				!       3D prikaz izracunatih potencijala - LabPlot
				!-----------------------------------------------------------------------
				open(7,file='3Dnapon_koraka.dat')
				do i = 1,Ntoc
					write(7,'(100e12.4)') (FF(j,i),j=1,Nprav)
				end do
				close(7)	
				deallocate(FF)
				
				write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D napona koraka snimljeni su u file: 3Dnaponkoraka.dat.")')
				write(*,'("Za graficki prikaz koristiti: Gnuplot. Za LabPlot: 3Dnapon_koraka.dat")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time8 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==8) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN ELEKTRICNOG POLJA U ODABRANOJ TOCKI
				! ---------------------------------------------------------------------------
				write(*,'(/,"Unesite koordinate tocke za koju zelite izracunati elektricno polje:")')
				write(*,'("Xt = ")',advance='no')
				read(*,*) xt
				write(*,'("Yt = ")',advance='no')
				read(*,*) yt
				write(*,'("Zt = ")',advance='no')
				read(*,*) zt
				write(*,'("Sloj u kojem je tocka: ")',advance='no')
				read(*,*) sloj

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti8 = 1

				! Poziv rutine za proracun Ex komponente elektricnog polja u odabranoj tocki
				! ..........................................................................
				call Epolje_Fix(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
				kapa,gama,Fix)
				call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
				gama,f,Ax,Ay,Az)
				! ..........................................................................
				Ex = -Fix - Ax
				! Poziv rutine za proracun Ey komponente elektricnog polja u odabranoj tocki
				! ..........................................................................
				call Epolje_Fiy(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
				kapa,gama,Fiy)
				call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
				gama,f,Ax,Ay,Az)
				! ..........................................................................
				Ey = -Fiy - Ay
				! Poziv rutine za proracun Ez komponente elektricnog polja u odabranoj tocki
				! ..........................................................................
				call Epolje_Fiz(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
				kapa,gama,Fiz)
				call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
				gama,f,Ax,Ay,Az)			
				! ..........................................................................
				Ez = -Fiz - Az

				! Tiskanje Ex:
				call zmk(Ex,modulx,kut)
				write(*,'(/,"Komponenta el. polja: Ex (modul/kut):",f16.6,f8.2," [V/m]")') modulx,kut
				! Tiskanje rezultata proracuna komponente Ex polja u odabranoj tocki
				write(2,'(/,"Elektricno polje u odabranoj tocki (modul/kut):")')
				write(2,'(67("-"))')
				write(2,'("Sloj: ",i2,"  ->  ","Ex: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V/m]")') sloj,xt,yt,zt,modulx,kut
				! Tiskanje Ey:
				call zmk(Ey,moduly,kut)
				write(*,'("Komponenta el. polja: Ey (modul/kut):",f16.6,f8.2," [V/m]")') moduly,kut
				! Tiskanje rezultata proracuna komponente Ey polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","Ey: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V/m]")') sloj,xt,yt,zt,moduly,kut				
				! Tiskanje Ez:
				call zmk(Ez,modulz,kut)
				write(*,'("Komponenta el. polja: Ez (modul/kut):",f16.6,f8.2," [V/m]")') modulz,kut
				! Tiskanje rezultata proracuna komponente Ex polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","Ez: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [V/m]")') sloj,xt,yt,zt,modulz,kut
				! Tiskanje rezultirajuceg elektricnog polja:
				modul = dsqrt(modulx**2+moduly**2+modulz**2)
				write(*,'("Rezultirajuce elektricno polje (modul):",f16.6" [V/m]")') modul
				write(2,'("Rezultirajuce elektricno polje (modul):",f16.6" [V/m]")') modul
				write(2,'(67("-"))')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time11 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==9) THEN
				! ---------------------------------------------------------------------------
				!        PRORACUN RASPODJELE ELEKTRICNOG POLJA DUZ ODABRANOG PRAVCA 
				! ---------------------------------------------------------------------------
				write(*,'(/,"Koordinate pocetne tocke pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Zp = ")',advance='no')
				read(*,*) zp1
				write(*,'("Koordinate krajnje tocke pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				write(*,'("Zk = ")',advance='no')
				read(*,*) zk1
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				write(*,'("Korak duz osi z: ")',advance='no')
				read(*,*) dz
!				write(*,'("Sloj u kojem je pravac: ")',advance='no')
!				read(*,*) sloj

				write(*,'(/,"Racunam ...")')
				
				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti9 = 1
				
				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Exx(Ntoc))
				allocate(Eyy(Ntoc))
				allocate(Ezz(Ntoc)) 
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					if (N==2) then
						if (zt>=0.d0) then
							sloj = 2
						else
							sloj = 1
						end if
					else
						!Odredjivanje polozaja tocke kod viseslojnog sredstva
						!da se moze uzeti pravac koji ide kroz vise slojeva
						sloj = 1
						do j = 1,N-1
							if (zt>HD(j)) sloj = sloj + 1
						end do
					end if
					! Poziv rutine za proracun Ex komponente elektricnog polja u odabranoj tocki
					! ..........................................................................
					call Epolje_Fix(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
					kapa,gama,Fix)
					call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
					gama,f,Ax,Ay,Az)
					! ..........................................................................
					Exx(i) = -Fix - Ax

					! Poziv rutine za proracun Ey komponente elektricnog polja u odabranoj tocki
					! ..........................................................................
					call Epolje_Fiy(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
					kapa,gama,Fiy)
					call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
					gama,f,Ax,Ay,Az)
					! ..........................................................................
					Eyy(i) = -Fiy - Ay

					! Poziv rutine za proracun Ez komponente elektricnog polja u odabranoj tocki
					! ..........................................................................
					call Epolje_Fiz(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
					kapa,gama,Fiz)
					call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
					gama,f,Ax,Ay,Az)			
					! ..........................................................................
					Ezz(i) = -Fiz - Az
									
					xt = xt + dx
					yt = yt + dy
					zt = zt + dz
				end do

				! .....................................................
				! PlotXY - graficki prikaz raspodjele elektricnog polja
				! .....................................................
				open(1,file='El_polje_pravac.adf')  
				write(1,*)
				write(1,'("x",15x,"Ex",15x,"Ey",15x,"Ez",15x,"Rezult.polje")')
				if (dx/=0.0) then
						t = xp1
				else if(dy/=0.0) then
						t = yp1
				else if(dz/=0.0) then
						t = zp1
				end if
				do i = 1,Ntoc
					modulx = dsqrt((dreal(Exx(i))**2+dimag(Exx(i))**2))
					moduly = dsqrt((dreal(Eyy(i))**2+dimag(Eyy(i))**2))
					modulz = dsqrt((dreal(Ezz(i))**2+dimag(Ezz(i))**2))
					modul = dsqrt(modulx**2+moduly**2+modulz**2)
					write(1,'(f10.3,f16.6,f16.6,f16.6,f16.6)') t,modulx,moduly,modulz,modul
					if (dx/=0.0) then
						t = t + dx
					else if(dy/=0.0) then
						t = t + dy
					else if(dz/=0.0) then
						t = t + dz
					end if
				end do
				close(1)
				deallocate(Exx,Eyy,Ezz)

				write(*,'(/,"Napomena:",/,9("-"),/,"Rezultati za graficki prikaz snimljeni su u fileove:")')
				write(*,'(4x,"Potencijal_pravac.adf - graf potencijala duz pravca u visesloju")')
				write(*,'(4x,"Potencijal_tlo.adf    - graf potencijala duz pravca na povrsini tla")')
				write(*,'(4x,"Napon_koraka.adf      - graf napona koraka duz odabranog pravca")')
				write(*,'(4x,"El_polje_pravac.adf   - graf raspodjele iznosa elektricnog polja duz pravca")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time12 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==10) THEN
				! ---------------------------------------------------------------------------
				! PRORACUN RASPODJELE ELEKTRICNOG POLJA DUZ VISE PARALELNIH PRAVACA PO BILO
				! KOJOJ PLOHI (UKLJUCUJUCI I POVRSINU TLA) - 3D RASPODJELA ELEKTRICNOG POLJA
				! ---------------------------------------------------------------------------
				write(*,'(/,"Ukupni broj paralelnih pravaca: ")',advance='no')
				read(*,*) Nprav
				write(*,'("Pravci paralelni koordinatnoj osi x ili y (x/y): ")',advance='no')
				read(*,*) paralelan
				write(*,'("Razmak medju paralelnim pravcima: ")',advance='no')
				read(*,*) razmak
				write(*,'(/,"Koordinate pocetne tocke prvog pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Koordinate krajnje tocke prvog pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				if ((paralelan=='x').or.(paralelan=='X')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dx
					dy = 0.d0
				else if ((paralelan=='y').or.(paralelan=='Y')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dy
					dx = 0.d0
				end if
				dz = 0.d0
				write(*,'(/,"Z koordinata plohe na kojoj leze pravci:")')
				write(*,'("   - negativne vrijednosti za plohu u zraku")')
				write(*,'("   - pozitivne vrijednosti za plohu u tlu")')
				write(*,'("   - vrijednost nula za povrsinu tla")')
				write(*,'("-------------------------------------------")')
				write(*,'("Z = ")',advance='no')
				read(*,*) zp1
				if (zp1<0.d0) then !ZRAK
					zp1 = zp1
					zk1 = zp1
					sloj = 1
				else if(zp1>0.d0) then !TLO
					zp1 = zp1
					zk1 = zp1
					write(*,'("Sloj u kojem je ploha: ")',advance='no')
					read(*,*) sloj
				else if(zp1==0.d0) then !POVRSINA TLA
					zp1 = 0.0
					zk1 = 0.0
					sloj = 2
				end if

				write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
				write(*,'("-------------------------------------------------")')
				write(*,'("(1) 3D raspodjela Ex komponente elektricnog polja")')
				write(*,'("(2) 3D raspodjela Ey komponente elektricnog polja")')
				write(*,'("(3) 3D raspodjela Ez komponente elektricnog polja")')
				write(*,'("(4) 3D raspodjela ukupnog/totalnog  elektr. polja")')
				write(*,'("-------------------------------------------------")')
				write(*,'("Izbor: ")',advance='no')
				read(*,*) izbor

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti10 = 1

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if

				! ................................................
				!      Proracun raspodjele elektricnog polja
				! ................................................
				write(*,'(/,"Racunam ...")')
				allocate(FF(Nprav,Ntoc))
				xtp = xp1
				ytp = yp1
				if (zp1<0.d0) then !ZRAK
					zt = zp1
				else if(zp1>0.d0) then !TLO
					zt = zp1
				else if(zp1==0.d0) then !POVRSINA TLA
					zt = 0.001 !1 mm
				end if
				do j = 1,Nprav
					allocate(Exx(Ntoc))
					allocate(Eyy(Ntoc))
					allocate(Ezz(Ntoc))
					xt = xtp
					yt = ytp
					do i = 1,Ntoc
						! Poziv rutine za proracun Ex komponente elektricnog polja u odabranoj tocki
						! ..........................................................................
						call Epolje_Fix(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
						kapa,gama,Fix)
						call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
						gama,f,Ax,Ay,Az)
						! ..........................................................................
						Exx(i) = -Fix - Ax

						! Poziv rutine za proracun Ey komponente elektricnog polja u odabranoj tocki
						! ..........................................................................
						call Epolje_Fiy(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
						kapa,gama,Fiy)
						call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
						gama,f,Ax,Ay,Az)
						! ..........................................................................
						Eyy(i) = -Fiy - Ay

						! Poziv rutine za proracun Ez komponente elektricnog polja u odabranoj tocki
						! ..........................................................................
						call Epolje_Fiz(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,r0,L,Ip,HD,h,&
						kapa,gama,Fiz)
						call Epolje_Ax_Ay_Az(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,&
						gama,f,Ax,Ay,Az)			
						! ..........................................................................
						Ezz(i) = -Fiz - Az
							
						xt = xt + dx
						yt = yt + dy
					end do

					select case(izbor)
						case(1)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Exx(i))**2+dimag(Exx(i))**2))
								FF(j,i) = modulx
							end do
							deallocate(Exx,Eyy,Ezz)
						case(2)
							do i = 1,Ntoc
								moduly = dsqrt((dreal(Eyy(i))**2+dimag(Eyy(i))**2))
								FF(j,i) = moduly
							end do
							deallocate(Exx,Eyy,Ezz)
						case(3)
							do i = 1,Ntoc
								modulz = dsqrt((dreal(Ezz(i))**2+dimag(Ezz(i))**2))
								FF(j,i) = modulz
							end do
							deallocate(Exx,Eyy,Ezz)
						case(4)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Exx(i))**2+dimag(Exx(i))**2))
								moduly = dsqrt((dreal(Eyy(i))**2+dimag(Eyy(i))**2))
								modulz = dsqrt((dreal(Ezz(i))**2+dimag(Ezz(i))**2))
								FF(j,i) = dsqrt(modulx**2+moduly**2+modulz**2)
							end do
							deallocate(Exx,Eyy,Ezz)
						case default
							write(*,'(/,"Pogresan odabir!!!")')
							goto 100
					end select

					if ((paralelan=='x').or.(paralelan=='X')) then
						ytp = ytp + razmak
					else if ((paralelan=='y').or.(paralelan=='Y')) then
						xtp = xtp + razmak
					end if
				end do

				!-----------------------------------------------------------------------
				!       3D prikaz izracunatog elektricnog polja - Gnuplot
				!-----------------------------------------------------------------------
				open(7,file='3D_el_polje.dat')
				do j = 1,Nprav
					do i = 1,Ntoc
						write(7,'(e12.4)') FF(j,i)
					end do
					write(7,*)
				end do
				close(7)	
				!-----------------------------------------------------------------------
				!       3D prikaz izracunatog elektricnog polja - LabPlot
				!-----------------------------------------------------------------------
				open(7,file='3Delek_polje.dat')
				do i = 1,Ntoc
					write(7,'(100e12.4)') (FF(j,i),j=1,Nprav)
				end do
				close(7)
				! ...............................................
				! Snimanje podataka za graficki prikaz : ParaView
				! ...............................................
				! Priprema koord. osi x i y za crtanje 3D raspodjele elektricnog polja
				if ((paralelan=='y').or.(paralelan=='Y')) then
					allocate(Yos(Ntoc))
					allocate(Xos(Nprav))
					yt = yp1
					do i = 1,Ntoc
						Yos(i) = yt
						yt = yt + dy
					end do
					xt = xp1
					do i = 1,Nprav
						Xos(i) = xt
						xt = xt + razmak
					end do
				else if((paralelan=='x').or.(paralelan=='X')) then
					allocate(Yos(Nprav))
					allocate(Xos(Ntoc))
					xt = xp1
					do i = 1,Ntoc
						Xos(i) = xt
						xt = xt + dx
					end do
					yt = yp1
					do i = 1,Nprav
						Yos(i) = yt
						yt = yt + razmak
					end do
				end if
				!ParaView format
				open(99,file='3DEpolje_para.vtk')
				write(99,'("# vtk DataFile Version 2.0")')
				write(99,'("3D prikaz raspodjele el. polja po plohi")')
				write(99,'("ASCII")')
				write(99,*)
				write(99,'("DATASET STRUCTURED_GRID")')
				write(99,'("DIMENSIONS",i5,i5,i5)') Ntoc,Nprav,1
				write(99,'("POINTS",i5," float")') Nprav*Ntoc
				do j = 1,Nprav
					do i = 1,Ntoc
						write(99,'(es10.3,2x,es10.3,2x,es10.3)') Xos(i),Yos(j),FF(j,i)
					end do
				end do
				write(99,*)
				write(99,'("POINT_DATA",i5)') Nprav*Ntoc
				write(99,*)
				write(99,'("SCALARS",2x,"epolje"," float"," 1")')
				write(99,'("LOOKUP_TABLE default")')
				do i = 1,Nprav
					do j = 1,Ntoc
						write(99,'(es10.3)') FF(i,j)
					end do
				end do
				close(99)

				deallocate(Xos)
				deallocate(Yos)
				deallocate(FF)
				
				write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele elektricnog polja snimljeni su u file:")')
				write(*,'("3D_el_polje.dat. Za graficki prikaz koristiti: Gnuplot. Za LabPlot: 3Delek_polje.dat")')
				write(*,'("Za ParaView koristiti file: 3DEpolje_para.vtk")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time13 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==11) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN MAGNETSKOG POLJA U ODABRANOJ TOCKI
				! ---------------------------------------------------------------------------
				write(*,'(/,"Unesite koordinate tocke za koju zelite izracunati magnetsko polje:")')
				write(*,'("Xt = ")',advance='no')
				read(*,*) xt
				write(*,'("Yt = ")',advance='no')
				read(*,*) yt
				write(*,'("Zt = ")',advance='no')
				read(*,*) zt
				write(*,'("Sloj u kojem je tocka: ")',advance='no')
				read(*,*) sloj

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti11 = 1

				! Poziv rutine za proracun Bx, By i Bz komponente magnetskog polja u odabranoj tocki
				! ..........................................................................
				call Hpolje_Bx(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bx)
				call Hpolje_By(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,By)
				call Hpolje_Bz(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bz)
				! ..........................................................................

				! Tiskanje Bx/Hx:
				call zmk(Bx,modulx,kut)
				write(*,'(/,"Komponenta mag. polja: Bx (modul/kut):",f16.6,f8.2," [mT]")') modulx*1.e3,kut
				! Tiskanje rezultata proracuna komponente Bx polja u odabranoj tocki
				write(2,'(/,"Magnetsko polje u odabranoj tocki (modul/kut):")')
				write(2,'(67("-"))')
				write(2,'("Sloj: ",i2,"  ->  ","Bx: T(",f4.1,";",f4.1,";",f4.1,") =",&
				e16.6,f8.2" [T]")') sloj,xt,yt,zt,modulx*1.e3,kut
				! Tiskanje By/Hy:
				call zmk(By,moduly,kut)
				write(*,'("Komponenta mag. polja: By (modul/kut):",f16.6,f8.2," [mT]")') moduly*1.e3,kut
				! Tiskanje rezultata proracuna komponente By polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","By: T(",f4.1,";",f4.1,";",f4.1,") =",&
				e16.6,f8.2" [T]")') sloj,xt,yt,zt,moduly*1.e3,kut				
				! Tiskanje Bz/Hz:
				call zmk(Bz,modulz,kut)
				write(*,'("Komponenta mag. polja: Bz (modul/kut):",f16.6,f8.2," [mT]")') modulz*1.e3,kut
				! Tiskanje rezultata proracuna komponente Bz polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","Bz: T(",f4.1,";",f4.1,";",f4.1,") =",&
				e16.6,f8.2" [T]")') sloj,xt,yt,zt,modulz*1.e3,kut
				! Tiskanje rezultirajuce magnetske indukcije:
				modul = dsqrt(modulx**2+moduly**2+modulz**2)
				write(*,'("Rezultirajuca magnet. indukcija (modul):",f16.6" [mT]")') modul
				write(2,'("Rezultirajuca magnet. indukcija (modul):",f16.6" [mT]")') modul
				write(2,'(67("-"))')
				
				mio = 4.d0*pi*1.e-7
				Bx = Bx/mio !Hx
				By = By/mio !Hy
				Bz = Bz/mio !Hz
			
				! Tiskanje Bx/Hx:
				call zmk(Bx,modulx,kut)
				write(*,'(/,"Komponenta mag. polja: Hx (modul/kut):",f16.6,f8.2," [A/m]")') modulx,kut
				! Tiskanje rezultata proracuna komponente Hx polja u odabranoj tocki
				write(2,'(/,"Magnetsko polje u odabranoj tocki (modul/kut):")')
				write(2,'(67("-"))')
				write(2,'("Sloj: ",i2,"  ->  ","Hx: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [A/m]")') sloj,xt,yt,zt,modulx,kut
				! Tiskanje By/Hy:
				call zmk(By,moduly,kut)
				write(*,'("Komponenta mag. polja: Hy (modul/kut):",f16.6,f8.2," [A/m]")') moduly,kut
				! Tiskanje rezultata proracuna komponente Hy polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","Hy: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [A/m]")') sloj,xt,yt,zt,moduly,kut				
				! Tiskanje Bz/Hz:
				call zmk(Bz,modulz,kut)
				write(*,'("Komponenta mag. polja: Hz (modul/kut):",f16.6,f8.2," [A/m]")') modulz,kut
				! Tiskanje rezultata proracuna komponente Hz polja u odabranoj tocki
				write(2,'("Sloj: ",i2,"  ->  ","Hz: T(",f4.1,";",f4.1,";",f4.1,") =",&
				f16.6,f8.2" [A/m]")') sloj,xt,yt,zt,modulz,kut
				! Tiskanje rezultirajuceg magnetskog polja:
				modul = dsqrt(modulx**2+moduly**2+modulz**2)
				write(*,'("Rezultirajuce magnetsko polje (modul):",f16.6" [A/m]")') modul
				write(2,'("Rezultirajuce magnetsko polje (modul):",f16.6" [A/m]")') modul
				write(2,'(67("-"))')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time14 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==12) THEN
				! ---------------------------------------------------------------------------
				!   PRORACUN RASPODJELE KOMPONENTI MAGNETSKOG POLJA DUZ VISE PARALELNIH 
				!        PRAVACA NA NEKOJ PLOHI (3D RASPODJELA MAGNETSKOG POLJA)
				! ---------------------------------------------------------------------------
				write(*,'(/,"Ukupni broj paralelnih pravaca: ")',advance='no')
				read(*,*) Nprav
				write(*,'("Pravci paralelni koordinatnoj osi x ili y (x/y): ")',advance='no')
				read(*,*) paralelan
				write(*,'("Razmak medju paralelnim pravcima: ")',advance='no')
				read(*,*) razmak
				write(*,'(/,"Koordinate pocetne tocke prvog pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Koordinate krajnje tocke prvog pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				if ((paralelan=='x').or.(paralelan=='X')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dx
					dy = 0.d0
				else if ((paralelan=='y').or.(paralelan=='Y')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dy
					dx = 0.d0
				end if
				dz = 0.d0
				write(*,'(/,"Z koordinata plohe na kojoj leze pravci:")')
				write(*,'("   - negativne vrijednosti za plohu u zraku")')
				write(*,'("   - pozitivne vrijednosti za plohu u tlu")')
				write(*,'("   - vrijednost nula za povrsinu tla")')
				write(*,'("-------------------------------------------")')
				write(*,'("Z = ")',advance='no')
				read(*,*) zp1
				if (zp1<0.d0) then !ZRAK
					zp1 = zp1
					zk1 = zp1
					sloj = 1
				else if(zp1>0.d0) then !TLO
					zp1 = zp1
					zk1 = zp1
					write(*,'("Sloj u kojem je ploha: ")',advance='no')
					read(*,*) sloj
				else if(zp1==0.d0) then !POVRSINA TLA
					zp1 = zp1
					zk1 = zp1
					sloj = 2
				end if

				write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
				write(*,'("------------------------------------------------")')
				write(*,'("(1) 3D raspodjela Hx komponente magnetskog polja")')
				write(*,'("(2) 3D raspodjela Hy komponente magnetskog polja")')
				write(*,'("(3) 3D raspodjela Hz komponente magnetskog polja")')
				write(*,'("(4) 3D raspodjela ukupnog/totalnog magnet. polja")')
				write(*,'("------------------------------------------------")')
				write(*,'("Izbor: ")',advance='no')
				read(*,*) izbor

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti12 = 1

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if

				! ................................................
				!      Proracun raspodjele magnetskog polja
				! ................................................
				write(*,'(/,"Racunam ...")')
				allocate(FF(Nprav,Ntoc))
				xtp = xp1
				ytp = yp1
				zt = zp1
				mio = 4.d0*pi*1.e-7
				do j = 1,Nprav
					allocate(Bxx(Ntoc))
					allocate(Byy(Ntoc))
					allocate(Bzz(Ntoc))
					xt = xtp
					yt = ytp
					do i = 1,Ntoc
						! Poziv rutine za proracun Bx, By i Bz komponente magnetskog polja u odabranoj tocki
						! ..........................................................................
						call Hpolje_Bx(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bx)
						call Hpolje_By(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,By)
						call Hpolje_Bz(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bz)
						! ..........................................................................
						Bx = Bx/mio !Hx
						By = By/mio !Hy
						Bz = Bz/mio !Hz
						! Ovdje sada Bxx,Byy i Bzz predstavljaju
						! Hx, Hy i Hz komponente magentskog polja
						Bxx(i) = Bx
						Byy(i) = By
						Bzz(i) = Bz
							
						xt = xt + dx
						yt = yt + dy
					end do

					select case(izbor)
						case(1)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Bxx(i))**2+dimag(Bxx(i))**2))
								FF(j,i) = modulx
							end do
							deallocate(Bxx,Byy,Bzz)
						case(2)
							do i = 1,Ntoc
								moduly = dsqrt((dreal(Byy(i))**2+dimag(Byy(i))**2))
								FF(j,i) = moduly
							end do
							deallocate(Bxx,Byy,Bzz)
						case(3)
							do i = 1,Ntoc
								modulz = dsqrt((dreal(Bzz(i))**2+dimag(Bzz(i))**2))
								FF(j,i) = modulz
							end do
							deallocate(Bxx,Byy,Bzz)
						case(4)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Bxx(i))**2+dimag(Bxx(i))**2))
								moduly = dsqrt((dreal(Byy(i))**2+dimag(Byy(i))**2))
								modulz = dsqrt((dreal(Bzz(i))**2+dimag(Bzz(i))**2))
								FF(j,i) = dsqrt(modulx**2+moduly**2+modulz**2)
							end do
							deallocate(Bxx,Byy,Bzz)
						case default
							write(*,'(/,"Pogresan odabir!!!")')
							goto 100
					end select
					
					if ((paralelan=='x').or.(paralelan=='X')) then
						ytp = ytp + razmak
					else if ((paralelan=='y').or.(paralelan=='Y')) then
						xtp = xtp + razmak
					end if
				end do

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time15 = end1 - start1

				!-----------------------------------------------------------------------
				!       3D prikaz izracunatog magnetskog polja - Gnuplot
				!-----------------------------------------------------------------------
				open(7,file='3D_mag_polje.dat')
				do j = 1,Nprav
					do i = 1,Ntoc
						write(7,'(e12.4)') FF(j,i)
					end do
					write(7,*)
				end do
				close(7)	
				!-----------------------------------------------------------------------
				!       3D prikaz izracunatog magnetskog polja - LabPlot
				!-----------------------------------------------------------------------
				open(7,file='3Dmagn_polje.dat')
				do i = 1,Ntoc
					write(7,'(100e12.4)') (FF(j,i),j=1,Nprav)
				end do
				close(7)		

				! ...............................................
				! Snimanje podataka za graficki prikaz : ParaView
				! ...............................................
				! Priprema koord. osi x i y za crtanje 3D raspodjele elektricnog polja
				if ((paralelan=='y').or.(paralelan=='Y')) then
					allocate(Yos(Ntoc))
					allocate(Xos(Nprav))
					yt = yp1
					do i = 1,Ntoc
						Yos(i) = yt
						yt = yt + dy
					end do
					xt = xp1
					do i = 1,Nprav
						Xos(i) = xt
						xt = xt + razmak
					end do
				else if((paralelan=='x').or.(paralelan=='X')) then
					allocate(Yos(Nprav))
					allocate(Xos(Ntoc))
					xt = xp1
					do i = 1,Ntoc
						Xos(i) = xt
						xt = xt + dx
					end do
					yt = yp1
					do i = 1,Nprav
						Yos(i) = yt
						yt = yt + razmak
					end do
				end if
				!ParaView format
				open(99,file='3DMpolje_para.vtk')
				write(99,'("# vtk DataFile Version 2.0")')
				write(99,'("3D prikaz raspodjele mag. polja po plohi")')
				write(99,'("ASCII")')
				write(99,*)
				write(99,'("DATASET STRUCTURED_GRID")')
				write(99,'("DIMENSIONS",i5,i5,i5)') Ntoc,Nprav,1
				write(99,'("POINTS",i5," float")') Nprav*Ntoc
				do j = 1,Nprav
					do i = 1,Ntoc
						write(99,'(es10.3,2x,es10.3,2x,es10.3)') Xos(i),Yos(j),FF(j,i)
					end do
				end do
				write(99,*)
				write(99,'("POINT_DATA",i5)') Nprav*Ntoc
				write(99,*)
				write(99,'("SCALARS",2x,"epolje"," float"," 1")')
				write(99,'("LOOKUP_TABLE default")')
				do i = 1,Nprav
					do j = 1,Ntoc
						write(99,'(es10.3)') FF(i,j)
					end do
				end do
				close(99)

				deallocate(Xos)
				deallocate(Yos)
				deallocate(FF)
				
				write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele magnetskog polja snimljeni su u file:")')
				write(*,'("3D_mag_polje.dat. Za graficki prikaz koristiti: Gnuplot. Za LabPlot: 3Dmagn_polje.dat")')
				write(*,'("Za ParaView koristiti file: 3DMpolje_para.vtk")')

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==13) THEN
				! ---------------------------------------------------------------------------
				!   PRORACUN RASPODJELE KOMPONENTI MAGNETSKE INDUKCIJE DUZ VISE PARALELNIH 
				!        PRAVACA NA NEKOJ PLOHI (3D RASPODJELA MAGNETSKOG POLJA)
				! ---------------------------------------------------------------------------
				write(*,'(/,"Ukupni broj paralelnih pravaca: ")',advance='no')
				read(*,*) Nprav
				write(*,'("Pravci paralelni koordinatnoj osi x ili y (x/y): ")',advance='no')
				read(*,*) paralelan
				write(*,'("Razmak medju paralelnim pravcima: ")',advance='no')
				read(*,*) razmak
				write(*,'(/,"Koordinate pocetne tocke prvog pravca:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Koordinate krajnje tocke prvog pravca:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				if ((paralelan=='x').or.(paralelan=='X')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dx
					dy = 0.d0
				else if ((paralelan=='y').or.(paralelan=='Y')) then
					write(*,'("Korak proracuna duz pravca: ")',advance='no')
					read(*,*) dy
					dx = 0.d0
				end if
				dz = 0.d0
				write(*,'(/,"Z koordinata plohe na kojoj leze pravci:")')
				write(*,'("   - negativne vrijednosti za plohu u zraku")')
				write(*,'("   - pozitivne vrijednosti za plohu u tlu")')
				write(*,'("   - vrijednost nula za povrsinu tla")')
				write(*,'("-------------------------------------------")')
				write(*,'("Z = ")',advance='no')
				read(*,*) zp1
				if (zp1<0.d0) then !ZRAK
					zp1 = zp1
					zk1 = zp1
					sloj = 1
				else if(zp1>0.d0) then !TLO
					zp1 = zp1
					zk1 = zp1
					write(*,'("Sloj u kojem je ploha: ")',advance='no')
					read(*,*) sloj
				else if(zp1==0.d0) then !POVRSINA TLA
					zp1 = zp1
					zk1 = zp1
					sloj = 2
				end if

				write(*,'(/,"Odaberite sto zelite graficki prikazati:")')
				write(*,'("------------------------------------------------")')
				write(*,'("(1) 3D raspodjela Bx komponente magnetske induk.")')
				write(*,'("(2) 3D raspodjela By komponente magnetske induk.")')
				write(*,'("(3) 3D raspodjela Bz komponente magnetske induk.")')
				write(*,'("(4) 3D raspodjela ukupne/totalne  magnet. induk.")')
				write(*,'("------------------------------------------------")')
				write(*,'("Izbor: ")',advance='no')
				read(*,*) izbor

				! Pocetak mjerenja CPU vremena trajanja proracuna
				call cpu_time(start1)
				ti13 = 1

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if

				! ................................................
				!      Proracun raspodjele magnetske indukcije
				! ................................................
				write(*,'(/,"Racunam ...")')
				allocate(FF(Nprav,Ntoc))
				xtp = xp1
				ytp = yp1
				zt = zp1
				mio = 4.d0*pi*1.e-7
				do j = 1,Nprav
					allocate(Bxx(Ntoc))
					allocate(Byy(Ntoc))
					allocate(Bzz(Ntoc))
					xt = xtp
					yt = ytp
					do i = 1,Ntoc
						! Poziv rutine za proracun Bx, By i Bz komponente magnetskog polja u odabranoj tocki
						! ..........................................................................
						call Hpolje_Bx(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bx)
						call Hpolje_By(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,By)
						call Hpolje_Bz(N,Ns,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,iso,r0,L,Iu,HD,gama,Bz)
						! ..........................................................................
						Bxx(i) = Bx
						Byy(i) = By
						Bzz(i) = Bz
							
						xt = xt + dx
						yt = yt + dy
					end do

					select case(izbor) !Vrijednosti se preracunavaju u mT (militesla)
						case(1)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Bxx(i))**2+dimag(Bxx(i))**2))
								FF(j,i) = modulx * 1.e3 !mT
							end do
							deallocate(Bxx,Byy,Bzz)
						case(2)
							do i = 1,Ntoc
								moduly = dsqrt((dreal(Byy(i))**2+dimag(Byy(i))**2))
								FF(j,i) = moduly * 1.e3 !mT
							end do
							deallocate(Bxx,Byy,Bzz)
						case(3)
							do i = 1,Ntoc
								modulz = dsqrt((dreal(Bzz(i))**2+dimag(Bzz(i))**2))
								FF(j,i) = modulz * 1.e3 !mT
							end do
							deallocate(Bxx,Byy,Bzz)
						case(4)
							do i = 1,Ntoc
								modulx = dsqrt((dreal(Bxx(i))**2+dimag(Bxx(i))**2))
								moduly = dsqrt((dreal(Byy(i))**2+dimag(Byy(i))**2))
								modulz = dsqrt((dreal(Bzz(i))**2+dimag(Bzz(i))**2))
								FF(j,i) = dsqrt(modulx**2+moduly**2+modulz**2) * 1.e3 !mT
							end do
							deallocate(Bxx,Byy,Bzz)
						case default
							write(*,'(/,"Pogresan odabir!!!")')
							goto 100
					end select
					
					if ((paralelan=='x').or.(paralelan=='X')) then
						ytp = ytp + razmak
					else if ((paralelan=='y').or.(paralelan=='Y')) then
						xtp = xtp + razmak
					end if
				end do

				!-----------------------------------------------------------------------
				!       3D prikaz izracunate magnetske indukcije - Gnuplot
				!-----------------------------------------------------------------------
				open(7,file='3D_mag_indukcija.dat')
				do j = 1,Nprav
					do i = 1,Ntoc
						write(7,'(e12.4)') FF(j,i)
					end do
					write(7,*)
				end do
				close(7)		
				!-----------------------------------------------------------------------
				!       3D prikaz izracunate magnetske indukcije - LabPlot
				!-----------------------------------------------------------------------
				open(7,file='3Dmagn_induk.dat')
				do i = 1,Ntoc
					write(7,'(100e12.4)') (FF(j,i),j=1,Nprav)
				end do
				close(7)	
				! ...............................................
				! Snimanje podataka za graficki prikaz : ParaView
				! ...............................................
				! Priprema koord. osi x i y za crtanje 3D raspodjele elektricnog polja
				if ((paralelan=='y').or.(paralelan=='Y')) then
					allocate(Yos(Ntoc))
					allocate(Xos(Nprav))
					yt = yp1
					do i = 1,Ntoc
						Yos(i) = yt
						yt = yt + dy
					end do
					xt = xp1
					do i = 1,Nprav
						Xos(i) = xt
						xt = xt + razmak
					end do
				else if((paralelan=='x').or.(paralelan=='X')) then
					allocate(Yos(Nprav))
					allocate(Xos(Ntoc))
					xt = xp1
					do i = 1,Ntoc
						Xos(i) = xt
						xt = xt + dx
					end do
					yt = yp1
					do i = 1,Nprav
						Yos(i) = yt
						yt = yt + razmak
					end do
				end if
				!ParaView format
				open(99,file='3DMpolje_para.vtk')
				write(99,'("# vtk DataFile Version 2.0")')
				write(99,'("3D prikaz raspodjele mag. polja po plohi")')
				write(99,'("ASCII")')
				write(99,*)
				write(99,'("DATASET STRUCTURED_GRID")')
				write(99,'("DIMENSIONS",i5,i5,i5)') Ntoc,Nprav,1
				write(99,'("POINTS",i5," float")') Nprav*Ntoc
				do j = 1,Nprav
					do i = 1,Ntoc
						write(99,'(es10.3,2x,es10.3,2x,es10.3)') Xos(i),Yos(j),FF(j,i)
					end do
				end do
				write(99,*)
				write(99,'("POINT_DATA",i5)') Nprav*Ntoc
				write(99,*)
				write(99,'("SCALARS",2x,"mpolje"," float"," 1")')
				write(99,'("LOOKUP_TABLE default")')
				do i = 1,Nprav
					do j = 1,Ntoc
						write(99,'(es10.3)') FF(i,j)
					end do
				end do
				close(99)
				deallocate(Xos)
				deallocate(Yos)
				deallocate(FF)
				
				write(*,'(/,"Napomena:",/,"Rezultati proracuna 3D raspodjele mag. indukcije snimljeni su u file:")')
				write(*,'("3D_mag_indukcija.dat. Za graficki prikaz koristiti: Gnuplot. Za LabPlot: 3Dmagn_induk.dat")')

				! Kraj mjerenja CPU vremena trajanja proracuna
				call cpu_time(end1)
				time16 = end1 - start1

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF (izbor==14) THEN
				! ---------------------------------------------------------------------------
				!            PRORACUN NAPONA DODIRA DUZ ODABRANOG PRAVCA 
				! ---------------------------------------------------------------------------
				! Raspodjela potencijala duz pravca na povrsini metala uzemljivackog sustava
				write(*,'(/,"Koordinate pocetne tocke pravca na metalu:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				write(*,'("Zp = ")',advance='no')
				read(*,*) zp1
				write(*,'("Koordinate krajnje tocke pravca na metalu:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				write(*,'("Zk = ")',advance='no')
				read(*,*) zk1
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				write(*,'("Korak duz osi z: ")',advance='no')
				read(*,*) dz
				write(*,'("Sloj u kojem je pravac: ")',advance='no')
				read(*,*) sloj

				write(*,'(/,"Racunam ...")')
				
				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Fi1(Ntoc)) 
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					! Poziv rutine za proracun potencijala u jednoj tocki
					! ..................................................................
					call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
					r0,L,Ip,HD,h,kapa,gama,potencijal)
					! ..................................................................
					Fi1(i) = potencijal
						
					xt = xt + dx
					yt = yt + dy
					zt = zt + dz
				end do
				
				! Raspodjela potencijala duz pravca na povrsini tla udaljenog 1 m od prijasnjeg
				write(*,'(/,"Koordinate pocetne tocke pravca na tlu 1m od prijasnjeg:")')
				write(*,'("Xp = ")',advance='no')
				read(*,*) xp1
				write(*,'("Yp = ")',advance='no')
				read(*,*) yp1
				zp1 = 0.d0
				write(*,'("Koordinate krajnje tocke pravca na tlu 1m od prijasnjeg:")')
				write(*,'("Xk = ")',advance='no')
				read(*,*) xk1
				write(*,'("Yk = ")',advance='no')
				read(*,*) yk1
				zk1 = 0.d0
				write(*,'("Korak duz osi x: ")',advance='no')
				read(*,*) dx
				write(*,'("Korak duz osi y: ")',advance='no')
				read(*,*) dy
				dz = 0.d0
				sloj = 2

				write(*,'(/,"Racunam ...")')

				if (dx/=0.0) then
					Ntoc = (xk1-xp1)/dx
				else if(dy/=0.0) then
					Ntoc = (yk1-yp1)/dy
				else if(dz/=0.0) then
					Ntoc = (zk1-zp1)/dz
				end if
				
				allocate(Fi2(Ntoc)) 
				xt = xp1
				yt = yp1
				zt = zp1
				do i = 1,Ntoc
					! Poziv rutine za proracun potencijala u jednoj tocki
					! ..................................................................
					call potencijal_tocke_freq(N,Ns,Nsi,sloj,xt,yt,zt,xp,yp,zp,xk,yk,zk,Fr,iso,&
					r0,L,Ip,HD,h,kapa,gama,potencijal)
					! ..................................................................
					Fi2(i) = potencijal
						
					xt = xt + dx
					yt = yt + dy
				end do
				
				! ................................................
				! PlotXY - graficki prikaz napona dodira
				! ................................................
				open(1,file='Napon_dodira.adf')  
				write(1,*)
				write(1,'("x/y",15x,"Nap_dodira")')
				if (dx/=0.0) then
						t = xp1
				else if(dy/=0.0) then
						t = yp1
				else if(dz/=0.0) then
						t = zp1
				end if
				do i = 1,Ntoc
					modulx = dsqrt((dreal(Fi1(i))**2+dimag(Fi1(i))**2))	!Recikliranje varijable
					moduly = dsqrt((dreal(Fi2(i))**2+dimag(Fi2(i))**2))	!Recikliranje varijeble
					modul = modulx - moduly	!Napon dodira (recikliranje varijable)

					write(1,'(f10.3,f16.6)') t,modul
					if (dx/=0.0) then
						t = t + dx
					else if(dy/=0.0) then
						t = t + dy
					else if(dz/=0.0) then
						t = t + dz
					end if
				end do
				close(1)
				deallocate(Fi1,Fi2)

				write(*,'(/,"Napomena:",/,9("-"),/,"Rezultati za graficki prikaz snimljeni su u file:")')
				write(*,'(4x,"Napon_dodira.adf      - graf napona dodira duz odabranog pravca")')

				! Jos proracuna na istom uzemljivacu
				write(*,'(/,"Da li zelite napraviti jos neki proracun (y/n): ")',advance='no')
				read(*,*) izb
				if ((izb=='Y').or.(izb=='y')) then
					ponovi = .true.
				else if ((izb=='N').or.(izb=='n')) then
					ponovi = .false.
				end if

			ELSE IF(izbor==0) THEN
				ponovi = .false.

			ELSE
				write(*,'(/,"Pogresan odabir - kraj programa!")')
				ponovi = .false.
			END IF
		END DO

		100 continue

		! Zatvaranje ulazno/izlaznog filea
		close(2)

		! Brisanje privremenog filea
		!res = DELFILESQQ('glob_pot.dat')

		! Oslobadjanje zauzete memorije
		deallocate(Iu,Ip)
		deallocate(L,r0)
		deallocate(xp,yp,zp)
		deallocate(xk,yk,zk)
		deallocate(kapa,gama)
		deallocate(HD,h)
		deallocate(Fr,iso)
	
	END IF

	! ---------------------------------------------------------------------------------
	!                  STATISTIKA CPU VREMENA TRAJANJA PRORACUNA
	! ---------------------------------------------------------------------------------
	write(*,'(//,"Statistika CPU vremena trajanja proracuna po fazama - [sekunde]:")')
	write(*,'(65("-"))')
	if (tip_pr==2) then
		write(*,'("1) Priprema, formiranje i rjesavanje globalnog sustava jednadzbi:",f8.2)') time1
		if (tt1 == 1) &
			write(*,'("2) Proracun tranzijentnog potenci. u jednoj tocki (FFT-MoM-IFFT):",f8.2)') time2
		if (tt2 == 1) &
			write(*,'("2) Proracun tranzijentnog  potencijala duz pravca (FFT-MoM-IFFT):",f8.2)') time22
		if (tt3 == 1) &
			write(*,'("2) Proracun tranzij. elektr. polja u jednoj tocki (FFT-MoM-IFFT):",f8.2)') time33
		if (tt4 == 1) &
			write(*,'("2) Proracun tranzij. magnet. polja u jednoj tocki (FFT-MoM-IFFT):",f8.2)') time44
		if (tt5 == 1) &
			write(*,'("2) Proracun 3D tranzij. magnet. polja  duz pravca (FFT-MoM-IFFT):",f8.2)') time55
		if (tt6 == 1) &
			write(*,'("2) Proracun 3D tranzij. elektr. polja  duz pravca (FFT-MoM-IFFT):",f8.2)') time66
		if (tt7==1) &
			write(*,'("2) Proracun  tranzijentne ulazne impedanc. uzemljivackog sustava:",f8.2)') time77			
 		
		time_uk = time1 + time2 + time22 + time33 + time44 + time55 + time66 + time77

	else if(tip_pr==1) then
		write(*,'("1) Priprema, formiranje i rjesavanje  globalnog sustava jednadzbi")')
		write(*,'("   te raspodjela struja segmenata uzemlj. (uzd. i poprecne str.):",f8.2)') time3
		if (ti1 == 1) &
			write(*,'("2) Proracun potencijala u samo jednoj tocki viseslojnog sredstva:",f8.2)') time9
		if (ti2 == 1) &
			write(*,'("2) Proracun potencijala  duz odabranog pravca u visesl. sredstvu:",f8.2)') time4
		if (ti3 == 1) &
			write(*,'("2) Proracun  potencijala duz  odabranog pravca na  povrsini  tla:",f8.2)') time5
		if (ti4 == 1) &
			write(*,'("2) Proracun napona koraka duz nekog proizvoljno odabranog pravca:",f8.2)') time6
		if (ti5 == 1) &
			write(*,'("2) Proracun raspodjele potenc. duz vise pravaca istovremeno (3D):",f8.2)') time7
		if (ti6 == 1) &
			write(*,'("2) Proracun  napona dodira  u jednoj proizvoljno odabranoj tocki:",f8.2)') time10
		if (ti7 == 1) &
			write(*,'("2) Proracun napona koraka duz vise paralel.  pravaca istovremeno:",f8.2)') time8
		if (ti8 == 1) &
			write(*,'("2) Proracun elektricnog polja  u jednoj tocki visesloj. sredstva:",f8.2)') time11
		if (ti9 == 1) &
			write(*,'("2) Proracun elektricnog  polja  duz pravca u viseslojn. sredstvu:",f8.2)') time12
		if (ti10 == 1) &
			write(*,'("2) Proracun 3D raspodjele elektricnog  polja po bilo kojoj plohi:",f8.2)') time13
		if (ti11 == 1) &
			write(*,'("2) Proracun magnetskog  polja  u jednoj tocki visesloj. sredstva:",f8.2)') time14
		if (ti12 == 1) &
			write(*,'("2) Proracun 3D raspodjele magnetskog  polja po  bilo kojoj plohi:",f8.2)') time15
		if (ti13 == 1) &
			write(*,'("2) Proracun 3D raspodjele  magnet. indukcije po bilo kojoj plohi:",f8.2)') time16

		time_uk = time3 + time4 + time5 + time6 + time7 + time8 + &
				  time9 + time10 + time11 + time12 + time13 + time14 + time15 + time16
	end if
	write(*,'(65("-"))')
	write(*,'(3x,"Ukupno CPU vrijeme trajanja proracuna:",f10.2," [sec]")') time_uk

	! Kraj programa
	write(*,'(//,"*** KRAJ PROGRAMA! ***")')
!	write(*,'(/,"Pritisnite Enter za izlaz.")')
!	read(*,*)
	
	stop
END PROGRAM
