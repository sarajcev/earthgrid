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

! Ispis naslova programskog paketa i osnovnih podataka o nacinu formiranja
! datoteke s ulaznim podacima, kao i ispis GPL licence.

subroutine naslov(ime_file)
	implicit none
	character(128) ime_file
	
	open(2,file='output.txt')
	write(2,'(/,tr10,"EartHGriD Linux Ver. 1.0",/)')
	write(2,'("PROGRAM ZA  NUMERICKI PRORACUN  HARMONICKOG POLJA  SUSTAVA")')
	write(2,'("UZEMLJIVACA POLOZENIH U  VISESLOJNO SREDSTVO - STACIONARNA")')
	write(2,'(8X,"I TRANZIJENTNA (FFT - MoM - IFFT) ANALIZA")')
	write(2,'(/,"Copyright (C) 2006. GPL - General Public Licence")')
	write(2,'("Author: Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)")')
	write(2,*)
	write(2,'("----------------------------------------------------------")')
	write(2,'(13x,"GNU General Public Licence:")')
	write(2,'("----------------------------------------------------------")')
    write(2,'("This program is free software; you can redistribute it and")')
    write(2,'("/or modify it under the terms of the  GNU  General  Public")')
    write(2,'("License (GPL) as published by the Free Software Foundation")')
    write(2,'("either version 2  of the License,  or (at your option) any")')
    write(2,'("later version.")')
	write(2,*)
    write(2,'("This program is  distributed in the  hope that it will  be")')
    write(2,'("useful, but WITHOUT ANY WARRANTY; without even the implied")')
    write(2,'("warranty of  MERCHANTABILITY  or  FITNESS FOR A PARTICULAR")')
    write(2,'("PURPOSE. See the attached  GNU General Public License  for")')
    write(2,'("more details.")')
	write(2,*)
    write(2,'("You should have received a copy of the GNU  General Public")')
    write(2,'("License along with this program; if not, write to the Free")')
    write(2,'("Software Foundation, Inc., 51 Franklin Street, Fifth Floor")')
	write(2,'(",Boston, MA02110-1301, USA.")')
	write(2,'("----------------------------------------------------------")')
	
	! ---------------------------------------------------------------------------
	!       CITANJE ULAZNIH PODATAKA IZ VANJSKOG FILE-A (INPUT)
	!       I TISKANJE ULAZNIH PODATAKA U DRUGI VANJSKI FILE (OUTPUT)
	! ---------------------------------------------------------------------------
	write(*,'(/,"----------------------------------------------------------")')		
	write(*,'("Slijedite jednu od ponudjenih  opcija programa: EartHGrid:")')
	write(*,'("----------------------------------------------------------")')	
	write(*,'("OPCIJA (1):")')
	write(*,'("Unesite ime file-a  koji sadrzi ulazne podatke za proracun")')
	write(*,'("(ukljucujuci i extenziju - .txt, .dat).  Doticni file mora")') 
	write(*,'("biti formatiran na propisan nacin  (ovo se ne provjerava).")')
	write(*,'("Automatski se unosom valjanog imena filea starta proracun.")')
	write(*,'("OPCIJA (2):")')
	write(*,'("Za pomoc o nacinu formatiranja  file-a s ulaznim podacima,")')
	write(*,'("upisite: help.  Prikazuje se potom forma i objasnjenje na-")')
	write(*,'("cina zadavanja  (formatiranja) file-a  s ulaznim podacima.")')
	write(*,'("OPCIJA (3):")')
	write(*,'("Za prikaz GNU  General Public Licence upisite: LIC.  Nakon")')
	write(*,'("ispisa teksta licence  (na engleskom) program  se prekida.")')
	write(*,'("OPCIJA (4):")')
	write(*,'("Za izlaz iz programa  bez ikakvih proracuna upisite: exit.")')
	write(*,'("----------------------------------------------------------")')	
	write(*,'("EartHGriD:> ")',advance='no')
	read(*,*) ime_file
	ime_file = trim(ime_file)
	
	if ((ime_file == 'exit').or.(ime_file == 'EXIT')) then
		write(*,'(/,"*** KRAJ PROGRAMA! ***")')	
		stop !Zavrsetak programa
	else if ((ime_file == 'lic').or.(ime_file == 'LIC')) then
		write(*,*)
		write(*,'("----------------------------------------------------------")')
		write(*,'(13x,"GNU General Public Licence:")')
		write(*,'("----------------------------------------------------------")')
		write(*,'("This program is free software; you can redistribute it and")')
		write(*,'("/or modify it under the terms of the  GNU  General  Public")')
		write(*,'("License (GPL) as published by the Free Software Foundation")')
		write(*,'("either version 2  of the License,  or (at your option) any")')
		write(*,'("later version.")')
		write(*,*)
		write(*,'("This program is  distributed in the  hope that it will  be")')
		write(*,'("useful, but WITHOUT ANY WARRANTY; without even the implied")')
		write(*,'("warranty of  MERCHANTABILITY  or  FITNESS FOR A PARTICULAR")')
		write(*,'("PURPOSE. See the attached  GNU General Public License  for")')
		write(*,'("more details.")')
		write(*,*)
		write(*,'("You should have received a copy of the GNU  General Public")')
		write(*,'("License along with this program; if not, write to the Free")')
		write(*,'("Software Foundation, Inc., 51 Franklin Street, Fifth Floor")')
		write(*,'(",Boston, MA02110-1301, USA.")')
		write(*,'("----------------------------------------------------------")')
		stop !Zavrsetak programa
	else if ((ime_file=='help').or.(ime_file=='HELP')) then
		write(*,*)
		write(*,'("U nastavku je prikazan jedan primjer forme ulazne datoteke")')
		write(*,'("na kojem je objasnjen nacin formatiranja ulaznih podataka.")')
		write(*,'("----------------------------------------------------------")')
		write(*,'("UlazniPodaci.txt")')
		write(*,'("----------------------------------------------------------")')
		write(*,'("4                     | Ukupni broj slojeva ukljuc. i zrak")')
		write(*,'("0.0 1.0 0.0           | Parametri zraka (uvijek: 0. 1. 0.)")')
		write(*,'("ro2 epsr2 h2          | specif. relativna otporn. (ohm*m),")')
		write(*,'("ro3 epsr3 h3          | specif. permitivno. i dubina sloj-")')
		write(*,'("ro4 epsr4 0.0         | eva tla-zadnji sloj je beskon.: 0.")')
		write(*,'("2                     | Ukupni broj tipova vodica  uzemlj.")')
		write(*,'("1 r0 mirv sigmavv     | Radijus (m), relat. specif. perme-")')
		write(*,'("2 r0 mirv sigmavv     | abilnost i vodljivost (S/m) vodica")')
		write(*,'("5 7 1 1               | Ukupni broj vodica i segmenata.***")')
		write(*,'("Xp Yp Zp Xk Yk Zk 1 2 | X, y i z koordinate pocetne tocke,")')
		write(*,'("Xp Yp Zp Xk Yk Zk 1 2 | te x, y i z koordinate krajnje to-")')
		write(*,'("Xp Yp Zp Xk Yk Zk 2 1 | cke svakog vodica (m), tip  vodica")')
		write(*,'("Xp Yp Zp Xk Yk Zk 2 1 | s obzirom na presjek, te broj pod-")')
		write(*,'("Xp Yp ZP XK YK ZK 1 1 | jela tog  vodica  na segmente ****")')
		write(*,'("tip_pr                | Tip proracuna koji ce se izvrsiti*")')
		write(*,'("Ni=1  f               | Broj narinutih str. i  frekv. str.")')
		write(*,'("xi yi zi Re Im        | Koord. narinute struje, te iznos**")')
		write(*,'("----------------------------------------------------------")')
		write(*,'(" * - na mjestu tip_pr dozvoljeno je upisati jednu od slje-")')
		write(*,'("     dece dvije vrijednosti: 1 ili 2. Ukoliko se  upise: 1")') 
		write(*,'("     rijec je o stacionar. proracunu uzemljivaca  za jednu")')
		write(*,'("     zadanu frekvenciju. Ako se upise: 2 rijec je o tranz-")')
		write(*,'("     ijentnoj analizi uzem. Pritom Ni=1 - mora biti jedan!")')
		write(*,'("** - zadaju se x, y i z koord. mjesta (m) gdje je narinuta")') 
		write(*,'("     struja (ili vise njih ako je rijec o stacion. prorac.")') 
		write(*,'("     biti ce onoliko redaka koliko je zadanih struja. Str-")')
		write(*,'("     uja se zadaje realnim i imaginarnim dijelom (Re, Im).")')
		write(*,'("*** -zadaju se respektivno ukupni broj svih vodica i segm-")')
		write(*,'("     enata sustav uzem., te potom broj izoliranih vodica i")')
		write(*,'("     segmenata.  Ake isti ne  postoje  postavlja se: 0 0 !")')
		write(*,'("****-prvo se zadaju geometrijski podaci svih  neizoliranih")')
		write(*,'("     (golih) segmenata, a potom se na kraju zadaju geomet-")')
		write(*,'("     rijski podaci izoliranih segmenata. Ovu  proceduru je")')
		write(*,'("     nuzno postivati zbog ispravnog proracuna!!!")')
		write(*,'("Napomena: Svi podaci u ulaznoj datoteci zadaju se u slobo-")')
		write(*,'("dnom formatu na gore propisan nacin. Ime datoteke je proi-")')
		write(*,'("zvoljno, kao i extenzija.  Ime ne smije sadrzavati razmak!")')
		write(*,*)
		stop !Zavrsetak programa
	end if

	return
end subroutine