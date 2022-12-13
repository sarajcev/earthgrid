!************************************************************************************
!                            PRORACUN Y(s)
!************************************************************************************
subroutine Y_single_frequency(N,ro,epsr,Ns,Nsi,Nc,Ni,Ii,mirv,sigmav,r0,f,L,HD,h,iso,&
    sloj,icv,iveze,xp,yp,zp,xk,yk,zk,xt,yt,zt,potenc_freq)
    use funkcije
    implicit none

    ! Input variables:
    integer N,Ns,Nsi,Nc,Ni,sloj
    integer,dimension(:) :: iveze,icv,iso
    real(8),dimension(:) :: ro,epsr,mirv,sigmav,L,HD,h
    real(8),dimension(:) :: xp,yp,zp,xk,yk,zk,r0
    real(8) xt,yt,zt,f
    complex(8),dimension(:) :: Ii
    ! Output variables:
    complex(8) potenc_freq
    ! Local variables:
    integer Nss
    complex(8),dimension(:),allocatable :: kapa,gama,Fr,zunjv,Fig,Fi,Iu,Ip
    complex(8),dimension(:,:),allocatable :: Yuv,Ypv,Yu,Ypp,Yg,B,Au,Ap,CC
    integer,dimension(:),allocatable :: ipiv
    complex(8) alpha,betha,potencijal
    character(1) transa,transb
    integer nrhs,info,ig,jg,i,j

    integer row_index
    complex(8),dimension(Nc-1) :: Y_2,Y_3
    complex(8),dimension(Nc-1,Nc-1) :: Y_4
    complex(8),dimension(Nc-1) :: ABC
    complex(8),dimension(Nc) :: swap
    integer NB,Lwork
    complex(8),dimension(:),allocatable :: Work
    complex(8) suma

    external zgetri

    ! ---------------------------------------------------------------------------
    !            PRORACUN FREKVENCIJSKI OVISNIH PARAMETARA PROGRAMA
    ! ---------------------------------------------------------------------------
    Nss = Ns - Nsi

    ! Proracun kompleksnih specificnih vodljivosti slojeva modela
    allocate(kapa(N))
    call proracun_kapa(N,f,ro,epsr,kapa)
    ! Proracun valnih konstanti slojeva modela
    allocate(gama(N))
    call proracun_gama(f,N,kapa,gama)
    ! Proracun faktora refleksije svih slojeva modela
    allocate(Fr(N))
    call vektor_F(N,kapa,Fr)
    ! Proracun jedinicnih unutarnjih impedancija svih segmenata
    allocate(zunjv(Ns))
    call proracun_zunjv(Ns,f,mirv,sigmav,r0,zunjv)

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
    deallocate(zunjv)

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
    !         Proracun redukcije matrice [Yg]
    ! ---------------------------------------------------------------------------
    outer: do j = 1,Ni
        inner: do i = 1,Nc
            if (icv(j)==i) then
                row_index = i    ! indeks retka i stupca
                exit outer
            end if
        end do inner
    end do outer
    !write(*,'("row_index =",i4)') row_index
    if (row_index /= 1) then
        ! PO RETCIMA
        swap = Yg(1,:)
        Yg(1,:) = Yg(row_index,:)
        Yg(row_index,:) = swap
        ! PO STUPCIMA
        swap = Yg(:,1)
        Yg(:,1) = Yg(:,row_index)
        Yg(:,row_index) = swap
    end if

    !Particioniranje globalne matrice
    Y_2(:) = Yg(1,2:Nc)
    Y_3(:) = Yg(2:Nc,1)
    Y_4(:,:) = Yg(2:Nc,2:Nc)

    !Invertiranje matrice Y_4
    allocate(ipiv(Nc-1))
    call zgetrf(Nc-1,Nc-1,Y_4,Nc-1,ipiv,info)
    NB = ILAENV( 1, 'ZGETRI', ' ', Nc-1, -1, -1, -1 )
    Lwork = (Nc-1)*NB
    allocate(work(Lwork))
    call zgetri(Nc-1,Y_4,Nc-1,ipiv,work,Lwork,info)
    deallocate(ipiv,work)

    !Proracun redukcije
    ABC = matmul(Y_4,Y_3)
    suma = cmplx(0.d0,0.d0)
    do i = 1,Nc-1
        suma = suma + Y_2(i)*ABC(i)
    end do

    potenc_freq = Yg(1,1) - suma

    deallocate(Yg)
    deallocate(kapa,gama)
    deallocate(Fr)

    return
end subroutine
