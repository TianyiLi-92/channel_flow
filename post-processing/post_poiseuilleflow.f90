program  post_poiseuilleflow
    use tecplot

    implicit none

    include 'mpif.h'

    !grid points nnx*nny*(nnz+1)
    integer, parameter :: nproc=64      ! the number of the processors
    integer, parameter :: nnx=128       !dimensions; NOTE:nny divided by 2 must be multiflies of nproc!!!
    integer, parameter :: nny=128
    integer, parameter :: nnz=128
    integer, parameter :: nx=nnx/nproc-1
    integer, parameter :: ny=nny-1
    integer, parameter :: nz=nnz

    real*8, parameter :: aphi=2.0d0,beta=1.0d0    ! the wave number of x and y direction

    integer :: istep,i,j,k,m

    real*4  ::  pi,x(0:nnx-1,0:ny,0:nz),y(0:nnx-1,0:ny,0:nz),z(0:nnx-1,0:ny,0:nz)

    character*70 outname1
    real*8,dimension(0:nx,0:ny,0:nz,3)   ::  ftmp
    real*8,dimension(0:nx,0:ny,0:nz,3)   ::  ftmp1
    real*8,dimension(0:nnx-1,0:ny,0:nz,3)   ::  ftmp2

    character(70)   ::  fnm1
    type(tecplot_time_file) ::  plt_file
    real*4,dimension(0:nnx-1,0:ny,0:nz,6) ::  field
    integer ::  locations(6), type_list(6), shared_list(6)

    real*8,dimension(0:nz/2)    ::  tmp_yp,tmp_up,yp,up
    real*8  ::  tmp
    real*8,dimension(0:nz)  ::  Rss,vs,tss,tmp_Rss,tmp_vs,tmp_tss
    real*8,dimension(0:nz)  ::  urms,vrms,wrms,tmp_urms,tmp_vrms,tmp_wrms

    integer id,ierr,p,nallgrp,ids,stat(mpi_status_size)

    !//////////////////////////////////////////////////////////////
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,id,ierr)
    call mpi_comm_size(mpi_comm_world,p,ierr)
    call mpi_barrier(mpi_comm_world,ierr)
    nallgrp = mpi_comm_world

    if(id==0) write(*,*)"MPI Init complete."
    !//////////////////////////////////change////////////////////////////


    pi=4.0*atan(1.0)
    do i=0,nnx-1
        x(i,:,:)=2.0*pi*real(i)/real(nnx)*aphi
    enddo

    do j=0,ny
        y(:,j,:)=2.0*pi*real(j)/real(ny+1)*beta
    enddo

    do k=0,nz
        z(:,:,k)=cos(pi*real(k)/real(nz))
    enddo

    do istep=300000,300000
        write(outname1,212)istep,id
212         format('../3Ddata/velocity',i8.8,'.',i3.3,'.dat')
        open(12,file=outname1,status='old',form='unformatted')
        read(12) ((((ftmp(i,j,k,m),i=0,nx),j=0,ny),k=0,nz),m=1,3)
        close(12)

        if(id.ne.0) call mpi_send(ftmp,(nx+1)*(ny+1)*(nz+1)*3,MPI_DOUBLE_PRECISION,0,id,nallgrp,ierr)
        if(id.eq.0)then
            ftmp2(0:nx,:,:,:)=ftmp
            do ids=1,nproc-1
                call mpi_recv(ftmp1,(nx+1)*(ny+1)*(nz+1)*3,MPI_DOUBLE_PRECISION,ids,ids,nallgrp,stat,ierr)
                ftmp2(ids*(nx+1):ids*(nx+1)+nx,:,:,:)=ftmp1
            enddo
        endif

        if(id==0)then
            field(:,:,:,1)=x
            field(:,:,:,2)=y
            field(:,:,:,3)=z
            field(:,:,:,4:6)=sngl(ftmp2)

            locations = 0
            type_list = 1       !  1 = float
            shared_list = -1    ! -1 = not shared

            write(fnm1,213)istep
213         format('3Ddata_p/velocity',i8.8,'.plt')
            call plt_file%init(fnm1,nnx,ny+1,nz+1,'file title','x,y,z,u,v,w')
            call plt_file%write_zone_header('zone title',real(istep,kind=4),0,locations)
            call plt_file%write_zone_data(type_list,shared_list,field)
            call plt_file%complete
        endif
    end do

!     yp=0.0
!     up=0.0
!     Rss=0.0
!     vs=0.0
!     tss=0.0
!     urms=0.0
!     vrms=0.0
!     wrms=0.0
!     do istep=200000,400000,1000
!         write(outname1,21116)istep
! 21116   format('../poiseuille_flow3/uplus/uplus',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz/2
!             read(12,*)tmp_yp(k),tmp_up(k)
!         end do
!         close(12)
!         yp=yp+1.0/tmp_yp
!         up=up+tmp_up

!         write(outname1,21113)istep
! 21113   format('../poiseuille_flow3/stress/Reynolds_shear_stress',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_Rss(k)
!         enddo
!         close(12)
!         Rss=Rss+tmp_Rss

!         write(outname1,21114)istep
! 21114   format('../poiseuille_flow3/stress/viscous_stress',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_vs(k)
!         enddo
!         close(12)
!         vs=vs+tmp_vs

!         write(outname1,2118)istep
! 2118    format('../poiseuille_flow3/stress/total_shear_stress',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_tss(k)
!         enddo
!         close(12)
!         tss=tss+tmp_tss

!         write(outname1,21119)istep
! 21119   format('../poiseuille_flow3/rms/urms',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_urms(k)
!         enddo
!         close(12)
!         urms=urms+tmp_urms

!         write(outname1,2115)istep
! 2115    format('../poiseuille_flow3/rms/vrms',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_vrms(k)
!         enddo
!         close(12)
!         vrms=vrms+tmp_vrms

!         write(outname1,2116)istep
! 2116    format('../poiseuille_flow3/rms/wrms',i6.6,'.','plt')
!         open(12,file=outname1,status='old')
!         do k=0,nz
!             read(12,*)tmp,tmp_wrms(k)
!         enddo
!         close(12)
!         wrms=wrms+tmp_wrms
!     end do
!     yp=201.0/yp
!     up=up/201.0
!     open(12,file='average/uplus.plt')
!     do k=0,nz/2
!         write(12,*)yp(k),up(k)
!     end do
!     close(12)

!     Rss=Rss/201.0
!     open(12,file='average/Reynolds_shear_stress.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),Rss(k)
!     enddo
!     close(12)

!     vs=vs/201.0
!     open(12,file='average/viscous_stress.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),vs(k)
!     enddo
!     close(12)

!     tss=tss/201.0
!     open(12,file='average/total_shear_stress.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),tss(k)
!     enddo
!     close(12)

!     urms=urms/201.0
!     open(12,file='average/urms.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),urms(k)
!     enddo
!     close(12)

!     vrms=vrms/201.0
!     open(12,file='average/vrms.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),vrms(k)
!     enddo
!     close(12)

!     wrms=wrms/201.0
!     open(12,file='average/wrms.plt')
!     do k=0,nz
!         write(12,*)z(0,0,k),wrms(k)
!     enddo
!     close(12)

    call mpi_finalize(ierr)
    stop
end program  post_poiseuilleflow