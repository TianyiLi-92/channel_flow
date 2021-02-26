program  post_poiseuilleflow
    implicit none

    include 'mpif.h'

    !grid points nnx*nny*(nnz+1)
    integer, parameter :: nproc=64      ! the number of the processors
    integer, parameter :: nnx=128       !dimensions; NOTE:nny divided by 2 must be multiplies of nproc!!!
    integer, parameter :: nny=128
    integer, parameter :: nnz=128
    integer, parameter :: nx=nnx/nproc-1
    integer, parameter :: ny=nny-1
    integer, parameter :: nz=nnz

    real*8, parameter :: aphi=2.0d0,beta=1.0d0    ! the wave number of x and y direction

    integer :: istep,i,j,k,m

    character*70 outname1
    real*8,dimension(0:nx,0:ny,0:nz,3)   ::  ftmp
    real*8,dimension(0:nx,0:ny,0:nz,3)   ::  ftmp1
    real*8,dimension(0:nnx-1,0:ny,0:nz,3)   ::  ftmp2

    character(70)   ::  fnm1

    integer id,ierr,p,nallgrp,ids,stat(mpi_status_size)

    real*8  ::  ua(0:nz)

    !//////////////////////////////////////////////////////////////
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,id,ierr)
    call mpi_comm_size(mpi_comm_world,p,ierr)
    call mpi_barrier(mpi_comm_world,ierr)
    nallgrp = mpi_comm_world

    if(id==0) write(*,*)"MPI Init complete."
    !//////////////////////////////////change////////////////////////////


    !if(id==0) open(999,file='images/us_y41',form='unformatted',position='append')
    !if(id==0) open(999,file='images/vel_y41',form='unformatted',position='append')
    !if(id==0) open(999,file='images/us_ver',form='unformatted',position='append')
    if(id==0) open(996,file='images/vel_y0.05391',form='unformatted',position='append')
    if(id==0) open(997,file='images/vel_y10.46',form='unformatted',position='append')
    if(id==0) open(998,file='images/vel_y30.17',form='unformatted',position='append')
    if(id==0) open(999,file='images/vel_y49.36',form='unformatted',position='append')

    !do istep=3000000,8000000,200
    do istep=500200,1200000,200
    !do istep=7204200,7204600,200
        write(outname1,212)istep,id
212         format('../3Ddata/velocity',i8.8,'.',i3.3,'.dat')
        open(12,file=outname1,status='old',form='unformatted')
        read(12) ((((ftmp(i,j,k,m),i=0,nx),j=0,ny),k=0,nz),m=1,3)
        close(12)

        if(id.ne.0) call mpi_send(ftmp,(nx+1)*(ny+1)*(nz+1)*3,mpi_real,0,id,nallgrp,ierr)
        if(id.eq.0)then
            ftmp2(0:nx,:,:,:)=ftmp
            do ids=1,nproc-1
                call mpi_recv(ftmp1,(nx+1)*(ny+1)*(nz+1)*3,mpi_real,ids,ids,nallgrp,stat,ierr)
                ftmp2(ids*(nx+1):ids*(nx+1)+nx,:,:,:)=ftmp1
            enddo
        
            !write(999) ((ftmp2(i,j,28,1),i=0,nnx-1),j=0,ny)
            !write(999) ((ftmp2(i,j,100,1),i=0,nnx-1),j=0,ny)
            !write(999) (((ftmp2(i,j,28,m),i=0,nnx-1),j=0,ny),m=1,3)
            !write(999) (((ftmp2(i,j,100,m),i=0,nnx-1),j=0,ny),m=1,3)
            !write(999) ((ftmp2(0,j,k,1),j=0,ny),k=0,nz)
            !write(999) ((ftmp2(64,j,k,1),j=0,ny),k=0,nz)
            write(996) (((ftmp2(i,j,1,m),i=0,nnx-1),j=0,ny),m=1,3)
            write(996) (((ftmp2(i,j,127,m),i=0,nnx-1),j=0,ny),m=1,3)

            write(997) (((ftmp2(i,j,14,m),i=0,nnx-1),j=0,ny),m=1,3)
            write(997) (((ftmp2(i,j,114,m),i=0,nnx-1),j=0,ny),m=1,3)

            write(998) (((ftmp2(i,j,24,m),i=0,nnx-1),j=0,ny),m=1,3)
            write(998) (((ftmp2(i,j,104,m),i=0,nnx-1),j=0,ny),m=1,3)

            write(999) (((ftmp2(i,j,31,m),i=0,nnx-1),j=0,ny),m=1,3)
            write(999) (((ftmp2(i,j,97,m),i=0,nnx-1),j=0,ny),m=1,3)

!             ua=0.0d0
!             do k=0,nz
!                 do j=0,ny
!                     do i=0,nnx-1
!                         ua(k)=ua(k)+ftmp2(i,j,k,1)
!                     enddo
!                 enddo
!             enddo
!             ua=ua/(nnx*(ny+1))
!             write(outname1,213)istep
! 213         format('meanu/average_velocity',i8.8,'.','plt')
!             open(12,file=outname1,status='replace')
!             do k=0,nz
!                 write(12,*)k,ua(k)
!             enddo
!             close(12)
        endif
    end do

    if(id==0) close(996)
    if(id==0) close(997)
    if(id==0) close(998)
    if(id==0) close(999)

    call mpi_finalize(ierr)
    stop
end program  post_poiseuilleflow