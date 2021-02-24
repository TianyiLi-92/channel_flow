    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    ! Pay attention to that the sequence of velocity and vorticty .            !
    ! This program is designed to compute channel flow.                        !
    ! u,v,w are the velocities;omegax,omegay,omegaz is the vorticities.        !
    ! This program employs cheybshev-tau method: in x,y directions, we use     !
    ! fourier series, in z direcition cheybshev polynomials are used.          !
    ! The numerical method is like what Luo used in his PC code(1998).         !
    ! please notice that in order to save the memory,in regard to the fourier  !
    ! coiefience conjugate,we use in the x direcition 0:(nx-1)/2 to save       !
    ! the real*8 part, and (nx+1)/2:nx to save the image part of the variable. !
    ! Give zero to the (nx+1)/2 mode,for it is small in general.               !
    ! Also, we must know that we use u,w to store the nolinear part of the     !
    ! equations. unp1 and wnp1, in computing v and omegay, represent the n-1   !
    ! time step nonlinear part; rhsvn and rhsomegayn represent the n time step !
    ! nonlinear part. "Single precision".(Oct 08, 2003)                        !
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    use paraandcons
    IMPLICIT NONE
    
    include 'mpif.h'
    include "fftw_f77.i"
    include 'mkl.fi'
    !IMPLICIT NONE
    !include 'channel.inc'
    !include 'channel1.inc'
    !include 'channel3.inc'
    !integer, parameter :: isign = 1
    INTEGER,parameter :: isign=-1 !-1:given by possuie flow; 1:given by file
    !real*8,dimension(0:nhx,0:ny) :: fr
    REAL(kind=8) ::  fr(0:nhx,0:ny)
    REAL(kind=8) ::  fi(0:nhx,0:ny)

    REAL(kind=8) ::  fr1(0:nx,0:ny)
    REAL(kind=8) ::  fi1(0:nx,0:ny)

    REAL(kind=8) ::  fr2(0:nx,0:ny)
    REAL(kind=8) ::  fi2(0:nx,0:ny)

    REAL(kind=8) ::  tvn(0:nz)
    REAL(kind=8) ::  tvni(0:nz)

    REAL(kind=8) ::  tvnm1(0:nz)
    REAL(kind=8) ::  tvnm1i(0:nz)

    REAL(kind=8) ::  tomegayn(0:nz)
    REAL(kind=8) ::  tomegayni(0:nz)


    REAL(kind=8) ::  x(0:nx)
    REAL(kind=8) ::  y(0:ny)
    REAL(kind=8) ::  z(0:nz)          ! the location of (x,y,z)

    REAL(kind=8) ::  rhsr(0:nz),rhsi(0:nz) ! the right hand of equations

    integer i,j,k,istep,kk,ka

    real*8 ibeta(0:ny),iapha(0:nhx) ! the tranform coeefience of x mode and y mode
    complex(kind=8) rhs(0:nz),a(0:nz,0:nz)    ! the coef matrix
    complex(kind=8) coef_e(0:nz,0:nz,0:nhx,0:ny)
    complex(kind=8) coef_a(0:nz,0:nz,0:nhx,0:ny)
    complex(kind=8) b(0:nz,0:nz)

    real*8 kmn,coef,pi,run_time
    real*8 bdv(4),bdomega(2)
    real*8 starttime,endtime,etime1,etime2,etime3,etime4
    real*8 s_eq_time,e_eq_time,t_eq_time
    real*8 bdvi(4),bdomegai(2),bdu(2),bdui(2)
    integer id,ierr,p,nallgrp,nstat
    character*50 inname

    integer :: ipiv(nz+1,0:nhx,0:ny), info
    complex(kind=8) work((nz+1)*64),mytemp(nz+1)

    !//////////////////////////////////////////////////////////////
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,id,ierr)
    call mpi_comm_size(mpi_comm_world,p,ierr)
    call mpi_barrier(mpi_comm_world,ierr)
    nallgrp = mpi_comm_world


    if(id==0) write(*,*)"MPI Init complete."
    !//////////////////////////////////change////////////////////////////

    allocate(t0(0:nz,0:nz),t2(0:nz,0:nz),t4(0:nz,0:nz),t1(0:nz,0:nz))
    allocate(u0(0:nz),du0(0:nz),ddu0(0:nz))

    allocate(rhsvn(0:nx,0:ny,0:nz),rhsomegayn(0:nx,0:ny,0:nz))
    allocate(unp1(0:nx,0:ny,0:nz),vnp1(0:nx,0:ny,0:nz))
    allocate(clamb1nm1(2,0:nz),clamb3nm1(2,0:nz),clamb1n(2,0:nz),clamb3n(2,0:nz))
    allocate(wnp1(0:nx,0:ny,0:nz))
    allocate(omegaxnp1(0:nx,0:ny,0:nz))
    allocate(omegaynp1(0:nx,0:ny,0:nz))
    allocate(omegaznp1(0:nx,0:ny,0:nz))
    allocate(au0(2,0:nz),w0(2,0:nz))
    if(id==0) write(*,*)"Matrix allocated."
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>change
    call fftw_f77_create_plan(plan1,2*nz,FFTW_FORWARD,FFTW_MEASURE)!**1
    call fftw_f77_create_plan(plan2,ny+1,FFTW_BACKWARD,FFTW_MEASURE)!**2
    call fftw_f77_create_plan(plan3,nnx,FFTW_BACKWARD,FFTW_MEASURE)!**3
    call fftw_f77_create_plan(plan4,nnx,FFTW_FORWARD,FFTW_MEASURE)!**4
    call fftw_f77_create_plan(plan5,ny+1,FFTW_FORWARD,FFTW_MEASURE)!**5
    call fftw_f77_create_plan(plan6,2*nz,FFTW_FORWARD,FFTW_MEASURE)!**6
    call fftw_f77_create_plan(plan7,2*nz,FFTW_FORWARD,FFTW_MEASURE)!**7
    call fftw_f77_create_plan(plan8,3*(ny+1)/2,FFTW_BACKWARD,FFTW_MEASURE)!**8
    call fftw_f77_create_plan(plan9,3*(nx+1)/2*nproc,FFTW_BACKWARD,FFTW_MEASURE)!**9
    call fftw_f77_create_plan(plan10,3*(nx+1)/2*nproc,FFTW_FORWARD,FFTW_MEASURE)!**10
    call fftw_f77_create_plan(plan11,3*(ny+1)/2,FFTW_FORWARD,FFTW_MEASURE)!**11
    call fftw_f77_create_plan(plan12,2*nz,FFTW_FORWARD,FFTW_MEASURE)!**12
    call fftw_f77_create_plan(plan13,4*nz,FFTW_FORWARD,FFTW_MEASURE)!**13
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call mpi_barrier(mpi_comm_world,ierr)
    if(id==0) write(*,*)"FFT plan created."
    !//////////////////////////////////////////////////////////////
    ! set output parameter nstat
    nstat=20
    ! start mpi enviroment
    pi=4.0d0*atan(1.d0)
    coef=1.0d0/(real(ny+1,kind=8))
    ! some paremeters
    do j=0,(ny+1)/2
        ibeta(j)=real(j,kind=8)/beta
        !        ibeta(j)=0.0d0
        if(j.le.(ny-1)/2-1)then
            !         ibeta(ny-j)=0.0d0
            ibeta(ny-j)=real(-j-1,kind=8)/beta
        endif
    enddo
    ! give ibeta value
    do i=0,nhx
        iapha(i)=(id*real(nx+1,kind=8)*0.5d0+real(i,kind=8))/aphi
    enddo
    ! give the iapha value
    ! give the aphi beta and pi value

    do i=0,nx
        x(i)=2.0d0*pi*(real(id*(nx+1)+i,kind=8))/real((nx+1)*nproc,kind=8)*aphi
    enddo

    do j=0,ny
        y(j)=2.0d0*pi*real(j,kind=8)/real(ny+1,kind=8)*beta
    enddo

    do k=0,nz
        z(k)=cos(pi*real(k,kind=8)/real(nz,kind=8))
    enddo
    ! give the location

    !      if(id.eq.1)then
    !        do i=0,nhx
    !          write(*,*)iapha(i)
    !        enddo
    !      endif
    call SCHEB(t0,t1,t2,t4,u0,du0,ddu0,nz)
    call initial(isign,x,y,z,tvn,tvni,iapha,ibeta,  &
        nallgrp,id,fr,fi)
    !      open(10,file='test.dat',status='replace')
    !      do j=0,ny
    !      do i=0,nx
    !      do k=0,nz
    !        if(abs(rhsvn(i,j,k)).ge.1e-6)then
    !          write(10,*)k,rhsvn(0,0,k),vnp1(0,0,k)
    !        endif
    !      enddo
    !      enddo
    !      enddo
    !      close(10)


    !     ----- begin iteration -----
    ! used by gaussian delimination process, ka = 1 after first time step.
    ! the biggest value of each wave number will be stored in subroutine
    ! after the first time of computation. which can save time.
    ka=0   

    run_time=0.0d0

    !compute coef matrix a, e and its LU decomposition
    do i=0,nhx
    do j=0,ny
        kmn=iapha(i)**2+ibeta(j)**2
        if(kmn.eq.0) cycle ! 'cycle' in Fortran means 'continue' in C
        call get_coef_v(coef_a(:,:,i,j),nz,kmn,iapha(i))
        call get_coef_omegay(coef_e(:,:,i,j),nz,kmn,iapha(i))
        ! if(id.eq.0 .and.  i.eq.0 .and. j.eq.1) then
        !     open(1997,file='matrix_2.txt',status='replace')
        !     write(1997,*) coef_a(:,:,i,j)
        !     write(1997,*) rhs
        !     close(1997)
        ! endif
        !call inverse_matrix(coef_a(:,:,i,j),nz+1)
        !call inverse_matrix(coef_e(:,:,i,j),nz+1)
        ! if(id.eq.0 .and. i.eq.0 .and. j.eq.1) then 
        !     do k=1,5
        !         write(*,*) coef_a(k,1:5,0,1)
        !     enddo
        ! endif
        ! ! get LU decomposition for coef matrix a and e
        call zgetrf(nz+1,nz+1,coef_a(:,:,i,j),nz+1,ipiv(:,i,j),info)   !
        call zgetri(nz+1, coef_a(:,:,i,j),nz+1, ipiv(:,i,j),work,(nz+1)*64,info)
        call zgetrf(nz+1,nz+1,coef_e(:,:,i,j),nz+1,ipiv(:,i,j),info)
        call zgetri(nz+1, coef_e(:,:,i,j),nz+1, ipiv(:,i,j),work,(nz+1)*64,info)
        
        call get_coef_v(a,nz,kmn,iapha(i))
        call matrix_multiply(a,coef_a(:,:,i,j),b,nz+1)
        if(id.eq.0 .and. i.eq.1 .and. j.eq.1) then
            do k=0,6
                write(*,'(100E16.8)') real(b(k,0:6))
            enddo
        endif
    enddo
    enddo

    if(id.eq.0)then !!*
        starttime=mpi_wtime() !!*
    endif !!*
    if(id.eq.0)write(*,*) 'Start main loop...',id
    do istep=istepbegin,nstep   !---------------->>> START MAIN LOOP

        if(istep.ne.istepbegin)then
            call statistics(istep,x,y,z,nstat,id,nallgrp,run_time,ierr)
            call getnonlinear(istep,tvn,tvni,fr,fi,id,nallgrp,iapha,ibeta)
            call output(istep,id)
        endif

        if(id.eq.0)then !!*
            etime1=mpi_wtime() !!*
        endif !!*
        t_eq_time = 0.0d0

        !      goto 1800
        !      if(id.eq.0)write(*,*)istep
        !//////////////////////////////////////////////////////////
        !      if(istep.le.dstep)then
        !       call disturb(istep,y,fr1,fi1,fr2,fi2,nx,ny,nallgrp,id,nproc)
        !      endif
        !//////////////////////////////////////////////////////////
        !
        !	continue
        do i=0,nhx
            do j=0,ny
                ! solve equation of wave number i,j

                kmn=iapha(i)**2+ibeta(j)**2
                if(j.eq.(ny+1)/2.or.kmn.eq.0.0d0)then
                    goto 1234
                endif
                !////////////////////////solve v equation 2.21d0////////////////////////////
                ! guess that unp1 means information of previous timeStep.
                do k=0,nz
                    rhsr(k)=(3.0d0*dt/2.0d0*rhsvn(2*i,j,k)-dt/2.0d0*  &
                        unp1(2*i,j,k))*1.0d0
                    rhsi(k)=(3.0d0*dt/2.0d0*rhsvn(2*i+1,j,k)-dt/2.0d0*  &
                        unp1(2*i+1,j,k))*1.0d0
                enddo

                do k=0,nz
                    tvn(k)=vnp1(2*i,j,k)
                    tvni(k)=vnp1(2*i+1,j,k)
                    tvnm1(k)=tvn(k)
                    tvnm1i(k)=tvni(k)
                enddo
                if(istep.le.dstep)then
                    bdv(1)=0.0d0
                    bdv(2)=0.0d0
                    bdv(3)=fr1(i,j)
                    bdv(4)=fr1(i,j)

                    bdvi(1)=0.0d0
                    bdvi(2)=0.0d0
                    bdvi(3)=fi1(i,j)
                    bdvi(4)=fi1(i,j)
                else

                    bdv(1)=0.0d0
                    bdv(2)=0.0d0
                    bdv(3)=0.0d0
                    bdv(4)=0.0d0

                    bdvi(1)=0.0d0
                    bdvi(2)=0.0d0
                    bdvi(3)=0.0d0
                    bdvi(4)=0.0d0
                endif

                ! Compute matrix Ampn in 2.21
                call get_coef_v(a,nz,kmn,iapha(i))
                ! a(p,q) = Aiqj(y_(p-1));

                ! Compute right hand side of 2.21
                call get_coef_rhs(rhsr,rhsi,tvn,tvni,rhs,nz,kmn,  &
                    iapha(i),bdv,bdvi)

                ! mytemp = rhs
                ! if(id.eq.0 .and. istep.eq.1 .and. i.eq.1 .and. j.eq.1) then
                !     open(1997,file='matrix.txt',status='replace')
                !     write(1997,*) a
                !     write(1997,*) coef_a(:,:,i,j)
                !     write(1997,*) rhs
                ! endif

                s_eq_time = mpi_wtime()
                ! Solve the linear equation 2.21 (a*x=rhs)
                ! call gsij(ka,i,j,a,rhs,nhx+1,ny,nz,id)
                ! call zgetrf(nz+1,nz+1,a,nz+1,ipiv,info)
                ! call zgetrs('N', nz+1, nz+1, coef_a(:,:,i,j), nz+1, &
                !            ipiv(:,i,j), rhs, nz+1, info)
                call a_multiply_x(coef_a(:,:,i,j),rhs,nz+1)

                ! call LUIJ(a,rhs,nz+1)
                ! The results will be stored in rhs(0:nz).
                e_eq_time = mpi_wtime()
                t_eq_time = t_eq_time + (e_eq_time - s_eq_time)
                ! if(id.eq.0 .and. istep.eq.1 .and. i.eq.1 .and. j.eq.1) then
                !     open(1997,file='matrix.txt',status='replace')
                !     write(1997,*) rhs
                !     write(1997,*) mytemp
                !     write(*,*) rhs(1:6)
                !     write(*,*) mytemp(1:6)
                !     close(1997) 
                ! endif                

                ! Copy current step data to last step array.
                do k=0,nz
                    tvn(k)=real(rhs(k),kind=8)
                    tvni(k)=dimag(rhs(k))
                    vnp1(2*i,j,k)=tvn(k)
                    vnp1(2*i+1,j,k)=tvni(k)
                enddo
                !//////////////// solve v equation 2.21 /////////////////
                !------------------------------------------------------------------------
                !///////////// solve omega_y equation 2.22 //////////////
                do k=0,nz
                    rhsr(k)=(3.0d0*dt/2.0d0*rhsomegayn(2*i,j,k)-dt/  &
                        2.0d0*wnp1(2*i,j,k))*1.0d0
                    rhsi(k)=(3.0d0*dt/2.0d0*rhsomegayn(2*i+1,j,k)-  &
                        dt/2.0d0*wnp1(2*i+1,j,k))*1.0d0
                enddo

                do k=0,nz
                    tomegayn(k)=omegaynp1(2*i,j,k)
                    tomegayni(k)=omegaynp1(2*i+1,j,k)
                enddo
                if(istep.le.dstep)then
                    bdomega(1)=fr2(i,j)
                    bdomega(2)=fr2(i,j)
                    bdomegai(1)=fi2(i,j)
                    bdomegai(2)=fi2(i,j)
                else
                    bdomega(1)=0.0d0
                    bdomega(2)=0.0d0
                    bdomegai(1)=0.0d0
                    bdomegai(2)=0.0d0
                endif

                ! compute Empn in 2.22
                call get_coef_omegay(a,nz,kmn,iapha(i))
                ! results stored in a

                ! compute right hand side of 2.22
                call get_rhs_omegay(rhsr,rhsi,rhs,tvn,tvni,tvnm1,tvnm1i, &
                    tomegayn,tomegayni,nz,kmn,iapha(i),ibeta(j), &
                    bdomega,bdomegai)
                ! results stored in rhs

                s_eq_time = mpi_wtime()
                ! solve 2.22 using gauss elimination (a * x = rhs)
                !call gsij(ka,i+nhx+1,j,a,rhs,nhx+1,ny,nz,id)
                ! call zgetrf(nz+1,nz+1,a,nz+1,ipiv,info)
                 ! call zgetrs('N', nz+1, nz+1, coef_e(:,:,i,j), nz+1, &
                 !            ipiv(:,i,j), rhs, nz+1, info)
                 call a_multiply_x(coef_e(:,:,i,j),rhs,nz+1)
                ! call LUIJ(a,rhs,nz+1)
                ! results stored in rhs.
                e_eq_time = mpi_wtime()
                t_eq_time = t_eq_time + (e_eq_time - s_eq_time)

                ! omegaynp1 means omega y, and is stored in module
                do k=0,nz
                    omegaynp1(2*i,j,k)=real(rhs(k),kind=8)
                    omegaynp1(2*i+1,j,k)=dimag(rhs(k))
                enddo

                ! compute u,v,w, omega x,y,z from omega_y and v;
                call getvelocityandvorticity(i,j,iapha,ibeta,tvn,  &
                    tvni,kmn)
                ! u,v,w stored in unp1,vnp1,wnp1
                ! omega x,y,z stored in omegaxnp1,omegaynp1,omegaznp1

1234            continue
            enddo
        enddo
        ! get the value of velocity and vorticity

        !//////////////////////////////////////////////////////////
        !     the following is to deal with the (0,0) wave number /
        !     so it is stored in the 0 computer                   /
        !//////////////////////////////////////////////////////////
        if(id.eq.0)then
            do k=0,nz
                rhsr(k)=-1.0d0*dt/2.0d0*(3.0d0*clamb1n(1,k)-clamb1nm1(1,k))
                rhsi(k)=-1.0d0*dt/2.0d0*(3.0d0*clamb1n(2,k)-clamb1nm1(2,k))

                tvn(k)=au0(1,k)
                tvni(k)=au0(2,k)
            enddo

            bdu(1)=0.0d0
            bdu(2)=0.0d0

            bdui(1)=0.0d0
            bdui(2)=0.0d0

            call get_coef_u0(a,nz)
            call get_rhs_u0(rhsr,rhsi,rhs,tvn,tvni,nz,bdu,bdui)
            call gsij(ka,2*(nhx+1),0,a,rhs,nhx+1,ny,nz,id)
            !-----------------------
            !        ka=1    !elimination information has been saved
            !        call mpi_bcast(ka,1,mpi_integer,0,nallgrp,ierr)
            !-----------------------
            do k=0,nz
                tvn(k)=real(rhs(k),kind=8)
                tvni(k)=dimag(rhs(k))
                au0(1,k)=tvn(k)
                au0(2,k)=tvni(k)
                unp1(0,0,k)=tvn(k)
                unp1(1,0,k)=tvni(k)
            enddo

            call dz(tvn,nz)
            call dz(tvni,nz)

            do k=0,nz
                omegaznp1(0,0,k)=-1.0d0*tvn(k)
                omegaznp1(1,0,k)=-1.0d0*tvni(k)
            enddo

            do k=0,nz
                rhsr(k)=-1.0d0*dt/2.0d0*(3.0d0*clamb3n(1,k)-clamb3nm1(1,k))
                rhsi(k)=-1.0d0*dt/2.0d0*(3.0d0*clamb3n(2,k)-clamb3nm1(2,k))
                tvn(k)=w0(1,k)
                tvni(k)=w0(2,k)
            enddo

            bdu(1)=0.0d0
            bdu(2)=0.0d0

            bdui(1)=0.0d0
            bdui(2)=0.0d0

            call get_coef_u0(a,nz)
            call get_rhs_u0(rhsr,rhsi,rhs,tvn,tvni,nz,bdu,bdui)
            call gsij(ka,2*(nhx+1),ny,a,rhs,nhx+1,ny,nz,id)

            do k=0,nz
                tvn(k)=real(rhs(k),kind=8)
                tvni(k)=dimag(rhs(k))
                w0(1,k)=tvn(k)
                w0(2,k)=tvni(k)
                wnp1(0,0,k)=tvn(k)
                wnp1(1,0,k)=tvni(k)
            enddo

            call dz(tvn,nz)
            call dz(tvni,nz)

            do k=0,nz
                omegaxnp1(0,0,k)=tvn(k)
                omegaxnp1(1,0,k)=tvni(k)
            enddo

        endif !if(id==0)then

        if(id==0)then !!*
            etime2=mpi_wtime()
            !     endtime=endtime-starttime
        endif !!*

        if(mod(istep,interupt).eq.0)then
            if(id==0)then !!*
                write(*,*)"solve eqns time=",etime2-etime1 !!*
                write(*,*)'solve eq group time=',t_eq_time
            endif !!*
        endif
        !      sumu=dsqrt(sum(vnp1**2))
        !      if(id.eq.0)then
        !      write(*,*)sumu,'#1'
        !      end if
1800    continue

        run_time=run_time+dt

        !call outputvelocityandvorticity(istep,run_time,nallgrp,id)
        !call statistics(istep,x,y,z,nstat,id,nallgrp,run_time,ierr)

        !if(mod(istep,interupt).eq.0.and.id==0)then
        !    etime3=mpi_wtime()!!*
        !    write(*,*)"statistics time=",etime3-etime2 !!*
        !endif
        ! putout the velocity and vorticity

        !call getnonlinear(istep,tvn,tvni,fr,fi,id,nallgrp,iapha,ibeta)
        ! the results are said to be stored in omegaxnp1,vn,omegaznp1

        !if(mod(istep,interupt).eq.0.and.id==0)then
        !    etime4=mpi_wtime()!!*
        !    write(*,*)"getnonlinear time=",etime4-etime3 !!*
        !endif
        ! get the new time lamb vector by psudo-spectral
        ! the result has been changed from this subroutine
        ! un,wn has been changed to the n-1 step of the nonlinear part
        !///////////////////////////////////////////////////////////////
        !      if(id.eq.0)then
        !      write(*,*)'#2'
        !      end if

        !call output(istep,id)

        !if(mod(istep,interupt).eq.0.and.id==0)then
        !    endtime=mpi_wtime()!!
        !    write(*,*)"output time=",endtime-etime4 !!*

        !    write(*,*)"from step=1 time=",endtime-starttime !!*
        !endif


        !      call realvelocityandvorticity(istep,x,y,z,id,nallgrp)

        ka=1
    end do   !<<<---------------------END MAIN LOOP
    call fftw_f77_destroy_plan(plan1)
    call fftw_f77_destroy_plan(plan2)
    call fftw_f77_destroy_plan(plan3)
    call fftw_f77_destroy_plan(plan4)
    call fftw_f77_destroy_plan(plan5)
    call fftw_f77_destroy_plan(plan6)
    call fftw_f77_destroy_plan(plan7)
    call fftw_f77_destroy_plan(plan8)
    call fftw_f77_destroy_plan(plan9)
    call fftw_f77_destroy_plan(plan10)
    call fftw_f77_destroy_plan(plan11)
    call fftw_f77_destroy_plan(plan12)
    call fftw_f77_destroy_plan(plan13)

    deallocate(t0,t2,t4,t1)
    deallocate(u0,du0,ddu0)

    deallocate(rhsvn,rhsomegayn)
    deallocate(unp1,vnp1)
    deallocate(clamb1nm1,clamb3nm1,clamb1n,clamb3n)
    deallocate(wnp1)
    deallocate(omegaxnp1)
    deallocate(omegaynp1)
    deallocate(omegaznp1)
    deallocate(au0,w0)

    call mpi_finalize(ierr)
    stop
    end


    !_______________________________________________________
    !-------------------------------------------------------
    !
    !     this program is to get the first order diffrention of y
    !     input v  is to be diffrentiated
    !     output v is the first order
    subroutine dz(v,n)
    real*8 v(0:n)
    real*8 dv(0:n)
    integer i,j
    do i=0,n
        dv(i)=0.0d0
    end do
    !//////////////////////////////////////////
    !	    do i=0,n
    !	     if(abs(v(i)).le.1.e-8)then
    !	         v(i)=0.d0
    !	     endif
    !	    enddo
    !//////////////////////////////////////////
    do i=0,n-1
        do j=i+1,n
            if(mod((i+j),2).ne.0)then
                if(i.eq.0)then
                    dv(i)=dv(i)+real(j,kind=8)*v(j)
                else
                    dv(i)=dv(i)+2.d0*real(j,kind=8)*v(j)
                endif
            endif
        enddo
    enddo
    do i=0,n
        v(i)=dv(i)
    enddo
    return
    end

    !////////////////////////////////////////////////////////////////////
    !     this program is to give the coeffience of Chebshev seres
    !     input(0:n) is the array which is to be transformed
    !     n is the number to be transfromed
    !     assist(0:2*n) is the assistant array for transformation
    !     the result is saved in the input array
    !     if the isign is equal 1,it represent the inverse transformation
    !     if the isign is equal -1,it represent the direct transformation
    !////////////////////////////////////////////////////////////////////
    subroutine chebyshev(insign,input,n)
    use paraandcons !!*
    include "fftw_f77.i"
    real*8 input(0:n)
    real*8 assist(0:2*n)
    complex(kind=8) in(2*n),out(2*n)
    integer insign,n,n2
    !      integer*8 plan
    real*8 dd,coef
    n2=2*n
    dd=1.0d0-0.5d0*real((1+insign)/2)
    coef=(0.5d0-1.0d0/real(4*n,kind=8))*real(insign,kind=8)+(0.5d0+1.0d0/real(4*n,kind=8))
    assist(0)=input(0)
    assist(n)=input(n)

    do j=1,n-1
        assist(j)=dd*input(j)
        assist(n2-j)=assist(j)
    end do
    do j=0,n2-1
        in(j+1)=assist(j)
    end do
    !     call fftw_f77_create_plan(plan,n2,FFTW_FORWARD,FFTW_ESTIMATE)!**1
    call fftw_f77_one(plan1,in,out)
    out=out*coef
    do j=0,n2-1
        assist(j)=out(j+1)
    end do
    input(0)=assist(0)
    input(n)=assist(n)
    do j=1,n-1
        input(j)=dd*(assist(j)+assist(n2-j))
    end do
    !      call fftw_f77_destroy_plan(plan)
    return
    end


    !-----------------------------------------------------------------------
    !     this program is to be uesed for parallel computering
    !     in pseudo spectral method,we need to collect the mode of all
    !     nodes,but for we expand the mode to 3/2 times mode,neccesary
    !     method is needed for communication
    !     if isign.eq.1
    !     distribute the modes that are collected from the nodes in order
    !     to guarantee the natural order of the mode
    !     if isign.eq.-1
    !     distribute the modes that have been transformed to the nodes
    !     in order to guarantee that every node has the right original
    !     node
    !     input u is the array that will be distributed
    !     mx,my,mz is the total number of each direction
    !     mproc is the total number of processors
    !------------------------------------------------------------------------

    subroutine equipment(isign,u,mx,my,mz,mproc)
    integer i,j,isgin,mp,k,mx,my,mz,moy,l
    real*8 u(0:mx,0:my,0:mz)
    real*8 au(0:mx,0:(my+1)*2/3-1)
    moy=(my+1)*2/3-1
    mp=(moy+1)/(mproc)
    if(isign.eq.1)then
        do l=0,mz
            do i=0,mx
                do k=0,mproc-1
                    do j=0,mp-1
                        au(i,k*(mp)+j)=u(i,k*(3*mp/2)+j,l)
                    enddo
                enddo
            enddo

            do i=0,mx
                do j=0,moy
                    u(i,j,l)=au(i,j)
                enddo
            enddo

        enddo
    endif

    if(isign.eq.-1)then
        do l=0,mz
            do i=0,mx
                do j=0,moy
                    au(i,j)=u(i,j,l)
                enddo
            enddo
            do i=0,mx
                do k=0,mproc-1
                    do j=0,mp-1
                        u(i,k*(3*mp/2)+j,l)=au(i,k*(mp)+j)
                    enddo
                enddo
            enddo
        enddo
    endif
    return
    end

    !-------------------------------------------------------------
    !   this program is to expand the original Fourier modes to
    !   3/2 modes
    !   if isign=1, the original to the 3/2
    !   if isign=-1,the 3/2 to original
    !
    !   fr is the real part of the original mode
    !   fi is the image part of the original mode
    !   fry is the real part of the 3/2 mode
    !   fiy is the image part of the 3/2 mode
    !-------------------------------------------------------------

    subroutine expand(isign,fr,fi,fry,fiy,mx,my)
    integer mx,my,amy,isign
    real*8 fr(0:mx,0:my),fi(0:mx,0:my)
    real*8 fry(0:mx,0:3*(my+1)/2-1)
    real*8 fiy(0:mx,0:3*(my+1)/2-1)
    amy=3*(my+1)/2-1
    if(isign.eq.1)then

        do i=0,mx
            fr(i,(my+1)/2)=0.0d0
            fi(i,(my+1)/2)=0.0d0
            do j=0,(my+1)/2
                fry(i,j)=fr(i,j)
                fiy(i,j)=fi(i,j)
                if(j.le.((my-1)/2-1))then
                    fry(i,(amy-j))=fr(i,(my-j))
                    fiy(i,(amy-j))=fi(i,(my-j))
                endif
                if(j.lt.(my+1)/4)then
                    fry(i,(my+1)/2+1+j)=0.0d0
                    fiy(i,(my+1)/2+1+j)=0.0d0
                    fry(i,(amy+1)/2+1+j)=0.0d0
                    fiy(i,(amy+1)/2+1+j)=0.0d0
                endif
            enddo
        enddo

    endif

    if(isign.eq.-1)then

        do i=0,mx
            do j=0,(my+1)/2
                fr(i,j)=fry(i,j)
                fi(i,j)=fiy(i,j)
                if(j.le.((my-1)/2-1))then
                    fr(i,(my-j))=fry(i,(amy-j))
                    fi(i,(my-j))=fiy(i,(amy-j))
                endif
            enddo
        enddo

    endif
    return
    end


    !----------------------------------------------------------------
    !     This subroutine is developed for fourier transformation  -!
    !     and inverse transformation: isign=1 for inverse trans-   -!
    !     formation and isign=-1 for direct transformation.        -!
    !----------------------------------------------------------------

    subroutine fourtrx_zcomp(isign,w,id,nallgrp,mx,my,mz,mproc)
    use paraandcons
    include "fftw_f77.i"
    integer i,j,k
    real*8 w(0:mx,0:my,0:mz),wp(0:(my+1)/mproc-1,0:(mx+1)*mproc-1,0:mz)
    !complex inx, outx,iny, outy
    COMPLEX(kind=8) ::  inx((mx+1)*mproc),outx((mx+1)*mproc), &
        iny(my+1),outy(my+1)
    real*8 fr(0:(mx-1)/2,0:my)
    real*8 fi(0:(mx-1)/2,0:my)
    real*8 fry(0:(my+1)/mproc-1,0:(mx+1)*mproc-1)
    real*8 fiy(0:(my+1)/mproc-1,0:(mx+1)*mproc-1)
    integer mx,my,mz,amz,nn,amx,amy,hmx,mproc
    !      integer*8 plan
    hmx=(mx-1)/2
    amx=(mx+1)*mproc
    amy=my+1

    if(isign.eq.1)then
        !/////////////////////////////////////////
        !     first    tranfrom in y direction
        !/////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,amy,FFTW_BACKWARD,FFTW_ESTIMATE)!**2

        do k=0,mz
            do j=0,my
                do i=0,(mx-1)/2
                    fr(i,j)=w(2*i,j,k)
                    fi(i,j)=w(2*i+1,j,k)
                enddo
            enddo

            do i=0,(mx-1)/2
                do j=0,my
                    iny(j+1)=dcmplx(fr(i,j),fi(i,j))
                end do

                iny((my+1)/2+1)=dcmplx(0.0d0,0.0d0)
                call fftw_f77_one(plan2,iny,outy)

                do j=0,my
                    w(2*i,j,k)=real(outy(j+1),kind=8)
                    w(2*i+1,j,k)=dimag(outy(j+1))
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)
        !/////////////////////////////////////////////////
        !     the following is to transpoze the matrix
        !/////////////////////////////////////////////////

        call transposefft_c_r(w,wp,id,nallgrp,mproc,mx,my,mz)
        !//////////////////////////////////////////////////
        !     the following is to transform in x direction
        !//////////////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,amx,FFTW_BACKWARD,FFTW_ESTIMATE)!**3

        do k=0,mz
            do j=0,((mx+1)*mproc-2)/2
                do i=0,(my+1)/mproc-1
                    fry(i,j)=wp(i,2*j,k)
                    fiy(i,j)=wp(i,2*j+1,k)
                enddo
            enddo
            call getall(fry,fiy,(my+1)/mproc-1,(mx+1)*mproc-1)

            do i=0,(my+1)/mproc-1
                do j=0,(mx+1)*mproc-1
                    inx(j+1)=dcmplx(fry(i,j),fiy(i,j))
                enddo
                inx(amx/2+1)=dcmplx(0.0d0,0.0d0)

                call fftw_f77_one(plan3,inx,outx)

                !	   call fftw_f77(plan,1,in,1,0,out,1,0)

                do j=0,(mx+1)*mproc-1
                    wp(i,j,k)=real(outx(j+1),kind=8)
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)
        call transposefft_r_c(w,wp,id,nallgrp,mproc,mx,my,mz)
    end if

    if(isign.eq.-1)then
        call transposefft_c_r(w,wp,id,nallgrp,mproc,mx,my,mz)
        !//////////////////////////////////////////////////
        !     the following is to transform in x direction from real*8 to fourier
        !//////////////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,amx,FFTW_FORWARD,FFTW_ESTIMATE)!**4
        do k=0,mz
            do i=0,(my+1)/mproc-1
                do j=0,(mx+1)*mproc-1
                    inx(j+1)=dcmplx(wp(i,j,k),0.0d0)
                enddo

                call fftw_f77_one(plan4,inx,outx)

                do j=0,(amx-2)/2
                    wp(i,2*j,k)=real(outx(j+1),kind=8)
                    wp(i,2*j+1,k)=dimag(outx(j+1))
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)
        !/////////////////////////////////////////////////
        !     the following is to transpoze the matrix
        !/////////////////////////////////////////////////
        call transposefft_r_c(w,wp,id,nallgrp,mproc,mx,my,mz)

        !//////////////////////////////////////////////////
        !     the following is to transform in x direction from real*8 to fourier
        !//////////////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,amy,FFTW_FORWARD,FFTW_ESTIMATE)!**5
        do  k=0,mz
            do i=0,hmx
                do j=0,my
                    iny(j+1)=dcmplx(w(2*i,j,k),w(2*i+1,j,k))
                enddo
                call fftw_f77_one(plan5,iny,outy)

                do j=0,my
                    fr(i,j)=real(outy(j+1),kind=8)
                    fi(i,j)=dimag(outy(j+1))
                enddo
            enddo

            do i=0,hmx
                do j=0,my
                    w(2*i,j,k)=fr(i,j)/real(amx*amy,kind=8)
                    w(2*i+1,j,k)=fi(i,j)/real(amx*amy,kind=8)
                enddo
            enddo
        enddo
        !     call fftw_f77_destroy_plan(plan)

    endif
    return
    end

    !---------------------------------------------------------
    !   isign=1: K-space --> Phys space; isign=-1: Phys space-
    !   --> K-space.                                         -
    !---------------------------------------------------------
    subroutine transform(isign,u,id,nallgrp,mx,my,mz,mproc)
    use paraandcons
    include "fftw_f77.i"
    !complex in, out
    complex(kind=8) ::  in(2*mz),out(2*mz)

    real*8 u(0:mx,0:my,0:mz)
    real*8 tmpu(0:mz)

    integer i,j,k,mx,my,mz,mproc,amx,amy,amz,mn
    !      integer*8 plan
    mn=2*mz
    if(isign.eq.1)then

        call fourtrx_zcomp(1,u,id,nallgrp,mx,my,mz,mproc)

        !      call fftw_f77_create_plan(plan,mn,FFTW_FORWARD,FFTW_ESTIMATE)!**6
        do ii=0,mx
            do j=0,my
                do k=0,mz
                    tmpu(k)=u(ii,j,k)
                enddo
                !///////////////////////////////////////////////////////
                !     the following is the transform	in z -direction
                !-------------------------------------------------------

                do i=0,mz
                    in(i+1)=dcmplx(tmpu(i),0)
                    if(i.lt.(mz-1))then
                        in(mn-i)=dcmplx(tmpu(i+1),0)
                    endif
                enddo
                in(1)=in(1)*2.0d0
                in(mz+1)=in(mz+1)*2.0d0

                call fftw_f77_one(plan6,in,out)

                do i=0,mz
                    tmpu(i)=real(out(i+1),kind=8)/2.0d0
                enddo

                do k=0,mz
                    u(ii,j,k)=tmpu(k)
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)

    end if

    if(isign.eq.-1)then   !transfer from physical space to Fourier space.

        !      call fftw_f77_create_plan(plan,mn,FFTW_FORWARD,FFTW_ESTIMATE)!**7
        do ii=0,mx
            do j=0,my
                do k=0,mz
                    tmpu(k)=u(ii,j,k)
                enddo
                !///////////////////////////////////////////////////////
                !     the following is the transform	in z -direction
                !-------------------------------------------------------

                do i=0,mz
                    in(i+1)=dcmplx(tmpu(i),0)
                    if(i.lt.(mz-1))then
                        in(mn-i)=dcmplx(tmpu(i+1),0)
                    endif
                enddo

                call fftw_f77_one(plan7,in,out)

                do i=0,mz
                    tmpu(i)=real(out(i+1),kind=8)/real(mz,kind=8)
                enddo
                tmpu(0)=tmpu(0)/2.0d0
                tmpu(mz)=tmpu(mz)/2.0d0
                do k=0,mz
                    u(ii,j,k)=tmpu(k)
                enddo
            enddo
        enddo
        !     call fftw_f77_destroy_plan(plan)

        call fourtrx_zcomp(-1,u,id,nallgrp,mx,my,mz,mproc)
    end if
    return
    end

    !----------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------
    !    this program is to transfrom  from fourier to real and real to fourier
    !    also this program is to deal with deleting alias error
    !    if isign  eq 1 transform from fourier to real equiped with deleting alias error
    !    if isign eq -1 transform from real to fourier equiped with deleting alias error
    !    if isign eq 0 transform from fourier to real without alias error
    !    input
    !    u is the orignal array
    !    up is the expanded array except z direction
    !    fr fi is the assistant array
    !    mx my mz is individually the node number of the x,y,z direction
    !    mproc is the the number of the computer
    !-----------------------------------------------------------------------------------

    subroutine ftrthreetotwo(isign,w,wpp,id,nallgrp,fr,fi  &
        ,mx,my,mz,mproc)
    use paraandcons
    include "fftw_f77.i"
    integer i,j,k
    real*8 w(0:mx,0:my,0:mz)
    real*8 wp(0:3*(mx+1)/2-1,0:3*(my+1)/2-1,0:mz), &
        wpp(0:3*(my+1)/2/mproc-1,0:3*(mx+1)/2*mproc-1,0:mz)
    !complex inx, outx,iny, outy
    complex(kind=8) ::  inx(3*(mx+1)*mproc/2),outx(3*(mx+1)*mproc/2), &
        iny(3*(my+1)/2),outy(3*(my+1)/2)
    real*8 fr(0:(mx-1)/2,0:my)
    real*8 fi(0:(mx-1)/2,0:my)
    real*8 fry(0:(mx-1)/2,0:3*(my+1)/2-1)
    real*8 fiy(0:(mx-1)/2,0:3*(my+1)/2-1)
    real*8 frx(0:3*(my+1)/2/mproc-1,0:(mx+1)*mproc-1)
    real*8 fix(0:3*(my+1)/2/mproc-1,0:(mx+1)*mproc-1)
    real*8 fra(0:3*(my+1)/2/mproc-1,0:3*(mx+1)*mproc/2-1)
    real*8 fia(0:3*(my+1)/2/mproc-1,0:3*(mx+1)*mproc/2-1)
    integer mx,my,mz,amz,nn,amx,amy,pmx,pmy,hmx,mproc
    !      integer*8 plan
    hmx=(mx-1)/2
    amx=3*(mx+1)/2-1
    amy=3*(my+1)/2-1

    pmx=3*(my+1)/2/mproc-1
    pmy=3*(mx+1)/2*mproc-1

    if(isign.eq.1)then
        !/////////////////////////////////////////
        !     first    tranfrom in y direction
        !/////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,amy+1,FFTW_BACKWARD,FFTW_ESTIMATE)!**8

        do k=0,mz
            do j=0,my
                do i=0,(mx-1)/2
                    fr(i,j)=w(2*i,j,k)
                    fi(i,j)=w(2*i+1,j,k)
                enddo
            enddo
            call expand(1,fr,fi,fry,fiy,hmx,my)
            do i=0,(mx-1)/2
                do j=0,amy
                    iny(j+1)=dcmplx(fry(i,j),fiy(i,j))
                end do

                iny((amy+1)/2+1)=dcmplx(0.0d0,0.0d0)

                call fftw_f77_one(plan8,iny,outy)

                !	   call fftw_f77(plan,1,in,1,0,out,1,0)

                do j=0,amy
                    wp(2*i,j,k)=real(outy(j+1),kind=8)
                    wp(2*i+1,j,k)=dimag(outy(j+1))
                enddo
            enddo
        enddo

        !      call fftw_f77_destroy_plan(plan)
        !/////////////////////////////////////////////////
        !     the following is to transpoze the matrix
        !/////////////////////////////////////////////////

        call transposefft_c_r(wp,wpp,id,nallgrp,mproc,amx,amy,mz)
        !//////////////////////////////////////////////////
        !     the following is to transform in x direction
        !//////////////////////////////////////////////////

        call equipment(1,wpp,pmx,pmy,mz,mproc)

        !      call fftw_f77_create_plan(plan,pmy+1,FFTW_BACKWARD,FFTW_ESTIMATE)!**9
        do k=0,mz
            do j=0,((mx+1)*mproc-2)/2
                do i=0,pmx
                    frx(i,j)=wpp(i,2*j,k)
                    fix(i,j)=wpp(i,2*j+1,k)
                enddo
            enddo
            call getall(frx,fix,pmx,(mx+1)*mproc-1)
            call expand(1,frx,fix,fra,fia,pmx,(mx+1)*mproc-1)
            do i=0,pmx
                do j=0,pmy
                    inx(j+1)=dcmplx(fra(i,j),fia(i,j))
                enddo
                inx((pmy+1)/2+1)=dcmplx(0.0d0,0.0d0)

                call fftw_f77_one(plan9,inx,outx)

                !	   call fftw_f77(plan,1,in,1,0,out,1,0)

                do j=0,pmy
                    wpp(i,j,k)=real(outx(j+1),kind=8)
                enddo
            enddo
        enddo
        !     call fftw_f77_destroy_plan(plan)

    end if
    if(isign.eq.-1)then

        !//////////////////////////////////////////////////
        !     the following is to transform in x direction from real*8 to fourier
        !//////////////////////////////////////////////////
        !      call fftw_f77_create_plan(plan,pmy+1,FFTW_FORWARD,FFTW_ESTIMATE)!**10
        do k=0,mz
            do i=0,pmx
                do j=0,pmy
                    inx(j+1)=dcmplx(wpp(i,j,k),0.0d0)
                enddo
                call fftw_f77_one(plan10,inx,outx)

                !          call fftw_f77(plan,1,in,1,0,out,1,0)

                do j=0,((mx+1)*mproc-2)/2
                    wpp(i,2*j,k)=real(outx(j+1),kind=8)
                    wpp(i,2*j+1,k)=dimag(outx(j+1))
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)

        call equipment(-1,wpp,pmx,pmy,mz,mproc)

        !/////////////////////////////////////////////////
        !     the following is to transpoze the matrix
        !/////////////////////////////////////////////////
        call transposefft_r_c(wp,wpp,id,nallgrp,mproc,amx,amy,mz)

        !      call fftw_f77_create_plan(plan,amy+1,FFTW_FORWARD,FFTW_ESTIMATE)!**11
        !//////////////////////////////////////////////////
        !     the following is to transform in x direction from real*8 to fourier
        !//////////////////////////////////////////////////
        do  k=0,mz
            do i=0,hmx
                do j=0,amy
                    iny(j+1)=dcmplx(wp(2*i,j,k),wp(2*i+1,j,k))
                enddo
                call fftw_f77_one(plan11,iny,outy)

                !          call fftw_f77(plan,1,in,1,0,out,1,0)

                do j=0,amy
                    fry(i,j)=real(outy(j+1),kind=8)
                    fiy(i,j)=dimag(outy(j+1))
                enddo
            enddo

            call expand(-1,fr,fi,fry,fiy,hmx,my)
            do i=0,hmx
                do j=0,my
                    w(2*i,j,k)=fr(i,j)/real((amy+1)*(pmy+1),kind=8)
                    w(2*i+1,j,k)=fi(i,j)/real((amy+1)*(pmy+1),kind=8)
                enddo
            enddo
        enddo
        !      call fftw_f77_destroy_plan(plan)

    endif
    return
    end

    !--------------------------------------------------------------------
    !     this program is to transpose the matrix in parallel computers
    !     complex_to_real: input
    !     ux is the array which is to be transposed
    !     id is the number of the processor
    !     nallgrp is the communication of the parallel compters
    !     nproc is the number of the total processors
    !     nx,ny,nz is the node number individual x,y,z directions
    !     output
    !     ux is the array that has been transposed
    !--------------------------------------------------------------------
    subroutine transposefft_c_r(ux,uxp,id,nallgrp,nproc,nx,ny,nz)
    include 'mpif.h'
    real*8 ux(0:nx,0:ny,0:nz), &
        uxp(0:(ny+1)/nproc-1,0:(nx+1)*nproc-1,0:nz)
    real*8 tmp1(nx+1,(ny+1)/nproc,nz+1), &
        tmp2(nx+1,(ny+1)/nproc,nz+1)
    real*8 tmp((nx+1)*nproc,(ny+1)/nproc,nz+1)
    integer isize,nxm,nxp,status(MPI_STATUS_SIZE,2),req(2)

    isize=(nx+1)*(ny+1)*(nz+1)/nproc
    lly=(ny+1)/nproc

    do i=1,nproc-1
        nxp=mod(id+i,nproc)
        nxm=mod(id-i+nproc,nproc)
        js=nxp*lly
        do k=1,nz+1
            do j=1,lly
                jl=js+j
                do ii=1,nx+1
                    tmp1(ii,j,k)=ux(ii-1,jl-1,k-1)
                end do
            end do
        end do

        call mpi_isend(tmp1,isize,MPI_DOUBLE_PRECISION,nxp,i,  &
            nallgrp,req(1),ierr)
        call mpi_irecv(tmp2,isize,MPI_DOUBLE_PRECISION,nxm,i,  &
            nallgrp,req(2),ierr)
        call MPi_Waitall(2,req,status,ierr)
        is=nxm*(nx+1)

        do ii=1,nx+1
            il=is+ii
            do j=1,lly
                do k=1,nz+1
                    tmp(il,j,k)=tmp2(ii,j,k)
                end do
            end do
        end do

    end do

    is=id*(nx+1)
    js=id*lly
    do i=1,nx+1
        il=is+i
        do j=1,(ny+1)/nproc
            jl=js+j
            do k=1,nz+1
                tmp(il,j,k)=ux(i-1,jl-1,k-1)
            end do
        end do
    end do

    do i=1,(ny+1)/nproc
        do j=1,(nx+1)*nproc
            do k=1,nz+1
                uxp(i-1,j-1,k-1)=tmp(j,i,k)
            end do
        end do
    end do
    return
    end

    !--------------------------------------------------------------------
    !     this program is to transpose the matrix in parallel computers
    !     real_to_complex: input
    !     ux is the array which is to be transposed
    !     id is the number of the processor
    !     nallgrp is the communication of the parallel compters
    !     nproc is the number of the total processors
    !     nx,ny,nz is the node number individual x,y,z directions
    !     output
    !     ux is the array that has been transposed
    !--------------------------------------------------------------------
    subroutine transposefft_r_c(ux,uxp,id,nallgrp,nproc,nx,ny,nz)
    include 'mpif.h'
    real*8 ux(0:nx,0:ny,0:nz), &
        uxp(0:(ny+1)/nproc-1,0:(nx+1)*nproc-1,0:nz)
    real*8 tmp1((ny+1)/nproc,nx+1,nz+1), &
        tmp2((ny+1)/nproc,nx+1,nz+1)
    real*8 tmp(ny+1,nx+1,nz+1)
    integer isize,nxm,nxp,status(MPI_STATUS_SIZE,2),req(2)

    isize=(nx+1)*(ny+1)*(nz+1)/nproc
    llx=nx+1
    do i=1,nproc-1
        nxp=mod(id+i,nproc)
        nxm=mod(id-i+nproc,nproc)
        js=nxp*llx
        do k=1,nz+1
            do j=1,llx
                jl=js+j
                do ii=1,(ny+1)/nproc
                    tmp1(ii,j,k)=uxp(ii-1,jl-1,k-1)
                end do
            end do
        end do

        call mpi_isend(tmp1,isize,MPI_DOUBLE_PRECISION,nxp,i,  &
            nallgrp,req(1),ierr)
        call mpi_irecv(tmp2,isize,MPI_DOUBLE_PRECISION,nxm,i,  &
            nallgrp,req(2),ierr)
        call MPi_Waitall(2,req,status,ierr)
        is=nxm*(ny+1)/nproc

        do ii=1,(ny+1)/nproc
            il=is+ii
            do j=1,llx
                do k=1,nz+1
                    tmp(il,j,k)=tmp2(ii,j,k)
                end do
            end do
        end do

    end do

    is=id*(ny+1)/nproc
    js=id*llx
    do i=1,(ny+1)/nproc
        il=is+i
        do j=1,nx+1
            jl=js+j
            do k=1,nz+1
                tmp(il,j,k)=uxp(i-1,jl-1,k-1)
            end do
        end do
    end do

    do i=1,nx+1
        do j=1,ny+1
            do k=1,nz+1
                ux(i-1,j-1,k-1)=tmp(j,i,k)
            end do
        end do
    end do

    return
    end

    !---------------------------------------------------------------------
    !     Because for saving memory,we fully use the symmetry of the array
    !     this program is to transform from the Fourier to real space
    !     input
    !     fr  the real part of the fourier series
    !     fi  the image part of the fourier series
    !     we must pay attention to that the input array contains just the
    !     wave numbers that are from 0 to ((my-1)/2)
    !     output
    !     fr and fi are the full wave numbers
    !---------------------------------------------------------------------
    subroutine getall(fr,fi,mx,my)
    real*8 fr(0:mx,0:my),fi(0:mx,0:my)
    integer i,j,mx,my
    do i=0,mx
        do j=0,(my-1)/2-1
            fr(i,my-j)=fr(i,j+1)
            fi(i,my-j)=-1.0d0*fi(i,j+1)
        end do
        fr(i,(my+1)/2)=0.0d0
        fi(i,(my+1)/2)=0.0d0
    end do
    return
    end

    !---------------------------------------------------------------------
    !     this program is to get the nonlinear part of the n time step
    !     input u,v,w,omegax,omegay,omegaz
    !     output
    !     in order to save memory in this program we use
    !     omegaxnp1,vn,omegaznp1 instead of lamb1 lamb2 and lamb3
    !     output the right hand of the equation v and omegay
    !     in here,we output it in unp1,wnp1,the reason is the same.
    !     output:
    !     rhsvn is the n time step result of the v equation
    !     rhsomegayn the n time step result of the omegay equation
    !     unp1 is the n-1 time step result of the v equation
    !     wnp1 is the n-1 time step result of the omegay equation
    !     also clamb1n clamb1nm1 clamb3n clamb3nm1 is what we want to
    !     solve the mod(0,0) velocity
    !-----------------------------------------------------------------------
    subroutine getnonlinear(istep,tvn,tvni,fr,fi,  &
        id,nallgrp,iapha,ibeta)
    use paraandcons
    include 'mpif.h'
    include "fftw_f77.i"
    !complex in, out
    complex(kind=8) ::  in(2*nz),out(2*nz)
    real*8 ibeta(0:ny),iapha(0:nhx)
    real*8 fr1(0:nx,0:ny)
    real*8 fi1(0:nx,0:ny)
    real*8 fr(0:(nx-1)/2,0:ny)
    real*8 fi(0:(nx-1)/2,0:ny)
    real*8 tvn(0:nz),tvni(0:nz),kmn,tmp
    real*8 lamb1(0:nx,0:ny,0:nz)
    real*8 lamb2(0:nx,0:ny,0:nz)
    integer i,j,k,amz
    !     integer*8 plan
    real*8 coriolist(0:nx,0:ny,0:nz) !changeRPCF
    real*8 coriolinm(0:nx,0:ny,0:nz) !changeRPCF
    real*8 coriolisp(0:nx,0:ny,0:nz) !changeRPCF
    real*8 unp1p(0:nx,0:ny,0:nz) !changeRPCF
    real stime,etime,ttime
    ttime = 0
    !====================================================
    ! if(rsp.eq.real(0,kind=8))then
    !     coriolist=0.0d0 !changeRPCF
    !     coriolinm=0.0d0 !changeRPCF
    !     coriolisp=0.0d0 !changeRPCF

    ! else
    !     unp1p=unp1
    !     ! transform unp1p from K-space --> Phys space
    !     call cpu_time(stime)
    !     call transform(1,unp1p,id,nallgrp,nx,ny,nz,nproc)
    !     call cpu_time(etime)
    !     ttime = ttime + (etime - stime)
    !     do k=0,nz
    !         unp1p(:,:,k)=unp1p(:,:,k)+u0(k)
    !     enddo
    !     ! transform back
    !     call cpu_time(stime)
    !     call transform(-1,unp1p,id,nallgrp,nx,ny,nz,nproc)
    !     call cpu_time(etime)
    !     ttime = ttime + (etime - stime)

    !     coriolist=rnm*wnp1-rsp*vnp1 !changeRPCF
    !     coriolinm=rsp*unp1p-rst*wnp1 !changeRPCF
    !     coriolisp=rst*vnp1-rnm*unp1p !changeRPCF
    ! endif
    !====================================================


    amz=2*nz

    ! u in omega y np1, omega in wnp1, results in lamb1
    ! call getlamb(1,omegaynp1,wnp1,lamb1,id,nallgrp,  &
    !     fr,fi,nx,ny,nz,nproc)
    ! mean pressure gradient(-1) considered here.
    call getlamb1(1,omegaynp1,wnp1,lamb1,id,nallgrp,  &
        fr,fi,nx,ny,nz,nproc)
    call getlamb(1,omegaznp1,vnp1,lamb2,id,nallgrp,  &
        fr,fi,nx,ny,nz,nproc)
    do k=0,nz
        do j=0,ny
            do i=0,nx
                lamb1(i,j,k)=lamb1(i,j,k)-lamb2(i,j,k)
            enddo
        enddo
    enddo

    call getlamb(1,omegaznp1,unp1,lamb2,id,nallgrp, &
        fr,fi,nx,ny,nz,nproc)
    call getlamb(1,omegaxnp1,wnp1,omegaznp1,id,nallgrp,  &
        fr,fi,nx,ny,nz,nproc)

    do k=0,nz
        do j=0,ny
            do i=0,nx
                lamb2(i,j,k)=lamb2(i,j,k)-omegaznp1(i,j,k)
            enddo
        enddo
    enddo

    call getlamb(1,omegaxnp1,vnp1,omegaznp1,id,nallgrp,  &
        fr,fi,nx,ny,nz,nproc)
    call getlamb(1,omegaynp1,unp1,omegaxnp1,id,nallgrp,  &
        fr,fi,nx,ny,nz,nproc)
    do k=0,nz
        do j=0,ny
            do i=0,nx
                omegaznp1(i,j,k)=omegaznp1(i,j,k)-omegaxnp1(i,j,k)
            enddo
        enddo
    enddo

    !============update: Lamb vector + Corioli Force===========================================================================
    ! lamb1=lamb1+coriolist !changeRPCF !streamwise direction
    ! lamb2=lamb2+coriolinm !changeRPCF !normal direction
    ! omegaznp1=omegaznp1+coriolisp !changeRPCF !spanwise direction
    !=======================================================================================


    do k=0,nz
        clamb1nm1(1,k)=lamb1(0,0,k)
        clamb3nm1(1,k)=omegaznp1(0,0,k)
        clamb1nm1(2,k)=lamb1(1,0,k)
        clamb3nm1(2,k)=omegaznp1(1,0,k)
    enddo

    do k=0,nz
        do j=0,ny
            do i=0,nx
                omegaxnp1(i,j,k)=lamb1(i,j,k)
                unp1(i,j,k)=lamb2(i,j,k)
            enddo
        enddo
    enddo

    !omegaxnp1(:,(ny+1)/2,:)=0.0d0
    !unp1(:,(ny+1)/2,:)=0.0d0
    !omegaznp1(:,(ny+1)/2,:)=0.0d0
    !     get the nonlinear part of the n time step.............change
    do i=0,nhx
        do j=0,ny
            !         do k=0,nz
            kmn=iapha(i)**2+ibeta(j)**2

            do kk=0,nz
                tvn(kk)=omegaxnp1(2*i,j,kk)
                tvni(kk)=omegaxnp1(2*i+1,j,kk)
            enddo

            call dz(tvn,nz)
            call dz(tvni,nz)

            do k=0,nz
                unp1(2*i,j,k)=kmn*unp1(2*i,j,k)-iapha(i)*tvni(k)
                unp1(2*i+1,j,k)=kmn*unp1(2*i+1,j,k)+iapha(i)*tvn(k)
            enddo

            do kk=0,nz
                tvn(kk)=omegaznp1(2*i,j,kk)
                tvni(kk)=omegaznp1(2*i+1,j,kk)
            enddo

            call dz(tvn,nz)
            call dz(tvni,nz)

            do k=0,nz
                unp1(2*i,j,k)=unp1(2*i,j,k)-ibeta(j)*tvni(k)
                unp1(2*i+1,j,k)=unp1(2*i+1,j,k)+ibeta(j)*tvn(k)

                wnp1(2*i,j,k)=-1.0d0*iapha(i)*omegaznp1(2*i+1,j,k)+  &
                    ibeta(j)*omegaxnp1(2*i+1,j,k)
                wnp1(2*i+1,j,k)=iapha(i)*omegaznp1(2*i,j,k)-  &
                    ibeta(j)*omegaxnp1(2*i,j,k)

            enddo

        enddo
    enddo

    do k=0,nz
        do j=0,ny
            do i=0,nx
                tmp=unp1(i,j,k)
                unp1(i,j,k)=rhsvn(i,j,k)
                rhsvn(i,j,k)=tmp
                tmp=wnp1(i,j,k)
                wnp1(i,j,k)=rhsomegayn(i,j,k)
                rhsomegayn(i,j,k)=tmp
            enddo
        enddo
    enddo
    do k=0,nz
        tmp=clamb1nm1(1,k)
        clamb1nm1(1,k)=clamb1n(1,k)
        clamb1n(1,k)=tmp

        tmp=clamb1nm1(2,k)
        clamb1nm1(2,k)=clamb1n(2,k)
        clamb1n(2,k)=tmp

        tmp=clamb3nm1(1,k)
        clamb3nm1(1,k)=clamb3n(1,k)
        clamb3n(1,k)=tmp

        tmp=clamb3nm1(2,k)
        clamb3nm1(2,k)=clamb3n(2,k)
        clamb3n(2,k)=tmp
    enddo

    !      call fftw_f77_create_plan(plan,amz,FFTW_FORWARD,FFTW_ESTIMATE)!**12
    do j=0,nx
        do k=0,ny
            do i=0,nz
                in(i+1)=dcmplx(rhsvn(j,k,i),rhsomegayn(j,k,i))
                if(i.lt.(nz-1))then
                    in(amz-i)=dcmplx(rhsvn(j,k,i+1),rhsomegayn(j,k,i+1))
                endif
            enddo
            in(1)=in(1)*2.0d0
            in(nz+1)=in(nz+1)*2.0d0
            call fftw_f77_one(plan12,in,out)

            do i=0,nz
                rhsvn(j,k,i)=real(out(i+1),kind=8)/2.0d0
                rhsomegayn(j,k,i)=dimag(out(i+1))/2.0d0
            enddo
        enddo
    enddo

    do k=1,2
        do i=0,nz
            in(i+1)=dcmplx(clamb1n(k,i),clamb3n(k,i))
            if(i.lt.(nz-1))then
                in(amz-i)=dcmplx(clamb1n(k,i+1),clamb3n(k,i+1))
            endif
        enddo
        in(1)=in(1)*2.0d0
        in(nz+1)=in(nz+1)*2.0d0

        call fftw_f77_one(plan12,in,out)

        do i=0,nz
            clamb1n(k,i)=real(out(i+1),kind=8)/2.0d0
            clamb3n(k,i)=dimag(out(i+1))/2.0d0
        enddo

    enddo
    !     call fftw_f77_destroy_plan(plan)
    return
    end

    !---------------------------------------------------------
    !     this program is to get the lamb vector
    !     input:
    !          u is the velocity  array
    !          omega is vorticity	 array
    !          fr is the assistant array
    !          fi is the same array as fr
    !     output:
    !          lamb is the u multiply omega
    !---------------------------------------------------------

    subroutine getlamb(isign,u,omega,lamb,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)
    use paraandcons
    include "fftw_f77.i"
    !complex(kind=8) in, out
    complex(kind=8) :: in(4*mz),out(4*mz)

    real*8 u(0:mx,0:my,0:mz)
    real*8 omega(0:mx,0:my,0:mz)
    real*8 lamb(0:mx,0:my,0:mz)
    real*8 up(0:3*(my+1)/2/mproc-1,0:3*(mx+1)/2*mproc-1,0:mz)
    real*8 omegap(0:3*(my+1)/2/mproc-1,0:3*(mx+1)/2*mproc-1,0:mz)
    real*8 tmpu(0:2*mz)

    real*8 tmpomega(0:2*mz)
    real*8 fr(0:(mx-1)/2,0:my)
    real*8 fi(0:(mx-1)/2,0:my)

    integer i,j,k,mx,my,mz,mproc,amx,amy,amz,mn
    !      integer*8 plan
    amx=3*(my+1)/2/mproc-1
    amy=3*(mx+1)/2*mproc-1
    amz=2*mz
    mn=2*amz

    call ftrthreetotwo(1,omega,omegap,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)

    call ftrthreetotwo(1,u,up,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)

    !     call fftw_f77_create_plan(plan,mn,FFTW_FORWARD,FFTW_ESTIMATE)!**13
    do ii=0,amx
        do j=0,amy
            do k=0,amz
                if(k.le.mz)then
                    tmpu(k)=up(ii,j,k)
                    tmpomega(k)=omegap(ii,j,k)
                else
                    tmpu(k)=0.0d0
                    tmpomega(k)=0.0d0
                endif
            enddo
            !///////////////////////////////////////////////////////
            !     the following is the transform	in z -direction
            !-------------------------------------------------------

            do i=0,amz
                in(i+1)=dcmplx(tmpu(i),tmpomega(i))
                if(i.lt.(amz-1))then
                    in(mn-i)=dcmplx(tmpu(i+1),tmpomega(i+1))
                endif
            enddo
            in(1)=in(1)*2.0d0
            in(amz+1)=in(amz+1)*2.0d0
            call fftw_f77_one(plan13,in,out)

            do i=0,amz
                tmpu(i)=real(out(i+1),kind=8)*dimag(out(i+1))/4.0d0
            enddo

            do i=0,amz
                in(i+1) = dcmplx(tmpu(i),0.0d0)

                if(i.lt.(amz-1))then
                    in(mn-i)=dcmplx(tmpu(i+1),0.0d0)
                endif
            enddo
            call fftw_f77_one(plan13,in,out)

            !          call fftw_f77(plan,1,in,1,0,out,1,0)
            do i=0,amz
                tmpu(i)=real(out(i+1),kind=8)/real(amz,kind=8)
            enddo
            tmpu(0)=tmpu(0)/2.0d0
            tmpu(amz)=tmpu(amz)/2.0d0

            do k=0,mz
                up(ii,j,k)=tmpu(k)
            enddo
        enddo
    enddo
    !      call fftw_f77_destroy_plan(plan)

    call ftrthreetotwo(-1,lamb,up,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)
    return
    end

    !---------------------------------------------------------
    !     this program is to get the lamb vector
    !     input:
    !          u is the velocity  array
    !          omega is vorticity	 array
    !          fr is the assistant array
    !          fi is the same array as fr
    !     output:
    !          lamb is the u multiply omega
    !     mean pressure gradient(-1) considered here.
    !---------------------------------------------------------

    subroutine getlamb1(isign,u,omega,lamb,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)
    use paraandcons
    include "fftw_f77.i"
    !complex(kind=8) in, out
    complex(kind=8) :: in(4*mz),out(4*mz)

    real*8 u(0:mx,0:my,0:mz)
    real*8 omega(0:mx,0:my,0:mz)
    real*8 lamb(0:mx,0:my,0:mz)
    real*8 up(0:3*(my+1)/2/mproc-1,0:3*(mx+1)/2*mproc-1,0:mz)
    real*8 omegap(0:3*(my+1)/2/mproc-1,0:3*(mx+1)/2*mproc-1,0:mz)
    real*8 tmpu(0:2*mz)

    real*8 tmpomega(0:2*mz)
    real*8 fr(0:(mx-1)/2,0:my)
    real*8 fi(0:(mx-1)/2,0:my)

    integer i,j,k,mx,my,mz,mproc,amx,amy,amz,mn
    !      integer*8 plan
    amx=3*(my+1)/2/mproc-1
    amy=3*(mx+1)/2*mproc-1
    amz=2*mz
    mn=2*amz

    call ftrthreetotwo(1,omega,omegap,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)

    call ftrthreetotwo(1,u,up,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)

    !     call fftw_f77_create_plan(plan,mn,FFTW_FORWARD,FFTW_ESTIMATE)!**13
    do ii=0,amx
        do j=0,amy
            do k=0,amz
                if(k.le.mz)then
                    tmpu(k)=up(ii,j,k)
                    tmpomega(k)=omegap(ii,j,k)
                else
                    tmpu(k)=0.0d0
                    tmpomega(k)=0.0d0
                endif
            enddo
            !///////////////////////////////////////////////////////
            !     the following is the transform	in z -direction
            !-------------------------------------------------------

            do i=0,amz
                in(i+1)=dcmplx(tmpu(i),tmpomega(i))
                if(i.lt.(amz-1))then
                    in(mn-i)=dcmplx(tmpu(i+1),tmpomega(i+1))
                endif
            enddo
            in(1)=in(1)*2.0d0
            in(amz+1)=in(amz+1)*2.0d0
            call fftw_f77_one(plan13,in,out)

            do i=0,amz
                tmpu(i)=real(out(i+1),kind=8)*dimag(out(i+1))/4.0d0
            enddo

            tmpu=tmpu-re_p  ! mean pressure gradient(-1) considered here.

            do i=0,amz
                in(i+1) = dcmplx(tmpu(i),0.0d0)

                if(i.lt.(amz-1))then
                    in(mn-i)=dcmplx(tmpu(i+1),0.0d0)
                endif
            enddo
            call fftw_f77_one(plan13,in,out)

            !          call fftw_f77(plan,1,in,1,0,out,1,0)
            do i=0,amz
                tmpu(i)=real(out(i+1),kind=8)/real(amz,kind=8)
            enddo
            tmpu(0)=tmpu(0)/2.0d0
            tmpu(amz)=tmpu(amz)/2.0d0

            do k=0,mz
                up(ii,j,k)=tmpu(k)
            enddo
        enddo
    enddo
    !      call fftw_f77_destroy_plan(plan)

    call ftrthreetotwo(-1,lamb,up,id,nallgrp,  &
        fr,fi,mx,my,mz,mproc)
    return
    end

    !----------------------------------------------------------------------------
    !     this program is to intial the velocity field
    !     if isign equal -1,it represent that the field is given by possuie fluid
    !     if isign equal 1, it represent that the field is given by file.
    !     if isign equal 0,it represent that we will change the time step.
    !----------------------------------------------------------------------------

    subroutine initial(isign,x,y,z,tvn,tvni,iapha,ibeta,  &
        nallgrp,id,fr,fi)

    !      subroutine initial(isign,x,z,ibeta,iapha,tvn,tvni,  &
    !     	                  nallgrp,id,rhsr,rhsi,fr,fi,a,a0)

    use paraandcons
    include 'mpif.h'
    real*8 ibeta(0:ny),iapha(0:nhx)
    real*8 fr(0:nhx,0:ny)
    real*8 fi(0:nhx,0:ny)
    real*8 fr1(0:ny,0:ny)
    real*8 fi1(0:ny,0:ny)
    real*8 tvn(0:nz)
    real*8 tvni(0:nz)
    !      integer is(nz+1),js(nz+1)

    real*8 x(0:nx),y(0:ny),z(0:nz)          ! the location of (x,y,z)
    !real(kind=4) buffer_type_1(2,0:nz)
    !real(kind=4) buffer_type_2(0:nx,0:ny,0:nz)

    integer i,j,k,isign,nallgrp,iii
    integer istat
    character*60 inname
    real*8 starttime,endtime
    integer*4 rec
    iii=2
    !/////////////////////////////////////////////////////
    !      istep=999
    !      coef=1.d0/(real((ny+1)))
    !-----give the number of file you want to load--------

    if(id.eq.0)then
        starttime=mpi_wtime()
    endif
    !     get the initial time

    if(isign.eq.1)then

!         write(inname,921)iii,id
! 921     format('data',i8.8,'_instant','.',i3.3,'.dat')

        ! 921      format('channel'  &
        !                 ,i6,'_instant','.',i3,'.dat')


        ! 921      format('/localhome/zhaoh/channel'  &
        !                 ,i6,'_instant','.',i3,'.dat')

        write(inname,9921)iii,id,5000000
9921    format('../poiseuille_flow5/data',i8.8,'_instant','.',i3.3,'_',i8.8,'.dat')

        open(10,file=inname,status='old',form='unformatted',iostat=istat)
        !open(10,file=inname,status='old',access='stream',iostat=istat)
        
        if(istat.ne.0) write(*,*) "File open failed in proc:",id

        ! temporarily parse float to double
        ! read(10) rec
        ! if(id==0) write(*,*) rec
        ! read(10) buffer_type_1 
        ! clamb1n = buffer_type_1
        ! read(10) buffer_type_1
        ! clamb1nm1  = buffer_type_1
        ! read(10) buffer_type_1
        ! clamb3n  = buffer_type_1
        ! read(10) buffer_type_1
        ! clamb3nm1  = buffer_type_1
        ! read(10) buffer_type_1
        ! au0 = buffer_type_1
        ! read(10) buffer_type_1
        ! w0 = buffer_type_1
        ! read(10) rec
        
        ! read(10) rec
        ! read(10) buffer_type_2
        ! vnp1 = buffer_type_2
        ! read(10) buffer_type_2
        ! omegaynp1 = buffer_type_2
        ! read(10) buffer_type_2
        ! rhsvn = buffer_type_2
        ! read(10) buffer_type_2
        ! rhsomegayn = buffer_type_2
        ! read(10) buffer_type_2
        ! unp1 = buffer_type_2
        ! read(10) buffer_type_2
        ! wnp1 = buffer_type_2
        ! read(10) rec
        !!!<<<<< Original reading codes >>>>>>!!!!!
        read(10)clamb1n,clamb1nm1,clamb3n,  &
          clamb3nm1,au0,w0

        read(10)vnp1,omegaynp1,rhsvn  &
           ,rhsomegayn,unp1,wnp1
        !!!>>>>>> Original reading codes <<<<<<!!!!!

        close(10)
        write(*,*)"read inital data from",inname
        return

    endif


    ! this file contain the first step lamb,if we want to change the
    ! time step,we need them.


    !///////////////////////////////////////////////////////////////////
    if(isign.eq.-1)then
        !////////////////////////////////////////////////////
        !      we will add perturbation which satisfies incompressibility
        !///////////////////////////////////////////////////
        call input(x,y,z,nallgrp,id,fr,fi,iapha,ibeta)   ! intial the velocity field
        if(id.eq.0)then
            do k=0,nz
                au0(1,k)=unp1(0,0,k)
                au0(2,k)=unp1(1,0,k)
                w0(1,k)=wnp1(0,0,k)
                w0(2,k)=wnp1(1,0,k)
            enddo
        endif

        call getnonlinear(1,tvn,tvni,fr,fi,id,nallgrp,iapha,ibeta)

        !/////////////////////////////////////////////////////
        if(id.eq.0)then
            endtime=mpi_wtime()
            write(*,*)'initialization time:',endtime-starttime
        endif
        return
    endif


    if(isign.eq.0)then !!!!!!!!!!!!!!!!!!!!!change
        !////////////////////////////////////////////////////
        !      we will add perturbation which satisfies incompressibility
        !///////////////////////////////////////////////////
        call inputbymyself(x,y,z,nallgrp,id,fr,fi,iapha,ibeta)   ! intial the velocity field
        if(id.eq.0)then
            do k=0,nz
                au0(1,k)=unp1(0,0,k)
                au0(2,k)=unp1(1,0,k)
                w0(1,k)=wnp1(0,0,k)
                w0(2,k)=wnp1(1,0,k)
            enddo
        endif

        call getnonlinear(istepbegin,tvn,tvni,fr,fi,id,nallgrp,iapha,ibeta)

        !/////////////////////////////////////////////////////
        if(id.eq.0)then
            endtime=mpi_wtime()
            write(*,*)'initialization time:',endtime-starttime
        endif
        return
    endif
    endsubroutine
    !-----------------------------------------------------------------

    subroutine inputbymyself(x,y,z,nallgrp,id,fr,fi,iapha,ibeta) ! intial the velocity field
    use paraandcons
    include 'mpif.h'
    real*8 fr(0:nhx,0:ny)
    real*8 fi(0:nhx,0:ny)
    real*8 x(0:nx),y(0:ny),z(0:nz)          ! the location of (x,y,z)
    real*8 ibeta(0:ny),iapha(0:nhx)
    real*8 upr(0:nz),upi(0:nz)
    real*8 kmn,pi,epsilon,lx,ly             !ly-->lz
    integer i,j,k,nallgrp,iseed(nproc),irr
    real*8 vorxa(0:ny,0:nz),summyz(0:ny,0:nz)
    real*8 array(0:nx,0:ny,0:nz)
    character*50 outname1

    pi=4.0d0*atan(1.0d0)
    lx=2.0d0*pi*aphi
    ly=2.0d0*pi*beta

    clamb1n=0.0d0
    clamb3n=0.0d0
    rhsvn=0.0d0
    rhsomegayn=0.0d0


    write(outname1,22212)istepbegin,id
22212 format('velocity',i8.8,'.',i3.3,'.dat')
    open(12,file=outname1,status='replace',form='unformatted')
    do k=0,nz
        do j=0,ny
            do i=0,nx
                read(12) unp1(i,j,k)!......................plus shear velocity
            enddo
        enddo
    enddo
    do k=0,nz
        do j=0,ny
            do i=0,nx
                read(12) wnp1(i,j,k)
            enddo
        enddo
    enddo
    do k=0,nz
        do j=0,ny
            do i=0,nx
                read(12) vnp1(i,j,k)
            enddo
        enddo
    enddo
    close(12)

    do k=0,nz
        unp1(:,:,k)=unp1(:,:,k)-u0(k)
    enddo

    unp1(:,:,0)=0.0d0
    vnp1(:,:,0)=0.0d0
    wnp1(:,:,0)=0.0d0
    omegaynp1(:,:,0)=0.0d0

    unp1(:,:,nz)=0.0d0
    vnp1(:,:,nz)=0.0d0
    wnp1(:,:,nz)=0.0d0
    omegaynp1(:,:,nz)=0.0d0

    !*********************************************omegaynp1=dunp1/dy-dwnp1/dx
    array=real(0,kind=8)
    do k=0,nz
        do i=0,nx
            call deri1d2ordernon(0,ny,unp1(i,:,k),array(i,:,k),y(:))
        enddo
    enddo
    omegaynp1=array
    array=real(0,kind=8)
    do k=0,nz
        do j=0,ny
            call derivexgather(0,nx,wnp1(:,j,k),array(:,j,k),2.0d0*pi*aphi/real(nnx,kind=8),nallgrp,nproc,ierr,id)
        enddo
    enddo

    omegaynp1=omegaynp1-array



    !*********************************************
    call transform(-1,vnp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,unp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,wnp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,omegaynp1,id,nallgrp,nx,ny,nz,nproc)


    do i=0,nhx
        do j=0,ny
            kmn=iapha(i)**2+ibeta(j)**2
            call getvelocityandvorticity(i,j,iapha,  &
                ibeta,upr,upi,kmn)

        enddo
    enddo
    return

    endsubroutine
    !-----------------------------------------------------------------

    subroutine input(x,y,z,nallgrp,id,fr,fi,iapha,ibeta) ! intial the velocity field
    use paraandcons
    include 'mpif.h'
    real*8 fr(0:nhx,0:ny)
    real*8 fi(0:nhx,0:ny)
    real*8 x(0:nx),y(0:ny),z(0:nz)          ! the location of (x,y,z)
    real*8 ibeta(0:ny),iapha(0:nhx)
    real*8 upr(0:nz),upi(0:nz)
    real*8 kmn,pi,epsilon,gamma,lx,ly             !ly-->lz
    integer i,j,k,nallgrp,iseed(nproc),irr
    real*8 vorxa(0:ny,0:nz),summyz(0:ny,0:nz)
    character*50 outname1
    real*8 tmp

    pi=4.0d0*atan(1.0d0)
    ! epsilon=0.005d0
    epsilon=1.0d0
    gamma=10.0d0
    lx=2.0d0*pi*aphi
    ly=2.0d0*pi*beta

    !      lx=pi/2.0d0
    !      ly=pi/2.0d0

    !      rhsvn=0.
    !      rhsomegayn=0.
    clamb1n=0.
    clamb3n=0.
    !      unp1=0.
    !      vnp1=0.
    !      wnp1=0.
    !      omegaxnp1=0.
    !      omegaynp1=0.
    !      omegaznp1=0.

    call random_seed()
    do i=0,nx
        do j=0,ny
            do k=0,nz
                rhsvn(i,j,k)=0.0d0
                rhsomegayn(i,j,k)=0.0d0
                ! vnp1(i,j,k)=-2.0d0/3.0d0*epsilon*(1.0d0+cos(1.5d0*pi*z(k)))*(sin(2.0d0*pi*x(i)/lx) &
                !     *sin(2.0d0*pi*y(j)/ly)+sin(4.0d0*pi*x(i)/lx) &
                !     *sin(2.0d0*pi*y(j)/ly)+sin(2.0d0*pi*x(i)/lx) &
                !     *sin(4.0d0*pi*y(j)/ly))
                ! omegaynp1(i,j,k)=epsilon*2.0d0*pi*sin(1.5d0*pi*z(k))*  &
                !     (lx/ly*(cos(2.0d0*pi*x(i)/lx)*cos(2.0d0*pi*y(j)/ly) &
                !     +0.5d0*cos(4.0d0*pi*x(i)/lx)*cos(2.0d0*pi*y(j)/ly) &
                !     +2.0d0*cos(2.0d0*pi*x(i)/lx)*cos(4.0d0*pi*y(j)/ly)) &
                !     +     &
                !     ly/lx*(0.5d0*cos(2.0d0*pi*x(i)/lx)*cos(2.0d0*pi*y(j)/ly) &
                !     +cos(4.0d0*pi*x(i)/lx)*cos(2.0d0*pi*y(j)/ly) &
                !     +0.25d0*cos(2.0d0*pi*x(i)/lx)*cos(4.0d0*pi*y(j)/ly)) )
                ! unp1(i,j,k)=epsilon*lx*sin(1.5d0*pi*z(k))*(cos(2.0d0*pi*x(i)/lx) &
                !     *sin(2.0d0*pi*y(j)/ly)+0.5d0*cos(4.0d0*pi*x(i)/lx) &
                !     *sin(2.0d0*pi*y(j)/ly)+cos(2.0d0*pi*x(i)/lx) &
                !     *sin(4.0d0*pi*y(j)/ly))
                ! wnp1(i,j,k)=-epsilon*ly*sin(1.5d0*pi*z(k))*(0.5d0*sin(2.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+0.5d0*sin(4.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+0.25d0*sin(2.0d0*pi*x(i)/lx) &
                !     *cos(4.0d0*pi*y(j)/ly))
                call random_number(tmp)
                vnp1(i,j,k)=-epsilon*(1.0+cos(pi*z(k)))*(sin(2.0*pi*x(i)/lx) &
                      *sin(2.0*pi*y(j)/ly)+sin(4.0*pi*x(i)/lx) &
                      *sin(2.0*pi*y(j)/ly)+sin(2.0*pi*x(i)/lx) &
                      *sin(4.0*pi*y(j)/ly))+0.2*tmp-0.1
                call random_number(tmp)
                omegaynp1(i,j,k)=epsilon*2.0*pi*sin(pi*z(k))*  &
                      (lx/ly*(cos(2.0*pi*x(i)/lx)*cos(2.0*pi*y(j)/ly) &
                      +0.5*cos(4.0*pi*x(i)/lx)*cos(2.0*pi*y(j)/ly) &
                      +2.0*cos(2.0*pi*x(i)/lx)*cos(4.0*pi*y(j)/ly)) &
                      +     &
                      ly/lx*(0.5*cos(2.0*pi*x(i)/lx)*cos(2.0*pi*y(j)/ly) &
                      +cos(4.0*pi*x(i)/lx)*cos(2.0*pi*y(j)/ly) &
                      +0.25*cos(2.0*pi*x(i)/lx)*cos(4.0*pi*y(j)/ly)) ) &
                      +0.2*tmp-0.1
                call random_number(tmp)
                unp1(i,j,k)=epsilon*lx*sin(pi*z(k))*(cos(2.0*pi*x(i)/lx) &
                      *sin(2.0*pi*y(j)/ly)+0.5*cos(4.0*pi*x(i)/lx) &
                      *sin(2.0*pi*y(j)/ly)+cos(2.0*pi*x(i)/lx) &
                      *sin(4.0*pi*y(j)/ly))+0.2*tmp-0.1
                call random_number(tmp)
                wnp1(i,j,k)=-epsilon*ly*sin(pi*z(k))*(0.5*sin(2.0*pi*x(i)/lx) &
                      *cos(2.0*pi*y(j)/ly)+0.5*sin(4.0*pi*x(i)/lx) &
                      *cos(2.0*pi*y(j)/ly)+0.25*sin(2.0*pi*x(i)/lx) &
                      *cos(4.0*pi*y(j)/ly))+0.2*tmp-0.1

                !--------------------------------------change----------------------------------
                !------------------------------------------------------------------------------
                ! omegaxnp1(i,j,k)=( -epsilon*ly*1.5d0*pi*cos(1.5d0*pi*z(k))*(0.5d0*sin(2.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+0.5d0*sin(4.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+0.25d0*sin(2.0d0*pi*x(i)/lx) &
                !     *cos(4.0d0*pi*y(j)/ly)) )&

                !     -( 2.0d0/3.0d0*epsilon*(1.0d0+cos(1.5d0*pi*z(k)))*4.0d0*pi/ly*(0.5d0*sin(2.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+0.5d0*sin(4.0d0*pi*x(i)/lx) &
                !     *cos(2.0d0*pi*y(j)/ly)+sin(2.0d0*pi*x(i)/lx) &
                !     *cos(4.0d0*pi*y(j)/ly)) )
                !------------------------------------------------------------------------------

            enddo
        enddo
    enddo

    unp1(:,:,0)=0.0d0
    vnp1(:,:,0)=0.0d0
    wnp1(:,:,0)=0.0d0
    omegaynp1(:,:,0)=0.0d0
    ! omegaxnp1(:,:,0)=0.0d0

    unp1(:,:,nz)=0.0d0
    vnp1(:,:,nz)=0.0d0
    wnp1(:,:,nz)=0.0d0
    omegaynp1(:,:,nz)=0.0d0
    ! omegaxnp1(:,:,nz)=0.0d0
    !----------------------------------------------------------------------------------
!     do k=0,nz
!         do j=0,ny
!             vorxa(j,k)=sum(omegaxnp1(:,j,k))/(real(nx+1,kind=8)*nproc)
!         enddo
!     enddo

!     call mpi_reduce(vorxa(0,0),summyz,(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,mpi_sum,0,nallgrp,ierr)

!     if(id.eq.0)then
!         vorxa=summyz

!         open(12,file='meanu/avervorx000000.plt',status='replace')
!         write(12,*)"variables=z,y,vorx"
!         write(12,"('zone',2X,'I=',I6,2X,'J=',I6,2X,'F=point')")ny+1,nz+1
!         do k=0,nz
!             do j=0,ny
!                 write(12,*)y(j),z(k),vorxa(j,k)
!             enddo
!         enddo
!         close(12)
!     endif

!     write(outname1,2212)0,id
! 2212 format('velocity',i8.8,'.',i3.3,'.dat')
!     open(12,file=outname1,status='replace',form='unformatted')
!     do k=0,nz
!         do j=0,ny
!             do i=0,nx
!                 write(12) unp1(i,j,k)+u0(k)!......................plus shear velocity
!             enddo
!         enddo
!     enddo
!     do k=0,nz
!         do j=0,ny
!             do i=0,nx
!                 write(12) wnp1(i,j,k)
!             enddo
!         enddo
!     enddo
!     do k=0,nz
!         do j=0,ny
!             do i=0,nx
!                 write(12) vnp1(i,j,k)
!             enddo
!         enddo
!     enddo
!     close(12)


    !/////////////////////////////////////////////////////
    !................................................................output instant physical vorticity..............
!     write(outname1,2213)0,id
! 2213 format('vorticity',i8.8,'.',i3.3,'.dat')
!     open(13,file=outname1,status='replace',form='unformatted')
!     do k=0,nz
!         do j=0,ny
!             do i=0,nx
!                 write(13)omegaxnp1(i,j,k)
!             enddo
!         enddo
!     enddo

!     omegaxnp1=0.0d0
    !*********************************************
    !      if(id.eq.8)then
    !      open(2,file='test1.dat',form='formatted')
    !      do k=0,nz
    !      write(2,*)z(k),unp1(1,20,k)
    !      end do
    !      close(2)
    !      end if


    !*********************************************
    call transform(-1,vnp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,unp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,wnp1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,omegaynp1,id,nallgrp,nx,ny,nz,nproc)

    !      if(id.eq.0)then
    !      vnp1(0,0,:)=0.
    !      vnp1(1,0,:)=0.
    !      omegaynp1(0,0,:)=0.
    !      omegaynp1(1,0,:)=0.
    !      end if

    do i=0,nhx
        do j=0,ny
            kmn=iapha(i)**2+ibeta(j)**2
            call getvelocityandvorticity(i,j,iapha,  &
                ibeta,upr,upi,kmn)

        enddo
    enddo

    !--------------------------------------------------!
    !--- we add Laminar solution(Possuie flow) here ---!
    !--------------------------------------------------!
    do k=0,nz
        rhsvn(:,:,k)=gamma*(1.0-z(k)*z(k)) ! +u0(k)
        rhsomegayn(:,:,k)=gamma*2.0*z(k)   ! -du0(k)
    end do
    rhsvn(:,:,0)=0.0
    rhsomegayn(:,:,0)=0.0

    rhsvn(:,:,nz)=0.0
    rhsomegayn(:,:,nz)=0.0

    call transform(-1,rhsvn,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,rhsomegayn,id,nallgrp,nx,ny,nz,nproc)

    unp1=unp1+rhsvn
    omegaznp1=omegaznp1+rhsomegayn

    rhsvn=0.0
    rhsomegayn=0.0
    !-----------------------------------------------------!
    !--- end of consideratoin of laminar flow solution ---!
    !-----------------------------------------------------!

    energy=0.0
    do k=0,nz
        do j=0,ny
            do i=0,nx
                energy=energy+0.5*(unp1(i,j,k)**2)+0.5*(vnp1(i,j,k)**2)  &
                    +0.5*(wnp1(i,j,k)**2)
            enddo
        enddo
    enddo

    call mpi_reduce(energy,e_t,1,mpi_real,mpi_sum, &
                    0,nallgrp,ierr)
    if(id.eq.0)then
        write(*,*)'Total energy= ',e_t
    endif

    return
    end

    !-------------------------------------------------------------------
    !     this program is get the u,w,omegax,omegaz except (0,0) mode
    !     input v,and omegay
    !     output u,w,omgeax,omegaz
    !-------------------------------------------------------------------

    subroutine getvelocityandvorticity(i,j,iapha,  &
        ibeta,dvr,dvi,kmn)
    use paraandcons
    real*8 dvr(0:nz)
    real*8 dvi(0:nz)
    integer i,j,k
    real*8 ibeta(0:ny),iapha(0:nhx)
    ! --- the tranform coeefience of x mode and y mode ---

    real*8 kmn
    do k=0,nz
        dvr(k)=vnp1(2*i,j,k)
        dvi(k)=vnp1(2*i+1,j,k)
    end do

    call dz(dvr,nz)
    call dz(dvi,nz)

    if(kmn.ne.0.0d0)then

        do k=0,nz
            unp1(2*i+1,j,k)=iapha(i)*dvr(k)/kmn  &
                -ibeta(j)*omegaynp1(2*i,j,k)/kmn
            unp1(2*i,j,k)=-1.0d0*iapha(i)*dvi(k)/kmn  &
                +ibeta(j)*omegaynp1(2*i+1,j,k)/kmn
            wnp1(2*i+1,j,k)=ibeta(j)*dvr(k)/kmn  &
                +iapha(i)*omegaynp1(2*i,j,k)/kmn
            wnp1(2*i,j,k)=-1.0d0*dvi(k)*ibeta(j)/kmn  &
                -1.0d0*omegaynp1(2*i+1,j,k)*iapha(i)/kmn

        end do

        do k=0,nz
            dvr(k)=wnp1(2*i,j,k)
            dvi(k)=wnp1(2*i+1,j,k)
        enddo

        call dz(dvr,nz)
        call dz(dvi,nz)

        do k=0,nz
            omegaxnp1(2*i,j,k)=dvr(k)+ibeta(j)*vnp1(2*i+1,j,k)
            omegaxnp1(2*i+1,j,k)=dvi(k)-ibeta(j)*vnp1(2*i,j,k)
        end do

        do k=0,nz
            dvr(k)=unp1(2*i,j,k)
            dvi(k)=unp1(2*i+1,j,k)
        end do

        call dz(dvr,nz)
        call dz(dvi,nz)

        do k=0,nz
            omegaznp1(2*i,j,k)=-dvr(k)-iapha(i)*vnp1(2*i+1,j,k)
            omegaznp1(2*i+1,j,k)=-dvi(k)+iapha(i)*vnp1(2*i,j,k)
        end do

    endif
    return
    end

    !---------------------------------------------------------------------
    !     this program is to output the n time velocity and vorticity
    !     the result is fourier and chebyshev spectral coeefiance
    !     in regard to the conjugate of the coeefiance,so
    !     in 0:(nx-1)/2 store the real*8 part (nx+1)/2:nx store the image part
    !---------------------------------------------------------------------

    subroutine outputvelocityandvorticity(istep,endtime,nallgrp,id)
    use paraandcons
    include 'mpif.h'
    integer i,j,k,istep
    real*8 energy,sum,endtime,tmp
    character*50 outname1

    !///////////////////////////////////////////////////////////
    energy=0.d0
    do k=0,nz
        do j=0,ny
            do i=0,nx
                energy=energy+0.5d0*(unp1(i,j,k)**2)+0.5d0*(vnp1(i,j,k)**2)  &
                    +0.5d0*(wnp1(i,j,k)**2)
            enddo
        enddo
    enddo

    call mpi_reduce(energy,sum,1,MPI_DOUBLE_PRECISION,mpi_sum, &
        0,nallgrp,ierr)
    if(id.eq.0.and.mod(istep,50).eq.0)then !.................change
        write(*,*)"istep=",istep
        write(*,121)istep,endtime,re,sum
    endif
121 format('istep=',I9,3x,'time=',f15.3,3x,'re=',f10.1,3x,'energy=' &
        ,e16.8)

    if(mod(istep,interupt).eq.0)then
        write(outname1,1111)istep,id

1111    format('data',i6,'_velocityandvorticity','.',i3,'.dat')
        open(10,file=outname1,status='replace',form='unformatted')

        write(10)unp1,wnp1,vnp1,  &
            omegaxnp1,omegaznp1,omegaynp1

        close(10)
    endif
    return
    end

    !-----------------------------------------------------------------------
    !     this subroutine to output the result that we need to run next time
    !-----------------------------------------------------------------------
    subroutine output(istep,id)

    use paraandcons
    integer i,j,k,istep,iii !..............change
    character*50 outname1,outname2
    iii=2
    if(mod(istep,interupt).eq.0)then
        write(outname1,2111)iii,id
2111    format('data',i8.8,'_instant','.',i3.3,'.dat')

        ! 2111   format('channel'  &
        !                 ,i6,'_instant','.',i3,'.dat')

        ! 2111   format('/localhome/zhaoh/channel'  &
        !                ,i6,'_instant','.',i3,'.dat')
        open(10,file=outname1,status='replace',form='unformatted')
        write(10)clamb1n,clamb1nm1,clamb3n,  &
            clamb3nm1,au0,w0
        write(10)vnp1,omegaynp1,rhsvn  &
            ,rhsomegayn,unp1,wnp1
        close(10)

    endif


    if(mod(istep,100000).eq.0)then
        write(outname2,22211)iii,id,istep
22211   format('data',i8.8,'_instant','.',i3.3,'_',i8.8,'.dat')

        ! 2111   format('channel'  &
        !                 ,i6,'_instant','.',i3,'.dat')

        ! 2111   format('/localhome/zhaoh/channel'  &
        !                ,i6,'_instant','.',i3,'.dat')
        open(10,file=outname2,status='replace',form='unformatted')
        write(10)clamb1n,clamb1nm1,clamb3n,  &
            clamb3nm1,au0,w0
        write(10)vnp1,omegaynp1,rhsvn  &
            ,rhsomegayn,unp1,wnp1
        close(10)

    endif

    return
    end

    !-----------------------------------------------------------------
    !     this program is to output the real velocity and vorticity
    !-----------------------------------------------------------------

    subroutine realvelocityandvorticity(istep,x,y,z,id,nallgrp)
    use paraandcons
    include 'mpif.h'
    !      include 'channel1.inc'
    !************************************
    !************************************
    real*8 fr(0:nhx,0:ny)
    real*8 fi(0:nhx,0:ny)
    real*8 x(0:nx),y(0:ny),z(0:nz)
    real*8 endtime
    !      double precision ua(0:nz),sumu(0:nz)
    character*50 outname1
    integer i,j,k,istep
    if(mod(istep,interupt).eq.0)then
        write(outname1,211)istep,id
211     format('data',i6,'_velocityandvorticity','.',i3,'.dat')
        open(10,file=outname1,status='replace',form='unformatted')
        read(10)unp1,wnp1,vnp1,  &
            omegaxnp1,omegaznp1,omegaynp1
        close(10)

        call transform(1,unp1,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,vnp1,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,wnp1,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,omegaxnp1,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,omegaynp1,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,omegaznp1,id,nallgrp,nx,ny,nz,nproc)
        !*************
        !      if(id.eq.0)then
        !      open(2,file='test2.dat',form='formatted')
        !      do k=0,nz
        !      write(2,*)z(k),vnp1(1,5,k)
        !      end do
        !      close(2)
        !      end if
        !*************
        !        call statistics(istep,x,y,z,unp1,vnp1,wnp1,id,nallgrp,endtime,ierr)
        call statistics(istep,x,y,z,nstat,id,nallgrp,endtime,ierr)
        goto 123
        !///////////////////////////////////////////////////////////
        !	output the result
        !///////////////////////////////////////////////////////////
        if(istep.eq.interupt)then
            write(outname1,212)istep,id
212         format('velocity',i6,'.',i3,'.dat')
            open(12,file=outname1,status='replace')
            write(12,741)
            write(12,742) nx+1,ny+1,nz+1
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)x(i)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)y(j)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)z(k)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)unp1(i,j,k)+(1-z(k)**2)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)wnp1(i,j,k)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)vnp1(i,j,k)
                    enddo
                enddo
            enddo
            close(12)
741         format(1x,'variables="X","Y","Z","U","V","W"')
742         format(1x,'zone i=',i4,2x,'j=',i4,2X, &
                'k=',i4,2x,'f=block')

            !/////////////////////////////////////////////////////
            write(outname1,213)istep,id
213         format('vorticity',i6,'.',i3,'.dat')
            open(12,file=outname1,status='replace')
            write(12,743)
            write(12,744) nx+1,ny+1,nz+1
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)x(i)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)y(j)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)z(k)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)omegaxnp1(i,j,k)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)omegaynp1(i,j,k)
                    enddo
                enddo
            enddo
            do k=0,nz
                do j=0,ny
                    do i=0,nx
                        write(12,*)omegaznp1(i,j,k)-2*z(k)
                    enddo
                enddo
            enddo
743         format(1x,'variables="X","Y","Z","omegax","omegay","omegaz"')
744         format(1x,'zone i=',i4,2x,'j=',i4,2X,  &
                'k=',i4,2x,'f=block')

        endif
123     continue
    endif
    return
    end

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------

    subroutine statistics(istep,x,y,z,nstat,id,nallgrp,endtime,ierr)
    use paraandcons
    include 'mpif.h'
    !     real*8,allocatable,dimension(:,:,:)::omegaxnp1p,omegaynp1p,omegaznp1p
    real*8 unp1p(0:nx,0:ny,0:nz),vnp1p(0:nx,0:ny,0:nz),wnp1p(0:nx,0:ny,0:nz),&
        omegaxnp1p(0:nx,0:ny,0:nz),pressure(0:nx,0:ny,0:nz)
    real*8 ua(0:nz),va(0:nz),wa(0:nz),pa(0:nz)
    real*8 uu(0:nz),uv(0:nz),uw(0:nz)
    real*8 vv(0:nz),vw(0:nz),ww(0:nz)
    real*8 summ(0:nz),sumn(0:nz),x(0:nx),y(0:ny),z(0:nz),yplus(0:nz)
    real*8 vorxa(0:ny,0:nz),summyz(0:ny,0:nz)
    real*8 coef,ubulk,energy,endtime,sume
    character*50 outname1,chara
    integer i,j,k,istep,ierr,kaa,id,nstat,m

    real*8 up(0:nproc-1),um(0:nproc-1),dy(1:nz-1),cfly(1:nz-1)
    real*8 pi,pi2,dx,ddy,ddz,amin,amax,umax,cfl_x,cfl_y,cfl_z
    real*8,allocatable,dimension(:,:,:,:) :: ftmp


    !     if(mod(istep,interupt/100).eq.0)then
    !.................................................change....
    !	  energy=0.d0
    !      do k=0,nz
    !        do j=0,ny
    !          do i=0,nx
    !            energy=energy+0.5d0*(unp1(i,j,k)**2)+0.5d0*(vnp1(i,j,k)**2)  &
    !                   +0.5d0*(wnp1(i,j,k)**2)
    !          enddo
    !        enddo
    !      enddo
    !
    !      call mpi_reduce(energy,sume,1,MPI_DOUBLE_PRECISION,mpi_sum, &
    !                      0,nallgrp,ierr)
    !     if(id.eq.0)then
    !     if(mod(istep,interupt/10).eq.0)then
    !            write(*,*)"istep=",istep
    !            write(*,121)istep,endtime,re,sume
    !            121  format('istep=',I8,3x,'time=',f15.3,3x,'re=',f10.1,3x,'engery=',e16.8)
    !     endif
    !     open(999,file='energy.plt',status='replace',form='formatted',position='append')
    !          WRITE(999,*) istep, sume
    !!		        call fseek(999,0,2)
    !			     write(999,*)"istep=",istep
    !			     write(999,*)"energy=",sume
    !     close(999)
    !
    !     endif !if (id ==0)


    !---------------------------------------------------------------------------------------------
    if(mod(istep,outputstep).eq.0)then
        !.......................................................................
        !	    allocate(omegaxnp1p(0:nx,0:ny,0:nz),omegaynp1p(0:nx,0:ny,0:nz),omegaznp1p(0:nx,0:ny,0:nz))
        unp1p=unp1
        vnp1p=vnp1
        wnp1p=wnp1

        !        omegaxnp1p=omegaxnp1
        !        omegaynp1p=omegaynp1
        !        omegaznp1p=omegaznp1



        call transform(1,unp1p,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,vnp1p,id,nallgrp,nx,ny,nz,nproc)
        call transform(1,wnp1p,id,nallgrp,nx,ny,nz,nproc)
        !        call transform(1,omegaxnp1p,id,nallgrp,nx,ny,nz,nproc)
        !		call transform(1,omegaynp1p,id,nallgrp,nx,ny,nz,nproc)
        !		call transform(1,omegaznp1p,id,nallgrp,nx,ny,nz,nproc)
        call mpi_barrier(nallgrp,ierr)
        !////////////////////////////////////change add by myself//////////////////////////////////////////
        if(istep==istepbegin)then
            kaa=0
        else
            kaa=1
        endif
        !---*** remove the part of solving pressure
        !		pressure=0.0d0
        !        call instantp(x,y,z,id,pressure,nallgrp,kaa,ierr,unp1p,vnp1p,wnp1p)
        !//////////////////////////////////////////////////////////////////////////////////////////////////
        !................................................................output instant physical velocity..............

        allocate(ftmp(0:nx,0:ny,0:nz,3))
        write(outname1,212)istep,id
212     format('3Ddata/velocity',i8.8,'.',i3.3,'.dat')
        open(12,file=outname1,status='replace',form='unformatted')
        do k=0,nz
            do j=0,ny
                do i=0,nx
                    ftmp(i,j,k,1) = unp1p(i,j,k)+u0(k) !......................plus shear velocity
                    ftmp(i,j,k,2) = wnp1p(i,j,k)
                    ftmp(i,j,k,3) = vnp1p(i,j,k)
                enddo
            enddo
        enddo
        !write(12) ((((sngl(ftmp(i,j,k,m)),i=0,nx),j=0,ny),k=0,nz),m=1,3)
        write(12) ((((ftmp(i,j,k,m),i=0,nx),j=0,ny),k=0,nz),m=1,3)
        close(12)

        energy = 0.d0
        do k=0,nz
            do j=0,ny
                do i=0,nx
                    energy=energy + 0.5d0*(ftmp(i,j,k,1)**2 + ftmp(i,j,k,2)**2  &
                        + ftmp(i,j,k,3)**2)
                enddo
            enddo
        enddo

        call mpi_reduce(energy,sume,1,MPI_DOUBLE_PRECISION,mpi_sum, &
            0,nallgrp,ierr)
        if(id.eq.0)then
            if(mod(istep,interupt/10).eq.0)then
                write(*,*)"istep=",istep
                write(*,121)istep,endtime,sume
121             format('istep=',I8,3x,'time=',f15.3,3x,'engery=',e16.8)
            endif
            ! open(999,file='energy.plt',status='replace',form='formatted',position='append')
            open(999,file='energy.plt',form='formatted',position='append')
            WRITE(999,*) istep, sume
            close(999)
            !
        endif !if (id ==0)


        deallocate(ftmp)

        !...................................................output instant physical vorticity..............
        !          write(outname1,22213)istep,id
        ! 22213      format('3Ddata/vorticity',i8.8,'.',i3.3,'.dat')
        !          open(153,file=outname1,status='replace',form='unformatted')
        !          do k=0,nz
        !            do j=0,ny
        !              do i=0,nx
        !                write(153)omegaxnp1p(i,j,k)
        !             enddo
        !            enddo
        !          enddo
        !         do k=0,nz
        !           do j=0,ny
        !             do i=0,nx
        !               write(13)omegaynp1p(i,j,k)
        !             enddo
        !           enddo
        !         enddo
        !         do k=0,nz
        !           do j=0,ny
        !             do i=0,nx
        !               write(13)omegaznp1p(i,j,k)+0.5d0
        !             enddo
        !           enddo
        !         enddo
        !		close(153)

        !          write(outname1,32213)istep,id
        !32213      format('3Ddata/pressure',i8.8,'.',i3.3,'.dat')
        !          open(153,file=outname1,status='replace',form='unformatted')
        !          do k=0,nz
        !            do j=0,ny
        !              do i=0,nx
        !                write(153)pressure(i,j,k)
        !             enddo
        !            enddo
        !          enddo
        !	     close(153)

        !		deallocate(omegaxnp1p,omegaynp1p,omegaznp1p)
        !----------------------------------------------------------------------------------
        !		do k=0,nz
        !		  do j=0,ny
        !		   vorxa(j,k)=sum(omegaxnp1p(:,j,k))/(real(nx+1,kind=8)*nproc)
        !		  enddo
        !		enddo
        !
        !        call mpi_reduce(vorxa(0,0),summyz,(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,mpi_sum,0,nallgrp,ierr)
        !
        !        if(id.eq.0)then
        !          vorxa=summyz
        !
        !!		write(chara,*)istep
        !!        outname1='meanu/avervorx'//trim(adjustl(chara))//'.plt'
        !          write(outname1,22113)istep
        !   22113     format('meanu/avervorx',i8.8,'.','plt')
        !          open(12,file=outname1,status='replace')
        ! 		  write(12,*)"variables=z,y,vorx"
        !          write(12,"('zone',2X,'I=',I6,2X,'J=',I6,2X,'F=point')")ny+1,nz+1
        !           do k=0,nz
        !            do j=0,ny
        !              write(12,*)y(j),z(k),vorxa(j,k)
        !			enddo
        !          enddo
        !          close(12)
        !
        !        end if

        !---------------------------------------------------------------------
        ! go to 20101
        if(mod(istep,interupt).eq.0)then
        ! if(mod(istep,10) .eq. 0)then
            coef=1.0d0/(real(nx+1,kind=8)*nproc*real(ny+1,kind=8))

            !.......................change CFL.............................
            pi=4.0d0*atan(1.0d0)
            pi2=2.0d0*pi
            dx=aphi*pi2/real(nnx,kind=8)
            ddz=beta*pi2/real(nny,kind=8)
            ddy=pi/real(nnz,kind=8)
            do k=1,nz/2
                dy(k)=-(cos(real(k,kind=8)*ddy)-cos(real(k-1,kind=8)*ddy))
                dy(nz-k)=dy(k)
            end do
            um=0.0d0
            up=0.0d0
            um(id)=minval(unp1p)    ! get max and min value of u from each process.
            up(id)=maxval(unp1p)
            call mpi_barrier(nallgrp,ierr)

            do i=0,nproc-1
                call mpi_bcast(um(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
                call mpi_bcast(up(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
            end do

            amin= minval(um)
            amax= maxval(up)
            umax= max(abs(amin),abs(amax))
            cfl_x=umax*dt/dx

            um=0.0d0
            up=0.0d0
            um(id)=minval(wnp1p)    ! get max and min value of w from each process.
            up(id)=maxval(wnp1p)
            call mpi_barrier(nallgrp,ierr)

            do i=0,nproc-1
                call mpi_bcast(um(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
                call mpi_bcast(up(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
            end do

            amin= minval(um)
            amax= maxval(up)
            umax= max(abs(amin),abs(amax))
            cfl_z=umax*dt/ddz

            do k=1,nz-1

                um=0.0d0
                up=0.0d0
                um(id)=minval(vnp1p(:,:,k))
                up(id)=maxval(vnp1p(:,:,k))
                call mpi_barrier(nallgrp,ierr)

                do i=0,nproc-1
                    call mpi_bcast(um(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
                    call mpi_bcast(up(i),1,MPI_DOUBLE_PRECISION, i,nallgrp,ierr)
                end do

                amin= minval(um)
                amax= maxval(up)
                umax= max(abs(amin),abs(amax))
                cfly(k)=umax*dt/dy(k)

            end do

            cfl_y=maxval(cfly)

            if(id.eq.0)then
                write(*,*)'CFL_X= ',cfl_x,'CFL_Y= ',cfl_y,'CFL_Z= ',cfl_z
            end if

            !............................................................
            !.......................................................................
            !  compute layer averaged u, v, w...

            do k=0,nz
                ua(k)=0.0d0
                va(k)=0.0d0
                wa(k)=0.0d0
                !		  pa(k)=0.0d0
                do j=0,ny
                    do i=0,nx
                        ua(k)=ua(k)+unp1p(i,j,k)
                        va(k)=va(k)+vnp1p(i,j,k)
                        wa(k)=wa(k)+wnp1p(i,j,k)
                        !			  pa(k)=pa(k)+pressure(i,j,k)
                    enddo
                enddo
                ua(k)=ua(k)*coef    
                va(k)=va(k)*coef
                wa(k)=wa(k)*coef
                !		  pa(k)=pa(k)*coef
            enddo

            call mpi_reduce(ua,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)

            if(id.eq.0)then
                do k=0,nz
                    ua(k)=summ(k)
                enddo
            end if
            call mpi_reduce(va,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    va(k)=summ(k)
                enddo
            endif
            call mpi_reduce(wa,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    wa(k)=summ(k)
                enddo
            endif

            !        call mpi_reduce(pa,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
            !                        0,nallgrp,ierr)
            !        if(id.eq.0)then
            !          do k=0,nz
            !            pa(k)=summ(k)
            !          enddo
            !        endif

            call mpi_bcast(ua,nz+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)
            call mpi_bcast(va,nz+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)
            call mpi_bcast(wa,nz+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)
            !        call mpi_bcast(pa,nz+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)

            !		if(id==0)then
            !		   write(*,*)"p(0)=",pa(0)
            !		endif
            ! =======^ ^ ^========= average =========^ ^ ^=============


            !//////////////  compute Reynolds-stress /////////////////
            do k=0,nz
                uu(k)=0.0d0
                uv(k)=0.0d0
                uw(k)=0.0d0
                vv(k)=0.0d0
                ww(k)=0.0d0
                vw(k)=0.0d0
                do j=0,ny
                    do i=0,nx
                        uu(k)=uu(k)+(unp1p(i,j,k)-ua(k))**2
                        vv(k)=vv(k)+(vnp1p(i,j,k)-va(k))**2
                        ww(k)=ww(k)+(wnp1p(i,j,k)-wa(k))**2
                        uv(k)=uv(k)+(unp1p(i,j,k)-ua(k))*(vnp1p(i,j,k)-va(k))
                        uw(k)=uw(k)+(unp1p(i,j,k)-ua(k))*(wnp1p(i,j,k)-wa(k))
                        vw(k)=vw(k)+(vnp1p(i,j,k)-va(k))*(wnp1p(i,j,k)-wa(k))
                    enddo
                enddo
                uu(k)=uu(k)*coef
                uv(k)=uv(k)*coef
                uw(k)=uw(k)*coef
                vv(k)=vv(k)*coef
                ww(k)=ww(k)*coef
                vw(k)=vw(k)*coef
            enddo

            call mpi_reduce(uu,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    uu(k)=summ(k)
                enddo
            endif
            call mpi_reduce(uv,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    uv(k)=summ(k)
                enddo
            endif
            call mpi_reduce(uw,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    uw(k)=summ(k)
                enddo
            endif
            call mpi_reduce(vv,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    vv(k)=summ(k)
                enddo
            endif
            call mpi_reduce(vw,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    vw(k)=summ(k)
                enddo
            endif
            call mpi_reduce(ww,summ,nz+1,MPI_DOUBLE_PRECISION,mpi_sum,  &
                0,nallgrp,ierr)
            if(id.eq.0)then
                do k=0,nz
                    ww(k)=summ(k)
                enddo
            endif

            ! ///////////////  Reynolds-Stress(z) evaluating //////////////

            ! ///////  output average velocity and reynolds stress ////////
            if(id==0)then
                write(outname1,2113)istep
2113            format('meanu/average_velocity',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
!                 open(12,file=outname1,position='append')
!                 write(12,2114) istep,nz+1
! 2114            format('ZONE T="istep=',I6,'" I=',I3,' DATAPACKING=POINT' )
                do k=0,nz
                    write(12,*)z(k),ua(k)+u0(k) !,pa(k)      !+1-z(k)**2..................change
                enddo
                close(12)
                !write(outname1,2117)istep
                !2117     format('uv',i6,'.','dat')
                summ=uv
                sumn=ua+u0   !+1.0d0-z**2................change

                call chebyshev(-1,ua,nz)
                call dz(ua,nz)
                call chebyshev(1,ua,nz)
                do k=0,nz
                    uv(k)=-uv(k)+ua(k)/re+du0(k)/re    !-2/re*z(k).............change
                enddo

                write(outname1,21113)istep

21113           format('stress/Reynolds_shear_stress',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),-summ(k)/uv(nz)
                enddo
                close(12)

                write(outname1,21114)istep
21114           format('stress/viscous_stress',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),(ua(k)/re+du0(k)/re)/uv(nz)
                enddo
                close(12)

                write(outname1,2118)istep
2118            format('stress/total_shear_stress',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),uv(k)/uv(nz)
                enddo
                close(12)

                write(outname1,21119)istep
21119           format('rms/urms',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),dsqrt(uu(k))/dsqrt(uv(nz))
                enddo
                close(12)
                write(outname1,2115)istep
2115            format('rms/vrms',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),dsqrt(vv(k))/dsqrt(uv(nz))
                enddo
                close(12)

                write(outname1,2116)istep
2116            format('rms/wrms',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz
                    write(12,*)z(k),dsqrt(ww(k))/dsqrt(uv(nz))
                enddo
                close(12)

                do k=0,nz/2
                    yplus(k)=re*((1.0d0-z(k))+(1.0d0+z(nz-k)))/2.0d0*dsqrt(uv(nz))    !y+
                end do


                write(outname1,21116)istep
21116           format('uplus/uplus',i8.8,'.','plt')
                open(12,file=outname1,status='replace')
                ! open(12,file=outname1,position='append')
                ! write(12,2114) istep,nz+1
                do k=0,nz/2
                    ! write(12,*)(yplus(k)+yplus(nz-k))/2.0d0, &
                    write(12,*)yplus(k), &
                        (sumn(k)+sumn(nz-k))/2.0d0/dsqrt(uv(nz))  !u+....................change
                enddo
                close(12)

                ubulk=0.0d0
                do k=nz/2,nz-1
                    ubulk=ubulk+sumn(k)*abs(z(k+1)-z(k))   !bulk velocity
                end do

                write(*,*)'istep=',istep
                write(*,*)'wall shear velocity: ',dsqrt(abs(uv(nz))),dsqrt(abs(uv(0)))
                write(*,*)'Mean pressure gradient[d<p>/dx]: ',-re_p, &
                    '& The shear-stress gradient[d\tau/dy=-\tau_w/\delta]: ',-uv(nz),uv(0)
                write(*,*)'Reynolds number R_t: ',dsqrt(abs(uv(nz)))*re,dsqrt(abs(uv(0)))*re
                write(*,*)'Reynolds number R_c: ',sumn(nz/2)*re
                write(*,*)'Reynolds number R_b: ',2.0d0*ubulk*re

                open(999,file='diagnostics.plt',form='formatted',position='append')
                write(999,*) istep, dsqrt(abs(uv(nz))), dsqrt(abs(uv(0))), -re_p, -uv(nz), uv(0), &
                    dsqrt(abs(uv(nz)))*re, dsqrt(abs(uv(0)))*re, sumn(nz/2)*re, 2.0d0*ubulk*re
                close(999)
            endif   !if id==0

        endif   !	  if(mod(istep,interupt).eq.0)then

! 20101   continue
    endif   !	  if(mod(istep,outputstep).eq.0)then

    return
    end


    !-----------------------------------------------------------------------
    !***********************************************************************
    subroutine disturb(istep,y,fr1,fi1,fr2,fi2,mx,my,nallgrp,id,mproc)
    include 'mpif.h'
    real*8 x(0:mx),y(0:my)
    real*8 fr1(0:mx,0:my),fi1(0:mx,0:my)
    real*8 fr2(0:mx,0:my),fi2(0:mx,0:my),pi
    integer i,j,nallgrp
    integer iseed(mproc),ierr
    real*8 t2,t1
    pi=4.0d0*atan(1.0d0)

    if(id.eq.0)then
        call srand(istep)
        do i=1,mproc
            iseed(i) = irand()
        end do
    end if
    !      call mpi_barrier(nallprc,ierr)
    !      call mpi_bcast(iseed,mproc,mpi_integer,0,nallgrp,ierr)
    irr = iseed(id+1)
    call srand(irr)

    do j=0,mx
        x(j)=2.0d0*pi*real(j,kind=8)/real(mx+1,kind=8)
    end do
    do i=0,mx
        do j=0,my
            fr1(i,j)=0.0d0
            fr2(i,j)=0.0d0
            !          fr1(i,j)=cos(x(i))
            !	   fr2(i,j)=cos(y(j))

        enddo
    enddo
    if(id.eq.0)then
        do k=0,mx
            do i=0,mx
                do j=0,mx
                    t1=rand()
                    t2=rand()
                    fr1(i,j)=fr1(i,j)+(cos(x(i)*(k+2.0d0*t1))  &
                        +cos(y(i)*(k-3.0d0*t2)))  &
                        +(sin(x(i)*(k+2.0d0*t1))  &
                        +sin(y(i)*(k-3.0d0*t2)))*mproc*0.01d0
                    fr2(i,j)=fr2(i,j)+(cos(x(i)*(k+2.0d0*t1))  &
                        -cos(y(i)*(k-3.0d0*t2)))  &
                        +(sin(x(i)*(k+2.0d0*t1))  &
                        -sin(y(i)*(k-3.0d0*t2)))*mproc*0.01d0
                enddo
            enddo
        enddo

        !        call fftx(1,fr1,fi1,mx,my)
        !        call ffty(1,fr1,fi1,1.d0,mx,my)
        !        call fftx(1,fr2,fi2,mx,my)
        !        call ffty(1,fr2,fi2,1.d0,mx,my)
    endif
    return
    end

    !-------------------------------------------------------------
    ! Compute chebyshev poly.s 
    ! Y0 : coordinates Y
    ! U0, du0, ddu0 : velocity, derivatives of velocity.
    ! T0(J,L) = cos(PI*I*L/N2)
    ! T1(J,L) = L*sin(PI*J*L/N2)/sin(PI*J/N2)
    ! T2(J,L) = L*(sin(PI*J*L/N2)*cos(PI*J/N2)-L*cos(PI*J*L/N2)/sin(PI*J/N2))/S**3
    ! T4(J,L) = ...
    !-------------------------------------------------------------

    SUBROUTINE SCHEB(t0,t1,t2,t4,u0,du0,ddu0,l2)
    IMPLICIT real*8(B-H,O-Z)
    PARAMETER(Z2=2.0d0,Z3=3.0d0,Z4=4.0d0,Z6=6.0d0,Z9=9.0d0,ZB=11.0d0)
    PARAMETER(L1=256)
    real*8 T0(0:L2,0:L2),T1(0:L2,0:L2),T2(0:L2,0:L2),T4(0:L2,0:L2)
    real*8 y0(0:l2),u0(0:l2),du0(0:l2),ddu0(0:l2)
    integer n2

    n2=l2
    !      H0=3.141592653589792d0/real(N2,kind=8)
    H0=4.0d0*atan(1.0d0)/real(N2,kind=8) !H0 = PI/N2
    ! U0(N2)=real(0,kind=8)!...................change
    ! U0(0)=real(1,kind=8)!....................change
    ! dU0(N2)=0.5d0 !...................change
    ! dU0(0)=0.5d0 !....................change
    ! ddU0(N2)=real(0,kind=8)!...................change
    ! ddU0(0)=real(0,kind=8)!....................change
    U0(N2)=real(0,kind=8)!...................change
    U0(0)=real(0,kind=8)!....................change
    dU0(N2)=0.0d0 !...................change
    dU0(0)=0.0d0 !....................change
    ddU0(N2)=real(0,kind=8)!...................change
    ddU0(0)=real(0,kind=8)!....................change

    DO L=0,N2
        T0(0,L)=real(1,kind=8)
        T0(N2,L)=real((-1)**L)
        T1(0,L)=real(L,kind=8)**2
        T1(N2,L)=(-1)**(L-1)*real(L,kind=8)**2
        T2(0,L)=real(L,kind=8)**2*(real(L,kind=8)**2-real(1,kind=8))/real(3,kind=8)
        T2(N2,L)=(-1)**L*(real(L,kind=8)**2)*(real(L,kind=8)**2-real(1,kind=8))/real(3,kind=8)
    ENDDO

    DO J=1,N2-1
        Q=H0*J    ! Q = PI*J/N2
        C1=cos(Q)
        C2=C1*C1
        S1=sin(Q)
        S2=S1*S1
        DS1=1.0d0/S1    !DS_ = 1/(S_)
        DS3=DS1/S2
        DS5=DS3/S2
        DS7=DS5/S2
        Y0(J)=C1
        ! U0(J)=0.5d0*(1.0d0+Y0(J))!.................change
        ! du0(j)=0.5d0 !.................change
        ! ddu0(j)=real(0,kind=8) !.................change
        U0(J)=0.0d0!.................change
        du0(j)=0.0d0 !.................change
        ddu0(j)=real(0,kind=8) !.................change
        !      U1(J)=-2.D0*Y0(J)*HTH
        !      U2(J)=-2.D0*HTH

        DO L=0,N2
            BJ=real(L,kind=8)
            QJ=Q*BJ
            SJ=sin(QJ)
            CJ=cos(QJ)
            BS=BJ*BJ*S2
            BC=BJ*CJ*S1
            T0(J,L)=CJ
            T1(J,L)=BJ*SJ*DS1
            T2(J,L)=BJ*(SJ*C1-BC)*DS3
            !	 T3(J,L)=BJ*(SJ*(1.D0+Z2*C2-BS)-Z3*BC*C1)*DS5
            T4(J,L)=BJ*(SJ*C1*(Z9+Z6*(C2-BS))-BC*(Z4+ZB*C2-BS))*DS7
        END DO
    END DO


    RETURN
    END


    !----------------------------------------------------------------
    ! compute Ampn in equation 2.21
    !----------------------------------------------------------------
    subroutine get_coef_v(a,n,kmn,alpha)
    use paraandcons
    complex(kind=8) a(0:n,0:n)
    real*8 alpha
    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 kmn,tmpr,tmpi
    integer i,j,k
    !.....................................................
    !      call scheb(t0,t2,t4,u0,du0,ddu0,n)
    do i=4,n
        !///////////////////////////////////////////////////////////
        !     the following give the coefence to the matrix A(i,j) in (2.21)
        do j=0,n
            tmpr=-t4(i-2,j)*dt*0.5d0/re+(1+kmn*dt/re)*t2(i-2,j)  &
                -kmn*(1+kmn*dt*0.5d0/re)*t0(i-2,j)
            ! tmpi=alpha*dt*u0(i-2)*0.5d0*t2(i-2,j)-  &
            !     (kmn*alpha*dt*0.5d0*u0(i-2)-alpha*dt*0.5d0*ddu0(i-2))*t0(i-2,j)
            tmpi=0.0
            a(i,j)=dcmplx(tmpr,tmpi)
        enddo
    enddo

    do j=0,n
        a(0,j)=(1.d0,0.d0)
        a(1,j)=dcmplx(real((-1)**(j)),0.d0)
        a(2,j)=dcmplx(real(j,kind=8)**2,0.d0)
        a(3,j)=cmplx((-1)**(j-1)*real(j,kind=8)**2,0.d0)
        !	 a(1,j)=(1.d0,0.d0)
        !	 a(2,j)=((-1)**(j-1),0.d0)
        !        a(3,j)=(real(j-1,kind=8)**2,0.d0)
        !        a(4,j)=((-1)**(j-2)*real(j-1,kind=8)**2,0.d0)
    enddo
    return
    end

    !--------------------------------------------------------------

    subroutine get_coef_rhs(rhsr,rhsi,tvn,tvni,rhs,n,kmn,alpha  &
        ,bdv,bdvi)
    use paraandcons
    complex(kind=8) rhs(0:n)
    real*8 rhsr(0:n),rhsi(0:n)
    real*8 tvn(0:n),tvni(0:n)
    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 kmn,bdv(4),bdvi(4),tmpr,tmpi,alpha
    integer i,j,k
    !.....................................................
    !      call scheb(t0,t2,t4,u0,du0,ddu0,n)
    do i=4,n
        rhs(i)=(0.0d0,0.0d0)
        do j=0,n
            !       tmpr=(t4(i-2,j)*dt*0.5d0/re+(1-kmn*dt/re)* &
            !       t2(i-2,j)-kmn*(1-kmn*dt*0.5d0/re)*t0(i-2,j))*tvn(j)-(-alpha*dt* &
            !       u0(i-2)*0.5d0*t2(i-2,j)+(kmn*alpha*dt*u0(i-2)*0.5d0-alpha*dt)* &
            !       t0(i-2,j))*tvni(j)
            !   	tmpi=(t4(i-2,j)*dt*0.5d0/re+(1-kmn*dt/re)* &
            !       t2(i-2,j)-kmn*(1-kmn*dt*0.5d0/re)*t0(i-2,j))*tvni(j)+(-alpha*dt* &
            !       u0(i-2)*0.5d0*t2(i-2,j)+(kmn*alpha*dt*u0(i-2)*0.5d0-alpha*dt)* &
            !       t0(i-2,j))*tvn(j)

            tmpr=(t4(i-2,j)*dt*0.5d0/re+(1-kmn*dt/re)*  &
                t2(i-2,j)-kmn*(1-kmn*dt*0.5d0/re)*t0(i-2,j))
            ! tmpi= (-alpha*dt*u0(i-2)*0.5d0*t2(i-2,j)+  &
            !     (kmn*alpha*dt*u0(i-2)*0.5d0-alpha*dt*0.5d0*ddu0(i-2))*t0(i-2,j))
            tmpi=0.0
            rhs(i)=rhs(i)+dcmplx(tmpr,tmpi)*dcmplx(tvn(j),tvni(j))
            !///////////////////////////////////////////////////////////
            !     the following give the coefence to the matrix
            !///////////////////////////////////////////////////////////
        enddo
        rhs(i)=dcmplx(rhsr(i-2),rhsi(i-2))+rhs(i)
    enddo
    rhs(0)=dcmplx(bdv(1),bdvi(1))
    rhs(1)=dcmplx(bdv(2),bdvi(2))
    rhs(2)=dcmplx(bdv(3),bdvi(3))
    rhs(3)=dcmplx(bdv(4),bdvi(4))
    return
    end

    !------------------------------------------------------------
    !------------------------------------------------------------
    !     this program is to get the coeffience of the omegay
    !------------------------------------------------------------

    subroutine get_coef_omegay(a,n,kmn,alpha)
    use paraandcons
    complex(kind=8) a(0:n,0:n)
    real*8 alpha
    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 kmn,tmpr,tmpi
    integer i,j
    !      call scheb(t0,t2,t4,u0,du0,ddu0,n)
    ! i : number of points in normal direction
    ! j : order of chebyshev poly.
    do i=2,n
        do j=0,n
            tmpr=-t2(i-1,j)*dt*0.5d0/re+(1+kmn*dt*0.5d0/re)*t0(i-1,j)
            ! tmpi=alpha*dt*0.5d0*u0(i-1)*t0(i-1,j)
            tmpi=0.0
            a(i,j)=dcmplx(tmpr,tmpi)
        enddo
    enddo
    do j=0,n
        a(0,j)=dcmplx(1.0d0,0.0d0)
        a(1,j)=dcmplx(real((-1)**(j)),0.0d0)
    enddo
    return
    end

    !------------------------------------------------------------
    subroutine get_rhs_omegay(rhsr,rhsi,rhs,tvn,tvni,tvnm1  &
        ,tvnm1i,tomegayn,tomegayni,n,kmn,alpha,betav, &
        bdomega,bdomegai)
    use paraandcons
    complex(kind=8) rhs(0:n)
    real*8 rhsr(0:n),rhsi(0:n)
    real*8 tvn(0:n),tvni(0:n)
    real*8 tvnm1(0:n),tvnm1i(0:n)
    real*8 tomegayn(0:n),tomegayni(0:n)
    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 kmn,bdomega(2),bdomegai(2)
    real*8 tmpr,tmpi,alpha,betav
    integer i,j,k
    !      call scheb(t0,t2,t4,u0,du0,ddu0,n)
    do i=2,n
        rhs(i)=(0.0d0,0.0d0)
        do j=0,n
            tmpr=t2(i-1,j)*dt*0.5d0/re+(1-kmn*dt*0.5d0/re)*t0(i-1,j)
            ! tmpi=-alpha*dt*0.5d0*u0(i-1)*t0(i-1,j)
            tmpi=0.0
            rhs(i)=rhs(i)+dcmplx(tmpr,tmpi)*dcmplx(tomegayn(j)  &
                ,tomegayni(j))
            ! tmpi=-0.5d0*du0(i-1)*dt*t0(i-1,j)*betav
            ! rhs(i)=rhs(i)+dcmplx(0.0d0,tmpi)*dcmplx(tvn(j)+tvnm1(j),  &
            !     tvni(j)+tvnm1i(j))
        enddo
        rhs(i)=rhs(i)+dcmplx(rhsr(i-1),rhsi(i-1))
    enddo
    rhs(0)=dcmplx(bdomega(1),bdomegai(1))
    rhs(1)=dcmplx(bdomega(2),bdomegai(2))
    return
    end

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine get_coef_u0(a,n)
    use paraandcons
    complex(kind=8) a(0:n,0:n)
    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 tmpr
    integer i,j
    !      call scheb(t0,t2,t4,u0,du0,ddu0,n)
    do i=2,n
        do j=0,n
            tmpr=-t2(i-1,j)*dt*0.5d0/re+t0(i-1,j)
            a(i,j)=dcmplx(tmpr,0.0d0)
        enddo
    enddo
    do j=0,n
        a(0,j)=(1.0d0,0.0d0)
        a(1,j)=dcmplx(real((-1)**(j)),0.0d0)
    enddo
    return
    end

    !--------------------------------------------------------------------

    subroutine get_rhs_u0(rhsr,rhsi,rhs,tvn,tvni,n,bdu,bdui)
    use paraandcons
    complex(kind=8) rhs(0:n)
    real*8 rhsr(0:n),rhsi(0:n)
    real*8 tvn(0:n),tvni(0:n),bdu(2),bdui(2)

    !      double precision T0(0:n,0:n),T2(0:n,0:n),T4(0:n,0:n)
    real*8 tmpr
    do i=2,n
        rhs(i)=(0.0d0,0.0d0)
        do j=0,n
            tmpr=t0(i-1,j)+dt*0.5d0/re*t2(i-1,j)
            rhs(i)=rhs(i)+dcmplx(tmpr,0)*dcmplx(tvn(j),tvni(j))
        enddo
        rhs(i)=rhs(i)+dcmplx(rhsr(i-1),rhsi(i-1))
    enddo
    rhs(0)=dcmplx(bdu(1),bdui(1))
    rhs(1)=dcmplx(bdu(2),bdui(2))
    return
    end

    !--------------------------------------------------------------------
    !------- Gauss Full delimination to Solve the system of eqn's--------
    !------- solving equation (aa * X = af), save result X in af --------
    !--------------------------------------------------------------------
    subroutine gsij(ka,k1,k3,aa,af,mx,my,mz,id)
    use paraandcons
    INTEGER ::  ik(0:2*(nhx+1),0:nz,0:ny),jk(0:2*(nhx+1),0:nz,0:ny)
    INTEGER ::  ki(0:mz),kj(0:mz)
    complex(kind=8) aa(0:mz,0:mz),af(0:mz)
    real*8 dd,dq
    complex(kind=8) add,as
    integer k1,k3,im,jm,ka,mx,my,mz,id,i,j,k
    integer :: ii,li,lj
    if(ka.eq.0) then
        ki=1
        kj=1
    end if
    do 10 k=0,mz
        if(ka.eq.0) then
            dd=0.0d0
            ! ---- find the max value aa(i,j)=dd, max's i,j = im, ij. ----
            do 15 i=0,mz
                if(ki(i).eq.1) then
                    do 16 j=0,mz
                        dq=dsqrt((real(aa(i,j),kind=8))**2+(dimag(aa(i,j)))**2)
                        if(kj(j).eq.1.and.dq.gt.dd) then
                            dd=dq
                            im=i
                            jm=j
                        end if
16                  continue
                end if
15          continue
            ! -^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^-
            ik(k1,k,k3)=im ! save i,j of max values for later step use.
            jk(k1,k,k3)=jm
            kj(jm)=0    ! mark the 'used' max value's i,j.
            ki(im)=0
            if(k.eq.mz) goto 10
            add=1.0d0/aa(im,jm)
            do 25 i=0,mz
                if(ki(i).eq.1) then
                    as=aa(i,jm)*add
                    do 20 j=0,mz
                        if(kj(j).eq.1) aa(i,j)=aa(i,j)-as*aa(im,j)
            ! Gauss elimination, use aa(im,jm) to eliminate aa(i,jm)
20                  continue
                    af(i)=af(i)-as*af(im)
                end if
25          continue
        else
            im=ik(k1,k,k3)
            jm=jk(k1,k,k3)
            if(k.eq.mz) goto 10
            add=1.0d0/aa(im,jm)
            do 65 li=k+1,mz
                i=ik(k1,li,k3)
                as=aa(i,jm)*add
                do 60 lj=k+1,mz
                    j=jk(k1,lj,k3)
                    aa(i,j)=aa(i,j)-as*aa(im,j)
60              continue
                af(i)=af(i)-as*af(im)
65          continue
        end if

10  continue
    aa(im,jm)=af(im)/aa(im,jm)
    do 30 k=mz-1,0,-1
        i=ik(k1,k,k3)
        j=jk(k1,k,k3)
        aa(im,j)=af(i)
        do 40 lj=k+1,mz
            jj=jk(k1,lj,k3)
            aa(im,j)=aa(im,j)-aa(i,jj)*aa(im,jj)
40      continue
        aa(im,j)=aa(im,j)/aa(i,j)
30  continue
    do 50 j=0,mz
        af(j)=aa(im,j)
50  continue
    return
    end
    !--------------------------------------------------------------------
    !------- Gauss Full delimination to Solve the system of eqn's--------
    !--------------------------------------------------------------------
    subroutine gsijsolvep(kaa,k11,k33,aax,aff,mx,my,mz,id)
    use paraandcons
    INTEGER ::  ik(0:2*(nhx+1),0:nz,0:ny),jk(0:2*(nhx+1),0:nz,0:ny)
    INTEGER ::  ki(0:mz),kj(0:mz)
    complex(kind=8) aax(0:mz,0:mz),aff(0:mz)
    real*8 dd,dq
    complex(kind=8) add,as
    integer k11,k33,im,jm,kaa,mx,my,mz,id,i,j,k

    if(kaa.eq.0) then
        ki=1
        kj=1
    end if
    do 10 k=0,mz
        if(kaa.eq.0) then
            dd=0.0d0
            do 15 i=0,mz
                if(ki(i).eq.1) then
                    do 16 j=0,mz
                        dq=dsqrt((real(aax(i,j),kind=8))**2+(dimag(aax(i,j)))**2)
                        if(kj(j).eq.1.and.dq.gt.dd) then
                            dd=dq
                            im=i
                            jm=j
                        end if
16                  continue
                end if
15          continue
            ik(k11,k,k33)=im
            jk(k11,k,k33)=jm
            kj(jm)=0
            ki(im)=0
            if(k.eq.mz) goto 10
            add=1.0d0/aax(im,jm)
            do 25 i=0,mz
                if(ki(i).eq.1) then
                    as=aax(i,jm)*add
                    do 20 j=0,mz
                        if(kj(j).eq.1) aax(i,j)=aax(i,j)-as*aax(im,j)
20                  continue
                    aff(i)=aff(i)-as*aff(im)
                end if
25          continue
        else
            im=ik(k11,k,k33)
            jm=jk(k11,k,k33)
            if(k.eq.mz) goto 10
            add=1.0d0/aax(im,jm)
            do 65 li=k+1,mz
                i=ik(k11,li,k33)
                as=aax(i,jm)*add
                do 60 lj=k+1,mz
                    j=jk(k11,lj,k33)
                    aax(i,j)=aax(i,j)-as*aax(im,j)
60              continue
                aff(i)=aff(i)-as*aff(im)
65          continue
        end if

10  continue

    aax(im,jm)=aff(im)/aax(im,jm)
    do 30 k=mz-1,0,-1
        i=ik(k11,k,k33)
        j=jk(k11,k,k33)
        aax(im,j)=aff(i)
        do 40 lj=k+1,mz
            jj=jk(k11,lj,k33)
            aax(im,j)=aax(im,j)-aax(i,jj)*aax(im,jj)
40      continue
        aax(im,j)=aax(im,j)/aax(i,j)
30  continue
    do 50 j=0,mz
        aff(j)=aax(im,j)
50  continue
    return
    end
    !================add by myself use to solve instant pressure==============================
    !unp1p,vnp1p,wnp1p are instant velocity in physical space without laminar flow
    subroutine instantp(x,y,z,id,pressure,nallgrp,kaa,ierr,unp1p,vnp1p,wnp1p)
    use paraandcons
    include 'mpif.h'
    real*8::convect1(0:nx,0:ny,0:nz),convect2(0:nx,0:ny,0:nz),convect3(0:nx,0:ny,0:nz),array(0:nx,0:ny,0:nz)
    real*8::pressure(0:nx,0:ny,0:nz),unp1p(0:nx,0:ny,0:nz),vnp1p(0:nx,0:ny,0:nz),wnp1p(0:nx,0:ny,0:nz)
    integer::i,j,k,p,id,nallgrp,kaa,ierr
    real*8::x(0:nx),y(0:ny),z(0:nz),derive(0:nz),rhsr(0:nz),rhsi(0:nz),du0dy(0:nz)
    complex(kind=8)::a(0:nz,0:nz),rhs(0:nz)
    real*8::ibeta(0:ny),iapha(0:nhx)
    real*8::pi
    real*8::coriolist(0:nx,0:ny,0:nz) !changeRPCF
    real*8::coriolinm(0:nx,0:ny,0:nz) !changeRPCF
    real*8::coriolisp(0:nx,0:ny,0:nz) !changeRPCF


    pi=4.0d0*atan(1.0d0)
    pressure=real(0,kind=8)
    du0dy=0.5d0 ! laminar flow couette du0/dy=0.5d0


    if(rsp.eq.real(0,kind=8))then
        coriolist=0.0d0 !changeRPCF
        coriolinm=0.0d0 !changeRPCF
        coriolisp=0.0d0 !changeRPCF

    else
        do k=0,nz
            unp1p(:,:,k)=unp1p(:,:,k)+u0(k)
        enddo

        coriolist=rnm*wnp1p-rsp*vnp1p !changeRPCF
        coriolinm=rsp*unp1p-rst*wnp1p !changeRPCF
        coriolisp=rst*vnp1p-rnm*unp1p !changeRPCF

        do k=0,nz
            unp1p(:,:,k)=unp1p(:,:,k)-u0(k)
        enddo
    endif
    !do k=0,nz
    !unp1p(:,:,k)=unp1p(:,:,k)-u0(k)
    !enddo

    !do i=0,nx
    !x(i)=2.0d0*pi*(real(id*(nx+1)+i))/real((nx+1)*nproc)*aphi
    !enddo
    !so delta x=2.0d0*pi*aphi/real(nnx,kind=8)
    !*******************************************************************************************
    !................................non (0,0) mod.............................................
    call getconvection(convect1,unp1p,x,y,z,id,nallgrp,ierr,2.0d0*pi*aphi/real(nnx,kind=8),unp1p,vnp1p,wnp1p)
    array=real(0,kind=8)
    do k=0,nz
        do j=0,ny
            !call deri1d2ordernon(0,nx,unp1p(:,j,k),array(:,j,k),x(:)) !dudx
            call derivexgather(0,nx,unp1p(:,j,k),array(:,j,k),2.0d0*pi*aphi/real(nnx,kind=8),nallgrp,nproc,ierr,id)
        enddo
    enddo
    do k=0,nz
        convect1(:,:,k)=convect1(:,:,k)+u0(k)*array(:,:,k)+du0dy(k)*vnp1p(:,:,k) !uo*dudx+v*du0dy
    enddo

    call getconvection(convect2,vnp1p,x,y,z,id,nallgrp,ierr,2.0d0*pi*aphi/real(nnx,kind=8),unp1p,vnp1p,wnp1p)
    array=real(0,kind=8)
    do k=0,nz
        do j=0,ny
            !call deri1d2ordernon(0,nx,vnp1p(:,j,k),array(:,j,k),x(:)) !dvdx
            call derivexgather(0,nx,vnp1p(:,j,k),array(:,j,k),2.0d0*pi*aphi/real(nnx,kind=8),nallgrp,nproc,ierr,id)
        enddo
    enddo
    do k=0,nz
        convect2(:,:,k)=convect2(:,:,k)+u0(k)*array(:,:,k) !uo*dvdx
    enddo

    call getconvection(convect3,wnp1p,x,y,z,id,nallgrp,ierr,2.0d0*pi*aphi/real(nnx,kind=8),unp1p,vnp1p,wnp1p)
    array=real(0,kind=8)
    do k=0,nz
        do j=0,ny
            !call deri1d2ordernon(0,nx,wnp1p(:,j,k),array(:,j,k),x(:)) !dwdx
            call derivexgather(0,nx,wnp1p(:,j,k),array(:,j,k),2.0d0*pi*aphi/real(nnx,kind=8),nallgrp,nproc,ierr,id)
        enddo
    enddo
    do k=0,nz
        convect3(:,:,k)=convect3(:,:,k)+u0(k)*array(:,:,k) !uo*dwdx
    enddo
    !===================Corioli force==================
    convect1=convect1+coriolist !changeRPCF
    convect2=convect2+coriolinm
    convect3=convect3+coriolisp
    !---------------------------------------------------------
    !   isign=1: K-space --> Phys space; isign=-1: Phys space-
    !   --> K-space.                                         -
    !---------------------------------------------------------
    call transform(-1,convect1,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,convect2,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,convect3,id,nallgrp,nx,ny,nz,nproc)
    call transform(-1,vnp1p,id,nallgrp,nx,ny,nz,nproc)
    !convect1(2*i,j,k)����real part
    !convect1(2*i+1,j,k)����image part



    do j=0,(ny+1)/2
        ibeta(j)=real(j,kind=8)/beta
        if(j.le.(ny-1)/2-1)then
            ibeta(ny-j)=real(-j-1,kind=8)/beta
        endif
    enddo
    do i=0,nhx
        iapha(i)=(real(id,kind=8)*real(nx+1,kind=8)*0.5d0+real(i,kind=8))/aphi
    enddo
    !for every (m,n) mod
    do i=0,nhx
        do j=0,ny
            if((j.ne.(ny+1)/2).and.((iapha(i)**2+ibeta(j)**2).ne.0.0d0))then
                do k=2,nz
                    rhsr(k)=real(0,kind=8)
                    rhsi(k)=real(0,kind=8)
                  do p=0,nz
                    rhsr(k)=rhsr(k)+iapha(i)*convect1(2*i+1,j,p)*t0(k-1,p) &
                    +ibeta(j)*convect3(2*i+1,j,p)*t0(k-1,p)-t1(k-1,p)*convect2(2*i,j,p)

                    rhsi(k)=rhsi(k)-iapha(i)*convect1(2*i,j,p)*t0(k-1,p) &
                    -ibeta(j)*convect3(2*i,j,p)*t0(k-1,p)-t1(k-1,p)*convect2(2*i+1,j,p)

                    enddo
                    rhs(k)=dcmplx(rhsr(k),rhsi(k))
                enddo
                rhsr(0)=real(0,kind=8)
                rhsi(0)=real(0,kind=8)
                rhsr(1)=real(0,kind=8)
                rhsi(1)=real(0,kind=8)
                do p=0,nz
                !rhsr(0)=rhsr(0)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(0,p)+t2(0,p))*vnp1p(2*i,j,p)/re+u0(0)*iapha(i)*t0(0,p)*vnp1p(2*i+1,j,p)!y=1
                !rhsi(0)=rhsi(0)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(0,p)+t2(0,p))*vnp1p(2*i+1,j,p)/re-u0(0)*iapha(i)*t0(0,p)*vnp1p(2*i,j,p)
                !rhsr(1)=rhsr(1)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(nz,p)+t2(nz,p))*vnp1p(2*i,j,p)/re+u0(nz)*iapha(i)*t0(nz,p)*vnp1p(2*i+1,j,p)!y=-1
                !rhsi(1)=rhsi(1)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(nz,p)+t2(nz,p))*vnp1p(2*i+1,j,p)/re-u0(nz)*iapha(i)*t0(nz,p)*vnp1p(2*i,j,p)
                rhsr(0)=rhsr(0)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(0,p) &
                    +t2(0,p))*vnp1p(2*i,j,p)/re-convect2(2*i,j,p)*t0(0,p)!y=1
                rhsi(0)=rhsi(0)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(0,p) &
                    +t2(0,p))*vnp1p(2*i+1,j,p)/re-convect2(2*i+1,j,p)*t0(0,p)
                rhsr(1)=rhsr(1)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(nz,p) &
                    +t2(nz,p))*vnp1p(2*i,j,p)/re-convect2(2*i,j,p)*t0(nz,p)!y=-1
                rhsi(1)=rhsi(1)+(-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(nz,p) &
                    +t2(nz,p))*vnp1p(2*i+1,j,p)/re-convect2(2*i+1,j,p)*t0(nz,p)
                enddo
                rhs(0)=dcmplx(rhsr(0),rhsi(0))
                rhs(1)=dcmplx(rhsr(1),rhsi(1))

                do k=2,nz
                    do p=0,nz
                        rhsr(p)=-1.0d0*(iapha(i)**2+ibeta(j)**2)*t0(k-1,p)+t2(k-1,p)
                        rhsi(p)=real(0,kind=8)
                        a(k,p)=dcmplx(rhsr(p),rhsi(p))
                    enddo
                enddo
                do p=0,nz
                    a(0,p)=dcmplx(t1(0,p),real(0,kind=8))
                    a(1,p)=dcmplx(t1(nz,p),real(0,kind=8))
                enddo


                call gsijsolvep(kaa,i,j,a,rhs,nhx+1,ny,nz,id)


                do k=0,nz
                    pressure(2*i+1,j,k)=dimag(rhs(k))
                    pressure(2*i,j,k)=real(rhs(k),kind=8)
                enddo
            endif !  if((j.ne.(ny+1)/2).and.((iapha(i)**2+ibeta(j)**2).ne.0.0d0))then
        enddo
    enddo
    !*************************************************************************************************
    !...............................................................(0,0) mod....
    if(id==0)then
        !for  (0,0) mod

        pressure(1,0,:)=real(0,kind=8)
        pressure(0,0,:)=real(0,kind=8)

    endif !if(id==0)then
    !............................................................................
    pressure(:,(ny+1)/2,:)=real(0,kind=8)


    call transform(1,pressure,id,nallgrp,nx,ny,nz,nproc)
    call transform(1,vnp1p,id,nallgrp,nx,ny,nz,nproc)

    !================================================================
    return

    endsubroutine


    !.................................................................................
    subroutine getconvection(convect1,a,x,y,z,id,nallgrp,ierr,dx,unp1p,vnp1p,wnp1p)
    use paraandcons
    include 'mpif.h'
    real*8::convect1(0:nx,0:ny,0:nz),a(0:nx,0:ny,0:nz),array(0:nx,0:ny,0:nz)
    real*8::unp1p(0:nx,0:ny,0:nz),vnp1p(0:nx,0:ny,0:nz),wnp1p(0:nx,0:ny,0:nz)
    integer::i,j,k,id,nallgrp,ierr
    real*8::x(0:nx),y(0:ny),z(0:nz),derive(0:nz)
    real*8::dx

    convect1=real(0,kind=8)
    array=real(0,kind=8)
    do k=0,nz
        do j=0,ny
            !call deri1d2ordernon(0,nx,a(:,j,k),array(:,j,k),x(:)) !dadx
            call derivexgather(0,nx,a(:,j,k),array(:,j,k),dx,nallgrp,nproc,ierr,id)
        enddo
    enddo
    convect1=convect1+unp1p*array   !u*dadx

    array=real(0,kind=8)
    do k=0,nz
        do i=0,nx
            call deri1d2ordernon(0,ny,a(i,:,k),array(i,:,k),y(:)) !dadz
        enddo
    enddo
    convect1=convect1+wnp1p*array !w*dadz

    array=real(0,kind=8)                                       !dady
    do j=0,ny
        do i=0,nx
            derive(:)=a(i,j,:)
            call chebyshev(-1,derive,nz)
            call dz(derive,nz)
            call chebyshev(1,derive,nz)
            array(i,j,:)=derive(:)
        enddo
    enddo
    convect1=convect1+vnp1p*array !v*dudy

    return
    endsubroutine
    !...................................................................................
    subroutine deri1d2ordernon(pbeg,pend,f,fd,yj)
    implicit none
    integer::j
    integer::pbeg,pend
    real*8,dimension(pbeg:pend)::f,fd,yj
    real*8 a

    fd=0.0d0

do j=pbeg+2,pend-2
    a=0.0d0
    a=f(j-2)*(yj(j)-yj(j-1))*(yj(j)-yj(j+1))*(yj(j)-yj(j+2))/(yj(j-2)-yj(j-1))/(yj(j-2)-yj(j))/(yj(j-2)-yj(j+1))/(yj(j-2)-yj(j+2))
    a=a+f(j-1)*(yj(j)-yj(j-2))*(yj(j)-yj(j+1))*(yj(j)-yj(j+2))/(yj(j-1)-yj(j-2))/(yj(j-1)-yj(j))/(yj(j-1)-yj(j+1))/(yj(j-1)-yj(j+2))
    a=a+f(j)*((yj(j)-yj(j-1))*(yj(j)-yj(j+1))*(yj(j)-yj(j+2))+(yj(j)-yj(j-2))*(yj(j)-yj(j+1))*(yj(j)-yj(j+2))+(yj(j)-yj(j-2))&
        *(yj(j)-yj(j-1))*(yj(j)-yj(j+2))+(yj(j)-yj(j-2))*(yj(j)-yj(j-1))*(yj(j)-yj(j+1)))/(yj(j)-yj(j-2))&
        /(yj(j)-yj(j-1))/(yj(j)-yj(j+1))/(yj(j)-yj(j+2))
    a=a+f(j+1)*(yj(j)-yj(j-2))*(yj(j)-yj(j-1))*(yj(j)-yj(j+2))/(yj(j+1)-yj(j-2))/(yj(j+1)-yj(j-1))/(yj(j+1)-yj(j))/(yj(j+1)-yj(j+2))
    a=a+f(j+2)*(yj(j)-yj(j-2))*(yj(j)-yj(j-1))*(yj(j)-yj(j+1))/(yj(j+2)-yj(j-2))/(yj(j+2)-yj(j-1))/(yj(j+2)-yj(j))/(yj(j+2)-yj(j+1))
    fd(j)=a
enddo
j=pbeg
fd(j)=(2.0d0*yj(j)-yj(j+1)-yj(j+2))/(yj(j)-yj(j+1))/(yj(j)-yj(j+2))*f(j)+(yj(j)-yj(j+2))/(yj(j+1)-yj(j))&
    /(yj(j+1)-yj(j+2))*f(j+1)+(yj(j)-yj(j+1))/(yj(j+2)-yj(j+1))/(yj(j+2)-yj(j))*f(j+2)
j=pbeg+1
fd(j)=(2.0d0*yj(j)-yj(j+1)-yj(j+2))/(yj(j)-yj(j+1))/(yj(j)-yj(j+2))*f(j)+(yj(j)-yj(j+2))/(yj(j+1)-yj(j))&
    /(yj(j+1)-yj(j+2))*f(j+1)+(yj(j)-yj(j+1))/(yj(j+2)-yj(j+1))/(yj(j+2)-yj(j))*f(j+2)
j=pend-3
fd(pend-1)=(2.0d0*yj(j)-yj(j+1)-yj(j+2))/(yj(j)-yj(j+1))/(yj(j)-yj(j+2))*f(j)+(yj(j)-yj(j+2))/(yj(j+1)-yj(j))&
    /(yj(j+1)-yj(j+2))*f(j+1)+(yj(j)-yj(j+1))/(yj(j+2)-yj(j+1))/(yj(j+2)-yj(j))*f(j+2)
j=pend-2
fd(pend)=(2.0d0*yj(j)-yj(j+1)-yj(j+2))/(yj(j)-yj(j+1))/(yj(j)-yj(j+2))*f(j)+(yj(j)-yj(j+2))/(yj(j+1)-yj(j))&
    /(yj(j+1)-yj(j+2))*f(j+1)+(yj(j)-yj(j+1))/(yj(j+2)-yj(j+1))/(yj(j+2)-yj(j))*f(j+2)
return

    endsubroutine
    !=====================================================================================
    subroutine derivexgather(pbeg,pend,f,fd,dx,nallgrp,nproc,ierr,id)
    include 'mpif.h'
    integer::pbeg,pend,pendall,nallgrp,nproc,ierr,id
    real*8,allocatable,dimension(:)::fall,fdall !,xall
    real*8::f(pbeg:pend),fd(pbeg:pend)
    real*8::dx

    pendall=(pend-pbeg+1)*nproc+pbeg-1

    allocate(fall(pbeg:pendall),fdall(pbeg:pendall)) !,xall(pbeg:pendall))

    call mpi_gather(f,pend-pbeg+1,MPI_DOUBLE_PRECISION,fall,pend-pbeg+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)
    !call mpi_gather(x,pend-pbeg+1,MPI_DOUBLE_PRECISION,xall,pend-pbeg+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)
    call mpi_barrier(nallgrp,ierr)
    if(id==0)then
        !call deri1d2ordernon(pbeg,pendall,fall,fdall,xall)
        call deri1d2order(pbeg,pendall,fall,fdall,dx)
    endif

    call mpi_scatter(fdall,pend-pbeg+1,MPI_DOUBLE_PRECISION,fd,pend-pbeg+1,MPI_DOUBLE_PRECISION,0,nallgrp,ierr)

    deallocate(fall,fdall)
    return
    endsubroutine
    !============================================================================================
    subroutine deri1d2order(pbeg,pend,f,fd,h)
    implicit none

    integer::pbeg,pend
    real*8,dimension(pbeg:pend)::f,fd
    real*8::h

    integer::j

    do j=pbeg+2,pend-2
        fd(j)=2.0d0/3.0d0/h*(f(j+1)-f(j-1))-1.0d0/12.0d0/h*(f(j+2)-f(j-2))
    enddo
    fd(pbeg)=(-11.d0*f(pbeg)+18.d0*f(pbeg+1)-9.d0*f(pbeg+2)+2.d0*f(pbeg+3))/(6.d0*h) !�߽���3�׾���
    fd(pbeg+1)=(-2.d0*f(pbeg)-3.d0*f(pbeg+1)+6.d0*f(pbeg+2)-f(pbeg+3))/(6.d0*h)  !�α߽�3�׾���
    fd(pend-1)=(f(pend-3)-6.d0*f(pend-2)+3.d0*f(pend-1) +2.d0*f(pend))/(6.d0*h) !3��
    fd(pend)=(-2.d0*f(pend-3)+9.d0*f(pend-2)-18.d0*f(pend-1)+11.d0*f(pend))/(6.d0*h) !3��

    return
    endsubroutine

    SUBROUTINE LUIJ(A,B,DIM)
    INTEGER DIM
    COMPLEX(KIND=8) A(DIM,DIM)
    COMPLEX(KIND=8) B(DIM)
    CALL LU(A,DIM)
    CALL solve_lu(A,B,DIM)
    !PRINT*,B
    END

    SUBROUTINE LU(A,DIM)
        INTEGER DIM
        COMPLEX(KIND=8) A(DIM,DIM)
        
        DO I=2,DIM
            DO J=I,DIM
                A(J,I-1) = A(J,I-1)*1.0/A(I-1,I-1)
                m = A(J,I-1)
                DO K=I,DIM
                    A(J,K) = A(J,K)-m*A(I-1,K)
                ENDDO
            ENDDO
        ENDDO
    END

    SUBROUTINE solve_lu(A,B,DIM)
        INTEGER:: DIM
        COMPLEX(KIND=8) A(DIM,DIM),B(DIM)
    !----------LY = B 
        DO I=2,DIM
            DO J=1,I-1
                B(I) = B(I) - A(I,J)*B(J)
            ENDDO
        ENDDO

        B(DIM) = B(DIM)/A(DIM,DIM)
    !---------UX=Y Y stores in B,same as X
        DO I=DIM-1,1,-1
            DO J=I+1,DIM
                B(I) = B(I) - A(I,J)*B(J)
            ENDDO
            B(I) = B(I)/A(I,I)
        ENDDO
    END

    subroutine inverse_matrix(a,n)
    IMPLICIT NONE
    integer :: n
    complex(kind=8) :: a(n,n)
    complex(kind=8),allocatable :: temp(:,:),swap(:)
    integer :: i,j,p,q
    complex(kind=8) :: ratio
    real(kind=8) max_abs
    integer :: max_i
    allocate(temp(n,n),swap(n))
    temp = dcmplx(0.0d0,0.0d0)
    do i=1,n
        temp(i,i) = 1.0d0
    enddo
    do i = 1,n-1
        max_abs = 0.0d0
        do j = i,n ! find max number between a(i,i) and a(i,n)
            if(abs(a(i,j))>max_abs) then
                max_abs = a(i,j)
                max_i = j
            endif
        enddo
        p = max_i
        if(a(p,i).ne.0) then
            swap = a(i,:)
            a(i,:) = a(p,:)
            a(p,:) = swap
            swap = temp(i,:)
            temp(i,:) = temp(p,:)
            temp(p,:) = swap
        else 
            write(*,*) 'error, matrix inversible'
            exit
        endif
        do j = i+1,n
            ratio = -a(j,i)/a(i,i)
            !write(*,*) i,j,ratio
            do p = 1,n
                a(j,p) = a(j,p)+ratio*a(i,p)
                temp(j,p) = temp(j,p)+ratio*temp(i,p)
            enddo
        enddo
    enddo

    do i = n,2,-1
        do j = i-1,1,-1
            ratio = -a(j,i)/a(i,i)
            do p = 1,n
                temp(j,p) = temp(j,p)+ratio*temp(i,p)
                a(j,p) = a(j,p)+ratio*a(i,p)
            enddo
        enddo
    enddo

    do i=1,n
        temp(i,:) = temp(i,:)/a(i,i)
        a(i,:) = a(i,:)/a(i,i)
    enddo
    a = temp
    deallocate(temp,swap)
    end subroutine

    subroutine a_multiply_x(a,x,n)
    IMPLICIT NONE
    integer :: n
    complex(kind=8) :: a(n,n)
    complex(kind=8) :: x(n)
    complex(kind=8),allocatable :: temp(:)
    integer :: i,j
    allocate(temp(n))
    do i=1,n
        temp(i) = 0.0d0
        do j=1,n
            temp(i) = temp(i) + a(i,j)*x(j)
        enddo
    enddo
    x = temp
    deallocate(temp)
    end subroutine

    subroutine matrix_multiply(a,b,c,n)
    IMPLICIT NONE
    integer :: n
    complex(kind=8) :: a(n,n),b(n,n),c(n,n)
    integer :: i,j,k
    do i=1,n
        do j=1,n
            c(i,j) = 0
            do k=1,n
                c(i,j) = c(i,j)+a(i,k)*b(k,j)
            enddo
        enddo
    enddo
    end subroutine
    !..............................................................................
    !-------------------------------------------------------------
    !--------------------END OF THE PROGRAM ----------------------
    !-------------------------------------------------------------
