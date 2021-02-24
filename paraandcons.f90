module paraandcons 
    implicit none !!!!!!!!!!!!!!!! channel.inc
    !grid points nnx*nny*(nnz+1)
    integer, parameter :: nproc=64      ! the number of the processors
    integer, parameter :: nnx=128       !dimensions; NOTE:nny divided by 2 must be multiflies of nproc!!!
    integer, parameter :: nny=128
    integer, parameter :: nnz=128
    integer, parameter :: nx=nnx/nproc-1
    integer, parameter :: ny=nny-1
    integer, parameter :: nz=nnz          
    integer, parameter :: nhx=(nx-1)/2    ! the mode that will be calculated
    integer, parameter :: ntime=2         ! time for input at the initial
    integer, parameter :: dstep=-1        ! steps which will disturb the bourdary conditions
        
    integer, parameter :: istepbegin=0   ! the loop will be stopped at 
    integer, parameter :: nstep=5000000    ! the loop will be stopped at 
    integer, parameter :: interupt=1000   ! output statistic meanu and so on at this times

    ! used in subroutine statistics
    integer, parameter :: outputstep=1000   ! Fourier inverse transform(spectra->physical) and output velocity of physical space at this times
        
    real*8, parameter :: aphi=2.0d0,beta=1.0d0    ! the wave number of x and y direction
    real*8, parameter :: dt=2.0e-4           ! time step
    real*8, parameter :: re=180.0d0             ! Reynolds number, which define viscosity
    real*8, parameter :: re_p=1.0d0             ! 
    real*8, parameter :: rst=0.0d0, rnm=0.0d0, rsp=0.0d0     ! rotate number(=2*ration vector/(U_w/h))in three (streamwise, spanwise and spanwise) directions
    
    real*8,save,allocatable,dimension(:,:) :: t0,t2,t4,t1 ! chebyshev poly. and its derivatives
    real*8,save,allocatable,dimension(:) :: u0,du0,ddu0   ! base flow velocity
    
    real*8,save,allocatable,dimension(:,:,:) :: rhsvn,rhsomegayn !rhs v n, rhs omega
    real*8,save,allocatable,dimension(:,:,:) :: unp1,vnp1
    real*8,save,allocatable,dimension(:,:) :: clamb1nm1,clamb3nm1,clamb1n,clamb3n ! coriolis lamb vector
    real*8,save,allocatable,dimension(:,:,:) :: wnp1
    real*8,save,allocatable,dimension(:,:,:) :: omegaxnp1
    real*8,save,allocatable,dimension(:,:,:) :: omegaynp1
    real*8,save,allocatable,dimension(:,:,:) :: omegaznp1
    real*8,save,allocatable,dimension(:,:) :: au0,w0
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer*8,save:: plan1,plan2,plan3,plan4,plan5,plan6,plan7,plan8,plan9,plan10,plan11,plan12,plan13
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    !real*8,save,dimension(:,:) :: t0(0:nz,0:nz),t2(0:nz,0:nz),t4(0:nz,0:nz)
    !real*8,save,dimension(:) :: u0(0:nz)
    !common/tt/ t0,t2,t4,u0
    
    !real*8,save,dimension(:,:,:) :: rhsvn(0:nx,0:ny,0:nz),rhsomegayn(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:,:) :: unp1(0:nx,0:ny,0:nz),vnp1(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:) :: clamb1nm1(2,0:nz),clamb3nm1(2,0:nz),clamb1n(2,0:nz),clamb3n(2,0:nz)
    !real*8,save,dimension(:,:,:) :: wnp1(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:,:) :: omegaxnp1(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:,:) :: omegaynp1(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:,:) :: omegaznp1(0:nx,0:ny,0:nz)
    !real*8,save,dimension(:,:) :: au0(2,0:nz),w0(2,0:nz)
!---------------------------------------
!      dimension dpx(2,0:nz)
!---------------------------------------
 !     common/rhs1/ rhsvn,rhsomegayn
 !     common/data1/ unp1,clamb1n
 !     common/data2/ vnp1,clamb3n
 !     common/clamb/ clamb1nm1,clamb3nm1,au0,w0
  !    common/data3/ wnp1
 !     common/cpomegax/ omegaxnp1
  !    common/cpomegay/ omegaynp1
  !    common/cpomegaz/ omegaznp1
 end module paraandcons
