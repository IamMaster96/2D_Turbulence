module mod_2flu
use,intrinsic :: iso_c_binding
implicit none
double precision :: scale
!mpi_variables
integer ::  ierr,nprocs,rank
!velocities
real(C_DOUBLE),pointer ::omega(:,:),ukx(:,:),uky(:,:),force(:,:)
complex(C_DOUBLE_COMPLEX),pointer :: fomega(:,:),fukx(:,:),fuky(:,:),psi(:,:),kforce(:,:),jac_old(:,:)
double precision,allocatable,dimension(:) :: E_omega(:),time_incr(:),E_omega_out
integer,allocatable,dimension(:) :: den_state
!array allocation variables(mpi_fftw)
type(C_PTR) :: rdata,cdata
integer(C_INTPTR_T) :: alloc_local,nlocal,local_offset
type(C_PTR) :: fplan,rplan
!fluid parameters
integer(C_INTPTR_T) :: n,nn,maxiter,navg,n1hf,n1h,t,kf,nalias,kasqr,ksqr_max,kdrag,kini,nshell
double precision :: l,vis,vis2,delta,pi,length,factor,mu,mu2
double precision :: famp,dx
character(100) :: fnum

end module
