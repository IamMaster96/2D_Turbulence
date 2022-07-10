subroutine gen_force
use mod_2flu
use fftw3
implicit none
integer :: i,j
double precision :: x

do i=1,nlocal
        do j=1,n
        x=j*dx
        force(j,i)=-famp*kf*dcos(kf*x)
        end do
end do

call fftw_mpi_execute_dft_r2c(fplan,force,kforce)
end subroutine gen_force
