subroutine vel_comm

use fftw3
use mpi
use mod_2flu
use mod_part

implicit none

call mpi_isend(ukx(1:n,1),n,mpi_double_precision,rank_l,1,mpi_comm_world,ireq1,ierr)
call mpi_isend(ukx(1:n,2),n,mpi_double_precision,rank_l,2,mpi_comm_world,ireq2,ierr)
call mpi_isend(ukx(1:n,nlocal),n,mpi_double_precision,rank_u,3,mpi_comm_world,ireq3,ierr)

call mpi_isend(uky(1:n,1),n,mpi_double_precision,rank_l,4,mpi_comm_world,ireq4,ierr)
call mpi_isend(uky(1:n,2),n,mpi_double_precision,rank_l,5,mpi_comm_world,ireq5,ierr)
call mpi_isend(uky(1:n,nlocal),n,mpi_double_precision,rank_u,6,mpi_comm_world,ireq6,ierr)


call mpi_recv(velx_recv_u_0(1),n,mpi_double_precision,rank_u,1,mpi_comm_world,istatus1,ierr)
call mpi_recv(velx_recv_u_1(1),n,mpi_double_precision,rank_u,2,mpi_comm_world,istatus2,ierr)
call mpi_recv(velx_recv_l(1),n,mpi_double_precision,rank_l,3,mpi_comm_world,istatus3,ierr)

call mpi_recv(vely_recv_u_0(1),n,mpi_double_precision,rank_u,4,mpi_comm_world,istatus4,ierr)
call mpi_recv(vely_recv_u_1(1),n,mpi_double_precision,rank_u,5,mpi_comm_world,istatus5,ierr)
call mpi_recv(vely_recv_l(1),n,mpi_double_precision,rank_l,6,mpi_comm_world,istatus6,ierr)

call mpi_wait(ireq1,istatus1,ierr)
call mpi_wait(ireq2,istatus2,ierr)
call mpi_wait(ireq3,istatus3,ierr)
call mpi_wait(ireq4,istatus4,ierr)
call mpi_wait(ireq5,istatus5,ierr)
call mpi_wait(ireq6,istatus6,ierr)

end subroutine vel_comm
