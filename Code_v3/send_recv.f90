subroutine send_recv
use mpi
use mod_2flu
use mod_part

implicit none

integer :: i,j,k
integer :: recv_l,recv_u
integer :: binx,biny

bin=0

send_ly=0.0d0
send_lx=0.0d0
send_uy=0.0d0
send_ux=0.0d0

!recv_ly=0.0d0
!recv_lx=0.0d0
!recv_uy=0.0d0
!recv_ux=0.0d0

j=0
k=0

do i=1,sizex
                xcord=x(i); ycord=y(i)
                binx=int(xcord*nbin)+1-int(nbox/nprocs)*rank
                biny=int(ycord*nbin)+1
                bin(0,binx,biny)=bin(0,binx,biny)+1
                bin(bin(0,binx,biny),binx,biny)=i
                        if(binx.le.1)then
                                j=j+1
                                send_ly(j)=ycord
                                send_lx(j)=xcord
                                send_ltrac(j)=part_trac(i)
                        else if(binx.ge.local_x_box)then
                                k=k+1
                                send_utrac(k)=part_trac(i)
                                send_uy(k)=ycord
                                send_ux(k)=xcord
                        end if
end do

call mpi_isend(j,1,mpi_int,rank_l,1,mpi_comm_world,ireq1,ierr)
call mpi_isend(k,1,mpi_int,rank_u,2,mpi_comm_world,ireq2,ierr)
call mpi_recv(recv_u,1,mpi_int,rank_u,1,mpi_comm_world,istatus1,ierr)
call mpi_recv(recv_l,1,mpi_int,rank_l,2,mpi_comm_world,istatus2,ierr)
call mpi_wait(ireq1,istatus1,ierr)
call mpi_wait(ireq2,istatus2,ierr)

call mpi_isend(send_ly(1),j,mpi_double_precision,rank_l,3,mpi_comm_world,ireqly,ierr)
call mpi_isend(send_lx(1),j,mpi_double_precision,rank_l,4,mpi_comm_world,ireqlx,ierr)
call mpi_isend(send_ltrac(1),j,mpi_integer,rank_l,9,mpi_comm_world,ireqltrac,ierr)

call mpi_isend(send_utrac(1),k,mpi_integer,rank_u,10,mpi_comm_world,irequtrac,ierr)
call mpi_isend(send_uy(1),k,mpi_double_precision,rank_u,6,mpi_comm_world,irequy,ierr)
call mpi_isend(send_ux(1),k,mpi_double_precision,rank_u,7,mpi_comm_world,irequx,ierr)

call mpi_recv(y(sizex+1),recv_u,mpi_double_precision,rank_u,3,mpi_comm_world,istatus3,ierr)
call mpi_recv(x(sizex+1),recv_u,mpi_double_precision,rank_u,4,mpi_comm_world,istatus4,ierr)
call mpi_recv(part_trac(sizex+1),recv_u,mpi_integer,rank_u,9,mpi_comm_world,istatus7,ierr)

call mpi_recv(y(sizex+1+recv_u),recv_l,mpi_double_precision,rank_l,6,mpi_comm_world,istatus5,ierr)
call mpi_recv(x(sizex+1+recv_u),recv_l,mpi_double_precision,rank_l,7,mpi_comm_world,istatus6,ierr)
call mpi_recv(part_trac(sizex+1+recv_u),recv_l,mpi_integer,rank_l,10,mpi_comm_world,istatus8,ierr)

call mpi_wait(ireqly,istatus3,ierr)
call mpi_wait(ireqlx,istatus4,ierr)
call mpi_wait(irequy,istatus5,ierr)
call mpi_wait(irequx,istatus6,ierr)
call mpi_wait(ireqltrac,istatus7,ierr)
call mpi_wait(irequtrac,istatus8,ierr)

do  i=sizex+1+recv_u,sizex+recv_u+recv_l       
        biny=int(y(i)*nbin)+1
        bin(0,0,biny)=bin(0,0,biny)+1
        bin(bin(0,0,biny),0,biny)=i
end do

do i=sizex+1,sizex+recv_u

        biny=int(y(i)*nbin)+1
        bin(0,local_x_box+1,biny)=bin(0,local_x_box+1,biny)+1
        bin(bin(0,local_x_box+1,biny),local_x_box+1,biny)=i
end do

end subroutine send_recv
