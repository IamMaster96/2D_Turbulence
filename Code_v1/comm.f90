subroutine send_recv_x
use mpi
use mod_part
use mod_2flu
implicit none

integer :: i,j,k,newsize
integer :: recv_l,recv_u

vxtemp=0.0d0
vytemp=0.0d0

xtemp=0.0d0
ytemp=0.0d0

send_ly=0.0d0
send_lx=0.0d0

send_ux=0.0d0
send_uy=0.0d0

j=0
k=0

do i=1,sizex
        if(x(i).lt.cord_l)then
                j=j+1
                send_ly(j)=y(i)
                send_lx(j)=x(i)
                send_ltrac(j)=part_trac(i)
                send_lvy(j)=vy(i)
                send_lvx(j)=vx(i)
        end if
        if(x(i).gt.cord_u)then
                k=k+1
                send_uy(k)=y(i)
                send_ux(k)=x(i)
                send_utrac(k)=part_trac(i)
                send_uvy(k)=vy(i)
                send_uvx(k)=vx(i)
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
call mpi_isend(send_uy(1),k,mpi_double_precision,rank_u,6,mpi_comm_world,irequy,ierr)
call mpi_isend(send_ux(1),k,mpi_double_precision,rank_u,7,mpi_comm_world,irequx,ierr)

call mpi_isend(send_utrac(1),k,mpi_integer,rank_u,8,mpi_comm_world,irequtrac,ierr)
call mpi_isend(send_ltrac(1),j,mpi_integer,rank_l,14,mpi_comm_world,ireqltrac,ierr)

call mpi_isend(send_lvy(1),j,mpi_double_precision,rank_l,9,mpi_comm_world,ireqvly,ierr)
call mpi_isend(send_lvx(1),j,mpi_double_precision,rank_l,10,mpi_comm_world,ireqvlx,ierr)
call mpi_isend(send_uvy(1),k,mpi_double_precision,rank_u,12,mpi_comm_world,ireqvuy,ierr)
call mpi_isend(send_uvx(1),k,mpi_double_precision,rank_u,13,mpi_comm_world,ireqvux,ierr)

call mpi_recv(y(sizex+1),recv_u,mpi_double_precision,rank_u,3,mpi_comm_world,istatus3,ierr)
call mpi_recv(x(sizex+1),recv_u,mpi_double_precision,rank_u,4,mpi_comm_world,istatus4,ierr)
call mpi_recv(y(sizex+recv_u+1),recv_l,mpi_double_precision,rank_l,6,mpi_comm_world,istatus5,ierr)
call mpi_recv(x(sizex+recv_u+1),recv_l,mpi_double_precision,rank_l,7,mpi_comm_world,istatus6,ierr)

call mpi_recv(part_trac(sizex+1),recv_u,mpi_integer,rank_u,14,mpi_comm_world,istatus11,ierr)
call mpi_recv(part_trac(sizex+recv_u+1),recv_l,mpi_integer,rank_l,8,mpi_comm_world,istatus12,ierr)

call mpi_recv(vy(sizex+1),recv_u,mpi_double_precision,rank_u,9,mpi_comm_world,istatus7,ierr)
call mpi_recv(vx(sizex+1),recv_u,mpi_double_precision,rank_u,10,mpi_comm_world,istatus8,ierr)
call mpi_recv(vy(sizex+recv_u+1),recv_l,mpi_double_precision,rank_l,12,mpi_comm_world,istatus9,ierr)
call mpi_recv(vx(sizex+recv_u+1),recv_l,mpi_double_precision,rank_l,13,mpi_comm_world,istatus10,ierr)

call mpi_wait(ireqly,istatus3,ierr)
call mpi_wait(ireqlx,istatus4,ierr)
call mpi_wait(irequy,istatus5,ierr)
call mpi_wait(irequx,istatus6,ierr)

call mpi_wait(irequtrac,istatus12,ierr)
call mpi_wait(ireqltrac,istatus11,ierr)

call mpi_wait(ireqvly,istatus7,ierr)
call mpi_wait(ireqvlx,istatus8,ierr)
call mpi_wait(ireqvuy,istatus9,ierr)
call mpi_wait(ireqvux,istatus10,ierr)

x=modulo(x,l)
y=modulo(y,l)

newsize=0

do i=1,sizex+recv_l+recv_u
        if(x(i).gt.cord_l.and.x(i).lt.cord_u)then
                newsize=newsize+1
                ytemp(newsize)=y(i)
                xtemp(newsize)=x(i)
                part_temp(newsize)=part_trac(i)
                vxtemp(newsize)=vx(i)
                vytemp(newsize)=vy(i)
        end if
end do

y=ytemp
x=xtemp
vx=vxtemp
vy=vytemp
part_trac=part_temp
sizex=sizex+recv_l+recv_u-j-k
end subroutine send_recv_x
