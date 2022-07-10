subroutine part_init
use mpi
use mod_part
use mod_2flu

implicit none

double precision,allocatable,dimension(:) :: yini,xini,vxini,vyini
double precision,allocatable,dimension(:,:) :: send_bufx,send_bufy,send_bufvx,send_bufvy
integer,allocatable,dimension(:,:) :: send_buftrac
integer :: i,j,k,send_rank,recv_count
double precision :: rankbin,eps,fac2,efmeand,efmeand2,npartgrid
integer :: istatus(mpi_status_size)
integer :: ipartinit,jpartinit,countpart,sqrtnpart
integer,allocatable,dimension(:) :: sed

rankbin=l/nprocs

if(rank.eq.0)then
k=10
allocate(yini(1:nparticle),xini(1:nparticle))
allocate(vyini(1:nparticle),vxini(1:nparticle))
allocate(sed(1:k))

eps=0.005
sqrtnpart=int(sqrt(1.0d0*nparticle))
npartgrid=l/(1.0d0*sqrtnpart)
if(nrunpart.eq.1)then
        sed=0
        countpart=0
        do ipartinit=1,sqrtnpart
                do jpartinit=1,sqrtnpart
                        countpart=countpart+1
                        xini(countpart)=eps+npartgrid*(jpartinit-1)
                        yini(countpart)=eps+npartgrid*(ipartinit-1)
                        print *, xini(countpart),yini(countpart)
                   end do
        end do
        print *,countpart

        call random_number(vxini(1:nparticle))
        call random_number(vyini(1:nparticle))          
else
open(unit=1000,file='xptrack.in',form='unformatted',access='stream')
open(unit=1001,file='yptrack.in',form='unformatted',access='stream')
open(unit=1002,file='vxptrack.in',form='unformatted',access='stream')
open(unit=1003,file='vyptrack.in',form='unformatted',access='stream')
read(1000)xini(1:nparticle)
read(1001)yini(1:nparticle)
read(1002)vxini(1:nparticle)
read(1003)vyini(1:nparticle)
close(1000)
close(1001)
close(1002)
close(1003)
end if
!call random_number(xini(1:nparticle))
!call random_number(yini(1:nparticle))
!xini=l*xini;
!yini=l*yini
!vxini=0.0d0
!vyini=0.0d0
allocate(send_bufx(0:2*local_nparticle,1:nprocs-1),send_bufy(0:2*local_nparticle,1:nprocs-1))
allocate(send_bufvx(0:2*local_nparticle,1:nprocs-1),send_bufvy(0:2*local_nparticle,1:nprocs-1),send_buftrac(0:2*local_nparticle,1:nprocs-1))

send_bufx=0.0d0
send_bufy=0.0d0
send_bufvx=0.0d0
send_bufvy=0.0d0

end if

if(rank.eq.0)then
                recv_count=0

                do i=1,nparticle
                send_rank=int(xini(i)/rankbin)
                if(send_rank.eq.0)then
                recv_count=recv_count+1
                y(recv_count)=yini(i)
                x(recv_count)=xini(i)
                vy(recv_count)=vyini(i)
                part_trac(recv_count)=i
                vx(recv_count)=vxini(i)
                      else
                send_bufy(0,send_rank)=send_bufy(0,send_rank)+1
                send_bufy(int(send_bufy(0,send_rank)),send_rank)=yini(i)
                send_bufx(int(send_bufy(0,send_rank)),send_rank)=xini(i)
                send_bufvy(int(send_bufy(0,send_rank)),send_rank)=vyini(i)
                send_bufvx(int(send_bufy(0,send_rank)),send_rank)=vxini(i)
                send_buftrac(int(send_bufy(0,send_rank)),send_rank)=i
                end if

                end do

                        do i=1,nprocs-1
                                call mpi_send(int(send_bufy(0,i)),1,mpi_int,i,1,mpi_comm_world,ierr)
                                call mpi_send(send_bufx(1,i),int(send_bufy(0,i)),mpi_double_precision,i,2,mpi_comm_world,ierr)
                                call mpi_send(send_bufy(1,i),int(send_bufy(0,i)),mpi_double_precision,i,3,mpi_comm_world,ierr)
                                call mpi_send(send_bufvx(1,i),int(send_bufy(0,i)),mpi_double_precision,i,5,mpi_comm_world,ierr)
                                call mpi_send(send_bufvy(1,i),int(send_bufy(0,i)),mpi_double_precision,i,6,mpi_comm_world,ierr)
                                call mpi_send(send_buftrac(1,i),int(send_bufy(0,i)),mpi_integer,i,7,mpi_comm_world,ierr)
                        end do
else
                call mpi_recv(recv_count,1,mpi_int,0,1,mpi_comm_world,istatus,ierr)
                call mpi_recv(x(1),recv_count,mpi_double_precision,0,2,mpi_comm_world,istatus,ierr)
                call mpi_recv(y(1),recv_count,mpi_double_precision,0,3,mpi_comm_world,istatus,ierr)
                call mpi_recv(vx(1),recv_count,mpi_double_precision,0,5,mpi_comm_world,istatus,ierr)
                call mpi_recv(vy(1),recv_count,mpi_double_precision,0,6,mpi_comm_world,istatus,ierr)
                call mpi_recv(part_trac(1),recv_count,mpi_integer,0,7,mpi_comm_world,istatus,ierr)
end if
sizex=recv_count
call mpi_barrier(mpi_comm_world,ierr)
print *,"initialized"
!if(rank==0)then
!deallocate(xini,yini,vxini,vyini,send_bufvx,send_bufvy,send_bufx,send_bufy,send_buftrac)
!end if
end subroutine part_init
