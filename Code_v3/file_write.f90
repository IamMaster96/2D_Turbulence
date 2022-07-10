subroutine file_write_in
use mpi
use mod_2flu
implicit none
integer :: vortfile
integer(kind=mpi_offset_kind) :: vort_disp
integer :: ispectra

!we write initial files here and  use second subroutine
!to write files in the run

vort_disp=n*nlocal*rank*8

call mpi_file_open(mpi_comm_world,'vortex.in',mpi_mode_wronly+mpi_mode_create,mpi_info_null,vortfile,ierr)
call mpi_file_set_view(vortfile,vort_disp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(vortfile,omega(1:n,1:nlocal),n*nlocal,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(vortfile,ierr)

call mpi_reduce(E_omega,E_omega_out,nshell,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(rank==0)then
open(unit=10,file='spectra.in')
do ispectra=0,nshell
write(10,*)dfloat(ispectra),E_omega_out(ispectra)
end do
end if
close(10)

end subroutine file_write_in


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine file_write
use mpi
use mod_part
use mod_2flu
implicit none
integer :: vortfile,ispectra,xfile,yfile,disp,tfile
integer(kind=mpi_offset_kind) :: vort_disp,posdisp

vort_disp=n*nlocal*rank*8

call mpi_file_open(mpi_comm_world,'vorticity/vortex'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,vortfile,ierr)
call mpi_file_set_view(vortfile,vort_disp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(vortfile,omega(1:n,1:nlocal),n*nlocal,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(vortfile,ierr)

if(particles)then

do istokes=1,nstokes

write(st,'(g8.0)')istokes
x=xpart(istokes,:)
y=ypart(istokes,:)
vx=vxpart(istokes,:)
vy=vypart(istokes,:)
sizex=partsize(istokes)
part_trac=xpart_trac(istokes,:)
disp=0

call mpi_scan(sizex,disp,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
disp=disp-sizex
posdisp=disp*8
call mpi_barrier(mpi_comm_world,ierr)
call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'//non_interacting/xptrack/xptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,xfile,ierr)
call mpi_file_set_view(xfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(xfile,x(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(xfile,ierr)


call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/non_interacting/yptrack/yptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,yfile,ierr)
call mpi_file_set_view(yfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(yfile,y(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(yfile,ierr)
posdisp=disp*4
!print *,part_trac(1:sizex)
call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/non_interacting/trac/trac'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,tfile,ierr)
call mpi_file_set_view(tfile,posdisp,mpi_integer,mpi_integer,'native',mpi_info_null,ierr)
call mpi_file_write_all(tfile,part_trac(1:sizex),sizex,mpi_integer,mpi_status_ignore,ierr)
call mpi_file_close(tfile,ierr)
posdisp=disp*8
call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'//non_interacting/vxptrack/vxptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,xfile,ierr)
call mpi_file_set_view(xfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(xfile,vx(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(xfile,ierr)


call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/non_interacting/vyptrack/vyptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,yfile,ierr)
call mpi_file_set_view(yfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(yfile,vy(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(yfile,ierr)
end do
do istokes=1,nstokes

write(st,'(g8.0)')istokes
x=sxpart(istokes,:)
y=sypart(istokes,:)
vx=svxpart(istokes,:)
vy=svypart(istokes,:)
sizex=spartsize(istokes)
part_trac=sxpart_trac(istokes,:)
disp=0

call mpi_scan(sizex,disp,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
disp=disp-sizex
posdisp=disp*8
call mpi_barrier(mpi_comm_world,ierr)
call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/interacting/xptrack/xptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,xfile,ierr)
call mpi_file_set_view(xfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(xfile,x(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(xfile,ierr)


call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/interacting/yptrack/yptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,yfile,ierr)
call mpi_file_set_view(yfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(yfile,y(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(yfile,ierr)

posdisp=disp*4
call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/interacting/trac/trac'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,tfile,ierr)
call mpi_file_set_view(tfile,posdisp,mpi_integer,mpi_integer,'native',mpi_info_null,ierr)
call mpi_file_write_all(tfile,part_trac(1:sizex),sizex,mpi_integer,mpi_status_ignore,ierr)
!print *,part_trac(1:sizex)
call mpi_file_close(tfile,ierr)
posdisp=disp*8

call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/interacting/vxptrack/vxptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,xfile,ierr)
call mpi_file_set_view(xfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(xfile,vx(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(xfile,ierr)


call mpi_file_open(mpi_comm_world,'tau'//trim(adjustl(st))//'/interacting/vyptrack/vyptrack'//trim(adjustl(fnum))//'.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,yfile,ierr)
call mpi_file_set_view(yfile,posdisp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
call mpi_file_write_all(yfile,vy(1:sizex),sizex,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_close(yfile,ierr)
end do
end if

call mpi_reduce(E_omega,E_omega_out,nshell,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(rank==0)then
open(unit=10,file='spectra/spectra'//trim(adjustl(fnum))//'.out')
do ispectra=0,nshell
write(10,*)dfloat(ispectra),E_omega_out(ispectra)
end do
end if
close(10)
end subroutine file_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine file_write_out
use fftw3
use mpi
use mod_2flu
implicit none
integer :: vortfile
integer(kind=mpi_offset_kind) :: vort_disp
integer :: ispectra

!we write initial files here and  use second subroutine
!to write files in the run

vort_disp=(n+2)*nlocal*rank*8

call fftw_mpi_execute_dft_r2c(fplan,omega,fomega)

call mpi_file_open(mpi_comm_world,'vortex.out',mpi_mode_wronly+mpi_mode_create,mpi_info_null,vortfile,ierr)
call mpi_file_set_view(vortfile,vort_disp,mpi_double_complex,mpi_double_complex,'native',mpi_info_null,ierr)
call mpi_file_write_all(vortfile,fomega(n/2+1,1:nlocal),(n/2+1)*nlocal,mpi_double_complex,mpi_status_ignore,ierr)
call mpi_file_close(vortfile,ierr)

call getpsivel
call energy_2d
call mpi_reduce(E_omega,E_omega_out,nshell,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(rank==0)then
open(unit=10,file='spectra.out')
do ispectra=0,nshell
write(10,*)dfloat(ispectra),E_omega_out(ispectra)
end do
end if
close(10)

end subroutine file_write_out
