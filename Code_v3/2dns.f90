program main
use mpi
use fftw3
use mod_2flu
use mod_part

implicit none
integer :: i,j,kx,ky,ksqr,i1
integer ::nrun, iiter,iiner
double precision :: t1,t2
integer :: vortfile,counter,part_t
integer(kind=mpi_offset_kind) :: vort_disp
character*100 :: fname1
!Read configuration file
!n=grid size,vis=viscosity
!we will be using exponential integrator with rk2 to solve in fourier space
!we start with a random distribution of velocities in fourier space



!initialize mpi

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

open(unit=10,file='para.in',status='old')
read(10,*)
read(10,*)n,vis,vis2,delta,maxiter,navg,nrun,nbox,bin_cap
read(10,*)
read(10,*)kini,kf,famp,mu,kdrag,mu2
read(10,*)
read(10,*)particles,nparticle,nstokes,meand,nrunpart,part_t
close(10)

!calculate all relevant variables
if(rank.eq.0)then
call system('mkdir vorticity spectra')
end if

if(rank.eq.0)then
if(particles) then
  	do i1 = 1,nstokes
   	 write(fname1,'(g8.0)')  i1
	 call system('mkdir tau'//trim(adjustl(fname1)))
         call system('mkdir tau'//trim(adjustl(fname1))//'/interacting')
         call system('mkdir tau'//trim(adjustl(fname1))//'/interacting/xptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/interacting/yptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/interacting/vxptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/interacting/vyptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/interacting/trac')
         call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting')
         call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting/xptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting/yptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting/vxptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting/vyptrack')
        call system('mkdir tau'//trim(adjustl(fname1))//'/non_interacting/trac')
        enddo
endif
end if
pi=4.0d0*atan(1.0d0)
l=2.0d0*pi
dx=l/dfloat(n)
nn=n*n
n1h=n/2
n1hf=n1h+1
nshell=int(1.414d0*n/2.0d0)+1
nalias=int(n/3)
kasqr=2*nalias**2
ksqr_max=2*(n/2)*(n/2)
scale=1.0d0/dfloat(nn)
factor=2.0d0*pi/l
meand2=2.0d0*meand
!start time

t1=mpi_wtime()

!array_allocation(full array avilible on each processor not used in fftw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if(particles)then
        nbin=(dfloat(nbox)/l)
        local_nparticle=int(nparticle/nprocs)
        local_x_box=int(nbox/nprocs)

        allocate(x(1:nparticle),y(1:nparticle),x0(nstokes,1:nparticle),y0(nstokes,1:nparticle),part_trac(1:nparticle))
        allocate(vx(1:nparticle),vy(1:nparticle),vx0(nstokes,1:nparticle),vy0(nstokes,1:nparticle),ux(1:nparticle),uy(1:nparticle))
        allocate(xtemp(1:nparticle),ytemp(1:nparticle))
        allocate(vxtemp(1:nparticle),vytemp(1:nparticle),part_temp(1:nparticle))
        allocate(fx(1:nparticle),fy(1:nparticle))
        allocate(ftemp(1:2))
        allocate(velx_recv_l(1:n),velx_recv_u_0(1:n),vely_recv_l(1:n),vely_recv_u_0(1:n))
        allocate(velx_recv_u_1(1:n),vely_recv_u_1(1:n))
        allocate(send_ly(1:nparticle/2),send_lx(1:nparticle/2))
        allocate(send_uy(1:nparticle/2),send_ux(1:nparticle/2))
        allocate(send_lvy(1:nparticle/2),send_lvx(1:nparticle/2),send_utrac(1:nparticle/2))
        allocate(send_uvy(1:nparticle/2),send_uvx(1:nparticle/2),send_ltrac(1:nparticle/2))
        allocate(vxpart(nstokes,1:nparticle),vypart(nstokes,1:nparticle),xpart_trac(nstokes,1:nparticle))
        allocate(xpart(nstokes,1:nparticle),ypart(nstokes,1:nparticle))
        allocate(sx0(nstokes,1:nparticle),sy0(nstokes,1:nparticle))
        allocate(svx0(nstokes,1:nparticle),svy0(nstokes,1:nparticle))
        allocate(bin(0:bin_cap,0:local_x_box+1,1:nbox))
        allocate(sxpart(nstokes,nparticle),sypart(nstokes,nparticle),svypart(nstokes,nparticle),svxpart(nstokes,nparticle),sxpart_trac(nstokes,1:nparticle))
        allocate(taup(1:nstokes),partsize(nstokes),spartsize(nstokes),trmsint(1:nstokes),trmsnint(1:nstokes),trmsint_out(1:nstokes),trmsnint_out(1:nstokes))
        open(unit=11,file='taup.in',status='old')
        do istokes=1,nstokes
                read(11,*)taup(istokes)
        end do
        close(11)
end if
!!!!!!!!!!!!!particle code!!!!!!!!!!!!
allocate(den_state(0:nshell))
allocate(E_omega(0:nshell),E_omega_out(0:nshell))
allocate(time_incr(0:ksqr_max))

!fftw_mpi array allocation

alloc_local=fftw_mpi_local_size_2d(n,n/2+1,mpi_comm_world,nlocal,local_offset)

rdata=fftw_alloc_real(2*alloc_local)
call c_f_pointer(rdata,omega,[n+2,nlocal])

rdata=fftw_alloc_real(2*alloc_local)
call c_f_pointer(rdata,ukx,[n+2,nlocal])

rdata=fftw_alloc_real(2*alloc_local)
call c_f_pointer(rdata,uky,[n+2,nlocal])

rdata=fftw_alloc_real(2*alloc_local)
call c_f_pointer(rdata,force,[n+2,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,jac_old,[n/2+1,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,kforce,[n/2+1,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,fomega,[n/2+1,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,fukx,[n/2+1,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,fuky,[n/2+1,nlocal])

cdata=fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata,psi,[n/2+1,nlocal])

!fftww_plan creation

fplan= fftw_mpi_plan_dft_r2c_2d(n,n,omega,fomega,mpi_comm_world,0)
rplan= fftw_mpi_plan_dft_c2r_2d(n,n,fomega,omega,mpi_comm_world,0)

!initialize arrays

psi=0.d0
omega=0.0d0
fomega=0.0d0
den_state=0
E_omega=0.0d0
E_omega_out=0.0d0
cord_l=(l/nprocs)*rank
cord_u=(l/nprocs)*(rank+1)

rank_l=rank-1
rank_u=rank+1

if(rank_l.eq.-1)then
        rank_l=nprocs-1
end if

if(rank_u.eq.nprocs)then
        rank_u=0
end if


!create global arrays and initialize the vortiity field and force field constant
!in time

if(particles)then
        call part_init
        do istokes=1,nstokes
                xpart(istokes,:)=x
                ypart(istokes,:)=y
                vxpart(istokes,:)=vx
                vypart(istokes,:)=vy
                partsize(istokes)=sizex
                xpart_trac(istokes,:)=part_trac
                        sxpart(istokes,:)=x
                        sypart(istokes,:)=y
                        svxpart(istokes,:)=vx
                        svypart(istokes,:)=vy
                        spartsize(istokes)=sizex
                        sxpart_trac(istokes,:)=part_trac
        end do
end if
call global2d
call gen_force
if(nrun.eq.1)then
        call iniconf2d
else
        vort_disp=n*nlocal*rank*8
        call mpi_file_open(mpi_comm_world,'vortex.in',mpi_mode_rdonly,mpi_info_null,vortfile,ierr)
        call mpi_file_set_view(vortfile,vort_disp,mpi_double_precision,mpi_double_precision,'native',mpi_info_null,ierr)
        call mpi_file_read_all(vortfile,omega(1:n,1:nlocal),n*nlocal,mpi_double_precision,mpi_status_ignore,ierr)
        close(vortfile)
        call fftw_mpi_execute_dft_r2c(fplan,omega,fomega)

do i=1,nlocal
        kx=i+local_offset
        kx=(kx-1)-n*(kx/(n1hf+1))
        do j=1,n1hf
        ky=j-1
        ksqr=kx**2+ky**2
        if(ksqr > kasqr)then
        fomega(j,i)=cmplx(0.0d0,0.0d0)
        end if
end do
end do
fomega=fomega*scale
end if

!write initial vorticity and spectra file

call fftw_mpi_execute_dft_c2r(rplan,fomega,omega)
call file_write_in
if(rank==0)then
if(nrunpart.eq.1)then
open(unit=200,file='energy.out')
else
open(unit=200,file='energy.out',position='append')
end if
end if

if(rank==0)then
        do istokes=1,nstokes
                write(st,'(g8.0)')istokes
                open(unit=istokes*312,file='tau'//trim(adjustl(st))//'/rms_dis.out')
        end do
end if

counter=1
do iiter=1,maxiter/navg
        do iiner=1,navg
                counter=counter+1+(part_t*navg)
                call time_march
                call getpsivel
                call energy_2d
                call mpi_reduce(E_omega,E_omega_out,nshell,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
                call mpi_reduce(trmsint,trmsint_out,nstokes,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
                call mpi_reduce(trmsnint,trmsnint_out,nstokes,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
                 if(rank==0)then
                        write(200,*)counter*delta,0.5d0*sum(E_omega_out(0:nshell))
                        do istokes=1,nstokes
                                write(istokes*312,*)trmsint_out(istokes)/nparticle,trmsnint_out(istokes)/nparticle
                        end do
                end if
        end do
                if(rank==0)then
                        print *, navg*delta*iiter
                end if
        write(fnum,'(g8.0)')iiter+part_t
!write file after every navg,total number of files=maxiter/navg
                call getpsivel
                call energy_2d
                call fftw_mpi_execute_dft_c2r(rplan,fomega,omega)
                call file_write
                
                !do istokes=1,nstokes
               
                !end do
                !if(rank==0)then
                        
                !end if
                
end do
if(rank==0)then
                        do istokes=1,nstokes
                        close(istokes*1000)
                        end do
                end if
call file_write_out

!writting output files

t2=mpi_wtime()

print *,"time taken=",t2-t1

if(particles)then
        deallocate(xtemp,ytemp)
        deallocate(vxtemp,vytemp)
        
        deallocate(fx,fy)
        deallocate(ftemp)
        deallocate(velx_recv_l,velx_recv_u_0,vely_recv_l,vely_recv_u_0)
        deallocate(velx_recv_u_1,vely_recv_u_1)
        deallocate(send_ly,send_lx)
        deallocate(send_uy,send_ux)
        deallocate(send_lvy,send_lvx)
        deallocate(send_uvy,send_uvx)
        deallocate(bin)
end if
call fftw_destroy_plan(fplan)
call fftw_destroy_plan(rplan)
call fftw_free(rdata)
call fftw_free(cdata)

call mpi_finalize(ierr)

end program main
