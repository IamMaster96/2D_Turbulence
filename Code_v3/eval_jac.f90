subroutine jacb(flag)
use mpi
use mod_2flu
use fftw3
use mod_part

implicit none
integer :: i,j,ksqr,ireal,iimag,il,flag
double precision :: k,kx,ky

call getpsivel

!get back to reality

call fftw_mpi_execute_dft_c2r(rplan,fomega,omega)
call fftw_mpi_execute_dft_c2r(rplan,fukx,ukx)
call fftw_mpi_execute_dft_c2r(rplan,fuky,uky)

ukx=ukx*scale
uky=uky*scale
omega=omega*scale

if(particles)then
        if(flag==1)then
                call vel_comm
                do istokes=1,nstokes
                        x=sxpart(istokes,:)
                        y=sypart(istokes,:)
                        vx=svxpart(istokes,:)
                        vy=svypart(istokes,:)
                        sizex=spartsize(istokes)
                        tau=taup(istokes)
                        part_trac=sxpart_trac(istokes,:)
                         

                               call linear_interp
                                call send_recv
       call soft_pt
                                        do i=1,sizex
                                                sx0(istokes,i)=x(i)
                                                svx0(istokes,i)=vx(i)
                                                sy0(istokes,i)=y(i)
                                                svy0(istokes,i)=vy(i)
                                                        x(i)=sx0(istokes,i)+delta*svx0(istokes,i)/2.0d0
                                                        y(i)=sy0(istokes,i)+delta*svy0(istokes,i)/2.0d0 
                                                                if(tau==0.0d0)then
                                                                        vx(i)=ux(i)+delta*fx(i)/2.0d0
                                                                        vy(i)=uy(i)+delta*fy(i)/2.0d0
                
                                                                else
                                                                        vx(i)=svx0(istokes,i)-delta*(svx0(istokes,i)-ux(i))/(2.0d0*tau)+delta*fx(i)/2.0d0
                                                                        vy(i)=svy0(istokes,i)-delta*(svy0(istokes,i)-uy(i))/(2.0d0*tau)+delta*fy(i)/2.0d0
                                                                end if
                                        end do
                        y=modulo(y,l)
                        sxpart(istokes,:)=x
                        sypart(istokes,:)=y
                        svxpart(istokes,:)=vx
                        svypart(istokes,:)=vy
                end do
                do istokes=1,nstokes
                        x=xpart(istokes,:)
                        y=ypart(istokes,:)
                        vx=vxpart(istokes,:)
                        vy=vypart(istokes,:)
                        sizex=partsize(istokes)
                        tau=taup(istokes)
                                call linear_interp
                                        do i=1,sizex
                                                x0(istokes,i)=x(i)
                                                vx0(istokes,i)=vx(i)
                                                y0(istokes,i)=y(i)
                                                vy0(istokes,i)=vy(i)
                                                        x(i)=x0(istokes,i)+delta*vx0(istokes,i)/2.0d0
                                                        y(i)=y0(istokes,i)+delta*vy0(istokes,i)/2.0d0 
                                                                if(tau==0.0d0)then
                                                                        vx(i)=ux(i)
                                                                        vy(i)=uy(i)
                
                                                                else
                                                                        vx(i)=vx0(istokes,i)-delta*(vx0(istokes,i)-ux(i))/(2.0d0*tau)
                                                                        vy(i)=vy0(istokes,i)-delta*(vy0(istokes,i)-uy(i))/(2.0d0*tau)
                                                                end if
                                        end do
                        y=modulo(y,l)
                        xpart(istokes,:)=x
                        ypart(istokes,:)=y
                        vxpart(istokes,:)=vx
                        vypart(istokes,:)=vy
                end do       
        else
                call vel_comm
                do istokes=1,nstokes
                        x=sxpart(istokes,:)
                        y=sypart(istokes,:)
                        vx=svxpart(istokes,:)
                        vy=svypart(istokes,:)
                        sizex=spartsize(istokes)
                        tau=taup(istokes)
                        part_trac=sxpart_trac(istokes,:)
                                call linear_interp
                                call send_recv
                                call soft_pt
                               trmsint(istokes)=0.0d0
                                        do i=1,sizex 
                                                x(i)=sx0(istokes,i)+delta*vx(i)
                                                y(i)=sy0(istokes,i)+delta*vy(i)
                                                        trmsint(istokes)=trmsint(istokes)+dsqrt((delta*vx(i))**2.0d0+(delta*vy(i))**2.0d0)

                                                        if(tau==0.0d0)then
                                                                vx(i)=ux(i)-delta*fx(i)
                                                                vy(i)=uy(i)-delta*fy(i)
                                                        else
                                                                vx(i)=svx0(istokes,i)-delta*(vx(i)-ux(i))/tau+delta*fx(i)
                                                                vy(i)=svy0(istokes,i)-delta*(vy(i)-uy(i))/tau+delta*fy(i)
                                                        end if
                                        end do

                        call send_recv_x
                        sxpart(istokes,:)=x
                        sypart(istokes,:)=y
                        svxpart(istokes,:)=vx
                        svypart(istokes,:)=vy
                        sxpart_trac(istokes,:)=part_trac
                        spartsize(istokes)=sizex
                        
                end do
                do istokes=1,nstokes
                        x=xpart(istokes,:)
                        y=ypart(istokes,:)
                        vx=vxpart(istokes,:)
                        vy=vypart(istokes,:)
                        sizex=partsize(istokes)
                        tau=taup(istokes)
                        part_trac=xpart_trac(istokes,:)
                                call linear_interp
                                trmsnint(istokes)=0.0d0
                                        do i=1,sizex 
                                                x(i)=x0(istokes,i)+delta*vx(i)
                                                y(i)=y0(istokes,i)+delta*vy(i)
                                                trmsnint(istokes)=trmsnint(istokes)+dsqrt((delta*vx(i))**2.0d0+(delta*vy(i))**2.0d0)
                                                        if(tau==0.0d0)then
                                                                vx(i)=ux(i)
                                                                vy(i)=uy(i)
                                                        else
                                                                vx(i)=vx0(istokes,i)-delta*(vx(i)-ux(i))/tau
                                                                vy(i)=vy0(istokes,i)-delta*(vy(i)-uy(i))/tau
                                                        end if
                                        end do

                        call send_recv_x
                        xpart(istokes,:)=x
                        ypart(istokes,:)=y
                        vxpart(istokes,:)=vx
                        vypart(istokes,:)=vy
                        xpart_trac(istokes,:)=part_trac
                        partsize(istokes)=sizex
                        
                end do

        end if
end if

do i=1,nlocal 
        do  j=1,n
        ukx(j,i)=ukx(j,i)*omega(j,i)
        uky(j,i)=uky(j,i)*omega(j,i)
        end do
end do

call fftw_mpi_execute_dft_r2c(fplan,omega,fomega)
call fftw_mpi_execute_dft_r2c(fplan,ukx,fukx)
call fftw_mpi_execute_dft_r2c(fplan,uky,fuky)

end subroutine jacb
