subroutine energy_2d
use mpi
use mod_2flu
implicit none
integer :: i,j,ireal,iimag,ksqr,mshl,kx,ky
double precision :: k,rk2,irank
double precision :: energy_norm,energy,one_by_nsqr

E_omega=0.0d0
one_by_nsqr=1.0d0/dfloat(nn)
energy_norm=one_by_nsqr**2
do i=1,nlocal
        kx=i+local_offset
        kx=(kx-1)-n*(kx/(n1hf+1))
        do j=1,n1hf
                ky=j-1
                ksqr=kx**2+ky**2
                rk2=factor*factor*dfloat(ksqr)
                k=dsqrt(rk2)                
                mshl=idnint(k)
                energy=energy_norm*(dble(fukx(j,i))**2+dimag(fukx(j,i))**2+dble(fukx(j,i))**2+dimag(fuky(j,i))**2)
                if(ky==0.or.ky==n1h)then
                E_omega(mshl)=E_omega(mshl)+energy
                else
                E_omega(mshl)=E_omega(mshl)+2.0d0*energy
                end if
        end do
end do


end subroutine energy_2d
