subroutine time_march
use mod_part
use mod_2flu
implicit none
double complex :: trc
double precision :: temph
double precision :: iniamp,omegak,rk2,k,k1,k2
integer :: i,j,kx,ky,ksqr,ireal,iimag,mshl

jac_old=fomega

call jacb(1)


do i=1,nlocal
        
        kx=i+local_offset


        kx=(kx-1)-n*(kx/(n1hf+1))
        
        do j=1,n1hf  
        ky=j-1
        ksqr=kx**2+ky**2
        if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
        fomega(j,i)=cmplx(0.0d0,0.0d0)
        else
        k1=factor*dfloat(kx)
        k2=factor*dfloat(ky)
        temph=time_incr(ksqr)
        trc=cmplx(k2*dimag(fuky(j,i))+k1*dimag(fukx(j,i)),-(k2*dble(fuky(j,i))+k1*dble(fukx(j,i))))
        fomega(j,i)=temph*(jac_old(j,i)+(delta/2.0d0)*(trc+kforce(j,i)))
end if
end do
end do

call jacb(2)

do i=1,nlocal
        
        kx=i+local_offset


        kx=(kx-1)-n*(kx/(n1hf+1))
        
        do j=1,n1hf  
        ky=j-1
        ksqr=kx**2+ky**2
        if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
        fomega(j,i)=cmplx(0.0d0,0.0d0)
        else
        k1=factor*dfloat(kx)
        k2=factor*dfloat(ky)
        temph=time_incr(ksqr)
        trc=cmplx(k2*dimag(fuky(j,i))+k1*dimag(fukx(j,i)),-(k2*dble(fuky(j,i))+k1*dble(fukx(j,i))))
        fomega(j,i)=temph*(temph*jac_old(j,i)+delta*(trc+kforce(j,i)))
end if
end do
end do

end subroutine time_march
