subroutine global2d
use mod_2flu
implicit none
double precision :: ran
double precision :: iniamp,omegak,rk2,k,nu
integer :: i,j,kx,ky,ksqr,ireal,iimag,mshl
do i=1,nlocal
        kx=i+local_offset
        kx=(kx-1)-n*(kx/(n1hf+1))
        
        do j=1,n1hf
        
        ky=j-1
        
        ksqr=kx**2+ky**2
        rk2=factor*factor*dfloat(ksqr)
        k=dsqrt(rk2)             
        mshl=nint(k)
        nu=vis+vis2*rk2*rk2*rk2
        if(mshl>kdrag)then
        time_incr(ksqr)= dexp(-delta*(nu*rk2+mu2)/2.0d0)
        else
        time_incr(ksqr)= dexp(-delta*(nu*rk2+mu)/2.0d0)
        end if
        if(ky.eq.0.or.ky.eq.n1h)then
        den_state(mshl)=den_state(mshl)+1
        else
        den_state(mshl)=den_state(mshl)+1
        end if

end do
end do
end subroutine global2d
