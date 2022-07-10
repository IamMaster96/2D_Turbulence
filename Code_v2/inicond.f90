subroutine iniconf2d
use mod_2flu
implicit none
double precision :: ran
double precision :: iniamp,omegak,rk2,k
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
               
        iniamp= (rk2**2)*dexp(-rk2**2)/den_state(mshl)
        omegak=dsqrt(iniamp)
        omegak=omegak*dfloat(nn)
        
        call random_number(ran)
        
        ran=-pi+2.0d0*pi*ran
        if(ksqr > kasqr)then
        fomega(j,i)=cmplx(0.0d0,0.0d0)
        else
        fomega(j,i)=cmplx(omegak*dcos(ran),omegak*dsin(ran))
        end if
end do
end do
end subroutine iniconf2d
