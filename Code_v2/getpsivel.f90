subroutine getpsivel
use mod_2flu
implicit none
integer :: i,j,ksqr,kx,ky
double precision :: rk2,k1,k2
do i=1,nlocal
       kx=i+local_offset 
       kx=(kx-1)-n*(kx/(n1hf+1))
        
        do j=1,n1hf
        
        ky=j-1
        ksqr=kx**2+ky**2
        rk2=factor*factor*dfloat(ksqr)
        k1=factor*dfloat(kx)
        k2=factor*dfloat(ky)
                if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
                psi(j,i)=0.0d0
                fukx(j,i)=0.0d0
                fuky(j,i)=0.0d0
                        else
                psi(j,i)=-fomega(j,i)/rk2       
                fukx(j,i)=cmplx(k2*dimag(psi(j,i)),-k2*dble(psi(j,i)))
                fuky(j,i)=cmplx(-k1*dimag(psi(j,i)),k1*dble(psi(j,i)))
       
       
                end if        
        
end do
end do

end subroutine getpsivel
