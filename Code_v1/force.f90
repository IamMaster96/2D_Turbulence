subroutine soft_pt
use mpi
use mod_part

implicit none

integer :: ix,iy,jx,jy,jxbin,jybin,i,j
double precision :: wi,wj,ri,rj,dx,dy,r,f,l

l=8.0d0*atan(1.0d0)

fx=0.0d0
fy=0.0d0
!colid=0
ftemp=0.0d0

do ix=0,local_x_box+1
do iy=1,nbox
                       
                do i=1,bin(0,ix,iy)
                        wi=x(bin(i,ix,iy))
                        wj=y(bin(i,ix,iy))
                                do jx=-1,1
                                do jy=-1,1
                                       jxbin=ix+jx
                                       jybin=iy+jy
                                        if(jybin.eq.0)then
                                                jybin=nbox
                                        else if(jybin.eq.nbox+1)then
                                                jybin=1
                                        end if
                                        
                                        if(jxbin.eq.-1)then
                                                jxbin=0
                                        else if(jxbin.eq.local_x_box+2)then
                                                jxbin=local_x_box+1
                                        end if
                                                if(jxbin.gt.local_x_box+1.or.jxbin.lt.0.or.jybin.lt.1.or.jybin.gt.nbox)then
                                                        print *,jxbin,jybin,x(bin(j,jxbin,jybin))
                                                end if
                                                do j=1,bin(0,jxbin,jybin)
                                                        if(bin(i,ix,iy).gt.bin(j,jxbin,jybin))then
                                                        ri=x(bin(j,jxbin,jybin))
                                                        rj=y(bin(j,jxbin,jybin))
                                                        dx=wi-ri
                                                        dy=wj-rj
                                                                if(abs(dx).gt.(l*0.5d0))then
                                                                    if(dx.lt.0.0d0)then    
                                                                    dx=l-abs(dx)
                                                                    else if(dx.gt.0.0d0)then
                                                                    dx=dx-l
                                                                    end if
                                                                end if
                                                                if(abs(dy).gt.(l*0.5d0))then
                                                                    if(dy.lt.0.0d0)then    
                                                                    dy=l-abs(dy)
                                                                    else if(dy.gt.0.0d0)then
                                                                    dy=dy-l
                                                                    end if
                                                                end if
                                                        r=dsqrt(dx**2+dy**2)
                                                                if(r.lt.meand)then
                                                                        f=1.0d0*(1.0d0-r/meand)/meand
                                                                        ftemp(1)=f*dx/r
                                                                        ftemp(2)=f*dy/r
                                                                        fx(bin(i,ix,iy))=fx(bin(i,ix,iy))+ftemp(1)
                                                                        fy(bin(i,ix,iy))=fy(bin(i,ix,iy))+ftemp(2)
                                                                        fx(bin(j,jxbin,jybin))=fx(bin(j,jxbin,jybin))-ftemp(1)
                                                                        fy(bin(j,jxbin,jybin))=fy(bin(j,jxbin,jybin))-ftemp(2)
                                                                    !    colid(bin(i,ix,iy))=colid(bin(i,ix,iy))+1
                                                                     !   colid(bin(j,jxbin,jybin))=colid(bin(j,jxbin,jybin))+1
                                                                end if
                                                        end if
                                                end do
                                end do
                                end do
                end do
end do
end do

end subroutine soft_pt
