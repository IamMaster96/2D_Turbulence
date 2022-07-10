subroutine linear_interp
use mpi
use mod_2flu
use mod_part

implicit none
integer :: intx,inty,i,intyp,intyp1
double precision ::xip1,yjp1,xi,yj

ux=0.0d0
uy=0.0d0

do i=1,sizex

xcord=x(i)
ycord=y(i)

intx=floor(xcord/dx)
inty=floor(ycord/dx)



xi=intx*dx
yj=inty*dx
xip1=xi+dx 
yjp1=yj+dx

intx=intx-local_offset
intx=intx+1
inty=inty+1

if(inty.eq.n)then
intyp=n
intyp1=1
else
intyp=inty
intyp1=inty+1

end if

if(intx.eq.nlocal)then
        uy(i)=(xip1-xcord)*(yjp1-ycord)*uky(intyp,intx)+(xip1-xcord)*(ycord-yj)*vely_recv_u_0(intyp)+(xcord-xi)*(yjp1-ycord)*uky(intyp1,intx)+(xcord-xi)*(ycord-yj)*vely_recv_u_0(intyp1)
        ux(i)=(xip1-xcord)*(yjp1-ycord)*ukx(intyp,intx)+(xip1-xcord)*(ycord-yj)*velx_recv_u_0(intyp)+(xcord-xi)*(yjp1-ycord)*ukx(intyp1,intx)+(xcord-xi)*(ycord-yj)*velx_recv_u_0(intyp1)

else if(intx.eq.nlocal+1)then
        uy(i)=(xip1-xcord)*(yjp1-ycord)*vely_recv_u_0(intyp)+(xip1-xcord)*(ycord-yj)*vely_recv_u_1(intyp)+(xcord-xi)*(yjp1-ycord)*vely_recv_u_0(intyp1)+(xcord-xi)*(ycord-yj)*vely_recv_u_1(intyp1)
        ux(i)=(xip1-xcord)*(yjp1-ycord)*velx_recv_u_0(intyp)+(xip1-xcord)*(ycord-yj)*velx_recv_u_1(intyp)+(xcord-xi)*(yjp1-ycord)*velx_recv_u_0(intyp1)+(xcord-xi)*(ycord-yj)*velx_recv_u_1(intyp1)
               

else if(intx.eq.0)then
        ux(i)=(xip1-xcord)*(yjp1-ycord)*velx_recv_l(intyp)+(xip1-xcord)*(ycord-yj)*ukx(intyp,intx+1)+(xcord-xi)*(yjp1-ycord)*velx_recv_l(intyp1)+(xcord-xi)*(ycord-yj)*ukx(intyp1,intx+1)
        uy(i)=(xip1-xcord)*(yjp1-ycord)*vely_recv_l(intyp)+(xip1-xcord)*(ycord-yj)*uky(intyp,intx+1)+(xcord-xi)*(yjp1-ycord)*vely_recv_l(intyp1)+(xcord-xi)*(ycord-yj)*uky(intyp1,intx+1)
else
        if(intx.ge.nlocal.or.intx.le.0)then
        print *, intx,xcord,rank
        end if
        uy(i)=(xip1-xcord)*(yjp1-ycord)*uky(intyp,intx)+(xip1-xcord)*(ycord-yj)*uky(intyp,intx+1)+(xcord-xi)*(yjp1-ycord)*uky(intyp1,intx)+(xcord-xi)*(ycord-yj)*uky(intyp1,intx+1)
        ux(i)=(xip1-xcord)*(yjp1-ycord)*ukx(intyp,intx)+(xip1-xcord)*(ycord-yj)*ukx(intyp,intx+1)+(xcord-xi)*(yjp1-ycord)*ukx(intyp1,intx)+(xcord-xi)*(ycord-yj)*ukx(intyp1,intx+1)

end if
ux(i)=ux(i)/dx**2
uy(i)=uy(i)/dx**2
end do


end subroutine linear_interp
