module mod_part
use mpi
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision,allocatable,dimension(:) :: x,y
double precision,allocatable,dimension(:) :: vx,vy,ux,uy,fx,fy,ftemp
double precision,allocatable,dimension(:) :: xtemp,ytemp
double precision,allocatable,dimension(:) :: vxtemp,vytemp
double precision,allocatable,dimension(:,:) ::xpart,ypart,vxpart,vypart,x0,y0,vx0,vy0,sxpart,sypart,svxpart,svypart,sx0,sy0,svx0,svy0
double precision,allocatable,dimension(:) :: taup,trmsint,trmsnint,trmsint_out,trmsnint_out
integer,allocatable,dimension(:) :: partsize,spartsize,part_trac,part_temp,send_utrac,send_ltrac
integer,allocatable,dimension(:,:):: xpart_trac,sxpart_trac
!!for interactions
integer,allocatable,dimension(:,:,:)::bin
!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: rank_l,rank_u
double precision :: cord_l,cord_u,xcord,ycord
double precision,allocatable,dimension(:) ::send_ly,send_lx,recv_ly,recv_lx
double precision,allocatable,dimension(:) ::send_uy,send_ux,recv_uy,recv_ux
double precision,allocatable,dimension(:) :: send_lvy,send_lvx,send_uvy,send_uvx
double precision,allocatable,dimension(:) :: velx_recv_u_0,velx_recv_u_1,velx_recv_l,vely_recv_u_0,vely_recv_u_1,vely_recv_l
integer :: istatus1(mpi_status_size),istatus2(mpi_status_size),istatus3(mpi_status_size),istatus4(mpi_status_size)
integer :: istatus8(mpi_status_size),istatus7(mpi_status_size),istatus6(mpi_status_size),istatus5(mpi_status_size)
integer :: istatus9(mpi_status_size),istatus10(mpi_status_size)
integer :: istatus12(mpi_status_size),istatus11(mpi_status_size)
integer :: ireqly,ireqlx,irequy,irequx,ireq1,ireq2,ireq3,ireq4,ireq5,ireq6,ireqltrac,irequtrac
integer :: ireqvly,ireqvlx,ireqvuy,ireqvux
!!!!!!!!!!!!!!!!!!!!!
integer ::nparticle,nstokes,istokes
double precision :: nbin,tau
integer :: sizex,local_nparticle
logical :: particles
double precision :: meand,meand2
integer :: nbox,local_x_box,bin_cap,nrunpart
character(100) :: st
end module
