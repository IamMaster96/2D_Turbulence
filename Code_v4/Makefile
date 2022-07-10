# Makefile for compiling the Serial Spectral Code Version 3.0
# This uses the Intel Fortran Compiler, with FFTW 
# ---------------------------------------------------------------------
compile =mpifort -ffree-line-length-none -g -fbacktrace -fcheck=bounds
link = -L/home/mohit.gupta/fftw_install/fftw3/lib -I/home/mohit.gupta/fftw_install/fftw3/include -lfftw3_mpi -lfftw3 -lm  
#---------------------------------------------------------------------
OBJ_rk =  2dns.o inicond.o eval_jac.o getpsivel.o energy.o file_write.o time_march.o gen_force.o global2d.o vel_com.o part_init.o interp.o comm.o force.o send_recv.o
MOD_FILE_1=fftw_mod.o
MOD_FILE_2=mod_2flu.o
MOD_FILE_3=mod_part.o
test : 
	echo $(OBJ_FILE) 
# -------------------------------------------------------------
nsrk.exe:$(OBJ_rk)  $(MOD_FILE_1) $(MOD_FILE_2) $(MOD_FILE_3)
	$(compile) $(OBJ_rk) $(MOD_FILE_2) $(MOD_FILE_1) $(MOD_FILE_3) -o nsrk.exe $(link)
#
2dns.o: 2dns.f90 $(MOD_FILE_1) $(MOD_FILE_2) $(MOD_FILE_3)
	$(compile) -c 2dns.f90 $(MOD_FILE_2) $(MOD_FILE_1) $(MOD_FILE_3)
#
fftw_mod.o: fftw_mod.f90
	$(compile) -c fftw_mod.f90 $(link)
#
energy.o: energy.f90
	$(compile) -c energy.f90 $(link)
#
part_init.o: part_init.f90
	$(compile) -c part_init.f90 $(link)
#
time_march.o: time_march.f90
	$(compile) -c time_march.f90 $(link)
#
comm.o: comm.f90
	$(compile) -c comm.f90 $(link)
#
interp.o: interp.f90
	$(compile) -c interp.f90 $(link)
#
vel_com.o: vel_com.f90
	$(compile) -c vel_com.f90 $(link)
#
file_write.o: file_write.f90
	$(compile) -c file_write.f90
#
force.o: force.f90
	$(compile) -c force.f90 $(link)
#
send_recv.o: send_recv.f90
	$(compile) -c send_recv.f90 $(link)
#
eval_jac.o: eval_jac.f90
	$(compile) -c eval_jac.f90 $(link)
#
getpsivel.o: getpsivel.f90
	$(compile) -c getpsivel.f90 $(link)
#
global2d.o: global2d.f90
	$(compile) -c global2d.f90 $(link)
#
gen_force.o: gen_force.f90
	$(compile) -c gen_force.f90 $(link)
#
inicond.o: inicond.f90
	$(compile) -c inicond.f90 $(link)
#
mod_part.o: mod_part.f90
	$(compile) -c mod_part.f90
#
mod_2flu.o: mod_2flu.f90
	$(compile) -c mod_2flu.f90
#
#
#
# ------------------------------------------------------------------------ 
clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe
data_clean:
	rm -f *.dat *.out data_psi/* spectras/* 
## --------------------------------------------------------------------##  
