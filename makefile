all: mm_2D sor_par
	
mm1: mm_2D
	mpirun -n 1 ./mm_2D
mm2: mm_2D
	mpirun -n 2 ./mm_2D
mm4: mm_2D
	mpirun -n 4 ./mm_2D
mm8: mm_2D
	mpirun -n 8 ./mm_2D
	
sor1: sor_par
	mpirun -n 1 ./sor_par -n 2048
sor2: sor_par
	mpirun -n 2 ./sor_par -n 2048
sor4: sor_par
	mpirun -n 4 ./sor_par -n 2048
sor8: sor_par
	mpirun -n 8 ./sor_par -n 2048
	
mm_2D:
	mpicc -o mm_2D matmul2D_mpi.c
	
sor_par:
	mpicc -o sor_par sor_par.c
clean:
	rm sor_par
	rm mm_2D