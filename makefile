all: 
	mpicc -o mm_2D matmul2D_mpi.c
	mpicc -o sor_par sor_par.c
	
clean:
	rm sor_par
	rm mm_2D