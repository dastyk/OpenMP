all: 
	gcc -o mm_seq matmul_seq.c
	mpicc -o mm_row matmulRow_mpi.c
	mpicc -o mm_2D matmul2D_mpi.c
	
	gcc -o sor_seq sor_seq.c
	mpicc -o sor_par sor_par.c
	
clean:
	rm mm_*
	rm sor_par
	rm sor_seq