all: 
	gcc -o mm_seq matmul_seq.c
	mpicc -o mm_row matmulRow_mpi.c
	mpicc -o mm_2D matmul2D_mpi.c