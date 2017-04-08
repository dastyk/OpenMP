
/*
 * Compile with:
 * mpicc -o mm matmul_mpi.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 1024	/* assumption: SIZE a multiple of number of nodes */
        /* SIZE should be 1024 in our measurements in the assignment */
        /* Hint: use small sizes when testing, e.g., SIZE 8 */
#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */
#define DEBUG	1	/* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
            /* Simple initialization, which enables us to easy check
             * the correct answer. Given SIZE size of the matrices, then
             * the output should be
             *     SIZE ... 2*SIZE ...
             *     ...
             *     2*SIZE ... 4*SIZE ...
             *     ...
             */
            a[i][j] = 1.0;
            if (i >= SIZE/2) a[i][j] = 2.0;
            b[i][j] = 1.0;
            if (j >= SIZE/2) b[i][j] = 2.0;
        }
}

static void
print_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
            printf(" %7.2f", c[i][j]);
        printf("\n");
    }
}

void SendBlock(void* data, int x, int y, int cols, int rows, int dest, int tag)
{
	int offset;
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Send(&a[x][y + offset], cols, MPI_INT, dest, tag, MPI_COMM_WORLD);
	}	
}

void RecvBlock(void* data, int x, int y, int cols, int rows, int src, int tag)
{
	int offset;
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Recv(&a[x][y + offset], cols, MPI_INT, src, tag, MPI_COMM_WORLD);
	}	
}




int main(int argc, char **argv)
{
	int i, j;
	int x,y;
	int myrank, numNodes;
	int dest, src, offset;
	MPI_Init(&argc, &argv);
	double start_time, end_time;
	int px, py, cx, cy;
	MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	px = numNodes / 2;
	py = 2;
	cx = SIZE / px;
	cy = SIZE / py;
	
	if(myrank == 0) // Master
	{
		printf("SIZE = %d, Number of columns: %d, Number of rows: %d", SIZE, px, py);
		init_matrix();
		start_time = MPI_Wtime();
		
		

		for(y = 0; y < py; y++)
		{
			for(x = 0; x < px; x++)
			{
				if(x*y != 0)
				{
					dest = y*cy + x;
#ifdef DEBUG
					printf("Sending %d colums and %d rows to node %d", cx,cy,dest);
#endif
					
					SendBlock(a, cx*x, cy*y, cx, cy, dest, FROM_MASTER);
				}
			}
		}
		
		
		
	}
	else
	{
		int a_l[cx][cy];
		RecvBlock(a, 0, 0, cx, cy, 0, FROM_MASTER);
#ifdef DEBUG
					printf("Node %d recvied...\n", myrank);
					for (i = 0; i < cx; i++) 
					{
						for (j = 0; j < cy; j++)
							printf(" %7.2f", c[i][j]);
						printf("\n");
					}
#endif		
		
		
		
	}
	
	
	
	
	MPI_Finalize();
	
	
	
    return 0;
}
