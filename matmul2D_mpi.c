
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
//#define DEBUG	

MPI_Status status;


static void
init_matrix(double* mat, int size)
{
    int x,y;

    for (y = 0; y < size; y++)
        for (x = 0; x < size; x++) {

            mat[y*size + x] = 1.0;
            if (x >= size/2) mat[y*size + x] = 2.0;

        }
}

static void
print_matrix(double* mat, int col, int row, int stride)
{
    int x,y;

    for (y = 0; y < row; y++){
        for (x = 0; x < col; x++) 
            printf(" %7.2f", mat[y*stride + x]);
        printf("\n");
    }
}
/*Send data into the data pointer, with an offset of x,y and a total of cols column and rows rows. Stride is the number if columns in the matrix.*/
void SendBlock(double* data, int x, int y, int cols, int rows, int stride, int dest, int tag)
{
	int offset;
	data += stride*y+x;
		#ifdef DEBUG
		printf("Sending %d colums and %d rows to node %d, with offsets %d, %d\n", cols, rows, dest, x,y); 
	//	print_matrix(data, cols, rows, stride); 
	#endif
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Send(data, cols, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);	
		data += stride;		
	}	
}


/*Recv data into the data pointer, with an offset of x,y and a total of cols column and rows rows. Stride is the number if columns in the matrix.*/
void RecvBlock(double* data, int x, int y, int cols, int rows, int stride, int src, int tag)
{
	int offset;
	#ifdef DEBUG
	double* temp = data;
	#endif
	data += stride*y+x;
	#ifdef DEBUG
		printf("Receiving %d colums and %d rows from node %d, with offsets %d, %d\n", cols, rows, src, x,y); 
	#endif
	
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Recv(data, cols, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);	
		data += stride;		
	}
	#ifdef DEBUG
	//	printf("Data received\n"); 
	//	print_matrix(temp, cols, rows, stride); 
	#endif
}


int main(int argc, char **argv)
{
	int i, j;
	int x,y;
	int myrank, numNodes;
	int dest, src, offset;
	
	
	double start_time, end_time;
	int px, py, cx, cy;
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	
	// Handle the first two cases a bit different to the rest.
	if(numNodes == 1)
		px = 1;
	else if(numNodes == 2)
		px = 1;
	else
		px = 2; // Then we divde the data once in the middle along the columns.
	py = numNodes  / px; // Number of processors in the y direction.

	cx = SIZE / px;	// Number of columns per process
	cy = SIZE / py; // Number of rows per process.
	
	if(myrank == 0) // Master
	{
		// Allocate the main matricies
		double* a = malloc(SIZE*SIZE*sizeof(double));
		init_matrix(a, SIZE);
		double* b = malloc(SIZE*SIZE*sizeof(double));
		init_matrix(b, SIZE);
		double* c = malloc(SIZE*SIZE*sizeof(double));
		

		printf("Num Nodes: %d\n", numNodes);

		
		start_time = MPI_Wtime();
		
		// Send one block to each worker
		for(y = 0; y < py; y++)
		{
			for(x = 0; x < px; x++)
			{
				dest = y*px + x;
				if(dest != 0) // Skip the master
				{			
					// SendBlock(/*matrix*/,/*x offset*/, /*y offset*/, /*cols*/, /*rows*/, /*Matrix stride*/...)
					SendBlock(a, 0, y*cy, SIZE, cy, SIZE, dest, FROM_MASTER); // Send the rows of A needed to calculate the area (x*cx, y*cy) - (x*cx +cx, y*cy + cy)
					SendBlock(b, x*cx, 0, cx, SIZE, SIZE, dest, FROM_MASTER); // Send the coloums of B needed. 
				}
			}
		}
		
		// Master calc
		for(y = 0; y < cy; y++)
		{
			for(x = 0; x < cx; x++)
			{
				double sum = 0.0f;
				for(i = 0; i < SIZE; i++)
				{
					sum += a[y*SIZE + i] * b[i*SIZE + x];
				}
				c[y*SIZE + x] = sum;
			}
		}
		
		#ifdef DEBUG
		printf("Node %d done\n", myrank);
		#endif
		
		//Collect the work from the children
		for(y = 0; y < py; y++)
		{
			for(x = 0; x < px; x++)
			{
				src = y*px + x;
				if(src != 0)
				{				
					RecvBlock(c, x*cx, y*cy, cx, cy, SIZE, src, FROM_WORKER);
				}
			}
		}
		
		end_time = MPI_Wtime();
		#ifdef DEBUG
		printf("\nA: \n");
		print_matrix(a, SIZE, SIZE, SIZE);
		printf("\nB: \n");
		print_matrix(b, SIZE, SIZE, SIZE);
		printf("\nC: \n");
		print_matrix(c, SIZE, SIZE, SIZE);
		#endif
		printf("\n\nExecution time on %2d nodes: %f\n", numNodes, end_time-start_time);
		
		free(a);
		free(b);
		free(c);
	}
	else // Workers
	{
		// Allocate local space for the worker.
		double* a = malloc(SIZE*cy*sizeof(double));
		double* b = malloc(cx*SIZE*sizeof(double));
		double* c = malloc(cx*cy*sizeof(double));
		
		// Collect work
		RecvBlock(a, 0, 0,SIZE,cy, SIZE, 0, FROM_MASTER);
		RecvBlock(b, 0, 0,cx,SIZE, cx, 0, FROM_MASTER);
	
		// Do calcs
		for(y = 0; y < cy; y++)
		{
			for(x = 0; x < cx; x++)
			{
				double sum = 0.0f;
				for(i = 0; i < SIZE; i++)
				{
					sum += a[y*SIZE + i] * b[i*cx + x];
				}
				c[y*cx + x] = sum;
			}
		}
		
		#ifdef DEBUG
		printf("Node %d done\n", myrank);
		#endif
		
		// Send local results to master
		SendBlock(c, 0 ,0, cx, cy, cx, 0, FROM_WORKER);
		
		
		free(a);
		free(b);
		free(c);
	}
	
	
	
	MPI_Finalize();
	
	
	
    return 0;
}
