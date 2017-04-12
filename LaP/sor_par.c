/*****************************************************
 *
 * S O R algorithm
 * ("Red-Black" solution to LaPlace approximation)
 *
 * sequential version
 *
 *****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#define MAX_SIZE 4096
#define EVEN_TURN 0 /* shall we calculate the 'red' or the 'black' elements */
#define ODD_TURN  1

#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */

typedef double matrix[MAX_SIZE+2][MAX_SIZE+2]; /* (+2) - boundary elements */

MPI_Status status;

struct Options {
    int		N;		/* matrix size		*/
    int		maxnum;		/* max number of element*/
    char	*Init;		/* matrix init type	*/
    double	difflimit;	/* stop condition	*/
    double	w;		/* relaxation factor	*/
    int		PRINT;		/* print switch		*/
    double*	A;		/* matrix A		*/
};

struct SlimOptions // Sent to worker nodes, as we don't need everything.
{
	int N;
	double difflimit;
	double w;
};

void SendOptions(struct Options* options)
{
	printf("Sending N=%d, difflimit=%f, w=%f\n", options->N, options->difflimit, options->w);
	
	MPI_Bcast(&options->N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&options->difflimit, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void RecvOptions(struct SlimOptions* options)
{
	
	MPI_Bcast(&options->N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&options->difflimit, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	printf("Recieved N=%d, difflimit=%f, w=%f\n", options->N, options->difflimit, options->w);
}
void SendBlock(double* data, int x, int y, int cols, int rows, int stride, int dest, int tag)
{
	int offset;
	data += stride*y+x;
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Send(data, cols, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);	
		data += stride;		
	}	
}

void RecvBlock(double* data, int x, int y, int cols, int rows, int stride, int src, int tag)
{
	int offset;
	data += stride*y+x;
	for(offset = 0; offset < rows; offset++)
	{
		MPI_Recv(data, cols, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);	
		data += stride;		
	}
}

/* forward declarations */
int work(int N, double w, double difflimit, double* A, int stride, int myrank, int numNodes);
void Init_Matrix(struct Options* options);
void Print_Matrix(struct Options* options);
void Init_Default(struct Options* options);
int Read_Options(int, char **, struct Options* options);
int Master(struct Options* options, int numNodes);
void Worker(int numNodes, int myrank);

int main(int argc, char **argv)
{
    int i, timestart, timeend, iter;
	int myrank, numNodes;
	int dest, src, offset;
    struct Options options;
	
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if(myrank == 0)
	{
		
		
		Init_Default(&options);		/* Init default values	*/
		Read_Options(argc,argv, &options);	/* Read arguments	*/
		Init_Matrix(&options);		/* Init the matrix	*/
		
		iter = Master(&options, numNodes);			/*Send work to each worker*/
		

		sleep(1);
		
		//iter = work(options);
		if (options.PRINT == 1)
			Print_Matrix(&options);
		printf("\nNumber of iterations = %d\n", iter);
	
	
	}
	else
	{
		Worker(numNodes, myrank);
		
	}
	
	MPI_Finalize();
}
int Master(struct Options* options, int numNodes)
{
	int i, iter,stride;
	int rowsPP = options->N / numNodes;
	stride = options->N + 2;
	SendOptions(options);

	for (i = 1; i < numNodes; i++)
	{
		SendBlock(options->A, 
		0, i*rowsPP + 1, // send the correct rows
		stride, rowsPP + (i == numNodes -1 ? 1 : 0), // One extra if last node
		stride, 
		i, FROM_MASTER);
		//sleep(1);
		
	}
	iter = work(options->N, options->w, options->difflimit, options->A, stride, 0, numNodes);
	
	// Get data from workers
	MPI_Gather(&options->A[1*stride], stride * rowsPP , MPI_DOUBLE, 
	&options->A[1*stride], rowsPP*stride, MPI_DOUBLE,
	0, MPI_COMM_WORLD);
	
	return iter;
	
}
void Worker(int numNodes, int myrank)
{
	int i;
	int rowsPP, stride;
	struct SlimOptions options;
	double* mat;

	RecvOptions(&options);
	
	rowsPP = options.N / numNodes;
	stride = options.N + 2; 	// Two extra cols for halo elements

	mat = malloc(stride*(rowsPP + 2)*sizeof(double));	// Two extra rows for halo elements
	
	
	RecvBlock(mat, 
	0, 1, /*Offset y by 1(the halo element)*/ 
	stride, rowsPP + (myrank == numNodes - 1 ? 1 : 0), /*one extra row if we are the last node*/
	stride, 
	0, FROM_MASTER);
  /*int x,y;

    for (y = 0; y < rowsPP + 2; y++){
        for (x = 0; x < stride; x++) 
            printf(" %7.2f", mat[y*stride + x]);
        printf("\n");
    }
	*/
	// Do the calcs
	work(options.N, options.w, options.difflimit, mat, stride, myrank, numNodes);
	
	
	// Send data to master
	MPI_Gather(&mat[1*stride], stride * rowsPP, MPI_DOUBLE, 
	mat, stride * rowsPP, MPI_DOUBLE,
	0, MPI_COMM_WORLD);
	
	
	free(mat);
	
}
int work(int N, double w, double difflimit, double* A, int stride, int myrank, int numNodes)
{
    double prevmax[2], maxi, sum, maxiall;
    int	m, n, i, rowsPP;
    int finished = 0;
    int turn = EVEN_TURN;
    int iteration = 0;

    prevmax[EVEN_TURN] = 0.0;
    prevmax[ODD_TURN] = 0.0;

    rowsPP = N / numNodes;
	
	
	/*printf("Node %d, Rows per node %d\n", myrank, rowsPP);
	   int x,y;

    for (y = 0; y < rowsPP + 2; y++){
        for (x = 0; x < N + 2; x++) 
            printf(" %7.2f", A[y*stride + x]);
        printf("\n");
    }*/
	
	
    while (!finished) {
	iteration++;
	
	
	
	//printf("Node %d, iter %d", myrank, iteration);
	
	// Send the halo elements, 
	if(numNodes > 1)
	{
		if(myrank == 0) // End node
		{
			SendBlock(A, 0, rowsPP, stride, 1, stride, 1, FROM_WORKER); // Send my halo elements to bottom neighbor 
			RecvBlock(A, 0, rowsPP + 1, stride, 1, stride, 1, FROM_WORKER); // Recv halo elements from bottom neighbor.
			
		}
		else if(myrank == numNodes - 1) // End node
		{
			RecvBlock(A, 0, 0, stride, 1, stride, myrank - 1, FROM_WORKER); // Recv from top neighbor
			SendBlock(A, 0, 1, stride, 1, stride, myrank - 1, FROM_WORKER); // Send to top neighbor	
		}
		else
		{
			RecvBlock(A, 0, 0, stride, 1, stride, myrank - 1, FROM_WORKER); // Recv from top neighbor
			SendBlock(A, 0, rowsPP, stride, 1, stride, myrank + 1, FROM_WORKER); // Send to bottom neighbor 
			RecvBlock(A, 0, rowsPP + 1, stride, 1, stride, myrank + 1, FROM_WORKER); // Recv from bottom neighbor
			SendBlock(A, 0, 1, stride, 1, stride, myrank - 1, FROM_WORKER); // Send to top neighbor			
		}
	}
	
	
	// Now that we have the halo elements we can do work.
	
	
	    /* CALCULATE */
	    for (m = 1; m < rowsPP + 1; m++) {
			for (n = 1; n < N+1; n++) {
				if (((m + n) % 2) == turn)
					A[m*stride + n] = (1 - w) * A[m*stride + n] 
						+ w * (A[(m-1)*stride + n] + A[(m+1)*stride + n] 
						+ A[m*stride + n-1] + A[m*stride + n+1]) / 4;
			}
	    }
	    /* Calculate the maximum sum of the elements */
	    maxi = -999999.0;
	    for (m = 1; m < rowsPP+1; m++) {
			sum = 0.0;
			for (n = 1; n < N+1; n++)
				sum += A[m*stride + n];
			if (sum > maxi)
				maxi = sum;
	    }
		
		printf("Node %d, maxi %f\n", myrank, maxi);
		
		// Now we need to share the maximum with all other nodes to see if we are finished
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&maxi, &maxiall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		printf("Node %d, maxiall %f\n", myrank, maxiall);
		
		

	   if(myrank == 0) // Only let master print this
	   {
		    if ((iteration%100) == 0)
			printf("Iteration: %d, maxiall = %f, prevmax[%s] = %f\n",
		       iteration, maxiall, turn ? "ODD_TURN" : "EVEN_TURN", prevmax[turn]);
	   }
	   
	   	/* Compare the sum with the prev sum, i.e., check wether 
	     * we are finished or not. */
	    if (fabs(maxiall - prevmax[turn]) <= difflimit)
		finished = 1;

	    prevmax[turn] = maxiall;
	    turn = (turn + 1) % 2;

	
	if (iteration > 100000) {
	    /* exit if we don't converge fast enough */
	    printf("Max number of iterations reached! Exit!\n");
	    finished = 1;
	}
    }

    return iteration;
}

/*--------------------------------------------------------------*/

void Init_Matrix(struct Options* options)
{
    int i, j, N, dmmy, stride;
 
    N = options->N;
	stride = N+ 2;
    printf("\nsize      = %dx%d ",N,N);
    printf("\nmaxnum    = %d \n",options->maxnum);
    printf("difflimit = %.7lf \n",options->difflimit);
    printf("Init	  = %s \n",options->Init);
    printf("w	  = %f \n\n",options->w);
    printf("Initializing matrix...");
 
 
	options->A = malloc(stride*stride*sizeof(double));
 
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i < N+2; i++) {
	for (j = 0; j < N+2; j++) {
	    options->A[i*stride + j] = 0.0;
	}
    }
    if (strcmp(options->Init,"count") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		options->A[i*stride + j] = (double)i/2;
	    }
	}
    }
    if (strcmp(options->Init,"rand") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		options->A[i*stride + j] = (rand() % options->maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(options->Init,"fast") == 0) {
	for (i = 1; i < N+1; i++){
	    dmmy++;
	    for (j = 1; j < N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    options->A[i*stride + j] = 1.0;
		else
		    options->A[i*stride + j] = 5.0;
	    }
	}
    }
	
	
	
    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    options->A[0* stride + 0] = options->A[1* stride + 1];
    options->A[0*stride + N+1] = options->A[1*stride + N];
    options->A[(N+1)*stride + 0] = options->A[N*stride + 1];
    options->A[(N+1)*stride + N+1] = options->A[N*stride + N];
    /* fix the top and bottom rows */
    for (i = 1; i < N+1; i++) {
	options->A[0*stride + i] = options->A[1* stride + i];
	options->A[(N+1)*stride + i] = options->A[N*stride + i];
    }
    /* fix the left and right columns */
    for (i = 1; i < N+1; i++) {
	options->A[i*stride + 0] = options->A[i*stride + 1];
	options->A[i* stride + N+1] = options->A[i*stride+ N];
    }

    printf("done \n\n");
    if (options->PRINT == 1)
		Print_Matrix(options);
}

void
Print_Matrix(struct Options* options)
{
    int i, j, N, stride;
 
    N = options->N;
	stride = N + 2;
    for (i=0; i<N+2 ;i++){
	for (j=0; j<N+2 ;j++){
	    printf(" %f",options->A[i*stride + j]);
	}
	printf("\n");
    }
    printf("\n\n");
}

void 
Init_Default(struct Options* options)
{
    options->N = 2048;
    options->difflimit = 0.00001*options->N;
    options->Init = "rand";
    options->maxnum = 15.0;
    options->w = 0.5;
    options->PRINT = 0;
}
 
int
Read_Options(int argc, char **argv, struct Options* options)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		options->N = atoi(*++argv);
		options->difflimit = 0.00001*options->N;
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-d difflimit] 0.1-0.000001 \n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand/count \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		printf("           [-w relaxation_factor] 1.0-0.1 \n\n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", options->N);
		printf("\n          difflimit = 0.0001 ");
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          w         = 0.5 \n");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		options->Init = *++argv;
		break;
	    case 'm':
		--argc;
		options->maxnum = atoi(*++argv);
		break;
	    case 'd':
		--argc;
		options->difflimit = atof(*++argv);
		break;
	    case 'w':
		--argc;
		options->w = atof(*++argv);
		break;
	    case 'P':
		--argc;
		options->PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
