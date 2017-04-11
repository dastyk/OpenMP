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

#define MAX_SIZE 4096
#define EVEN_TURN 0 /* shall we calculate the 'red' or the 'black' elements */
#define ODD_TURN  1

typedef double matrix[MAX_SIZE+2][MAX_SIZE+2]; /* (+2) - boundary elements */

struct Options {
    int		N;		/* matrix size		*/
    int		maxnum;		/* max number of element*/
    char	*Init;		/* matrix init type	*/
    double	difflimit;	/* stop condition	*/
    double	w;		/* relaxation factor	*/
    int		PRINT;		/* print switch		*/
    matrix	A;		/* matrix A		*/
};

/* forward declarations */
int work(struct Options* options);
void Init_Matrix(struct Options* options);
void Print_Matrix(struct Options* options);
void Init_Default(struct Options* options);
int Read_Options(int, char **, struct Options* options);

int main(int argc, char **argv)
{
    int i, timestart, timeend, iter;
 
    struct Options* options = malloc(sizeof(struct Options));

    Init_Default(options);		/* Init default values	*/
    Read_Options(argc,argv, options);	/* Read arguments	*/
    Init_Matrix(options);		/* Init the matrix	*/
    iter = work(options);
    if (options.PRINT == 1)
		Print_Matrix();
    printf("\nNumber of iterations = %d\n", iter);
	
	
	free(options);
	
}

int work(struct Options* options)
{
    double prevmax_even, prevmax_odd, maxi, sum, w;
    int	m, n, N, i;
    int finished = 0;
    int turn = EVEN_TURN;
    int iteration = 0;

    prevmax_even = 0.0;
    prevmax_odd = 0.0;
    N = options->N;
    w = options->w;
    
    while (!finished) {
	iteration++;
	if (turn == EVEN_TURN) {
	    /* CALCULATE part A - even elements */
	    for (m = 1; m < N+1; m++) {
		for (n = 1; n < N+1; n++) {
		    if (((m + n) % 2) == 0)
			options->A[m][n] = (1 - w) * options->A[m][n] 
			    + w * (options->A[m-1][n] + options->A[m+1][n] 
				   + options->A[m][n-1] + options->A[m][n+1]) / 4;
		}
	    }
	    /* Calculate the maximum sum of the elements */
	    maxi = -999999.0;
	    for (m = 1; m < N+1; m++) {
		sum = 0.0;
		for (n = 1; n < N+1; n++)
		    sum += options->A[m][n];
		if (sum > maxi)
		    maxi = sum;
	    }
	    /* Compare the sum with the prev sum, i.e., check wether 
	     * we are finished or not. */
	    if (fabs(maxi - prevmax_even) <= options->difflimit)
		finished = 1;
	    if ((iteration%100) == 0)
		printf("Iteration: %d, maxi = %f, prevmax_even = %f\n",
		       iteration, maxi, prevmax_even);
	    prevmax_even = maxi;
	    turn = ODD_TURN;

	} else if (turn == ODD_TURN) {
	    /* CALCULATE part B - odd elements*/
	    for (m = 1; m < N+1; m++) {
		for (n = 1; n < N+1; n++) {
		    if (((m + n) % 2) == 1)
			options->A[m][n] = (1 - w) * options->A[m][n] 
			    + w * (options->A[m-1][n] + options->A[m+1][n] 
				   + options->A[m][n-1] + options->A[m][n+1]) / 4;
		}
	    }
	    /* Calculate the maximum sum of the elements */
	    maxi = -999999.0;
	    for (m = 1; m < N+1; m++) {
		sum = 0.0;
		for (n = 1; n < N+1; n++)
		    sum += options->A[m][n];	
		if (sum > maxi)			
		    maxi = sum;
	    }
	    /* Compare the sum with the prev sum, i.e., check wether 
	     * we are finished or not. */
	    if (fabs(maxi - prevmax_odd) <= options->difflimit)
		finished = 1;
	    if ((iteration%100) == 0)
		printf("Iteration: %d, maxi = %f, prevmax_odd = %f\n",
		       iteration, maxi, prevmax_odd);
	    prevmax_odd = maxi;
	    turn = EVEN_TURN;
	} else {
	    /* something is very wrong... */
	    printf("PANIC: Something is really wrong!!!\n");
	    exit(-1);
	}
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
    int i, j, N, dmmy;
 
    N = options->N;
    printf("\nsize      = %dx%d ",N,N);
    printf("\nmaxnum    = %d \n",options->maxnum);
    printf("difflimit = %.7lf \n",options->difflimit);
    printf("Init	  = %s \n",options->Init);
    printf("w	  = %f \n\n",options->w);
    printf("Initializing matrix...");
 
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i < N+2; i++) {
	for (j = 0; j < N+2; j++) {
	    options->A[i][j] = 0.0;
	}
    }
    if (strcmp(options->Init,"count") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		options->A[i][j] = (double)i/2;
	    }
	}
    }
    if (strcmp(options->Init,"rand") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		options->A[i][j] = (rand() % options->maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(options->Init,"fast") == 0) {
	for (i = 1; i < N+1; i++){
	    dmmy++;
	    for (j = 1; j < N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    options->A[i][j] = 1.0;
		else
		    options->A[i][j] = 5.0;
	    }
	}
    }
	
	
	
    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    options->A[0][0] = options->A[1][1];
    options->A[0][N+1] = options->A[1][N];
    options->A[N+1][0] = options->A[N][1];
    options->A[N+1][N+1] = options->A[N][N];
    /* fix the top and bottom rows */
    for (i = 1; i < N+1; i++) {
	options->A[0][i] = options->A[1][i];
	options->A[N+1][i] = options->A[N][i];
    }
    /* fix the left and right columns */
    for (i = 1; i < N+1; i++) {
	options->A[i][0] = options->A[i][1];
	options->A[i][N+1] = options->A[i][N];
    }

    printf("done \n\n");
    if (options->PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix(struct Options* options)
{
    int i, j, N;
 
    N = options->N;
    for (i=0; i<N+2 ;i++){
	for (j=0; j<N+2 ;j++){
	    printf(" %f",options->A[i][j]);
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
