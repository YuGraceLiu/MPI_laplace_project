/*************************************************
 * Laplace Serial C Version
 *
 * Temperature is initially 0.0
 * Boundaries are as follows:
 *
 *      0         T         0
 *   0  +-------------------+  0
 *      |                   |
 *      |                   |
 *      |                   |
 *   T  |                   |  T
 *      |                   |
 *      |                   |
 *      |                   |
 *   0  +-------------------+ 100
 *      0         T        100
 *
 *  John Urbanic, PSC 2014
 *
 ************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

// size of plate
#define COLUMNS    1000
#define ROWS_GLOBAL   1000
#define chunk (ROWS_GLOBAL/4)

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01


//   helper routines
void initialize(int rank_index, int num_rank);
void track_progress(int iter);
double Temperature[chunk+2][COLUMNS+2];      // temperature grid
double Temperature_last[chunk+2][COLUMNS+2]; // temperature grid from last iteration
void output(int rand_index, int iteration); 


int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    int num_rank, rank_index;
    double dt_global=10000.;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_index);

    if ( num_rank !=4 && rank_index == 0)  {
	printf("num of processes should be 4.\n");
	MPI_Finalize();
	return 1;
    }


    if (rank_index == 0){
    	printf("Maximum iterations [100-4000]?\n");
	fflush(stdout);
    	scanf("%d", &max_iterations);
    }

    MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if ( rank_index == 0)  gettimeofday(&start_time,NULL); // Unix timer }

    initialize(rank_index, num_rank);                   // initialize Temp_last including boundary conditions



    // do until error is minimal or until max steps
    while ( dt_global > MAX_TEMP_ERROR && iteration <= max_iterations ) {

        // main calculation: average my four neighbors
        for(i = 1; i <= chunk; i++) {
            for(j = 1; j <= COLUMNS; j++) {
                Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] + Temperature_last[i][j+1] + Temperature_last[i][j-1]);
            }
        }
        
	// send first and last columns between different ranks
        if (rank_index != 0 && rank_index != num_rank-1) {
            MPI_Send(&Temperature[1][1], COLUMNS, MPI_DOUBLE, rank_index-1, 1, MPI_COMM_WORLD);
            MPI_Recv(&Temperature_last[chunk+1][1], COLUMNS, MPI_DOUBLE, rank_index+1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&Temperature[chunk][1], COLUMNS, MPI_DOUBLE, rank_index+1, 0, MPI_COMM_WORLD);
            MPI_Recv(&Temperature_last[0][1], COLUMNS, MPI_DOUBLE, rank_index-1, 0, MPI_COMM_WORLD, &status);
        }
        else if (rank_index == 0){
            MPI_Recv(&Temperature_last[chunk+1][1], COLUMNS, MPI_DOUBLE, rank_index+1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&Temperature[chunk][1], COLUMNS, MPI_DOUBLE, rank_index+1, 0, MPI_COMM_WORLD);
        }
        else{
            MPI_Send(&Temperature[1][1], COLUMNS, MPI_DOUBLE, rank_index-1, 1, MPI_COMM_WORLD);
            MPI_Recv(&Temperature_last[0][1], COLUMNS, MPI_DOUBLE, rank_index-1, 0, MPI_COMM_WORLD, &status);
        }
	//if (iteration == 1 || iteration ==2) { 
	//    MPI_Barrier(MPI_COMM_WORLD);
	//    for (int ind=0; ind<=3; ind++) { 
	//    	MPI_Barrier(MPI_COMM_WORLD);
	//	if (rank_index == ind) {
	//    	    printf("iter: %d, rank_index %d, temp     [chunk-1]:          ", iteration, rank_index);
	//    	    fflush(stdout);
	//    	    for (i=COLUMNS-5; i<=COLUMNS; i++){
	//    	    	printf("%f ", Temperature[chunk-1][i]);
	//	    }
	//    	    printf("\n");
	//    	    fflush(stdout);
	//    	    printf("iter: %d, rank_index %d, temp     [chunk]:          ", iteration, rank_index);
	//    	    fflush(stdout);
	//    	    for (i=COLUMNS-5; i<=COLUMNS; i++){
	//    	    	printf("%f ", Temperature[chunk][i]);
	//	    }
	//    	    printf("\n");
	//    	    fflush(stdout);
	//	}
	//    }
	//    fflush(stdout);
	//    MPI_Barrier(MPI_COMM_WORLD);
	//    for (int ind=0; ind<=3; ind++) { 
	//        MPI_Barrier(MPI_COMM_WORLD);
	//	if (rank_index == ind) {
	//    	    printf("iter: %d, rank_index %d, temp_last[0]:    ", iteration, rank_index);
	//    	    fflush(stdout);
	//    	    for (i=COLUMNS-5; i<=COLUMNS; i++){
	//    	    	printf("%f ", Temperature_last[0][i]);
	//    	    }
	//    	    printf("\n");
	//    	    fflush(stdout);
	//    	    printf("iter: %d, rank_index %d, temp_last[1]:    ", iteration, rank_index);
	//    	    fflush(stdout);
	//    	    for (i=COLUMNS-5; i<=COLUMNS; i++){
	//    	    	printf("%f ", Temperature_last[1][i]);
	//    	    }
	//    	    printf("\n");
	//    	    fflush(stdout);
	//	}
	//    }
	//    MPI_Barrier(MPI_COMM_WORLD);
	//}
	
	if (iteration <5) 
	    output(rank_index, iteration);

        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration and find latest dt
        for(i = 1; i <= chunk; i++){
            for(j = 1; j <= COLUMNS; j++){
	      dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
	      Temperature_last[i][j] = Temperature[i][j];
            }
        }
        MPI_Reduce(&dt, &dt_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Bcast(&dt_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // periodically print test values
        if((iteration % 100) == 0 && iteration< 2000) {
	    if (rank_index == 3 )
 	        track_progress(iteration);
        }

	iteration++;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank_index == 0 ) { 
	gettimeofday(&stop_time,NULL);
        timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
        printf("\nMax error at iteration %d was %f\n", iteration-1, dt_global);
        printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    }

    MPI_Finalize();

}


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(int rank_index, int num_rank) {

    int i,j;
    
    for(i = 0; i <= chunk+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
    }

    // these boundary conditions never change throughout run

    // set  right side to a linear increase
    for (i=0; i <= chunk+1; i++) {
    	Temperature_last[i][COLUMNS+1] = (100.0/ROWS_GLOBAL)*(i + rank_index*chunk) ;
    }
    // set bottom to linear increase
    if (rank_index == num_rank -1) {
	for (j = 0; j <= COLUMNS+1; j++) 
    	    Temperature_last[chunk+1][j] = (100.0/COLUMNS)*j ;
    }
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    printf("---------- Iteration number: %d ------------\n", iteration);
    //printf("global position [750,900]: %5.2f  ", Temperature_last[750%chunk][900]);
    //printf("global position [750,900]: %5.2f  \n", Temperature[750%chunk][900]);
    for(int i = 1; i <= 5; i++) {
        printf("[%d,900]: last:%5.2f  \n", i, i, Temperature_last[i][900]);
    }
    //for(int i = 0; i <= 5; i++) {
    //    printf("[i,988]: last:%5.2f ", i, Temperature_last[i][988]);
    //}
    printf("\n");
}

void output(int rank_index, int iteration) {
    FILE* fp;
    char filename[50];

    sprintf(filename, "output%d.txt", iteration);

    for (int pe=0; pe<4; pe++){
	if (rank_index == pe) {
	    fp = fopen(filename, "a");
	    for (int y=1; y<=chunk; y++) {
		for (int x=1;x<=COLUMNS;x++) {
		    fprintf(fp, "%5.2f ", Temperature[y][x]);
		}
		fprintf(fp, "\n");
	    }
	    fflush(fp);
	    fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }
}
