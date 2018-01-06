#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "Utils/barrier.h"

#define NUMTHREADS 2

int n;
double *mat, *matOriginal;
pthread_barrier_t mybarrier;

struct arg_struct {
    int rank;
};

void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}


void cholesky(int rank) {

    printf("Cholesky\n");
	double start, end;
	int i, j, k;
	double *copy = calloc(n, sizeof(double));



	for(j = 0; j < n; j++) {

		/*
		 * Step 0:
		 * Replace the entries above the diagonal with zeroes
		 */
		if (rank == 0) {
			for (i = 0; i < j; i++) {
                //L[i][j] = 0.0;
				mat[i * n + j] = 0;
			}
		}

		/*
		 * Step 1:
		 * Update the diagonal element
		 */

		if (j % NUMTHREADS == rank) {

			for (k = 0; k < j; k++) {
				//L[j][j] = L[j][j] - L[j][k] * L[j][k];
                mat[j * n + j] = mat[j * n + j] - mat[j * n + k] * mat[j * n + k];
			}

			//L[j][j] = sqrt(L[j][j]);
            mat[j * n + j] = sqrt(mat[j * n + j]);
		}

        pthread_barrier_wait(&mybarrier);

		/*
		 * Step 2:
		 * Update the elements below the diagonal element
		 */

		// Divide the rest of the work
		for (i = j + 1; i < n; i++) {
			if (i % NUMTHREADS == rank) {
				for (k = 0; k < j; k++) {
					//L[i][j] = L[i][j] - L[i][k] * L[j][k];
                    mat[i * n + j] = mat[i * n + j] - mat[i * n + k] * mat[j * n + k];
				}

				//L[i][j] = L[i][j] / L[j][j];
                mat[i * n + k] = mat[i * n + k] / mat[j * n + j];
			}
		}
	}

    /*
	if (rank == 0){
		end = MPI_Wtime();
		printf("Testing OpenMpi implementation Output: \n");
		printf("Runtime = %lf\n", end-start);
		printf("Testing MPI implementation Output: ");
		testBasicOutput(A,L,n);
        // Test
        /*double ** LLT = matrixMultiply(L, transpose(L, n), n);
        printf("L*L^T = \n");
        print(LLT, n);

	}
    */

    //printf("matIN\n");
    //show_matrix(mat, n);

}

void *choleskyThread(void *arguments) {
    struct arg_struct *args = arguments;
    cholesky(args->rank);
}

void verifyCholesky(double *mat, double *a, int n) {
    printf("am intrat\n");
    int i, j, k;
    double rez[n * n], sum;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum = 0;
            for(k = 0; k < n; k++) {
                sum += a[i * n + k] * a[j * n + k];
            }
            rez[i*n+j] = sum;
        }
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(round(mat[i*n+j]) != round(rez[i*n+j])) {
                printf("%d %d\n", i, j);
                printf("%lf -- %lf\n", mat[i*n+j], rez[i*n+j]);
                printf("sugi o ceapa\n");
                return;
            }
        }
    }

    printf("sugi un cartof\n");
}



int main(int argc, char *argv[]) {
    int i, j, iret;
    FILE *f = NULL;

    pthread_t thread[NUMTHREADS];
    struct arg_struct args_s[NUMTHREADS];

    f = fopen("testFile4.txt", "r");
    fscanf(f, "%d", &n);
    printf("%d\n", n);


    mat = calloc(n*n, sizeof(double));
    matOriginal = calloc(n*n, sizeof(double));

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            fscanf(f, "%lf", &mat[i * n + j]);
        }
    }
    printf("barabula\n");
    fclose(f);

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            matOriginal[i * n + j] = mat[i * n + j];
        }
    }



	// for(i = 0; i < n; i++) {
	// 	for(j = 0; j < n; j++) {
	// 		printf("%lf ", mat[i * n + j]);
	// 	}
	// 	printf("\n");
	// }

    //show_matrix(mat, n);
    //show_matrix(matOriginal, n);

    for (i = 0; i < NUMTHREADS; ++i) {
        args_s[i].rank = i;
    }

    pthread_barrier_init(&mybarrier, NULL, NUMTHREADS);

    // create threads
    for (i = 0; i < NUMTHREADS; ++i) {
        iret = pthread_create(&thread[i], NULL, choleskyThread, (void *)&args_s[i]);
        if(iret) {
            fprintf(stderr,"Error - pthread_create() return code: %d\n",iret);
            exit(EXIT_FAILURE);
        }
    }

    // join threads
    for (i = 0; i < NUMTHREADS; ++i) {
        pthread_join( thread[i], NULL);
    }

    pthread_barrier_destroy(&mybarrier);

    //show_matrix(mat, n);
    //show_matrix(matOriginal, n);

    //cholesky(mat, n, rank, nProcesses);

    // for(i = 0; i < n; i++) {
	// 	for(j = 0; j < n; j++) {
	// 		printf("%lf ", mat[i * n + j]);
	// 	}
	// 	printf("\n");
	// }

	printf("dupa cholesky\n");

    verifyCholesky(matOriginal, mat, n);


    free(mat);
    free(matOriginal);


    return 0;
}
