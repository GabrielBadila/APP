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

void cholesky(int rank) {
	int i, j, k;

	for(j = 0; j < n; j++) {
        /* Replace the entries above the diagonal with zeroes */
		if (rank == 0) {
			for (i = 0; i < j; i++) {
				mat[i * n + j] = 0;
			}
		}

        /* Update the diagonal element */
		if (j % NUMTHREADS == rank) {
			for (k = 0; k < j; k++) {
                mat[j * n + j] = mat[j * n + j] - mat[j * n + k] * mat[j * n + k];
			}

            mat[j * n + j] = sqrt(mat[j * n + j]);
		}

        /* Wait for all threads to update the diagonal */
        pthread_barrier_wait(&mybarrier);

		/* Divide the rest of the work and update the elements below the diagonal */
		for (i = j + 1; i < n; i++) {
			if (i % NUMTHREADS == rank) {
				for (k = 0; k < j; k++) {
                    mat[i * n + j] = mat[i * n + j] - mat[i * n + k] * mat[j * n + k];
				}

                mat[i * n + k] = mat[i * n + k] / mat[j * n + j];
			}
		}
	}
}

void *choleskyThread(void *arguments) {
    struct arg_struct *args = arguments;
    cholesky(args->rank);
}

void show_matrix(double *mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", mat[i * n + j]);
        printf("\n");
    }
}

void verifyCholesky(double *mat, double *L, int n) {
    int i, j, k;
    double sum;
    double *rez = (double*) malloc(n * n * sizeof(double));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum = 0;
            for(k = 0; k < n; k++) {
                sum += L[i * n + k] * L[j * n + k];
            }
            rez[i*n+j] = sum;
        }
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(round(mat[i*n+j]) != round(rez[i*n+j])) {
                printf("Wrong matrix!\n");
                printf("Position: (%d, %d)\n", i, j);
                return;
            }
        }
    }
    printf("Correct matrix!\n");
}

int main(int argc, char *argv[]) {
    int i, j, iret;
    FILE *f = NULL;

    pthread_t thread[NUMTHREADS];
    struct arg_struct args_s[NUMTHREADS];

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <in_file>\n", argv[0]);
        exit(1);
    }

    f = fopen(argv[1], "r");
    fscanf(f, "%d", &n);
    printf("Size: %dx%d\n", n, n);


    mat = (double*) malloc(n * n * sizeof(double));
    matOriginal = (double*) malloc(n * n * sizeof(double));

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            fscanf(f, "%lf", &mat[i * n + j]);
        }
    }
    fclose(f);

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            matOriginal[i * n + j] = mat[i * n + j];
        }
    }

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
    verifyCholesky(matOriginal, mat, n);

    free(mat);
    free(matOriginal);

    return 0;
}
