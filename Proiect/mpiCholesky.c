#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void cholesky(double *a, int n, int rank, int nProcesses) {
	int i, j, k;
	double *copy = (double*) malloc(n * sizeof(double));

    for(j = 0; j < n; j++) {
        /* Replace the entries above the diagonal with zeroes */
        if (rank == 0) {
			for (i = 0; i < j; i++) {
				a[i * n + j] = 0;
            }
		}

        /* Update the diagonal element */
        if (j % nProcesses == rank) {
            for (k = 0; k < j; k++) {
                a[j * n + j] = a[j * n + j] - a[j * n + k] * a[j * n + k];
            }

            a[j * n + j] = sqrt(a[j * n + j]);
        }

        /* Broadcast row with new values to other processes */
		for(i = 0; i < n; i++) {
            copy[i] = a[j * n + i];
		}

		MPI_Bcast(copy, n, MPI_DOUBLE, j % nProcesses, MPI_COMM_WORLD);

        for(i = 0; i < n; i++) {
			a[j * n + i] = copy[i];
		}

        /* Divide the rest of the work and update the elements below the diagonal */
        for (i = j + 1; i < n; i++) {
			if (i % nProcesses == rank) {
				for (k = 0; k < j; k++) {
                    a[i * n + j] = a[i * n + j] - a[i * n + k] * a[j * n + k];
				}

                a[i * n + k] = a[i * n + k] / a[j * n + j];
			}
		}
	}
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
    int n, i, j, rank, nProcesses;
    double *mat, *matOriginal;
	FILE * f = NULL;

	if (argc < 2) {
        fprintf(stderr, "Usage: mpirun --oversubscribe -np <num_procs> %s <in_file>\n", argv[0]);
        exit(1);
    }

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
	printf("Hello from %i/%i\n", rank, nProcesses);

    if(rank == 0) {
		f = fopen(argv[1], "r");

	    fscanf(f, "%d", &n);
	    printf("Size: %dx%d\n", n, n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mat = (double*) malloc(n * n * sizeof(double));
    matOriginal = (double*) malloc(n * n * sizeof(double));


    if (rank == 0) {
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

    }

    MPI_Bcast(mat, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    cholesky(mat, n, rank, nProcesses);

    if (rank == 0) {
        verifyCholesky(matOriginal, mat, n);
    }

    free(mat);
    free(matOriginal);

    MPI_Finalize();
    return 0;
}
