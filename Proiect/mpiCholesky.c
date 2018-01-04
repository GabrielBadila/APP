#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void cholesky(double *a, int n, int rank, int nProcesses) {
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
				a[i * n + j] = 0;
			}
		}

		/*
		 * Step 1:
		 * Update the diagonal element
		 */

		if (j % nProcesses == rank) {

			for (k = 0; k < j; k++) {
				//L[j][j] = L[j][j] - L[j][k] * L[j][k];
                a[j * n + j] = a[j * n + j] - a[j * n + k] * a[j * n + k];
			}

			//L[j][j] = sqrt(L[j][j]);
            a[j * n + j] = sqrt(a[j * n + j]);
		}

		// Broadcast row with new values to other processes
		//MPI_Bcast(L[j], n, MPI_DOUBLE, j % nProcesses, MPI_COMM_WORLD);

		for(i=0;i<n;i++) {
			copy[i] = a[j * n + i];
		}

		MPI_Bcast(copy, n, MPI_DOUBLE, j % nProcesses, MPI_COMM_WORLD);

        for(i=0;i<n;i++) {
			a[j * n + i] = copy[i];
		}

		/*
		 * Step 2:
		 * Update the elements below the diagonal element
		 */

		// Divide the rest of the work
		for (i = j + 1; i < n; i++) {
			if (i % nProcesses == rank) {
				for (k = 0; k < j; k++) {
					//L[i][j] = L[i][j] - L[i][k] * L[j][k];
                    a[i * n + j] = a[i * n + j] - a[i * n + k] * a[j * n + k];
				}

				//L[i][j] = L[i][j] / L[j][j];
                a[i * n + k] = a[i * n + k] / a[j * n + j];
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
    //MPI_Barrier(MPI_COMM_WORLD);


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

void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int n, rank, nProcesses;
    int i, j, tid;
    FILE *f = NULL;
    double *mat, *matOriginal;

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
	printf("Hello from %i/%i\n", rank, nProcesses);

    if(rank == 0) {
        f = fopen("testFile7.txt", "r");

        fscanf(f, "%d", &n);
        printf("%d\n", n);
    }


    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mat = calloc(n*n, sizeof(double));
    matOriginal = calloc(n*n, sizeof(double));


    if (rank == 0) {
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

    }


    MPI_Bcast(mat, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// for(i = 0; i < n; i++) {
	// 	for(j = 0; j < n; j++) {
	// 		printf("%lf ", mat[i * n + j]);
	// 	}
	// 	printf("\n");
	// }

    cholesky(mat, n, rank, nProcesses);

    // for(i = 0; i < n; i++) {
	// 	for(j = 0; j < n; j++) {
	// 		printf("%lf ", mat[i * n + j]);
	// 	}
	// 	printf("\n");
	// }

	printf("dupa cholesky\n");

    if (rank == 0) {
        verifyCholesky(matOriginal, mat, n);
    }

    free(mat);
    free(matOriginal);



    MPI_Finalize();

    return 0;
}
