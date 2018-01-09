#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* The algorithm for Cholesky decomposition */
void cholesky(double *mat, int n) {
    int i, j, k;

    if (mat == NULL)
        exit(EXIT_FAILURE);

    // main Cholesky algorithm
    for(k = 0; k < n; k++) {
        mat[k * n + k] = sqrt(mat[k * n + k]);

        for(i = k + 1; i < n; i++) {
            mat[i * n + k] = mat[i * n + k] / mat[k * n + k];
        }

        for(j = k + 1; j < n; j++) {
            for(i = j; i < n; i++) {
                mat[i * n + j] = mat[i * n + j] - mat[i * n + k] * mat[j * n + k];
            }
        }
    }

    // zero all elements above the diagonal
    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            mat[i*n + j] = 0;
        }
    }
}

/* Function to print a matrix */
void show_matrix(double *mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", mat[i * n + j]);
        printf("\n");
    }
}

/* Function to verify the accuracy of the lower triangular matrix obtained */
void verifyCholesky(double *mat, double *L, int n) {
    int i, j, k;
    double sum;
    double *rez = (double*) malloc(n * n * sizeof(double));

    // multiply L*L'
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum = 0;
            for(k = 0; k < n; k++) {
                sum += L[i * n + k] * L[j * n + k];
            }
            rez[i*n+j] = sum;
        }
    }

    // compare rez with original matrix
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(round(mat[i*n+j]) != round(rez[i*n+j])) {
                printf("Wrong matrix!\n");
                printf("Position: (%d, %d)\n", i, j);
                return;
            }
        }
    }
    free(rez);
    printf("Correct matrix!\n");
}

/* Driver program to test above functions */
int main(int argc, char *argv[]) {
    int n, i, j;
    FILE * f = NULL;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <in_file>\n", argv[0]);
        exit(1);
    }

    f = fopen(argv[1], "r");

    fscanf(f, "%d", &n);
    printf("Size: %dx%d\n", n, n);

    double *mat = (double*) malloc(n * n * sizeof(double));
    double *matOriginal = (double*) malloc(n * n * sizeof(double));

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            fscanf(f, "%lf", &mat[i * n + j]);
        }
    }

    // save original matrix into a copy because it will be modified
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            matOriginal[i * n + j] = mat[i * n + j];
        }
    }

    cholesky(mat, n);
    verifyCholesky(matOriginal, mat, n);

    fclose(f);
    free(mat);
    free(matOriginal);

    return 0;
}
