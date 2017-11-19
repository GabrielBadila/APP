#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void cholesky(double *a, int n) {
    int i, j, k;
    double *L = (double*)calloc(n * n, sizeof(double));
    if (a == NULL)
        exit(EXIT_FAILURE);

    for(k = 0; k < n; k++) {
        a[k * n + k] = sqrt(a[k * n + k]);

        #pragma omp parallel for
        for(i = k + 1; i < n; i++) {
            a[i * n + k] = a[i * n + k] / a[k * n + k];
        }

        #pragma omp parallel for collapse(2)
        for(j = k + 1; j < n; j++) {
            for(i = 0; i < n; i++) {
                if (i >= j) {
                    a[i * n + j] = a[i * n + j] - a[i * n + k] * a[j * n + k];
                }
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j > i) {
                a[i * n + j] = 0;
            }
        }
    }
}

void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}

void verifyCholesky(double *mat, double *a, int n) {
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

    //show_matrix(rez, n);
    //printf("\n");
    //show_matrix(mat, n);
    //printf("\n");


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

int main() {
    int n, i, j;
    FILE *f = fopen("testFile5.txt", "r");

    fscanf(f, "%d", &n);

    printf("%d\n", n);

    double mat[n * n];

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            fscanf(f, "%lf", &mat[i * n + j]);
        }
    }

    double matOriginal[n * n];

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            matOriginal[i * n + j] = mat[i * n + j];
        }
    }

    cholesky(mat, n);
    // show_matrix(mat, n);
    // printf("\n");

    verifyCholesky(matOriginal, mat, n);

    fclose(f);

    return 0;
}
