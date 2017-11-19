#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cholesky(double *a, int n) {
    int i, j, k;
    double *L = (double*)calloc(n * n, sizeof(double));
    if (a == NULL)
        exit(EXIT_FAILURE);

    for(k = 0; k < n; k++) {
        a[k * n + k] = sqrt(a[k * n + k]);

        for(i = k + 1; i < n; i++) {
            a[i * n + k] = a[i * n + k] / a[k * n + k];
        }

        for(j = k + 1; j < n; j++) {
            for(i = j; i < n; i++) {
                a[i * n + j] = a[i * n + j] - a[i * n + k] * a[j * n + k];
            }
        }
    }

    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            a[i*n + j] = 0;
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

int main() {
    int n = 3;
    double m1[] = {25, 15, -5,
                   15, 18,  0,
                   -5,  0, 11};
    cholesky(m1, n);
    show_matrix(m1, n);
    printf("\n");

    n = 4;
    double m2[] = {18, 22,  54,  42,
                   22, 70,  86,  62,
                   54, 86, 174, 134,
                   42, 62, 134, 106};
    cholesky(m2, n);
    show_matrix(m2, n);

    return 0;
}
