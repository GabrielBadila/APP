#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *cholesky(double *a, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (a == NULL)
        exit(EXIT_FAILURE);

    for(int k=0;k<n;k++) {
      a[k*n+k] = sqrt(a[k*n+k]);

      for(int i=k+1;i<n;i++) {
        a[i*n+k] = a[i*n+k]/a[k*n+k];
      }

      for(int j=k+1;j<n;j++) {
        for(int i=j;i<n;i++) {
          a[i*n+j] = a[i*n+j] - a[i*n+k]*a[j*n+k];
        }
      }
    }

    return L;
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
    double *c1 = cholesky(m1, n);
    show_matrix(m1, n);
    printf("\n");
    free(c1);

    n = 4;
    double m2[] = {18, 22,  54,  42,
                   22, 70,  86,  62,
                   54, 86, 174, 134,
                   42, 62, 134, 106};
    double *c2 = cholesky(m2, n);
    show_matrix(m2, n);
    free(c2);

    return 0;
}
