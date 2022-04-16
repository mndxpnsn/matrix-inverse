/*
 * main.cpp
 *
 *  Created on: Apr 16, 2022
 *      Author: d-w-h
 */

#include <stdio.h>
#include <math.h>

double ** mat2D(int n) {

    double ** mat = new double * [n];

    for(int i = 0; i < n; ++i)
        mat[i] = new double[n];

    return mat;
}

void free_mat2D(double ** mat, int n) {

    for(int i = 0; i < n; ++i)
        delete [] mat[i];

    delete [] mat;
}


double determinant(double ** A, int n) {
    double det = 0;

    if(n == 1) {
        return A[0][0];
    }

    if(n == 2) {
        return A[0][0] * A[1][1] - A[1][0] * A[0][1];
    }

    if(n > 2) {
        for(int c = 0; c < n; ++c) {

            double ** M = mat2D(n - 1);

            for(int i = 1; i < n; ++i) {
                int j_m = 0;
                for(int j = 0; j < n; ++j) {
                    if(j != c) {
                        M[i - 1][j_m] = A[i][j];
                        j_m++;
                    }

                }
            }

            double fac = pow(-1, c + 2);

            det = det + A[0][c] * fac * determinant(M, n - 1);

            free_mat2D(M, n - 1);

        }
    }

    return det;
}


double co_factor(double ** A, int n, int i, int j) {

    double fac = 0;

    double ** M = mat2D(n - 1);

    int i_m = 0;
    for(int r = 0; r < n; ++r) {
        int j_m = 0;
        if(r != i) {
            for(int c = 0; c < n; ++c) {
                if(c != j) {
                    M[i_m][j_m] = A[r][c];
                    j_m++;
                }
            }
            i_m++;
        }
    }

    fac = pow(-1, i + j + 2) * determinant(M, n - 1);

    free_mat2D(M, n - 1);

    return fac;
}

void adj(double ** A, int n, double ** adj_mat) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            adj_mat[i][j] = co_factor(A, n, i, j);
        }
    }
}

void mat_inverse(double ** A, int n, double ** mat_inv) {

    double ** adj_mat = mat2D(n);

    adj(A, n, adj_mat);

    double det = determinant(A, n);

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            mat_inv[j][i] = 1.0 / det * adj_mat[i][j];
        }
    }

    free_mat2D(adj_mat, n);
}

void mat_mult_sq(double ** A, double ** A_inv, int n, double ** mat_res) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0;

            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + A[i][k] * A_inv[k][j];
            }

            mat_res[i][j] = sum_loc;
        }
    }
}

void print_mat(double ** mat, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            printf("%.2f ", mat[i][j]);
        }
        printf("\n");
    }

    printf("\n");
}

int main(int argc, char * argv[]) {

    // Size of matrix
    int n = 4;

    // Allocate space for matrices
    double ** A = mat2D(n);
    double ** A_inv = mat2D(n);
    double ** mat_prod = mat2D(n);

    // Populate matrix A with some data
    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = -3;
    A[0][3] = 4;

    A[1][0] = -4;
    A[1][1] = 2;
    A[1][2] = 1;
    A[1][3] = 3;

    A[2][0] = 3;
    A[2][1] = 0;
    A[2][2] = 0;
    A[2][3] = -3;

    A[3][0] = 2;
    A[3][1] = 0;
    A[3][2] = -2;
    A[3][3] = 3;

    // Compute the inverse of A
    mat_inverse(A, n, A_inv);

    // Compute the product A and its inverse A_inv
    mat_mult_sq(A, A_inv, n,  mat_prod);

    // Print A
    print_mat(A, n);

    // Print the inverse of A
    print_mat(A_inv, n);

    // Print product of A with its inverse
    print_mat(mat_prod, n);

    // Free allocated space
    free_mat2D(A, n);
    free_mat2D(A_inv, n);
    free_mat2D(mat_prod, n);

    return 0;
}



