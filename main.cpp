//
//  main.cpp
//  mat-inv
//
//  Created by mndx on 16/04/2022.
//

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


double determinant(double ** mat, int n) {
    double det = 0;

    if(n == 1) {
        return mat[0][0];
    }

    if(n == 2) {
        return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    }

    if(n > 2) {
        for(int c = 0; c < n; ++c) {

            double ** mat_red = mat2D(n - 1);

            for(int i = 1; i < n; ++i) {
                int j_m = 0;
                for(int j = 0; j < n; ++j) {
                    if(j != c) {
                        mat_red[i - 1][j_m] = mat[i][j];
                        j_m++;
                    }

                }
            }

            double fac = pow(-1, c + 2);

            det = det + mat[0][c] * fac * determinant(mat_red, n - 1);

            free_mat2D(mat_red, n - 1);

        }
    }

    return det;
}


double co_factor(double ** mat, int n, int i, int j) {

    double fac = 0;

    double ** mat_red = mat2D(n - 1);

    int i_m = 0;
    for(int r = 0; r < n; ++r) {
        int j_m = 0;
        if(r != i) {
            for(int c = 0; c < n; ++c) {
                if(c != j) {
                    mat_red[i_m][j_m] = mat[r][c];
                    j_m++;
                }
            }
            i_m++;
        }
    }

    fac = pow(-1, i + j + 2) * determinant(mat_red, n - 1);

    free_mat2D(mat_red, n - 1);

    return fac;
}

void adjoint_mat(double ** mat, int n, double ** adj_mat) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            adj_mat[i][j] = co_factor(mat, n, i, j);
        }
    }
}

void mat_inverse(double ** mat, int n, double ** mat_inv) {

    double ** adj_mat = mat2D(n);

    adjoint_mat(mat, n, adj_mat);

    double det = determinant(mat, n);

    // Check if matrix is singular
    if(det == 0)
        printf("Matrix is singular\n");
    
    // Matrix is nonsingular
    if(det != 0) {
        
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                mat_inv[j][i] = 1.0 / det * adj_mat[i][j];
            }
        }
        
    }

    free_mat2D(adj_mat, n);
}

void mat_mult_sq(double ** mat, double ** mat_inv, int n, double ** mat_res) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0;

            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + mat[i][k] * mat_inv[k][j];
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

double rand_num(double min, double max) {
    
    double val = (double) rand() / (RAND_MAX + 1.0);
    
    return val * (max - min) - (max - min) / 2.0;
}

void init_mat(int n, double ** mat) {
    
    double max_range = 50;
    double min_range = -50;
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            mat[i][j] = rand_num(min_range, max_range);
        }
    }
}

int main(int argc, char * argv[]) {

    // Size of matrix
    int n = 4;

    // Allocate space for matrices
    double ** mat = mat2D(n);
    double ** mat_inv = mat2D(n);
    double ** mat_prod = mat2D(n);

    // Populate matrix mat with some data
    init_mat(n, mat);

    // Compute the inverse of mat
    mat_inverse(mat, n, mat_inv);

    // Compute the product of mat and its inverse mat_inv
    mat_mult_sq(mat, mat_inv, n,  mat_prod);

    // Print mat
    print_mat(mat, n);

    // Print the inverse of mat
    print_mat(mat_inv, n);

    // Print product of mat with its inverse
    print_mat(mat_prod, n);

    // Free allocated space
    free_mat2D(mat, n);
    free_mat2D(mat_inv, n);
    free_mat2D(mat_prod, n);

    return 0;
}



