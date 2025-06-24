// laplace_fem.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5  // Number of nodes per side (change for finer mesh)
#define NN (N*N) // Total number of nodes
#define NE (2*(N-1)*(N-1)) // Number of elements (2 triangles per square)
#define IDX(i,j) ((i)*N+(j))

// Dirichlet boundary condition function
double boundary(double x, double y) {
    // Example: u = x + y on the boundary
    return x + y;
}

// Mesh node coordinates
void generate_nodes(double *x, double *y) {
    int i, j, k=0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            x[k] = (double)j/(N-1);
            y[k] = (double)i/(N-1);
            k++;
        }
    }
}

// Element connectivity (triangles)
void generate_elements(int elem[NE][3]) {
    int e = 0;
    for (int i=0; i<N-1; i++) {
        for (int j=0; j<N-1; j++) {
            int n0 = IDX(i,j);
            int n1 = IDX(i,j+1);
            int n2 = IDX(i+1,j);
            int n3 = IDX(i+1,j+1);
            // Lower triangle
            elem[e][0] = n0; elem[e][1] = n1; elem[e][2] = n2; e++;
            // Upper triangle
            elem[e][0] = n1; elem[e][1] = n3; elem[e][2] = n2; e++;
        }
    }
}

// Assemble global stiffness matrix and load vector
void assemble(double *x, double *y, int elem[NE][3], double K[NN][NN], double F[NN]) {
    for (int e=0; e<NE; e++) {
        int *n = elem[e];
        double x0=x[n[0]], y0=y[n[0]];
        double x1=x[n[1]], y1=y[n[1]];
        double x2=x[n[2]], y2=y[n[2]];
        // Area of triangle
        double area = fabs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))/2.0;
        // Gradients of shape functions
        double b[3], c[3];
        b[0] = y1 - y2; b[1] = y2 - y0; b[2] = y0 - y1;
        c[0] = x2 - x1; c[1] = x0 - x2; c[2] = x1 - x0;
        // Local stiffness matrix
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                double kij = (b[i]*b[j] + c[i]*c[j])/(4.0*area)*area;
                K[n[i]][n[j]] += kij;
            }
        }
        // No source term, so F unchanged
    }
}

// Apply Dirichlet boundary conditions
void apply_bc(double *x, double *y, double K[NN][NN], double F[NN], double U[NN]) {
    for (int i=0; i<NN; i++) {
        if (x[i]==0 || x[i]==1 || y[i]==0 || y[i]==1) {
            U[i] = boundary(x[i], y[i]);
            for (int j=0; j<NN; j++) {
                K[i][j] = 0.0;
                K[j][i] = 0.0;
            }
            K[i][i] = 1.0;
            F[i] = U[i];
        }
    }
}

// Simple Gaussian elimination solver
void solve(double K[NN][NN], double F[NN], double U[NN]) {
    int i, j, k;
    double factor;
    // Forward elimination
    for (k=0; k<NN-1; k++) {
        for (i=k+1; i<NN; i++) {
            factor = K[i][k]/K[k][k];
            for (j=k; j<NN; j++)
                K[i][j] -= factor*K[k][j];
            F[i] -= factor*F[k];
        }
    }
    // Back substitution
    for (i=NN-1; i>=0; i--) {
        U[i] = F[i];
        for (j=i+1; j<NN; j++)
            U[i] -= K[i][j]*U[j];
        U[i] /= K[i][i];
    }
}

int main() {
    double x[NN], y[NN], K[NN][NN]={0}, F[NN]={0}, U[NN]={0};
    int elem[NE][3];
    generate_nodes(x, y);
    generate_elements(elem);
    assemble(x, y, elem, K, F);
    apply_bc(x, y, K, F, U);
    solve(K, F, U);

    printf("Node (x, y): Solution u\n");
    for (int i=0; i<NN; i++)
        printf("(%.2f, %.2f): %.6f\n", x[i], y[i], U[i]);
    return 0;
}
