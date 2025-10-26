// compile: mpicc -o testpetsc testpetsc.c -I${PETSC_DIR}/include -L${PETSC_DIR}/lib -Wl,-rpath=${PETSC_DIR}/lib -lpetsc
// run: mpirun -np 4 ./testpetsc -ksp_type gmres -pc_type hypre -ksp_rtol 1e-6

#include <petscksp.h>


int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    /* Problem parameters */
    PetscInt    N = 10;          /* System size */
    PetscScalar v;               /* Matrix value */
    PetscInt    Ii, J, Istart, Iend;
    Vec         x, b;            /* Solution and RHS vectors */
    Mat         A;               /* Sparse matrix */
    KSP         ksp;             /* Linear solver context */
    PC          pc;              /* Preconditioner context */
    PetscReal   norm;            /* Norm of solution */
    PetscInt    its;             /* Iteration count */

    /* Create vectors */
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    /* Create matrix */
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N);
    MatSetFromOptions(A);
    MatSetUp(A);

    /* Assemble matrix */
    MatGetOwnershipRange(A, &Istart, &Iend);
    for (Ii = Istart; Ii < Iend; Ii++) {
        v = -1.0; J = Ii - 1;
        if (J >= 0) MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES);

        v = 2.0; J = Ii;
        MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES);

        v = -1.0; J = Ii + 1;
        if (J < N) MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES);

    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    /* Set RHS vector */
    VecSet(b, 1.0);  /* All entries = 1.0 */

    /* Create solver context */
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);

    /* Set solver options at runtime */
    KSPSetFromOptions(ksp);

    /* Solve the system */
    KSPSolve(ksp, b, x);

    /* View solver info */
    KSPGetIterationNumber(ksp, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Iterations %D\n", its);

    /* Check solution */
    VecNorm(x, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD, "Solution norm: %g\n", (double)norm);

    /* Clean up */
    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

    PetscFinalize();
    return 0;

}
