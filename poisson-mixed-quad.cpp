// ================================================
// solve 
// \vec{u} + \nabla p = 0
// \nabla \cdot \vec{u} = f
// using RT0/DG0 BDFM1/DG0 on quad meshes.
// ================================================
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <petscksp.h>
#include <stdio.h>

#define N 2
#define h (1.0/N)
#define area (h*h)
#define First 0
#define Second 1
using namespace std;
using RealMat = Eigen::MatrixXd;
using RealVec = Eigen::VectorXd;
typedef double Coord[2];

class Quad1d {
    public:
        double points[2] = {(1-1/sqrt(3))/2, (1+1/sqrt(3))/2};
        double weights[2] = {0.5, 0.5};
};

class Quad2d {
    public:
        Quad1d quad1d;
        double points[4][2] = {
            {quad1d.points[0], quad1d.points[0]},
            {quad1d.points[1], quad1d.points[0]},
            {quad1d.points[0], quad1d.points[1]},
            {quad1d.points[1], quad1d.points[1]}
        };
        double weights[4] = {0.25, 0.25, 0.25, 0.25};
};

class Elem {
    public:
        Elem(int i, int j) : x_left(i*h), x_right((i+1)*h), y_down(j*h), y_up((j+1)*h) {
            x_center = (x_left + x_right) / 2;
            y_center = (y_down + y_up) / 2;
        }
        int i, j;
        double x_left;
        double x_right;
        double x_center;
        double y_down;
        double y_up;
        double y_center;
};

int map_u(int i, int j, int ibas) {
    int idof_global;
    switch (ibas) {
        case 0:
            idof_global = N*(N+1) + j*N + i;
            break;
        case 1:
            idof_global = j*(N+1) + (i+1);
            break;
        case 2:
            idof_global = N*(N+1) + (j+1)*N + i;
            break;
        case 3:
            idof_global = j*(N+1) + i;
            break;
    }
    return idof_global;
}

int map_p(int i, int j) {
    return (2*N*(N+1) + j*N + i);
}

class FERT0 {
    public:
        FERT0();
        void CalcShapeRef(const double *ip, RealMat &shape) const;
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(const double *ip, RealVec &divshape) const;
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Coord *nodes = nullptr; 
        Eigen::PartialPivLU<RealMat> Ti;
        const int nbas = 4;
        const double nk[8] =
        { 0., -1.,
          1.,  0.,
          0.,  1.,
         -1.,  0. };
        vector<int> dof2nk;
};

FERT0::FERT0() {
    nodes = (Coord *) calloc (nbas, sizeof(Coord));
    dof2nk.resize(nbas);
    int o = 0;
    nodes[o][0] = 1./2.; nodes[o][1] = 0.0;   dof2nk[o++] = 0;
    nodes[o][0] = 1.;    nodes[o][1] = 1./2.; dof2nk[o++] = 1;
    nodes[o][0] = 1./2.; nodes[o][1] = 1.0;   dof2nk[o++] = 2;
    nodes[o][0] = 0.;    nodes[o][1] = 1./2.; dof2nk[o++] = 3;
    RealMat T(nbas, nbas);
    for (int k = 0; k < nbas; k++) {
        const double *ip = nodes[k];
        const double *n_k = nk + 2*dof2nk[k];

        o = 0;
        T(o++,k) = n_k[0];
        T(o++,k) = n_k[1];
        T(o++,k) = ip[0] * n_k[0];
        T(o++,k) = ip[1] * n_k[1];
    }
    Ti = Eigen::PartialPivLU<RealMat> (T);
}

void FERT0::CalcShapeRef(const double *ip, RealMat &shape) const {
    RealMat u(nbas,2);
    int o = 0;
    u(o,0) = 1; u(o,1) = 0;
    o++;
    u(o,0) = 0; u(o,1) = 1;
    o++;
    u(o,0) = ip[0]; u(o,1) = 0;
    o++;
    u(o,0) = 0; u(o,1) = ip[1];
    o++;
    assert(o == nbas);

    shape = Ti.solve(u);
}

void FERT0::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(ip, shape);
}

void FERT0::CalcDivShapeRef(const double *ip, RealVec &divshape) const {
    int o = 0; 
    RealVec divu(nbas);
    divu(o++) = 0; 
    divu(o++) = 0; 
    divu(o++) = 1; 
    divu(o++) = 1; 
    assert(o == nbas);

    divshape = Ti.solve(divu);
}

void FERT0::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(ip, divshape);
}

class FEDG0 {
    public:
        FEDG0();
        void CalcShapeRef(const double *ip, RealVec &shape) const;
        void CalcShape(const Elem &e, const double *ip, RealVec &shape);
        void CalcGradShapeRef(const double *ip, RealMat &gradshape) const;
        void CalcGradShape(const Elem &e, const double *ip, RealMat &gradshape);
        Coord *nodes = nullptr;
};

FEDG0::FEDG0() {
    nodes = (Coord *) calloc (1, sizeof(Coord));
    nodes[First][0] = 1./2.;
    nodes[First][1] = 1./2.;
}

void FEDG0::CalcShapeRef(const double *ip, RealVec &shape) const {
    assert(shape.size() == 1);
    shape(First) = 1.;
}

void FEDG0::CalcShape(const Elem &e, const double *ip, RealVec &shape) {
    assert(shape.size() == 1);
    shape(First) = 1.;
}

void FEDG0::CalcGradShapeRef(const double *ip, RealMat &gradshape) const {
    assert(gradshape.rows() == 1);
    gradshape(First, 0) = 0.;
    gradshape(First, 1) = 0.;
}

void FEDG0::CalcGradShape(const Elem &e, const double *ip, RealMat &gradshape) {
    assert(gradshape.rows() == 1);
    gradshape(First, 0) = 0.;
    gradshape(First, 1) = 0.;
}



// --------------------------------- Solver using Eigen ---------------------------------
class Solver {
    public:
        Solver() {
            u = new FERT0();
            p = new FEDG0();
        }
        int total_dofs = 2*N*(N+1) + N*N;

        void init();
        void build_mat();
        void build_rhs();
        void solve();
        void step();
        FERT0 *u;
        FEDG0 *p;
        Quad2d quad2d;
        Quad1d quad1d;
        RealMat mat;
        RealVec rhs;
        RealVec x;
};


// --------------------------------- Solver using Petsc ---------------------------------
/* class Solver { */
/*     public: */
/*         Solver() { */
/*             u = new FERT0(); */
/*             p = new FEDG0(); */
/*         } */
/*         ~Solver() { */
/*             delete u; */
/*             delete p; */
/*             MatDestroy(&mat); */
/*             VecDestroy(&rhs); */
/*             VecDestroy(&x); */
/*             KSPDestroy(&ksp); */
/*         } */
/*         void init(); */
/*         void build_mat(); */
/*         void build_rhs(); */
/*         void solve(); */
/*         void step(); */
/*         Mat mat = nullptr; */
/*         Vec rhs = nullptr; */
/*         KSP ksp = nullptr; */
/*         Vec x = nullptr; */
/*         Vec x0 = nullptr; */
/*         FERT0 *u; */
/*         FEDG0 *p; */
/*         Quad2d quad2d; */
/* }; */

void Solver::init() {
    cout << "Solver Init\n";
    cout << "u_space: RT0\n";
    cout << "p_space: DG0\n";
    
    // Initialize PETSc objects
    mat.resize(total_dofs, total_dofs);
    mat.setZero();
    rhs.resize(total_dofs);
    rhs.setZero();
    x.resize(total_dofs);
    x.setZero();
}


void Solver::build_mat() {
    cout << "Build Mat\n";
    int Nu = 4;
    int Np = 1;
    int Ndof = Nu + Np;
    int Id[Ndof];
    int *Iu = Id, *Ip = Iu + Nu;

    RealMat A(Ndof,Ndof);
    RealMat u_shape(Nu,2), dp_shape(Np,2);
    RealVec p_shape(Np), divu_shape(Nu);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Elem e(i,j);
            for (int k = 0; k < Nu; k++) {
                Iu[k] = map_u(i,j,k);
            }
            Ip[0] = map_p(i,j);

            A.setZero();
            for (int k = 0; k < 4; k++) {
                u->CalcShape(e, quad2d.points[k], u_shape);
                u->CalcDivShape(e, quad2d.points[k], divu_shape);
                p->CalcShape(e, quad2d.points[k], p_shape);
                p->CalcGradShape(e, quad2d.points[k], dp_shape);
                
                for (int m = 0; m < Nu; m++) {
                    // A(u,v)
                    for (int n = 0; n < Nu; n++) {
                        A(m,n) += quad2d.weights[k] * area * (u_shape(m,0)*u_shape(n,0) + u_shape(m,1)*u_shape(n,1));
                    }
                    // A(p,v)
                    for (int n = 0; n < Np; n++) {
                        A(m,n+Nu) += quad2d.weights[k] * area * (-divu_shape(m) * p_shape(n));
                    }
                }
                
                // A(q,u)
                for (int m = 0; m < Np; m++) {
                    for (int n = 0; n < Nu; n++) {
                        A(m + Nu, n) += quad2d.weights[k] * area * (p_shape(m) * divu_shape(n));
                    }
                }
            }

            // Apply boundary conditions (Dirichlet)
            if (i == 0) {
                A(3,3) = 1; 
            }
            if (i == N-1) {
                A(1,1) = 1;
            }
            if (j == 0) {
                A(0,0) = 1;
            }
            if (j == N-1) {
                A(2,2) = 1;
            }

            // Add element matrix to global matrix
            for (int m = 0; m < Ndof; m++) {
                for (int n = 0; n < Ndof; n++) {
                    mat(Id[m], Id[n]) += A(m,n);
                }
            }
        }
    }

    /* MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); */
    /* MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY); */
    
    /* PetscViewer viewer; */
    /* PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A.m", &viewer); */
    /* PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); */
    /* PetscObjectSetName((PetscObject)mat, "A"); */
    /* MatView(mat, viewer); */
    /* PetscViewerPopFormat(viewer); */
    /* PetscViewerDestroy(&viewer); */

    cout << "Build Mat Done\n";

    /* for (int i = 0; i < total_dofs; i++) { */
    /*     for (int j = 0; j < total_dofs; j++) */ 
    /*         printf("%4.2f\t", mat(i,j)); */
    /*         /1* cout << "%4.2f" << mat(i,j) << "\t"; *1/ */
    /*     cout << endl; */
    /* } */
}

void Solver::build_rhs() {
    cout << "Build Rhs\n";
    /* VecSet(rhs, 0.0); */
    
    int Nu = 4;
    int Np = 1;
    int Ndof = Nu + Np;
    int Id[Ndof];
    int *Iu = Id, *Ip = Iu + Nu;
    
    RealVec F(Np);
    RealVec p_shape(Np), divu_shape(Nu);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Elem e(i,j);
            for (int k = 0; k < Nu; k++) {
                Iu[k] = map_u(i,j,k);
            }
            Ip[0] = map_p(i,j);
            
            F.setZero();
            for (int k = 0; k < 4; k++) {
                p->CalcShape(e, quad2d.points[k], p_shape);
                /* u->CalcDivShape(e, quad2d.points[k], divu_shape); */
                
                // Define source term f (example: f = 1)
                double f = 1.0;
                
                // F(q)
                F(0) += quad2d.weights[k] * area * (p_shape(0) * f);
            }
            
            /* // Apply boundary conditions (Dirichlet) */
            /* if (i == 0) { */
            /*     F(3) = 0.0; // Value for Dirichlet BC */
            /* } */
            /* if (i == N-1) { */
            /*     F(1) = 0.0; */
            /* } */
            /* if (j == 0) { */
            /*     F(0) = 0.0; */
            /* } */
            /* if (j == N-1) { */
            /*     F(2) = 0.0; */
            /* } */
            
            /* VecSetValues(rhs, Ndof, Id, &F[0], ADD_VALUES); */
        }
    }
    
    /* VecAssemblyBegin(rhs); */
    /* VecAssemblyEnd(rhs); */
}

/* void Solver::solve() { */
/*     cout << "Solve Use Petsc\n"; */
    
/*     KSPSetOperators(ksp, mat, mat); */
/*     KSPSetFromOptions(ksp); */
    
/*     KSPSolve(ksp, rhs, x); */
    
/*     // View solution */
/*     PetscViewer viewer; */
/*     PetscViewerASCIIOpen(PETSC_COMM_WORLD, "solution.m", &viewer); */
/*     PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); */
/*     PetscObjectSetName((PetscObject)x, "x"); */
/*     VecView(x, viewer); */
/*     PetscViewerPopFormat(viewer); */
/*     PetscViewerDestroy(&viewer); */
/* } */

void Solver::step() {
    build_mat();
    build_rhs();
    /* solve(); */
}

int main(int argc, char* argv[]) {
    /* PetscInitialize(&argc, &argv, NULL, NULL); */
    
    Solver solver;
    solver.init();
    solver.build_mat();
    /* solver.step(); */
    
    /* PetscFinalize(); */
    return 0;
}
