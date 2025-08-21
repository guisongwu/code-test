// ================================================
//  solve 
//  \vec{u} + \nabla p = 0
//  \nabla \cdot \vec{u} + beta p = f
//  using RT1/DGQ1 on quad meshes with mass-lumping.
// ================================================
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdio.h>
/* #include "fe-rt.hpp" */
/* #include "fe-dg.hpp" */
/* #include "elem.hpp" */
/* #include "header.hpp" */


#define N 8
#define beta 1
#define h (1.0/N)
#define area (h*h)
#define len h
static int udof_global = 4*N*(2*N+1);

using namespace std;
using RealMat = Eigen::MatrixXd;
using RealVec = Eigen::VectorXd;
typedef double Coord[2];



class Test {
    public:
        double func_p(const Coord &phy_coord) {
            double x = phy_coord[0];
            double y = phy_coord[1];
            return x*x*y*y;
        }
        double func_p(double x, double y) {
            return x*x*y*y;
        }

        void func_u(const Coord &phy_coord, double *u) {
            double x = phy_coord[0];
            double y = phy_coord[1];
            u[0] = -2*x*y*y;
            u[1] = -2*x*x*y;
        }
        void func_u(double x, double y, double *u) {
            u[0] = -2*x*y*y;
            u[1] = -2*x*x*y;
        }

        double func_f(const Coord &phy_coord) {
            double x = phy_coord[0];
            double y = phy_coord[1];
            return -2*(x*x + y*y) + beta * func_p(x,y);
        }
        double func_f(double x, double y) {
            return -2*(x*x + y*y) + beta * func_p(x,y);
        }
};


/* class Test { */
/*     public: */
/*         double func_p(const Coord &phy_coord) { */
/*             double x = phy_coord[0]; */
/*             double y = phy_coord[1]; */
/*             return exp(x)*sin(y); */
/*         } */
/*         double func_p(double x, double y) { */
/*             return exp(x)*sin(y); */
/*         } */

/*         void func_u(const Coord &phy_coord, double *u) { */
/*             double x = phy_coord[0]; */
/*             double y = phy_coord[1]; */
/*             u[0] = -exp(x)*sin(y); */
/*             u[1] = -exp(x)*cos(y); */
/*         } */
/*         void func_u(double x, double y, double *u) { */
/*             u[0] = -exp(x)*sin(y); */
/*             u[1] = -exp(x)*cos(y); */
/*         } */

/*         double func_f(const Coord &phy_coord) { */
/*             return beta * func_p(x,y); */
/*         } */
/*         double func_f(double x, double y) { */
/*             return beta * func_p(x,y); */
/*         } */
/* }; */

Test test;


/* class Quad1d { */
/*     public: */
/*         double points[2] = {(1-1/sqrt(3))/2, (1+1/sqrt(3))/2}; */
/*         double weights[2] = {0.5, 0.5}; */
/* }; */
/* class Quad2d { */
/*     public: */
/*         Quad1d quad1d; */
/*         double points[4][2] = { */
/*             {quad1d.points[0], quad1d.points[0]}, */
/*             {quad1d.points[1], quad1d.points[0]}, */
/*             {quad1d.points[0], quad1d.points[1]}, */
/*             {quad1d.points[1], quad1d.points[1]} */
/*         }; */
/*         double weights[4] = {0.25, 0.25, 0.25, 0.25}; */
/* }; */

class Quad1d {
    public:
        const int npoints = 3;
        double points[3] = {0.112701665379258, 0.5, 0.887298334620742};
        double weights[3] = {0.277777777777778, 0.444444444444444, 0.277777777777778};
};
class Quad2d {
    public:
        Quad1d quad1d;
        const int npoints = 9;
        double points[9][2] = {
            {quad1d.points[0], quad1d.points[0]},
            {quad1d.points[1], quad1d.points[0]},
            {quad1d.points[2], quad1d.points[0]},
            {quad1d.points[0], quad1d.points[1]},
            {quad1d.points[1], quad1d.points[1]},
            {quad1d.points[2], quad1d.points[1]},
            {quad1d.points[0], quad1d.points[2]},
            {quad1d.points[1], quad1d.points[2]},
            {quad1d.points[2], quad1d.points[2]}
        };
        double weights[9] = {
            quad1d.weights[0]*quad1d.weights[0],
            quad1d.weights[1]*quad1d.weights[0],
            quad1d.weights[2]*quad1d.weights[0],
            quad1d.weights[0]*quad1d.weights[1],
            quad1d.weights[1]*quad1d.weights[1],
            quad1d.weights[2]*quad1d.weights[1],
            quad1d.weights[0]*quad1d.weights[2],
            quad1d.weights[1]*quad1d.weights[2],
            quad1d.weights[2]*quad1d.weights[2]
        };
};






class Elem {
    public:
        Elem(int i, int j) : ix(i), jy(j), x_left(i*h), y_down(j*h) {
            if (ix > 0) {
                edge_sign[6] = -1; 
                edge_sign[7] = -1; 
            }
            if (jy > 0) {
                edge_sign[0] = -1; 
                edge_sign[1] = -1;
            }
        }
        
        void ref2phy(const Coord ref_coord, Coord &phy_coord) {
            phy_coord[0] = ref_coord[0]*h + x_left; 
            phy_coord[1] = ref_coord[1]*h + y_down; 
        }

        int ix, jy;
        int edge_sign[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
        double x_left;
        double y_down;
};


int map_u(int i, int j, int ibas) {
    int idof_global;
    switch (ibas) {
        case 0:
        case 1:
            idof_global = i*N*8 + j*8 + ibas;
            break;
        case 2:
            if (i < N-1)
                idof_global = map_u(i+1, j, 7);
            else
                idof_global = N*N*8 + 2*j;
            break;
        case 3:
            if (i < N-1)
                idof_global = map_u(i+1, j, 6);
            else
                idof_global = N*N*8 + 2*j + 1;
            break;
        case 4:
            if (j < N-1)
                idof_global = map_u(i, j+1, 1);
            else
                idof_global = N*N*8 + 2*N + 2*i + 1;
            break;
        case 5:
            if (j < N-1)
                idof_global = map_u(i, j+1, 0);
            else
                idof_global = N*N*8 + 2*N + 2*i;
            break;
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
            idof_global = i*N*8 + j*8 + ibas-4;
            break;
    }
    return idof_global;
}


int map_p(int i, int j, int ibas) {
    return udof_global + 4*i*N + 4*j + ibas;
}


class FERT1 {
    public:
        FERT1();
        void CalcShapeRef(const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(const double *ip, RealVec &divshape);
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Coord *nodes; 
        Eigen::PartialPivLU<RealMat> Ti;
        const int nbas = 12;
        const double nk[8] =
        { 0., -1.,
          1.,  0.,
          0.,  1.,
         -1.,  0. };
        vector<int> dof2nk;
};


// ------------------------------ RT1 ---------------------------------
FERT1::FERT1() {
    nodes = (Coord *) calloc (nbas, sizeof(Coord));
    dof2nk.resize(nbas);
    int o = 0;
    double gauss1 = (1-1/sqrt(3))/2;
    double gauss2 = (1+1/sqrt(3))/2;

    // edge
    nodes[o][0] = gauss1; nodes[o][1] = 0.0;    dof2nk[o++] = 0;
    nodes[o][0] = gauss2; nodes[o][1] = 0.0;    dof2nk[o++] = 0;
    nodes[o][0] = 1.;     nodes[o][1] = gauss1; dof2nk[o++] = 1;
    nodes[o][0] = 1.;     nodes[o][1] = gauss2; dof2nk[o++] = 1;
    nodes[o][0] = gauss2; nodes[o][1] = 1.0;    dof2nk[o++] = 2;
    nodes[o][0] = gauss1; nodes[o][1] = 1.0;    dof2nk[o++] = 2;
    nodes[o][0] = 0.;     nodes[o][1] = gauss2; dof2nk[o++] = 3;
    nodes[o][0] = 0.;     nodes[o][1] = gauss1; dof2nk[o++] = 3;
    // interior
    nodes[o][0] = 0.5;    nodes[o][1] = gauss1; dof2nk[o++] = 1;
    nodes[o][0] = gauss2; nodes[o][1] = 0.5;    dof2nk[o++] = 2;
    nodes[o][0] = 0.5;    nodes[o][1] = gauss2; dof2nk[o++] = 1;
    nodes[o][0] = gauss1; nodes[o][1] = 0.5;    dof2nk[o++] = 2;

    RealMat T(nbas, nbas);
    for (int k = 0; k < nbas; k++) {
        const double *ip = nodes[k];
        const double *n_k = nk + 2*dof2nk[k];

        o = 0;
        double x = ip[0];
        double y = ip[1];
        double n1 = n_k[0];
        double n2 = n_k[1];

        T(o++,k) = n1;
        T(o++,k) = n2;
        T(o++,k) = x*n1;
        T(o++,k) = x*n2;
        T(o++,k) = y*n1;
        T(o++,k) = y*n2;

        T(o++,k) = x*x*n1;
        T(o++,k) = y*y*n2;
        T(o++,k) = x*y*n1;
        T(o++,k) = x*y*n2;

        T(o++,k) = x*x*y*n1;
        T(o++,k) = x*y*y*n2;
        assert(o == nbas);
    }
    Ti = Eigen::PartialPivLU<RealMat> (T);
}


void FERT1::CalcShapeRef(const double *ip, RealMat &shape) {
    RealMat u(nbas,2);
    int o = 0;

    double x = ip[0];
    double y = ip[1];

    // order 1
    u(o,0) = 1; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = 1; o++;
    u(o,0) = x; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = x; o++;
    u(o,0) = y; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = y; o++;
    // order 2 (not complete)
    u(o,0) = x*x; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = y*y; o++;
    u(o,0) = x*y; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = x*y; o++;
    // order 3 (not complete)
    u(o,0) = x*x*y; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = x*y*y; o++;
    assert(o == nbas);

    shape = Ti.solve(u);
}

void FERT1::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(ip, shape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}

void FERT1::CalcDivShapeRef(const double *ip, RealVec &divshape) {
    int o = 0; 
    RealVec divu(nbas);
    double x = ip[0];
    double y = ip[1];

    divu(o++) = 0; 
    divu(o++) = 0; 
    divu(o++) = 1; 
    divu(o++) = 0; 
    divu(o++) = 0; 
    divu(o++) = 1; 

    divu(o++) = 2*x; 
    divu(o++) = 2*y; 
    divu(o++) = y; 
    divu(o++) = x; 

    divu(o++) = 2*x*y; 
    divu(o++) = 2*x*y; 
    assert(o == nbas);

    divshape = Ti.solve(divu);
}

void FERT1::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(ip, divshape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        divshape.row(i) *= e.edge_sign[i];
    }
}


class FEDGQ1 {
    public:
        FEDGQ1();
        void CalcShapeRef(const double *ip, RealVec &shape) const;
        void CalcShape(const Elem &e, const double *ip, RealVec &shape);
        void CalcGradShapeRef(const double *ip, RealMat &gradshape) const;
        void CalcGradShape(const Elem &e, const double *ip, RealMat &gradshape);
        Coord *nodes = nullptr;
        const int nbas = 4;
        Eigen::PartialPivLU<RealMat> Ti;
};



// ----------------------------- DGQ1 --------------------------------
FEDGQ1::FEDGQ1() {
    nodes = (Coord *) calloc (4, sizeof(Coord));
    double gauss1 = (1-1/sqrt(3))/2;
    double gauss2 = (1+1/sqrt(3))/2;

    int o = 0;
    nodes[o][0] = gauss1; nodes[o][1] = gauss1; o++; 
    nodes[o][0] = gauss2; nodes[o][1] = gauss1; o++; 
    nodes[o][0] = gauss1; nodes[o][1] = gauss2; o++; 
    nodes[o][0] = gauss2; nodes[o][1] = gauss2; o++; 
    assert(o == nbas);

    RealMat T(nbas, nbas);
    for (int k = 0; k < nbas; k++) {
        const double *ip = nodes[k];
        o = 0;
        double x = ip[0];
        double y = ip[1];

        T(o++,k) = 1;
        T(o++,k) = x;
        T(o++,k) = y;
        T(o++,k) = x*y;
        assert(o == nbas);
    }
    Ti = Eigen::PartialPivLU<RealMat> (T);
}



void FEDGQ1::CalcShapeRef(const double *ip, RealVec &shape) const {
    assert(shape.size() == nbas);
    RealVec p(nbas);
    double x = ip[0];
    double y = ip[1];

    int o = 0;
    p(o++) = 1;
    p(o++) = x;
    p(o++) = y;
    p(o++) = x*y;
    assert(o == nbas);
    shape = Ti.solve(p);
}


void FEDGQ1::CalcShape(const Elem &e, const double *ip, RealVec &shape) {
    assert(shape.size() == nbas);
    CalcShapeRef(ip, shape);
}

void FEDGQ1::CalcGradShapeRef(const double *ip, RealMat &gradshape) const {
    assert(gradshape.rows() == nbas);
    RealMat dp(nbas,2);
    double x = ip[0];
    double y = ip[1];
    int o = 0;
    dp(o,0) = 0; dp(o,1) = 0; o++;
    dp(o,0) = 1; dp(o,1) = 0; o++;
    dp(o,0) = 0; dp(o,1) = 1; o++;
    dp(o,0) = y; dp(o,1) = x; o++;
    gradshape = Ti.solve(dp);
}


void FEDGQ1::CalcGradShape(const Elem &e, const double *ip, RealMat &gradshape) {
    assert(gradshape.rows() == nbas);
    CalcGradShapeRef(ip, gradshape);
}



// --------------------------------- Solver using Eigen ---------------------------------
class Solver {
    public:
        Solver(); 
        int total_dofs = 4*N*(2*N+1) + 4*N*N;

        void build_mat();
        void build_rhs();
        void solve();
        void error();
        
        FERT1 *u;
        FEDGQ1 *p;
        /* Elem *e; */
        Quad2d quad2d;
        Quad1d quad1d;
        double u_error_L2 = 0;
        double p_error_L2 = 0;
        RealMat mat;
        RealVec rhs;
        RealVec x;
        /* double geterror(Elem &e, int ibas, RealVec &x); */
};


Solver::Solver() {
    u = new FERT1();
    p = new FEDGQ1();

    cout << "Solver Init\n";
    cout << "u_space: RT1\n";
    cout << "p_space: DGQ1\n";
    
    mat.resize(total_dofs, total_dofs);
    mat.setZero();
    rhs.resize(total_dofs);
    rhs.setZero();
    x.resize(total_dofs);
    x.setZero();
    /* cout << total_dofs << endl; */ 
}




void Solver::build_mat() {
    cout << "Building Mat\n";
    int Nu = 12;
    int Np = 4;
    int Ndof = Nu + Np;
    int Id[Ndof];
    int *Iu = Id, *Ip = Iu + Nu;

    RealMat A(Ndof,Ndof);
    RealMat u_shape(Nu,2), dp_shape(Np,2);
    RealVec p_shape(Np), divu_shape(Nu);
    
    for (int i = 0; i < N; i++) { // iter elem
        for (int j = 0; j < N; j++) { // iter elem
            /* e = new Elem(i, j); // dynamically allocates memory for an Elem object on the heap. */
            Elem e(i,j);           // created an Elem object on the stack.

            /* cout << i << "\t" << j << endl; */
            /* for (int k = 0; k < 4; k++) { */
            /*     cout << e.edge_sign[k] << "\t"; */
            /* } */
            /* cout << endl; */

            for (int k = 0; k < Nu; k++) {
                Iu[k] = map_u(i,j,k);
            }
            for (int k = 0; k < Np; k++) {
                Ip[k] = map_p(i,j,k);
            }

            A.setZero();
            for (int k = 0; k < 9; k++) {
                u->CalcShape(e, quad2d.points[k], u_shape);
                u->CalcDivShape(e, quad2d.points[k], divu_shape);
                p->CalcShape(e, quad2d.points[k], p_shape);
                p->CalcGradShape(e, quad2d.points[k], dp_shape);

                /* if (i == 0 && j == 0) { */
                /*     cout << "u_shape\n" << u_shape << endl; */ 
                /*     cout << "p_shape\n" << p_shape << endl; */ 
                /* } */
                /* if (i == 1 && j == 0) { */
                /*     cout << "u_shape\n" << u_shape << endl; */ 
                /*     cout << "p_shape\n" << p_shape << endl; */ 
                /* } */
                
                for (int m = 0; m < Nu; m++) {
                    // A(u,v)
                    for (int n = 0; n < Nu; n++) {
                        A(m,n) += quad2d.weights[k] * area * (u_shape(m,0)*u_shape(n,0) + u_shape(m,1)*u_shape(n,1));
                    }
                    // A(p,v)
                    for (int n = 0; n < Np; n++) {
                        A(m,n+Nu) += quad2d.weights[k] * area / h * (-divu_shape(m) * p_shape(n));
                    }
                }
                
                for (int m = 0; m < Np; m++) {
                    // A(q,u)
                    for (int n = 0; n < Nu; n++) {
                        A(m + Nu, n) += quad2d.weights[k] * area / h * (p_shape(m) * divu_shape(n));
                    }
                    // A(p,q)
                    for (int n = 0; n < Np; n++) {
                        A(m + Nu, n + Nu) += beta * quad2d.weights[k] * area * (p_shape(m) * p_shape(n));
                    }
                }
            }

            // Apply boundary conditions (Dirichlet)
            /* if (i == 0) { */
            /*     A.row(3).setZero(); */
            /*     A(3,3) = 1; */ 
            /* } */
            /* if (i == N-1) { */
            /*     A.row(1).setZero(); */
            /*     A(1,1) = 1; */
            /* } */
            /* if (j == 0) { */
            /*     A.row(0).setZero(); */
            /*     A(0,0) = 1; */
            /* } */
            /* if (j == N-1) { */
            /*     A.row(2).setZero(); */
            /*     A(2,2) = 1; */
            /* } */


            /* cout << A << endl << endl; */

            // Add element matrix to global matrix
            for (int m = 0; m < Ndof; m++) {
                for (int n = 0; n < Ndof; n++) {
                    mat(Id[m], Id[n]) += A(m,n);
                }
            }
        } // end iter elem
    } // end iter elem
    cout << "Build Mat Done\n";

    // output full matrix
    /* for (int i = 0; i < total_dofs; i++) { */
    /*     printf("%d\t", i); */
    /*     for (int j = 0; j < total_dofs; j++) */ 
    /*         printf("%5.4f\t", mat(i,j)); */
    /*         /1* cout << "%4.2f" << mat(i,j) << "\t"; *1/ */
    /*     cout << endl; */
    /* } */

    // output only u-part matrix
    /* for (int i = 0; i < 2*(N)*(N+1); i++) { */
    /*     printf("%d\t", i); */
    /*     for (int j = 0; j < 2*(N)*(N+1); j++) */ 
    /*         printf("%4.2f\t", mat(i,j)); */
    /*         /1* cout << "%4.2f" << mat(i,j) << "\t"; *1/ */
    /*     cout << endl; */
    /* } */
}


void Solver::build_rhs() {
    cout << "Building Rhs\n";
    /* VecSet(rhs, 0.0); */
    
    int Nu = 12;
    int Np = 4;
    int Ndof = Nu + Np;
    int Id[Ndof];
    int *Iu = Id, *Ip = Iu + Nu;
    
    RealVec F(Ndof);
    RealMat u_shape(Nu, 2);
    RealVec p_shape(Np), divu_shape(Nu);
    double ubc[2];
    double vf;
    
    for (int i = 0; i < N; i++) { // iter elem
        for (int j = 0; j < N; j++) { // iter elem
                                     
            Elem e(i,j);

            for (int k = 0; k < Nu; k++) {
                Iu[k] = map_u(i,j,k);
            }
            for (int k = 0; k < Np; k++) {
                Ip[k] = map_p(i,j,k);
            }
            
            F.setZero();
            for (int k = 0; k < 9; k++) {
                /* u->CalcShape(e, quad2d.points[k], u_shape); */
                u->CalcDivShape(e, quad2d.points[k], divu_shape);
                p->CalcShape(e, quad2d.points[k], p_shape);
                
                Coord ref_coord = { quad2d.points[k][0], quad2d.points[k][1] }; 
                Coord phy_coord;
                e.ref2phy(ref_coord, phy_coord);
                vf = test.func_f(phy_coord);
                
                // F(q)
                for (int m = 0; m < Np; m++) {
                    F(m + Nu) += quad2d.weights[k] * area * (p_shape(m) * vf);
                }
            }
            
            // -(p, v·n)_{\partial \Omega} and dirichlet boundary condition
            if (i == 0) {
                double n[2] = {-1,0};
                for (int k = 0; k < 3; k++) {
                    Coord ref_coord = {0, quad1d.points[k]}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);

                    /* if (j == 1) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */

                    double pbc = test.func_p(phy_coord);
                    u->CalcShape(e, ref_coord, u_shape);
                    for (int m = 0; m < Nu; m++) {
                        F(m) += len * quad1d.weights[k] * (-1) * pbc * (u_shape(m,0) * n[0] + u_shape(m,1) * n[1]);
                    }
                }
                /* test.func_u(0, e.y_center, ubc); */
                /* F(3) = ubc[0] * n[0] + ubc[1] * n[1]; */
            }

            if (i == N-1) {
                double n[2] = {1,0};
                for (int k = 0; k < 3; k++) {
                    Coord ref_coord = {1, quad1d.points[k]}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);

                    /* if (j == 0) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */

                    double pbc = test.func_p(phy_coord);
                    u->CalcShape(e, ref_coord, u_shape);
                    for (int m = 0; m < Nu; m++) {
                        F(m) += len * quad1d.weights[k] * (-1) * pbc * (u_shape(m,0) * n[0] + u_shape(m,1) * n[1]);
                    }
                }
                /* test.func_u(1, e.y_center, ubc); */
                /* F(1) = ubc[0] * n[0] + ubc[1] * n[1]; */
                /* cout << ubc[0] << "\t" << ubc[1] << endl; */
                /* cout << F(1) << endl; */
            }

            if (j == 0) {
                double n[2] = {0,-1};
                for (int k = 0; k < 3; k++) {
                    Coord ref_coord = {quad1d.points[k], 0}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);
                    /* if (i == 1) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */
                    double pbc = test.func_p(phy_coord);
                    u->CalcShape(e, ref_coord, u_shape);
                    for (int m = 0; m < Nu; m++) {
                        F(m) += len * quad1d.weights[k] * (-1) * pbc * (u_shape(m,0) * n[0] + u_shape(m,1) * n[1]);
                    }
                }
                /* test.func_u(e.x_center, 0, ubc); */
                /* F(0) = ubc[0] * n[0] + ubc[1] * n[1]; */
            }

            if (j == N-1) {
                double n[2] = {0,1};
                for (int k = 0; k < 3; k++) {
                    Coord ref_coord = {quad1d.points[k], 1}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);
                    /* if (i == 0) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */
                    double pbc = test.func_p(phy_coord);
                    u->CalcShape(e, ref_coord, u_shape);
                    for (int m = 0; m < Nu; m++) {
                        F(m) += len * quad1d.weights[k] * (-1) * pbc * (u_shape(m,0) * n[0] + u_shape(m,1) * n[1]);
                    }
                }
                /* test.func_u(e.x_center, 1, ubc); */
                /* F(2) = ubc[0] * n[0] + ubc[1] * n[1]; */
            }

            // dirichlet boundary condition
            /* if (i == 0) { */  
            /*     double n[2] = {-1,0}; */
            /*     test.func_u(0, e.y_center, ubc); */
            /*     F(3) = ubc[0] * n[0] + ubc[1] * n[1]; */
            /* } */
            /* if (i == N-1) { */
            /*     double n[2] = {1,0}; */
            /*     test.func_u(1, e.y_center, ubc); */
            /*     F(1) = ubc[0] * n[0] + ubc[1] * n[1]; */
            /* } */
            /* if (j == 0) { */
            /*     double n[2] = {0,-1}; */
            /*     test.func_u(e.x_center, 0, ubc); */
            /*     F(0) = ubc[0] * n[0] + ubc[1] * n[1]; */
            /* } */
            /* if (j == N-1) { */
            /*     double n[2] = {0,1}; */
            /*     test.func_u(e.x_center, 1, ubc); */
            /*     F(2) = ubc[0] * n[0] + ubc[1] * n[1]; */
            /* } */


            for (int m = 0; m < Ndof; m++) {
                rhs(Id[m]) += F(m);
                /* for (int n = 0; n < Ndof; n++) { */
                /*     mat(Id[m], Id[n]) += A(m,n); */
                /* } */
            }

        } // end iter elem
    } // end iter elem
 
    // output rhs vector
    /* for (int m = 0; m < total_dofs; m++) { */
    /*     printf("%d\t%5.3f\n", m, rhs(m)); */
    /* } */
    cout << "Build Rhs Done\n";
}



void Solver::solve() {
    /* x = mat.lu().solve(rhs); */
    x = mat.partialPivLu().solve(rhs);
    /* cout << x << endl; */

    error();
    cout << "u L2 error: " << u_error_L2 << endl;
    cout << "p L2 error: " << p_error_L2 << endl;
}





void Solver::error() {
    double exact_sol;
    double num_sol;
    /* RealVec u_value(2); */
    double u_value[2];
    double Iu[12];
    double Ip[4];
    RealMat u_shape(12,2);
    RealVec p_shape(4);


    // ------------------------------ u L2 error -------------------------------------
    for (int i = 0; i < N; i++) { // iter elem
        for (int j = 0; j < N; j++) { // iter elem
            Elem e(i,j);

            for (int k = 0; k < 12; k++) {
                Iu[k] = map_u(i,j,k);
            }

            for (int k = 0; k < 9; k++) {
                u->CalcShape(e, quad2d.points[k], u_shape);

                Coord ref_coord = {quad2d.points[k][0], quad2d.points[k][1]}; 
                Coord phy_coord;
                e.ref2phy(ref_coord, phy_coord);
                test.func_u(phy_coord, u_value);

                double diff1 = 0, diff2 = 0;
                for (int l = 0; l < 12; l++) {
                    diff1 += x[Iu[l]] * u_shape(l,0);
                    diff2 += x[Iu[l]] * u_shape(l,1);
                }
                diff1 -= u_value[0];
                diff2 -= u_value[1];
                
                u_error_L2 += quad2d.weights[k] * area * (diff1 * diff1 + diff2 * diff2);
            }
        }
    }

    u_error_L2 = sqrt(u_error_L2);




    // ------------------------------ p L2 error -----------------------------------
    for (int i = 0; i < N; i++) { // iter elem
        for (int j = 0; j < N; j++) { // iter elem
            Elem e(i,j);
            
            for (int k = 0; k < 4; k++) {
                Ip[k] = map_p(i,j,k);
            }

            for (int k = 0; k < 9; k++) {
                /* u->CalcShape(e, quad2d.points[k], u_shape); */
                Coord ref_coord = {quad2d.points[k][0], quad2d.points[k][1]}; 
                Coord phy_coord;
                e.ref2phy(ref_coord, phy_coord);
                exact_sol = test.func_p(phy_coord);
                p->CalcShape(e, quad2d.points[k], p_shape);

                double diff = 0;
                for (int l = 0; l < 4; l++) {
                    diff += x[Ip[l]]*p_shape(l); 
                }
                diff -= exact_sol;
                p_error_L2 += quad2d.weights[k] * area * (diff * diff);
            }
        }
    }

    p_error_L2 = sqrt(p_error_L2);
}



int main(int argc, char* argv[]) {
    // ---------------------------- solve ------------------------------ 
    Solver solver;
    solver.build_mat();
    solver.build_rhs();
    solver.solve();
    
    // ---------------------------- FEDGQ1 ------------------------------
    /* FEDGQ1 *p = new FEDGQ1(); */
    /* Elem e(1,0); */
    /* RealVec p_shape(4); */
    /* for (int i = 0; i < 4; i++) { */
    /*     p->CalcShape(e, p->nodes[i], p_shape); */
    /*     cout << "i=" << i << endl << p_shape << endl; */
    /* } */
    /* for (int i = 0; i < p->nbas; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << p->nodes[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */

    // ----------------------------- FERT1 ------------------------------
    /* FERT1 *u = new FERT1(); */
    /* Elem e(1,0); */
    /* RealMat u_shape(12,2); */
    /* for (int i = 0; i < 12; i++) { */
    /*     u->CalcShape(e, u->nodes[i], u_shape); */
    /*     cout << "i=" << i << endl << u_shape << endl; */
    /* } */

    /* FERT1 *u = new FERT1(); */
    /* for (int i = 0; i < 12; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << u->nodes[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */

    // ----------------------------- map_p -------------------------------
    /* for (int i = 0; i < 2; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         for (int ibas = 0; ibas < 4; ibas++) { */
    /*             cout << map_p(i,j,ibas) << endl; */
    /*         } */
    /*     } */
    /* } */
    /* cout << udof_global; */

    // ----------------------------- map_u ----------------------------------
    /* for (int i = 0; i < 2; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         for (int ibas = 0; ibas < 12; ibas++) { */
    /*             cout << map_u(i,j,ibas) << endl; */
    /*         } */
    /*     } */
    /* } */

    // ------------------------------- quadrature ---------------------------------
    /* Quad1d quad1d; */
    /* double integral; */
    /* for (int i = 0; i < quad1d.npoints; i++) { */
    /*     integral += fun1(quad1d.points[i]) * quad1d.weights[i]; */
    /* } */
    /* cout << integral << endl; */
    
    // ------------------------------- nodes ------------------------------
    /* FERT0 u; */
    /* for (int i = 0; i < 4; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << u.nodes[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */

    // --------------------------- test basis function ---------------------------------
    /* FERT0 *u_test = new FERT0(); */
    /* double point[2] = {0.5, 0}; */
    /* RealMat value(4,2); */
    /* u_test->CalcShapeRef(point, value); */
    /* cout << value << endl; */

    /* point[0] = 1.0; */
    /* point[1] = 0.5; */
    /* u_test->CalcShapeRef(point, value); */
    /* cout << value << endl; */

    // ----------------------------- test dof map ------------------------------------
    /* for (int i = 0; i < N; i++) { */
    /*     for (int j = 0; j < N; j++) { */
    /*         for (int ibas = 0; ibas < 4; ibas++) { */
    /*             cout << map_u(i, j, ibas) << "\t"; */
    /*         } */
    /*         cout << map_p(i, j); */
    /*         cout << endl; */
    /*     } */
    /* } */
    
    return 0;
}
