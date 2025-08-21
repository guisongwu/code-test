// ====================================================================
//  solve 
//  \vec{u} + \nabla p = 0
//  \nabla \cdot \vec{u} + \beta p = f
//  using RT0/DG0 on longitude-latitude meshes and with mass-lumping.
// ====================================================================
#include <iostream>
#include <vector>
#include <Eigen/Dense>
/* #include <petscksp.h> */
#include <stdio.h>

#define N 8
#define beta 1
#define h (1.0/N)
#define area (h*h)
#define len h
#define First 0
#define Second 1

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
        Elem(int i, int j) : ix(i), jy(j), x_left(i*h), y_down(j*h) {
            if (i > 0) edge_sign[3] = -1;
            if (j > 0) edge_sign[0] = -1;
        }

        void ref2phy(const Coord ref_coord, Coord &phy_coord) {
            phy_coord[0] = ref_coord[0]*h + x_left; 
            phy_coord[1] = ref_coord[1]*h + y_down; 
        }

        int ix, jy;
        int edge_sign[4] = { 1,1,1,1 };
        double x_left;
        double y_down;
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
        Coord *nodes; 
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
    double x = ip[0];
    double y = ip[1];
    int o = 0;
    u(o,0) = 1; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = 1; o++;
    u(o,0) = x; u(o,1) = 0; o++;
    u(o,0) = 0; u(o,1) = y; o++;
    assert(o == nbas);

    shape = Ti.solve(u);
}


void FERT0::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(ip, shape);
    // apply edge signs
    for (int i = 0; i < 4; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
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
    // apply edge signs
    for (int i = 0; i < 4; i++) {
        divshape.row(i) *= e.edge_sign[i];
    }
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
        Solver(); 
        int total_dofs = 2*N*(N+1) + N*N;

        void build_mat();
        void build_rhs();
        void solve();
        /* void step(); */
        void error();
        
        FERT0 *u;
        FEDG0 *p;
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
    u = new FERT0();
    p = new FEDG0();

    cout << "Solver Init\n";
    cout << "u_space: RT0\n";
    cout << "p_space: DG0\n";
    
    mat.resize(total_dofs, total_dofs);
    mat.setZero();
    rhs.resize(total_dofs);
    rhs.setZero();
    x.resize(total_dofs);
    x.setZero();
}




void Solver::build_mat() {
    cout << "Building Mat\n";
    int Nu = 4;
    int Np = 1;
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
            Ip[0] = map_p(i,j);

            A.setZero();
            for (int k = 0; k < 4; k++) {
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
    
    int Nu = 4;
    int Np = 1;
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
            /* if (i > 0) e.edge_sign[3] = -1; */
            /* if (j > 0) e.edge_sign[0] = -1; */

            for (int k = 0; k < Nu; k++) {
                Iu[k] = map_u(i, j, k);
            }
            Ip[0] = map_p(i, j);
            
            F.setZero();
            for (int k = 0; k < 4; k++) {
                /* u->CalcShape(e, quad2d.points[k], u_shape); */
                u->CalcDivShape(e, quad2d.points[k], divu_shape);
                p->CalcShape(e, quad2d.points[k], p_shape);
                
                Coord ref_coord = {quad2d.points[k][0], quad2d.points[k][1]}; 
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
                for (int k = 0; k < 2; k++) {
                    Coord ref_coord = {0, quad1d.points[k]}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);

                    /* if (j == 1) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */

                    double pbc = test.func_p(phy_coord);
                    /* u->CalcShape(e, phy_coord, u_shape); */
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
                for (int k = 0; k < 2; k++) {
                    Coord ref_coord = {1, quad1d.points[k]}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);

                    /* if (j == 0) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */

                    double pbc = test.func_p(phy_coord);
                    /* u->CalcShape(e, phy_coord, u_shape); */
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
                for (int k = 0; k < 2; k++) {
                    Coord ref_coord = {quad1d.points[k], 0}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);
                    /* if (i == 1) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */
                    double pbc = test.func_p(phy_coord);
                    /* u->CalcShape(e, phy_coord, u_shape); */
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
                for (int k = 0; k < 2; k++) {
                    Coord ref_coord = {quad1d.points[k], 1}; 
                    Coord phy_coord;
                    e.ref2phy(ref_coord, phy_coord);
                    /* if (i == 0) */
                    /*     cout << "real physical coord:\n" << phy_coord[0] << endl << phy_coord[1] << endl; */
                    double pbc = test.func_p(phy_coord);
                    /* u->CalcShape(e, phy_coord, u_shape); */
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
    double Id[4];
    RealMat u_shape(4,2);


    // ------------------------------ u L2 error -------------------------------------
    for (int i = 0; i < N; i++) { // iter elem
        for (int j = 0; j < N; j++) { // iter elem
            Elem e(i,j);
            /* if (i > 0) e.edge_sign[3] = -1; */
            /* if (j > 0) e.edge_sign[0] = -1; */

            for (int k = 0; k < 4; k++) {
                Id[k] = map_u(i, j, k);
            }

            for (int k = 0; k < 4; k++) {
                u->CalcShape(e, quad2d.points[k], u_shape);

                Coord ref_coord = {quad2d.points[k][0], quad2d.points[k][1]}; 
                Coord phy_coord;
                e.ref2phy(ref_coord, phy_coord);
                test.func_u(phy_coord, u_value);

                double diff1 = x[Id[0]]*u_shape(0,0) + x[Id[1]]*u_shape(1,0) + x[Id[2]]*u_shape(2,0) + x[Id[3]]*u_shape(3,0) - u_value[0];
                double diff2 = x[Id[0]]*u_shape(0,1) + x[Id[1]]*u_shape(1,1) + x[Id[2]]*u_shape(2,1) + x[Id[3]]*u_shape(3,1) - u_value[1];
                
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
                /* u->CalcShape(e, quad2d.points[k], u_shape); */
                Coord ref_coord = {quad2d.points[k][0], quad2d.points[k][1]}; 
                Coord phy_coord;
                e.ref2phy(ref_coord, phy_coord);
                exact_sol = test.func_p(phy_coord);

                p_error_L2 += quad2d.weights[k] * area * (exact_sol-x[map_p(i,j)]) * (exact_sol-x[map_p(i,j)]);
            }
        }
    }

    p_error_L2 = sqrt(p_error_L2);
}



int main(int argc, char* argv[]) {
    // ---------------------------- solve ------------------------------ 
    Solver solver;
    /* solver.init(); */
    solver.build_mat();
    solver.build_rhs();
    solver.solve();
    
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
