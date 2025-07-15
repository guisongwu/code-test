#include "fe-dg.hpp"

// --------------------- DG0 ----------------------

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


