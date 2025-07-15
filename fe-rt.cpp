#include "fe-rt.hpp"


// ------------------------------ RT0 ---------------------------------
FERT0::FERT0(Elem *e) {
    int o = 0;

    if (elem->elem_type == 'T') {
        nbas = 3;
        nodes = (Coord *) calloc (nbas, sizeof(Coord));
        dof2nk.resize(nbas);
        nodes[o][0] = 1./2.; nodes[o][1] = 0.0;   dof2nk[o++] = 0;
        nodes[o][0] = 1./2.; nodes[o][1] = 1./2.; dof2nk[o++] = 1;
        nodes[o][0] = 0.0;   nodes[o][1] = 1./2.; dof2nk[o++] = 2;
        RealMat T(nbas, nbas);
        for (int k = 0; k < nbas; k++) {
            const double *ip = nodes[k];
            const double *n_k = nk_tria + 2*dof2nk[k];

            o = 0;
            T(o++,k) = n_k[0];
            T(o++,k) = n_k[1];
            T(o++,k) = ip[0] * n_k[0] + ip[1] * n_k[1];
        }
    } else if (elem->elem_type == 'Q') {
        nbas = 4;
        nodes = (Coord *) calloc (nbas, sizeof(Coord));
        dof2nk.resize(nbas);
        nodes[o][0] = 1./2.; nodes[o][1] = 0.0;   dof2nk[o++] = 0;
        nodes[o][0] = 1.;    nodes[o][1] = 1./2.; dof2nk[o++] = 1;
        nodes[o][0] = 1./2.; nodes[o][1] = 1.0;   dof2nk[o++] = 2;
        nodes[o][0] = 0.;    nodes[o][1] = 1./2.; dof2nk[o++] = 3;
        RealMat T(nbas, nbas);
        for (int k = 0; k < nbas; k++) {
            const double *ip = nodes[k];
            const double *n_k = nk_quad + 2*dof2nk[k];

            o = 0;
            T(o++,k) = n_k[0];
            T(o++,k) = n_k[1];
            T(o++,k) = ip[0] * n_k[0];
            T(o++,k) = ip[1] * n_k[1];
        }
    } 
    
    Ti = Eigen::PartialPivLU<RealMat> (T);
}



void FERT0::CalcShapeRef(const double *ip, RealMat &shape) {
    RealMat u(nbas,2);
    int o = 0;
    double x = ip[0];
    double y = ip[1];

    if (elem->elem_type == 'T') {
        u(o,0) = 1; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = 1; o++;
        u(o,0) = x; u(o,1) = y; o++;
    } else if (elem->elem_type == 'Q') {
        u(o,0) = 1; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = 1; o++;
        u(o,0) = x; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = y; o++;
    }

    assert(o == nbas);
    shape = Ti.solve(u);
}


// TODO
void FERT0::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(ip, shape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}



void FERT0::CalcDivShapeRef(const double *ip, RealVec &divshape) {
    int o = 0; 
    RealVec divu(nbas);
    if (elem->elem_type == 'T') {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 2; 
    } else if (elem->elem_type == 'Q') {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 
        divu(o++) = 1; 
    }
    assert(o == nbas);
    divshape = Ti.solve(divu);
}


// TODO
void FERT0::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(ip, divshape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}




// ------------------------------ RT1 ---------------------------------
FERT1::FERT1(Elem *e) {
    int o = 0;
    double gauss1 = (1-1/sqrt(3))/2;
    double gauss2 = (1+1/sqrt(3))/2;

    if (elem->elem_type == 'T') {
        nbas = 8;
        nodes = (Coord *) calloc (nbas, sizeof(Coord));
        dof2nk.resize(nbas);

        // edge
        nodes[o][0] = gauss1; nodes[o][1] = 0.0;    dof2nk[o++] = 0;
        nodes[o][0] = gauss2; nodes[o][1] = 0.0;    dof2nk[o++] = 0;
        nodes[o][0] = gauss2; nodes[o][1] = gauss1; dof2nk[o++] = 1;
        nodes[o][0] = gauss1; nodes[o][1] = gauss2; dof2nk[o++] = 1;
        nodes[o][0] = 0.0;    nodes[o][1] = gauss2; dof2nk[o++] = 2;
        nodes[o][0] = 0.0;    nodes[o][1] = gauss1; dof2nk[o++] = 2;
        // interior
        nodes[o][0] = 1./3.;  nodes[o][1] = 1./3.;  dof2nk[o++] = 0;
        nodes[o][0] = 1./3.;  nodes[o][1] = 1./3.;  dof2nk[o++] = 2;

        RealMat T(nbas, nbas);
        for (int k = 0; k < nbas; k++) {
            const double *ip = nodes[k];
            const double *n_k = nk_tria + 2*dof2nk[k];

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

            T(o++,k) = x*x*n1 + x*y*n2;
            T(o++,k) = x*y*n1 + y*y*n2;

            assert(o == nbas);
        }
    } else if (elem->elem_type == 'Q') {
        nbas = 12;
        nodes = (Coord *) calloc (nbas, sizeof(Coord));
        dof2nk.resize(nbas);

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
            const double *n_k = nk_quad + 2*dof2nk[k];

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
    }

    Ti = Eigen::PartialPivLU<RealMat> (T);
}



void FERT1::CalcShapeRef(const double *ip, RealMat &shape) {
    RealMat u(nbas,2);
    int o = 0;

    double x = ip[0];
    double y = ip[1];

    if (elem->elem_type == 'T') {
        // order 1
        u(o,0) = 1; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = 1; o++;
        u(o,0) = x; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = x; o++;
        u(o,0) = y; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = y; o++;
        // order 2 (not complete)
        u(o,0) = x*x; u(o,1) = x*y; o++;
        u(o,0) = x*y; u(o,1) = y*y; o++;

    } else if (elem->elem_type == 'Q') {
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
    }
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

    if (elem->elem_type == 'T') {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 

        divu(o++) = 2*x + x; 
        divu(o++) = y + 2*y; 

    } else if (elem->elem_type == 'Q') {
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
    }
    assert(o == nbas);
    divshape = Ti.solve(divu);
}



void FERT1::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(ip, divshape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}




int main(int argc, char *argv[]) {

    return 0;
}
