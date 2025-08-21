#include "fe-rt.hpp"


// ------------------------------ RT0 ---------------------------------
FERT0::FERT0() {
    int o = 0;
    RealMat T;
    nodes = (Coord *) calloc (4, sizeof(Coord)); // enough to put either tria or quad nodes.

    // -------------------------- Tria ----------------------------
    int nbas = nbas_[ELEM_TYPE::TRIA];
    dof2nk.resize(nbas);
    T.resize(nbas, nbas);

    nodes[o][0] = 1./2.; nodes[o][1] = 0.0;   dof2nk[o++] = 0;
    nodes[o][0] = 1./2.; nodes[o][1] = 1./2.; dof2nk[o++] = 1;
    nodes[o][0] = 0.0;   nodes[o][1] = 1./2.; dof2nk[o++] = 2;

    for (int k = 0; k < nbas; k++) {
        const double *ip = nodes[k];
        const double *n_k = nk_tria + 2*dof2nk[k];

        o = 0;
        T(o++,k) = n_k[0];
        T(o++,k) = n_k[1];
        T(o++,k) = ip[0] * n_k[0] + ip[1] * n_k[1];
    }
    TiTria = Eigen::PartialPivLU<RealMat> (T);

    // -------------------------- Quad ----------------------------
    nbas = nbas_[ELEM_TYPE::QUAD];
    /* nodes = (Coord *) calloc (nbas, sizeof(Coord)); */
    dof2nk.resize(nbas);
    T.resize(nbas, nbas);

    nodes[o][0] = 1./2.; nodes[o][1] = 0.0;   dof2nk[o++] = 0;
    nodes[o][0] = 1.;    nodes[o][1] = 1./2.; dof2nk[o++] = 1;
    nodes[o][0] = 1./2.; nodes[o][1] = 1.0;   dof2nk[o++] = 2;
    nodes[o][0] = 0.;    nodes[o][1] = 1./2.; dof2nk[o++] = 3;

    for (int k = 0; k < nbas; k++) {
        const double *ip = nodes[k];
        const double *n_k = nk_quad + 2*dof2nk[k];

        o = 0;
        T(o++,k) = n_k[0];
        T(o++,k) = n_k[1];
        T(o++,k) = ip[0] * n_k[0];
        T(o++,k) = ip[1] * n_k[1];
    }
    
    TiQuad = Eigen::PartialPivLU<RealMat> (T);
}




void FERT0::CalcShapeRef(Btype type, const double *ip, RealMat &shape) {
    int nbas = nbas_[type];
    RealMat u(nbas,2);
    int o = 0;
    double x = ip[0];
    double y = ip[1];

    if (type == ELEM_TYPE::TRIA) {
        u(o,0) = 1; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = 1; o++;
        u(o,0) = x; u(o,1) = y; o++;

        /* assert(o == nbas[type]); */
        shape = TiTria.solve(u);
    } else if (type == ELEM_TYPE::QUAD) {
        u(o,0) = 1; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = 1; o++;
        u(o,0) = x; u(o,1) = 0; o++;
        u(o,0) = 0; u(o,1) = y; o++;

        /* assert(o == nbas[type]); */
        shape = TiQuad.solve(u);
    }
}


// TODO
void FERT0::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(e.elem_type, ip, shape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}





void FERT0::CalcDivShapeRef(Btype type, const double *ip, RealVec &divshape) {
    int nbas = nbas_[type];
    RealVec divu(nbas);
    int o = 0;
    double x = ip[0];
    double y = ip[1];

    if (type == ELEM_TYPE::TRIA) {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 2; 

        /* assert(o == nbas[type]); */
        divshape = TiTria.solve(divu);
    } else if (type == ELEM_TYPE::QUAD) {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 
        divu(o++) = 1; 

        /* assert(o == nbas[type]); */
        divshape = TiQuad.solve(divu);
    }
}




// TODO
void FERT0::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(e.elem_type, ip, divshape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}




// ------------------------------ RT1 ---------------------------------
FERT1::FERT1() {
    int o = 0;
    double gauss1 = (1-1/sqrt(3))/2;
    double gauss2 = (1+1/sqrt(3))/2;
    RealMat T;
    nodes = (Coord *) calloc (12, sizeof(Coord)); // enough to put either tria or quad nodes.

    // -------------------------- Tria ----------------------------
    int nbas = nbas_[ELEM_TYPE::TRIA];
    dof2nk.resize(nbas);
    T.resize(nbas, nbas);

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
    TiTria = Eigen::PartialPivLU<RealMat> (T);

    // -------------------------- Quad ----------------------------
    nbas = nbas_[ELEM_TYPE::QUAD];
    /* nodes = (Coord *) calloc (nbas, sizeof(Coord)); */
    dof2nk.resize(nbas);
    T.resize(nbas, nbas);

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
    TiQuad = Eigen::PartialPivLU<RealMat> (T);
}



void FERT1::CalcShapeRef(Btype type, const double *ip, RealMat &shape) {
    int nbas = nbas_[type];
    RealMat u(nbas,2);
    int o = 0;
    double x = ip[0];
    double y = ip[1];

    if (type == ELEM_TYPE::TRIA) {
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

        /* assert(o == nbas[type]); */
        shape = TiTria.solve(u);
    } else if (type == ELEM_TYPE::QUAD) {
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

        /* assert(o == nbas[type]); */
        shape = TiQuad.solve(u);
    }

}


// TODO
void FERT1::CalcShape(const Elem &e, const double *ip, RealMat &shape) {
    CalcShapeRef(e.elem_type, ip, shape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}


void FERT1::CalcDivShapeRef(Btype type, const double *ip, RealVec &divshape) {
    int nbas = nbas_[type];
    RealVec divu(nbas);
    int o = 0;
    double x = ip[0];
    double y = ip[1];

    if (type == ELEM_TYPE::TRIA) {
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 
        divu(o++) = 0; 
        divu(o++) = 0; 
        divu(o++) = 1; 

        divu(o++) = 2*x + x; 
        divu(o++) = y + 2*y; 

        /* assert(o == nbas[type]); */
        divshape = TiTria.solve(divu);
    } else if (type == ELEM_TYPE::QUAD) {
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

        /* assert(o == nbas[type]); */
        divshape = TiQuad.solve(divu);
    }
}



// TODO
void FERT1::CalcDivShape(const Elem &e, const double *ip, RealVec &divshape) {
    CalcDivShapeRef(e.elem_type, ip, divshape);
    // apply edge signs
    for (int i = 0; i < nbas; i++) {
        shape.row(i) *= e.edge_sign[i];
    }
}




/* int main(int argc, char *argv[]) { */

/*     return 0; */
/* } */
