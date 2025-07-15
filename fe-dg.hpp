#ifndef FE_DG_HPP


#include "header.hpp"
// -------------------------- DG0 ---------------------------------
class FEDG0 {
    public:
        FEDG0();
        void CalcShapeRef(const double *ip, RealVec &shape) const;
        void CalcShape(const Elem &e, const double *ip, RealVec &shape);
        void CalcGradShapeRef(const double *ip, RealMat &gradshape) const;
        void CalcGradShape(const Elem &e, const double *ip, RealMat &gradshape);
        Coord *nodes = nullptr;
        const int nbas = 1;
};


// -------------------------- DGQ1 ---------------------------------
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


#define FE_DG_HPP
#endif
