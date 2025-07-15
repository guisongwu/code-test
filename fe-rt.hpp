#ifndef FE_RT_HPP

#include "header.hpp"

// -------------------------- RT0 ---------------------------------
class FERT0 {
    public:
        FERT0(Elem *e) : elem(e);
        void CalcShapeRef(const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(const double *ip, RealVec &divshape);
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Elem *elem;
        Coord *nodes; 
        Eigen::PartialPivLU<RealMat> Ti;
        int nbas; // 3 for tria
                  // 4 for quad
        const double nk_tria[6] =
        { 0., -1.,
          1.,  1.,
         -1.,  0. };
        const double nk_quad[8] =
        { 0., -1.,
          1.,  0.,
          0.,  1.,
         -1.,  0. };
        vector<int> dof2nk;
};




// -------------------------- RT1 ---------------------------------
class FERT1 {
    public:
        FERT1(Elem *e) : elem(e);
        void CalcShapeRef(const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(const double *ip, RealVec &divshape);
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Coord *nodes; 
        Elem *elem;
        Eigen::PartialPivLU<RealMat> Ti;
        int nbas; // 8  for tria
                  // 12 for quad 
        const double nk_tria[6] =
        { 0., -1.,
          1.,  1.,
         -1.,  0. };
        const double nk_quad[8] =
        { 0., -1.,
          1.,  0.,
          0.,  1.,
         -1.,  0. };
        vector<int> dof2nk;
};


#define FE_RT_HPP
#endif
