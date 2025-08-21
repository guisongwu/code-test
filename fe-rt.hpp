#ifndef FE_RT_HPP
#include "header.hpp"

// -------------------------- RT0 ---------------------------------
namespace ELEM_TYPE {
	enum{
		TRIA, // 0
		QUAD  // 1
	};
}



class FERT0 {
    public:
        FERT0();
        void CalcShapeRef(Btype type, const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(Btype type, const double *ip, RealVec &divshape);
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Coord *nodes; // common use

        Eigen::PartialPivLU<RealMat> TiTria;
        Eigen::PartialPivLU<RealMat> TiQuad;

        int nbas_[2] = { 3, 4 }; // 3 for tria
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
        FERT1();
        void CalcShapeRef(Btype type, const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(Btype type, const double *ip, RealVec &divshape);
        void CalcDivShape(const Elem &e, const double *ip, RealVec &divshape);
        Coord *nodes; // common use
        Eigen::PartialPivLU<RealMat> TiTria;
        Eigen::PartialPivLU<RealMat> TiQuad;
        int nbas_[2] = { 8, 12 }; // 8  for tria
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
