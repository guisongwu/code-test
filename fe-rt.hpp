#include <Eigen/Dense>
using RealMat = Eigen::MatrixXd;
using RealVec = Eigen::VectorXd;
typedef double Coord[2];



class FERT0 {
    public:
        FERT0();
        void CalcShapeRef(const double *ip, RealMat &shape);
        void CalcShape(const Elem &e, const double *ip, RealMat &shape);
        void CalcDivShapeRef(const double *ip, RealVec &divshape);
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



