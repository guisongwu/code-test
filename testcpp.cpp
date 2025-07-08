#include <iostream>
#include <cstdio>
#include <cmath>
#include <Eigen/Dense>
using namespace std;

using RealMat = Eigen::MatrixXd;
using RealVec = Eigen::VectorXd;

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

double func1d(double x) {
    return x*x;
}
double func2d(double x, double y) {
    return (x*x+2)*(y*y+2*y+3); 
}

typedef double Coord[2];
/* void ref2phy(const Coord *ref_coord, Coord &phy_coord) { */
/*     double h = 1.0/2.0; */
/*     phy_coord[0] = (*ref_coord)[0]*h + 0.5; */ 
/*     phy_coord[1] = (*ref_coord)[1]*h + 0.0; */ 
/* } */
void ref2phy(const Coord &ref_coord, Coord &phy_coord) {
    double h = 1.0/2.0;
    phy_coord[0] = ref_coord[0]*h + 0.5; 
    phy_coord[1] = ref_coord[1]*h + 0.0; 
}


int main(int argc, char *argv[]) {
    // ---------------------------------- unsigned int ----------------------------------------
    unsigned int i = 10;
    fprintf(stdout, "%d\n", i);
    fprintf(stdout, "%u\n", i);
    // --------------------------------- ref2phy ----------------------------------------------
    /* double points[4][2] = { */
    /*     {1,2}, */
    /*     {2,1}, */
    /*     {0.5, 0.5}, */
    /*     {4,5} */
    /* }; */
    /* Coord ref_coord = {points[2][0], points[2][1]}; */
    /* Coord phy_coord; */
    /* ref2phy(ref_coord, phy_coord); */
    /* cout << phy_coord[0] << "\t" << phy_coord[1] << endl; */

    
    // ------------------------- COORD *verts and COORD *verts[2] -----------------------------
    /* typedef double COORD[3]; */
    /* COORD *vert; */
    /* vert = (COORD *) calloc (2, sizeof(COORD)); */
    /* for (int i = 0; i < 2; i++) { */
    /*     for (int j = 0; j < 3; j++) { */
    /*         vert[i][j] = (double)i*j; */
    /*         printf("%lf\t", vert[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    /* cout << endl; */

    /* COORD vert1 = {1.0, 2.0, 3.0}; */
    /* COORD vert2 = {4.0, 5.0, 6.0}; */
    /* COORD *verts[2] = {&vert1, &vert2}; */
    /* for (int i = 0; i < 2; i++) { */
    /*     for (int j = 0; j < 3; j++) { */
    /*         printf("%lf\t", (*verts[i])[j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    /* cout << endl; */

    /* vert = new COORD[10](); */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 3; j++) { */
    /*         vert[i][j] = (double)i*j; */
    /*         printf("%lf\t", vert[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

    // ----------------------------- eigen solve -------------------------------
    /* RealMat mat(3,3); */
    /* mat << 1.0, 0, 0, */
    /*     0, 1.0, 0, */
    /*     0, 0, 1.0; */
    /* RealVec rhs(3); */
    /* rhs << 4, 5, 6; */
    /* RealVec x = mat.lu().solve(rhs); */
    /* cout << x << endl; */
    
    // --------------------------- numerical integral ----------------------------
    /* Quad1d quad1d; */
    /* Quad2d quad2d; */
    /* double sum = 0; */
    /* for (int i = 0; i < 2; i++) { */
    /*     sum += quad1d.weights[i] * func1d(quad1d.points[i]); */
    /* } */
    /* for (int i = 0; i < 4; i++) { */
    /*     sum += quad2d.weights[i] * func2d(quad2d.points[i][0], quad2d.points[i][1]); */
    /* } */
    /* cout << sum << endl; */

    // --------------------------- array and pointer -----------------------------
    /* double *weights = {1.0, 1.0};      ------------ wrong! */
    /* double weights[2] = {1.0, 1.0};    ------------ right! */
    // --------------------------- stdout and stderr ------------------------------- 
    /* std::cout << "This goes to stdout" << std::endl; */
    /* std::cerr << "Error: Something went wrong!" << std::endl; */
    return 0;
}
