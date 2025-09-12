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

namespace BDRY_MASK {
    enum {
        INTERIOR = 1,
        OWNER = 2,
        REMOTE = 4,
        DIRICHLET = 8,
        NEUMANN = 16
    };
}

typedef int BdryEdge[3];
typedef int Edge[2];

int compare_bdry_edge(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    BdryEdge *e0 = (BdryEdge *) p0;
    BdryEdge *e1 = (BdryEdge *) p1;

    i = (*e0)[0] - (*e1)[0];

    if (i < 0) {
		return -1;
    } else if (i > 0) {
		return 1;
    } else {
		i = (*e0)[1] - (*e1)[1];
		return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }

    /* abort(); */
}


int compare_edge(const void *p0, const void *p1)
/* compares two "INT's (used in qsort and bsearch) */
{
    int i;
    Edge *e0 = (Edge *) p0;
    Edge *e1 = (Edge *) p1;

    i = (*e0)[0] - (*e1)[0];

    if (i < 0) {
		return -1;
    } else if (i > 0) {
		return 1;
    } else {
		i = (*e0)[1] - (*e1)[1];
		return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }
}


class Quad2d_tria {
    public:
        double points[3][2] = {
            {1.0/6.0, 1.0/6.0},
            {2.0/3.0, 1.0/6.0},
            {1.0/6.0, 2.0/3.0}
        };
        double weights[3] = {1.0/6.0, 1.0/6.0, 1.0/6.0};
};


double func2d(double x, double y) {
    /* return (x*x+2)*(y*y+2*y+3); */ 
    /* return x*x; */ 
    /* return x*y; */ 
    /* return x; */ 
    return 1; 
}

namespace ELEM_TYPE {
    enum {
        TRIA,
        QUAD
    };
}

typedef unsigned short Btype;


void CalcLegendre(const int p, const double x, RealVec &u)
{
	// use the recursive definition for [-1,1]:
	// (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
	double z;
	u[0] = 1.;
	if (p == 0) { return; }
	u[1] = z = 2.*x - 1.;
	for (int n = 1; n < p; n++)
    {
        u[n+1] = ((2*n + 1)*z*u[n] - n*u[n-1])/(n + 1);
    }
}

void CalcLegendre(const int p, const double x, RealVec &u, RealVec &d)
{
	// use the recursive definition for [-1,1]:
	// (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
	// for the derivative use, z in [-1,1]:
	// P'_{n+1}(z) = (2*n+1)*P_n(z)+P'_{n-1}(z)
	double z;
	u[0] = 1.;
	d[0] = 0.;
	if (p == 0) { return; }
	u[1] = z = 2.*x - 1.;
	d[1] = 2.;
	for (int n = 1; n < p; n++)
    {
        u[n+1] = ((2*n + 1)*z*u[n] - n*u[n-1])/(n + 1);
        d[n+1] = (4*n + 2)*u[n] + d[n-1];
    }
}


void CalcChebyshev(const int p, const double x, RealVec &u)
{
	// recursive definition, z in [-1,1]
	// T_0(z) = 1,  T_1(z) = z
	// T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
	double z;
	u[0] = 1.;
	if (p == 0) { return; }
	u[1] = z = 2.*x - 1.;
	for (int n = 1; n < p; n++)
    {
        u[n+1] = 2*z*u[n] - u[n-1];
    }
}

void CalcChebyshev(const int p, const double x, RealVec &u, RealVec &d, RealVec &dd)
{
	// recursive definition, z in [-1,1]
	// T_0(z) = 1,  T_1(z) = z
	// T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
	// T'_n(z) = n*U_{n-1}(z)
	// U_0(z) = 1  U_1(z) = 2*z
	// U_{n+1}(z) = 2*z*U_n(z) - U_{n-1}(z)
	// U_n(z) = z*U_{n-1}(z) + T_n(z) = z*T'_n(z)/n + T_n(z)
	// T'_{n+1}(z) = (n + 1)*(z*T'_n(z)/n + T_n(z))
	// T''_{n+1}(z) = (n + 1)*(2*(n + 1)*T'_n(z) + z*T''_n(z)) / n
	double z;
	u[0] = 1.;
	d[0] = 0.;
	dd[0]= 0.;
	if (p == 0) { return; }
	u[1] = z = 2.*x - 1.;
	d[1] = 2.;
	dd[1] = 0;
	for (int n = 1; n < p; n++)
    {
        u[n+1] = 2*z*u[n] - u[n-1];
        d[n+1] = (n + 1)*(z*d[n]/n + 2*u[n]);
        dd[n+1] = (n + 1)*(2.*(n + 1)*d[n] + z*dd[n])/n;
    }
}

int main(int argc, char *argv[]) {
    // ------------------------------------ Chebyshev --------------------------------------------
    RealVec u(10);
    RealVec d(10);
    RealVec dd(10);
    CalcChebyshev(9, 0.2, u, d, dd);
    for (int i = 0; i < 10; i++) {
        printf("%10.8f \t %10.8f \t %10.8f\n", u[i], d[i], dd[i]);
        /* cout << u[i] << "\t" << d[i] << '\t' << dd[i] << endl; */
    }

    // ------------------------------------ Legendre ---------------------------------------------
    /* RealVec u(10); */
    /* RealVec d(10); */
    /* CalcLegendre(9, 0.5, u); */
    /* cout << u << endl << endl; */
    /* CalcLegendre(9, 0.5, u, d); */
    /* cout << u << endl << endl; */
    /* cout << d << endl; */

    // ------------------------------------ sin --------------------------------------------------
    /* cout << M_PI << endl; */
    /* double x = sin(M_PI); */
    /* cout << x << endl; */

    // --------------------------------------- std::vector ---------------------------------------
    /* vector<double> p(10); // 动态数组，自动管理内存 */
    /* for (unsigned int i = 0; i < 10; i++) { */
    /*     p[i] = (double)i; */
    /* } */
    /* // C++11 引入的 range-based for 循环（基于范围的 for 循环） */
    /* for (auto val : p) {     // iterate every item */
    /*     cout << val << " ";  // output this item */
    /* } */
    /* cout << endl; */

    // ---------------------------------- stack memory -------------------------------------------
    /* double p[10]; // 直接在栈上分配数组 */
    /* for (unsigned int i = 0; i < 10; i++) { */
    /*     p[i] = (double)i; */
    /* } */
    /* for (unsigned int i = 0; i < 10; i++) { */
    /*     cout << p[i] << " "; */
    /* } */
    /* cout << endl; */

    // ------------------------------ new for dynamic memory allocation --------------------------
    /* double *p = new double[10]; // 动态分配10个double的内存 */
    /* for (unsigned int i = 0; i < 10; i++) { */
    /*     p[i] = (double)i; // 等价于 *(p+i) = i; */
    /* } */

    /* for (unsigned int i = 0; i < 10; i++) { */
    /*     cout << p[i] << " "; */
    /* } */
    /* cout << endl; */
    /* delete[] p; // 释放内存 */

    // ----------------------------------- command line parameter -------------------------------- 
    /* cout << "the number of parameters: " << argc - 1 << endl; */
    /* cout << "program name: " << argv[0] << endl; */
    /* for (int i = 1; i < argc; ++i) { */
    /*     cout << "parameter " << i << ": " << argv[i] << endl; */
    /* } */

    // ------------------------------------ ELEM_TYPE --------------------------------------------
    /* Btype elem_type = ELEM_TYPE::TRIA; */
    /* cout << (elem_type == 0 ? "TRIA" : "QUAD") << endl; */
    /* elem_type = ELEM_TYPE::QUAD; */
    /* /1* cout << elem_type == 0 ? "TRIA" : "QUAD" << endl; *1/ // Wrong! */
    /* cout << (elem_type == 0 ? "TRIA" : "QUAD") << endl; */

    // -------------------------- resize of RealVec and RealMat ----------------------------------
    /* RealVec vec; */
    /* vec.resize(4); */
    /* vec << 1.0, 2.0, 3.0, 4.0; */
    /* for (int i = 0; i < 4; i++) */
    /*     cout << vec(i) << endl; */
    /* vec.resize(7); */
    /* vec << 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0; */
    /* for (int i = 0; i < 7; i++) */
    /*     cout << vec(i) << endl; */

    // ----------------------------------- tria quad ---------------------------------------------
    /* Quad2d_tria quad2d_tria; */
    /* double sum = 0; */
    /* for (int i = 0; i < 3; i++) { */
    /*     sum += quad2d_tria.weights[i] * func2d(quad2d_tria.points[i][0], quad2d_tria.points[i][1]); */
    /* } */
    /* cout << sum << endl; */


    // --------------------------------------- const ---------------------------------------------
    /* const char a = 'T'; // 'const' variables can be only assigned while declaration. */ 
    /* // a = 'T'; // Totally Wrong! */
    /* cout << a << endl; */

    // ----------------------------------- nullptr in if -----------------------------------------
    /* RealVec *ptr = nullptr; */
    /* if (!ptr) */
    /*     cout << "ptr is nullptr\n"; */

    // ----------------------------------- RealVec indexing --------------------------------------
    /* RealVec a(3); */
    /* a << 1.0, 2.0, 3.0; */
    /* for (int i = 0; i < 3; i++) { */
    /*     cout << a(i) << endl; */
    /* } */
    /* for (int i = 0; i < 3; i++) { */
    /*     cout << a[i] << endl; */
    /* } */

    // ------------------------------- function read_mesh edge -----------------------------------
    /* Edge *edges = new Edge[10](); */
    /* int nedge = 10; */
    /* for (int i = 0; i < 10; i++) { */
    /*     if (i < 3) { */
    /*         edges[i][0] = 5; */
    /*         edges[i][1] = 7; */
    /*     } else if (i < 6) { */
    /*         edges[i][0] = 4; */
    /*         edges[i][1] = 9; */
    /*     } else if (i < 9) { */
    /*         edges[i][0] = 2; */
    /*         edges[i][1] = 4; */
    /*     } else { */
    /*         edges[i][0] = 9; */
    /*         edges[i][1] = 10; */
    /*     } */
    /* } */
    /* cout << "before sort:\n"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << edges[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */
    /* qsort(edges, 10, sizeof(Edge), compare_edge); */
    /* cout << "after sort:\n"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << edges[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */
    /* int count = 0; */ 
    /* { // inside the brace are local variables which won't affect the outside variables */
		/* int i0 = 0; */
		/* for (unsigned int i = i0+1; i < 10; i++) { */
			/* int cmp = compare_edge(edges + i0, edges + i); */
			/* if (cmp < 0) { */
				/* i0++; */
				/* memmove(edges + i0, */
						/* edges + i, */
						/* sizeof(Edge)); */
			/* } */
		/* } */
		/* count = i0 + 1; */
    /* } */
    /* nedge = count; */
    /* cout << "nedge = " << nedge << endl; */
    /* cout << "after remove the duplicates:\n"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         cout << edges[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */

    // ------------------------------------- memmove -----------------------------------------
    // void* memmove(void* dest, const void* src, size_t count); // dest: destination memory location
                                                                 // src: source memory location
                                                                 // count: byte numbers to be copied.
    /* int src[] = {1, 2, 3, 4, 5}; */
    /* int dest[5]; */
    /* memmove(dest, src, sizeof(src)); */
    /* for (int i = 0; i < 5; i++) { */
    /*     cout << dest[i] << endl; // 1 2 3 4 5 */
    /* } */
    /* int data[] = {1, 2, 3, 4, 5}; */

    /* memmove(data + 1, data, 3 * sizeof(int)); // move data[0..2] to data[1..3]（memory overlap） */
    /* for (int i = 0; i < 5; i++) { */
    /*     cout << data[i] << endl; // 1 1 2 3 5 */
    /* } */

    // -------------------------------------- qsort ------------------------------------------
    /* BdryEdge *bdry_edges = new BdryEdge[10](); */ 
    /* int nbdry_edge = 10; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         bdry_edges[i][j] = (10-i) * (10-j); */
    /*         /1* (*(bdry_edges+i))[j] = (10-i) * (10-j); // same as above *1/ */
    /*     } */
    /*     bdry_edges[i][2] = 1; */
    /* } */
    /* cout << "before sort:\n"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 3; j++) { */
    /*         cout << bdry_edges[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */
    /* qsort(bdry_edges, nbdry_edge, sizeof(BdryEdge), compare_bdry_edge); */
    /* cout << "after sort:\n"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     for (int j = 0; j < 3; j++) { */
    /*         cout << bdry_edges[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */

    // ------------------------------- NULL and nullptr --------------------------------------
    /* int *p1 = NULL; // easy to be confused with integer(0 or (void*)0) */
    /* int *p2 = nullptr; // more safer after C++11 standard */
    /* if (p1 == nullptr) { */
    /*     cout << "p2 is nullptr\n" << endl; */
    /* } */

    // --------------------------------- abort -----------------------------------------------
    /* bool success = false; */
    /* if (!success) { */
    /*     cerr << "Fatal error encountered!\n"; */
    /*     abort(); */
    /* } */
    /* cout << "This won't be printed.\n"; */

    // --------------------------- petsc: See testpetsc.c -------------------------------------
    /* PetscInitialize(&argc, &argv, NULL, NULL); */
    /* Mat A; */
    /* PetscInt m = 1000, n = 1000; // Matrix dimension */
    /* MatCreate(PETSC_COMM_WORLD, &A); */
    /* MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, n); */
    /* MatSetType(A, MATAIJ); */
    /* MatSetUp(A); */

    // ------------------------------------- enum ---------------------------------------------
    // user-defined types that consist of a set of named integral constants
    /* enum Color { */
    /*     RED,    // 0 */
    /*     BLUE,   // 1 */
    /*     YELLOW, // 2 */
    /*     GREEN = 5 */
    /* }; */
    /* cout << Color::RED << endl; // this is legal */
    /* Color g = GREEN; */ 
    /* cout << g << endl; */
    /* int color = 2; */
    /* if (color == Color::YELLOW) */ 
    /*     cout << "color is " << color << endl; */

    // --------------------------------- enum in namespace ------------------------------------
    /* cout << BDRY_MASK::OWNER << endl; // namespace must defined as global. */

    // ---------------------------------- unsigned int ----------------------------------------
    /* unsigned int i = 10; */
    /* fprintf(stdout, "%d\n", i); */
    /* fprintf(stdout, "%u\n", i); */

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
