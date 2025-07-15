#include <iostream>
#include <cstring>
#include <cassert>
#include <Eigen/Dense>

using namespace std;
typedef unsigned short Btype;
typedef Btype (*BC_FUNCTION)(int mark);

typedef int Edge[2];
typedef int BdryEdge[3];
typedef double Coord[3];

using Real3 = Eigen::Vector3d;
typedef double Real;

#define SortIndex(v0,v1) {                 \
    if (v0 > v1) {                         \
        int _tmp = v0; v0 = v1; v1 = _tmp; \
    }                                      \
}
int nelem = 0;


namespace BDRY_MASK {
	enum{
		INTERIOR  = 1,
		OWNER     = 2,
		REMOTE    = 4,
		DIRICHLET = 8,
		NEUMANN   = 16,
		BDRY_USER1 = 32,
		BDRY_USER2 = 64,
		BDRY_USER3 = 128,
		BDRY_USER4 = 256,
		BDRY_USER5 = 512   //max 65535
	};
}


namespace ELEM_TYPE {
    enum {
        TRIA,
        QUAD
    };
}

class Elem {
    public:
        // Note: call after updateJacobian().
        short elem_type; 
                            
        Real getArea(Real3 *face_normal = nullptr,
                     Real3 *face_df1 = nullptr,
                     Real3 *face_df2 = nullptr);
        Real getEdgeLength(int iside, Real3 *edge_normal = nullptr);

        int index;
        int verts[4];
        int edges[4];
        int neigh[4];  // >=0: elem index
                       // -1:  boundry
                       // -2:  remote neigh
                       // int rneigh[NEdge];

        char ordering[4]; // Ordering of vertices:

        Btype btypes[4];
        int edges_sgn[4];

        // Geom 
        Real area;
        Real edge_length[4];
        Real edge_normal[4][3];  // edge outer normal
        Real face_normal[3];	 // face outer normal

        // Update for each quad point
        // Eigen::Matrix2d dFdx;
        // Eigen::Matrix2d Jac;
        Eigen::Matrix<Real, 3, 2> dFdx;
        Eigen::Matrix<Real, 2, 3> Jac;
        Real Det;
        int sgn;
};



/* class TriaElem { */
/*     public: */
/*         // Note: call after updateJacobian(). */
/*         Real getArea(Real3 *face_normal = nullptr, */
/*                      Real3 *face_df1 = nullptr, */
/*                      Real3 *face_df2 = nullptr); */
/*         Real getEdgeLength(int iside, Real3 *edge_normal = nullptr); */

/*         int index; */
/*         int verts[3]; */
/*         int edges[3]; */
/*         int neigh[3];  // >=0: elem index */
/*                        // -1:  boundry */
/*                        // -2:  remote neigh */
/*                        // int rneigh[NEdge]; */

/*         char ordering[3]; // Ordering of vertices: */

/*         Btype btypes[3]; */
/*         int edges_sgn[3]; */

/*         // Geom */ 
/*         Real area; */
/*         Real edge_length[3]; */
/*         Real edge_normal[3][3];  // edge outer normal */
/*         Real face_normal[3];		 // face outer normal */

/*         // Update for each quad point */
/*         // Eigen::Matrix2d dFdx; */
/*         // Eigen::Matrix2d Jac; */
/*         Eigen::Matrix<Real, 3, 2> dFdx; */
/*         Eigen::Matrix<Real, 2, 3> Jac; */
/*         Real Det; */
/*         int sgn; */
/* }; */



/* class QuadElem { */
/*     public: */
/*         // Note: call after updateJacobian(). */
/*         Real getArea(Real3 *face_normal = nullptr, */
/*                      Real3 *face_df1 = nullptr, */
/*                      Real3 *face_df2 = nullptr); */
/*         Real getEdgeLength(int iside, Real3 *edge_normal = nullptr); */

/*         int index; */
/*         int verts[4]; */
/*         int edges[4]; */
/*         int neigh[4];  // >=0: elem index */
/*                        // -1:  boundry */
/*                        // -2:  remote neigh */
/*                        // int rneigh[NEdge]; */

/*         char ordering[4]; // Ordering of vertices: */

/*         Btype btypes[4]; */
/*         int edges_sgn[4]; */

/*         // Geom */ 
/*         Real area; */
/*         Real edge_length[4]; */
/*         Real edge_normal[4][3];  // edge outer normal */
/*         Real face_normal[3];		 // face outer normal */

/*         // Update for each quad point */
/*         // Eigen::Matrix2d dFdx; */
/*         // Eigen::Matrix2d Jac; */
/*         Eigen::Matrix<Real, 3, 2> dFdx; */
/*         Eigen::Matrix<Real, 2, 3> Jac; */
/*         Real Det; */
/*         int sgn; */
/* }; */




class Grid {
    public:
        /* Grid() : rank(phgRank), nprocs(phgNProcs), comm(MPI_COMM_WORLD), verb(1) {} */
        Grid() : verb(1) {} 

        //
        // Initial grid:
        //    1. coord, elems
        //    2. edges
        //    3. types_xxxx
        //    4. btypes and neighs
        //    5. partition
        //    6. local_grid
        //    7. gref
        //    8. geom
        //    9. statistics
        //
        void read_mesh(const char *name, BC_FUNCTION user_bc_map = nullptr, const char *part_file_name = nullptr);

        virtual void partition(const char *part_file_name = nullptr) {}
        virtual void get_local_grid() {}

        void init_geom(); // isop_map, element area, length, normal, etc
        virtual void init_gref() {} // vert_to_elem, edge_to_elem, edge sign
        void statistic(); // statitics

        virtual void check_conformal() {}
        void check_conformal2();
        void output_gmsh(const char *file_name, bool second_order = false);

        // Geom
        void lambda2xyz(const Elem &e, Real xi, Real eta, Real &x, Real &y, Real &z) const;
        void lambda2xyzDirect(const Elem &e, Real xi, Real eta, Real &x, Real &y, Real &z) const;
        void xyz2lambda(const Elem &e, Real x, Real y, Real &xi, Real &eta) const;
        void xyz2lambda(const Elem &e, Real x, Real y, Real z, Real &xi, Real &eta, Real &yeta) const;
        void updateJacobian(Elem &e, const double *lambda);
        void updateJacobian(Elem &e);
        // void getFaceLambda(int v0, int v1, int v2, const Real *p, Real *lambda) const;
        bool getFaceLambda(Elem &e, int iside, const Real *p, Real *lambda, bool do_swap = true) const;

        // Utils
        static int compare_edge(const void *p0, const void *p1);

        virtual int nvert_global() { return nvert; }
        virtual int nedge_global() { return nedge; }
        virtual int nelem_global() { return nelem; }

        unsigned int until_elem(unsigned int ielem, Btype btype = BDRY_MASK::OWNER );

        // nums
        int nvert;
        int nedge;
        int nelem;
        int ntria_elem;
        int nquad_elem;

        Coord *verts = nullptr;
        Edge *edges = nullptr;
        Elem *elems = nullptr;
        /* TriaElem *tria_elems = nullptr; */
        /* QuadElem *quad_elems = nullptr; */

        Btype *types_vert = nullptr;
        Btype *types_edge = nullptr;
        Btype *types_elem = nullptr;

        BC_FUNCTION bc_map;
        struct gRef {
            unsigned int eind;			// elem index
            unsigned int gind;			// geom index
        };
        std::vector<std::vector<gRef>> vert_to_elem;
        std::vector<std::vector<gRef>> edge_to_elem;


        // Isop
        /* Dof *isop_map = nullptr; */
        /* Real proj_sphere_radius = 1.; */

        // layered
        /* int nlayer; */

        /* int rank; */
        /* int nprocs; */
        /* MPI_Comm comm; */

        int verb;
    private:
};


#define ForElements(g, ielem)						                \
    for (unsigned int ielem = g->until_elem(0); ielem < g->nelem;	\
	 ielem = g->until_elem(ielem + 1))


#define SQUARE(x) ((x)*(x))
#define INNER_PRODUCT(p, q)                     \
    (*(p  ) * *(q  ) +                          \
     *(p+1) * *(q+1) +                          \
     *(p+2) * *(q+2))
#define CROSS_PRODUCT(a, b, n) {				\
        n[0] =  (a[1] * b[2] - b[1] * a[2]);	\
        n[1] = -(a[0] * b[2] - b[0] * a[2]);	\
        n[2] =  (a[0] * b[1] - b[0] * a[1]);	\
    }
#define NORMALIZE(a) {							\
		Real l_ = sqrt(a[0]*a[0]				\
					   + a[1]*a[1]				\
					   + a[2]*a[2]);			\
		a[0] /= l_;								\
		a[1] /= l_;								\
		a[2] /= l_;								\
	}


// ------------------------------------- parallel grid ----------------------------------------
/* class ParallelGrid: public Grid { */
/*     public: */
/*         ParallelGrid() : Grid() {}; */ 

/*         void init_gref();  // override edge sign */
/*                            // */
/*                            // Setup partition: */
/*                            //   part, */
/*                            //   assign owner ship */
/*                            // */ 
/*         void partition(const char *part_file_name = nullptr); */
/*         void save_part(const char *part_file_name); */

/*         // */ 
/*         // Setup local grid: */
/*         //   assign owner index */
/*         // */      
/*         // */ 
/*         void get_local_grid(); */
/*         void check_conformal(); */

/*         std::vector<unsigned int> L2Gmap_vert; */
/*         std::vector<unsigned int> L2Gmap_edge; */
/*         std::vector<unsigned int> L2Gmap_elem; */

/*         std::vector<int > G2Lmap_vert; */
/*         std::vector<int > G2Lmap_edge; */
/*         std::vector<int > G2Lmap_elem; */


/*         std::vector<int > owner_index_vert; */
/*         std::vector<int > owner_rank_vert; */

/*         std::vector<int > owner_index_edge; */
/*         std::vector<int > owner_rank_edge; */

/*         std::vector<int > owner_index_elem; */
/*         std::vector<int > owner_rank_elem; */

/*         int nvert_owned; */
/*         int nedge_owned; */
/*         int nelem_owned; */

/*         std::vector<int > flux_edges; */

/*         class Global { */
/*             public: */

/*                 int nvert; */
/*                 int nedge; */
/*                 int nelem; */

/*                 Coord *verts = nullptr; */
/*                 Edge *edges = nullptr; */
/*                 Elem *elems = nullptr; */

/*                 Btype *types_vert = nullptr; */
/*                 Btype *types_edge = nullptr; */
/*                 Btype *types_elem = nullptr; */

/*                 std::vector<int > owner_index_vert; */
/*                 std::vector<int > owner_rank_vert; */

/*                 std::vector<int > owner_index_edge; */
/*                 std::vector<int > owner_rank_edge; */

/*                 std::vector<int > owner_index_elem; */
/*                 std::vector<int > owner_rank_elem; */

/*                 std::vector<std::vector<gRef>> vert_to_elem; */
/*                 std::vector<std::vector<gRef>> edge_to_elem; */
/*         }; */

/*         Global global; */

/*         // Partition */
/*         std::vector<unsigned int> part_id; */

/*         int nvert_global() {return global.nvert;}; */
/*         int nedge_global() {return global.nedge;}; */
/*         int nelem_global() {return global.nelem;}; */

/* }; */


