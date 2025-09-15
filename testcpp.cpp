#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
using namespace std;

using RealMat = Eigen::MatrixXd;
using RealVec = Eigen::VectorXd;
typedef double Real;

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


typedef double Coord[3];
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

void phgInfo(int verbose_level, const char *fmt, ...) 
{
    // 检查详细级别是否满足输出条件
    if (verbose_level < PHG_MIN_VERBOSE_LEVEL) {
        return;
    }

    va_list ap;
    va_start(ap, fmt);

    // 直接格式化输出，避免缓冲区溢出风险
    /* if (log_file != NULL) { */
    /*     vfprintf(log_file, fmt, ap); */
    /*     fflush(log_file); */
    /* } else { */
    /*     vfprintf(stdout, fmt, ap); */
    /*     fflush(stdout); */
    /* } */
    vfprintf(stdout, fmt, ap);
    fflush(stdout);

    va_end(ap);
}

class Elem {
  public:

	// Note: call after updateJacobian().
    Real getArea(Real3 *face_normal = nullptr,
				 Real3 *face_df1 = nullptr,
				 Real3 *face_df2 = nullptr);
    Real getEdgeLength(int iside, Real3 *edge_normal = nullptr);

    int index;
    int verts[NVert];
    int edges[NEdge];
    int neigh[NEdge];  // >=0: elem index
		       // -1:  boundry
		       // -2:  remote neigh
    //int rneigh[NEdge];

    char ordering[NVert]; // Ordering of vertices:

    
    Btype btypes[NEdge];
    //char ordering[NVert];
    int edges_sgn[NEdge];

    // Geom 
    Real area;
    Real edge_length[NEdge];
    Real edge_normal[NEdge][3];  // edge outer normal
    Real face_normal[3];		 // face outer normal
	
    
    // Update for each quad point
    // Eigen::Matrix2d dFdx;
    // Eigen::Matrix2d Jac;
    Eigen::Matrix<Real, 3, 2> dFdx;
    Eigen::Matrix<Real, 2, 3> Jac;
    Real Det;
    int  sgn;

	// // Linearized 
    // Eigen::Matrix<Real, 3, 2> dFdxLn;
    // Eigen::Matrix<Real, 2, 3> JacLn;
    // Real DetLn;

  private:
};


class Elem {
    public:
        Elem();
        double verts[4][3];
        RealMat dFdx;
        RealMat Jac;
        double Det;
        int sgn;
};


Elem::Elem() : dFdx(3,2), Jac(2,3)
{
    double external_data[4][3] = {
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 4.0},
        {2.0, 2.0, 4.0},
        {0.0, 2.0, 0.0}
    };
    for(int i=0; i<4; ++i) {
        for(int j=0; j<3; ++j) {
            verts[i][j] = external_data[i][j];
        }
    }
}


void linear_proj(Real xi, Real eta, RealVec &shape) {
    shape << (1.-xi)*(1.-eta), 
             xi*(1.-eta), 
             xi*eta, 
             (1.-xi)*eta;
}


void updateJacobian(Elem &e, const double *lambda)
//
// Update elem's dFdx[3][2], Jac[2][3], det
// with simple linear projection
//
{
    int N = 4;
    Real xi = lambda[0];
    Real eta = lambda[1];
    RealVec shape(N), xyz(3);
    linear_proj(xi, eta, shape);

    Real x0 = e.verts[0][0];
    Real x1 = e.verts[1][0];
    Real x2 = e.verts[2][0];
    Real x3 = e.verts[3][0];
    
    Real y0 = e.verts[0][1];
    Real y1 = e.verts[1][1];
    Real y2 = e.verts[2][1];
    Real y3 = e.verts[3][1];

    Real z0 = e.verts[0][2];
    Real z1 = e.verts[1][2];
    Real z2 = e.verts[2][2];
    Real z3 = e.verts[3][2];

    e.dFdx << (1.-eta)*(x1-x0) + eta*(x2-x3), (1.-xi)*(x3-x0) + xi*(x2-x1),
              (1.-eta)*(y1-y0) + eta*(y2-y3), (1.-xi)*(y3-y0) + xi*(y2-y1),
              (1.-eta)*(z1-z0) + eta*(z2-z3), (1.-xi)*(z3-z0) + xi*(z2-z1);

    xyz(0) = x0*shape(0) + x1*shape(1) + x2*shape(2) + x3*shape(3);
    xyz(1) = y0*shape(0) + y1*shape(1) + y2*shape(2) + y3*shape(3);
    xyz(2) = z0*shape(0) + z1*shape(1) + z2*shape(2) + z3*shape(3);


    //      Mult(PointMat, dshape, dFdx);
    //           3 x nbas  nbas x 2    3x2
    //
	// Det = || dF1 x dF2 ||
	// Jac = (A^T A)^{-1} A^T
	// 

	/* if (Params->use_proj_sphere) { */
	/* 	const Real R = proj_sphere_radius; */
	/* 	Real x = xyz(0); */
	/* 	Real y = xyz(1); */
	/* 	Real z = xyz(2); */
	/* 	Real r = sqrt(x*x + y*y + z*z); */
			
	/* 	RealMat	Jx(3, 3); */

	/* 	Jx << y*y + z*z, -x*y	 , -x*z, */
	/* 	    -x*y	  , x*x + z*z, -y*z, */
	/* 	    -x*z	  , -y*z	 , x*x + y*y; */

	/* 	Jx *= R / (r*r*r); */

	/* 	e.dFdx = Jx * e.dFdx; */
	/* } */

	
    {
		//
		//
		// Physical 
		// 
		// 
		Real d[6] = {e.dFdx(0,0), e.dFdx(1,0), e.dFdx(2,0),  // DF^T (2x3)
					 e.dFdx(0,1), e.dFdx(1,1), e.dFdx(2,1)}; 
		Real ad[6];
		Real e_, g_, f_;
		// DF^T DF = [e f]
		//           [f g]
		// adj(...) = [g, -f]          
		//            [-f, e]
		e_ = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
		g_ = d[3]*d[3] + d[4]*d[4] + d[5]*d[5];
		f_ = d[0]*d[3] + d[1]*d[4] + d[2]*d[5];

		//  adj(DF^T x DF) x DF^T
		ad[0] = d[0]*g_ - d[3]*f_;	
		ad[1] = d[3]*e_ - d[0]*f_;
		ad[2] = d[1]*g_ - d[4]*f_;
		ad[3] = d[4]*e_ - d[1]*f_;
		ad[4] = d[2]*g_ - d[5]*f_;
		ad[5] = d[5]*e_ - d[2]*f_;

		e.Det = e_ * g_ - f_ * f_;	// det^2
		e.Jac << ad[0], ad[2], ad[4], ad[1], ad[3], ad[5];
		e.Jac /= e.Det; // inv(DF^T x DF) x DF
		e.Det = sqrt(e.Det);
		e.sgn = 1;
	}

    return;
}

void lambda2xyzDirect(const Elem &e, Real xi, Real eta,
				  Real &x, Real &y, Real &z) 
//
// lambda -> xyz using simple linear projection
// 
{
    Real N0 = (1.-xi) * (1.-eta);
    Real N1 = xi * (1.-eta);
    Real N2 = xi * eta;
    Real N3 = (1.-xi) * eta;

    // Use simple linear
    x = e.verts[0][0] * N0
	    + e.verts[1][0] * N1
	    + e.verts[2][0] * N2
	    + e.verts[3][0] * N3;
    y = e.verts[0][1] * N0
	    + e.verts[1][1] * N1
	    + e.verts[2][1] * N2
	    + e.verts[3][1] * N3;
    z = e.verts[0][2] * N0
	    + e.verts[1][2] * N1
	    + e.verts[2][2] * N2
	    + e.verts[3][2] * N3;
    return;
}


void xyz2lambdaDirect(const Elem &e, Real x, Real y, Real z,
				 Real &xi, Real &eta)
//
//  xyz -------------> lambda
//        Affine map
//        
{
	//
	// x - x0
	// 
	x -= e.verts[0][0];
	y -= e.verts[0][1];
	z -= e.verts[0][2];

    xi   = e.Jac(0, 0) * x + e.Jac(0, 1) * y + e.Jac(0, 2) * z;
    eta  = e.Jac(1, 0) * x + e.Jac(1, 1) * y + e.Jac(1, 2) * z;
	return;
}

typedef Btype (*BC_FUNCTION)(int mark);


static Btype default_bc_map(int mark)
{
    return BDRY_MASK::DIRICHLET;
}

static int get_token_ret = 0;
bool get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (true) {
        memset(token, 0, 100);
        if (fscanf(fp, "%s", token) != 1)  // everytime 'fscanf' is called, 
                                           // the pointer will move to the end of this string.
            return (get_token_ret = false);
        if (token[0] != '#')  
            break;
        /* skip to newline */
        do {
            if ((c = fgetc(fp)) == EOF)
                return (get_token_ret = false);
        } while (c != '\n');
    }
    if ((p = strchr(token, '#')) != nullptr)
        *p = '\0';
    return (get_token_ret = true);
}


#define MAX_TOKEN_LEN 1024
/* void read_mesh(const char *mesh_file_name, BC_FUNCTION user_bc_map, const char *part_file_name) */
void read_mesh(const char *mesh_file_name, BC_FUNCTION user_bc_map)
// ------------------------------------------------------------
//
//
// Read in 2D mesh
//
//
// ------------------------------------------------------------
{
    int nbdry_edge;
    BdryEdge *bdry_edges;
    BC_FUNCTION bc_map;

    if (user_bc_map == NULL)
		bc_map = default_bc_map;
    else
		bc_map = user_bc_map;

    /* if (verb >=1 && phgRank == 0) */
    cout << "Reading Mesh " << mesh_file_name << endl;

    Coord *verts = NULL;
    int nvert;
    // int *bdry_edges = NULL;
    // int *faces = NULL;
    
#define READ_NUMBER													\
    if (!get_token(fp, token)) strcpy(token, "End");				\
    if (isalpha((int)(token[0]))) {									\
		fprintf(stderr, "fewer entries (%d) than expected.\n", i);	\
		break;														\
    }
#undef ERROR
#define ERROR	{n = __LINE__; goto error;}
    
    FILE *fp ;
    char token[MAX_TOKEN_LEN]; // 1024
    int n;
    if ((fp = fopen(mesh_file_name, "r")) == NULL) {
		printf("can't open mesh file \"%s\"!\n", mesh_file_name);
		/* phgAbort(0); */
		return;
    }

    token[0] = '\0';
    if (!get_token(fp, token) || strcasecmp(token, "$MeshFormat") ||
		!get_token(fp, token) || strcmp(token, "2.2") ||
		!get_token(fp, token) || strcmp(token, "0") ||
		!get_token(fp, token) || strcmp(token, "8") ||
		!get_token(fp, token) || strcasecmp(token, "$EndMeshFormat")) {
		n = __LINE__;
	error:
		fprintf(stderr, "(%s:%d) invalid input file: ret=%s, token=%s\n",
				__FILE__, n, get_token_ret ? "TRUE" : "FALSE", token);
    }

    while (true) {
		if (!get_token(fp, token))
			break;
		//next_token:
		if (!strcasecmp(token, "$Nodes")) {
			if (!get_token(fp, token))
				ERROR;
			n = atoi(token);
			/* if (verb >= 1 && phgRank == 0) */
            fprintf(stderr, "number of vertices: %d\n", n);
			verts = new Coord[n]();

			unsigned int i;
			for (i = 0; i < n; i++) {
				READ_NUMBER;	// Unused
				READ_NUMBER;
				verts[i][0] = atof(token);
				if (!get_token(fp, token))
					ERROR;
				verts[i][1] = atof(token);
				//phgInfo(0, "vert: %d %lf %lf\n", i, verts[i][0],verts[i][1]);
				if (!get_token(fp, token))
					ERROR;
				verts[i][2] = atof(token); 
			}
			assert(i == n);
			// if (i < n)
			// 	goto next_token;
			nvert = n;

			if (!get_token(fp, token))
				ERROR;
			assert(!strcasecmp(token, "$EndNodes"));
		}
		else if (!strcasecmp(token, "$Elements")) {
			if (!get_token(fp, token))
				ERROR;
			n = atoi(token);
			/* if (verb >= 1 && phgRank == 0) */
            fprintf(stderr, "number of elements: %d\n", n);

			bdry_edges = new BdryEdge[n]();
			elems = new Elem[n]();
			int nedge_read = 0, nface_read = 0;
			int elem_type = 0;

			unsigned int i;
			for (i = 0; i < n; i++) {
				int phy_id = 0, entity_id = 0;
				READ_NUMBER;	// id 
				assert(atoi(token) > 0);

				READ_NUMBER;	// # element type
				elem_type = atoi(token);
				// edge or quadrilateral
				assert(elem_type == 1 || elem_type == 3);

				READ_NUMBER;	// # of tags
				int ntag = atoi(token);
				int itag = 0;
				if (itag < ntag) {
					READ_NUMBER;
					entity_id = atoi(token);
					itag++;
				}
				if (itag < ntag) {
					READ_NUMBER;
					phy_id = atoi(token);
					itag++;
				}
				for (; itag < ntag; itag++) {
					READ_NUMBER;  //  rest tag not used
				}

				if (elem_type == 1) {
					// 
					// edges
					//
					int v0, v1;
					READ_NUMBER;
					v0 = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					v1 = atoi(token) - 1;
					SortIndex(v0, v1);

					bdry_edges[nedge_read][0] = v0;
					bdry_edges[nedge_read][1] = v1;
					bdry_edges[nedge_read][2] = phy_id;  // edge mark
					nedge_read++;
				}
				else if (elem_type == 3) {
					//
					// faces
					// 
					READ_NUMBER;
					Elem *e = &elems[nface_read];
					e->verts[0] = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					e->verts[1] = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					e->verts[2] = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					e->verts[3] = atoi(token) - 1;
					e->index = nface_read;
					//e->mark = phy_id;

					// if (phy_id == 0) {
					// 	phgInfo(2, "Discard elem\n");
					// }
					// else {
					nface_read++;
					//}
				}
				else {
					ERROR;
				}
			}
			nelem = nface_read;
			nbdry_edge = nedge_read;
			assert(i == n);
			// if (i < n)
			// 	goto next_token;
	    
			if (!get_token(fp, token))
				ERROR;
			assert(!strcasecmp(token, "$EndElements"));
		}
    }
    qsort(bdry_edges, nbdry_edge, sizeof(BdryEdge),
		  compare_bdry_edge);



    // ------------------------------------------------------------
    //
    // Assign edge indices
    // 
    // ------------------------------------------------------------

    edges	= new Edge[nelem * 4]();
    for (unsigned int i = 0; i < nelem; i++) {
		int v0, v1;
		for (unsigned int j = 0; j < 4; j++) {
			v0 = elems[i].verts[GetEdgeVert(j, 0)];
			v1 = elems[i].verts[GetEdgeVert(j, 1)];

			SortIndex(v0, v1);
	    
			edges[i*4+j][0] = v0;
			edges[i*4+j][1] = v1;
		}

		// Ordering
		{
			int V0, V1, V2, V3, v0, v1, v2, v3;
			V0 = elems[i].verts[0];
			V1 = elems[i].verts[1];
			V2 = elems[i].verts[2];
			V3 = elems[i].verts[3];

			v0 = v1 = v2 = v3 = 0;
			(V0 > V1) ? v0++ : v1++;
			(V0 > V2) ? v0++ : v2++;
			(V0 > V3) ? v0++ : v3++;
			(V1 > V2) ? v1++ : v2++;
			(V1 > V3) ? v1++ : v3++;
			(V2 > V3) ? v2++ : v3++;

			elems[i].ordering[0] = v0;
			elems[i].ordering[1] = v1;
			elems[i].ordering[2] = v2;
			elems[i].ordering[3] = v3;
		}
    }
    // for (i = 0; i < 3*nelem; i++)
    // 	phgInfo(2, "edge %d: %d %d\n", i, edges[i][0], edges[i][1]);

    qsort(edges, nelem*4, sizeof(Edge), compare_edge);

    int count = 0; 
    {
		int i0 = 0;
		for (unsigned int i = i0+1; i < 4*nelem; i++) {
			int cmp = compare_edge(edges + i0,
								   edges + i);
			if (cmp < 0) {
				i0++;
				memmove(edges + i0,
						edges + i,
						sizeof(Edge));
			}
		}
		count = i0 + 1;
    }
    nedge = count;


    // for (unsigned int i = 0; i < nedge; i++) 
    // 	phgInfo(2, "edge: %d %d %d\n", i,
    // 	       edges[i][0],edges[i][1]);


#define SortIndex(v0, v1)  {			\
	if (v0 > v1) {				\
	    int _tmp = v0; v0 = v1; v1 = _tmp;	\
	}					\
    }

    // ------------------------------------------------------------
    //
    // Boundary makers
    //
    // ------------------------------------------------------------
    types_vert = new Btype[nvert](); 
    types_edge = new Btype[nedge]();
    types_elem = new Btype[nelem]();

    Edge *edge2elem = new Edge[nedge];
    for (unsigned int i = 0; i < nedge; i++) 
		edge2elem[i][0] = edge2elem[i][1] = -1;

    for (unsigned int ielem = 0; ielem < nelem; ielem++) {
		Elem *e = &elems[ielem];

		for (unsigned int j = 0; j < NEdge; j++) {
			int v0 = e->verts[GetEdgeVert(j, 0)];
			int v1 = e->verts[GetEdgeVert(j, 1)];

			SortIndex(v0, v1);

			Edge edge = {v0, v1}, *p; // Edge *p;

			p = (Edge *) bsearch(&edge, edges,
								 nedge, sizeof(Edge),
								 compare_edge);
			assert (p != NULL);
			int iedge = p - edges;
			e->edges[j] = iedge;

			if (edge2elem[iedge][0] == -1)
				edge2elem[iedge][0] = ielem;
			else if (edge2elem[iedge][1] == -1)
				edge2elem[iedge][1] = ielem;
			else
				/* phgAbort(0); */
                abort();
                

			e->neigh[j] = -1;
		}
    }


    

    //
    // make neighs
    //
    for (unsigned int iedge = 0; iedge < nedge; iedge++) {
		//
		// Search boundary info in mesh
		// 
		BdryEdge edge, *p;
		int mark = -1;	// defualt mark -1

		int v0_ = edges[iedge][0];
		int v1_ = edges[iedge][1];
		SortIndex(v0_, v1_);

		edge[0] = v0_;
		edge[1] = v1_;
		edge[2] = -1;
		//phgInfo(2, "search: %d %d\n", edge[0], edge[1]);

		p = (BdryEdge *) bsearch(&edge, bdry_edges,
								 nbdry_edge, sizeof(BdryEdge),
								 compare_bdry_edge);
		if (p != NULL) {
			mark = (*p)[2];
		}

		Btype btype = (mark >= 0) ? bc_map(mark) : 0;
		types_edge[iedge] = btype;
		// phgInfo(0, "edge %d %d %d\n", iedge, mark, btype);
		// if (btype & BDRY_MASK::DIRICHLET)
		//     phgInfo(0, "(D %d)", iedge);
		//     //phgInfo(2, "(D %d-%d )", e0->index, j0);


		//
		// Update elem edges info
		// 
		int ielem0 = edge2elem[iedge][0];
		int ielem1 = edge2elem[iedge][1];

		if (ielem0 >= 0 && ielem1 >= 0) {
			// Interior edge
			Elem *e0 = &elems[ielem0];
			Elem *e1 = &elems[ielem1];
			int j0, j1;
	    
			for (j0 = 0; j0 < NEdge; j0++)
				if (e0->edges[j0] == iedge)
					break;
			assert(j0 < NEdge);

			for (j1 = 0; j1 < NEdge; j1++)
				if (e1->edges[j1] == iedge)
					break;
			assert(j1 < NEdge);

			types_edge[iedge] |= BDRY_MASK::INTERIOR;	    
			e0->btypes[j0] = BDRY_MASK::INTERIOR | btype;
			e1->btypes[j1] = BDRY_MASK::INTERIOR | btype;
	    
			e0->neigh[j0] = ielem1;
			e1->neigh[j1] = ielem0;
	    
			// to vert
			int v0 = e0->verts[GetEdgeVert(j0, 0)];  // fixeme: use v0_
			int v1 = e0->verts[GetEdgeVert(j0, 1)];
			types_vert[v0] |= btype;
			types_vert[v1] |= btype;
		}
		else if (ielem0 >= 0) {
			// Boundary edge
			Elem *e0 = &elems[ielem0];
			int j0;
	    
			for (j0 = 0; j0 < NEdge; j0++)
				if (e0->edges[j0] == iedge)
					break;
			assert(j0 < NEdge);

			e0->btypes[j0] = btype;
			e0->neigh[j0] = -1;

			// to vert
			int v0 = e0->verts[GetEdgeVert(j0, 0)];  // fixeme: use v0_
			int v1 = e0->verts[GetEdgeVert(j0, 1)];
			types_vert[v0] |= btype;
			types_vert[v1] |= btype;
		}
		else {
			/* phgAbort(0);  // edge with no element */
            abort();
		}
    }  

    // // Dump neighs
    // for (unsigned int ielem = 0; ielem < nelem; ielem++) {
    // 	Elem *e = &elems[ielem];
    // 	phgInfo(2, "e neigh %3d: %3d %3d %3d\n", ielem,
    // 	       e->neigh[0], e->neigh[1], e->neigh[2]
    // 	       );
    // }

    // ------------------------ not tested ----------------------------
	/* partition(part_file_name); */
    /* get_local_grid(); */

    /* init_gref(); */
    /* init_geom(); */

    /* statistic(); */

    return;
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

    // -------------------------------------- read_mesh ------------------------------------------
    /* read_mesh("square_quad.msh", NULL); */

    // ------------------------------------ dFdx and Jac -----------------------------------------
    /* Elem e; */
    /* double lambda[2] = {0,0}; */
    /* double x,y,z; */
    /* updateJacobian(e, lambda); */
    /* cout << e.dFdx << endl; */
    /* cout << e.Jac << endl; */
    /* lambda2xyzDirect(e, 0.5, 0.5, x, y, z); */
    /* cout << "x = " << x << endl; */
    /* cout << "y = " << y << endl; */
    /* cout << "z = " << z << endl; */
    /* double xi, eta; */
    /* xyz2lambdaDirect(e, 2, 2, 4, xi, eta); */
    /* cout << "xi = " << xi << endl; */
    /* cout << "eta = " << eta << endl; */

    // ----------------------------------- const char * ------------------------------------------
    // const 是指不能通过这个指针修改字符串内容，但是可以通过别的指向这个内存的指针修改
    /* char str[10] = "hello"; */
    /* const char *pt = str; */
    /* cout << str << endl; */
    /* cout << pt << endl; */
    /* cout << *pt << endl; */
    /* /1* *pt = "world"; *1/ */
    /* strcpy(str, "world"); */ 
    /* cout << str << endl; */
    /* cout << *pt << endl; */
    
    // ------------------------------- Operation of RealMat --------------------------------------
    /* RealMat A(3,3); */
    /* A << 1,1,1,1,1,1,1,1,1; */
    /* A = A * 4; */
    /* cout << A << endl; */

    // ---------------------------------------- auto ---------------------------------------------
    // auto 自动类型推导，必须进行初始化
    /* auto x = 1; */
    /* auto y = 1.5; */
    /* cout << "bytes of x = " << sizeof(x) << " x = " <<  x << endl; */
    /* cout << "bytes of y = " << sizeof(y) << " y = " <<  y << endl; */

    // -------------------------------------- std::vector ----------------------------------------
    // std::vector 是 C++ 中最常用的动态数组容器，可以自动管理内存，支持动态扩容、随机访问和多种高效操作。以下是 详细用法指南，包含基础操作、高级技巧和性能优化建议。    
    // Initialization
    /* std::vector<int> vec1; // 空向量 */
    /* std::vector<int> vec2(5); // 指定初始大小和默认值（5个0） */
    /* std::vector<int> vec3 = {1, 2, 3, 4, 5}; // 初始化列表 */
    /* std::vector<int> vec4(vec3); // 复制另一个向量 */
    /* std::vector<int> vec5(5, 42); // 指定大小和初始值（5个42） */

    /* // 访问元素 */
    /* std::vector<int> vec = {10, 20, 30}; */
    /* int a = vec[0];  // 10 下标访问（不检查越界） */
    /* int b = vec.at(1); // 20 at() 访问（越界抛出异常） */
    /* int first = vec.front(); // 10 首元素 */
    /* int last = vec.back();   // 30 末元素 */
    /* for (auto it = vec.begin(); it != vec.end(); ++it) { // 迭代器访问 */
    /*     std::cout << *it << std::endl; */
    /* } */

    /* // 添加元素 */
    /* vec.push_back(1); // {10, 20, 30, 1} */
    /* vec.push_back(2); // {10, 20, 30, 1, 2} 尾部插入 */
    /* vec.insert(vec.begin() + 1, 99); // {10, 99, 20, 30, 1, 2} 插入到指定位置 */
    /* vec.insert(vec.end(), {3, 4, 5}); // {10, 99, 20, 30, 1, 2, 3, 4, 5} 批量插入 */
    /* for (auto it = vec.begin(); it != vec.end(); ++it) { */
    /*     std::cout << *it << std::endl; */
    /* } */

    /* // 删除元素 */
    /* vec.pop_back(); // {10, 99, 20, 30, 1, 2, 3, 4} 删除尾部元素 */
    /* vec.erase(vec.begin() + 1); // {10, 20, 30, 1, 2, 3, 4} 删除指定位置元素 */
    /* vec.clear(); // {} 删除所有元素 */
    /* vec.resize(10); */
    /* for (auto it = vec.begin(); it != vec.end(); ++it) { */
    /*     std::cout << *it << std::endl; */
    /* } */

    // --------------------------------------- phgInfo -------------------------------------------
    /* FILE* log_file = fopen("app.log", "w"); */
    /* phgInfo(1, "System initialized. Version: %s\n", "1.0.0"); */
    /* phgInfo(2, "Debug info: %d, %f\n", 42, 3.14); */

    // ---------------------------------------- %04d ---------------------------------------------
    /* char s[200]; */
    /* char phgLogFilename[] = "log"; */
    /* for (int i = 0; i < 10; i++) { */
    /*     sprintf(s, "%s.%04d" "", phgLogFilename, (int)i); */
    /*     cout << s << endl; */
    /* } */
    // %04d 中的 4 表示最小宽度为 4, 0 表示不足 4 位时用 0 填充。例如：
    // •5 → 0005
    // •123 → 0123
    // •4567 → 4567

    // ----------------------------------------- sprintf -----------------------------------------
    /* char buffer[100];  // 确保足够大 */
    /* int num = 42; */
    /* double pi = 3.14159; */
    /* sprintf(buffer, "Number: %d, Pi: %.2f", num, pi); // 在 C++ 中，sprintf 是一个 C 标准库函数，用于将格式化数据写入字符串缓冲区。它的功能类似于 printf，但不是输出到标准输出（stdout），而是写入指定的字符数组（C 风格字符串） */
    /* std::cout << buffer << std::endl;  // 输出: "Number: 42, Pi: 3.14" */

    // -------------------------------------- PATH_MAX -------------------------------------------
    /* cout << PATH_MAX << endl; // PATH_MAX 是一个宏，用于表示系统支持的文件路径的最大长度(以字节为单位) = 4096 in Linux */

    // ----------------------------------------- sin ---------------------------------------------
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
