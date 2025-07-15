#include "mesh.hpp"

namespace REF_ELEM {
	int EdgeVertTria[3][3] = {
		{0, 1, 2}, 
		{1, 2, 0}, 
		{2, 0, 1}
	};
	int EdgeVertQuad[4][3] = {
		{0, 1, 2}, 
		{1, 2, 3}, 
		{2, 3, 0},
        {3, 0, 1}
	};
}

#define GetEdgeVertTria(i,j) REF_ELEM::EdgeVertTria[i][j]
#define GetEdgeVertQuad(i,j) REF_ELEM::EdgeVertQuad[i][j]

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



static Btype default_bc_map(int mark)
{
    return BDRY_MASK::DIRICHLET;
}


int
compare_bdry_edge(const void *p0, const void *p1)
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

    abort();
}



void
Grid::read_mesh(const char *mesh_file_name, BC_FUNCTION user_bc_map, const char *part_file_name)
{
    int nbdry_edge;
    BdryEdge *bdry_edges;

    if (user_bc_map == NULL)
		bc_map = default_bc_map;
    else
		bc_map = user_bc_map;


	cout << "---------- Reading Mesh: " << mesh_file_name << " ----------" << endl;

#define READ_NUMBER													\
    if (!get_token(fp, token)) strcpy(token, "End");				\
    if (isalpha((int)(token[0]))) {									\
		fprintf(stderr, "fewer entries (%u) than expected.\n", i);	\
		break;														\
    }
#undef ERROR
#define ERROR {n = __LINE__; goto error;}
    
    FILE *fp ;
    char token[100]; 
    int n;
    if ((fp = fopen(mesh_file_name, "r")) == NULL) {
		printf("can't open mesh file \"%s\"!\n", mesh_file_name);
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
		// next_token:
		if (!strcasecmp(token, "$Nodes")) {
			if (!get_token(fp, token))
				ERROR;
			n = atoi(token);
			fprintf(stderr, "number of vertices: %d\n", n);
			verts = new Coord[n](); // verts is a member viriable of class Grid:
                                    // Coord *verts;

			int i;
			for (i = 0; i < n; i++) {
				READ_NUMBER;	// Unused
				READ_NUMBER;
				verts[i][0] = atof(token);
				if (!get_token(fp, token))
					ERROR;
				verts[i][1] = atof(token);
				if (!get_token(fp, token))
					ERROR;
				verts[i][2] = atof(token); 
			}
			assert(i == n);
			nvert = n;
            /* cout << verts[0][0] << "\t" << verts[0][1] << "\t" << verts[0][2] << endl; */
            /* cout << verts[1][0] << "\t" << verts[1][1] << "\t" << verts[1][2] << endl; */

			if (!get_token(fp, token))
				ERROR;
			assert(!strcasecmp(token, "$EndNodes"));
		}
		else if (!strcasecmp(token, "$Elements")) {
			if (!get_token(fp, token))
				ERROR;
			n = atoi(token);
            fprintf(stderr, "number of elements(edges + triangles + quadrilaterals): %d\n", n);

			bdry_edges = new BdryEdge[n]();
			elems = new Elem[n]();
			/* tria_elems = new TriaElem[n](); */
			/* quad_elems = new QuadElem[n](); */

			int nedge_read = 0, nread = 0, ntria_read = 0, nquad_read = 0;
			int elem_type_msh = 0;

			int i;
			for (i = 0; i < n; i++) {
				int phy_id = 0, entity_id = 0;
				READ_NUMBER;	// id 
				assert(atoi(token) > 0);

				READ_NUMBER;	// # element type
				elem_type_msh = atoi(token);
				// edge or triangle or quadrilateral
				assert(elem_type_msh == 1 || elem_type_msh == 2 || elem_type_msh == 3);

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

				if (elem_type_msh == 1) {
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
				} else if (elem_type_msh == 2) {
					//
					// triangle
					// 
					READ_NUMBER;
					/* TriaElem *e = &tria_elems[ntria_read]; */
					Elem *e = &elems[nread];
                    e->elem_type = 'T';
					e->verts[0] = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					e->verts[1] = atoi(token) - 1;
					if (!get_token(fp, token))
						ERROR;
					e->verts[2] = atoi(token) - 1;
					e->index = nread;

                    nread++;
					ntria_read++;
				}
                else if (elem_type_msh == 3) {
					//
					// quadrilateral
					// 
					READ_NUMBER;
					/* QuadElem *e = &quad_elems[nquad_read]; */
					Elem *e = &elems[nread];
                    e->elem_type = 'Q';
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
					e->index = nread; // tria and quad's index are seperate
                    
                    nread++;
					nquad_read++;
				}
				else {
					ERROR;
				}
			}
            fprintf(stderr, "number of edges in #Elements: %d\n", nedge_read);
            fprintf(stderr, "number of triangles in #Elements: %d\n", ntria_read);
            fprintf(stderr, "number of quadrilaterals in #Elements: %d\n", nquad_read);
            ntria_elem = ntria_read;
            nquad_elem = nquad_read;
			nelem = nread;
            /* cout << tria_elems[0].verts[0] << "\t" << tria_elems[0].verts[1] << "\t" << tria_elems[0].verts[2] << endl; */
            /* cout << quad_elems[0].verts[0] << "\t" << quad_elems[0].verts[1] << "\t" << quad_elems[0].verts[2] << "\t" << quad_elems[0].verts[3] << endl; */
            
            /* cout << nelem << endl; */
			nbdry_edge = nedge_read;
			assert(i == n);
			// if (i < n)
			// 	goto next_token;
	    
			if (!get_token(fp, token))
				ERROR;
			assert(!strcasecmp(token, "$EndElements"));
		}
    }

    qsort(bdry_edges, nbdry_edge, sizeof(BdryEdge), compare_bdry_edge);




    // ------------------------------------------------------------
    //
    // Assign edge indices (above is bdry_edges, which is not
    //                      one of Grid's member)
    // 
    // ------------------------------------------------------------
    /* cout << endl << nelem << endl; */
    /* cout << endl << ntria_elem << endl << nquad_elem << endl; */
#if 0
    for (int i = 0; i < nelem; i++) {
        Elem e = elems[i];
        cout << e.elem_type << endl;
    }
#endif

    edges = new Edge[ntria_elem * 3 + nquad_elem * 4]();
    int iter = 0;
    for (unsigned int i = 0; i < nelem; i++) {
        /* if (i == 0) { */
        /*     cout << tria_elems[i].verts[0] << "\t" << tria_elems[i].verts[1] << "\t" << tria_elems[i].verts[2] << endl; */
        /* } */
        int v0, v1;
        
        Elem e = elems[i]; 
        if (e.elem_type == 'T') {
            for (unsigned int j = 0; j < 3; j++) {
                v0 = e.verts[GetEdgeVertTria(j, 0)];
                v1 = e.verts[GetEdgeVertTria(j, 1)];
#if 0
                if (i == 8) {
                    cout << endl << "v0 = " << v0 << endl;
                    cout << endl << "v1 = " << v1 << endl;
                }
#endif

                SortIndex(v0, v1);

                edges[iter][0] = v0;
                edges[iter][1] = v1;
                iter++;
            }

            // Ordering
            {
                int V0, V1, V2, v0, v1, v2;
                V0 = e.verts[0];
                V1 = e.verts[1];
                V2 = e.verts[2];

                v0 = v1 = v2 = 0;
                (V0 > V1) ? v0++ : v1++;
                (V0 > V2) ? v0++ : v2++;
                (V1 > V2) ? v1++ : v2++;

                e.ordering[0] = v0;
                e.ordering[1] = v1;
                e.ordering[2] = v2;
            }
        } else if (e.elem_type == 'Q') {
            for (unsigned int j = 0; j < 4; j++) {
                v0 = e.verts[GetEdgeVertQuad(j, 0)];
                v1 = e.verts[GetEdgeVertQuad(j, 1)];
#if 0
                if (i == 0) {
                    cout << "v0 = " << v0 << endl;
                    cout << "v1 = " << v1 << endl;
                }
#endif
                SortIndex(v0, v1);
                edges[iter][0] = v0;
                edges[iter][1] = v1;
                iter++;
            }

            // Ordering
            {
                int V0, V1, V2, V3, v0, v1, v2, v3;
                V0 = e.verts[0];
                V1 = e.verts[1];
                V2 = e.verts[2];
                V3 = e.verts[3];

                v0 = v1 = v2 = v3 = 0;
                (V0 > V1) ? v0++ : v1++;
                (V0 > V2) ? v0++ : v2++;
                (V0 > V3) ? v0++ : v3++;
                (V1 > V2) ? v1++ : v2++;
                (V1 > V3) ? v1++ : v3++;
                (V2 > V3) ? v2++ : v3++;

                e.ordering[0] = v0;
                e.ordering[1] = v1;
                e.ordering[2] = v2;
                e.ordering[3] = v3;
            }
        } // end if
    }
    /* cout << "iter = " << iter << endl; */ 
    /* cout << "num of edge = " << ntria_elem*3 + nquad_elem*4 << endl; */ 

    qsort(edges, ntria_elem*3+nquad_elem*4, sizeof(Edge), compare_edge);

    int count = 0; 
    {
		int i0 = 0;
		for (unsigned int i = i0+1; i < ntria_elem*3+nquad_elem*4; i++) {
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
    cout << "number of edges(iter elems): " << nedge << endl;


    // ------------------------------------------------------------
    //
    // Boudnary makers
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
        if (e->elem_type == 'T') {
            for (unsigned int j = 0; j < 3; j++) {
                int v0 = e->verts[GetEdgeVertTria(j, 0)];
                int v1 = e->verts[GetEdgeVertTria(j, 1)];

                SortIndex(v0, v1);

                Edge edge = {v0, v1}, *p;

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
        } else if (e->elem_type == 'Q') {
            for (unsigned int j = 0; j < 4; j++) {
                int v0 = e->verts[GetEdgeVertQuad(j, 0)];
                int v1 = e->verts[GetEdgeVertQuad(j, 1)];

                SortIndex(v0, v1);
                Edge edge = {v0, v1}, *p;

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
		// phgInfo(2, "search: %d %d\n", edge[0], edge[1]);

		p = (BdryEdge *) bsearch(&edge, bdry_edges,
								 nbdry_edge, sizeof(BdryEdge),
								 compare_bdry_edge);
		if (p != NULL) {
			mark = (*p)[2];
		}

		Btype btype = (mark >= 0) ? bc_map(mark) : 0;
		types_edge[iedge] = btype;


		//
		// Update elem edges info
		// 
		int ielem0 = edge2elem[iedge][0];
		int ielem1 = edge2elem[iedge][1];

        if (ielem0 >= 0 && ielem1 >= 0) {
            // Interior edge
            int j0, j1, v0, v1;
            Elem *e0 = &elems[ielem0];
            Elem *e1 = &elems[ielem1];

            if (e0->elem_type == 'T') { // tria elem
                for (j0 = 0; j0 < 3; j0++)
                    if (e0->edges[j0] == iedge)
                        break;
                assert(j0 < 3);
                v0 = e0->verts[GetEdgeVertTria(j0, 0)];
                v1 = e0->verts[GetEdgeVertTria(j0, 1)];
                e0->btypes[j0] = BDRY_MASK::INTERIOR | btype;
                e0->neigh[j0] = ielem1;
            } else if (e0->elem_type == 'Q'){ // quad elem
                for (j0 = 0; j0 < 4; j0++)
                    if (e0->edges[j0] == iedge)
                        break;
                assert(j0 < 4);
                v0 = e0->verts[GetEdgeVertQuad(j0, 0)];
                v1 = e0->verts[GetEdgeVertQuad(j0, 1)];
                e0->btypes[j0] = BDRY_MASK::INTERIOR | btype;
                e0->neigh[j0] = ielem1;
            }

            if (e1->elem_type == 'T') { // tria elem
                for (j1 = 0; j1 < 3; j1++)
                    if (e1->edges[j1] == iedge)
                        break;
                assert(j1 < 3);
                e1->btypes[j1] = BDRY_MASK::INTERIOR | btype;
                e1->neigh[j1] = ielem0;
            } else if (e1->elem_type == 'Q') { // quad elem
                for (j1 = 0; j1 < 4; j1++)
                    if (e1->edges[j1] == iedge)
                        break;
                assert(j1 < 4);
                e1->btypes[j1] = BDRY_MASK::INTERIOR | btype;
                e1->neigh[j1] = ielem0;
            }


            types_edge[iedge] |= BDRY_MASK::INTERIOR;	    

            types_vert[v0] |= btype;
            types_vert[v1] |= btype;
        } else if (ielem0 >= 0) {
            // Boundary edge
            int j0, v0, v1;

            Elem *e0 = &elems[ielem0];
            if (e0->elem_type == 'T') {
                for (j0 = 0; j0 < 3; j0++)
                    if (e0->edges[j0] == iedge)
                        break;
                assert(j0 < 3);
                v0 = e0->verts[GetEdgeVertTria(j0, 0)];
                v1 = e0->verts[GetEdgeVertTria(j0, 1)];
                e0->btypes[j0] = btype;
                e0->neigh[j0] = -1;
            } else if (e0->elem_type == 'Q') {
                for (j0 = 0; j0 < 4; j0++)
                    if (e0->edges[j0] == iedge)
                        break;
                assert(j0 < 4);
                v0 = e0->verts[GetEdgeVertQuad(j0, 0)];
                v1 = e0->verts[GetEdgeVertQuad(j0, 1)];
                e0->btypes[j0] = btype;
                e0->neigh[j0] = -1;
            }

            types_vert[v0] |= btype;
            types_vert[v1] |= btype;
        }
        else {
            abort();
        }
    }

		    
	cout << "-------------------- Read Mesh Done ---------------------" << endl;
    return;
}



/* bool */
/* Grid::getFaceLambda(Elem &e, int iside, const Real *p, Real *lambda, bool do_swap) const */
/* { */
/* 	bool swaped = false; */
		
/*     int v0 = GetEdgeVert(iside, 0); */
/*     int v1 = GetEdgeVert(iside, 1); */
/*     int v2 = GetEdgeVert(iside, 2);	 // op vert */

/*     int V0 = e.verts[v0]; */
/*     int V1 = e.verts[v1]; */
/*     int V2 = e.verts[v2];  // op vert */

/*     // Coresponde to above EdgeVert */
/*     lambda[0] = lambda[1] = lambda[2] = 0.; */

/*     Real pp[2] = {p[0], p[1]}; */
/*     if (do_swap && V0 > V1) { */
/* 		// swap */
/* 		// order changed so that when computing flux, the quad points matches */
/* 		pp[0] = p[1]; */
/* 		pp[1] = p[0]; */
/* 		swaped = true; */
/*     } */

    
/*     if (v0 >= 1) lambda[v0 - 1] = pp[0]; //*(p++); */
/*     if (v1 >= 1) lambda[v1 - 1] = pp[1]; //*(p++); */

/*     //phgInfo(2, "%f %f %f\n", lambda[0], lambda[1], lambda[2]); */

/* 	/* */
/*      * */     
/*      *     2 */ 
/*      *     |\ */ 
/*      *     |  \ */ 
/*      *     |    \ */ 
/*      * */    
/*      *     |        \ */ 
/*      *     |          \ */ 
/*      *     |            \ */ 
/*      *     0 ---    ---- 1 */
/*      * */
/*      * v0 >= 1 ---> v0 = {1, 2} */    
/*      * v1 >= 1 ---> v1 = {1, 2} */
/*      * */
/*      * Order is not important */
/*      * */ 
/* 	 * case 1: v0 = 1, v1 = 0, (ξ↑, 0) */
/* 	 * case 1: v0 = 2, v1 = 0, (0, η↑)  -> */
/* 	 * case 2: v0 = 0, v1 = 1, (ξ↑, 0)  -> */
/* 	 * case 2: v0 = 0, v1 = 2, (0, η↑) */  
/* 	 * case 3: v0 = 1, v1 = 2, (ξ↑, η↓)  -> */
/* 	 * case 4: v0 = 2, v1 = 1, (ξ↓, η↑) */
/* 	 * This is indp to swap */
/*      * */    
/* 	 *1/ */
    
/*     return swaped; */
/* } */


/* void */
/* Grid::check_conformal2() */
/* // */ 
/* // Check edge to elem and vert to elem */
/* // */ 
/* { */

/*     for (unsigned int ielem = 0; ielem < nelem; ielem++) { */
/* 		Elem &e = elems[ielem]; */

/* 		// */
/* 		// check edge_to_elem */
/* 		// */
/* 		for (int i = 0; i < NEdge; i++) { */
/* 			int iedge = e.edges[i]; */

/* 			int ii; */
/* 			for (ii = 0; ii < edge_to_elem[iedge].size(); ii++) */ 
/* 				if (edge_to_elem[iedge][ii].eind == e.index) */
/* 					break; */
/* 			assert(ii < edge_to_elem[iedge].size()); */
/* 			assert(edge_to_elem[iedge][ii].gind == i); */
/* 		} */

/* 		// */
/* 		// check vert_to_elem */
/* 		// */ 
/* 		for (int i = 0; i < NVert; i++) { */
/* 			int ivert = e.verts[i]; */

/* 			int ii; */
/* 			for (ii = 0; ii < vert_to_elem[ivert].size(); ii++) */ 
/* 				if (vert_to_elem[ivert][ii].eind == e.index) */
/* 					break; */
/* 			assert(ii < vert_to_elem[ivert].size()); */
/* 			assert(vert_to_elem[ivert][ii].gind == i); */
/* 		} */
/*     } */    

/*     // TODO: add btype check */
    
/*     return; */
/* } */




/* void */
/* Grid::output_gmsh(const char *prefix, bool second_order) */
/* { */
/*     const char ofs = 1; */
/*     char fname[1000]; */
/*     FILE *fp = NULL; */

/*     if (phgNProcs == 1) { */
/* 		sprintf(fname, "%s.msh", prefix); */
/*     } */
/*     else { */
/* 		sprintf(fname, "%s.p%02d.msh", prefix, phgRank); */
/*     } */
/*     fp = fopen(fname, "w"); */
/*     assert(fp != NULL); */
    
/*     fprintf(fp, "$MeshFormat\n"); */
/*     fprintf(fp, "2.2 0 8\n"); */
/*     fprintf(fp, "$EndMeshFormat\n"); */

/*     if (second_order) { */
/* 		// */
/* 		// */
/* 		// Second order mesh */
/* 		// */
/* 		// */ 
/* 		Coord *edge_verts = new Coord[nedge]{}; */
/* 		for (unsigned int ielem = 0; ielem < nelem; ielem++) { */
/* 			Elem &e = elems[ielem]; */

/* 			// Real x0, y0; */
/* 			// lambda2xyz(e, 1./3., 1./3., x0, y0); */

/* 			for (int i = 0; i < NEdge; i++) { */
/* 				int iedge = e.edges[i]; */
/* 				int ivert0 = e.verts[GetEdgeVert(i, 0)]; */
/* 				int ivert1 = e.verts[GetEdgeVert(i, 1)]; */
/* 				Real x = 0, y = 0; */

/* 				x = .5 * (verts[ivert0][X_DIR] + verts[ivert1][X_DIR]); */
/* 				y = .5 * (verts[ivert0][Y_DIR] + verts[ivert1][Y_DIR]); */
/* 				// if (_dbg_shrink_) { */
/* 				//     x = (x - x0) * 0.85 + x0; */
/* 				//     y = (y - y0) * 0.85 + y0; */
/* 				// } */
/* 				edge_verts[iedge][X_DIR] = x; */
/* 				edge_verts[iedge][Y_DIR] = y; */
/* 			} */
/* 		} */
	
	
/* 		fprintf(fp, "$Nodes\n"); */
/* 		fprintf(fp, "%d\n", nvert + nedge); */
/* 		for (int i = 0; i < nvert; i++) { */
/* 			fprintf(fp, "%d %30.18e %30.18e %30.18e\n", */ 
/* 					i + ofs, verts[i][0], verts[i][1], verts[i][2]); */
/* 		} */
/* 		for (int i = 0; i < nedge; i++) { */
/* 			fprintf(fp, "%d %30.18e %30.18e %30.18e\n", */ 
/* 					nvert + i + ofs, edge_verts[i][0], edge_verts[i][1], edge_verts[i][2]); */
/* 		} */
/* 		fprintf(fp, "$EndNodes\n"); */

/* 		fprintf(fp, "$Elements\n"); */
/* 		fprintf(fp, "%d\n", nvert + nedge + nelem ); */

/* 		// Vert */
/* 		for (int i = 0; i < nvert; i++) { */
/* 			fprintf(fp, "%d 15 2 1 %d  %d\n", */ 
/* 					i + ofs, (int) types_vert[i], */
/* 					i + ofs); */
/* 		} */

/* 		// // Edge as node */
/* 		// for (int i = 0; i < nedge; i++) { */
/* 		//     fprintf(fp, "%d 15 2 1 %d  %d\n", */ 
/* 		// 	    i + nvert + ofs, (int) types_edge[i], */
/* 		// 	    nvert + i + ofs); */
/* 		// } */

/* 		for (int i = 0; i < nedge; i++) { */
/* 			fprintf(fp, "%d 8 2 1 %d  %d %d %d\n", */ 
/* 					i + nvert + ofs, (int) types_edge[i], */
/* 					edges[i][0] + ofs, edges[i][1] + ofs, */
/* 					nvert + i + ofs); */
/* 		} */

/* 		// Elem */
/* 		for (int i = 0; i < nelem; i++) { */
/* 			fprintf(fp, "%d 9 2 1 %d  %d %d %d %d %d %d\n", */ 
/* 					i + nvert + nedge + ofs, (int) types_elem[i], */
/* 					elems[i].verts[0] + ofs, */
/* 					elems[i].verts[1] + ofs, */
/* 					elems[i].verts[2] + ofs, */
/* 					nvert + elems[i].edges[0] + ofs, */
/* 					nvert + elems[i].edges[1] + ofs, */
/* 					nvert + elems[i].edges[2] + ofs */
/* 					); */
/* 		} */
/* 		fprintf(fp,  "$EndElements\n"); */

/* 		fprintf(fp, "$Comments\n"); */
/* 		fprintf(fp, "0 %d %d %d\n", nvert, nvert + nedge, nvert + nedge + nelem); */
/* 		fprintf(fp, "$EndComments\n"); */
/*     } */
/*     else { */
/* 		// */
/* 		// */
/* 		// Linear mesh */
/* 		// */
/* 		// */ 
/* 		fprintf(fp, "$Nodes\n"); */
/* 		fprintf(fp, "%d\n", nvert); */
/* 		for (int i = 0; i < nvert; i++) { */
/* 			fprintf(fp, "%d %30.18e %30.18e %30.18e\n", */ 
/* 					i + ofs, verts[i][0], verts[i][1], 0.); //verts[i][2]); */
/* 		} */
/* 		fprintf(fp, "$EndNodes\n"); */


/* 		fprintf(fp, "$Elements\n"); */
/* 		fprintf(fp, "%d\n", nvert + nedge + nelem ); */

/* 		// Vert */
/* 		for (int i = 0; i < nvert; i++) { */
/* 			fprintf(fp, "%d 15 2 1 %d  %d\n", */ 
/* 					i + ofs, (int) types_vert[i], */
/* 					i + ofs); */
/* 		} */

/* 		// Edge */
/* 		for (int i = 0; i < nedge; i++) { */
/* 			fprintf(fp, "%d 1 2 1 %d  %d %d\n", */ 
/* 					i + nvert + ofs, (int) types_edge[i], */
/* 					edges[i][0] + ofs, edges[i][1] + ofs); */
/* 		} */

/* 		// Elem */
/* 		for (int i = 0; i < nelem; i++) { */
/* 			fprintf(fp, "%d 2 2 1 %d  %d %d %d\n", */ 
/* 					i + nvert + nedge + ofs, (int) types_elem[i], */
/* 					elems[i].verts[0] + ofs, */ 
/* 					elems[i].verts[1] + ofs, */ 
/* 					elems[i].verts[2] + ofs); */
/* 		} */
/* 		fprintf(fp,  "$EndElements\n"); */

/* 		fprintf(fp, "$Comments\n"); */
/* 		fprintf(fp, "0 %d %d %d\n", nvert, nvert + nedge, nvert + nedge + nelem); */
/* 		fprintf(fp, "$EndComments\n"); */
/*     } */

/*     fclose(fp); */
/*     return; */
/* } */


int
Grid::compare_edge(const void *p0, const void *p1)
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





/* void */
/* Grid::check_conformal() */
/* // */ 
/* // Check element's edges and verts */
/* // */ 
/* { */

/*     for (unsigned int ielem = 0; ielem < nelem; ielem++) { */
/* 		Elem &e = elems[ielem]; */

/* 		for (int i = 0; i < NEdge; i++) { */
/* 			int iedge = e.edges[i]; */
/* 			int ivert0 = e.verts[GetEdgeVert(i, 0)]; */
/* 			int ivert1 = e.verts[GetEdgeVert(i, 1)]; */

/* 			if (ivert0 > ivert1) { */
/* 				assert(edges[iedge][0] == ivert1); */
/* 				assert(edges[iedge][1] == ivert0); */
/* 			} */
/* 			else { */
/* 				assert(edges[iedge][0] == ivert0); */
/* 				assert(edges[iedge][1] == ivert1); */
/* 			} */
/* 		} */
/*     } */    

/*     // TODO: add btype check */
    
/*     return; */
/* } */



/* void */
/* Grid::init_geom() */
/* { */
/*     // build isop map */
/*     // if (isop_map) { */
/*     phgInfo(0, "Geom init\n"); */
/*     isop_map = NULL;  // TODO: free dof */
/*     // } */
/*     // else { */
/*     // 	phgInfo(0, "Geom reinit\n"); */
/*     // } */

/*     FESpace *isop_space = new FESpace(this, FEM::getLagrangePn(1), SDim, "isop space P1"); */
/*     isop_map = new Dof(isop_space, "isop dof"); */
/* 	if (Params->use_proj_sphere) { */
/* 		assert(Params->proj_sphere_radius > 0.); */
/* 		proj_sphere_radius = Params->proj_sphere_radius; */
/* 		phgInfo(0, "Grid: proj sphere radius %30.15f\n", Params->proj_sphere_radius); */
/* 	} */
    
/*     Scalar *data = isop_map->coefData(); */
/*     for (unsigned int ivert = 0; ivert < nvert; ivert++) { */
/* 		*(data++) = verts[ivert][X_DIR]; */
/* 		*(data++) = verts[ivert][Y_DIR]; */
/* 		*(data++) = verts[ivert][Z_DIR]; */
/*     } */

/*     //isop_map->dump(); */

/*     for (unsigned int ielem = 0; ielem < nelem; ielem++) { */
/* 		Elem *e = &elems[ielem]; */


/* #if 0 */
/* 		// jacobian */
/* 		// */
/* 		/1* computes the Jabobian */
/* 		 *	J = D\lambda / Dx = [ D\lambda_0/Dx, D\lambda_0/Dy, D\lambda_0/Dz, c0; */
/* 		 *			      D\lambda_1/Dx, D\lambda_1/Dy, D\lambda_1/Dz, c1; */
/* 		 *			      D\lambda_2/Dx, D\lambda_2/Dy, D\lambda_2/Dz, c2; */
/* 		 *			      D\lambda_3/Dx, D\lambda_3/Dy, D\lambda_3/Dz, c3 ] */
/* 		 * of the element, by inverting the matrix */
/* 		 * */
/* 		 *	[x0, x1, x2, x3; y0, y1, y2, y3; z0, z1, z2, z3; 1 1 1 1]. */
/* 		 * */
/* 		 * Note: lambda = J[0..Dim][0..Dim-1] * x + J[0..Dim][Dim] */
/* 		 *1/ */
/* 		{ */
/* 			Real a[3][3]; */
/* 			for (unsigned int i = 0; i < 2; i++) */
/* 				for (unsigned int j = 0; j < 3; j++) */
/* 					a[i][j] = verts[e->verts[j]][i]; // vert[e->verts[j]][i] */

/* 			for (unsigned int j = 0; j < 3; j++) */
/* 				a[2][j] = 1.; */

/* 			Real JJ[3][3] = { */
/* 				{a[1][1]*a[2][2] - a[1][2]*a[2][1], a[0][2]*a[2][1] - a[0][1]*a[2][2], a[0][1]*a[1][2] - a[0][2]*a[1][1]}, */
/* 				{a[1][2]*a[2][0] - a[1][0]*a[2][2], a[0][0]*a[2][2] - a[0][2]*a[2][0], a[0][2]*a[1][0] - a[0][0]*a[1][2]}, */
/* 				{a[1][0]*a[2][1] - a[1][1]*a[2][0], a[0][1]*a[2][0] - a[0][0]*a[2][1], a[0][0]*a[1][1] - a[0][1]*a[1][0]}}; */

/* 			Real det = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] */
/* 				+ a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0]; */

/* 			for (unsigned int i = 0; i < 3; i++) */
/* 				for (unsigned int j = 0; j < 3; j++) */
/* 					e->Jac[i][j] = JJ[i][j] / det; */
/* 		} */
/* #else */
/* #endif */


/* 		// area */
/* 		Real nf[3];				   // face normal */
/* 		{ */
/* 			const Real *p0 = verts[e->verts[0]]; */
/* 			const Real *p1 = verts[e->verts[1]]; */
/* 			const Real *p2 = verts[e->verts[2]]; */

/* 			Real d, a[3], b[3]; */
/* 			a[0] = p1[0] - p0[0]; a[1] = p1[1] - p0[1]; a[2] = p1[2] - p0[2]; */
/* 			b[0] = p2[0] - p0[0]; b[1] = p2[1] - p0[1]; b[2] = p2[2] - p0[2]; */

/* 			CROSS_PRODUCT(a, b, nf); */
/* 			NORMALIZE(nf); */
/* 			e->face_normal[X_DIR] = nf[X_DIR]; */
/* 			e->face_normal[Y_DIR] = nf[Y_DIR]; */
/* 			e->face_normal[Z_DIR] = nf[Z_DIR]; */

/* 			/1* compute the area (Heron's formula) and diameter */
/* 			   (http://mathworld.wolfram.com/Circumradius.html) *1/ */
/* 			a[0] = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); */
/* 			a[1] = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]); */
/* 			b[0] = p2[0] - p1[0]; b[1] = p2[1] - p1[1]; b[2] = p2[2] - p1[2]; */
/* 			a[2] = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]); */
/* 			d = (a[0] + a[1] + a[2]) * 0.5; */

/* 			Real c = d * (d - a[0]) * (d - a[1]) * (d - a[2]); */
/* 			assert(c > 0.); */
/* 			e->area = sqrt(c); */
/* 		} */



/* 		// edge length */
/* 		//   TODO: use isop */
/* 		Coord center; */
/* 		center[X_DIR] = (verts[e->verts[0]][X_DIR] */
/* 						 + verts[e->verts[1]][X_DIR] */
/* 						 + verts[e->verts[2]][X_DIR]) / 3.; */
/* 		center[Y_DIR] = (verts[e->verts[0]][Y_DIR] */
/* 						 + verts[e->verts[1]][Y_DIR] */
/* 						 + verts[e->verts[2]][Y_DIR]) / 3.; */
/* 		center[Z_DIR] = (verts[e->verts[0]][Z_DIR] */
/* 						 + verts[e->verts[1]][Z_DIR] */
/* 						 + verts[e->verts[2]][Z_DIR]) / 3.; */

/* 		for (unsigned int iside = 0; iside < NEdge; iside++) { */
/* 			Real *p0 = verts[e->verts[GetEdgeVert(iside, 0)]]; */
/* 			Real *p1 = verts[e->verts[GetEdgeVert(iside, 1)]]; */
/* 			Real dx = p1[X_DIR] - p0[X_DIR]; */
/* 			Real dy = p1[Y_DIR] - p0[Y_DIR]; */
/* 			Real dz = p1[Z_DIR] - p0[Z_DIR]; */
/* 			Real len = sqrt(dx*dx + dy*dy + dz*dz); */

/* #if 0 */			
/* 			Real nx = dy / len; */
/* 			Real ny = - dx / len; */
/* #else */
/* 			Real te[3] = {dx, dy, dz}; // edge tangential */
/* 			Real ne[3];				   // edge normal */

/* 			CROSS_PRODUCT(nf, te, ne); */
/* 			NORMALIZE(ne); */
/* #endif */			


/* 			Real midx = (p0[X_DIR] + p1[X_DIR]) / 2.; */
/* 			Real midy = (p0[Y_DIR] + p1[Y_DIR]) / 2.; */
/* 			Real midz = (p0[Z_DIR] + p1[Z_DIR]) / 2.; */

/* 			Real rx = center[X_DIR] - midx; */
/* 			Real ry = center[Y_DIR] - midy; */
/* 			Real rz = center[Z_DIR] - midz; */

/* 			int sign = 1; */
/* 			if (rx * ne[X_DIR] + ry * ne[Y_DIR] + rz * ne[Z_DIR] > 0) */
/* 				sign = -1; */

/* 			e->edge_length[iside] = len; */
/* 			e->edge_normal[iside][X_DIR] = sign * ne[X_DIR]; */
/* 			e->edge_normal[iside][Y_DIR] = sign * ne[Y_DIR]; */
/* 			e->edge_normal[iside][Z_DIR] = sign * ne[Z_DIR]; */

/* 			// if (e->btypes[iside] & BDRY_MASK::DIRICHLET) { */
/* 			// 	phgInfo(2, "dirich len %f\n", len); */
/* 			// } */
/* 		} */
/*     } */

/*     return; */
/* } */



/* void */
/* Grid::statistic() */
/* { */
/*     Real area = 0; */
/*     Real length[10] = {0};  // all zero */
/* 	Real max_length = 0, min_length = 1e10; */

/*     if (verb < 1) { */
/* 		return; */
/*     } */

/*     const Qrule *qrule = quad2d.GetRule(4); */
	

/*     // */
/*     // Domain info */
/*     // */
/* 	ForElements(this, ielem) { */
/* 		Elem *e = &elems[ielem]; */

/* 		ForQuadPointsInElem(qrule) { */
/* 			updateJacobian(*e, qp); */
/* 			area += (*qw) * e->getArea(); */
/* 			if (e->index == 0) */
/* 				phgInfo(0, "area: %e\n", e->getArea()); */
/* 		} */
			
/* 		// area += e->area; */

/* 		for (unsigned int iside = 0; iside < NEdge; iside++) { */
/* 			Real len = e->edge_length[iside]; */

/* 			if (len < min_length) */
/* 				min_length = len; */
/* 			if (len > max_length) */
/* 				max_length = len; */
			

/* 			if (!(types_edge[e->edges[iside]] & BDRY_MASK::OWNER)) */
/* 				continue; */
			
/* 			// if (edge_to_elem[e->edges[iside]][0].eind != e->index) */
/* 			// 	continue; */
/* 			if (e->edges_sgn[iside] != 1) */
/* 				continue; */
			
/* 			if (e->btypes[iside] & BDRY_MASK::DIRICHLET) { */
/* 				length[0] += len; */
/* 				//phgInfo(3, "dirich len %f\n", len); */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::NEUMANN) { */
/* 				length[1] += len; */
/* 				//phgInfo(3, "neuman len %f\n", len); */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::BDRY_USER1) { */
/* 				length[2] += len; */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::BDRY_USER2) { */
/* 				length[3] += len; */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::BDRY_USER3) { */
/* 				length[4] += len; */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::BDRY_USER4) { */
/* 				length[5] += len; */
/* 			} */
/* 			if (e->btypes[iside] & BDRY_MASK::BDRY_USER5) { */
/* 				length[6] += len; */
/* 			} */
/* 		} */
/*     } */

/*     if (phgNProcs > 1) { */
/* 		Real sum_length[10] = {0}; */
/* 		length[7] = area; */
/* 		MPI_Allreduce(length, sum_length, */
/* 					  sizeof(length)/sizeof(length[0]), */
/* 					  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/* 		memcpy(length, sum_length, sizeof(length)); */
/* 		area = length[7]; */

/* 		Real length_local[2] = {-min_length, max_length}, length_global[2]; */
/* 		MPI_Allreduce(length_local, length_global, 2, */
/* 					  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); */
/* 		min_length = -length_global[0]; */
/* 		max_length = length_global[1]; */
/*     } */

/*     phgPrintf("\n   Domain area: %30.15f\n", area); */
/*     phgPrintf("   Length dirich  : %30.15f\n", length[0]); */
/*     phgPrintf("   Length neumann : %30.15f\n", length[1]); */
/*     phgPrintf("   Length user1   : %30.15f\n", length[2]); */
/*     phgPrintf("   Length user2   : %30.15f\n", length[3]); */
/*     phgPrintf("   Length user3   : %30.15f\n", length[4]); */
/*     phgPrintf("   Length user4   : %30.15f\n", length[5]); */
/*     phgPrintf("   Length user5   : %30.15f\n\n", length[6]); */
/*     phgPrintf("   Max length     : %30.15f\n", max_length); */
/*     phgPrintf("   Min length     : %30.15f\n", min_length); */
/*     phgPrintf("\n\n"); */


/*     // */
/*     // */
/*     //  Mesh info */
/*     // */
/*     // */
/*     phgPrintf("---\nMesh info:\n"); */
/*     phgPrintf("   nvert: %d\n", nvert_global()); */
/*     phgPrintf("   nedge: %d\n", nedge_global()); */
/*     phgPrintf("   nface: %d\n", nelem_global()); */

/*     int nedge_type[10] = {}; */
/*     for (unsigned int i = 0; i < nedge; i++) { */
/* 		Btype btype = types_edge[i]; */

/* 		if (!(btype & BDRY_MASK::OWNER)) */
/* 			continue; */
		
/* 		if (btype & BDRY_MASK::INTERIOR) */
/* 			nedge_type[0]++; */
/* 		if (btype & BDRY_MASK::DIRICHLET) */
/* 			nedge_type[1]++; */
/* 		if (btype & BDRY_MASK::NEUMANN) */
/* 			nedge_type[2]++; */
/* 		if (btype & BDRY_MASK::BDRY_USER1) */
/* 			nedge_type[3]++; */
/* 		if (btype & BDRY_MASK::BDRY_USER2) */
/* 			nedge_type[4]++; */
/* 		if (btype & BDRY_MASK::BDRY_USER3) */
/* 			nedge_type[5]++; */
/* 		if (btype & BDRY_MASK::BDRY_USER4) */
/* 			nedge_type[6]++; */
/* 		if (btype & BDRY_MASK::BDRY_USER5) */
/* 			nedge_type[7]++; */
/*     } */


/*     if (phgNProcs > 1) { */
/* 		int sum_edge_type[10] = {}; */
/* 		MPI_Allreduce(nedge_type, sum_edge_type, */
/* 					  sizeof(nedge_type)/sizeof(nedge_type[0]), */
/* 					  MPI_INT, MPI_SUM, */
/* 					  MPI_COMM_WORLD); */
/* 		memcpy(nedge_type, sum_edge_type, sizeof(nedge_type)); */
/*     } */

/*     phgPrintf("\n   nedge interior: %d\n", nedge_type[0]); */
/*     phgPrintf("   nedge dirich  : %d\n", nedge_type[1]); */
/*     phgPrintf("   nedge neumann : %d\n", nedge_type[2]); */
/*     phgPrintf("   nedge user1   : %d\n", nedge_type[3]); */
/*     phgPrintf("   nedge user2   : %d\n", nedge_type[4]); */
/*     phgPrintf("   nedge user3   : %d\n", nedge_type[5]); */
/*     phgPrintf("   nedge user4   : %d\n", nedge_type[6]); */
/*     phgPrintf("   nedge user5   : %d\n", nedge_type[7]); */
/*     phgPrintf("\n\n"); */
/* } */



/* void */
/* Grid::updateJacobian(Elem &e, const double *lambda) */
/* // */
/* // Update elem's dFdx[3][2], Jac[2][3], Det */
/* // */
/* { */
/*     const FEBasis *fe_basis = isop_map->fe_basis; */
/*     int N = fe_basis->nbas; */
/*     RealMat PointMat(N, SDim); */
/*     RealMat dshape(N, 2); */
/*     RealVec shape(N), xyz(SDim); */
/*     fe_basis->CalcDShapeRef(lambda, dshape); */
/*     fe_basis->CalcShapeRef(lambda, shape); */

/*     //      Mult(PointMat, dshape, dFdx); */
/*     //           3 x nbas  nbas x 2    3x2 */

/*     for (unsigned int i = 0; i < N; i++) { */
/* 		Real x, y; */
/* 		lambda2xyzDirect(e, */
/* 						 fe_basis->nodes[i][0], */
/* 						 fe_basis->nodes[i][1], */
/* 						 PointMat(i, 0), */
/* 						 PointMat(i, 1), */
/* 						 PointMat(i, 2) */
/* 						 ); */
/*     } */

/*     e.dFdx.setZero(); */
/*     e.dFdx = PointMat.transpose() * dshape; */
/* 	xyz = PointMat.transpose() * shape; */

/*     // std::cout << "dFdx" << std::endl; */
/*     // std::cout << e.dFdx << std::endl; */

/* #if 0 */
/*     e.Jac = e.dFdx.inverse(); */
/*     e.Det = e.dFdx.determinant(); */
/*     e.sgn = e.Det > 0 ? 1:-1; */
/* #else */
/* 	// */
/* 	// Det = || dF1 x dF2 || */
/* 	// Jac = (A^T A)^{-1} A^T */
/* 	// */ 

/* 	// TODO: Stabilized ??? */

/* 	//e.Jac = e.dFdx.completeOrthogonalDecomposition().pseudoInverse(); */

/* 	if (Params->use_proj_sphere) { */
/* 		const Real R = proj_sphere_radius; */
/* 		Real x = xyz(0); */
/* 		Real y = xyz(1); */
/* 		Real z = xyz(2); */
/* 		Real r = sqrt(x*x + y*y + z*z); */
			
/* 		RealMat	Jx(3, 3); */

/* 		Jx << y*y + z*z, -x*y	 , -x*z, */
/* 		    -x*y	  , x*x + z*z, -y*z, */
/* 		    -x*z	  , -y*z	 , x*x + y*y; */

/* 		Jx *= R / (r*r*r); */

/* 		// Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]"); */
/* 		// std::cout << Jx.format(HeavyFmt) << std::endl; */
		
/* 		// e.dFdxLn = e.dFdx; */
/* 		e.dFdx = Jx * e.dFdx; */
/* 	} */

	
/* 	 { */
/* 		// */
/* 		// */
/* 		// Physical */ 
/* 		// */ 
/* 		// */ 
/* 		Real d[6] = {e.dFdx(0,0), e.dFdx(1,0), e.dFdx(2,0),  // DF^T (2x3) */
/* 					 e.dFdx(0,1), e.dFdx(1,1), e.dFdx(2,1)}; */ 
/* 		Real ad[6]; */
/* 		Real e_, g_, f_; */
/* 		// DF^T DF = [e f] */
/* 		//           [f g] */
/* 		// adj(...) = [g, -f] */          
/* 		//            [-f, e] */
/* 		e_ = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; */
/* 		g_ = d[3]*d[3] + d[4]*d[4] + d[5]*d[5]; */
/* 		f_ = d[0]*d[3] + d[1]*d[4] + d[2]*d[5]; */

/* 		//  adj(DF^T x DF) x DF^T */
/* 		ad[0] = d[0]*g_ - d[3]*f_; */	
/* 		ad[1] = d[3]*e_ - d[0]*f_; */
/* 		ad[2] = d[1]*g_ - d[4]*f_; */
/* 		ad[3] = d[4]*e_ - d[1]*f_; */
/* 		ad[4] = d[2]*g_ - d[5]*f_; */
/* 		ad[5] = d[5]*e_ - d[2]*f_; */

/* 		e.Det = e_ * g_ - f_ * f_;	// det^2 */
/* 		e.Jac << ad[0], ad[2], ad[4], ad[1], ad[3], ad[5]; */
/* 		e.Jac /= e.Det; // inv(DF^T x DF) x DF */
/* 		e.Det = sqrt(e.Det); */
/* 		e.sgn = 1; */
/* 	} */

/* 	// { */
/* 	// 	// */
/* 	// 	// */
/* 	// 	// Linearized (without Sphere) */
/* 	// 	// */
/* 	// 	// */ 
/* 	// 	Real d[6] = {e.dFdxLn(0,0), e.dFdxLn(1,0), e.dFdxLn(2,0),  // DF */
/* 	// 				 e.dFdxLn(0,1), e.dFdxLn(1,1), e.dFdxLn(2,1)}; */ 
/* 	// 	Real ad[6]; */
/* 	// 	Real e_, g_, f_;			// DF^T DF = [e f] */
/* 	// 	//           [f g] */
/* 	// 	// */
/* 	// 	e_ = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; */
/* 	// 	g_ = d[3]*d[3] + d[4]*d[4] + d[5]*d[5]; */
/* 	// 	f_ = d[0]*d[3] + d[1]*d[4] + d[2]*d[5]; */

/* 	// 	//  adj(DF^T x DF) x DF */
/* 	// 	ad[0] = d[0]*g_ - d[3]*f_; */	
/* 	// 	ad[1] = d[3]*e_ - d[0]*f_; */
/* 	// 	ad[2] = d[1]*g_ - d[4]*f_; */
/* 	// 	ad[3] = d[4]*e_ - d[1]*f_; */
/* 	// 	ad[4] = d[2]*g_ - d[5]*f_; */
/* 	// 	ad[5] = d[5]*e_ - d[2]*f_; */

/* 	// 	e.DetLn = e_ * g_ - f_ * f_;	// det^2 */
/* 	// 	e.JacLn << ad[0], ad[2], ad[4], ad[1], ad[3], ad[5]; */
/* 	// 	e.JacLn /= e.DetLn; // inv(DF^T x DF) x DF */
/* 	// 	//e.Det = sqrt(e.Det); */
/* 	// 	//e.sgn = 1; */
/* 	// } */	
/* #endif */	

/*     // std::cout << "Jac: " */ 
/* 	// 		  << e.Jac << std::endl; */
/*     // std::cout << "Det: " */ 
/* 	// 		  << e.Det << std::endl; */

/*     // std::cout << "Det2: " << e.Det/2 << std::endl; */
/*     // std::cout << "Area: " << e.area << std::endl; */

/*     return; */
/* } */


/* void */
/* Grid::updateJacobian(Elem &e) */
/* { */
/*     const Coord lambda0 = {1./3., 1./3., 1./3.}; */
/*     updateJacobian(e, &lambda0[0]); */
/* } */


/* Real */
/* Elem::getArea(Real3 *normal, Real3 *df1, Real3 *df2) */
/* // */
/* // Area on point, depend on Det */
/* // */ 
/* { */
/* 	RealVec DF1 = dFdx.col(0); // dxi */
/* 	RealVec DF2 = dFdx.col(1); // deta */

/* 	if (normal) { */
/* 		CROSS_PRODUCT(DF1, DF2, (*normal)); */
/* 		NORMALIZE((*normal)); */
/* 	} */

/* 	if (df1) { */
/* 		*df1 = DF1; */
/* 		NORMALIZE((*df1)); */
/* 	} */

/* 	if (df2) { */
/* 		*df2 = DF2; */
/* 		NORMALIZE((*df2)); */
/* 	} */
	
/* 	return .5 * fabs(Det); */
/* } */


/* Real */
/* Elem::getEdgeLength(int iside, Real3 *gn) */
/* // See mfem.org for details */
/* { */
/* 	RealVec DF1 = dFdx.col(0); // dxi */
/* 	RealVec DF2 = dFdx.col(1); // deta */
/* 	RealVec DF3 = (-DF1 + DF2); // side */

/*     /* */
/*      * */
/*      *     + */
/*      *     |\ */ 
/*      *     |  \ _ */
/*      *     |   |\ */ 
/*      *     ^     DF3 */ 
/*      *     |        \ */ 
/*      *    DF2         \ */ 
/*      *     |            \ */ 
/*      *     + -- DF1 -->-- + */
/*      * */  
/*      *1/ */ 

/* 	Real length = 0; */
/* 	Real3 ne; */
/* 	getArea(&ne); */

/* 	if (iside == 0) { */
/* 		length = DF1.norm(); */
/* 		if (gn) { */
/* 			NORMALIZE(DF1); */
/* 			CROSS_PRODUCT(DF1, ne, (*gn)); */
/* 			// phgInfo(0, "kn: %30.15e %30.15e %30.15e\n", */ 
/* 			// 		ne(0), ne(1), ne(2)); */ 
/* 			// phgInfo(0, "df: %30.15e %30.15e %30.15e\n", */ 
/* 			// 		DF1(0), DF1(1), DF1(2)); */ 
/* 		} */
/* 	} */
/* 	else if (iside == 1) { */
/* 		length = DF3.norm(); */
/* 		if (gn) { */
/* 			NORMALIZE(DF3); */
/* 			CROSS_PRODUCT(DF3, ne, (*gn)); */
/* 		} */
/* 	} */
/* 	else if (iside == 2) { */
/* 		length = DF2.norm(); */
/* 		if (gn) { */
/* 			NORMALIZE(DF2); */
/* 			//CROSS_PRODUCT(DF2, ne, gn); */
/* 			// gn[0] *= -1; */
/* 			// gn[1] *= -1; */
/* 			// gn[2] *= -1; */
/* 			CROSS_PRODUCT(ne, DF2, (*gn)); */
/* 		} */
/* 	} */

/* 	if (0) { */
/* 		// */
/* 		// Compare to linear */
/* 		// */
/* 		Real len0 = edge_length[iside]; */
/* 		const Real *n0 = edge_normal[iside]; */

/* 		phgInfo(0, "elem %3d edge %d, %12.6f %12.6f %12.6e\n", */
/* 				index, iside, */
/* 				length, len0, fabs(length - len0)); */
/* 		phgInfo(0, "%12.6f %12.6f %12.6e\n", (*gn)[0], (*gn)[1], (*gn)[2]); */ 
/* 		phgInfo(0, "%12.6f %12.6f %12.6e\n", n0[0], n0[1], n0[2]); */ 
/* 	} */
	
	
/* 	return length; */
/* } */



/* void */
/* Grid::lambda2xyzDirect(const TriaElem &e, Real xi, Real eta, Real &x, Real &y, Real &z) const */
/* // */
/* // lambda -> xyz using Isop map */
/* // */ 
/* { */
/*     Real yeta = 1. - xi - eta; */

/* #if 0 */
/*     // Use simple linear */
/*     x = verts[e.verts[0]][0] * yeta */
/* 	    + verts[e.verts[1]][0] * xi */
/* 	    + verts[e.verts[2]][0] * eta; */
/*     y = verts[e.verts[0]][1] * yeta */
/* 	    + verts[e.verts[1]][1] * xi */
/* 	    + verts[e.verts[2]][1] * eta; */
/* #else */
/*     // Use isop */
/*     Coord lambda = {xi, eta, yeta}; */
/*     ScalarVec coord(isop_map->dim); */
/*     isop_map->eval(e, lambda, coord); */
/*     x = std::real(coord[X_DIR]); /1* extract the real part of a complex number *1/ */
/*     y = std::real(coord[Y_DIR]); */
/*     z = std::real(coord[Z_DIR]); */

/*     return; */
/* #endif */
/* } */



/* void */
/* Grid::lambda2xyzDirect(const QuadElem &e, Real xi, Real eta, Real &x, Real &y, Real &z) const */
/* // */
/* // lambda -> xyz using Isop map */
/* // */ 
/* { */
/*     Real yeta = 1. - xi - eta; */

/* #if 1 */
/*     // Use simple linear */
/*     x = verts[e.verts[0]][0] * yeta */
/* 	    + verts[e.verts[1]][0] * xi */
/* 	    + verts[e.verts[2]][0] * eta; */
/*     y = verts[e.verts[0]][1] * yeta */
/* 	    + verts[e.verts[1]][1] * xi */
/* 	    + verts[e.verts[2]][1] * eta; */
/* #else */
/*     // Use isop */
/*     Coord lambda = {xi, eta, yeta}; */
/*     ScalarVec coord(isop_map->dim); */
/*     isop_map->eval(e, lambda, coord); */
/*     x = std::real(coord[X_DIR]); /1* extract the real part of a complex number *1/ */
/*     y = std::real(coord[Y_DIR]); */
/*     z = std::real(coord[Z_DIR]); */

/*     return; */
/* #endif */

/* } */



/* void */
/* Grid::lambda2xyz(const Elem &e, Real xi, Real eta, */
/* 				 Real &x, Real &y, Real &z) const */
/* // */
/* // lambda ------> xyz' -------> xyz */ 
/* //         Isop         Proj sphere */
/* { */
/* 	lambda2xyzDirect(e, xi, eta, x, y, z); */
/* 	/1* phgInfo(3, "X %f %f %f\n", x, y, z); *1/ */

/* 	/1* if (Params->use_proj_sphere) { *1/ */
/* 	if (0) { */
/* 		const Real R = proj_sphere_radius; */
/* 		Real r = sqrt(x*x + y*y + z*z); */
/* 		x = x / r * R; */
/* 		y = y / r * R; */
/* 		z = z / r * R; */

/* 		phgInfo(3, "R %f %f %f\n", x, y, z); */
/* 	} */
/* } */



/* void */
/* Grid::xyz2lambda(const Elem &e, Real x, Real y, Real z, */
/* 				 Real &xi, Real &eta, Real &yeta) const */
/* // */
/* // xyz ------> xyz' -------> lambda */
/* //    Proj plane         Isop */          
/* { */
/* 	// assert(!Params->use_proj_sphere); */
/* 	// assert not use isop */

/* #if 0 */	
/* 	phgInfo(0, "%30.15e %30.15e %30.15e\n", */
/* 			verts[e.verts[0]][0], */
/* 			verts[e.verts[0]][1], */
/* 			verts[e.verts[0]][2] */
/* 			); */
/* 	phgInfo(0, "%30.15e %30.15e %30.15e\n", */
/* 			verts[e.verts[1]][0], */
/* 			verts[e.verts[1]][1], */
/* 			verts[e.verts[1]][2] */
/* 			); */
/* 	phgInfo(0, "%30.15e %30.15e %30.15e\n", */
/* 			verts[e.verts[2]][0], */
/* 			verts[e.verts[2]][1], */
/* 			verts[e.verts[2]][2] */
/* 			); */
/* 	phgInfo(0, "%30.15e %30.15e %30.15e\n\n", */
/* 			x, y, z */
/* 			); */

/* 	// DF */
/* 	phgInfo(0, "%30.15e %30.15e \n" */
/* 			"%30.15e %30.15e \n" */
/* 			"%30.15e %30.15e \n\n", */
/* 			e.dFdx(0,0), e.dFdx(0,1), */
/* 			e.dFdx(1,0), e.dFdx(1,1), */
/* 			e.dFdx(2,0), e.dFdx(2,1)); */

/* 	// Jac */
/* 	phgInfo(0, "%30.15e %30.15e %30.15e\n" */
/* 			"%30.15e %30.15e %30.15e\n\n", */
/* 			e.Jac(0,0), e.Jac(0,1), e.Jac(0,2), */
/* 			e.Jac(1,0), e.Jac(1,1), e.Jac(1,2)) ; */
/* #endif */	
	

/* 	if (Params->use_proj_sphere) { */
/* 		// */
/* 		// Proj to triangle */
/* 		// */
/* 		// TODO: cache */
/* 		Real x1[3] = { */
/* 			verts[e.verts[1]][0] - verts[e.verts[0]][0], */
/* 			verts[e.verts[1]][1] - verts[e.verts[0]][1], */
/* 			verts[e.verts[1]][2] - verts[e.verts[0]][2] */
/* 		}; */
/* 		Real x2[3] = { */
/* 			verts[e.verts[2]][0] - verts[e.verts[0]][0], */
/* 			verts[e.verts[2]][1] - verts[e.verts[0]][1], */
/* 			verts[e.verts[2]][2] - verts[e.verts[0]][2] */
/* 		}; */
/* 		Real n[3]; */
/* 		CROSS_PRODUCT(x1, x2, n); */
/* 		NORMALIZE(n); */

/* 		const Real *x0 = verts[e.verts[0]]; */
/* 		Real x0n = INNER_PRODUCT(x0, n); */
/* 		Real xn  = x*n[0] + y*n[1] + z*n[2]; */
/* 		Real t   = x0n / xn; */
	
/* 		x *= t; */
/* 		y *= t; */
/* 		z *= t; */
/* 	} */

/* 	// */
/* 	// x - x0 */
/* 	// */ 
/* 	x -= verts[e.verts[0]][0]; */
/* 	y -= verts[e.verts[0]][1]; */
/* 	z -= verts[e.verts[0]][2]; */

	
/* 	// phgInfo(0, "%30.15e %30.15e %30.15e\n", */
/* 	// 		x + verts[e.verts[0]][0], */
/* 	// 		y + verts[e.verts[0]][1], */
/* 	// 		z + verts[e.verts[0]][2] */
/* 	// 		); */

/*     xi   = e.Jac(0, 0) * x + e.Jac(0, 1) * y + e.Jac(0, 2) * z; */
/*     eta  = e.Jac(1, 0) * x + e.Jac(1, 1) * y + e.Jac(1, 2) * z; */
/*     yeta = 1. - xi - eta; */

	
/* 	return; */
/* } */

/* void */
/* Grid::xyz2lambda(const Elem &e, Real x, Real y, */
/* 				 Real &xi, Real &eta) const */
/* { */
/* 	/1* assert(!Params->use_proj_sphere); *1/ */

/* 	x -= verts[e.verts[0]][0]; */
/* 	y -= verts[e.verts[0]][1]; */
	
/*     xi   = e.Jac(0, 0) * x + e.Jac(0, 1) * y + e.Jac(0, 2) * 0; */
/*     eta  = e.Jac(1, 0) * x + e.Jac(1, 1) * y + e.Jac(1, 2) * 0; */

/* 	return; */
/* } */


/* void */
/* Grid::init_gref() */
/* { */

/*     /1* phgInfo(0, "Serial or local gref\n"); *1/ */
/*     cout << "Serial or local gref\n"; */

/*     // */
/*     // Vert to elem list, edge to elem list */
/*     // The element in the list in ordered so that the element indices increase. */
/*     // In case of crack mesh, the order is not garenteed for copied edges. */
/*     // */
/*     vert_to_elem.resize(nvert); */
/*     edge_to_elem.resize(nedge); */

/*     for (unsigned int ielem = 0; ielem < nelem; ielem++) { */
/* 		const Elem &e = elems[ielem]; */

/* 		for (unsigned int i = 0; i < NVert; i++) { */
/* 			int vid = e.verts[i]; */
/* 			gRef gr = {ielem, i}; */

/* 			vert_to_elem[vid].push_back(gr); */
/* 		} */

/* 		for (unsigned int i = 0; i < NEdge; i++) { */
/* 			int eid = e.edges[i]; */
/* 			gRef gr = {ielem, i}; */
/* 			std::vector<gRef > &e2E = edge_to_elem[eid]; */

/* 			// Sort edge to elem list */
/* 			if (e2E.size() == 1 */
/* 				&& e2E[0].eind > ielem ) { */
/* 				e2E.insert(e2E.begin(), gr); */
/* 			} */
/* 			else { */
/* 				e2E.push_back(gr); */
/* 			} */
/* 		} */
/*     } */

/*     // Sort vert to elem list */
/*     for (unsigned int ivert = 0; ivert < nvert; ivert++) { */
/* 		sort(begin(vert_to_elem[ivert]), */
/* 			 end(vert_to_elem[ivert]), */
/* 			 [](gRef a, gRef b) {return a.eind > b.eind; }); */

/* 		// for (auto item : vert_to_elem[ivert]) */
/* 		//     printf("%d(%d) ", item.eind, item.gind); */
/* 		// printf("\n"); */
/* 		// if (vert_to_elem[ivert].size() == 1) */
/* 		//     printf("%f %f\n", verts[ivert][0], verts[ivert][1]); */
/*     } */

/*     // assign element edge signs */
/*     for (unsigned int iedge = 0; iedge < nedge; iedge++) { */
/* 		int ne = edge_to_elem[iedge].size(); */
/* 		assert(ne == 1 || ne == 2); */

/* 		// First elem has POISITIVE edge sgn. */
/* 		// In FESpace::interp, the value is interpolated from first element. */
/* 		gRef *gr = &edge_to_elem[iedge][0]; */
/* 		elems[gr->eind].edges_sgn[gr->gind] = 1; */

/* 		if (ne == 2) { */
/* 			// Second elem has NEGETIVE edge sgn */
/* 			gr = &edge_to_elem[iedge][1]; */
/* 			elems[gr->eind].edges_sgn[gr->gind] = -1; */
/* 		} */
/*     } */


/* #if 0 */
/*     // Check */
/*     for (unsigned int i = 0; i < nvert; i++) { */
/* 		phgInfo(2, "Vert[%5d]: ", i); */
/* 		for (unsigned int j = 0; j < vert_to_elem[i].size(); j++) { */
/* 			phgInfo(2, "(%4d, %d), ",  vert_to_elem[i][j].eind, */
/* 					vert_to_elem[i][j].gind */
/* 					); */
/* 		} */
/* 		phgInfo(2, "\n"); */
/*     } */
/* #endif */
/* #if 0 */	
/*     for (unsigned int i = 0; i < nedge; i++) { */
/* 		phgInfo(2, "Edge[%5d]: ", i); */
/* 		for (unsigned int j = 0; j < edge_to_elem[i].size(); j++) { */
/* 			phgInfo(2, "(%4d, %d), ",  edge_to_elem[i][j].eind, */
/* 					edge_to_elem[i][j].gind */
/* 					); */
/* 		} */
/* 		phgInfo(2, "\n"); */
/*     } */
/* #endif */
/*     return; */
/* } */



/* /1* void *1/ */
/* /1* ParallelGrid::init_gref() *1/ */
/* /1* { *1/ */

/* /1*     Grid::init_gref(); *1/ */

/* /1*     if (phgNProcs == 1) *1/ */
/* /1* 		return; *1/ */

/* /1*     phgInfo(0, "Override edge sign.\n"); *1/ */

/* /1*     // Overide element edge signs *1/ */
/* /1*     //    using element global index *1/ */
/* /1*     for (unsigned int iedge = 0; iedge < nedge; iedge++) { *1/ */
/* /1* 		unsigned int Iedge = L2Gmap_edge[iedge]; *1/ */

/* /1* 		int ne = global.edge_to_elem[Iedge].size(); *1/ */
/* /1* 		assert(ne == 1 || ne == 2); *1/ */

/* /1* 		// First elem has POISITIVE edge sgn. *1/ */
/* /1* 		// In FESpace::interp, the value is interpolated from first element. *1/ */
/* /1* 		gRef *gr = &global.edge_to_elem[Iedge][0]; *1/ */

/* /1* 		int ielem1 = G2Lmap_elem[gr->eind]; *1/ */
/* /1* 		if (ielem1 >= 0) *1/ */
/* /1* 			elems[ielem1].edges_sgn[gr->gind] = 1; *1/ */

/* /1* 		if (ne == 2) { *1/ */
/* /1* 			// Second elem has NEGETIVE edge sgn *1/ */
/* /1* 			gr = &global.edge_to_elem[Iedge][1]; *1/ */

/* /1* 			int ielem2 = G2Lmap_elem[gr->eind]; *1/ */
/* /1* 			if (ielem2 >= 0) *1/ */
/* /1* 				elems[ielem2].edges_sgn[gr->gind] = -1; *1/ */
/* /1* 		} *1/ */
/* /1*     } *1/ */



/* /1*     return; *1/ */
/* /1* } *1/ */











/* unsigned int */
/* Grid::until_elem(unsigned int ielem, Btype btype) */
/* // */
/* // Step utils elem of btype */
/* // */
/* { */
/*     while (ielem < nelem && !(types_elem[ielem] & btype) ) */
/* 		ielem++; */

/*     return ielem; */
/* } */





int main(int argc, char *argv[]) {
    Grid *g = new Grid();
    g->read_mesh("mixedR1.msh");

    for (int i = 0; i < g->nelem; i++) {
        Elem e = g->elems[i];
        cout << e.elem_type << endl; 
    }

    // ------------------------------ edge test ------------------------------
    /* cout << "number of edges: " << g->nedge << endl; */
    /* cout << "edges\' index and nodes:\n"; */
    /* for (int i = 0; i < g->nedge; i++) { */
    /*     cout << i << "\t" << g->edges[i][0] << "\t" << g->edges[i][1] << endl; */
    /* } */
    // ------ quadrilateral ------
    /* cout << endl << "quadrilateral elements\' edge index:" << endl; */
    /* for (int i = 0; i < g->nquad_elem; i++) { */
    /*     cout << i << "\t"; */
    /*     QuadElem e = g->quad_elems[i]; */
    /*     for (int j = 0; j < 4; j++) { */
    /*         cout << e.edges[j] << "\t"; */ 
    /*     } */
    /*     cout << endl; */
    /* } */
    // ------ triangle ------
    /* cout << endl << "triangle elements\' edge index:" << endl; */
    /* for (int i = 0; i < g->ntria_elem; i++) { */
    /*     cout << i << "\t"; */
    /*     TriaElem e = g->tria_elems[i]; */
    /*     for (int j = 0; j < 3; j++) { */
    /*         cout << e.edges[j] << "\t"; */ 
    /*     } */
    /*     cout << endl; */
    /* } */

    /* for (int i = 0; i < nvert; i++) { */
    /*     cout << i << "\t"; */
    /*     for (int j = 0; j < 3; j++) { */
    /*         cout << verts[i][j] << "\t"; */
    /*     } */
    /*     cout << endl; */
    /* } */
    return 0;
}
