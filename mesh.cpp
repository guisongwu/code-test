#include <iostream>
#include <cstring>

typedef unsigned short Btype;
typedef Btype (*BC_FUNCTION)(int mark);


bool get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (true) {
        memset(token, 0, 100);
        if (fscanf(fp, "%s", token) != 1)  // everytime 'fscanf' is called, 
                                           // the pointer will move to the end of this string.
            return false;
        if (token[0] != '#')  
            break;
        /* skip to newline */
        do {
            if ((c = fgetc(fp)) == EOF)
                return false;
        } while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
        *p = '\0';
    return true;
}



void read_mesh(const char *mesh_file_name, BC_FUNCTION user_bc_map, const char *part_file_name)
{
    int nbdry_edge;
    BdryEdge *bdry_edges;


    if (user_bc_map == NULL)
		bc_map = default_bc_map;
    else
		bc_map = user_bc_map;

    if (verb >=1 && phgRank == 0)
		cout << "Read Mesh " << mesh_file_name << endl;

#define READ_NUMBER													\
    if (!get_token(fp, token)) strcpy(token, "End");				\
    if (isalpha((int)(token[0]))) {									\
		fprintf(stderr, "fewer entries (%d) than expected.\n", i);	\
		break;														\
    }
#undef ERROR
#define ERROR {n = __LINE__; goto error;}
    
    FILE *fp ;
    char token[MAX_TOKEN_LEN]; 
    int n;
    if ((fp = fopen(mesh_file_name, "r")) == NULL) {
		phgPrintf("can't open mesh file \"%s\"!\n", mesh_file_name);
		phgAbort(0);
		return;
    }

    token[0] = '\0';
    /* if ( */
    /*         !get_token(fp, token) ||          // Read token (fail if EOF/error) */
    /*         strcasecmp(token, "$MeshFormat") || // Must be "$MeshFormat" (case-insensitive) */
    /*         !get_token(fp, token) ||           // Read next token */
    /*         strcmp(token, "2.2") ||            // Must be "2.2" (exact match) */
    /*         !get_token(fp, token) ||           // Read next token */
    /*         strcmp(token, "0") ||              // Must be "0" */
    /*         !get_token(fp, token) ||           // Read next token */
    /*         strcmp(token, "8") ||              // Must be "8" */
    /*         !get_token(fp, token) ||           // Read next token */
    /*         strcasecmp(token, "$EndMeshFormat") // Must be "$EndMeshFormat" (case-insensitive) */
    /*    ) */
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

    while (true) { // end at line 316
		if (!get_token(fp, token))
			break;
		//next_token:
		if (!strcasecmp(token, "$Nodes")) {
			if (!get_token(fp, token))
				ERROR;
			n = atoi(token);
			if (verb >= 1 && phgRank == 0)
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
			if (verb >= 1 && phgRank == 0)
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
				// edge or triangle
				assert(elem_type == 1 || elem_type == 2);

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
				else if (elem_type == 2) {
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

    edges	= new Edge[nelem * 3]();
    for (unsigned int i = 0; i < nelem; i++) {
		int v0, v1;
		for (unsigned int j = 0; j < 3; j++) {
			v0 = elems[i].verts[GetEdgeVert(j, 0)];
			v1 = elems[i].verts[GetEdgeVert(j, 1)];

			SortIndex(v0, v1);
	    
			edges[i*3+j][0] = v0;
			edges[i*3+j][1] = v1;
		}

		// Ordering
		{
			int V0, V1, V2, v0, v1, v2;
			V0 = elems[i].verts[0];
			V1 = elems[i].verts[1];
			V2 = elems[i].verts[2];

			v0 = v1 = v2 = 0;
			(V0 > V1) ? v0++ : v1++;
			(V0 > V2) ? v0++ : v2++;
			(V1 > V2) ? v1++ : v2++;

			elems[i].ordering[0] = v0;
			elems[i].ordering[1] = v1;
			elems[i].ordering[2] = v2;
		}
    }
    // for (i = 0; i < 3*nelem; i++)
    // 	phgInfo(2, "edge %d: %d %d\n", i, edges[i][0], edges[i][1]);

    qsort(edges, nelem*3, sizeof(Edge), compare_edge);

    int count = 0; 
    {
		int i0 = 0;
		for (unsigned int i = i0+1; i < 3*nelem; i++) {
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

		for (unsigned int j = 0; j < NEdge; j++) {
			int v0 = e->verts[GetEdgeVert(j, 0)];
			int v1 = e->verts[GetEdgeVert(j, 1)];

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
				phgAbort(0);

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
			phgAbort(0);  // edge with no element
		}
    }  

    // // Dump neighs
    // for (unsigned int ielem = 0; ielem < nelem; ielem++) {
    // 	Elem *e = &elems[ielem];
    // 	phgInfo(2, "e neigh %3d: %3d %3d %3d\n", ielem,
    // 	       e->neigh[0], e->neigh[1], e->neigh[2]
    // 	       );
    // }

	partition(part_file_name);
    get_local_grid();

    
    init_gref();
    init_geom();


    // phgInfo(0, "debug exit\n");
    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);
    statistic();

    return;
}

int main() {


}
