extern "C"
{
#include "fsls.h"
}
#include "rcm.hpp"

int RCM( fsls_CSRMatrix *A, fsls_XVector *f, fsls_XVector *u, int nt )
{
	using namespace boost;
	using namespace std;
	int ierr = 0, row = 0, col = 0, i = 0, j = 0;
	int num_rows = fsls_CSRMatrixNumRows(A);
	int num_cols = fsls_CSRMatrixNumCols(A);
	int nnz = fsls_CSRMatrixNumNonzeros(A);
	if (num_rows != num_cols)
		printf("is there anything wrong with input matrix A?\n");
	int ngrid = num_rows;
	int *matrix_i = fsls_CSRMatrixI(A);
	int *matrix_j = fsls_CSRMatrixJ(A);
	double *matrix_data = fsls_CSRMatrixData(A);

	int *new_i = fsls_CTAlloc(int, num_rows + 1);
	new_i[num_rows] = matrix_i[num_rows];
	int *new_j = fsls_CTAlloc(int, nnz);
	double *new_data = fsls_CTAlloc(double, nnz);
	double *new_f, *new_u;
	if (nt == 0)
	{
		new_f    = fsls_CTAlloc(double, ngrid);
		new_u    = fsls_CTAlloc(double, ngrid);
	}
	else
	{
		new_f    = fsls_CTAlloc(double, ngrid*nt);
		new_u    = fsls_CTAlloc(double, ngrid*nt);
	}
	typedef adjacency_list<vecS, vecS, undirectedS, 
					property<vertex_color_t, default_color_type,
					property<vertex_degree_t,int> > > Graph;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	typedef graph_traits<Graph>::vertices_size_type size_type;

	typedef std::pair<std::size_t, std::size_t> Pair;
	Pair edges[nnz];
//	Pair new_edges[nnz];
	for (i = 0; i < nnz; ++i)
	{
		if (i==matrix_i[row+1])
			row = row+1;
		col = matrix_j[i];
		{
			edges[j].first  = row;
			edges[j].second = col;
			j++;
		}
	}

	Graph G(ngrid);
//	Graph new_G(ngrid);
	for (i = 0; i < nnz; ++i)
	{
		add_edge(edges[i].first, edges[i].second, G);
	}
	graph_traits<Graph>::vertex_iterator ui, ui_end;

	property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
	for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
		deg[*ui] = degree(*ui, G);

	property_map<Graph, vertex_index_t>::type
		index_map = get(vertex_index, G);

	std::cout << "original bandwidth: " << bandwidth(G) << std::endl;

	std::vector<Vertex> new_idx(num_vertices(G));
	std::vector<size_type> perm(num_vertices(G));
	{
		//reverse cuthill_mckee_ordering
		cuthill_mckee_ordering(G, new_idx.rbegin(), get(vertex_color, G),
					make_degree_map(G));

		for (size_type c = 0; c != new_idx.size(); ++c)
		{
			perm[index_map[new_idx[c]]] = c;
		}

		int flag = 0;
		for (i = 0; i < ngrid; ++i)
		{
			new_i[i] = flag;
			row = matrix_i[new_idx[i]];
			if (nt == 0)
			{
				new_f[i] = f->data[new_idx[i]];
				new_u[i] = u->data[new_idx[i]];
			}
			else
			{
				for (j = 0; j < nt; ++j)
				{
					new_f[i+ngrid*j] = f->data[new_idx[i]+ngrid*j];
					new_u[i+ngrid*j] = u->data[new_idx[i]+ngrid*j];
				}
			}
			while ( edges[row].first == new_idx[i])
			{
				new_j[flag]    = perm[edges[row].second];
				new_data[flag] = matrix_data[row];
				row++;
				flag ++;
			}
		}
		fsls_CSRMatrixI(A) = new_i;
		fsls_CSRMatrixJ(A) = new_j;
		fsls_CSRMatrixData(A) = new_data;
		fsls_XVectorData(f) = new_f;
		fsls_XVectorData(u) = new_u;
//		row = 0, j = 0;
//		for (i = 0; i < nnz; ++i)
//		{
//			if (i==new_i[row+1])
//				row = row+1;
//			col = new_j[i];
//			{
//				new_edges[j].first  = row;
//				new_edges[j].second = col;
//				j++;
//			}
//		}
//		Graph new_G(ngrid);
//		for (i = 0; i < nnz; ++i)
//		{
//			add_edge(new_edges[i].first, new_edges[i].second, new_G);
//		}

//		std::cout << "new original bandwidth: " << bandwidth(new_G) << std::endl;

		std::cout << "new(RCM) bandwidth: " 
			<< bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
			<< std::endl;
	}
	return ierr;
}
