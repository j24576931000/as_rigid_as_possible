#pragma once
#include "TriangulationBase.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
class TriangulationCgal : public TriangulationBase
{
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K,Tds> Triangulation;

	typedef CGAL::Delaunay_mesh_size_criteria_2<Triangulation> Criteria;
	typedef CGAL::Delaunay_mesher_2<Triangulation, Criteria> Mesher;

	typedef Triangulation::Face_handle Face_handle;
	typedef Triangulation::Finite_faces_iterator Finite_faces_iterator;
	typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
	typedef Triangulation::Point  Point;
	typedef Triangulation::Vertex_handle Vertex_handle;

	TriangulationCgal(void);
	~TriangulationCgal(void);
	virtual void AddPoint(float x, float y);
	virtual void Compute();
	virtual void DelaunayMesher2();


	Vector2s MeshPointSet();
	int getVertexID(float x, float y);

	std::vector<Vertex_handle> ContourPoint;
	Triangulation	m_Triangulation;
};

