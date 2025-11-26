#pragma once
#include <iostream>
#include <vector>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>
#include<Eigen/Dense> 
#include <CGAL/Surface_mesh.h>
#include <glm/glm.hpp>

#include <CGAL/convex_hull_3.h> 
#include <CGAL/Polygon_mesh_processing/corefinement.h> 



#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>


#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Polygon_mesh_slicer.h>


#include "Runcuda.h"

#include "trimesh.h"

#include "earcut.hpp"


#define PI 3.1415926


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Facet_const_handle Facet_const_handle;
typedef Polyhedron::Point_3 Point3;

typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef K::Point_2 Point2;


typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef CGAL::AABB_halfedge_graph_segment_primitive<SurfaceMesh>     HGSP;
typedef CGAL::AABB_tree<CGAL::AABB_traits<K, HGSP>>           AABB_tree;

typedef CGAL::Alpha_shape_vertex_base_3<K>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb, Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>       Triangulation;
typedef CGAL::Alpha_shape_3<Triangulation>         Alpha_shape;
typedef Alpha_shape::Alpha_iterator                Alpha_iterator;

typedef K::FT                                               FT;
typedef K::Weighted_point_2                                 Weighted_point_2;
typedef K::Segment_2                                        Segment;

typedef CGAL::Regular_triangulation_vertex_base_2<K>        Rvb_2;
typedef CGAL::Alpha_shape_vertex_base_2<K, Rvb_2>              Vb_2;
typedef CGAL::Regular_triangulation_face_base_2<K>          Rf_2;
typedef CGAL::Alpha_shape_face_base_2<K, Rf_2>                 Fb_2;
typedef CGAL::Triangulation_data_structure_2<Vb_2, Fb_2>         Tds_2;
typedef CGAL::Regular_triangulation_2<K, Tds_2>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_iterator                Alpha_iterator_2;



namespace PMP = CGAL::Polygon_mesh_processing;


void Polyhedron2Trimesh(const Polyhedron& poly, trimesh::TriMesh* tri);

struct VoxelSize {
	double x;
	double y;
	double z;
    
    VoxelSize() : x(0), y(0), z(0) {}
    VoxelSize(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    VoxelSize(const double& edge) : x(edge), y(edge), z(edge) {}
};

struct VoxelColor {
    double r;
    double g;
    double b;
    VoxelColor(double r = 0, double g = 0, double b = 0) {
        this->r = r;
		this->g = g;
		this->b = b;
    }
};


class SAMPLE_ON_BALL
{
public:
    void EnumerationOritation();
    std::vector<Vector> sample_points;
    
    int num_ori_sample;
};


inline bool compare_bbx_volume_trimesh(trimesh::TriMesh& p1, trimesh::TriMesh& p2) {
    p1.need_bbox();
    p2.need_bbox();

    double volume1 = (p1.bbox.max.x - p1.bbox.min.x) * (p1.bbox.max.y - p1.bbox.min.y) * (p1.bbox.max.z - p1.bbox.min.z);
    double volume2 = (p2.bbox.max.x - p2.bbox.min.x) * (p2.bbox.max.y - p2.bbox.min.y) * (p2.bbox.max.z - p2.bbox.min.z);
	
    return volume1 > volume2;   
}

inline bool compare_volume(const Polyhedron& p1, const Polyhedron& p2) {
    double volume1 = CGAL::Polygon_mesh_processing::volume(p1);
    double volume2 = CGAL::Polygon_mesh_processing::volume(p2);
    return volume1 > volume2;   
}

inline bool compare_bbx_volume(const Polyhedron& p1, const Polyhedron& p2) {
    Tree tree1(faces(p1).first, faces(p1).second, p1);
    tree1.accelerate_distance_queries();

    CGAL::Bbox_3 bbox1 = tree1.bbox();

    Tree tree2(faces(p2).first, faces(p2).second, p2);
    tree2.accelerate_distance_queries();

    CGAL::Bbox_3 bbox2 = tree2.bbox();

    double volume1 = (bbox1.xmax() - bbox1.xmin()) * (bbox1.ymax() - bbox1.ymin()) * (bbox1.zmax() - bbox1.zmin());
    double volume2 = (bbox2.xmax() - bbox2.xmin()) * (bbox2.ymax() - bbox2.ymin()) * (bbox2.zmax() - bbox2.zmin());
    return volume1 > volume2;   
}



void voxelizeMesh(
    const Polyhedron& mesh,
    const size_t gridX,
    const size_t gridY,
    const size_t gridZ,
    const VoxelSize& voxelS,
    std::vector<bool>& voxelGrid);


void Export_Voxel(std::ofstream& export_file_output, int& export_index,
    std::string s_name, const Point& point, const VoxelColor& voxelcolor, const VoxelSize& voxelsize);
void Export_Line(std::ofstream& export_file_output, int& export_index, std::string s_name,
    const Point& start, const Point& end, const double radius, const VoxelColor& voxelcolor);

void Bool2Voxel(const std::string& filename, const std::vector<bool>& h_, const VoxelColor& voxelcolor, const VoxelSize& voxelsize);

void Bool2Voxel(const std::string& filename, const std::vector<bool>& h_, const VoxelColor& voxelcolor, const VoxelSize& voxelsize, const std::vector<int>& bbox);



void Int2VoxelUseColorBar(const std::string& filename, const std::vector<float>& h_, const VoxelSize& voxelsize);
void Int2VoxelUseColorBar(const std::string& filename, const std::vector<int>& h_, const VoxelSize& voxelsize);

void translateMesh(Polyhedron& mesh, const Vector& translation);
void translateMesh(trimesh::TriMesh& mesh, const Vector& translation);



void rotateMesh(Polyhedron& mesh, const Vector& ori);

void rotateMesh(trimesh::TriMesh& mesh, const Vector& ori);

void scaleMesh(Polyhedron& mesh, const double& scale);
void scaleMesh(Polyhedron& mesh, const Vector& scale);


void scaleMesh(trimesh::TriMesh& mesh, const double& scale);


void compute_intersection(std::vector<Polyhedron>& meshes, Polyhedron& intersection_result);
double computeDistance(const Point& p1, const Point& p2);

std::vector<std::pair<int, int>> computeEdges(const std::vector<Point>& points, double epsilon);


Point computeCentroid(const std::vector<Point>& points);

std::vector<Point> generateOffset(const std::vector<Point>& points, const Vector& offsetDirection, double offsetDistance);

void computeBoundary(const std::vector<Point>& vecs, const std::vector<std::pair<int, int>> edges, std::set<int>& boundaries);


void generateMesh(const std::vector<Point>& original, const Vector& offsetDirection, const double offsetLength, Polyhedron& mesh);

VoxelColor getVoxelColor(double value);
