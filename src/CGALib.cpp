#include "CGALib.h"




void SAMPLE_ON_BALL::EnumerationOritation()
{
    sample_points.clear();
    int number_points = num_ori_sample;    

    std::vector<Vector> ori_vec;
    if (NZ == 1) {

        ori_vec.push_back(Vector(1, 0, 0));       
        ori_vec.push_back(Vector(-1, 0, 0));      
        ori_vec.push_back(Vector(0, 1, 0));       
        ori_vec.push_back(Vector(0, -1, 0));      

        float inv_sqrt2 = 1.0 / sqrt(2.0);
        ori_vec.push_back(Vector(inv_sqrt2, inv_sqrt2, 0));      
        ori_vec.push_back(Vector(inv_sqrt2, -inv_sqrt2, 0));     
        ori_vec.push_back(Vector(-inv_sqrt2, inv_sqrt2, 0));     
        ori_vec.push_back(Vector(-inv_sqrt2, -inv_sqrt2, 0));    
    }
    else {  
        ori_vec.push_back(Vector(0, 1, 0));
        ori_vec.push_back(Vector(0, -1, 0));
        ori_vec.push_back(Vector(0, 0, 1));
        ori_vec.push_back(Vector(0, 0, -1));
        ori_vec.push_back(Vector(1, 0, 0));
        ori_vec.push_back(Vector(-1, 0, 0));

        float inv_sqrt3 = 1.0 / sqrt(3.0);
        ori_vec.push_back(Vector(inv_sqrt3, inv_sqrt3, inv_sqrt3));        
        ori_vec.push_back(Vector(inv_sqrt3, inv_sqrt3, -inv_sqrt3));       
        ori_vec.push_back(Vector(inv_sqrt3, -inv_sqrt3, inv_sqrt3));       
        ori_vec.push_back(Vector(inv_sqrt3, -inv_sqrt3, -inv_sqrt3));      
        ori_vec.push_back(Vector(-inv_sqrt3, inv_sqrt3, inv_sqrt3));       
        ori_vec.push_back(Vector(-inv_sqrt3, inv_sqrt3, -inv_sqrt3));      
        ori_vec.push_back(Vector(-inv_sqrt3, -inv_sqrt3, inv_sqrt3));      
        ori_vec.push_back(Vector(-inv_sqrt3, -inv_sqrt3, -inv_sqrt3));     
    }
    if (number_points == 1)
    {
        sample_points.push_back(ori_vec[3]);
    }
    else if (number_points == 2) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[4]);

    }
    else if (number_points == 3) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[2]);
        sample_points.push_back(ori_vec[4]);
    }
    else if (number_points == 6 || number_points == 14) {
        for (int i = 0; i < number_points; i++) {
            sample_points.push_back(ori_vec[i]);
        }
    }
    else if (number_points == 8) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[1]);
        sample_points.push_back(ori_vec[4]);
        sample_points.push_back(ori_vec[5]);
        sample_points.push_back(ori_vec[7]);
        sample_points.push_back(ori_vec[9]);
        sample_points.push_back(ori_vec[11]);
        sample_points.push_back(ori_vec[13]);
    }
    else if (number_points == 9) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[1]);
        sample_points.push_back(ori_vec[3]);
        sample_points.push_back(ori_vec[4]);
        sample_points.push_back(ori_vec[5]);
        sample_points.push_back(ori_vec[7]);
        sample_points.push_back(ori_vec[9]);
        sample_points.push_back(ori_vec[11]);
        sample_points.push_back(ori_vec[13]);
    }
    else if (number_points == 10) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[1]);
        sample_points.push_back(ori_vec[2]);
        sample_points.push_back(ori_vec[3]);
        sample_points.push_back(ori_vec[4]);
        sample_points.push_back(ori_vec[5]);
        sample_points.push_back(ori_vec[7]);
        sample_points.push_back(ori_vec[9]);
        sample_points.push_back(ori_vec[11]);
        sample_points.push_back(ori_vec[13]);
    }
    else if (number_points == 18) {
        for (int i = 0; i < 14; i++)
            sample_points.push_back(ori_vec[i]);

        float inv_sqrt2 = 1.0 / sqrt(2.0);
        sample_points.push_back(Vector(inv_sqrt2, inv_sqrt2, 0));        
        sample_points.push_back(Vector(-inv_sqrt2, -inv_sqrt2, 0));       
        sample_points.push_back(Vector(inv_sqrt2, -inv_sqrt2, 0));       
        sample_points.push_back(Vector(inv_sqrt2, -inv_sqrt2, 0));      
    }
    else if (number_points == 3) {
        sample_points.push_back(ori_vec[0]);
        sample_points.push_back(ori_vec[2]);
        sample_points.push_back(ori_vec[4]);
    }
}


void Polyhedron2Trimesh(const Polyhedron& poly, trimesh::TriMesh* tri)
{
    std::vector<trimesh::point> trimeshVertices;
    std::vector<trimesh::TriMesh::Face> trimeshFaces;

    for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
        const Point& p = v->point();
        trimeshVertices.emplace_back(p.x(), p.y(), p.z());
    }

    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        auto h = f->facet_begin();
        std::vector<int> indices;
        do {
            indices.push_back(std::distance(poly.vertices_begin(), h->vertex()));
            ++h;
        } while (h != f->facet_begin());

        if (indices.size() == 3) {     
            trimeshFaces.emplace_back(indices[0], indices[1], indices[2]);
        }
    }

    tri->vertices = trimeshVertices;
    tri->faces = trimeshFaces;
}

void voxelizeMesh(
    const Polyhedron& mesh,
    const size_t gridX,
    const size_t gridY,
    const size_t gridZ,
    const VoxelSize &voxelS,
    std::vector<bool>& voxelGrid)
{
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);
    tree.accelerate_distance_queries();

    CGAL::Bbox_3 bbox = tree.bbox();
    Point minPoint(bbox.xmin(), bbox.ymin(), bbox.zmin());
    Point maxPoint(bbox.xmax(), bbox.ymax(), bbox.zmax());

    std::cout << "BoundingBox span: " << bbox.x_span() << " " << bbox.y_span() << " " << bbox.z_span() << std::endl;

    double voxelSizeX = voxelS.x;
    double voxelSizeY = voxelS.y;
    double voxelSizeZ = voxelS.z;

    VoxelColor voxelcolor;
    voxelcolor.r = voxelcolor.g = voxelcolor.b = 1;

    int index = 0;

    for (size_t x = 0; x < (bbox.x_span() / voxelSizeX + 1) && x < gridX; ++x) {
        for (size_t y = 0; y < (bbox.y_span() / voxelSizeY + 1) && y < gridY; ++y) {
            for (size_t z = 0; z < (bbox.z_span() / voxelSizeZ + 1) && z < gridZ; ++z) {
                Point voxelCenter(
                    minPoint.x() + (x + 0.5) * voxelSizeX,
                    minPoint.y() + (y + 0.5) * voxelSizeY,
                    minPoint.z() + (z + 0.5) * voxelSizeZ
                );
                
                Point RayPoint(
                    minPoint.x() + (x + 1) * voxelSizeX,
                    minPoint.y() + (y + 0.5) * voxelSizeY,
                    minPoint.z() + (z + 0.5) * voxelSizeZ
                );

                Ray Ray_query(voxelCenter, RayPoint);

                if (tree.number_of_intersected_primitives(Ray_query) % 2 == 1) {
					int idx = x * gridY * gridZ + y * gridZ + z;
                    voxelGrid[idx] = true;  
                    index++;
                }
                else {
                    int idx = x * gridY * gridZ + y * gridZ + z;
                    voxelGrid[idx] = false;  
                }
            }
        }
    }
}


void Export_Voxel(std::ofstream& export_file_output, int& export_index,
    std::string s_name, const Point &point, const VoxelColor &voxelcolor, const VoxelSize &voxelsize)
{	


    std::vector<Point> vecs;
    vecs.push_back(Point(0.5 * voxelsize.x, 0.5 * voxelsize.y, 0.5 * voxelsize.z));
    vecs.push_back(Point(-0.5 * voxelsize.x, 0.5 * voxelsize.y, 0.5 * voxelsize.z));
    vecs.push_back(Point(-0.5 * voxelsize.x, 0.5 * voxelsize.y, -0.5 * voxelsize.z));
    vecs.push_back(Point(0.5 * voxelsize.x, 0.5 * voxelsize.y, -0.5 * voxelsize.z));

    vecs.push_back(Point(0.5 * voxelsize.x, -0.5 * voxelsize.y, 0.5 * voxelsize.z));
    vecs.push_back(Point(-0.5 * voxelsize.x, -0.5 * voxelsize.y, 0.5 * voxelsize.z));
    vecs.push_back(Point(-0.5 * voxelsize.x, -0.5 * voxelsize.y, -0.5 * voxelsize.z));
    vecs.push_back(Point(0.5 * voxelsize.x, -0.5 * voxelsize.y, -0.5 * voxelsize.z));

    std::vector<std::vector<int>> faces;

    int face_index_0[4] = { 0, 1, 2, 3 };
    int face_index_1[4] = { 5, 1, 0, 4 };
    int face_index_2[4] = { 4, 0, 3, 7 };
    int face_index_3[4] = { 5, 4, 7, 6 };
    int face_index_4[4] = { 7, 3, 2, 6 };
    int face_index_5[4] = { 6, 2, 1, 5 };

    faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
    faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
    faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
    faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
    faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
    faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

    export_file_output << "# " + s_name << std::endl;

    for (int i = 0; i < vecs.size(); i++)
    {
        Point p(point[0] + vecs[i][0], point[1] + vecs[i][1], point[2] + vecs[i][2]);

        export_file_output << "v " << p[0] << " " << p[1] << " " << p[2] << " " << voxelcolor.r << " " << voxelcolor.g << " " << voxelcolor.b << std::endl;
    }
    for (int i = 0; i < faces.size(); i++)
    {
        export_file_output << "f ";

        for (int j = faces[i].size() - 1; j >= 0; j--)
        {
            export_file_output << faces[i][j] + export_index << " ";
        }
        export_file_output << "" << std::endl;
    }
    export_index += 8;
}

inline Vector RotationAxis(Vector p, double angle, Vector n)
{
    glm::mat4 inputMatrix(0.0);
    inputMatrix[0][0] = p[0];
    inputMatrix[1][0] = p[1];
    inputMatrix[2][0] = p[2];
    inputMatrix[3][0] = 1.0;
    double u = n[0];
    double v = n[1];
    double w = n[2];

    glm::mat4  rotationMatrix;

    double L = (u * u + v * v + w * w);

    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;

    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][3] = 0.0;

    rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][3] = 0.0;

    rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix[2][3] = 0.0;

    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;

    double outputMatrix[4][1] = { 0.0, 0.0, 0.0, 0.0 };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 1; j++) {
            outputMatrix[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
    return Vector(outputMatrix[0][0], outputMatrix[0][1], outputMatrix[0][2]);
}


void Export_Line(std::ofstream& export_file_output, int& export_index, std::string s_name,
    const Point& start, const Point& end,const double radius, const VoxelColor& voxelcolor)
{	
    Vector normal = end - start;
    Plane plane(start, normal);
    Vector base_1 = plane.base1();

    double length_base_1 = sqrt(base_1[0] * base_1[0] + base_1[1] * base_1[1] + base_1[2] * base_1[2]);  

    base_1 = base_1 / length_base_1 * (radius * sqrt(2) / 2);

    std::vector<Point> vecs;

    for (int i = 0; i < 4; i++)
    {
        double angle = (i + 0.5) * 2 * PI / 4;
        Vector v = RotationAxis(normal + base_1, angle, normal);
        vecs.push_back(start + v);
    }
    for (int i = 0; i < 4; i++)
    {
        vecs.push_back(vecs[i] - normal);
    }

    std::vector<std::vector<int>> faces;

    int face_index_0[4] = { 0, 1, 2, 3 };
    int face_index_1[4] = { 5, 1, 0, 4 };
    int face_index_2[4] = { 4, 0, 3, 7 };
    int face_index_3[4] = { 5, 4, 7, 6 };
    int face_index_4[4] = { 7, 3, 2, 6 };
    int face_index_5[4] = { 6, 2, 1, 5 };

    faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
    faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
    faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
    faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
    faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
    faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

    export_file_output << "g " + s_name << std::endl;

    for (int i = 0; i < vecs.size(); i++)
    {
        export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << voxelcolor.r << " " << voxelcolor.g << " " << voxelcolor.b << std::endl;
    }

    for (int i = 0; i < faces.size(); i++)
    {
        export_file_output << "f ";

        for (int j = 0; j < faces[i].size(); j++)
        {
            export_file_output << faces[i][j] + export_index << " ";
        }
        export_file_output << "" << std::endl;
    }

    export_index += 8;
}

void Bool2Voxel(const std::string& filename, const std::vector<bool>& h_, const VoxelColor& voxelcolor, const VoxelSize& voxelsize)
{
    std::ofstream file(filename);
    
    int index = 1;
    for (int ix = 0; ix < NX; ++ix)
        for (int iy = 0; iy < NY; ++iy)
            for (int iz = 0; iz < NZ; ++iz) {

                int idx = ix * NY * NZ + iy * NZ + iz;
                if (h_[idx]) {
                    Point voxelCenter(
                        (ix + 0.5) * voxelsize.x,
                        (iy + 0.5) * voxelsize.y,
                        (iz + 0.5) * voxelsize.z
                    );
                    Export_Voxel(file, index, "voxel" + std::to_string(ix) + "_" + std::to_string(iy) + "_" + std::to_string(iz), voxelCenter, voxelcolor, voxelsize);
                }
            }
    file.clear();
    file.close();
}

void Bool2Voxel(const std::string& filename, const std::vector<bool>& h_, const VoxelColor& voxelcolor, const VoxelSize& voxelsize, const std::vector<int>& bbox)
{
    std::ofstream file(filename);

    int index = 1;
    for (int ix = 0; ix < bbox[0]; ++ix)
        for (int iy = 0; iy < bbox[1]; ++iy)
            for (int iz = 0; iz < bbox[2]; ++iz) {

                int idx = ix * bbox[1] * bbox[2] + iy * bbox[2] + iz;
                if (h_[idx]) {
                    Point voxelCenter(
                        (ix + 0.5) * voxelsize.x,
                        (iy + 0.5) * voxelsize.y,
                        (iz + 0.5) * voxelsize.z
                    );
                    Export_Voxel(file, index, "voxel" + std::to_string(ix) + "_" + std::to_string(iy) + "_" + std::to_string(iz), voxelCenter, voxelcolor, voxelsize);
                }
            }
    file.clear();
    file.close();
}


VoxelColor getVoxelColor(double value) {
    const int color_index = 255 * value;
    double r, g, b;
    if (color_index < 32) {
        r = 0;
        g = 0;
        b = 0.5156 + 0.0156 * color_index;
    }
    else if (color_index < 96) {
        r = 0;
        g = 0.0156 + 0.9844 * (color_index - 32.0) / 64;
        b = 1;
    }
    else if (color_index < 158) {
        r = 0.0156 + (color_index - 96.0) / 64;
        g = 1;
        b = 0.9844 - (color_index - 96.0) / 64;
    }
    else if (color_index < 223) {
        r = 1;
        g = 1 - (color_index - 158.0) / 65;
        b = 0;
    }
    else {
        r = (2 - (color_index - 223.0) / 32) / 2.0;
        g = 0;
        b = 0;
    }
    return VoxelColor(static_cast<double>(r * 255),
        static_cast<double>(g * 255),
        static_cast<double>(b * 255));
}


void Int2VoxelUseColorBar(const std::string& filename, const std::vector<float>& h_, const VoxelSize& voxelsize)
{
    std::ofstream file(filename);
    float max_value = h_[0];
    float min_value = h_[0];
    for (int i = 0; i < TX; i++) {
        if (h_[i] > max_value && h_[i] > 0) {
            max_value = h_[i];
        }
        if (h_[i] < min_value && h_[i] > 0) {
            min_value = h_[i];
        }
    }

    VoxelColor voxelcolor;
    voxelcolor.r = voxelcolor.g = voxelcolor.b = 0;


    int index = 1;
    for (int ix = 0; ix < NX; ++ix)
        for (int iy = 0; iy < NY; ++iy)
            for (int iz = 0; iz < NZ; ++iz) {
                int idx = ix * NY * NZ + iy * NZ + iz;
                if (h_[idx] <= 0) {
                    continue;
                }
                const double normalized_value = (double)h_[idx] / (double)(max_value - min_value);
                voxelcolor = getVoxelColor(normalized_value);
                Point voxelCenter(
                    (ix + 0.5) * voxelsize.x,
                    (iy + 0.5) * voxelsize.y,
                    (iz + 0.5) * voxelsize.z
                );
                Export_Voxel(file, index, "voxel" + std::to_string(ix) + "_" + std::to_string(iy) + "_" + std::to_string(iz), voxelCenter, voxelcolor, voxelsize);
            }
    file.clear();
    file.close();
}


void Int2VoxelUseColorBar(const std::string& filename, const std::vector<int>& h_, const VoxelSize& voxelsize)
{
    std::ofstream file(filename);
    int max_value = h_[0];
    int min_value = h_[0];
    for (int i = 0; i < TX; i++) {
        if (h_[i] > max_value && h_[i]!= 0) {
            max_value = h_[i];
        }
        if (h_[i] < min_value && h_[i]!= 0) {
            min_value = h_[i];
        }
    } 

    VoxelColor voxelcolor;
    voxelcolor.r = voxelcolor.g = voxelcolor.b = 0;


    int index = 1;
    for (int ix = 0; ix < NX; ++ix)
        for (int iy = 0; iy < NY; ++iy)
            for (int iz = 0; iz < NZ; ++iz) {
                int idx = ix * NY * NZ + iy * NZ + iz;
                if (h_[idx] == 0) {
                    continue;
                }
                const double normalized_value = (double)h_[idx] / (double)(max_value - min_value);
                voxelcolor = getVoxelColor(normalized_value);
                Point voxelCenter(
                    (ix + 0.5) * voxelsize.x,
                    (iy + 0.5) * voxelsize.y,
                    (iz + 0.5) * voxelsize.z
                );
                Export_Voxel(file, index, "voxel" + std::to_string(ix) + "_" + std::to_string(iy) + "_" + std::to_string(iz), voxelCenter, voxelcolor, voxelsize);
            }
    file.clear();
    file.close();
}


void translateMesh(Polyhedron& mesh, const Vector& translation)
{
    std::cout << translation.x() << " " << translation.y() << " " << translation.z() << std::endl;
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);
    tree.accelerate_distance_queries();

    CGAL::Bbox_3 bbox = tree.bbox();

    Vector v_min(bbox.xmin(), bbox.ymin(), bbox.zmin());

    for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
        Point p = v->point();
        p = Point(p.x() - v_min.x() + translation.x() * VOXELSIZE,
            p.y() - v_min.y() + translation.y() * VOXELSIZE,
            p.z() - v_min.z() + translation.z() * VOXELSIZE);
        v->point() = p;
    }
}

void translateMesh(trimesh::TriMesh& mesh, const Vector& translation)
{
    mesh.need_bbox();

    for (auto v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        trimesh::point vertex = *v;
        Point p(vertex.x, vertex.y, vertex.z);
        p = Point(p.x() - mesh.bbox.min.x + translation.x() * VOXELSIZE,
            p.y() - mesh.bbox.min.y + translation.y() * VOXELSIZE,
            p.z() - mesh.bbox.min.z + translation.z() * VOXELSIZE);
        vertex.x = p.x();
        vertex.y = p.y();
        vertex.z = p.z();
        *v = vertex;
    }
}


void rotateMesh(Polyhedron& mesh, const Vector& ori)
{
    Eigen::Vector3d vector_Before(ori.x(), ori.y(), ori.z());
    Eigen::Vector3d vector_After(0, 1, 0);
    Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_Before, vector_After).toRotationMatrix();
     
    for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
        Point p = v->point();
        Eigen::MatrixXd temp_V;
        temp_V.resize(3, 1);
        temp_V(0, 0) = p.x();
        temp_V(1, 0) = p.y();
        temp_V(2, 0) = p.z();

        temp_V = rotMatrix.inverse() * temp_V;

        p = Point(temp_V(0, 0), temp_V(1, 0), temp_V(2, 0));
        v->point() = p;
    }
}

void rotateMesh(trimesh::TriMesh& mesh, const Vector& ori)
{
    Eigen::Vector3d vector_Before(ori.x(), ori.y(), ori.z());
    Eigen::Vector3d vector_After(0, 1, 0);
    Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_Before, vector_After).toRotationMatrix();

    for (auto v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        trimesh::point vertex = *v;
        Point p(vertex.x, vertex.y, vertex.z);
        Eigen::MatrixXd temp_V;
        temp_V.resize(3, 1);
        temp_V(0, 0) = p.x();
        temp_V(1, 0) = p.y();
        temp_V(2, 0) = p.z();

        temp_V = rotMatrix.inverse() * temp_V;

        p = Point(temp_V(0, 0), temp_V(1, 0), temp_V(2, 0));
        vertex.x = p.x();
        vertex.y = p.y();
        vertex.z = p.z();
        *v = vertex;
    }
}

void scaleMesh(Polyhedron& mesh, const double& scale)
{
    for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
        Point p = v->point();
        p = Point(p.x()*scale, p.y()*scale, p.z()*scale);
        v->point() = p;
    }
}

void scaleMesh(Polyhedron& mesh, const Vector& scale)
{
    for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
        Point p = v->point();
        p = Point(p.x() * scale.x(), p.y() * scale.y(), p.z() * scale.z());
        v->point() = p;
    }
}

void scaleMesh(trimesh::TriMesh& mesh, const double& scale)
{
    for (auto v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        trimesh::point vertex = *v;
        Point p(vertex.x, vertex.y, vertex.z);
       
        p = Point(p.x() * scale, p.y() * scale, p.z() * scale);
        vertex.x = p.x();
        vertex.y = p.y();
        vertex.z = p.z();
        *v = vertex;
    }
}



void compute_intersection(std::vector<Polyhedron>& meshes, Polyhedron& intersection_result) {
    if (meshes.empty()) return;

    intersection_result = meshes[0];  

    for (size_t i = 1; i < meshes.size(); ++i) {
        Polyhedron temp_result;
        PMP::corefine_and_compute_intersection(intersection_result, meshes[i], temp_result);
        intersection_result = temp_result;  

    }
}



double computeDistance(const Point& p1, const Point& p2) {
    return std::sqrt(CGAL::squared_distance(p1, p2));
}

std::vector<std::pair<int, int>> computeEdges(const std::vector<Point>& points, double epsilon) {
    size_t n = points.size();
    std::vector<std::pair<int, int>> edges;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (computeDistance(points[i], points[j]) < epsilon) {
                edges.push_back({ i, j });
            }
        }
    }

    return edges;
}


Point computeCentroid(const std::vector<Point>& points) {
    double cx = 0, cy = 0, cz = 0;
    for (const auto& p : points) {
        cx += p.x();
        cy += p.y();
        cz += p.z();
    }
    double n = static_cast<double>(points.size());
    return Point(cx / n, cy / n, cz / n);
}

std::vector<Point> generateOffset(const std::vector<Point>& points, const Vector& offsetDirection, double offsetDistance) {
    std::vector<Point> offsetPoints;
    Vector normalizedOffset = offsetDirection / std::sqrt(offsetDirection.squared_length());
    Vector offsetVector = normalizedOffset * offsetDistance;

    for (const auto& p : points) {
        offsetPoints.push_back(p + offsetVector);
    }

    return offsetPoints;
}

void updateNeighbor(int index_0, int index_1, std::vector<std::vector<int>>& vecs_neigbor, std::vector<std::vector<int>>& vecs_neigbor_lable) {
    int search_0 = -1;
    for (int j = 0; j < vecs_neigbor[index_0].size() && search_0 < 0; j++) {
        if (vecs_neigbor[index_0][j] == index_1) {
            search_0 = j;
            vecs_neigbor_lable[index_0][j]++;
        }
    }
    if (search_0 < 0) {
        vecs_neigbor[index_0].push_back(index_1);
        vecs_neigbor_lable[index_0].push_back(1);
    }
}

void computeBoundary(const std::vector<Point>& vecs, const std::vector<std::pair<int,int>> edges, std::set<int>& boundaries)
{

    std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
    std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());

    for (const auto& edge : edges) {
        updateNeighbor(edge.first, edge.second, vecs_neigbor, vecs_neigbor_lable);
        updateNeighbor(edge.second, edge.first, vecs_neigbor, vecs_neigbor_lable);
    }

    std::set<std::pair<int, int>> deleted_edges;

    for (int i = 0; i < vecs.size(); i++) {
        for (int j = 0; j < vecs_neigbor_lable[i].size(); j++) {
            if (vecs_neigbor_lable[i][j] == 1) {     
                int index_1 = vecs_neigbor[i][j];
                if (deleted_edges.count({ i, index_1 }) == 0) {
                    boundaries.insert(i);
                    boundaries.insert(index_1);
                    deleted_edges.insert({ i, index_1 });
                    deleted_edges.insert({ index_1, i });
                }
            }
        }
    }
    
}

void samplePointsAround(std::vector<Weighted_point_2>& sampled_points, const Weighted_point_2& center, double radius, int num_samples) {
    double angle_increment = 2 * PI / num_samples;  
    for (int i = 0; i < num_samples; ++i) {
        double angle = i * angle_increment;  

        double x = center.x() + radius * std::cos(angle);   
        double y = center.y() + radius * std::sin(angle);   
        Weighted_point_2 p(x, y);
        sampled_points.push_back(p);
    }
}


void generateMesh(const std::vector<Point>& original, const Vector& offsetDirection,const double offsetLength, Polyhedron& mesh) {
   
    std::vector<Point> all_points_;
    for(const auto& p : original)
        all_points_.push_back(p);

    Eigen::Vector3d normal(offsetDirection.x(), offsetDirection.y(), offsetDirection.z());
    normal.normalize();
    Eigen::Vector3d vectorBefore(normal.x(), normal.y(), normal.z());
    Eigen::Vector3d vectorAfter(0, 0, 1);
    Eigen::Matrix3d rotMatrix3;
    rotMatrix3 = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();


    std::vector<Eigen::MatrixXd> temp_V2;
    temp_V2.resize(all_points_.size());
    for (int i = 0; i < all_points_.size(); i++) {
        temp_V2[i].resize(3, 1);
        temp_V2[i](0, 0) = all_points_[i].x();
        temp_V2[i](1, 0) = all_points_[i].y();
        temp_V2[i](2, 0) = all_points_[i].z();
    }
    for (int i = 0; i < all_points_.size(); i++)
        temp_V2[i] = rotMatrix3.inverse() * temp_V2[i];

    std::vector<Weighted_point_2> initial_points;
    for (int i = 0; i < all_points_.size(); i++)
        initial_points.push_back(Weighted_point_2(temp_V2[i](0, 0), temp_V2[i](1, 0)));

    std::vector<Weighted_point_2> all_sampled_points;

    for (int i = 0; i < initial_points.size(); i++) {
        std::vector<Weighted_point_2> sampled_points;
        samplePointsAround(sampled_points, initial_points[i], 0.5 * VOXELSIZE, 4);
        for (int j = 0; j < sampled_points.size(); j++) {
            all_sampled_points.push_back(sampled_points[j]);
        }
    }
    for (int i = 0; i < all_sampled_points.size(); i++) {
        initial_points.push_back(all_sampled_points[i]);
    }



    Alpha_shape_2 alpha_shape(initial_points.begin(), initial_points.end());    

    Alpha_iterator_2 opt = alpha_shape.find_optimal_alpha(1);
    std::cout << "Optimal alpha value to get one connected component is " << *opt << std::endl;

    alpha_shape.set_alpha((*opt));

    assert(alpha_shape.number_of_solid_components() == 1);

    std::vector<Point2> points;

    for (auto eit = alpha_shape.alpha_shape_edges_begin(); eit != alpha_shape.alpha_shape_edges_end(); ++eit) {
        auto seg = alpha_shape.segment(*eit);
        Point2 p1(seg.source().x(), seg.source().y());
        Point2 p2(seg.target().x(), seg.target().y());
        points.push_back(p1);
        points.push_back(p2);
    }
    
    std::cout << "Number of points in the alpha shape: " << points.size() << std::endl;
    double sumX = 0, sumY = 0;
    for (const auto& p : points) {
        sumX += p.x();
        sumY += p.y();
    }
    Point2 centroid(sumX / points.size(), sumY / points.size());
    std::sort(points.begin(), points.end(), [&centroid](const Point2& p1, const Point2& p2) {
        return std::atan2(p1.y() - centroid.y(), p1.x() - centroid.x()) < std::atan2(p2.y() - centroid.y(), p2.x() - centroid.x());
        });
    double tmp_V2_val = temp_V2[0](2, 0);
    std::vector<Point> points3D, points3D_2;
    for (int i = 0; i < points.size(); i++) {
        points3D.push_back(Point(points[i].x(), points[i].y(), tmp_V2_val));
        points3D_2.push_back(Point(points[i].x(), points[i].y(), tmp_V2_val + offsetLength));
    }

    temp_V2.resize(points.size());
    rotMatrix3 = Eigen::Quaterniond::FromTwoVectors(vectorAfter, vectorBefore).toRotationMatrix();
    for (int i = 0; i < points3D.size(); i++) {
        temp_V2[i].resize(3, 1);
        temp_V2[i](0, 0) = points3D[i].x();
        temp_V2[i](1, 0) = points3D[i].y();
        temp_V2[i](2, 0) = points3D[i].z();
    }
    for (int i = 0; i < points3D.size(); i++)
        temp_V2[i] = rotMatrix3.inverse() * temp_V2[i];
    for (int i = 0; i < points3D.size(); i++)
        points3D[i] = Point(temp_V2[i](0, 0), temp_V2[i](1, 0), temp_V2[i](2, 0));

    for (int i = 0; i < points3D_2.size(); i++) {
        temp_V2[i].resize(3, 1);
        temp_V2[i](0, 0) = points3D_2[i].x();
        temp_V2[i](1, 0) = points3D_2[i].y();
        temp_V2[i](2, 0) = points3D_2[i].z();
    }
    for (int i = 0; i < points3D.size(); i++)
        temp_V2[i] = rotMatrix3.inverse() * temp_V2[i];
    for (int i = 0; i < points3D_2.size(); i++)
        points3D_2[i] = Point(temp_V2[i](0, 0), temp_V2[i](1, 0), temp_V2[i](2, 0));


    std::string filename = "output/RESULTS/bunny/";
    std::ofstream file(filename + "aa.obj");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
    }

    int n = points3D.size();   
    std::vector<Point> vertices;
    std::vector<std::vector<int>> faces;

    vertices.insert(vertices.end(), points3D.begin(), points3D.end());
    vertices.insert(vertices.end(), points3D_2.begin(), points3D_2.end());

    for (const auto& vertex : vertices) {
        file << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
    }

    using Coord = double;
    using NN = uint32_t;
    using PPoint = std::array<Coord, 2>;
    std::vector<std::vector<PPoint>> polygon(1);
    std::vector<NN> indices = mapbox::earcut<NN>(polygon);
    polygon[0].clear();
    for (int i = 0; i < points.size(); i++) {
        polygon[0].push_back({ points[i].x(), points[i].y()});
    }

    std::vector<NN> indicess = mapbox::earcut<NN>(polygon);
    for (int j = indicess.size()-1; j >=0;) {
        file << "f";
        for (int k = 0; k < 3; k++) {
            file << " " << indicess[j] + 1;
            j--;
        }
        file << std::endl;
    }
    for (int j = 0; j < indicess.size();) {
        file << "f";
        for (int k = 0; k < 3; k++) {
            file << " " << indicess[j] + 1 + points.size();
            j++;
        }
        file << std::endl;
    }

    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        faces.push_back({ i, next, next + n, i + n });   
    }
    
    for (const auto& face : faces) {
        file << "f";
        for (int idx : face) {
            file << " " << idx + 1;   
        }
        file << "\n";
    }
    file.close();
}
