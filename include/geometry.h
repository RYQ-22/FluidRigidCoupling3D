#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "util.h"

namespace backend {

struct Face {
    int v0;
    int v1;	
    int v2;
    Face() {}
    Face(int v0_, int v1_, int v2_) :v0(v0_), v1(v1_), v2(v2_) {}
};

Field3d computeSDF_cuboid(const int& n1, const int& n2, const int& n3, const Vector3i& left_corner, const Vector3i& cuboid_size);
double computeSDF_cuboid(const Vector3d& pos, const Vector3d& cuboid_size);
Field3d computeSDF_sphere(const int n1, const int& n2, const int& n3, const Vector3d& center, const double& r);

void writePLY(const std::string& ply_path, const std::vector<Vector3d>& vertices, const std::vector<Face>& faces);
void writePLY(const std::string& ply_path, const std::vector<Vector3d>& particles);
void writeNPY(const std::string& npy_path, const Field3d& phi);
void writeVideo(const std::string& file_name, const int& fps, const int& t);

Field3d objToSDF(int N1, int N2, int N3, int size, double l, std::string obj_path);
bool objToSDF(int N1, int N2, int N3, double l, std::string obj_path, std::vector<double>& phi, double scale = 1.0, double translate_x = 0.0, double translate_y = 0.0, double translate_z = 0.0, bool inverse = false);

}

#endif