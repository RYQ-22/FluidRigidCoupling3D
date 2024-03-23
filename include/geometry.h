#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>

inline const std::string projectPath = "/home/ricky/Documents/MyProject/FluidRigidCoupling3D";

bool objToSDF(int N1, int N2, int N3, int size, double l, std::string obj_path, std::vector<double>& phi, bool inverse = false);
bool objToSDF(int N1, int N2, int N3, double l, std::string obj_path, std::vector<double>& phi, double scale = 1.0, double translate_x = 0.0, double translate_y = 0.0, double translate_z = 0.0, bool inverse = false);

#endif