#ifndef CONFIG_H
#define CONFIG_H

// Commonly used std headers.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>
//#include <omp.h>

// Eigen headers.
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "unsupported/Eigen/Polynomials"
#include "unsupported/Eigen/MatrixFunctions"

// openCV
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// tinyply
#include "tinyply.h"

// to write .npy
#include "cnpy.h"

namespace backend {

// Define constant
const double PI = 3.1415926535897932384626433832795;
const double G = 9.8;

// Alias
using VectorXd = Eigen::VectorXd;
using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using Matrix3d = Eigen::Matrix3d;
using Vector8d = Eigen::Matrix<double, 8, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Quaterniond = Eigen::Quaterniond;

using FuncPtr = double (*)(const Vector3d&);

// Path
const std::string projectPath = std::string("/home/ricky/Documents/MyProject/FluidRigidCoupling3D");
const std::string pythonPath = std::string("/home/ricky/anaconda3/envs/myenv/bin/python");
const std::string pbrtPath = std::string("/home/ricky/Documents/MyProject/pbrt-v4/build");

}

#endif