#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "Eigen/Dense"
#include "field.h"
#include "util.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

using VectorXd = Eigen::VectorXd;
using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using Field3d = Field3<double>;
using Field3i = Field3<int>;

using appconstants::PI;
using appconstants::G;

class FluidSim{
private:
    double l;
    int n1, n2, n3;  
    // velocity at the cell's center
    Field3d u, v, w;
    Field3d u_new, v_new, w_new;
    // level set
    Field3d phi;
    // boundary
    Field3d phi_solid;
    Field3i valid_u, valid_v, valid_w;
    Field3i valid_u_old, valid_v_old, valid_w_old;
    // for projection
    Field3d Adiag, Aplusi, Aplusj, Aplusk;
    Field3d d, p;
    Field3d z, s; // z: auxiliary vetor, s: search vector

    // particles
    double particle_radius;
    std::vector<Vector3d> particles;
    
    double computeVolume(const int& i, const int& j, const int& k, const Field3d& phi, const int& id);
    // advance
    // 1. add force
    void applyForce(const double& dt);
    // 2. advect
    void advect(const double& dt);
    void advectParticles(const double& dt);
    Vector3d traceRk2(const Vector3d& position, const double& dt);
    void semiLagrangian(const Field3d& field, Field3d& field_new, const double& dt, const int& id);
    // 3. project
    void project();
    void solve(const int& maxIterations);
    void applyA(const Field3d& x, Field3d& ans);
    // some details
    void computePhi();
    void extrapolate();
    void extrapolate(Field3d& field, Field3d field_new, Field3i& valid, Field3i& valid_old);
    void constrain();


public:
    // init
    FluidSim(const int& n1_init, const int& n2_init, const int& n3_init, const double& l_init, const std::vector<double>& phi_init, const std::vector<double>& phi_solid_init);
    void advance(const double& dt);
    void setVelocity(const Vector3d& v);
    bool valid(const int& i, const int& j, const int& k);
    Vector3d getVelocity(const Vector3d& position);
};

#endif