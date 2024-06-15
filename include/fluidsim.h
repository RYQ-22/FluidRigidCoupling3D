#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "util.h"
#include "rigidbody.h"

namespace backend {

class FluidSim{
private:
    double l;
    int n1, n2, n3;  
    int n_liquid = 0;
    // velocity at the cell's center
    Field3d u, v, w;
    Field3d u_new, v_new, w_new;
    // level set
    Field3d phi;
    // boundary
    Field3d phi_solid;
    Field3i valid_u, valid_v, valid_w;
    Field3i valid_u_old, valid_v_old, valid_w_old;

    Field3d u_solid, v_solid, w_solid;
    // for projection
    Field3d weights_u, weights_v, weights_w;
    Field3d Adiag, Aplusi, Aplusj, Aplusk;
    Field3d d, p;
    Field3d z, s; // z: auxiliary vetor, s: search vector

    // particles
    double particle_radius, particle_output_radius;
    std::vector<Vector3d> particles;
    std::vector<Vector3d> particles_output;
    Field3i particles_num;

    // rigidbody
    bool add_rigidbody = false;
    RigidBody rigidbody;
    double rigidbody_density;
    double fluid_density;

    Field3d phi_rigidbody;
    Field3d phi_solid_rigidbody;

    Field3d J_x, J_y, J_z, J_rot_x, J_rot_y, J_rot_z;

    Field3d weights_rigid_u, weights_rigid_v, weights_rigid_w;

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
    void computeWeights();        
    void solve(int maxIterations);
    void applyA(const Field3d& x, Field3d& ans);
    Field3d applyA(const Field3d& x);
    // some details
    void setParticlesOutput();
    void computePhi();
    void extrapolate();
    void extrapolate(Field3d& field, Field3d field_new, Field3i& valid, Field3i& valid_old);
    void constrain();
    void computeN();
    // rigidbody
    void updateRigidBodyGrids();// update phi_rigidbody, phi_solid_rigidbody, weights_rigid
    void recomputeRigidBodyMass();
    void recomputeSolidVelocity();

public:
    // init
    FluidSim(const int& n1_init, const int& n2_init, const int& n3_init, const double& l_init, const Field3d& phi_init, const Field3d& phi_solid_init);
    void advance(const double& time_step);
    void setVelocity(const Vector3d& vec);

    Vector3d getVelocity(const Vector3d& pos);
    Vector3d getSolidVelocity(const Vector3d& pos);
    double getPhi(const Vector3d& pos);
    double getPhiSolid(const Vector3d& pos);
    double getPhiRigidBody(const Vector3d& pos);
    double getPhiSolidRigidBody(const Vector3d& pos);

    // rigid body
    void addRigidBody(
        const std::vector<Vector3d>& vertices_init,
        const std::vector<Face>& faces_init,
        const double& m_init, 
        const Matrix3d& I_init, 
        const Vector3d& c_init, 
        const Quaterniond& q_init, 
        const Vector3d& velocity_init, 
        const Vector3d& omega_init,
        const FuncPtr& get_sdf_func_init
    );

    void run(double time_step);// no UI    
    void outputPLY(int fps, int t, double dt, int n);
    void outputVedio(const std::string& file_name, int fps, int t);
};

}

#endif