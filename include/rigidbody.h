#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "geometry.h"

namespace backend {

class RigidBody {
private:
    std::vector<Vector3d> vertices;// (at theta = 0, in S_c)
    std::vector<Vector3d> vertices_temp;
    std::vector<Face> faces;    
    FuncPtr get_sdf_func;// pointer of the function to compute sdf in S_c
    
    // translation
    Matrix3d M;
    Matrix3d M_inv;
    Vector3d c;    
    Vector3d velocity;
    // rotation
    Matrix3d I_rel;
    Matrix3d I_rel_inv;
    Quaterniond q;
    Vector3d omega;
    Vector3d L;

    Vector3d Force();
    Vector3d Torque();

    // for collision
    double e = 0.3f;
    double mu = 0.8f;

    void handleCollision(const Field3d& phi_solid, const double& l);

public:
    RigidBody();
    RigidBody(
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

    Matrix3d getRotationMatrix();
    Matrix3d getMass();
    Matrix3d getI();
    Matrix3d getI_inv();
    Vector3d getC();
    Vector3d getVelocity();
    Vector3d getVelocity(const Vector3d& pos);
    Vector3d getOmega();
    Vector3d getVertexPosition(const int& idx);    
    std::vector<Vector3d> getVertices();
    std::vector<Face> getFaces();
    double getSDF(const Vector3d& pos);

    void setVelocity(const Vector3d& velocity_);
    void setOmega(const Vector3d& omega_);
    void setMass(const double& m_u, const double& m_v, const double& m_w);

    void advance(const double& time_step, const Field3d& phi_solid, const double& l);    

};
    
}


#endif