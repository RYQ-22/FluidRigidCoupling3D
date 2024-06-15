#include "rigidbody.h"

namespace backend{

Matrix3d ToCrossProduct(const Vector3d& vec) {
    Matrix3d vec_cross; vec_cross.setZero();
    vec_cross(0, 1) = -vec(2); vec_cross(0, 2) = vec(1);
    vec_cross(1, 0) = vec(2); vec_cross(1, 2) = -vec(0);
    vec_cross(2, 0) = -vec(1); vec_cross(2, 1) = vec(0);
    return vec_cross;
}

RigidBody::RigidBody() {}

RigidBody::RigidBody(
        const std::vector<Vector3d>& vertices_init,
        const std::vector<Face>& faces_init,
        const double& m_init, 
        const Matrix3d& I_init, 
        const Vector3d& c_init, 
        const Quaterniond& q_init, 
        const Vector3d& velocity_init, 
        const Vector3d& omega_init,
        const FuncPtr& get_sdf_func_init
    ) : vertices(vertices_init), 
        faces(faces_init), 
        c(c_init),
        velocity(velocity_init),
        I_rel(I_init),
        I_rel_inv(I_init.inverse()),
        q(q_init),
        omega(omega_init),
        get_sdf_func(get_sdf_func_init) {
            M.setZero(); M(0, 0) = m_init; M(1, 1) = m_init; M(2, 2) = m_init;
            M_inv = M.inverse();
            L = getI() * omega;
            vertices_temp.resize(vertices.size());            
            for (int idx = 0; idx < vertices.size(); idx++) {
                vertices_temp[idx] = getVertexPosition(idx);
            }

    }

Matrix3d RigidBody::getRotationMatrix() {
    return q.toRotationMatrix();
}
    
Matrix3d RigidBody::getMass() {
    return M;
}
    
Matrix3d RigidBody::getI() {
    return getRotationMatrix() * I_rel * getRotationMatrix().transpose();
}

Matrix3d RigidBody::getI_inv() {
    return getRotationMatrix() * I_rel_inv * getRotationMatrix().transpose();
}
    
Vector3d RigidBody::getC() {
    return c;
}

Vector3d RigidBody::RigidBody::getVelocity() {
    return velocity;
}
    
Vector3d RigidBody::getVelocity(const Vector3d& pos) {
    Vector3d r_rel = pos - c;
    return velocity + omega.cross(r_rel);
}
    
Vector3d RigidBody::getOmega() {
    return omega;
}
    
Vector3d RigidBody::getVertexPosition(const int& idx) {
    return c + getRotationMatrix() * vertices[idx];
}

std::vector<Vector3d> RigidBody::getVertices() {
    return vertices_temp;
}

std::vector<Face> RigidBody::getFaces() {
    return faces;
}

Vector3d RigidBody::Force() {
    return Vector3d(0, -G * M(1, 1), 0) -0.1 * M * velocity;
}

Vector3d RigidBody::Torque() {
    return Vector3d::Zero();
}

double RigidBody::getSDF(const Vector3d& pos) {
    Vector3d pos_rel = getRotationMatrix().transpose() * (pos - c);
    return get_sdf_func(pos_rel);
}

void RigidBody::setVelocity(const Vector3d& velocity_) {
    velocity = velocity_;
    return;
}

void RigidBody::setOmega(const Vector3d& omega_) {
    omega = omega_;
    return;
}

void RigidBody::setMass(const double& m_u, const double& m_v, const double& m_w) {
    M.setZero();
    M(0, 0) = m_u;
    M(1, 1) = m_v;
    M(2, 2) = m_w;
    return;
}

void RigidBody::handleCollision(const Field3d& phi_solid, const double& l) {
    // handle collision
    Vector3d pos, grad;
    Vector3d vi, vi_n, vi_t;
    Matrix3d K;
    Vector3d r_rel, j, vi_new;
    for (int idx = 0; idx < vertices.size(); idx++) {        
        pos(0) = clamp(getVertexPosition(idx)(0), 0., (phi_solid.getN1()-1-1e-6) * l);
        pos(1) = clamp(getVertexPosition(idx)(1), 0., (phi_solid.getN2()-1-1e-6) * l);
        pos(2) = clamp(getVertexPosition(idx)(2), 0., (phi_solid.getN3()-1-1e-6) * l);
        if (interpolate_value(pos / l, phi_solid) < 0) {
            interpolate_gradient(grad, pos / l, phi_solid);
            grad.normalize();
            vi = getVelocity(pos);
            if (vi.dot(grad) < 0) {// v_n < 0
                vi_n = vi.dot(grad) * grad;
                vi_t = vi - vi_n;
                r_rel = getRotationMatrix() * vertices[idx];
                K = M_inv - ToCrossProduct(r_rel) * getI_inv() * ToCrossProduct(r_rel);
                vi_new = -e * vi_n + max(0., 1. - mu * (1. + e) * vi_n.norm() / vi_t.norm()) * vi_t;
                j = K.inverse() * (vi_new - vi);

                // update velocity and omega
                velocity += M_inv * j;
                L += r_rel.cross(j);
                omega = getI_inv() * L;
            }
        }
    }

    // project rigidbody out of solid
    double depth;
    for (int idx = 0; idx < vertices.size(); idx++) {
        pos(0) = clamp(getVertexPosition(idx)(0), 0., (phi_solid.getN1()-1-1e-6) * l);
        pos(1) = clamp(getVertexPosition(idx)(1), 0., (phi_solid.getN2()-1-1e-6) * l);
        pos(2) = clamp(getVertexPosition(idx)(2), 0., (phi_solid.getN3()-1-1e-6) * l);
        depth = interpolate_value(pos / l, phi_solid) * l;
        if (depth < 0) {
            interpolate_gradient(grad, pos / l, phi_solid);
            grad.normalize();
            c(0) += -depth * grad(0);
            c(1) += -depth * grad(1);
            c(2) += -depth * grad(2);
        }
    }

    return;
}

void RigidBody::advance(const double& dt, const Field3d& phi_solid, const double& l) {
    velocity += M_inv * Force() * dt;
    L += Torque() * dt;
    c += velocity * dt;
    q = Quaterniond(1., 0.5 * omega(0) * dt, 0.5 * omega(1) * dt, 0.5 * omega(2) * dt) * q;
    omega = getI_inv() * L;

    handleCollision(phi_solid, l);

    // update vertices position
    for (int idx = 0; idx < vertices.size(); idx++) {
        vertices_temp[idx] = getVertexPosition(idx);
    }

    return;
}



}
