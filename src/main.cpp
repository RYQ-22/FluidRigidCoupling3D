#include "fluidsim.h"
#include "geometry.h"

int n1 = 60;
int n2 = 60;
int n3 = 120;
double l = 1. / n3;

double rigidbodySDF(const backend::Vector3d& pos) {
    return backend::computeSDF_cuboid(pos, backend::Vector3d(static_cast<double>(n1) / 4, static_cast<double>(n1) / 4, static_cast<double>(n1) / 4) * l);
}

int main(){   
    backend::Field3d phi_fluid;
    backend::Field3d phi_solid;
    //std::string cube_path = backend::projectPath + "/OBJ/test.obj";
    //std::string cuboid_path = backend::projectPath + "/OBJ/cuboid.obj";
    //std::string boundary_path = backend::projectPath + "/OBJ/test2.obj";
    //backend::objToSDF(n1+1, n2+1, n3+1, 37, l, boundary_path, phi_solid, true);
    //backend::objToSDF(n1, n2, n3, l, cuboid_path, phi, 0.16f, 0.2f, 0.2f, 0.64, false);
    phi_solid = -backend::computeSDF_cuboid(n1, n2, n3, backend::Vector3i(2, 2, 2), backend::Vector3i(n1 - 4, n2 - 4, n3 - 4)) * l;
    //phi_fluid = backend::computeSDF_sphere(n1, n2, n3, backend::Vector3d((double)n1 / 2, (double)n2 / 2, (double)n2 / 2), 20) * l;
    //phi_fluid = backend::computeSDF_cuboid(n1, n2, n3, backend::Vector3i(2, 2, 2), backend::Vector3i(n1 - 4, 22, n3 - 4)) * l;
    //phi_fluid = backend::computeSDF_cuboid(n1, n2, n3, backend::Vector3i(10, 5, 15), backend::Vector3i(30, 48, 30)) * l;
    phi_fluid = backend::objToSDF(n1+1, n2+1, n3+1, 52, l, backend::projectPath + "/OBJ/bunny.obj");
    backend::FluidSim sim(n1, n2, n3, l, phi_fluid, phi_solid);

    // add rigidbody
    double a = n1 * l / 4;
    std::vector<backend::Vector3d> vertices;
    std::vector<backend::Face> faces;    
    vertices.push_back(backend::Vector3d(-a / 2, -a / 2, -a / 2));
    vertices.push_back(backend::Vector3d(-a / 2, -a / 2, a / 2));
    vertices.push_back(backend::Vector3d(-a / 2, a / 2, -a / 2));
    vertices.push_back(backend::Vector3d(-a / 2, a / 2, a / 2));
    vertices.push_back(backend::Vector3d(a / 2, a / 2, -a / 2));
    vertices.push_back(backend::Vector3d(a / 2, a / 2, a / 2));
    vertices.push_back(backend::Vector3d(a / 2, -a / 2, -a / 2));
    vertices.push_back(backend::Vector3d(a / 2, -a / 2, a / 2));
    
    faces.push_back(backend::Face(0, 1, 3));
    faces.push_back(backend::Face(0, 3, 2));
    faces.push_back(backend::Face(2, 5, 4));
    faces.push_back(backend::Face(2, 3, 5));
    faces.push_back(backend::Face(4, 5, 7));
    faces.push_back(backend::Face(4, 7, 6));
    faces.push_back(backend::Face(0, 7, 1));
    faces.push_back(backend::Face(0, 6, 7));
    faces.push_back(backend::Face(1, 5, 3));
    faces.push_back(backend::Face(1, 7, 5));
    faces.push_back(backend::Face(0, 4, 6));
    faces.push_back(backend::Face(0, 2, 4));

    double m_init = 1.0;
    backend::Matrix3d I_init; I_init.setZero();
    I_init(0, 0) = m_init * a * a / 3;
    I_init(1, 1) = m_init * a * a / 3;
    I_init(2, 2) = m_init * a * a / 3;
    backend::Vector3d c_init(n1 * l / 2, n2 * l * 0.75, n3 * l * 0.71);
    backend::Quaterniond q_init(cos(0.3), sin(0.3)/sqrt(3), sin(0.3)/sqrt(3), sin(0.3)/sqrt(3));
    backend::Vector3d velocity_init(0., -0.4, -0.);
    backend::Vector3d omega_init; omega_init.setZero();

    sim.addRigidBody(
        vertices,
        faces,
        m_init,
        I_init,
        c_init,
        q_init,
        velocity_init,
        omega_init,
        rigidbodySDF
    );

    //sim.run(0.005);
    int fps = 30;
    int t = 3;
    sim.outputPLY(fps, t, 0.01, 1);
    sim.outputVedio("output10", fps, t);
    return 0;
}