#include "fluidsim.h"
#include "geometry.h"

#include <iostream>

int main(){
    int n1 = 40;
    int n2 = 40;
    int n3 = 80;
    double l = 0.01f;
    std::vector<double> phi;
    std::vector<double> phi_solid;
    std::string cube_path = projectPath + "/OBJ/test.obj";
    std::string cuboid_path = projectPath + "/OBJ/cuboid.obj";
    std::string boundary_path = projectPath + "/OBJ/test2.obj";
    objToSDF(n1+1, n2+1, n3+1, 36, l, boundary_path, phi_solid, true);
    objToSDF(n1, n2, n3, l, cuboid_path, phi, 0.16f, 0.2f, 0.2f, 0.64, false);
    FluidSim sim(n1, n2, n3, l, phi, phi_solid);
    sim.run(0.001, 1);

    return 0;
}