#include "fluidsim.h"

namespace backend {

FluidSim::FluidSim(
    const int& n1_init,
    const int& n2_init,
    const int& n3_init,
    const double& l_init,
    const Field3d& phi_init,
    const Field3d& phi_solid_init) {
        Assert(n1_init+1 == phi_init.getN1() && n2_init+1 == phi_init.getN2() && n3_init+1 == phi_init.getN3(),
        "FluidSim::FluidSim", "The shape of phi_init is not correct.");
        Assert(n1_init+1 == phi_solid_init.getN1() && n2_init+1 == phi_solid_init.getN2() && n3_init+1 == phi_solid_init.getN3(),
        "FluidSim::FluidSim", "The shape of phi_solid_init is not correct.");    
        n1 = n1_init;
        n2 = n2_init;
        n3 = n3_init;
        l = l_init;  
        u.resize3D(n1+1, n2, n3); u_new.resize3D(n1+1, n2, n3); valid_u.resize3D(n1+1, n2, n3); valid_u_old.resize3D(n1+1, n2, n3); weights_u.resize3D(n1+1, n2, n3); u_solid.resize3D(n1+1, n2, n3);
        u.setZero(); u_new.setZero(); valid_u.setZero(); valid_u_old.setZero(); weights_u.setZero(); u_solid.setZero();
        v.resize3D(n1, n2+1, n3); v_new.resize3D(n1, n2+1, n3); valid_v.resize3D(n1, n2+1, n3); valid_v_old.resize3D(n1, n2+1, n3); weights_v.resize3D(n1, n2+1, n3); v_solid.resize3D(n1, n2+1, n3);
        v.setZero(); v_new.setZero(); valid_v.setZero(); valid_v_old.setZero(); weights_v.setZero(); v_solid.setZero();
        w.resize3D(n1, n2, n3+1); w_new.resize3D(n1, n2, n3+1); valid_w.resize3D(n1, n2, n3+1); valid_w_old.resize3D(n1, n2, n3+1); weights_w.resize3D(n1, n2, n3+1); w_solid.resize3D(n1, n2, n3+1);
        w.setZero(); w_new.setZero(); valid_w.setZero(); valid_w_old.setZero(); weights_w.setZero(); w_solid.setZero();
        Adiag.resize3D(n1, n2, n3); Aplusi.resize3D(n1, n2, n3); Aplusj.resize3D(n1, n2, n3); Aplusk.resize3D(n1, n2, n3);
        d.resize3D(n1, n2, n3); p.resize3D(n1, n2, n3); z.resize3D(n1, n2, n3); s.resize3D(n1, n2, n3);
        phi_solid = phi_solid_init;
        phi.resize3D(n1, n2, n3);
        n_liquid = 0;
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
            phi(i, j, k) = interpolate_value(i+0.5, j+0.5, k+0.5, phi_init);
            if (phi(i, j, k) < 0) {
                n_liquid++;
            }
        }       
        particles_num.resize3D(n1, n2, n3);
        
    // make the particles large enough so they always appear on the grid
    int resolution = 1;
    particle_radius = (l * 1.01 * sqrt(3.0) / 2.0) / resolution * 0.9;
    // init particles
    int seed = 0;
    Vector3d pos;
	for (int n = 0; n < 64; n++) {
		for (int i = resolution; i < resolution * (n1 - 1); i++)
        for (int j = resolution; j < resolution * (n2 - 1); j++)
        for (int k = resolution; k < resolution * (n3 - 1); k++) {
			double a = randhashd(seed++, -.5 / resolution, .5 / resolution);
            double b = randhashd(seed++, -.5 / resolution, .5 / resolution);
            double c = randhashd(seed++, -.5 / resolution, .5 / resolution);
            pos = Vector3d((i + 0.5) / resolution + a, (j + 0.5) / resolution + b, (k + 0.5) / resolution + c);			
			if (interpolate_value(pos, phi_init) <= -particle_radius) {
                Assert(pos(0) >= 0 && pos(1) >= 0 && pos(2) >= 0, "FluidSim::FluidSim", "x_, y_, z_ should >= 0.");
                Assert(pos(0) < phi_solid.getN1()-1 && pos(1) < phi_solid.getN2()-1 && pos(2) < phi_solid.getN3()-1, "FluidSim::FluidSim", "out of index.");
				if (getPhiSolid(pos * l) > 0)
					particles.push_back(pos * l);
			}
		}
	}
    // init particles_output
    int resolution_output = 1;
    particle_output_radius = (l * 1.01 * sqrt(3.0) / 2.0) / resolution_output;
    seed = 0;    
	for (int n = 0; n < 1; n++) {
		for (int i = resolution_output; i < resolution_output * (n1 - 1); i++)
        for (int j = resolution_output; j < resolution_output * (n2 - 1); j++)
        for (int k = resolution_output; k < resolution_output * (n3 - 1); k++) {
			double a = randhashd(seed++, -.5 / resolution_output, .5 / resolution_output);
            double b = randhashd(seed++, -.5 / resolution_output, .5 / resolution_output);
            double c = randhashd(seed++, -.5 / resolution_output, .5 / resolution_output);
            pos = Vector3d((i + 0.5) / resolution_output + a, (j + 0.5) / resolution_output + b, (k + 0.5) / resolution_output + c);			
			if (interpolate_value(pos, phi_init) <= -particle_output_radius) {
                Assert(pos(0) >= 0 && pos(1) >= 0 && pos(2) >= 0, "FluidSim::FluidSim", "x_, y_, z_ should >= 0.");
                Assert(pos(0) < phi_solid.getN1()-1 && pos(1) < phi_solid.getN2()-1 && pos(2) < phi_solid.getN3()-1, "FluidSim::FluidSim", "out of index.");
				if (getPhiSolid(pos * l) > 0)
					particles_output.push_back(pos * l);
			}
		}
	}
    computePhi();

    std::cout << "FluidSim is initialized." << std::endl;
}

void FluidSim::advance(const double& time_step){
    // debug
    Assert(time_step > 0, "FluidSim::advance", "Non-positive time step.");
    computeN();
    
    std::cout 
    << "n_liquid: " << n_liquid 
    << " V_liquid: " << max(u.cwiseAbs().maxCoeff(), v.cwiseAbs().maxCoeff(), w.cwiseAbs().maxCoeff()) 
    << " p_max: " << p.cwiseAbs().maxCoeff() 
    << " d_max: " << d.cwiseAbs().maxCoeff() 
    << " rigidbody v_y: " << rigidbody.getVelocity()(1)
    << std::endl;
    
    double t = 0;
    bool done = false;
    while (!done) {
        // 5-dx rule
        double dt = 5 * l / (u.cwiseAbs().maxCoeff() + v.cwiseAbs().maxCoeff() + w.cwiseAbs().maxCoeff());
        if (t + dt >= time_step) {
            dt = time_step - t;
            done = true;
        }

        // core steps
        if (add_rigidbody) {
            rigidbody.advance(dt, phi_solid, l);
            updateRigidBodyGrids();
            recomputeRigidBodyMass();
        }
        advectParticles(dt);        
        computePhi();
        computeWeights();
        advect(dt);
        applyForce(dt);
        //extrapolate();
        if (add_rigidbody) {
            recomputeSolidVelocity();
        }
        constrain();
        project();
        extrapolate();        
        constrain();
        t += dt;
    } 
    return;
}

double FluidSim::getPhi(const Vector3d& pos) {  
    Vector3d position = pos;
    position(0) = clamp(position(0), l, (n1-1) * l);
    position(1) = clamp(position(1), l, (n2-1) * l); 
    position(2) = clamp(position(2), l, (n3-1) * l); 
    return interpolate_value(position / l - Vector3d(0.5, 0.5, 0.5), phi);
}

double FluidSim::getPhiSolid(const Vector3d& pos) {    
    Assert(pos(0) >= 0 && pos(1) >= 0 && pos(2) >= 0, "FluidSim::getPhiSolid", "pos(0), pos(1) and pos(2) should >= 0"); 
    Assert(pos(0) / l < phi_solid.getN1()-1 && pos(1) / l < phi_solid.getN2()-1 && pos(2) / l < phi_solid.getN3()-1, "FluidSim::getPhiSolid", "out of index");  
    return interpolate_value(pos / l, phi_solid);
}

double FluidSim::getPhiRigidBody(const Vector3d& pos) {
    Assert(add_rigidbody, "FluidSim::getPhiRigidBody", "rigidbody not defined.");
    Assert(pos(0) >= 0 && pos(1) >= 0 && pos(2) >= 0, "FluidSim::getPhiRigidBody", "pos(0), pos(1) and pos(2) should >= 0"); 
    Assert(pos(0) / l < phi_solid.getN1()-1 && pos(1) / l < phi_solid.getN2()-1 && pos(2) / l < phi_solid.getN3()-1, "FluidSim::getPhiRigidBody", "out of index");  
    return interpolate_value(pos / l, phi_rigidbody);
}

double FluidSim::getPhiSolidRigidBody(const Vector3d& pos) {
    Assert(add_rigidbody, "FluidSim::getPhiSolidRigidBody", "rigidbody not defined.");
    Assert(pos(0) >= 0 && pos(1) >= 0 && pos(2) >= 0, "FluidSim::getPhiSolidRigidBody", "pos(0), pos(1) and pos(2) should >= 0"); 
    Assert(pos(0) / l < phi_solid.getN1()-1 && pos(1) / l < phi_solid.getN2()-1 && pos(2) / l < phi_solid.getN3()-1, "FluidSim::getPhiSolidRigidBody", "out of index");    
    return interpolate_value(pos / l, phi_rigidbody);
}

Vector3d FluidSim::getVelocity(const Vector3d& position){
    Vector3d lattice_position = position / l;    
    // project_to_interior
    lattice_position(0) = clamp(lattice_position(0), 1., n1-1.);
    lattice_position(1) = clamp(lattice_position(1), 1., n2-1.);
    lattice_position(2) = clamp(lattice_position(2), 1., n3-1.);

    double u_ave = interpolate_value(lattice_position(0), lattice_position(1) - 0.5, lattice_position(2) - 0.5, u);
    double v_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1), lattice_position(2) - 0.5, v);
    double w_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1) - 0.5, lattice_position(2), w);

    return Vector3d(u_ave, v_ave, w_ave);
}

Vector3d FluidSim::getSolidVelocity(const Vector3d& position){
    Vector3d lattice_position = position / l;    
    // project_to_interior
    lattice_position(0) = clamp(lattice_position(0), 1., n1-1.);
    lattice_position(1) = clamp(lattice_position(1), 1., n2-1.);
    lattice_position(2) = clamp(lattice_position(2), 1., n3-1.);

    double u_ave = interpolate_value(lattice_position(0), lattice_position(1) - 0.5, lattice_position(2) - 0.5, u_solid);
    double v_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1), lattice_position(2) - 0.5, v_solid);
    double w_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1) - 0.5, lattice_position(2), w_solid);
    
    return Vector3d(u_ave, v_ave, w_ave);
}

Vector3d FluidSim::traceRk2(const Vector3d& position, const double& dt){
    Vector3d input = position;
	Vector3d velocity = getVelocity(input);
	velocity = getVelocity(input + 0.5f * dt * velocity);
	input += dt * velocity;
	return input;
}

Field3d FluidSim::applyA(const Field3d& x){
    Field3d ans(x.getN1(), x.getN2(), x.getN3());
    ans.setZero();
    //#pragma omp parallel for
    for (int i = 0; i < x.getN1(); i++) for (int j = 0; j < x.getN2(); j++) for (int k = 0; k < x.getN3(); k++) {
		ans(i, j, k) = Adiag(i, j, k) * x(i, j, k);
        if(i < x.getN1()-1)
		    ans(i, j, k) += Aplusi(i, j, k) * x(i + 1, j, k);
        if(j < x.getN2()-1)    
		    ans(i, j, k) += Aplusj(i, j, k) * x(i, j + 1, k);
        if(k < x.getN3()-1)
		    ans(i, j, k) += Aplusk(i, j, k) * x(i, j, k + 1);
        if(i > 0)
		    ans(i, j, k) += Aplusi(i - 1, j, k) * x(i - 1, j, k);
        if(j > 0)
		    ans(i, j, k) += Aplusj(i, j - 1, k) * x(i, j - 1, k);
        if(k > 0)
		    ans(i, j, k) += Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
    if (add_rigidbody) {
        ans += fluid_density * (J_x.reshaped()).dot(x.reshaped()) * J_x / rigidbody.getMass()(0, 0);
        ans += fluid_density * (J_y.reshaped()).dot(x.reshaped()) * J_y / rigidbody.getMass()(1, 1);
        ans += fluid_density * (J_z.reshaped()).dot(x.reshaped()) * J_z / rigidbody.getMass()(2, 2);
        Vector3d J_rot_dot_p((J_rot_x.reshaped()).dot(x.reshaped()), (J_rot_y.reshaped()).dot(x.reshaped()), (J_rot_z.reshaped()).dot(x.reshaped()));
        J_rot_dot_p = rigidbody.getI_inv() * J_rot_dot_p;
        ans += fluid_density * J_rot_x * J_rot_dot_p(0);
        ans += fluid_density * J_rot_y * J_rot_dot_p(1);
        ans += fluid_density * J_rot_z * J_rot_dot_p(2);
    }
    return ans;
}

void FluidSim::applyA(const Field3d& x, Field3d& ans){
    //#pragma omp parallel for
    for (int i = 0; i < x.getN1(); i++) for (int j = 0; j < x.getN2(); j++) for (int k = 0; k < x.getN3(); k++) {
		ans(i, j, k) = Adiag(i, j, k) * x(i, j, k);
        if(i < x.getN1()-1)
		    ans(i, j, k) += Aplusi(i, j, k) * x(i + 1, j, k);
        if(j < x.getN2()-1)    
		    ans(i, j, k) += Aplusj(i, j, k) * x(i, j + 1, k);
        if(k < x.getN3()-1)
		    ans(i, j, k) += Aplusk(i, j, k) * x(i, j, k + 1);
        if(i > 0)
		    ans(i, j, k) += Aplusi(i - 1, j, k) * x(i - 1, j, k);
        if(j > 0)
		    ans(i, j, k) += Aplusj(i, j - 1, k) * x(i, j - 1, k);
        if(k > 0)
		    ans(i, j, k) += Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
    if (add_rigidbody) {
        ans += fluid_density * (J_x.reshaped()).dot(x.reshaped()) * J_x / rigidbody.getMass()(0, 0);
        ans += fluid_density * (J_y.reshaped()).dot(x.reshaped()) * J_y / rigidbody.getMass()(1, 1);
        ans += fluid_density * (J_z.reshaped()).dot(x.reshaped()) * J_z / rigidbody.getMass()(2, 2);
        Vector3d J_rot_dot_p((J_rot_x.reshaped()).dot(x.reshaped()), (J_rot_y.reshaped()).dot(x.reshaped()), (J_rot_z.reshaped()).dot(x.reshaped()));
        J_rot_dot_p = rigidbody.getI_inv() * J_rot_dot_p;
        ans += fluid_density * J_rot_x * J_rot_dot_p(0);
        ans += fluid_density * J_rot_y * J_rot_dot_p(1);
        ans += fluid_density * J_rot_z * J_rot_dot_p(2);
    }
    return;
}

void FluidSim::extrapolate(){
    valid_u.setZero();
    valid_v.setZero();
    valid_w.setZero();
    // reset valid_u
    //#pragma omp parallel for
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_u(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i - 1, j, k) < 0)) {
			valid_u(i, j, k) = 1;
		}		
    }
    // reset valid_v
    //#pragma omp parallel for
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_v(i, j, k) > 0 && (phi(i, j, k) < 0 ||  phi(i, j - 1, k) < 0)) {
			valid_v(i, j, k) = 1;
		}		
    }
    // reset valid_w
    //#pragma omp parallel for
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 1; k < n3; k++){
        if (weights_w(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i, j, k - 1) < 0)) {
			valid_w(i, j, k) = 1;
		}		
    }
    for (int idx = 0; idx < u.size(); idx++){
        if(valid_u[idx] == 0){
            u[idx] = 0;
        }
    }
    for (int idx = 0; idx < v.size(); idx++){
        if(valid_v[idx] == 0){
            v[idx] = 0;
        }
    }
    for (int idx = 0; idx < w.size(); idx++){
        if(valid_w[idx] == 0){
            w[idx] = 0;
        }
    }

	// Apply several iterations of a very simple propagation of valid velocity data in all directions
	extrapolate(u, u_new, valid_u, valid_u_old);
	extrapolate(v, v_new, valid_v, valid_v_old);
	extrapolate(w, w_new, valid_w, valid_w_old);
    return;
}

void FluidSim::extrapolate(Field3d& field, Field3d field_new, Field3i& valid, Field3i& valid_old){
    field_new = field;
	for (int num = 0; num < 10; num++) {
		double sum = 0;
		int count = 0;
		valid_old = valid;
        //#pragma omp parallel for schedule(dynamic, 3)
		for (int i = 1; i < field.getN1() - 1; i++) for (int j = 1; j < field.getN2() - 1; j++) for (int k = 1; k < field.getN3() - 1; k++) {
			sum = 0;
			count = 0;
			if (!valid_old(i, j, k)) {

				if (valid_old(i + 1, j, k)) {
					sum += field(i + 1, j, k);
					count++;
				}
				if (valid_old(i - 1, j, k)) {
					sum += field(i - 1, j, k);
					count++;
				}
				if (valid_old(i, j + 1, k)) {
					sum += field(i, j + 1, k);
					count++;
				}
				if (valid_old(i, j - 1, k)) {
					sum += field(i, j - 1, k);
					count++;
				}
				if (valid_old(i, j, k + 1)) {
					sum += field(i, j, k + 1);
					count++;
				}
				if (valid_old(i, j, k - 1)) {
					sum += field(i, j, k - 1);
					count++;
				}

				if (count > 0) {
					field_new(i, j, k) = sum / count;
					valid(i, j, k) = 1;
				}
			}
		}
		field = field_new;
	}
    return;
}

void FluidSim::advectParticles(const double& dt){
    // advect particles
    //#pragma omp parallel for schedule(dynamic, 10)
    for(int i=0; i<particles.size(); i++){
        particles[i] = traceRk2(particles[i], dt);
        particles[i](0) = clamp(particles[i](0), l, (n1-1) * l);
        particles[i](1) = clamp(particles[i](1), l, (n2-1) * l);
        particles[i](2) = clamp(particles[i](2), l, (n3-1) * l);
        // solid boundary
        double phi_val = getPhiSolid(particles[i]);
        Vector3d grad;        
        if(phi_val < 0){
            interpolate_gradient(grad, particles[i] / l, phi_solid);            
            grad.normalize();
            particles[i] -= phi_val * grad;
        }        
        if (add_rigidbody) {
            particles[i](0) = clamp(particles[i](0), l, (n1-1) * l);
            particles[i](1) = clamp(particles[i](1), l, (n2-1) * l);
            particles[i](2) = clamp(particles[i](2), l, (n3-1) * l);
            phi_val = getPhiRigidBody(particles[i]);  
            if (phi_val < 0) {
                interpolate_gradient(grad, particles[i] / l, phi_rigidbody);
                grad.normalize();
                particles[i] -= phi_val * grad;                
            }
        }
    }
    /*
    // advect particles_output
    for(int i=0; i<particles_output.size(); i++){
        particles_output[i] = traceRk2(particles_output[i], dt);
        particles_output[i](0) = clamp(particles_output[i](0), l, (n1-1) * l);
        particles_output[i](1) = clamp(particles_output[i](1), l, (n2-1) * l);
        particles_output[i](2) = clamp(particles_output[i](2), l, (n3-1) * l);
        // solid boundary
        double phi_val = getPhiSolid(particles_output[i]);
        Vector3d grad;
        if(phi_val < 0){
            interpolate_gradient(grad, particles_output[i] / l, phi_solid);            
            grad.normalize();
            particles_output[i] -= phi_val * grad;
        }
        if (add_rigidbody) {
            particles_output[i](0) = clamp(particles_output[i](0), l, (n1-1) * l);
            particles_output[i](1) = clamp(particles_output[i](1), l, (n2-1) * l);
            particles_output[i](2) = clamp(particles_output[i](2), l, (n3-1) * l);
            phi_val = getPhiRigidBody(particles_output[i]);  
            if (phi_val < 0) {
                interpolate_gradient(grad, particles_output[i] / l, phi_rigidbody);
                grad.normalize();
                particles_output[i] -= phi_val * grad;                
            }
        }
    }    
    */
    return;
}

void FluidSim::setParticlesOutput() {
    int max_num = 20;
    particles_num.setZero();
    particles_output.clear();
    int i, j, k;
    for (int idx = 0; idx < particles.size(); idx++) {
        i = static_cast<int>(particles[idx](0) / l);
        j = static_cast<int>(particles[idx](1) / l);
        k = static_cast<int>(particles[idx](2) / l);
        if (particles_num(i, j, k) > max_num) {
            continue;
        }
        else {
            particles_output.push_back(particles[idx]);
            particles_num(i, j, k) += 1;
        }
    }
    return;
}

void FluidSim::computePhi(){
    phi.setConstant(3*l);
    //#pragma omp parallel for 
    for(int p=0; p<particles.size(); p++){
        Vector3i cell_ind = (particles[p]/l).cast<int>();
        for (int k = max(0, cell_ind[2] - 1); k <= min(cell_ind[2] + 1, n3 - 1); k++) 
        for (int j = max(0, cell_ind[1] - 1); j <= min(cell_ind[1] + 1, n2 - 1); j++) 
        for (int i = max(0, cell_ind[0] - 1); i <= min(cell_ind[0] + 1, n1 - 1); i++) {
            Vector3d sample_pos(i*l, j*l, k*l);                    
            double test_val = (sample_pos - particles[p]).norm() - particle_radius;
            if (test_val < phi(i, j, k)) {
                phi(i, j, k) = test_val;
            }
        }            		
    }

    // extend phi into solid
    //#pragma omp parallel for
    /*
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
		if (phi(i, j, k) < 0.5 * l && getPhiSolid(Vector3d(i+0.5, j+0.5, k+0.5) * l) < 0) {
			phi(i, j, k) = -0.5 * l;
		}
	}
    */
    return;
}

void FluidSim::advect(const double& dt){
    semiLagrangian(u, u_new, dt, 0);
    semiLagrangian(v, v_new, dt, 1);
    semiLagrangian(w, w_new, dt, 2);
    u = u_new;
    v = v_new;
    w = w_new;
    return;
}

void FluidSim::semiLagrangian(const Field3d& field, Field3d& field_new, const double& dt, const int& id){
    assert(id == 0 || id == 1 || id == 2);
    Vector3d offset;
    if(id == 0){
        offset = Vector3d(0, 0.5, 0.5);
    }
    else if(id == 1){
        offset = Vector3d(0.5, 0, 0.5);
    }
    else if(id == 2){
        offset = Vector3d(0.5, 0.5, 0);
    }
    //#pragma omp parallel for schedule(dynamic, 5)
    for (int i = 1; i < field.getN1()-1; i++) for (int j = 1; j < field.getN2()-1; j++) for (int k = 1; k < field.getN3()-1; k++) {
		Vector3d pos = Vector3d(i * l, j * l, k * l) + offset * l;        
        field_new(i, j, k) = getVelocity(traceRk2(pos, -dt))(id);        
	}
    return;
}

void FluidSim::applyForce(const double& dt){
    for (int idx = 0; idx < v.size(); idx++) {
        v[idx] += -G * dt;
    }

    return;
}

void FluidSim::computeWeights(){
    // compute weights_u
    //#pragma omp parallel for schedule(dynamic, 5)
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++) {
        weights_u(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i, j+1, k), phi_solid(i, j+1, k+1), phi_solid(i, j, k+1));
        if (add_rigidbody) {
            weights_u(i, j, k) -= weights_rigid_u(i, j, k);
        }    
        weights_u(i, j, k) = clamp(weights_u(i, j, k), 0., 1.);
    }
    // compute weights_v
    //#pragma omp parallel for schedule(dynamic, 5)
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3; k++){
        weights_v(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i+1, j, k), phi_solid(i+1, j, k+1), phi_solid(i, j, k+1));  
        if (add_rigidbody) {
            weights_v(i, j, k) -= weights_rigid_v(i, j, k);
        }        
        weights_v(i, j, k) = clamp(weights_v(i, j, k), 0., 1.);
    }
    // compute weights_w
    //#pragma omp parallel for schedule(dynamic, 5)
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3+1; k++){
        weights_w(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i, j+1, k), phi_solid(i+1, j+1, k), phi_solid(i+1, j, k));     
        if (add_rigidbody) {
            weights_w(i, j, k) -= weights_rigid_w(i, j, k);
        }     
        weights_w(i, j, k) = clamp(weights_w(i, j, k), 0., 1.);
    }
    return;
}

void FluidSim::project() {
    // compute J
    if (add_rigidbody) {
        J_x.setZero(); J_y.setZero(); J_rot_x.setZero(); J_rot_y.setZero(); J_rot_z.setZero();     
        //#pragma omp parallel for schedule(dynamic, 2)   
        for(int i = 0; i < n1; i++) for(int j = 0; j < n2; j++) for(int k = 0; k < n3; k++) {
            if (phi(i, j, k) < 0 && getPhiRigidBody(Vector3d((i+0.5) * l, (j+0.5) * l, (k+0.5) * l)) >= 0) {
                J_x(i, j, k) = weights_rigid_u(i+1, j, k) - weights_rigid_u(i, j, k);
                J_y(i, j, k) = weights_rigid_v(i, j+1, k) - weights_rigid_v(i, j, k);
                J_z(i, j, k) = weights_rigid_w(i, j, k+1) - weights_rigid_w(i, j, k);
                Vector3d pos((i+0.5) * l, (j+0.5) * l, (k+0.5) * l);
                pos = pos - rigidbody.getC();
                J_rot_x(i, j, k) = pos(1) * J_z(i, j, k) - pos(2) * J_y(i, j, k);
                J_rot_y(i, j, k) = pos(2) * J_x(i, j, k) - pos(0) * J_z(i, j, k);
                J_rot_z(i, j, k) = pos(0) * J_y(i, j, k) - pos(1) * J_x(i, j, k);
            }            
        }         
    }

    // compute A and d
    Adiag.setZero(); Aplusi.setZero(); Aplusj.setZero(); Aplusk.setZero(); d.setZero();
    //#pragma omp parallel for 
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        double theta = 0;
        if(phi(i, j, k) < 0){
            // right
            if(weights_u(i+1, j, k) > 0){
                if(phi(i+1, j, k) < 0){
                    Adiag(i, j, k) += weights_u(i+1, j, k);
                    Aplusi(i, j, k) += -weights_u(i+1, j, k);
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i+1, j, k));
                    if(theta < 0.1) theta = 0.1;                    
                    Adiag(i, j, k) += weights_u(i+1, j, k) / theta;
                }
                d(i, j, k) += -weights_u(i+1, j, k) * u(i+1, j, k);
            }
            // left            
            if(weights_u(i, j, k) > 0){
                if(phi(i-1, j, k) < 0){
                    Adiag(i, j, k) += weights_u(i, j, k);                    
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i-1, j, k));
                    if(theta < 0.1) theta = 0.1;
                    Adiag(i, j, k) += weights_u(i, j, k) / theta;
                }
                d(i, j, k) += weights_u(i, j, k) * u(i, j, k);
            }
            // top
            if(weights_v(i, j+1, k) > 0){
                if(phi(i, j+1, k) < 0){
                    Adiag(i, j, k) += weights_v(i, j+1, k);
                    Aplusj(i, j, k) += -weights_v(i, j+1, k);
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i, j+1, k));
                    if(theta < 0.1) theta = 0.1;                    
                    Adiag(i, j, k) += weights_v(i, j+1, k) / theta;
                }
                d(i, j, k) += -weights_v(i, j+1, k) * v(i, j+1, k);
            }
            // bottom            
            if(weights_v(i, j, k) > 0){
                if(phi(i, j-1, k) < 0){
                    Adiag(i, j, k) += weights_v(i, j, k);                
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i, j-1, k));
                    if(theta < 0.1) theta = 0.1;                    
                    Adiag(i, j, k) += weights_v(i, j, k) / theta;
                }
                d(i, j, k) += weights_v(i, j, k) * v(i, j, k);
            }
            // far
            if(weights_w(i, j, k+1) > 0){
                if(phi(i, j, k+1) < 0){
                    Adiag(i, j, k) += weights_w(i, j, k+1);
                    Aplusk(i, j, k) += -weights_w(i, j, k+1);
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i, j, k+1));
                    if(theta < 0.1) theta = 0.1;                    
                    Adiag(i, j, k) += weights_w(i, j, k+1) / theta;
                }
                d(i, j, k) += -weights_w(i, j, k+1) * w(i, j, k+1);
            }
            // near
            if(weights_w(i, j, k) > 0){
                if(phi(i, j, k-1) < 0){
                    Adiag(i, j, k) += weights_w(i, j, k);                    
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i, j, k-1));
                    if(theta < 0.1) theta = 0.1;                    
                    Adiag(i, j, k) += weights_w(i, j, k) / theta;
                }
                d(i, j, k) += weights_w(i, j, k) * w(i, j, k);
            }            
        }        
    }

    // add rigidbody terms
    if (add_rigidbody) {
        // rhs
        d -= J_x * rigidbody.getVelocity()(0);
        d -= J_y * rigidbody.getVelocity()(1);
        d -= J_z * rigidbody.getVelocity()(2);
        d -= J_rot_x * rigidbody.getOmega()(0);
        d -= J_rot_y * rigidbody.getOmega()(1);
        d -= J_rot_z * rigidbody.getOmega()(2);
    }

    // ensure A is positive definite
    //#pragma omp parallel for 
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if(Adiag(i, j, k)==0 && Aplusi(i, j, k)==0 && Aplusj(i, j, k)==0 && Aplusk(i, j, k)==0 && 
        (i==0 || Aplusi(i-1, j, k)==0) && (j==0 || Aplusj(i, j-1, k)==0) && (k==0 || Aplusk(i, j, k-1)==0)){
            Adiag(i, j, k) = 1;
            d(i, j, k) = 0;
        }		
	}
    // TODO: remove the 1D null space
    // solve A * p = d
    solve(4000);

    // update u
    //#pragma omp parallel for schedule(dynamic, 3)   
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        double theta = 0;
        if (weights_u(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i - 1, j, k) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i - 1, j, k) > 0) theta = computeFraction(phi(i, j, k), phi(i - 1, j, k));
            if (theta < 0.1) theta = 0.1;            
            u(i, j, k) += -((p(i, j, k) - p(i - 1, j, k))) / theta;
		}		
    }
    // update v
    //#pragma omp parallel for schedule(dynamic, 3)   
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) for (int k = 0; k < n3; k++){
        double theta = 0;
        if (weights_v(i, j, k) > 0 && (phi(i, j, k) < 0 ||  phi(i, j - 1, k) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i, j - 1, k) > 0) theta = computeFraction(phi(i, j, k), phi(i, j - 1, k));
            if (theta < 0.1) theta = 0.1;            
            v(i, j, k) += -((p(i, j, k) - p(i, j - 1, k))) / theta;
		}		
    }
    // update w
    //#pragma omp parallel for schedule(dynamic, 3)   
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 1; k < n3; k++){
        double theta = 0;
        if (weights_w(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i, j, k - 1) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i, j, k - 1) > 0) theta = computeFraction(phi(i, j, k), phi(i, j, k - 1));
            if (theta < 0.1) theta = 0.1;            
            w(i, j, k) += -((p(i, j, k) - p(i, j, k - 1))) / theta;
		}		
    }

    // update rigidbody
    if (add_rigidbody) {        
        Vector3d velocity_temp = rigidbody.getVelocity();        
        Vector3d L_temp = rigidbody.getI() * rigidbody.getOmega();
        double velocity_temp_0, velocity_temp_1, velocity_temp_2;
        double L_temp_0, L_temp_1, L_temp_2;
        velocity_temp_0 = velocity_temp(0); velocity_temp_1 = velocity_temp(1); velocity_temp_2 = velocity_temp(2);
        L_temp_0 = L_temp(0); L_temp_1 = L_temp(1); L_temp_2 = L_temp(2);
        //#pragma omp parallel for reduction(+:velocity_temp_0) reduction(+:velocity_temp_1) reduction(+:velocity_temp_2) reduction(+:L_temp_0) reduction(+:L_temp_1) reduction(+:L_temp_2)
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
            if (phi(i, j, k) < 0 && getPhiRigidBody(Vector3d((i+0.5) * l, (j+0.5) * l, (k+0.5) * l)) >= 0) {
                velocity_temp_0 += fluid_density * J_x(i, j, k) * p(i, j, k) / rigidbody.getMass()(0, 0);
                velocity_temp_1 += fluid_density * J_y(i, j, k) * p(i, j, k) / rigidbody.getMass()(1, 1);
                velocity_temp_2 += fluid_density * J_z(i, j, k) * p(i, j, k) / rigidbody.getMass()(2, 2);
                L_temp_0 += fluid_density * J_rot_x(i, j, k) * p(i, j, k);
                L_temp_1 += fluid_density * J_rot_y(i, j, k) * p(i, j, k);
                L_temp_2 += fluid_density * J_rot_z(i, j, k) * p(i, j, k);
            }
        }
        rigidbody.setVelocity(Vector3d(velocity_temp_0, velocity_temp_1, velocity_temp_2));
        rigidbody.setOmega(rigidbody.getI_inv() * Vector3d(L_temp_0, L_temp_1, L_temp_2));        
    }

    return;
}

void FluidSim::solve(int maxIterations){
    // set p = 0
	//p.setConstant(0);
	if (d.cwiseAbs().maxCoeff() == 0) {		
		p.setConstant(0);
        return;
	}
    d = d - applyA(p);
	// TODO
	// applyPrecon();
    // ######
    z = d;
	s = z;
	double sigma = z.dot(d);
	double sigma_new;
	double alpha;
    while (d.cwiseAbs().maxCoeff() > 1e-3) {
        Assert(maxIterations >= 0, "FluidSim::Solve", "CG failed.");
        applyA(s, z);
		alpha = sigma / z.dot(s);
		p += s * alpha;
		d -= z * alpha;		
        // TODO
		//applyPrecon();
        // ######		
		z = d;
		sigma_new = z.dot(d);
		alpha = sigma_new / sigma;
		sigma = sigma_new;        		
        s = z + alpha * s;
        maxIterations--;
    }	
	return;
}

void FluidSim::constrain(){
    u_new = u;
    v_new = v;
    w_new = w;

    // constrain u
    //#pragma omp parallel for schedule(dynamic, 3)
    for(int i=1; i<n1; i++) for(int j=1; j<n2-1; j++) for(int k=1; k<n3-1; k++){
        Vector3d pos, vec, vec_solid, normal_solid;
        double vec_n, vec_solid_n;
        if(weights_u(i, j, k) == 0){
            pos = Vector3d(i, j+0.5, k+0.5) * l;
            vec = getVelocity(pos);
            vec_solid = getSolidVelocity(pos);
            if (add_rigidbody) {
                interpolate_gradient(normal_solid, pos / l, phi_solid_rigidbody);
            }
            else {
                interpolate_gradient(normal_solid, pos / l, phi_solid);
            }            
            normal_solid.normalize();
            vec_n = vec.dot(normal_solid);
            vec_solid_n = vec_solid.dot(normal_solid);            
            vec -= vec_n * normal_solid;
            vec += vec_solid_n * normal_solid;            
            u_new(i, j, k) = vec(0);
        }
    }
    // constrain v
    //#pragma omp parallel for schedule(dynamic, 3)
    for(int i=1; i<n1-1; i++) for(int j=1; j<n2; j++) for(int k=1; k<n3-1; k++){
        Vector3d pos, vec, vec_solid, normal_solid;
        double vec_n, vec_solid_n;
        if(weights_v(i, j, k) == 0){
            pos = Vector3d(i+0.5, j, k+0.5) * l;
            vec = getVelocity(pos);
            vec_solid = getSolidVelocity(pos);
            if (add_rigidbody) {
                interpolate_gradient(normal_solid, pos / l, phi_solid_rigidbody);
            }
            else {
                interpolate_gradient(normal_solid, pos / l, phi_solid);
            }
            normal_solid.normalize();
            vec_n = vec.dot(normal_solid);
            vec_solid_n = vec_solid.dot(normal_solid);
            vec -= vec_n * normal_solid;
            vec += vec_solid_n * normal_solid;
            v_new(i, j, k) = vec(1);
        }
    }
    // constrain w
    //#pragma omp parallel for schedule(dynamic, 3)
    for(int i=1; i<n1-1; i++) for(int j=1; j<n2-1; j++) for(int k=1; k<n3; k++){
        Vector3d pos, vec, vec_solid, normal_solid;
        double vec_n, vec_solid_n;
        if(weights_w(i, j, k) == 0){
            pos = Vector3d(i+0.5, j+0.5, k) * l;
            vec = getVelocity(pos);
            vec_solid = getSolidVelocity(pos);
            if (add_rigidbody) {
                interpolate_gradient(normal_solid, pos / l, phi_solid_rigidbody);
            }
            else {
                interpolate_gradient(normal_solid, pos / l, phi_solid);
            }
            normal_solid.normalize();
            vec_n = vec.dot(normal_solid);
            vec_solid_n = vec_solid.dot(normal_solid);
            vec -= vec_n * normal_solid;
            vec += vec_solid_n * normal_solid;
            w_new(i, j, k) = vec(2);
        }
    }
    // boundary
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++) if(i==0 || i==n1) 
        u_new(i, j, k) = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3; k++) if(j==0 || j==n2) 
        v_new(i, j, k) = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3+1; k++) if(k==0 || k==n3) 
        w_new(i, j, k) = 0;

    u = u_new;
    v = v_new;
    w = w_new;
    return;
}

void FluidSim::computeN(){
    n_liquid = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++) {
        if(phi(i, j, k) < 0) n_liquid++; 
    }
}

void FluidSim::setVelocity(const Vector3d& vec){
    u.setConstant(vec(0));
    v.setConstant(vec(1));
    w.setConstant(vec(2));
    return;
}

void FluidSim::addRigidBody(
        const std::vector<Vector3d>& vertices_init,
        const std::vector<Face>& faces_init,
        const double& m_init, 
        const Matrix3d& I_init, 
        const Vector3d& c_init, 
        const Quaterniond& q_init, 
        const Vector3d& velocity_init, 
        const Vector3d& omega_init,
        const FuncPtr& get_sdf_func_init
    ) {
        Assert(!add_rigidbody, "FluidSim::addRigidBody", "Only one rigid body is supported.");
        std::cout << "Rigidbody is added." << std::endl;
        add_rigidbody = true;

        // init
        rigidbody = RigidBody(vertices_init, faces_init, m_init, I_init, c_init, q_init, velocity_init, omega_init, get_sdf_func_init);
        J_x.resize3D(n1, n2, n3); J_y.resize3D(n1, n2, n3); J_z.resize3D(n1, n2, n3);
        J_rot_x.resize3D(n1, n2, n3); J_rot_y.resize3D(n1, n2, n3); J_rot_z.resize3D(n1, n2, n3);
        phi_rigidbody.resize3D(n1+1, n2+1, n3+1); phi_solid_rigidbody.resize3D(n1+1, n2+1, n3+1);
        weights_rigid_u.resize3D(n1+1, n2, n3); weights_rigid_v.resize3D(n1, n2+1, n3); weights_rigid_w.resize3D(n1, n2, n3+1);

        updateRigidBodyGrids();
        //compute rigidbody_density
        rigidbody_density = (rigidbody.getMass()(0, 0) + rigidbody.getMass()(1, 1) + rigidbody.getMass()(2, 2)) 
        / (weights_rigid_u.array().abs().sum() + weights_rigid_v.array().abs().sum() + weights_rigid_w.array().abs().sum());
        fluid_density = rigidbody_density * 3.8;

        return;
}

void FluidSim::updateRigidBodyGrids() {
    // update phi_rigidbody
    //#pragma omp parallel for 
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) for (int k = 0; k < n3+1; k++) {
        Vector3d pos(i * l, j * l, k * l);
        phi_rigidbody(i, j, k) = rigidbody.getSDF(pos);
    }
    // update phi_solid_rigidbody (no intersection)
    //#pragma omp parallel for 
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) for (int k = 0; k < n3+1; k++) {
        if (abs(phi_rigidbody(i, j, k)) > abs(phi_solid(i, j, k))) {
            phi_solid_rigidbody(i, j, k) = phi_solid(i, j, k);
        }
        else {
            phi_solid_rigidbody(i, j, k) = phi_rigidbody(i, j, k);
        }
    }    
    // compute weights_rigid_u
    //#pragma omp parallel for 
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for (int k = 0; k < n3; k++) {
        weights_rigid_u(i, j, k) = computeFraction(phi_rigidbody(i, j, k), phi_rigidbody(i, j+1, k), phi_rigidbody(i, j+1, k+1), phi_rigidbody(i, j, k+1));
        weights_rigid_u(i, j, k) = clamp(weights_rigid_u(i, j, k), 0., 1.);
    }
    // compute weights_rigid_v
    //#pragma omp parallel for 
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for (int k = 0; k < n3; k++) {
        weights_rigid_v(i, j, k) = computeFraction(phi_rigidbody(i, j, k), phi_rigidbody(i+1, j, k), phi_rigidbody(i+1, j, k+1), phi_rigidbody(i, j, k+1)); 
        weights_rigid_v(i, j, k) = clamp(weights_rigid_v(i, j, k), 0., 1.);
    }
    // compute weights_rigid_w
    //#pragma omp parallel for 
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for (int k = 0; k < n3+1; k++) {
        weights_rigid_w(i, j, k) = computeFraction(phi_rigidbody(i, j, k), phi_rigidbody(i, j+1, k), phi_rigidbody(i+1, j+1, k), phi_rigidbody(i+1, j, k));
        weights_rigid_w(i, j, k) = clamp(weights_rigid_w(i, j, k), 0., 1.);
    }
    return;
}

void FluidSim::recomputeRigidBodyMass() {
    // recompute effective mass
    rigidbody.setMass(
        rigidbody_density * weights_rigid_u.array().abs().sum(),
        rigidbody_density * weights_rigid_v.array().abs().sum(),
        rigidbody_density * weights_rigid_w.array().abs().sum()
    );
}

void FluidSim::recomputeSolidVelocity() {
    // compute u_solid
    u_solid.setZero();
    //#pragma omp parallel for schedule(dynamic, 2)
    for(int i=0; i < n1+1; i++) for(int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
        Vector3d pos(i * l, (j+0.5) * l, (k+0.5) * l);        
        pos(0) = clamp(pos(0), 0., (n1-1e-6) * l);
        pos(1) = clamp(pos(1), 0., (n2-1-1e-6) * l);
        pos(2) = clamp(pos(2), 0., (n3-1-1e-6) * l);
        if (getPhiSolid(pos) < 0) {
            u_solid(i, j, k) = 0.;
        }
        else if (getPhiRigidBody(pos) < 0) {
            u_solid(i, j, k) = rigidbody.getVelocity(pos)(0);
        }        
    }
    
    // compute v_solid
    v_solid.setZero();
    //#pragma omp parallel for schedule(dynamic, 2)
    for(int i=0; i < n1; i++) for(int j = 0; j < n2+1; j++) for (int k = 0; k < n3; k++) {
        Vector3d pos((i+0.5) * l, j * l, (k+0.5) * l);
        pos(0) = clamp(pos(0), 0., (n1-1-1e-6) * l);
        pos(1) = clamp(pos(1), 0., (n2-1e-6) * l); 
        pos(2) = clamp(pos(2), 0., (n3-1-1e-6) * l);
        if (getPhiSolid(pos) < 0) {
            v_solid(i, j, k) = 0.;            
        }
        else if (getPhiRigidBody(pos) < 0) {
            v_solid(i, j, k) = rigidbody.getVelocity(pos)(1);
        }
    }

    // compute w_solid
    w_solid.setZero();
    //#pragma omp parallel for schedule(dynamic, 2)
    for(int i=0; i < n1; i++) for(int j = 0; j < n2; j++) for (int k = 0; k < n3+1; k++) {
        Vector3d pos((i+0.5) * l, (j+0.5) * l, k * l);
        pos(0) = clamp(pos(0), 0., (n1-1-1e-6) * l);
        pos(1) = clamp(pos(1), 0., (n2-1-1e-6) * l); 
        pos(2) = clamp(pos(2), 0., (n3-1e-6) * l);
        if (getPhiSolid(pos) < 0) {
            w_solid(i, j, k) = 0.;            
        }
        else if (getPhiRigidBody(pos) < 0) {
            w_solid(i, j, k) = rigidbody.getVelocity(pos)(2);
        }
    }
    return;    
}

void FluidSim::run(double time_step){
    long long count = 0;
    std::chrono::duration<double, std::milli> elapsed;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    while (1) {
        count++;
        if (count == 1) {
            count = 0;
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "fps: " << 1. / (elapsed.count() * 1e-3) << std::endl;
            start = end;
        }
		advance(time_step);
	}	
    return;
}

void FluidSim::outputPLY(int fps, int t, double dt, int n) {
	// delete .ply files
	std::string folderPath_ = projectPath + "/python/ply";
    if (!std::filesystem::exists(folderPath_)) {
        std::filesystem::create_directories(folderPath_);
    }
	try {
		for (const auto& entry : std::filesystem::directory_iterator(folderPath_)) {
			std::filesystem::remove(entry.path());
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}

    // sdf
    std::string sdfPath_ = projectPath + "/python/sdf";
    if (!std::filesystem::exists(sdfPath_)) {
        std::filesystem::create_directories(sdfPath_);
    }
	
	for (int frame_num = 0; frame_num < fps * t; frame_num++) {		
		std::cout << "frame: " << frame_num << std::endl;
        // output rigidBody .ply
		if (add_rigidbody) {
			std::string rigidbody_path = folderPath_ + "/rigidbody_" + std::to_string(frame_num) + ".ply";
			writePLY(rigidbody_path, rigidbody.getVertices(), rigidbody.getFaces());
		}        
		// output sdf .npy
        std::string sdf_path = sdfPath_ + "/sdf_" + std::to_string(frame_num) + ".npy";    
        writeNPY(sdf_path, phi);
        // sdf to ply
        std::string command = pythonPath + " " + projectPath + "/python/sdf2ply.py"
        + " -frame " + std::to_string(frame_num)
		+ " -l " + std::to_string(l)
		+ " -project_path " + projectPath;
	    system(command.c_str());

		// update system
		for (int i = 0; i < n; i++) {
			advance(dt);
		}
	}	

    try {
		for (const auto& entry : std::filesystem::directory_iterator(sdfPath_)) {
			std::filesystem::remove(entry.path());
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}

	return;
}

void FluidSim::outputVedio(const std::string& file_name, int fps, int t) {
    Assert(add_rigidbody, "FluidSim::outputVedio", "rigidbody does't exist.");
    // delete .png files
	std::string png_folder = projectPath + "/python/png";
    if (!std::filesystem::exists(png_folder)) {
        std::filesystem::create_directories(png_folder);
    }
	try {
		for (const auto& entry : std::filesystem::directory_iterator(png_folder)) {
			std::filesystem::remove(entry.path());
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}
   
    // ply to png
    std::string command = pythonPath + " " + projectPath + "/python/pbrt_renderer.py"
        + " -frame " + std::to_string(fps * t)
        + " -pbrt_folder " + pbrtPath
        + " -project_path " + projectPath;
    system(command.c_str());
    
    // png to mp4
    command = pythonPath + " " + projectPath + "/python/png2mp4.py"
        + " -fps " + std::to_string(fps)
        + " -output_file " + file_name + std::string(".mp4")
        + " -project_path " + projectPath;
    system(command.c_str());

    // delete .pbrt
    try {
		for (const auto& entry : std::filesystem::directory_iterator(projectPath + "/python/scene")) {
			if (entry.is_regular_file() && entry.path().extension() == ".pbrt") {
                std::filesystem::remove(entry.path());        
            }
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}

    return;
}

}