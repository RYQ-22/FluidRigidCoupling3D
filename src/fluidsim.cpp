#include "fluidsim.h"

FluidSim::FluidSim(
    const int& n1_init,
    const int& n2_init,
    const int& n3_init,
    const double& l_init,
    const std::vector<double>& phi_init,
    const std::vector<double>& phi_solid_init){
        assert(phi_init.size() == n1_init*n2_init*n3_init);
        assert(phi_solid_init.size() == (n1_init+1)*(n2_init+1)*(n3_init+1));        
        n1 = n1_init;
        n2 = n2_init;
        n3 = n3_init;
        l = l_init;       
        u.resize3D(n1+1, n2, n3); u_new.resize3D(n1+1, n2, n3); valid_u.resize3D(n1+1, n2, n3); valid_u_old.resize3D(n1+1, n2, n3); weights_u.resize3D(n1+1, n2, n3);
        v.resize3D(n1, n2+1, n3); v_new.resize3D(n1, n2+1, n3); valid_v.resize3D(n1, n2+1, n3); valid_v_old.resize3D(n1, n2+1, n3); weights_v.resize3D(n1, n2+1, n3);
        w.resize3D(n1, n2, n3+1); w_new.resize3D(n1, n2, n3+1); valid_w.resize3D(n1, n2, n3+1); valid_w_old.resize3D(n1, n2, n3+1); weights_w.resize3D(n1, n2, n3+1);
        Adiag.resize3D(n1, n2, n3); Aplusi.resize3D(n1, n2, n3); Aplusj.resize3D(n1, n2, n3); Aplusk.resize3D(n1, n2, n3);
        d.resize3D(n1, n2, n3); p.resize3D(n1, n2, n3); z.resize3D(n1, n2, n3); s.resize3D(n1, n2, n3);
        phi.resize3D(n1, n2, n3);
        phi_solid.resize3D(n1+1, n2+1, n3+1);
        for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){        
            phi(i, j, k) = phi_init[i*n2*n3+j*n3+k];
        }
        for(int i=0; i<n1+1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3+1; k++){
            phi_solid(i, j, k) = phi_solid_init[i*(n2+1)*(n3+1)+j*(n3+1)+k];
        }

    // make the particles large enough so they always appear on the grid
    particle_radius = (l * 1.01 * sqrt(3.0) / 2.0);
    // init particles
    int seed = 0;
	for (int n = 0; n < 64; n++) {
		for (int i = 1; i < n1-1; i++) for (int j = 1; j < n2-1; j++) for (int k = 1; k < n3-1; k++) {
			double a = randhashd(seed++, -1, 1); double b = randhashd(seed++, -1, 1); double c = randhashd(seed++, -1, 1);
			double x_ = i + a, y_ = j + b, z_ = k + c;
			if (interpolate_value(x_, y_, z_, phi) <= -particle_radius) {				
				if (phi_solid_ave(x_, y_, z_) > 0)
					particles.push_back(Vector3d(x_ * l, y_ * l, z_ * l));
			}
		}
	}
}

void FluidSim::advance(const double& dt){
    // debug
    std::cout 
    << "V_liquid: " << max(u.cwiseAbs().maxCoeff(), v.cwiseAbs().maxCoeff(), w.cwiseAbs().maxCoeff()) 
    << " p_max: " << p.cwiseAbs().maxCoeff() 
    << " d_max: " << d.cwiseAbs().maxCoeff() 
    << std::endl;

    advectParticles(dt);
    computePhi();
    computeWeights();
    advect(dt);
    applyForce(dt);
    //project();
    extrapolate();
    constrain();
    return;
}

double FluidSim::phi_solid_ave(const double& i, const double& j, const double& k){
    return interpolate_value(i+0.5, j+0.5, k+0.5, phi_solid);
}

double FluidSim::phi_solid_ave(const Vector3d& pos){
    return interpolate_value(pos.x()+0.5, pos.y()+0.5, pos.z()+0.5, phi_solid);
}

Vector3d FluidSim::getVelocity(const Vector3d& position){
    double u_ave = interpolate_value(position.x()/l+0.5, position.y()/l, position.z()/l, u);
    double v_ave = interpolate_value(position.x()/l, position.y()/l+0.5, position.z()/l, v);
    double w_ave = interpolate_value(position.x()/l, position.y()/l, position.z()/l+0.5, w);
    return Vector3d(u_ave, v_ave, w_ave);
}

Vector3d FluidSim::traceRk2(const Vector3d& position, const double& dt){
    Vector3d input = position;
	Vector3d velocity = getVelocity(input);
	velocity = getVelocity(input + 0.5f * dt * velocity);
	input += dt * velocity;
	return input;
}

double FluidSim::computeFraction(const double& phi0, const double& phi1, const double& phi2, const double& phi3){
    return (computeTriangleArea(phi0, phi1, phi2) + computeTriangleArea(phi2, phi3, phi0)) / 2;
}

double FluidSim::computeFraction(const double& phi0, const double& phi1){
    assert((phi0 < 0 && phi1 > 0) || (phi0 > 0 && phi1 < 0));
    if(phi0 < 0) return -phi0/(phi1-phi0);
    if(phi1 < 0) return -phi1/(phi0-phi1);
}

void FluidSim::applyA(const Field3d& x, Field3d& ans){
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
    return;
}


void FluidSim::extrapolate(){
    valid_u.setConstant(0);
    valid_v.setConstant(0);
    valid_w.setConstant(0);
    // reset valid_u
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_u(i, j, k) > 0 && (phi(i, j, k) < 0.0 || phi(i - 1, j, k) < 0.0)) {
			valid_u(i, j, k) = 1;
		}		
    }
    // reset valid_v
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_v(i, j, k) > 0 && (phi(i, j, k) < 0 ||  phi(i, j - 1, k) < 0)) {
			valid_v(i, j, k) = 1;
		}		
    }
    // reset valid_w
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 1; k < n3; k++){
        if (weights_w(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i, j, k - 1) < 0)) {
			valid_w(i, j, k) = 1;
		}		
    }

    for (int idx = 0; idx < u.size(); idx++){
        if(valid_u[idx] == 0){
            u[idx] == 0;
        }
    }
    for (int idx = 0; idx < v.size(); idx++){
        if(valid_v[idx] == 0){
            v[idx] == 0;
        }
    }
    for (int idx = 0; idx < w.size(); idx++){
        if(valid_w[idx] == 0){
            w[idx] == 0;
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
    for(int i=0; i<particles.size(); i++){
        particles[i] = traceRk2(particles[i], dt);
        // solid boundary
        double phi_val = phi_solid_ave(particles[i]/l + Vector3d(0.5, 0.5, 0.5));
        if(phi_val < 0){
            Vector3d grad;
            interpolate_gradient(grad, particles[i]/l + Vector3d(0.5, 0.5, 0.5), phi_solid);
            grad.normalize();
            particles[i] -= phi_val * grad;
        }
    }
    return;
}

void FluidSim::computePhi(){
    phi.setConstant(3*l);
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
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
		if (phi(i, j, k) < 0.5 * l && phi_solid_ave(i+0.5,j+0.5,k+0.5) < 0) {
			phi(i, j, k) = -0.5 * l;
		}
	}
    return;
}

void FluidSim::advect(const double& dt){
    semiLagrangian(u, u_new, dt, 1);
    semiLagrangian(v, v_new, dt, 2);
    semiLagrangian(w, w_new, dt, 3);
    u = u_new;
    v = v_new;
    w = w_new;
    return;
}

void FluidSim::semiLagrangian(const Field3d& field, Field3d& field_new, const double& dt, const int& id){
    assert(id == 1 || id == 2 || id == 3);
    Vector3d offset;
    if(id == 1){
        offset = Vector3d(0.5, 0, 0);
    }
    else if(id == 2){
        offset = Vector3d(0, 0.5, 0);
    }
    else if(id == 3){
        offset = Vector3d(0, 0, 0.5);
    }
    for (int i = 1; i < field.getN1()-1; i++) for (int j = 1; j < field.getN2()-1; j++) for (int k = 1; k < field.getN3()-1; k++) {
		Vector3d pos(i * l, j * l, k * l);
        Vector3d pos_field = traceRk2(pos - offset * l, -dt) / l + offset; 
        field_new(i, j, k) = interpolate_value(pos_field, field);        
	}
    return;
}

void FluidSim::applyForce(const double& dt){
    u.array() += -G * dt;
    return;
}

void FluidSim::computeWeights(){
    // compute weights_u
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        weights_u(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i, j+1, k), phi_solid(i, j+1, k+1), phi_solid(i, j, k+1));        
    }
    // compute weights_v
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3; k++){
        weights_v(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i+1, j, k), phi_solid(i+1, j, k+1), phi_solid(i, j, k+1));        
    }
    // compute weights_w
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3+1; k++){
        weights_w(i, j, k) = 1 - computeFraction(phi_solid(i, j, k), phi_solid(i, j+1, k), phi_solid(i+1, j+1, k), phi_solid(i+1, j, k));        
    }
    return;
}


void FluidSim::project(){
    // compute A and d
    Adiag.setZero(); Aplusi.setZero(); Aplusj.setZero(); Aplusk.setZero(); d.setZero();
    double theta = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        if(phi(i, j, k) < 0){
            // right
            if(weights_u(i+1, j, k) > 0){
                if(phi(i+1, j, k) < 0){
                    Adiag(i, j, k) += weights_u(i+1, j, k);
                    Aplusi(i, j, k) += -weights_u(i+1, j, k);
                }
                else{
                    theta = computeFraction(phi(i, j, k), phi(i+1, j, k));
                    if(theta < 0.01) theta = 0.01;
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
                    if(theta < 0.01) theta = 0.01;
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
                    if(theta < 0.01) theta = 0.01;
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
                    if(theta < 0.01) theta = 0.01;
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
                    if(theta < 0.01) theta = 0.01;
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
                    if(theta < 0.01) theta = 0.01;
                    Adiag(i, j, k) += weights_w(i, j, k) / theta;
                }
                d(i, j, k) += weights_w(i, j, k) * w(i, j, k);
            }            
        }        
    }
    // ensure A is positive definite
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if(Adiag(i, j, k)==0 && Aplusi(i, j, k)==0 && Aplusj(i, j, k)==0 && Aplusk(i, j, k)==0 && 
        (i==0 || Aplusi(i-1, j, k)==0) && (j==0 || Aplusj(i, j-1, k)==0) && (k==0 || Aplusk(i, j, k-1)==0)){
            Adiag(i, j, k) = 1;
            d(i, j, k) = 0;
        }		
	}
    // TODO: remove the 1D null space
    // solve A * p = d
    solve(10);
    // update u
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_u(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i - 1, j, k) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i - 1, j, k) > 0) theta = computeFraction(phi(i, j, k), phi(i - 1, j, k));
            if (theta < 0.01) theta = 0.01;
            u(i, j, k) += -((p(i, j, k) - p(i - 1, j, k))) / theta;
		}		
    }
    // update v
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) for (int k = 0; k < n3; k++){
        if (weights_v(i, j, k) > 0 && (phi(i, j, k) < 0 ||  phi(i, j - 1, k) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i, j - 1, k) > 0) theta = computeFraction(phi(i, j, k), phi(i, j - 1, k));
            if (theta < 0.01) theta = 0.01;
            v(i, j, k) += -((p(i, j, k) - p(i, j - 1, k))) / theta;
		}		
    }
    // update w
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 1; k < n3; k++){
        if (weights_w(i, j, k) > 0 && (phi(i, j, k) < 0 || phi(i, j, k - 1) < 0)) {
			theta = 1;
            if (phi(i, j, k) > 0 || phi(i, j, k - 1) > 0) theta = computeFraction(phi(i, j, k), phi(i, j, k - 1));
            if (theta < 0.01) theta = 0.01;
            w(i, j, k) += -((p(i, j, k) - p(i, j, k - 1))) / theta;
		}		
    }
    return;
}

void FluidSim::solve(const int& maxIterations){
    // set p = 0
	p.setConstant(0);
	if (d.cwiseAbs().maxCoeff() == 0) {		
		return;
	}
	// PCG
	//applyPrecon();
	z = d;
	s = z;
	double sigma = z.dot(d);
	double sigma_new;
	double alpha;
	for (int i = 0; i < maxIterations; i++) {
		applyA(s, z);
		alpha = sigma / z.dot(s);
		p += s * alpha;
		d -= z * alpha;
		if (d.cwiseAbs().maxCoeff() <= 1e-6)
			return;
		//applyPrecon();		
		z = d;
		sigma_new = z.dot(d);
		alpha = sigma_new / sigma;
		sigma = sigma_new;        		
        s = z + alpha * s;
	}
	
	return;
}

void FluidSim::constrain(){
    u_new = u;
    v_new = v;
    w_new = w;

    Vector3d pos, vec, normal_solid;
    // constrain u
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        if(weights_u(i, j, k) == 0){
            pos = Vector3d(i-0.5, j, k)*l;
            vec = getVelocity(pos);
            interpolate_gradient(normal_solid, pos/l + Vector3d(0.5, 0.5, 0.5), phi_solid);
            normal_solid.normalize();
            vec -= vec.dot(normal_solid) * normal_solid;
            u_new(i, j, k) = vec(0);
        }
    }
    // constrain v
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3; k++){
        if(weights_v(i, j, k) == 0){
            pos = Vector3d(i, j-0.5, k)*l;
            vec = getVelocity(pos);
            interpolate_gradient(normal_solid, pos/l + Vector3d(0.5, 0.5, 0.5), phi_solid);
            normal_solid.normalize();
            vec -= vec.dot(normal_solid) * normal_solid;
            v_new(i, j, k) = vec(1);
        }
    }
    // constrain w
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3+1; k++){
        if(weights_w(i, j, k) == 0){
            pos = Vector3d(i, j, k-0.5)*l;
            vec = getVelocity(pos);
            interpolate_gradient(normal_solid, pos/l + Vector3d(0.5, 0.5, 0.5), phi_solid);
            normal_solid.normalize();
            vec -= vec.dot(normal_solid) * normal_solid;
            w_new(i, j, k) = vec(2);
        }
    }

    u = u_new;
    v = v_new;
    w = w_new;
    return;
}

void FluidSim::setVelocity(const Vector3d& vec){
    u.setConstant(vec.x());
    v.setConstant(vec.y());
    w.setConstant(vec.z());
    return;
}

void FluidSim::run(double dt, int n){
    while (1) {				
		for (int i = 0; i < n; i++) {
			advance(dt);
		}
	}	
    return;
}
