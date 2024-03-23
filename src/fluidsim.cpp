#include "fluidsim.h"

FluidSim::FluidSim(
    const int& n1_init,
    const int& n2_init,
    const int& n3_init,
    const double& l_init,
    const std::vector<double>& phi_init,
    const std::vector<double>& phi_solid_init){
        assert(phi_init.size() == n1_init*n2_init*n3_init);
        assert(phi_solid_init.size() == n1_init*n2_init*n3_init);        
        n1 = n1_init;
        n2 = n2_init;
        n3 = n3_init;
        l = l_init;       
        u.resize3D(n1+1, n2, n3);
        v.resize3D(n1, n2+1, n3);
        w.resize3D(n1, n2, n3+1);
        u_new.resize3D(n1+1, n2, n3);
        v_new.resize3D(n1, n2+1, n3);
        w_new.resize3D(n1, n2, n3+1);
        valid_u.resize3D(n1+1, n2, n3);
        valid_v.resize3D(n1, n2+1, n3);
        valid_w.resize3D(n1, n2, n3+1);
        valid_u_old.resize3D(n1+1, n2, n3);
        valid_v_old.resize3D(n1, n2+1, n3);
        valid_w_old.resize3D(n1, n2, n3+1);
        Adiag.resize3D(n1, n2, n3);
	    Aplusi.resize3D(n1, n2, n3);
	    Aplusj.resize3D(n1, n2, n3);
	    Aplusk.resize3D(n1, n2, n3);
        d.resize3D(n1, n2, n3);
        p.resize3D(n1, n2, n3);
        z.resize3D(n1, n2, n3);
        s.resize3D(n1, n2, n3);
        
        phi.resize3D(n1, n2, n3);
        phi_solid.resize3D(n1, n2, n3);
        for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
            // init phi and phi_solid
            phi(i, j, k) = phi_init[i*n2*n3+j*n3+k];
            phi_solid(i, j, k) = phi_solid_init[i*n2*n3+j*n3+k];            
        }

    // make the particles large enough so they always appear on the grid
    particle_radius = (l * 1.01 * sqrt(3.0) / 2.0);
    // init particles
    int seed = 0;
	for (int n = 0; n < 64*4; n++) {
		for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
			double a = randhashd(seed++); double b = randhashd(seed++); double c = randhashd(seed++);
			double x_ = i + a, y_ = j + b, z_ = k + c;
			if (interpolate_value(x_, y_, z_, phi) <= -particle_radius) {				
				if (interpolate_value(x_, y_, z_, phi_solid) > 0)
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

    extrapolate();
    advectParticles(dt);
    computePhi();
    advect(dt);
    applyForce(dt);
    project();
    constrain();
    return;
}

bool FluidSim::valid(const int& i, const int& j, const int& k){
    if (i >= 0 && i < n1 && j >= 0 && j < n2 && k >= 0 && k < n3 && phi_solid(i, j, k) > 0 && phi(i, j, k) < 0.0f) {
			return true;
    }
    return false;
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

double FluidSim::computeVolume(const int& i, const int& j, const int& k, const Field3d& phi, const int& id){
    /* id = 
	1: +x
	2: +y
	3: +z
	4: -x
	5: -y
	6: -z
	*/
    //--------------------------------------------------------------------- to be check ------------------------------------------------------------------
	double phi0, phi1, phi2, phi3;
	if (id == 1) {
		phi0 = interpolate_value(i + 0.5, j + 0.5, k + 0.5, phi);
		phi1 = interpolate_value(i + 0.5, j - 0.5, k + 0.5, phi);
		phi2 = interpolate_value(i + 0.5, j + 0.5, k - 0.5, phi);
		phi3 = interpolate_value(i + 0.5, j - 0.5, k - 0.5, phi);
	}
	else if (id == 2) {
		phi0 = interpolate_value(i + 0.5, j + 0.5, k + 0.5, phi);
		phi1 = interpolate_value(i - 0.5, j + 0.5, k + 0.5, phi);
		phi2 = interpolate_value(i + 0.5, j + 0.5, k - 0.5, phi);
		phi3 = interpolate_value(i - 0.5, j + 0.5, k - 0.5, phi);
	}
	else if (id == 3) {
		phi0 = interpolate_value(i + 0.5, j + 0.5, k + 0.5, phi);
		phi1 = interpolate_value(i - 0.5, j + 0.5, k + 0.5, phi);
		phi2 = interpolate_value(i + 0.5, j - 0.5, k + 0.5, phi);
		phi3 = interpolate_value(i - 0.5, j - 0.5, k + 0.5, phi);
	}
	else if (id == 4) {
		phi0 = interpolate_value(i - 0.5, j + 0.5, k + 0.5, phi);
		phi1 = interpolate_value(i - 0.5, j - 0.5, k + 0.5, phi);
		phi2 = interpolate_value(i - 0.5, j + 0.5, k - 0.5, phi);
		phi3 = interpolate_value(i - 0.5, j - 0.5, k - 0.5, phi);
	}
	else if (id == 5) {
		phi0 = interpolate_value(i + 0.5, j - 0.5, k + 0.5, phi);
		phi1 = interpolate_value(i - 0.5, j - 0.5, k + 0.5, phi);
		phi2 = interpolate_value(i + 0.5, j - 0.5, k - 0.5, phi);
		phi3 = interpolate_value(i - 0.5, j - 0.5, k - 0.5, phi);
	}
	else if (id == 6) {
		phi0 = interpolate_value(i + 0.5, j + 0.5, k - 0.5, phi);
		phi1 = interpolate_value(i - 0.5, j + 0.5, k - 0.5, phi);
		phi2 = interpolate_value(i + 0.5, j - 0.5, k - 0.5, phi);
		phi3 = interpolate_value(i - 0.5, j - 0.5, k - 0.5, phi);
	}
	else {
		std::cerr << "invalid input id" << std::endl;
		abort();
	}	
	// compute s1: 0 1 3
	double s1 = computeTriangleArea(phi0, phi1, phi3);
	// compute s2: 0 2 3
	double s2 = computeTriangleArea(phi0, phi2, phi3);
	return (s1 + s2) / 2;
}

void FluidSim::applyA(const Field3d& x, Field3d& ans){
    for (int i = 0; i < x.getN1(); i++) for (int j = 0; j < x.getN2(); j++) for (int k = 0; k < x.getN3(); k++) {
		ans(i, j, k) = Adiag(i, j, k) * x(i, j, k);		
		ans(i, j, k) += Aplusi(i, j, k) * x(i + 1, j, k);
		ans(i, j, k) += Aplusj(i, j, k) * x(i, j + 1, k);
		ans(i, j, k) += Aplusk(i, j, k) * x(i, j, k + 1);
		ans(i, j, k) += Aplusi(i - 1, j, k) * x(i - 1, j, k);
		ans(i, j, k) += Aplusj(i, j - 1, k) * x(i, j - 1, k);
		ans(i, j, k) += Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
    return;
}


void FluidSim::extrapolate(){
    // reset valid_u
    for (int i = 0; i < n1 + 1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if ((i < n1 && phi_solid(i, j, k) > 0 && phi(i, j, k) < 0.0) || (i - 1 >= 0 && phi_solid(i - 1, j, k) > 0 && phi(i - 1, j, k) < 0.0f)) {
			valid_u(i, j, k) = 1;
		}
		else {
			valid_u(i, j, k) = 0;
		}
    }
    // reset valid_v
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2 + 1; j++) for (int k = 0; k < n3; k++){
        if ((j < n2 && phi_solid(i, j, k) > 0 && phi(i, j, k) < 0.0) || (j - 1 >= 0 && phi_solid(i, j - 1, k) > 0 && phi(i, j - 1, k) < 0.0)) {
			valid_v(i, j, k) = 1;
		}
		else {
			valid_v(i, j, k) = 0;
		}
    }
    // reset valid_w
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3 + 1; k++){
        if ((k < n3 && phi_solid(i, j, k) > 0 && phi(i, j, k) < 0.0f) || (k - 1 >= 0 && phi_solid(i, j, k - 1) < 0 && phi(i, j, k - 1) < 0.0)) {
			valid_w(i, j, k) = 1;
		}
		else {
			valid_w(i, j, k) = 0;
		}
    }

	// Apply several iterations of a very simple propagation of valid velocity data in all directions
	extrapolate(u, u_new, valid_u, valid_u_old);
	extrapolate(v, v_new, valid_v, valid_v_old);
	extrapolate(w, w_new, valid_w, valid_w_old);
    return;
}

void FluidSim::extrapolate(Field3d& field, Field3d field_new, Field3i& valid, Field3i& valid_old){
    v_new = v;
	for (int num = 0; num < 10; num++) {
		float sum = 0;
		int count = 0;
		valid_old = valid;
		for (int i = 0; i < v.getN1(); i++) for (int j = 0; j < v.getN2(); j++) for (int k = 0; k < v.getN3(); k++) {
			sum = 0;
			count = 0;
			if (!valid_old(i, j, k)) {

				if (valid_old(i + 1, j, k)) {
					sum += v(i + 1, j, k);
					count++;
				}
				if (valid_old(i - 1, j, k)) {
					sum += v(i - 1, j, k);
					count++;
				}
				if (valid_old(i, j + 1, k)) {
					sum += v(i, j + 1, k);
					count++;
				}
				if (valid_old(i, j - 1, k)) {
					sum += v(i, j - 1, k);
					count++;
				}
				if (valid_old(i, j, k + 1)) {
					sum += v(i, j, k + 1);
					count++;
				}
				if (valid_old(i, j, k - 1)) {
					sum += v(i, j, k - 1);
					count++;
				}

				if (count > 0) {
					v_new(i, j, k) = sum / count;
					valid(i, j, k) = 1;
				}
			}
		}
		v = v_new;
	}
    return;
}

void FluidSim::advectParticles(const double& dt){
    for(int i=0; i<particles.size(); i++){
        particles[i] = traceRk2(particles[i], dt);
        // solid boundary
        double phi_val = interpolate_value(particles[i]/l, phi_solid);
        if(phi_val < 0){
            Vector3d grad;
            interpolate_gradient(grad, particles[i]/l, phi_solid);
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
		if (phi(i, j, k) < 0.5 * l && phi_solid(i,j,k) < 0) {
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
    else{
        abort();
    }
    for (int i = 0; i < field.getN1(); i++) for (int j = 0; j < field.getN2(); j++) for (int k = 0; k < field.getN3(); k++) {
		Vector3d pos(i * l, j * l, k * l);
        field_new(i, j, k) = interpolate_value(traceRk2(pos - offset * l, -dt) / l + offset, field);
	}
    return;
}

void FluidSim::applyForce(const double& dt){
    u.array() += -G * dt;
}

void FluidSim::project(){
    // compute A and d
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        Adiag(i, j, k) = computeVolume(i, j, k, phi, 1) + computeVolume(i, j, k, phi, 2) + computeVolume(i, j, k, phi, 3)
                       + computeVolume(i-1, j, k, phi, 1) + computeVolume(i, j-1, k, phi, 2) + computeVolume(i, j, k-1, phi, 3);
        Aplusi(i, j, k) = -computeVolume(i, j, k, phi, 1);
        Aplusj(i, j, k) = -computeVolume(i, j, k, phi, 2);
        Aplusk(i, j, k) = -computeVolume(i, j, k, phi, 3);
        d(i, j, k) = -computeVolume(i, j, k, phi, 1) * u(i+1, j, k) + computeVolume(i-1, j, k, phi, 1) * u(i, j, k)
                     -computeVolume(i, j, k, phi, 2) * v(i, j+1, k) + computeVolume(i, j-1, k, phi, 2) * v(i, j, k)
                     -computeVolume(i, j, k, phi, 3) * w(i, j, k+1) + computeVolume(i, j, k-1, phi, 3) * w(i, j, k);
    }
    // ensure A is positive definite
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++){
        if(Adiag(i, j, k)==0 && Aplusi(i, j, k)==0 && Aplusj(i, j, k)==0 && Aplusk(i, j, k)==0 && Aplusi(i-1, j, k)==0 && Aplusj(i, j-1, k)==0 && Aplusk(i, j, k-1)==0){
            Adiag(i, j, k) = 1;
            d(i, j, k) = 0;
        }		
	}
    // TODO: remove the 1D null space
    // solve A * p = d
    solve(10);
    // update u
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3; k++){
        if (!(i < n1 && phi_solid(i, j, k) < 0) && !(i != 0 && phi_solid(i - 1, j, k) < 0)) {// not solid
            if (i < n1 )// p+
                u(i, j, k) += -p(i, j, k);
            if (i != 0 )// p-
                u(i, j, k) += p(i - 1, j, k);
        }
    }
    // update v
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) for(int k=0; k<n3; k++){
        if (!(j < n2 && phi_solid(i, j, k) < 0) && !(j != 0 && phi_solid(i, j - 1, k) < 0)) {// not solid
            if (j < n2 )// p+
                v(i, j, k) += -p(i, j, k);
            if (j != 0 )// p-
                v(i, j, k) += p(i, j - 1, k);
        }
    }
    // update w
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) for(int k=0; k<n3+1; k++){
        if (!(k < n3 && phi_solid(i, j, k) < 0) && !(k != 0 && phi_solid(i, j, k - 1) < 0)) {// not solid
            if (k < n3 )// p+
                w(i, j, k) += -p(i, j, k);
            if (i != 0 )// p-
                w(i, j, k) += p(i, j, k - 1);
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
    // TODO    
}

void FluidSim::run(double dt, int n){
    while (1) {				
		for (int i = 0; i < n; i++) {
			advance(dt);
		}
	}	
}
