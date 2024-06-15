#include "geometry.h"

namespace backend {

Field3d computeSDF_cuboid(const int& n1, const int& n2, const int& n3, const Vector3i& left_corner, const Vector3i& cuboid_size) {
    Vector3i r[8];
    r[0] = left_corner;
    r[1] = r[0] + Vector3i(1, 0, 0) * cuboid_size(0);
    r[2] = r[1] + Vector3i(0, 1, 0) * cuboid_size(1);
    r[3] = r[0] + Vector3i(0, 1, 0) * cuboid_size(1);
    r[4] = r[0] + Vector3i(0, 0, 1) * cuboid_size(2);
    r[5] = r[4] + Vector3i(1, 0, 0) * cuboid_size(0);
    r[6] = r[5] + Vector3i(0, 1, 0) * cuboid_size(1);
    r[7] = r[4] + Vector3i(0, 1, 0) * cuboid_size(1);
    Vector3i pos;
    Vector3i d[8];
    Vector8d l;
    Vector6d h;
    Field3d phi(n1+1, n2+1, n3+1);

    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) for (int k = 0; k < n3+1; k++) {
        pos = Vector3i(i, j, k);
        for (int idx = 0; idx < 8; idx++) {
            d[idx] = pos - r[idx];
            l(idx) = (d[idx].cast<double>()).norm();     
        }
        h(0) = static_cast<double>(d[0](0));
        h(1) = static_cast<double>(d[1](0));
        h(2) = static_cast<double>(d[0](1));
        h(3) = static_cast<double>(d[3](1));
        h(4) = static_cast<double>(d[0](2));
        h(5) = static_cast<double>(d[4](2));        
        if (d[0](0) > 0 && d[1](0) < 0 && d[0](1) > 0 && d[3](1) < 0 && d[0](2) > 0 && d[4](2) < 0) {// inside
            phi(i, j, k) = -min(l.cwiseAbs().minCoeff(), h.cwiseAbs().minCoeff());
        }
        else {
            phi(i, j, k) = min(l.cwiseAbs().minCoeff(), h.cwiseAbs().minCoeff());
        }
    }
    
    return phi;
}

double computeSDF_cuboid(const Vector3d& pos, const Vector3d& cuboid_size) {
    Vector3d r[8];
    r[0] = -cuboid_size * 0.5;
    r[1] = r[0] + Vector3d(1, 0, 0) * cuboid_size(0);
    r[2] = r[1] + Vector3d(0, 1, 0) * cuboid_size(1);
    r[3] = r[0] + Vector3d(0, 1, 0) * cuboid_size(1);
    r[4] = r[0] + Vector3d(0, 0, 1) * cuboid_size(2);
    r[5] = r[4] + Vector3d(1, 0, 0) * cuboid_size(0);
    r[6] = r[5] + Vector3d(0, 1, 0) * cuboid_size(1);
    r[7] = r[4] + Vector3d(0, 1, 0) * cuboid_size(1);    
    Vector3d d[8];
    Vector8d l;
    Vector6d h;
    double phi;
    
    for (int idx = 0; idx < 8; idx++) {
        d[idx] = pos - r[idx];
        l(idx) = (d[idx].cast<double>()).norm();     
    }
    h(0) = d[0](0);
    h(1) = d[1](0);
    h(2) = d[0](1);
    h(3) = d[3](1);
    h(4) = d[0](2);
    h(5) = d[4](2);        
    if (d[0](0) > 0 && d[1](0) < 0 && d[0](1) > 0 && d[3](1) < 0 && d[0](2) > 0 && d[4](2) < 0) {// inside
        phi = -min(l.cwiseAbs().minCoeff(), h.cwiseAbs().minCoeff());
    }
    else {
        phi = min(l.cwiseAbs().minCoeff(), h.cwiseAbs().minCoeff());
    }    
    
    return phi;
}

Field3d computeSDF_sphere(const int n1, const int& n2, const int& n3, const Vector3d& center, const double& r) {
    Field3d phi(n1+1, n2+1, n3+1);
    Vector3d pos;
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) for (int k = 0; k < n3+1; k++) {
        pos = Vector3d(i, j, k);
        phi(i, j, k) = (pos - center).norm() - r;
    }
    return phi;
}

void writePLY(const std::string& ply_path, const std::vector<Vector3d>& vertices, const std::vector<Face>& faces) {
    tinyply::PlyFile plyFile;
	std::vector<float> vertexData;
	for (int i = 0; i < vertices.size(); i++) {
		vertexData.push_back(static_cast<float>(vertices[i](0)));
		vertexData.push_back(static_cast<float>(vertices[i](1)));
		vertexData.push_back(static_cast<float>(vertices[i](2)));
	}
	plyFile.add_properties_to_element("vertex", { "x", "y", "z" }, tinyply::Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertexData.data()), tinyply::Type::INVALID, 0);
	std::vector<int32_t> faceData;
	faceData.reserve(faces.size() * 3);
	for (int i = 0; i < faces.size(); i++) {		
		faceData.push_back(faces[i].v0);
		faceData.push_back(faces[i].v1);
		faceData.push_back(faces[i].v2);
	}
	plyFile.add_properties_to_element("face", { "vertex_indices" }, tinyply::Type::INT32, faces.size(), reinterpret_cast<uint8_t*>(faceData.data()), tinyply::Type::UINT8, 3);
	std::filebuf fb;
	fb.open(ply_path, std::ios::out | std::ios::binary);
	std::ostream outputStream(&fb);
	if (outputStream.fail())
		throw std::runtime_error("failed to open " + ply_path);
	plyFile.write(outputStream, true);

	return;
}

void writePLY(const std::string& ply_path, const std::vector<Vector3d>& particles) {
    tinyply::PlyFile plyFile;
	std::vector<float> vertexData;
	for (int i = 0; i < particles.size(); i++) {
		vertexData.push_back(static_cast<float>(particles[i](0)));
		vertexData.push_back(static_cast<float>(particles[i](1)));
		vertexData.push_back(static_cast<float>(particles[i](2)));
	}
	plyFile.add_properties_to_element("vertex", { "x", "y", "z" }, tinyply::Type::FLOAT32, particles.size(), reinterpret_cast<uint8_t*>(vertexData.data()), tinyply::Type::INVALID, 0);
	std::filebuf fb;
	fb.open(ply_path, std::ios::out | std::ios::binary);
	std::ostream outputStream(&fb);
	if (outputStream.fail())
		throw std::runtime_error("failed to open " + ply_path);
	plyFile.write(outputStream, true);

	return;
}

void writeNPY(const std::string& npy_path, const Field3d& phi) {
    cnpy::npy_save(npy_path, phi.data(), {(size_t)phi.getN1(), (size_t)phi.getN2(), (size_t)phi.getN3()});
    return;
}

void writeVideo(const std::string& file_name, const int& fps, const int& t) {
    std::string file_path = file_name + ".avi";
    std::string png_folder = projectPath + "/python/png";
    int codec = cv::VideoWriter::fourcc('M', 'J', 'P', 'G');
    cv::Size videoSize(1080, 1920);
    cv::VideoWriter video(file_name, codec, (double)fps, videoSize);
    Assert(video.isOpened(), "writeVideo", "can't write video");
    for (int i = 0; i < fps * t; i++) {
        std::string png_path = png_folder + "/liquid_" + std::to_string(i) + ".png";
        cv::Mat frame = cv::imread(png_path);
        if (!frame.empty()) {
            cv::resize(frame, frame, videoSize);
            video.write(frame);
        }
    }
    video.release();
    std::cout << "Video has be written." << std::endl;
    return;
}


 Field3d objToSDF(int N1, int N2, int N3, int size, double l, std::string obj_path) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = projectPath + "/python/sdf/temp.txt";
	// call python script	
	std::string command = pythonPath + " " + projectPath + "/python/obj2sdf.py"
		+ (" -N1 " + std::to_string(N1))
		+ " -N2 " + std::to_string(N2)
		+ " -N3 " + std::to_string(N3)
		+ " -size " + std::to_string(size)
		+ " -voxel_size " + std::to_string(l)
		+ " -input_file " + obj_path
		+ " -output_file " + file_path;
	system(command.c_str());
	// read sdf data
	std::ifstream file(file_path);
    Assert(file.is_open(), "objToSDF", "Can open the file.");
	file >> N1 >> N2 >> N3 >> l;
    Field3d phi(N1, N2, N3);
	double value;
	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
		if (!(file >> value)) {
			abort();
		}
		phi(i, j, k) = value;
	}
	file.close();
	std::cout << "Model is loaded.\n";
	// delet temp.txt	
	try {
		if (std::filesystem::exists(file_path)) { // check whether file exists
			std::filesystem::remove(file_path); // delete file
		}
		else {
			std::cout << "File does not exist." << std::endl;
		}
	}
	catch (const std::filesystem::filesystem_error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return phi;
}

bool objToSDF(int N1, int N2, int N3, double l, std::string obj_path, std::vector<double>& phi, double scale, double translate_x, double translate_y, double translate_z, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = projectPath + "/python/sdf/temp.txt";
	// call python script
	std::string command = pythonPath + " " + projectPath + "/python/obj2sdf.py"
		+ (" -N1 " + std::to_string(N1))
		+ " -N2 " + std::to_string(N2)
		+ " -N3 " + std::to_string(N3)
		+ " -scale " + std::to_string(scale)
		+ " -translate_x " + std::to_string(translate_x)
		+ " -translate_y " + std::to_string(translate_y)
		+ " -translate_z " + std::to_string(translate_z)
		+ " -voxel_size " + std::to_string(l)
		+ " -input_file " + obj_path
		+ " -output_file " + file_path;
	system(command.c_str());
	// read sdf data
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_path << std::endl;
		return false;
	}
	file >> N1 >> N2 >> N3 >> l;
	phi.resize(N1 * N2 * N3);
	double value;
	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
		if (!(file >> value)) {
			abort();
		}
		phi[i * N2 * N3 + j * N3 + k] = value;
		if (inverse) phi[i * N2 * N3 + j * N3 + k] *= -1;
	}
	file.close();
	std::cout << "Model is loaded.\n";
	// delet temp.txt	
	try {
		if (std::filesystem::exists(file_path)) { // check whether file exists
			std::filesystem::remove(file_path); // delete file
		}
		else {
			std::cout << "File does not exist." << std::endl;
		}
	}
	catch (const std::filesystem::filesystem_error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return true;
}

}