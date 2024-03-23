#include "geometry.h"

bool objToSDF(int N1, int N2, int N3, int size, double l, std::string obj_path, std::vector<double>& phi, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = projectPath + "/python/sdf/temp.txt";
	// call python script	
	std::string command = "/home/ricky/anaconda3/envs/myenv/bin/python " + projectPath + "/python/obj2sdf.py"
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

bool objToSDF(int N1, int N2, int N3, double l, std::string obj_path, std::vector<double>& phi, double scale, double translate_x, double translate_y, double translate_z, bool inverse) {
	std::cout << "Loading " + obj_path << "\n";
	std::string file_path = projectPath + "/python/sdf/temp.txt";
	// call python script
	std::string command = "/home/ricky/anaconda3/envs/myenv/bin/python " + projectPath + "/python/obj2sdf.py"
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