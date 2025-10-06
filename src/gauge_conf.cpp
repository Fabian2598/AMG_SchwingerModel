#include "gauge_conf.h"

std::complex<double> RandomU1() {
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	std::complex<double> z(cos(theta), sin(theta));
	return z;
}

void GaugeConf::initialize() {
	for (int i = 0; i < Ntot; i++) {
		for (int mu = 0; mu < 2; mu++) {
			Conf[i][mu] = RandomU1(); //Conf[Ns x Nt][mu in {0,1}]
		}
	}
}

void GaugeConf::read_conf(const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    int x, t, mu;
    double re, im; 
    while (infile >> x >> t >> mu >> re >> im) {
        Conf[Coords[x][t]][mu] = c_double(re, im); 
    }
    infile.close();
    std::cout << "Conf read from " << name << std::endl;
    
}

void GaugeConf::readBinary(const std::string& name){
    using namespace LV;
    std::ifstream infile(name, std::ios::binary);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    
    for (int x = 0; x < Nx; x++) {
    for (int t = 0; t < Nt; t++) {
        int n = x * Nx + t;
        for (int mu = 0; mu < 2; mu++) {
            int x_read, t_read, mu_read;
            double re, im;
            infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&re), sizeof(double));
            infile.read(reinterpret_cast<char*>(&im), sizeof(double));
                Conf[n][mu_read] = c_double(re, im);
           
        }
    }
    }
    infile.close();
      
}
