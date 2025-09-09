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

