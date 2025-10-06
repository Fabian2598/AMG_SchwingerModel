#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "variables.h"
#include <fstream>

/*
Generate a random U(1) variable
*/
std::complex<double> RandomU1(); 


class GaugeConf {
public:
	/*
	Nspace: number of lattice points in the space direction
	Ntime: number of lattice points in the time direction
	*/
	GaugeConf(const int& Nspace, const int& Ntime) : Nx(Nspace), Nt(Ntime), Ntot(Nspace* Ntime) {
		Conf = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //Gauge configurationion copy
		
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) : Nx(GConfig.getNx()), Nt(GConfig.getNt()), Ntot(Nx*Nt) {
		Conf = GConfig.Conf; 
	}
	~GaugeConf() {}; 

	/*
	Random initialization of the gauge configuration
	*/
	void initialize(); 

	/*
	Set the gauge configuration
	*/
	void setGconf(const std::vector<std::vector<std::complex<double>>>& CONF) {Conf = CONF;}

	int getNx() const { return Nx; }
	int getNt() const { return Nt; }

	std::vector<std::vector<std::complex<double>>> Conf; //Conf[Nx Nt][2] 	

	void read_conf(const std::string& name); //read gauge configuration from .txt file
	void readBinary(const std::string& name);
private:
	int Nx, Nt, Ntot;
};




#endif
