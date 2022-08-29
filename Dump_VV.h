#pragma once
#ifndef _DUMP_VV_H_
#define _DUMP_VV_H_
#include "Particle_VV.h"

class Particle_VV;
class Dump_VV
{
private:
	Particle_VV* ParticleMassive;
	int particle_COUNT;
	ofstream f_Data_rVector;
	ofstream f_Data_vVector;
	ofstream f_Data_e_calculation;
	ofstream f_Data_t_calculation;
	ofstream f_Data_p_calculation;
	ofstream f_Data_borderline;

	ofstream f_dx;
public:
	Dump_VV();
	Dump_VV(const Particle_VV&, int);
	Dump_VV(Particle_VV*, int);
	~Dump_VV();
	void set_ParticleMassive(Particle_VV&);
	void set_ParticleMassive(Particle_VV*);
	void set_particle_COUNT(int);
	void print_rVector();
	void print_vVector();
	void print_e_Calculation(double, double, double);
	void print_t_Calculation(double);
	void print_p_Calculation(double);
	void print_Borderline(double);
	void print_Bline(double, double);

	void print_dx(double, double, double, double, double);

	double random(double, double);
};
#endif // !_DUMP_VV_H_#pragma once
