#pragma once
#ifndef _DUMP_H_
#define _DUMP_H_
#include "Particle.h"
#include "Particle_VV.h"

class Particle;
class Dump
{
private:
	Particle* ParticleMassive;
	int particle_COUNT;
	ofstream f_Data_rVector;
	ofstream f_Data_vVector;
	ofstream f_Data_e_calculation;
	ofstream f_Data_t_calculation;
	ofstream f_Data_borderline;
public:
	Dump();
	Dump(const Particle&, int);
	Dump(Particle*, int);
	~Dump();
	void set_ParticleMassive(Particle&);
	void set_ParticleMassive(Particle*);
	void set_particle_COUNT(int);
	void print_rVector();
	void print_vVector();
	void print_t_Calculation(double);
	void print_e_Calculation(double, double, double);
	void print_Borderline(double, double);
};
#endif // !_DUMP_H_