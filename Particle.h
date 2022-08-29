#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include "Vector3D.h"

class Particle
{
protected:
	int particle_COUNT;
	int particle_ID;
	double particle_Mass;
	double particle_Charge;
	Vector3D<double>* radius_Vector;
	Vector3D<double>* velocity_Vector;
	Vector3D<double>* acceleration_Vector;

	Vector3D<double>* previous_rVector;			// радиус-вектор на i-1 шаге итерации
	Vector3D<double>* previous_vVector;			// вектор скоростей на i-1 шаге итерации
public:
	Particle();
	~Particle();
	int get_Particle_ID();
	void set_Particle_ID(int);
	
	double get_Particle_Mass();
	void set_Particle_Mass(double);

	double get_Particle_Charge();
	void set_Particle_Charge(double);

	Vector3D<double> get_Radius_Vector();
	void set_Radius_Vector(double, double);
	void set_Radius_Vector(double, double, double);
	void set_Radius_Vector(Vector3D<double>);
	//void set_rVector(double, double, double);

	Vector3D<double> get_Velocity_Vector();
	void set_Velocity_Vector(double, double);
	void set_Velocity_Vector(Vector3D<double>);
	void set_Velocity_Vector(double, double, double);

	Vector3D<double> get_previous_rVector();
	void set_previous_rVector(Vector3D<double>);
	void set_previous_rVector(double, double, double);
	void set_previous_vVector(double, double, double);

	Vector3D<double> get_Acceleration_Vector();
	void set_Acceleration_Vector(double, double);
	void set_aVector(double, double, double);

	void set_Particle_Count(const int);
};
#endif // !_PARTICLE_H_