#pragma once
#ifndef _PARTICLE_VV_H_
#define _PARTICLE_VV_H_
#include "Vector3D.h"

class Particle_VV
{
protected:
	int particle_COUNT;
	int particle_ID;
	double particle_Mass;
	double particle_Charge;
	Vector3D<double>* radius_Vector;
	Vector3D<double>* velocity_Vector;
	Vector3D<double>* acceleration_Vector;
	Vector3D<double>* velocity_Vector_05t;

	int flag;
public:
	Particle_VV();
	~Particle_VV();
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
	void set_Radius_Vector_Circle(double, double, double);

	Vector3D<double> get_Velocity_Vector();
	void set_Velocity_Vector(double, double);
	void set_Velocity_Vector(Vector3D<double>);
	void set_Velocity_Vector(double, double, double);

	Vector3D<double> get_Acceleration_Vector();
	void set_Acceleration_Vector(double, double);
	void set_Acceleration_Vector(double, double, double);

	Vector3D<double> get_Velocity_Vector_05t();
	void set_Velocity_Vector_05t(double, double);
	void set_Velocity_Vector_05t(Vector3D<double>);
	void set_Velocity_Vector_05t(double, double, double);

	int get_Flag();
	void set_Flag(int);

	void set_Particle_Count(const int);
};
#endif // !_PARTICLE_VV_H_