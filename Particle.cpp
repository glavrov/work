#include "Particle.h"

Particle::Particle()
{
	this->particle_ID = 0;
	this->particle_Mass = 1.;
	this->particle_Charge = 0.;
	this->particle_COUNT = 0;

	this->radius_Vector = new Vector3D<double>(3);
	this->velocity_Vector = new Vector3D<double>(3);
	this->acceleration_Vector = new Vector3D<double>(3);

	this->previous_rVector = new Vector3D<double>(3);
	this->previous_vVector = new Vector3D<double>(3);
}

Particle::~Particle()
{}
// метод get дл€ пол€ particle_ID
int Particle::get_Particle_ID()
{
	return this->particle_ID;
}
// метод set дл€ пол€ particle_ID
void Particle::set_Particle_ID(int _particle_ID)
{
	this->particle_ID = _particle_ID;
}
// метод get дл€ пол€ particle_Mass
double Particle::get_Particle_Mass()
{
	return this->particle_Mass;
}
// метод set дл€ пол€ particle_Mass
void Particle::set_Particle_Mass(double _particle_Mass)
{
	this->particle_Mass = _particle_Mass;
}
// метод get дл€ пол€ particle_Charge
double Particle::get_Particle_Charge()
{
	return particle_Charge;
}
// метод set дл€ пол€ particle_Charge
void Particle::set_Particle_Charge(double _particle_Charge)
{
	this->particle_Charge = _particle_Charge;
}
// метод get дл€ указател€ radius_Vector на эклемпл€р класса Vector3D<double>
Vector3D<double> Particle::get_Radius_Vector()
{
	return *this->radius_Vector;
}
// метод set дл€ указател€ radius_Vector на эклемпл€р класса Vector3D<double>
void Particle::set_Radius_Vector(double _coord_MIN, double _coord_MAX)
{
	this->radius_Vector->set_Coord_MIN(_coord_MIN);
	this->radius_Vector->set_Coord_MAX(_coord_MAX);
	this->radius_Vector->rand_Vector();
}
// метод set дл€ указател€ radius_Vector на эклемпл€р класса Vector3D<double>
void Particle::set_Radius_Vector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->radius_Vector->set_Vector(_coord_X, _coord_Y, _coord_Z);
}
// метод set дл€ указател€ radius_Vector на эклемпл€р класса Vector3D<double>
void Particle::set_Radius_Vector(Vector3D<double> _radius_Vector)
{
	*this->radius_Vector = _radius_Vector;
}
// метод get дл€ указател€ velocity_Vector на эклемпл€р класса Vector3D<double>
Vector3D<double> Particle::get_Velocity_Vector()
{
	return *this->velocity_Vector;
}
// метод set дл€ указател€ velocity_Vector на эклемпл€р класса Vector3D<double>
void Particle::set_Velocity_Vector(double _velocity_MIN, double _velocity_MAX)
{
	this->velocity_Vector->set_Coord_MIN(_velocity_MIN);
	this->velocity_Vector->set_Coord_MAX(_velocity_MAX);
	this->velocity_Vector->rand_Vector();
}
// метод set дл€ указател€ velocity_Vector на эклемпл€р класса Vector3D<double>
void Particle::set_Velocity_Vector(Vector3D<double> _velocity_Vector)
{
	*this->velocity_Vector = _velocity_Vector;
}
// метод set дл€ изменени€ координат вектора скорости
void Particle::set_Velocity_Vector(double _velocity_X, double _velocity_Y, double _velocity_Z)
{
	this->velocity_Vector->set_Vector(_velocity_X, _velocity_Y, _velocity_Z);
}
// метод get дл€ пол€ previous_rVector
Vector3D<double> Particle::get_previous_rVector()
{
	return *this->previous_rVector;
}
// метод set дл€ указател€ previous_rVector на эклемпл€р класса Vector3D<double>
void Particle::set_previous_rVector(Vector3D<double> _previuos_rVector)
{
	*this->previous_rVector = _previuos_rVector;
}
// метод set дл€ изменени€ координат радиус-вектора
void Particle::set_previous_rVector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->previous_rVector->set_Vector(_coord_X, _coord_Y, _coord_Z);
}
// метод set дл€ изменени€ координат вектора скорости
void Particle::set_previous_vVector(double _velocity_X, double _velocity_Y, double _velocity_Z)
{
	this->previous_vVector->set_Vector(_velocity_X, _velocity_Y, _velocity_Z);
}
// метод set дл€ пол€ particle_COUNT
void Particle::set_Particle_Count(const int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
}
// метод get дл€ пол€ acceleration_Vector
Vector3D<double> Particle::get_Acceleration_Vector()
{
	return *this->acceleration_Vector;
}
//
void Particle::set_Acceleration_Vector(double _acceleration_MIN, double _acceleration_MAX)
{
	this->acceleration_Vector->set_Coord_MIN(_acceleration_MIN);
	this->acceleration_Vector->set_Coord_MAX(_acceleration_MAX);
	this->acceleration_Vector->rand_Vector();
}
// метод set дл€ изменени€ координат вектора ускорени€
void Particle::set_aVector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->acceleration_Vector->set_Vector(_coord_X, _coord_Y, _coord_Z);
}