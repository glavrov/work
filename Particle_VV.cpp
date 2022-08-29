#include "Particle_VV.h"
// конструктор
Particle_VV::Particle_VV()
{
	this->particle_ID = 0;
	this->particle_Mass = 1.;
	this->particle_Charge = 0.;
	this->particle_COUNT = 0;
	this->flag = 0;

	this->radius_Vector = new Vector3D<double>(3);
	this->velocity_Vector = new Vector3D<double>(3);
	this->acceleration_Vector = new Vector3D<double>(3);
	this->velocity_Vector_05t = new Vector3D<double>(3);
}
// деструктор
Particle_VV::~Particle_VV()
{}
// метод get для поля Particle_ID
int Particle_VV::get_Particle_ID()
{
	return this->particle_ID;
}
// метод set для поля particle_ID
void Particle_VV::set_Particle_ID(int _particle_ID)
{
	this->particle_ID = _particle_ID;
}
// метод get для поля particle_Mass
double Particle_VV::get_Particle_Mass()
{
	return this->particle_Mass;
}
// метод set для поля particle_Mass
void Particle_VV::set_Particle_Mass(double _particle_Mass)
{
	this->particle_Mass = _particle_Mass;
}
// метод get для поля particle_Charge
double Particle_VV::get_Particle_Charge()
{
	return this->particle_Charge;
}
// метод set для поля particle_Charge
void Particle_VV::set_Particle_Charge(double _particle_Charge)
{
	this->particle_Charge = _particle_Charge;
}
// метод get для указателя radius_Vector на эклемпляр класса Vector3D<double>
Vector3D<double> Particle_VV::get_Radius_Vector()
{
	return *this->radius_Vector;
}
// метод set для указателя radius_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Radius_Vector(double _coord_MIN, double _coord_MAX)
{
	this->radius_Vector->set_Coord_MIN(_coord_MIN);
	this->radius_Vector->set_Coord_MAX(_coord_MAX);
	this->radius_Vector->rand_Vector();
}
// метод set для указателя radius_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Radius_Vector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->radius_Vector->set_Vector(_coord_X, _coord_Y, _coord_Z);
}
// метод set для указателя radius_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Radius_Vector(Vector3D<double> _radius_Vector)
{
	*this->radius_Vector = _radius_Vector;
}
//
void Particle_VV::set_Radius_Vector_Circle(double _center_X, double _center_Y, double _R)
{
	this->radius_Vector->rand_Vector(_center_X, _center_Y, _R);
}
// метод get для указателя velocity_Vector на эклемпляр класса Vector3D<double>
Vector3D<double> Particle_VV::get_Velocity_Vector()
{
	return *this->velocity_Vector;
}
// метод set для указателя velocity_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Velocity_Vector(double _velocity_MIN, double _velocity_MAX)
{
	this->velocity_Vector->set_Coord_MIN(_velocity_MIN);
	this->velocity_Vector->set_Coord_MAX(_velocity_MAX);
	this->velocity_Vector->rand_Vector();
}
// метод set для указателя velocity_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Velocity_Vector(Vector3D<double> _velocity_Vector)
{
	*this->velocity_Vector = _velocity_Vector;
}
// метод set для изменения координат вектора скорости
void Particle_VV::set_Velocity_Vector(double _velocity_X, double _velocity_Y, double _velocity_Z)
{
	this->velocity_Vector->set_Vector(_velocity_X, _velocity_Y, _velocity_Z);
}
// метод get для поля acceleration_Vector
Vector3D<double> Particle_VV::get_Acceleration_Vector()
{
	return *this->acceleration_Vector;
}
//
void Particle_VV::set_Acceleration_Vector(double _acceleration_MIN, double _acceleration_MAX)
{
	this->acceleration_Vector->set_Coord_MIN(_acceleration_MIN);
	this->acceleration_Vector->set_Coord_MAX(_acceleration_MAX);
	this->acceleration_Vector->rand_Vector();
}
// метод set для изменения координат вектора ускорения
void Particle_VV::set_Acceleration_Vector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->acceleration_Vector->set_Vector(_coord_X, _coord_Y, _coord_Z);
}
// метод get для указателя velocity_Vector_05t на экземпляр класса Vector3D<double>
Vector3D<double> Particle_VV::get_Velocity_Vector_05t()
{
	return *this->velocity_Vector_05t;
}
// метод set для указателя velocity_Vector_05t на экземпляр класса Vector3D<double>
void Particle_VV::set_Velocity_Vector_05t(double _velocity_MIN, double _velocity_MAX)
{
	this->velocity_Vector_05t->set_Coord_MIN(_velocity_MIN);
	this->velocity_Vector_05t->set_Coord_MAX(_velocity_MAX);
	this->velocity_Vector_05t->rand_Vector();
}
// метод set для указателя velocity_Vector на эклемпляр класса Vector3D<double>
void Particle_VV::set_Velocity_Vector_05t(Vector3D<double> _velocity_Vector_05t)
{
	*this->velocity_Vector_05t = _velocity_Vector_05t;
}
// метод set для изменения координат вектора скорости
void Particle_VV::set_Velocity_Vector_05t(double _velocity_X, double _velocity_Y, double _velocity_Z)
{
	this->velocity_Vector_05t->set_Vector(_velocity_X, _velocity_Y, _velocity_Z);
}
// метод get для поля flag
int Particle_VV::get_Flag()
{
	return this->flag;
}
// метод set для изменения поля flag
void Particle_VV::set_Flag(int _flag)
{
	this->flag = _flag;
}
// метод set для поля particle_COUNT
void Particle_VV::set_Particle_Count(const int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
}