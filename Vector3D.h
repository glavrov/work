#pragma once
#define _USE_MATH_DEFINES
#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_
#include <iomanip>
#include <math.h>
#include "VectorMD.h"

template <typename T>
class Vector3D : public VectorMD<T>
{
private:
	double coord_MIN, coord_MAX;
public:
	Vector3D<T>();
	Vector3D<T>(int);
	~Vector3D<T>();
	void set_Coord_MIN(double);
	void set_Coord_MAX(double);
	T get_Coord_MIN();
	T get_Coord_MAX();
	void set_Vector(double, double, double);
	void rand_Vector();
	void rand_Vector(double, double, double);
	double random(double, double);
};
// конструкторы
template <typename T>
Vector3D<T>::Vector3D() : VectorMD<T>::VectorMD()
{
	this->set_Coord_MIN(0.);
	this->set_Coord_MAX(0.);
}

template <typename T>
Vector3D<T>::Vector3D(int _size) : VectorMD<T>::VectorMD(_size)
{
	this->set_Coord_MIN(0.);
	this->set_Coord_MAX(0.);
}
// деструктор
template <typename T>
Vector3D<T>::~Vector3D()
{}
// метод set для левой границы диапазона координат вектора
template <typename T>
void Vector3D<T>::set_Coord_MIN(double _coord_MIN)
{
	this->coord_MIN = _coord_MIN;
}
// метод set для правой границы диапазона координат вектора
template <typename T>
void Vector3D<T>::set_Coord_MAX(double _coord_MAX)
{
	this->coord_MAX = _coord_MAX;
}
// метод get для получения левой границы диапазона координат вектора
template <typename T>
T Vector3D<T>::get_Coord_MIN()
{
	return this->coord_MIN;
}
// метод get для получения правой границы диапазона координат вектора
template <typename T>
T Vector3D<T>::get_Coord_MAX()
{
	return this->coord_MAX;
}
// метод set изменения координат вектора VectorMD<T>::Vector
template <typename T>
void Vector3D<T>::set_Vector(double _coord_X, double _coord_Y, double _coord_Z)
{
	this->Vector.at(0) = _coord_X;
	this->Vector.at(1) = _coord_Y;
	this->Vector.at(2) = _coord_Z;
}
// метод генерации равномерно распределенных случайных чисел типа <double> и заполнение VectorMD<T>::Vector ими (для случая квадратной области)
template <typename T>
void Vector3D<T>::rand_Vector()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> urd(this->get_Coord_MIN(), this->get_Coord_MAX());
//#pragma omp parallel for
	for (int i = 0; i < int(this->Vector.size()); i++)
	{
		if (i != 2)
			this->Vector.at(i) = urd(gen);
		else
			this->Vector.at(i) = 0.;
	}
}
// метод генерации равномерно распределенных случайных чисел типа <double> и заполнение VectorMD<T>::Vector ими (для случая круговой области)
template <typename T>
void Vector3D<T>::rand_Vector(double _center_X, double _center_Y, double _R)
{
	/* double R = _R * this->random(0., 1.);
	//double R = _R * sqrt(this->random(0., 1.));
	double theta =  2. * M_PI * this->random(0., 1.);
	//this->Vector[0] = _center_X + R * cos(theta);
	//this->Vector[1] = _center_Y + R * sin(theta);
	this->Vector.at(0) = _center_X + R * cos(theta);
	this->Vector.at(1) = _center_Y + R * sin(theta);
	this->Vector.at(2) = 0.; */
	
	double u = this->random(0., 1.) + this->random(0., 1.);
	double theta = 2. * M_PI * this->random(0., 1.);
	double R;
	if (u > 1) R = _R * (2. - u);
	else R = _R * u;
	this->Vector.at(0) = _center_X + R * cos(theta);
	this->Vector.at(1) = _center_Y + R * sin(theta);
	this->Vector.at(2) = 0.;
}
// метод генерации равномерно распределенных случайных чисел типа <double> для получения полярного радиуса и угла (для случая круговой области)
template <typename T>
double Vector3D<T>::random(double _min, double _max)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> urd(_min, _max);

	double random = urd(gen);
	return random;
}
#endif // !_VECTOR3D_H_