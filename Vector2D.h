#pragma once
#ifndef _VECTOR2D_H_
#define _VECTOR2D_H_
#include "VectorMD.h"

template <typename T>
class Vector2D : public VectorMD<T>
{
private:
	double coord_MIN, coord_MAX;
public:
	Vector2D<T>();
	Vector2D<T>(int);
	~Vector2D<T>();
	void set_Coord_MIN(double);
	void set_Coord_MAX(double);
	void set_Vector(double, double);
	void rand_Vector();
};

template <typename T>
Vector2D<T>::Vector2D() : VectorMD<T>::VectorMD()
{
	this->coord_MIN = 0.;
	this->coord_MAX = 0.;
}

template <typename T>
Vector2D<T>::Vector2D(int _size) : VectorMD<T>::VectorMD(_size)
{
	this->coord_MIN = 0.;
	this->coord_MAX = 0.;
}

template <typename T>
Vector2D<T>::~Vector2D()
{}
// метод set для левой границы диапазона координат вектора
template <typename T>
void Vector2D<T>::set_Coord_MIN(double _coord_MIN)
{
	this->coord_MIN = _coord_MIN;
}
// метод set для правой границы диапазона координат вектора
template <typename T>
void Vector2D<T>::set_Coord_MAX(double _coord_MAX)
{
	this->coord_MAX = _coord_MAX;
}
// метод set изменения координат вектора VectorMD<T>::Vector
template <typename T>
void Vector2D<T>::set_Vector(double _coord_X, double _coord_Y)
{
	this->Vector[0] = _coord_X;
	this->Vector[1] = _coord_Y;
}
// метод генерации равномерно распределенных случайных чисел типа <double> и заполнение VectorMD<T>::Vector ими
template <typename T>
void Vector2D<T>::rand_Vector()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> urd(this->coord_MIN, this->coord_MAX);
	for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
		*it = urd(gen);
}
#endif // !_VECTOR2D_H_