#pragma once
#ifndef _VECTORMD_H_
#define _VECTORMD_H_

#include <cstdlib>
#include <vector>
#include <cmath>
#include <math.h>
#include <ctime>
#include <random>
#include <iterator>
#include <omp.h>
#include "Exception.h"

using namespace std;

template <typename T> class VectorMD;
template <typename T> ostream& operator<< (ostream&, const VectorMD<T>&);
template <typename T> istream& operator>> (istream&, VectorMD<T>&);

template <typename T>
class VectorMD
{
protected:
    vector<T> Vector;
    int size;
    double coord_MIN, coord_MAX;
public:
    VectorMD<T>();
    VectorMD<T>(int);
    ~VectorMD<T>();
    vector<T> get_Vector();
    //void set_Vector(double, double);
    const int get_Size() const;
    void set_Size(int);
    void init_Vector(int);
    void rand_Vector();
    VectorMD<T>(const VectorMD&);
    VectorMD<T> operator=(const VectorMD&);
    VectorMD<T> operator+(const VectorMD&);
    VectorMD<T> operator-(const VectorMD&);

    friend VectorMD operator*(const VectorMD&, T);

    VectorMD<T> operator/(T);
    VectorMD<T> operator*(T);
    T operator*(const VectorMD&);
    VectorMD<T> cross(const VectorMD&);
    VectorMD<T> operator+=(const VectorMD&);
    VectorMD<T> operator-=(const VectorMD&);
    VectorMD<T> operator*=(T);
    VectorMD<T> operator/=(T);
    friend ostream& operator<< <T>(ostream&, const VectorMD<T>&);
    friend istream& operator>> <T>(istream&, VectorMD<T>&);
    T distance(const VectorMD&);
    T module();
    void print_Vector();
};

template <typename T>
VectorMD<T>::VectorMD()
{
    this->set_Size(3);
    this->coord_MIN = 0.;
    this->coord_MAX = 0.;
}

template <typename T>
VectorMD<T>::VectorMD(int _size)
{
    this->size = _size;
    this->Vector.reserve(this->size);
    this->Vector.resize(this->size);
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        *it = (T)0;
}

template <typename T>
VectorMD<T>::~VectorMD()
{}

template <typename T>
vector<T> VectorMD<T>::get_Vector()
{ 
    return this->Vector; 
}

template <typename T>
const int VectorMD<T>::get_Size() const 
{ 
    return this->size; 
}

template <typename T>
void VectorMD<T>::set_Size(int _size)
{
    this->size = _size;
}

template <typename T>
void VectorMD<T>::init_Vector(int _size)
{
    this->set_Size(_size);
    this->Vector.reserve(this->size);
    this->Vector.resize(this->size);
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        *it = (T)0;
}

template <typename T>
void VectorMD<T>::rand_Vector()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> urd(0., 10.);
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        *it = urd(gen);
}

template <typename T>
VectorMD<T>::VectorMD(const VectorMD& _vector)
{
    this->set_Size(int(_vector.Vector.size()));
    this->Vector.resize(this->get_Size());
    auto _it = _vector.Vector.begin();
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
    {
        *it = *_it;
        ++_it;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator= (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.Vector.size())
            throw BadVectorDimensionException();
        auto _it = _vector.Vector.begin();
        for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        {
            *it = *_it;
            ++_it;
        }

        return *this;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator+ (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.get_Vector().size())
            throw BadVectorDimensionException();
        VectorMD<T>* temp = new VectorMD<T>(_vector.get_Vector().size());
        auto it = this->Vector.begin();
        auto _it = _vector.Vector.begin();
        for (auto it_temp = temp->Vector.begin(), end = temp->Vector.end(); it_temp != end; ++it_temp)
        {
            *it_temp = *it + *_it;
            ++it; ++_it;
        }

        return *temp;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator- (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.get_Vector().size())
            throw BadVectorDimensionException();
        VectorMD<T>* temp = new VectorMD<T>(_vector.get_Vector().size());
        auto it = this->Vector.begin();
        auto _it = _vector.Vector.begin();
        for (auto it_temp = temp->Vector.begin(), end = temp->Vector.end(); it_temp != end; ++it_temp)
        {
            *it_temp = *it - *_it;
            ++it; ++_it;
        }

        return *temp;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator* (T _multiply)
{
    VectorMD<T>* temp = new VectorMD<T>(this->get_Vector().size());
    auto _it = this->Vector.begin();
    for (auto it = temp->Vector.begin(), end = temp->Vector.end; it != end; ++it)
    {
        *it = *_it * _multiply;
        ++_it;
    }

    return *temp;
}

template <typename T>
VectorMD<T> VectorMD<T>::operator/ (T _division)
{
    try
    {
        if (_division == 0 || _division == 0.0)
            throw ZeroDivideException();
        VectorMD<T>* temp = new VectorMD<T>(this->get_Vector().size());
        auto _it = this->Vector.begin();
        for (auto it = temp->Vector.begin(), end = temp->Vector.end(); it != end; ++it)
        {
            *it = *_it / _division;
            ++_it;
        }

        return *temp;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
T VectorMD<T>::operator* (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.get_Vector().size())
            throw BadVectorDimensionException();
        T temp = 0;
        auto _it = _vector.Vector.begin();
        for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        {
            temp += *it * *_it;
            ++_it;
        }

        return temp;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::cross(const VectorMD& _vector)
{
    // рассмотреть отдельно
}

template <typename T>
VectorMD<T> VectorMD<T>::operator+= (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.get_Vector().size())
            throw BadVectorDimensionException();
        auto _it = _vector.Vector.begin();
        for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        {
            *it += *_it;
            ++_it;
        }

        return *this;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator-= (const VectorMD& _vector)
{
    try
    {
        if (this->get_Vector().size() != _vector.get_Vector().size())
            throw BadVectorDimensionException();
        auto _it = _vector.Vector.begin();
        for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        {
            *it -= *_it;
            ++_it;
        }

        return *this;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
VectorMD<T> VectorMD<T>::operator*= (T multiply)
{
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        *it *= multiply;

    return *this;
}

template <typename T>
VectorMD<T> VectorMD<T>::operator/= (T _division)
{
    try
    {
        if (_division == 0 || _division == 0.)
            throw ZeroDivideException();
        for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
            *it /= _division;

        return *this;
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 0;
    }
}

template <typename T>
ostream& operator<< (ostream& out, const VectorMD<T>& _vector)
{
    out << "( " << _vector.Vector.at(0);
    for (auto it = _vector.Vector.begin() + 1, end = _vector.Vector.end() - 1; it != end; ++it) /* */
        out << "; " << *it;
    out << "; " << _vector.Vector.at(_vector.Vector.size() - 1) << " )\n";
    return out;
}

template <typename T>
istream& operator>> (istream& in, VectorMD<T>& _vector)
{
    for (auto it = _vector.Vector.begin(), end = _vector.Vector.end(); it != end; ++it)
        in >> *it;
    return in;
}

template <typename T>
T VectorMD<T>::distance(const VectorMD& _vector)
{
    T temp = NULL;
    auto _it = _vector.Vector.begin();
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
    {
        temp += (*it - *_it) * (*it - *_it);
        ++_it;
    }

    return (T)(sqrt(temp));
}

template <typename T>
T VectorMD<T>::module()
{
    T temp = NULL;
    for (auto it = this->Vector.begin(), end = this->Vector.end(); it != end; ++it)
        temp += *it * *it;

    return (T)(sqrt(temp));
}

template <typename T>
void VectorMD<T>::print_Vector()
{
    cout << "( " << this->Vector[0];
    for (auto it = this->Vector.begin() + 1, end = this->Vector.end() - 1; it != end; ++it)
        cout << "; " << *it;
    cout << "; " << this->Vector.at(this->Vector.size() - 1) << " )\n";
}

template <typename T>
VectorMD<T> operator*(const VectorMD<T>& _vector, T _multiply)
{
    VectorMD<T>* temp = new VectorMD<T>(_vector.Vector.size());
    auto _it = _vector.Vector.begin();
    for (auto it = temp->Vector.begin(), end = temp->Vector.end(); it != end; ++it)
    {
        *it = *_it * _multiply;
        ++_it;
    }

    return *temp;
}
#endif // !_VECTORMD_H_