//-----------------------------------------------------------------------------------------------------
// Моделирование движения частиц в ограниченной квадратной (кубической) области [coord_MIN; coord_MAX],
// границы которой являются "прозрачными" для частиц, и действуют периодические граничные условия.
// Моделирование осуществляется классическим методом Верле.
//-----------------------------------------------------------------------------------------------------

#pragma once
#ifndef _MD_H_
#define _MD_H_
#include <fstream>
#include <iomanip>
#include "Particle.h"
#include "Dump.h"

class MD
{
private:
	Particle* ParticleMassive;					// массив указателей на объект типа Particle
	Dump* dump_FILE;							// указатель на объект типа Dump (для печати результатов моделирования в файлы)
	Vector3D<double>* SumV;						// массив указателей на объект типа Vector3D<double> (для расчета средней скорости)
	Vector3D<double>* SumV2;					// массив указателей на объект типа Vector3D<double> (для расчета среднеквадратической скорости)
	Vector3D<double>* velocity_Coefficient;		// массив указателей на объект типа Vector3D<double> (для расчета коэффициента нормализации)
	int particle_COUNT;							// количество частиц
	double particle_Temperature;				// температура
	double coord_MIN;							// минимальная координата
	double coord_MAX;							// максимальная координата
	double velocity_MIN;						// минимальная скорость
	double velocity_MAX;						// максимальная скорость
	double delta_T;								// шаг по времени
	double k_B = 1.;							// постоянная Больцмана
	int coeff_Ax;
	int coeff_Ay;

	double E_kinetic, E_potential, E_full;		// кинетическая, потенциальная и полная энергии
public:
	MD();										// конструкторы
	MD(int);
	MD(int, double);
	MD(int, double, double, double, double);
	~MD();										// деструктор

	void Initialization();						// метод инициализации массива частиц
	void velocity_Normalization();				// метод нормализации скоростей частиц
	void Move_СV(int);							// метод моделирования движения частиц
	void Evolution_СV(int, int);				// метод моделирования эволюции системы
	
	Vector3D<double> force(Vector3D<double>&);	// метод расчета сил, действующих на частицу
	double LJ_potential(Vector3D<double>&);		// метод расчета потенциала Леннард-Джонса
	void Next_Radius_Vector(int);				// метод расчета радиус-векторов частиц в момент t+Δt
	void Next_Velocity_Vector();				// метод расчета скоростей частиц в момент t+Δt
	void compute_Acceleration();				// метод расчета ускорений частиц
	void v_Normalization();						// метод нормализации скоростей частиц

	double Average_velocity(int);				// метод расчета средней скорости частиц
	double Average_velocity2(int);				// метод расчета среднеквадратической скорости частиц

	void E_calculation();						// метод расчета кинетической, потенциальной и полной энергий
	double lambda_Thermostat(double, double);				// метод расчета параметра лямбда для термостата Берендсена
	double lambda_Barostat(double, double);					// метод расчета параметра лямбда для баростата Берендсена
	double lambda_Proportional_Thermostat(double, double);	// метод расчета параметра лямбда для пропорционального термостата 
};
#endif // !_MD_H_