//-----------------------------------------------------------------------------------------------------
// Моделирование движения частиц в ограниченной квадратной (кубической) области [coord_MIN; coord_MAX],
// границы которой являются "прозрачными" для частиц, и действуют периодические граничные условия,
// либо границы являются "непрозрачными", и действует закон отражения.
// Моделирование осуществляется скоростным методом Верле.
//-----------------------------------------------------------------------------------------------------

#pragma once
#ifndef _MD_VV_H_
#define _MD_VV_H_
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Particle_VV.h"
#include "Dump_VV.h"

class MD_VV
{
private:
	Particle_VV* ParticleMassive;							// массив указателей на объект типа Particle
	Dump_VV* dump_FILE;										// указатель на объект типа Dump
	Vector3D<double>* SumV;									// массив указателей на объект типа Vector3D<double> (расчет средней координаты)
	Vector3D<double>* SumV2;								// массив указателей на объект типа Vector3D<double> (расчет среднеквадратичной координаты)
	Vector3D<double>* velocity_Coefficient;					// массив указателей на объект типа Vector3D<double> (расчет коэффициента нормализации)
	int particle_COUNT;										// количество частиц
	double system_Temperature = 40.;						// температура
	double coord_MIN = 0.;									// минимальная координата
	double coord_MAX = 10.;									// максимальная координата
	double velocity_MIN = -1.;								// минимальная скорость
	double velocity_MAX = 1.;								// максимальная скорость
	double delta_T = 0.;									// шаг по времени
	double k_B = 1.;										// постоянная Больцмана
	double particle_RADIUS = .65;							// 1. // sqrt(3); // Rmin // pow(2., (1. / 6.));
	int space_TYPE = 0;										// тип пространства
	
	double r_Cut = (this->coord_MAX - this->coord_MIN) / 2.;
	double displ = .0001;
	int m_Size = 50000; //(int)(r_Cut / displ);
	double* table_pot = new double[m_Size];

	double E_kinetic, E_potential, E_full;					// кинетическая, потенциальная, полная энергии
public:
	MD_VV();												// конструкторы
	MD_VV(int);
	MD_VV(int, double);
	MD_VV(int, double, double, double, double);
	~MD_VV();												// деструктор

	void Initialization();									// метод инициализации массива частиц, квадратное пространство
	void Initialization(double);							// метод инициализации массива частиц, круглое пространство
	void Velocity_Normalization();							// метод нормализации скоростей частиц
	void Move_VV(int, int);									// метод моделирования движения частиц
	void Evolution_VV(int, int, int);						// метод моделирования эволюции системы

	Vector3D<double> force(Vector3D<double>&);				// метод расчета сил, действующих на частицу
	double LJ_potential(Vector3D<double>&);					// метод расчета потенциала Леннард-Джонса
	void Next_Radius_Vector();								// метод расчета радиус-векторов частиц в момент t+Δt
	void Next_Radius_Vector(int, int);						// метод расчета радиус-векторов частиц в момент t+Δt
	void Next_Velocity_Vector_05t();						// метод расчета скоростей частиц в момент t+Δt/2
	void Next_Velocity_Vector();							// метод расчета скоростей частиц в момент t+Δt
	void Compute_Acceleration();							// метод расчета ускорений частиц
	void v_Normalization();									// метод нормализации скоростей частиц

	double Average_velocity(int);							// метод расчета средней скорости частиц
	double Average_velocity2(int);							// метод расчета среднеквадратической скорости частиц

	void E_calculation_VV();								// метод расчета кинетической, потенциальной и полной энергий
	double lambda_Thermostat(double, double);				// метод расчета параметра лямбда для термостата Берендсена
	double lambda_Barostat(double, double);					// метод расчета параметра лямбда для баростата Берендсена
	double lambda_Proportional_Thermostat(double, double);	// метод расчета параметра лямбда для пропорционального термостата

	void load_LJ_forces();
	void load_forces();
};
#endif // !_MD_VV_H_