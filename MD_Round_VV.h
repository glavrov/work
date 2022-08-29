//-----------------------------------------------------------------------------------------------------
// Моделирование движения частиц в ограниченной круглой (шарообразной) области,
// границы которой являются "прозрачными" для частиц, и действуют периодические граничные условия,
// либо границы являются "непрозрачными", и действует закон отражения.
// Моделирование осуществляется скоростным методом Верле.
//-----------------------------------------------------------------------------------------------------

#pragma once
#ifndef _MD_VV_ROUND_H_
#define _MD_VV_ROUND_H_
#include <fstream>
#include <iomanip>
#include "Particle_VV.h"
#include "Dump_VV.h"

class MD_Round_VV
{
private:
	Particle_VV* ParticleMassive;							// массив указателей на объект типа Particle
	Dump_VV* dump_FILE;										// указатель на объект типа Dump
	Vector3D<double>* SumV;									// массив указателей на объект типа Vector3D<double> (расчет средней координаты)
	Vector3D<double>* SumV2;								// массив указателей на объект типа Vector3D<double> (расчет среднеквадратичной координаты)
	Vector3D<double>* velocity_Coefficient;					// массив указателей на объект типа Vector3D<double> (расчет коэффициента нормализации)
	int particle_COUNT;										// количество частиц
	double system_Temperature = 10.;						// температура
	double R;												// радиус области
	double velocity_MIN;									// минимальная скорость
	double velocity_MAX;									// максимальная скорость
	double delta_T;											// шаг по времени
	double k_B = 1.;										// постоянная Больцмана
	double T0 = .01;
	double particle_RADIUS = 0.;							// .65;

	double E_kinetic, E_potential, E_full;					// кинетическая, потенциальная, полная энергии
public:
	MD_Round_VV();											// конструкторы
	MD_Round_VV(int);
	MD_Round_VV(int, double);
	MD_Round_VV(int, double, double);
	MD_Round_VV(int, double, double, double);
	~MD_Round_VV();											// деструктор

	void Initialization();									// метод инициализации системы частиц (массив), круглое пространство
	void Velocity_Normalization();							// метод нормализации скоростей системы частиц
	void v_Normalization();									// метод нормализации скоростей системы частиц
	void Move_VV();											// метод моделирования движения системы частиц
	void Evolution_VV(int);									// метод моделирования эволюции системы

	Vector3D<double> force(Vector3D<double>&);				// метод расчета сил, действующих на частицу
	double LJ_potential(Vector3D<double>&);					// метод расчета потенциала Леннард-Джонса
	void Next_Radius_Vector();								// метод расчета радиус-векторов системы частиц в момент t+Δt
	void Next_Velocity_Vector_05t();						// метод расчета скоростей системы частиц в момент t+Δt/2
	void Next_Velocity_Vector();							// метод расчета скоростей системы частиц в момент t+Δt
	void Compute_Acceleration();							// метод расчета ускорений системы частиц
	void Deceleration();									// метод торможения системы частиц

	double Average_velocity(int);							// метод расчета средней скорости системы частиц
	double Average_velocity2(int);							// метод расчета среднеквадратической скорости системы частиц

	void E_calculation_VV();								// метод расчета кинетической, потенциальной и полной энергий системы частиц
	double lambda_Thermostat(double);						// метод расчета параметра лямбда для термостата Берендсена
	double lambda_Barostat(double, double);					// метод расчета параметра лямбда для баростата Берендсена
	double lambda_Proportional_Thermostat(double);			// метод расчета параметра лямбда для пропорционального термостата
};
#endif // !_MD_VV_ROUND_H_