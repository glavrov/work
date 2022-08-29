#include <fstream>
#include <iomanip>
#include "Dump.h"

Dump::Dump()
{
	// очистка файлов
	fstream f_open;
	f_open.open("Data/borderline_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/electrons_square_CV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/velocity_n_CV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/e_calculation_CV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/tt_calculation_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/e_calculation_CV.xyz", ios::app);
	f_open << setw(15) << right << "Кинетическая:";
	f_open << setw(15) << right << "Потенциальная:";
	f_open << setw(13) << right << "Полная:" << endl << endl;
	f_open.close();
	// открытие файлов для вывода данных
	this->f_Data_borderline.open("Data/borderline_VV.xyz", ios::app);
	this->f_Data_rVector.open("Data/electrons_square_CV.xyz", ios::app);
	this->f_Data_vVector.open("Data/velocity_n_CV.xyz", ios::app);
	this->f_Data_e_calculation.open("Data/e_calculation_CV.xyz", ios::app);
	this->f_Data_t_calculation.open("Data/tt_calculation_VV.xyz", ios::app);

	this->particle_COUNT = 0;
	this->ParticleMassive = new Particle[this->particle_COUNT];
}

Dump::Dump(const Particle& _ParticleMassive, int _particle_COUNT)
{
	
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	*this->ParticleMassive = _ParticleMassive;
}

Dump::Dump(Particle _ParticleMassive[], int _particle_COUNT)				// разобрать момент с указателем
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->ParticleMassive = _ParticleMassive;
}

Dump::~Dump()
{
	this->f_Data_borderline.close();
	this->f_Data_rVector.close();
	this->f_Data_vVector.close();
	this->f_Data_e_calculation.close();
	this->f_Data_t_calculation.close();
}

void Dump::set_ParticleMassive(Particle& _ParticleMassive)
{
	*this->ParticleMassive = _ParticleMassive;
}

void Dump::set_ParticleMassive(Particle _ParticleMassive[])
{
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->ParticleMassive = _ParticleMassive;
}

void Dump::set_particle_COUNT(int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
}
// метод печати радиус-векторов в файл
void Dump::print_rVector()
{
	this->f_Data_rVector << this->particle_COUNT << endl << endl;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->f_Data_rVector << setprecision(6);
		this->f_Data_rVector << setiosflags(ios::showpoint);
		this->f_Data_rVector << left << "h" << setw(4) << this->ParticleMassive[i].get_Particle_ID();
		this->f_Data_rVector << setw(12) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0];
		this->f_Data_rVector << setw(12) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1];
		this->f_Data_rVector << setw(12) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2];
		this->f_Data_rVector << endl;
	}
}
// метод печати векторов скоростей в файл
void Dump::print_vVector()
{
	this->f_Data_vVector << this->particle_COUNT << endl << endl;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->f_Data_vVector << setprecision(6);
		this->f_Data_vVector << setiosflags(ios::showpoint);
		this->f_Data_vVector << left << "h" << setw(4) << this->ParticleMassive[i].get_Particle_ID();
		this->f_Data_vVector << setw(12) << right << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0];
		this->f_Data_vVector << setw(12) << right << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1];
		this->f_Data_vVector << setw(12) << right << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2];
		this->f_Data_vVector << endl;
	}
}
// метод печати кинетической, потенциальной и полной механической энегрий в файл
void Dump::print_e_Calculation(double _E_kinetic, double _E_potential, double _E_full)
{
	this->f_Data_e_calculation << setw(15) << right << _E_kinetic;
	this->f_Data_e_calculation << setw(15) << right << _E_potential;
	this->f_Data_e_calculation << setw(15) << right << _E_full << endl;
}
// метод печати температуры системы в файл
void Dump::print_t_Calculation(double _T)
{
	this->f_Data_t_calculation << setw(15) << right << _T << endl;
}
//
void Dump::print_Borderline(double _coord_MIN, double _coord_MAX)
{
	int p_Count = 40000;
	double x = 0., y = 0.;
	double delta = (_coord_MAX - _coord_MIN) / 10000.;
	
	this->f_Data_borderline << p_Count << endl << endl;
	for (int i = 0; i < p_Count / 2; i++)
	{
		this->f_Data_borderline << setprecision(6);
		this->f_Data_borderline << setiosflags(ios::showpoint);
		this->f_Data_borderline << left << "h" << setw(4) << i;
		this->f_Data_borderline << setw(20) << right << x;
		this->f_Data_borderline << setw(20) << right << y;
		this->f_Data_borderline << setw(20) << right << 0.;
		this->f_Data_borderline << endl;

		y += _coord_MAX;

		this->f_Data_borderline << setprecision(6);
		this->f_Data_borderline << setiosflags(ios::showpoint);
		this->f_Data_borderline << left << "h" << setw(4) << i;
		this->f_Data_borderline << setw(20) << right << x;
		this->f_Data_borderline << setw(20) << right << y;
		this->f_Data_borderline << setw(20) << right << 0.;
		this->f_Data_borderline << endl;

		x += delta;
		y = 0.;
	}
	x = 0.; y = 0.;
	for (int i = 0; i < p_Count / 2; i++)
	{
		this->f_Data_borderline << setprecision(6);
		this->f_Data_borderline << setiosflags(ios::showpoint);
		this->f_Data_borderline << left << "h" << setw(4) << i;
		this->f_Data_borderline << setw(20) << right << x;
		this->f_Data_borderline << setw(20) << right << y;
		this->f_Data_borderline << setw(20) << right << 0.;
		this->f_Data_borderline << endl;

		x += _coord_MAX;

		this->f_Data_borderline << setprecision(6);
		this->f_Data_borderline << setiosflags(ios::showpoint);
		this->f_Data_borderline << left << "h" << setw(4) << i;
		this->f_Data_borderline << setw(20) << right << x;
		this->f_Data_borderline << setw(20) << right << y;
		this->f_Data_borderline << setw(20) << right << 0.;
		this->f_Data_borderline << endl;

		x -= _coord_MAX;
		y += delta;
	}
}