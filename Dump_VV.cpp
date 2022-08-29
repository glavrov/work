#include <fstream>
#include <iomanip>
#include "Dump_VV.h"
#include "Vector3D.h"

Dump_VV::Dump_VV()
{
	// очистка файлов
	fstream f_open;
	f_open.open("Data/borderline_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/electrons_square_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/velocity_n_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/p_calculation_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/tt_calculation_VV.xyz", ios::out);
	f_open.close();
	f_open.open("Data/e_calculation_VV.xyz", ios::out);
	f_open << setw(15) << right << "Кинетическая:";
	f_open << setw(15) << right << "Потенциальная:";
	f_open << setw(15) << right << "Полная:" << endl << endl;
	f_open.close();

	// временный вывод в файл для поиска ошибок по ускорению
	f_open.open("Data/dx.xyz", ios::out);
	f_open << setprecision(6);
	f_open << setiosflags(ios::showpoint);
	f_open << left << setw(15) << "Pid";
	f_open << setw(15) << left << "dx";
	f_open << setw(15) << left << "dy";
	f_open << setw(15) << left << "dz";
	f_open << setw(15) << left << "forceX";
	f_open << setw(15) << left << "forceY";
	f_open << setw(15) << left << "Ax";
	f_open << setw(15) << left << "Ay";
	f_open << endl << endl;
	f_open.close();

	// открытие файлов для вывода данных
	this->f_Data_borderline.open("Data/borderline_VV.xyz", ios::app);
	this->f_Data_rVector.open("Data/electrons_square_VV.xyz", ios::app);
	this->f_Data_vVector.open("Data/velocity_n_VV.xyz", ios::app);
	this->f_Data_e_calculation.open("Data/e_calculation_VV.xyz", ios::app);
	this->f_Data_t_calculation.open("Data/tt_calculation_VV.xyz", ios::app);
	this->f_Data_p_calculation.open("Data/p_calculation_VV.xyz", ios::app);
	
	// временный вывод в файл для поиска ошибок по ускорению
	this->f_dx.open("Data/dx.xyz", ios::app);

	this->particle_COUNT = 0;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
}

Dump_VV::Dump_VV(const Particle_VV& _ParticleMassive, int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	*this->ParticleMassive = _ParticleMassive;
}

Dump_VV::Dump_VV(Particle_VV _ParticleMassive[], int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->ParticleMassive = _ParticleMassive;
}

Dump_VV::~Dump_VV()
{
	this->f_Data_rVector.close();
	this->f_Data_vVector.close();
	this->f_Data_e_calculation.close();
	this->f_Data_p_calculation.close();
	this->f_Data_t_calculation.close();
	this->f_Data_borderline.close();
}

void Dump_VV::set_ParticleMassive(Particle_VV& _ParticleMassive)
{
	this->ParticleMassive = nullptr;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	*this->ParticleMassive = _ParticleMassive;
}

void Dump_VV::set_ParticleMassive(Particle_VV _ParticleMassive[])
{
	this->ParticleMassive = nullptr;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->ParticleMassive = _ParticleMassive;
	//this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	//this->ParticleMassive = _ParticleMassive;
}

void Dump_VV::set_particle_COUNT(int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
}
// метод печати радиус-векторов в файл
void Dump_VV::print_rVector()
{
	this->f_Data_rVector << this->particle_COUNT << endl << endl;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->f_Data_rVector << setprecision(6);
		this->f_Data_rVector << setiosflags(ios::showpoint);
		this->f_Data_rVector << left  << "h" << setw(5) << this->ParticleMassive[i].get_Particle_ID();
		this->f_Data_rVector << setw(15) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0];
		this->f_Data_rVector << setw(15) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1];
		this->f_Data_rVector << setw(15) << right << this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2];
		this->f_Data_rVector << endl;
	}
}
// метод печати векторов скоростей в файл
void Dump_VV::print_vVector()
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
void Dump_VV::print_e_Calculation(double _E_kinetic, double _E_potential, double _E_full)
{
	this->f_Data_e_calculation << setw(15) << right << _E_kinetic;
	this->f_Data_e_calculation << setw(15) << right << _E_potential;
	this->f_Data_e_calculation << setw(15) << right << _E_full << endl;
}
// метод печати давления системы в файл
void Dump_VV::print_p_Calculation(double _P)
{
	this->f_Data_p_calculation << setw(15) << right << _P << endl;
}
// метод печати температуры системы в файл
void Dump_VV::print_t_Calculation(double _T)
{
	this->f_Data_t_calculation << setw(15) << right << _T << endl;
}
// методо печати границ области в файл
void Dump_VV::print_Borderline(double _R)
{
	int p_Count = 10000;
	double x = 0., y = 0.;
	double theta = 0.;
	
	this->f_Data_borderline << p_Count << endl << endl;
	for (int i = 0; i < p_Count; i++)
	{
		theta = this->random(0, 2. * M_PI);
		x = _R * cos(theta);
		y = _R * sin(theta);
		this->f_Data_borderline << setprecision(6);
		this->f_Data_borderline << setiosflags(ios::showpoint);
		this->f_Data_borderline << left << "h" << setw(4) << i;
		this->f_Data_borderline << setw(20) << right << x;
		this->f_Data_borderline << setw(20) << right << y;
		this->f_Data_borderline << setw(20) << right << 0.;
		this->f_Data_borderline << endl;
	}
}
// метод генерации равномерно распределенных случайных чисел типа <double> для получения полярного угла
double Dump_VV::random(double _min, double _max)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> urd(_min, _max);

	double random = urd(gen);
	return random;
}
//
void Dump_VV::print_dx(double _dx, double _dy, double _dz, double _forceX, double _forceY)
{
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->f_dx << setprecision(6);
		this->f_dx << setiosflags(ios::showpoint);
		this->f_dx << left <<  "h" << setw(14) << this->ParticleMassive[i].get_Particle_ID();
		this->f_dx << setw(15) << left << _dx;
		this->f_dx << setw(15) << left << _dy;
		this->f_dx << setw(15) << left << _dz;
		this->f_dx << setw(15) << left << _forceX;
		this->f_dx << setw(15) << left << _forceY;
		this->f_dx << setw(15) << left << this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0];
		this->f_dx << setw(15) << left << this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1];
		this->f_dx << endl;
	}
	this->f_dx << endl;
}
//
void Dump_VV::print_Bline(double _coord_MIN, double _coord_MAX)
{
	int p_Count = 40000;
	double x = 0., y = 0.;
	double delta = (_coord_MAX - _coord_MIN) / 10000.;

	this->f_Data_borderline << p_Count << endl << endl;
	for (int i = 0; i < p_Count / 2; i++)
	{
		if (i < 10000)
		{
			this->f_Data_borderline << setprecision(5);
			this->f_Data_borderline << setiosflags(ios::showpoint);
			this->f_Data_borderline << left << "h" << setw(4) << i;
			this->f_Data_borderline << setw(15) << right << x;
			this->f_Data_borderline << setw(15) << right << y;
			this->f_Data_borderline << setw(15) << right << 0.;
			this->f_Data_borderline << endl;
			x += delta;
		}
		if (i == 10000)
		{
			x = 0.;
			y += _coord_MAX;
		}
		if (i >= 10000)
		{
			this->f_Data_borderline << setprecision(5);
			this->f_Data_borderline << setiosflags(ios::showpoint);
			this->f_Data_borderline << left << "h" << setw(4) << i;
			this->f_Data_borderline << setw(15) << right << x;
			this->f_Data_borderline << setw(15) << right << y;
			this->f_Data_borderline << setw(15) << right << 0.;
			this->f_Data_borderline << endl;
			x += delta;
		}
	}
	x = 0.; y = 0.;
	for (int i = 20000; i < p_Count; i++)
	{
		if (i < 30000)
		{
			this->f_Data_borderline << setprecision(5);
			this->f_Data_borderline << setiosflags(ios::showpoint);
			this->f_Data_borderline << left << "h" << setw(4) << i;
			this->f_Data_borderline << setw(15) << right << x;
			this->f_Data_borderline << setw(15) << right << y;
			this->f_Data_borderline << setw(15) << right << 0.;
			this->f_Data_borderline << endl;
			y += delta;
		}
		if (i == 30000)
		{
			x += _coord_MAX;
			y = 0.;
		}
		if (i >= 30000)
		{
			this->f_Data_borderline << setprecision(5);
			this->f_Data_borderline << setiosflags(ios::showpoint);
			this->f_Data_borderline << left << "h" << setw(4) << i;
			this->f_Data_borderline << setw(15) << right << x;
			this->f_Data_borderline << setw(15) << right << y;
			this->f_Data_borderline << setw(15) << right << 0.;
			this->f_Data_borderline << endl;
			y += delta;
		}
	}
}