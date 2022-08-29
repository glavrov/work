#include "MD_VV.h"
#include <omp.h>
// контрукторы
MD_VV::MD_VV()
{
	this->particle_COUNT = 20;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_VV::MD_VV(int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_VV::MD_VV(int _particle_COUNT, double _delta_T)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->delta_T = _delta_T;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_VV::MD_VV(int _particle_COUNT, double _coord_MIN, double _coord_MAX, double _velocity_MIN, double _velocity_MAX)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->coord_MIN = _coord_MIN;
	this->coord_MAX = _coord_MAX;
	this->velocity_MIN = _velocity_MIN;
	this->velocity_MAX = _velocity_MAX;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// деструктор
MD_VV::~MD_VV()
{}
// метод инициализации массива частиц, сгенерированных случайным образом (равномерное распределение)
void MD_VV::Initialization()
{
	//this->load_LJ_forces();
	cout << "инициализация" << endl;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_Particle_ID(i);
		this->ParticleMassive[i].set_Radius_Vector(this->coord_MIN, this->coord_MAX);
		this->ParticleMassive[i].set_Velocity_Vector(this->velocity_MIN, this->velocity_MAX);
	}
	// расчет векторов средних скоростей и среднеквадратичных скоростей
	this->SumV->set_Vector(this->Average_velocity(0), this->Average_velocity(1), this->Average_velocity(2));
	this->SumV2->set_Vector(this->Average_velocity2(0), this->Average_velocity2(1), this->Average_velocity2(2));
	// вывод радиус-веторов и векторов скоростей частиц после генерации до вызова метода "нормализации"
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
	this->dump_FILE->print_vVector();
	this->dump_FILE->print_Bline(this->coord_MIN, this->coord_MAX);
	for (int i = 0; i < particle_COUNT; i++)
	{
		cout << i << ": " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] << "; " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] << endl;
	}
}
// метод инициализации массива частиц, сгенерированных случайным образом (равномерное распределение), внутри окружности радиуса R
void MD_VV::Initialization(double _R)
{
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_Particle_ID(i);
		this->ParticleMassive[i].set_Radius_Vector_Circle(_R, _R, _R);
		this->ParticleMassive[i].set_Velocity_Vector(this->velocity_MIN, this->velocity_MAX);
	}
	// расчет векторов средних скоростей и среднеквадратичных скоростей
	this->SumV->set_Vector(this->Average_velocity(0), this->Average_velocity(1), this->Average_velocity(2));
	this->SumV2->set_Vector(this->Average_velocity2(0), this->Average_velocity2(1), this->Average_velocity2(2));
	// вывод радиус-веторов и векторов скоростей частиц после генерации до вызова метода "нормализации"
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
	this->dump_FILE->print_vVector();
	for (int i = 0; i < particle_COUNT; i++)
		cout << i << ": " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] << "; " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] << endl;	
}
// метод нормализации векторов скоростей частиц
void MD_VV::Velocity_Normalization()
{
	cout << "нормализация" << endl;
	// коэффициент нормализации
	this->velocity_Coefficient->set_Vector(sqrt((2. * this->system_Temperature) / this->SumV2->get_Vector()[0]),
										   sqrt((2. * this->system_Temperature) / this->SumV2->get_Vector()[1]),
										   0.); // sqrt((2. * this->particle_Temperature) / this->SumV2->get_Vector()[2]));
	// координаты после нормализации
	double velocity_X_NEW = 0., velocity_Y_NEW = 0., velocity_Z_NEW = 0.;
	// нормализация векторов скоростей
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		velocity_X_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] - this->Average_velocity(0)) * this->velocity_Coefficient->get_Vector()[0];
		velocity_Y_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] - this->Average_velocity(1)) * this->velocity_Coefficient->get_Vector()[1];
		//velocity_Z_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] - this->Average_velocity(2)) * this->velocity_Coefficient->get_Vector()[2];
		this->ParticleMassive[i].set_Velocity_Vector(velocity_X_NEW, velocity_Y_NEW, velocity_Z_NEW);
	}
	// печать векторов скоростей в файл
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_vVector();
	for (int i = 0; i < particle_COUNT; i++)
	{
		cout << i << ": " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] << "; " << this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] << endl;
	}
	// расчет кинетической, потенциальной и полной энергий в момент времени t
	this->dump_FILE->print_t_Calculation(this->system_Temperature);
	this->E_calculation_VV();
}
// метод для моделирования эволюции молекулярной динамики на основе скоростного метода Верле
void MD_VV::Evolution_VV(int _step_COUNT, int _space_TYPE, int _calculation_METHOD)
{
	this->space_TYPE = _space_TYPE;
	cout << "[0%";
	for (int i = 0; i < _step_COUNT; i++)
	{
		//if ((i + 1) % 100 == 0)
		//	this->v_Normalization();
		this->Move_VV(_space_TYPE, _calculation_METHOD);
		if (i % 1000 == 0)
			cout << ".";
	}
	cout << "100%]";
}
// метод движения частиц - реализация скоростного метода Верле
void MD_VV::Move_VV(int _space_TYPE, int _calculation_METHOD)
{
	if (_space_TYPE == 1)
	{
		if (_calculation_METHOD == 1)
		{
			// расчет скоростей частиц в момент времени t+Δt/2
			this->Next_Velocity_Vector_05t();
			// расчет радиус-векторов в момент времени t+Δt
			this->Next_Radius_Vector(_space_TYPE, _calculation_METHOD);
			// расчет ускорений в момент времени t+Δt
			this->Compute_Acceleration(); //this->load_forces(); 
			// расчет скоростей частиц в момент времени t+Δt
			this->Next_Velocity_Vector();
			// расчет кинетической, потенциальной и полной энергий в момент времени t+Δt
			//this->E_calculation_VV();
			// печать в файл
			//this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
			this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
			this->dump_FILE->print_rVector();

			for (int i = 0; i < this->particle_COUNT; i++)
				this->ParticleMassive[i].set_Flag(0);
		}
		else if (_calculation_METHOD == 2)
		{
			// расчет радиус-векторов в момент времени t+Δt
			this->Next_Radius_Vector(_space_TYPE, _calculation_METHOD);
			// расчет скоростей частиц в момент времени t+Δt/2
			this->Next_Velocity_Vector_05t();
			// расчет ускорений в момент времени t+Δt
			this->Compute_Acceleration(); //this->load_forces(); 
			// расчет скоростей частиц в момент времени t+Δt
			this->Next_Velocity_Vector();
			// расчет кинетической, потенциальной и полной энергий в момент времени t+Δt
			//this->E_calculation_VV();
			// печать в файл
			//this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
			this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
			this->dump_FILE->print_rVector();

			for (int i = 0; i < this->particle_COUNT; i++)
				this->ParticleMassive[i].set_Flag(0);
		}
	}
	else if (_space_TYPE == 2)
	{
		if (_calculation_METHOD == 1)
		{
			// расчет скоростей частиц в момент времени t+Δt/2
			this->Next_Velocity_Vector_05t();
			// расчет радиус-векторов в момент времени t+Δt
			this->Next_Radius_Vector();
			//this->Next_Radius_Vector(_space_TYPE, _calculation_METHOD);
			// расчет ускорений в момент времени t+Δt
			this->Compute_Acceleration();
			// расчет скоростей частиц в момент времени t+Δt
			this->Next_Velocity_Vector();
			// расчет кинетической, потенциальной и полной энергий в момент времени t+Δt
			//this->E_calculation_VV();
			// печать в файл
			//this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
			this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
			this->dump_FILE->print_rVector();

			for (int i = 0; i < this->particle_COUNT; i++)
				this->ParticleMassive[i].set_Flag(0);
		}
		else if (_calculation_METHOD == 2)
		{
			// расчет радиус-векторов в момент времени t+Δt
			this->Next_Radius_Vector(_space_TYPE, _calculation_METHOD);
			// расчет скоростей частиц в момент времени t+Δt/2
			this->Next_Velocity_Vector_05t();
			// расчет ускорений в момент времени t+Δt
			this->Compute_Acceleration();
			// расчет скоростей частиц в момент времени t+Δt
			this->Next_Velocity_Vector();
			// расчет кинетической, потенциальной и полной энергий в момент времени t+Δt
			//this->E_calculation_VV();
			// печать в файл
			//this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
			this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
			this->dump_FILE->print_rVector();

			for (int i = 0; i < this->particle_COUNT; i++)
				this->ParticleMassive[i].set_Flag(0);
		}
	}
}
// расчет радиус-векторов в момент времени t+Δt
void MD_VV::Next_Radius_Vector(int _space_TYPE, int _calculation_METHOD)
{
	double x_NEW = 0., y_NEW = 0., z_NEW = 0.;		// новые координаты радиус-вектора для периодических граничный условий
	double lambda = 0.;

	if (_space_TYPE == 1)
	{// незамкнутое пространство (с периодическими граничными условиями)
		if (_calculation_METHOD == 1) // R(t+Δt) = R(t) + v(t+Δt/2)*Δt, незамкнутое пространство (с периодическими граничными условиями)
		{
			for (int i = 0; i < this->particle_COUNT; i++)
			{
				x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0];
				y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1];
				//z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2];
				if (x_NEW > this->coord_MAX)
					x_NEW -= (this->coord_MAX - this->coord_MIN);
				else if (x_NEW < this->coord_MIN)
					x_NEW += (this->coord_MAX - this->coord_MIN);
				if (y_NEW > this->coord_MAX)
					y_NEW -= (this->coord_MAX - this->coord_MIN);
				else if (y_NEW < this->coord_MIN)
					y_NEW += (this->coord_MAX - this->coord_MIN);
				/*if (z_NEW > this->coord_MAX)
					z_NEW -= (this->coord_MAX - this->coord_MIN);
				else if (z_NEW < this->coord_MIN)
					z_NEW += (this->coord_MAX - this->coord_MIN);*/
				this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
			}
		}
		else if (_calculation_METHOD == 2) // R(t+Δt) = R(t) + v(t)*Δt + aΔt*Δt/2, незамкнутое пространство (с периодическими граничными условиями)
		{
			for (int i = 0; i < this->particle_COUNT; i++)
			{
				x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T * this->delta_T / 2.;
				y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T * this->delta_T / 2.;
				//z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T * this->delta_T / 2.;
				if (x_NEW > this->coord_MAX)
					x_NEW -= (this->coord_MAX - this->coord_MIN);
				else if (x_NEW < this->coord_MIN)
					x_NEW += (this->coord_MAX - this->coord_MIN);
				if (y_NEW > this->coord_MAX)
					y_NEW -= (this->coord_MAX - this->coord_MIN);
				else if (y_NEW < this->coord_MIN)
					y_NEW += (this->coord_MAX - this->coord_MIN);
				/*if (z_NEW > this->coord_MAX)
					z_NEW -= this->coord_MAX;
				else if (z_NEW < this->coord_MIN)
					z_NEW += this->coord_MAX;*/
				this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
			}
		}
	}
	else if (_space_TYPE == 2)
	{// замкнутое пространство
		if (_calculation_METHOD == 1) // R(t+Δt) = R(t) + v(t+Δt/2)*Δt, замкнутое пространство
		{
			for (int i = 0; i < this->particle_COUNT; i++)
			{
				x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0];
				y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1];
				//z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2];
				if (x_NEW >= this->coord_MAX)
				{
					x_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (x_NEW <= this->coord_MIN)
				{
					x_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				if (y_NEW >= this->coord_MAX)
				{
					y_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (y_NEW <= this->coord_MIN)
				{
					y_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				if (z_NEW >= this->coord_MAX)
				{
					z_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (z_NEW <= this->coord_MIN)
				{
					z_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
			}
			// лямбда(t)
			//lambda = this->lambda_Barostat(.001, 1.2);
			//for (int i = 0; i < this->particle_COUNT; i++)
			//{
			//	x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] * lambda;
			//	y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] * lambda;
			//	z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] * lambda;
			//	this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
			//}
		}
		else if (_calculation_METHOD == 2) // R(t+Δt) = R(t) + v(t)*Δt + aΔt*Δt/2, замкнутое пространство
		{
			for (int i = 0; i < this->particle_COUNT; i++)
			{
				x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T * this->delta_T / 2.;
				y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T * this->delta_T / 2.;
				//z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] + this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] * this->delta_T + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T * this->delta_T / 2.;
				if (x_NEW >= this->coord_MAX)
				{
					x_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (x_NEW <= this->coord_MIN)
				{
					x_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				if (y_NEW >= this->coord_MAX)
				{
					y_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (y_NEW <= this->coord_MIN)
				{
					y_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				if (z_NEW >= this->coord_MAX)
				{
					z_NEW = this->coord_MAX;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				else if (z_NEW <= this->coord_MIN)
				{
					z_NEW = this->coord_MIN;
					this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
					this->ParticleMassive[i].set_Velocity_Vector_05t(this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0],
						this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1],
						-this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2]);
				}
				this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
			}
		}
	}
}
// расчет радиус-векторов в момент времени t+Δt
void MD_VV::Next_Radius_Vector()
{
	double x_NEW = 0., y_NEW = 0., z_NEW = 0.;		// новые координаты радиус-вектора для периодических граничный условий
	double lambda = 0.;
	double k = 0., b = 0.;
	double x = 0., y = 0.;
	double Vx = 0., Vy = 0., Vz = 0.;
	double x_CROSS = 0., y_CROSS = 0.;
	double alpha = 0.;
	double a = 0., c = 0., bb = 0.;

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		x = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0];
		y = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1];

		x_NEW = x + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0];
		y_NEW = y + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1];

		Vx = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0];
		Vy = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1];

		if (x_NEW > this->coord_MAX)
		{
			this->ParticleMassive[i].set_Flag(1);
			k = (y - y_NEW) / (x - x_NEW);
			b = (x * y_NEW - x_NEW * y) / (x - x_NEW);

			x_CROSS = this->coord_MAX;
			y_CROSS = k * x_CROSS + b;

			a = x_NEW - this->coord_MAX;
			
			x_NEW -= 2. * a;
			
			this->ParticleMassive[i].set_Velocity_Vector(-Vx, Vy, Vz);
		}
		else if (x_NEW < this->coord_MIN)
		{
			this->ParticleMassive[i].set_Flag(1);
			k = (y - y_NEW) / (x - x_NEW);
			b = (x * y_NEW - x_NEW * y) / (x - x_NEW);

			x_CROSS = this->coord_MIN;
			y_CROSS = k * x_CROSS + b;

			a = x_NEW - this->coord_MIN;

			x_NEW -= 2. * a;

			this->ParticleMassive[i].set_Velocity_Vector(-Vx, Vy, Vz);
		}

		if (y_NEW > this->coord_MAX)
		{
			this->ParticleMassive[i].set_Flag(1);
			k = (y - y_NEW) / (x - x_NEW);
			b = (x * y_NEW - x_NEW * y) / (x - x_NEW);

			x_CROSS = this->coord_MIN;
			y_CROSS = k * x_CROSS + b;

			a = y_NEW - this->coord_MAX;

			y_NEW -= 2. * a;

			this->ParticleMassive[i].set_Velocity_Vector(Vx, -Vy, Vz);
		}
		else if (y_NEW < this->coord_MIN)
		{
			this->ParticleMassive[i].set_Flag(1);
			k = (y - y_NEW) / (x - x_NEW);
			b = (x * y_NEW - x_NEW * y) / (x - x_NEW);

			y_CROSS = this->coord_MIN;
			x_CROSS = (y_CROSS - b) / k;

			a = y_NEW - this->coord_MIN;

			y_NEW -= 2. * a;

			this->ParticleMassive[i].set_Velocity_Vector(Vx, -Vy, Vz);
		}
		this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
	}
}
// расчет скоростей частиц в момент времени t+Δt/2
void MD_VV::Next_Velocity_Vector_05t()
{
	double V_05X = 0., V_05Y = 0., V_05Z = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		V_05X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T / 2.;
		V_05Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T / 2.;
		//V_05Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T / 2.;
		
		this->ParticleMassive[i].set_Velocity_Vector_05t(V_05X, V_05Y, V_05Z);
	}
}
// расчет скоростей частиц в момент времени t+Δt
void MD_VV::Next_Velocity_Vector()
{
	double V_X = 0., V_Y = 0., V_Z = 0.;
	double lambda = 1.;

	//lambda = this->lambda_Thermostat(.1, .01);
	lambda = this->lambda_Proportional_Thermostat(1., 2.);
	 
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		if (this->ParticleMassive[i].get_Flag() == 1)
			continue;
		else
		{
			V_X = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[0] * lambda + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T / 2.;
			V_Y = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[1] * lambda + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T / 2.;
			//V_Z = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector()[2] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T / 2.;
			this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
		}
	}
	// лямбда(t)
	//lambda = this->lambda_Proportional_Thermostat(.1, 2.);
	//lambda = this->lambda_Thermostat(.1, .01);
	//for (int i = 0; i < this->particle_COUNT; i++)
	//{
	//	V_X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] * lambda;
	//	V_Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] * lambda;
	//	//V_Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] * lambda;
	//	this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
	//}
}
// метод нормализации скоростей частиц
void MD_VV::v_Normalization()
{
	// координаты после нормализации
	double velocity_X_NEW = 0., velocity_Y_NEW = 0., velocity_Z_NEW = 0.;
	double averageX = this->Average_velocity(0);
	double avarageY = this->Average_velocity(1);
	// нормализация векторов скоростей
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		velocity_X_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] - averageX;
		velocity_Y_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] - avarageY;
		//velocity_Z_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] - this->Average_velocity(2)) * this->velocity_Coefficient->get_Vector()[2];
		this->ParticleMassive[i].set_Velocity_Vector(velocity_X_NEW, velocity_Y_NEW, velocity_Z_NEW);
	}
}
// расчет средней скорости
double MD_VV::Average_velocity(int coordinate_Number)
{
	double average_Velocity = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		average_Velocity += this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[coordinate_Number];
	average_Velocity /= this->particle_COUNT;
	return average_Velocity;
}
// расчет среднеквадратической скорости
double MD_VV::Average_velocity2(int coordiname_Number)
{
	double average_Velocity2 = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		average_Velocity2 += pow(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[coordiname_Number], 2.);
	average_Velocity2 /= this->particle_COUNT;
	return average_Velocity2;
}
// расчет ускорений для метода Верле
void MD_VV::Compute_Acceleration()
{
	Vector3D<double>* dR = new Vector3D<double>(3);			// расстояние между i-ой и j-ой частицами
	Vector3D<double>* Rcut = new Vector3D<double>(3);		// радиус "обрезания"
	double aXi = 0., aYi = 0., aZi = 0.;					// ускорение i-ой частицы
	double aXj = 0., aYj = 0., aZj = 0.;					// ускорение j-ой частицы
	double particle_mass_i = 0., particle_mass_j = 0.;		// масса частиц

	double dx = 0., dy = 0., dz = 0.;						// компоненты dR
	double r = 0;											// |dR|
	double k = 0., b = 0.;									// коэффициенты уравнения прямой
	double x0 = 0., y0 = 0.;								// координаты i-ой частицы
	double xj = 0., yj = 0.;								// координаты j-ой частицы
	double x1 = 0., x2 = 0.;								// корни квадратного уравнения
	double y1 = 0., y2 = 0.;								// значения функции в точках пересечения прямой с окружностью

	double X_razn = 0., Y_razn = 0.;
	double x_razn1 = 0., y_razn1 = 0.;
	double x_razn2 = 0., y_razn2 = 0.;

	for (int i = 0; i < this->particle_COUNT; i++)
		this->ParticleMassive[i].set_Acceleration_Vector(0., 0., 0.);

	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		particle_mass_i = this->ParticleMassive[i].get_Particle_Mass();

		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			particle_mass_j = this->ParticleMassive[j].get_Particle_Mass();

			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];

			r = sqrt(dx * dx + dy * dy);// +dz * dz);

			if (r < this->particle_RADIUS)
			{
				x0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0];
				y0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1];
				xj = this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
				yj = this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];

				k = dy / dx;
				b = (x0 * yj - xj * y0) / dx;
				
				double D = 4. * pow(k * b - x0 - y0 * k, 2.) - 4 * (k * k + 1) * (x0 * x0 + b * b - 2. * y0 * b + y0 * y0 - this->particle_RADIUS * this->particle_RADIUS); // дискриминант кв. ур-ния (x-x0)^2 + (y-y0)^2 = pR^2
				if (D < 0.) { cout << "D < 0"; }
				else
				{
					x1 = (-k * b + x0 + k * y0 - sqrt(D)) / (k * k + 1.);
					x2 = (-k * b + x0 + k * y0 + sqrt(D)) / (k * k + 1.);
				}

				y1 = k * x1 + b;
				y2 = k * x2 + b;

				x_razn1 = x1 - xj;
				y_razn1 = y1 - yj;
				x_razn2 = x2 - xj;
				y_razn2 = y2 - yj;
				
				double u1 = 0., u2 = 0.;
				if (sqrt(x_razn1 * x_razn1 + y_razn1 * y_razn1) <= sqrt(x_razn2 * x_razn2 + y_razn2 * y_razn2))
				{
					u1 = x1; u2 = y1;
					if (this->space_TYPE == 1)
					{
						if (u1 < this->coord_MIN)
							u1 += this->coord_MAX;
						else if (u1 > this->coord_MAX)
							u1 -= this->coord_MAX;

						if (u2 < this->coord_MIN)
							u2 += this->coord_MAX;
						else if (u2 > this->coord_MAX)
							u2 -= this->coord_MAX;
					}
					
					this->ParticleMassive[j].set_Radius_Vector(u1, u2, 0.);

					dx = x0 - u1;
					dy = y0 - u2;
					/*if (x0 > xj)
						dx = x0 - u1;
					else
						dx = u1 - x0;

					if (y0 > yj)
						dy = y0 - u2;
					else
						dy = u2 - y0;*/
				}
				else
				{
					u1 = x2; u2 = y2;

					if (this->space_TYPE == 1)
					{
						if (u1 < this->coord_MIN)
							u1 += this->coord_MAX;
						else if (u1 > this->coord_MAX)
							u1 -= this->coord_MAX;

						if (u2 < this->coord_MIN)
							u2 += this->coord_MAX;
						else if (u2 > this->coord_MAX)
							u2 -= this->coord_MAX;
					}

					this->ParticleMassive[j].set_Radius_Vector(u1, u2, 0.);

					dx = x0 - u1;
					dy = y0 - u2;
					/*if (x0 > xj)
						dx = x0 - u1;
					else
						dx = u1 - x0;

					if (y0 > yj)
						dy = y0 - u2;
					else
						dy = u2 - y0;*/
				}
			}

			dR->set_Vector(dx, dy, dz);

			this->ParticleMassive[i].set_Acceleration_Vector(this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] + this->force(*dR).get_Vector()[0] / particle_mass_i,	// aXi
															 this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] + this->force(*dR).get_Vector()[1] / particle_mass_i,	// aYi
															 aZi);	// aZi //this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] + this->force(*dR).get_Vector()[2] / this->ParticleMassive[i].get_Particle_Mass());
			this->ParticleMassive[j].set_Acceleration_Vector(this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[0] - this->force(*dR).get_Vector()[0] / particle_mass_j,	// aXj
															 this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[1] - this->force(*dR).get_Vector()[1] / particle_mass_j,	// aYj
															 aZj);	// aZj //this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[2] - this->force(*dR).get_Vector()[2] / this->ParticleMassive[j].get_Particle_Mass());

			this->E_potential += this->LJ_potential(*dR);
		}
	}

	for (int i = 0; i < this->particle_COUNT; i++)
		this->E_kinetic += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	this->E_kinetic /= 2.;
	// полная энергия
	this->E_full = this->E_kinetic + this->E_potential;
	// вывод расчитанных энергий в файл
	this->dump_FILE->print_e_Calculation(this->E_kinetic, this->E_potential, this->E_full);
	// обнуление полей, отвечающих за значения энергий
	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// расчет кинетической, потенциальной и полной энергий
void MD_VV::E_calculation_VV()
{
	Vector3D<double>* dr = new Vector3D<double>(3);
	double dx = 0., dy = 0., dz = 0.;
	// кинетическая энергия
	for (int i = 0; i < this->particle_COUNT; i++)
		this->E_kinetic += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	this->E_kinetic /= 2.;
	// потенциальная энергия
	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		for (int j = i + 1 ; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];

			double r = sqrt(dx * dx + dy * dy + dz * dz);

			if (r < this->particle_RADIUS)
			{
				double x0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0];
				double y0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1];
				double xj = this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
				double yj = this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];

				double k = dy / dx;
				double b = (x0 * yj - xj * y0) / dx;

				double D = 4. * pow(k * b - x0 - y0 * k, 2.) - 4 * (k * k + 1) * (x0 * x0 + b * b - 2. * y0 * b + y0 * y0 - this->particle_RADIUS * this->particle_RADIUS);
				double x1 = 0., x2 = 0.;
				if (D < 0.) { cout << "D < 0"; }
				else
				{
					x1 = (-k * b + x0 + k * y0 - sqrt(D)) / (k * k + 1.);
					x2 = (-k * b + x0 + k * y0 + sqrt(D)) / (k * k + 1.);
				}

				double y1 = k * x1 + b;
				double y2 = k * x2 + b;

				double x_razn1 = x1 - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
				double y_razn1 = y1 - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
				double x_razn2 = x2 - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
				double y_razn2 = y2 - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];

				double X_razn = 0., Y_razn = 0.;
				double u1 = 0., u2 = 0.;
				
				if (sqrt(x_razn1 * x_razn1 + y_razn1 * y_razn1) <= sqrt(x_razn2 * x_razn2 + y_razn2 * y_razn2))
				{
					u1 = x1; u2 = y1;

					dx = x0 - u1;
					dy = y0 - u2;
				}
				else
				{
					u1 = x2; u2 = y2;

					dx = x0 - u1;
					dy = y0 - u2;
				}
			}

			dr->set_Vector(dx, dy, dz);

			this->E_potential += this->LJ_potential(*dr);
		}
	}
	// полная энергия
	this->E_full = this->E_kinetic + this->E_potential;
	// вывод расчитанных энергий в файл
	this->dump_FILE->print_e_Calculation(this->E_kinetic, this->E_potential, this->E_full);
	// обнуление полей, отвечающих за значения энергий
	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// расчет силы, действующей на частицы, для метода Верле
Vector3D<double> MD_VV::force(Vector3D<double>& _dR)
{
	double rl = _dR.module();
	double r6 = pow(rl, 6.);
	double fr = ((48. / r6) * (1. / r6 - .5)) / pow(rl, 2.);
	Vector3D<double>* force = new Vector3D<double>(3);
	force->set_Vector(fr * _dR.get_Vector()[0], fr * _dR.get_Vector()[1], 0.); // fr* _dR.get_Vector()[2]);
	return *force;
}
// вычисление потенциала Леннард-Джонса
double MD_VV::LJ_potential(Vector3D<double>& _dr)
{
	double rm = _dr.module();
	double rm6 = pow(rm, 6);
	double ur = (4. * ((1. / rm6) - 1.)) / rm6;
	return ur;
}
// метод расчета параметра лямбда для термостата Берендсена
double MD_VV::lambda_Thermostat(double _T0, double _tau)
{
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	double T = E / (2. * this->k_B * this->particle_COUNT);
	double lambda = sqrt(1. + this->delta_T / _tau * ((_T0 / T) - 1.));
	this->dump_FILE->print_t_Calculation(T);
	return lambda;
}
// метод расчета параметра лямбда для баростата Берендсена
double MD_VV::lambda_Barostat(double _P0, double _tau)
{
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Radius_Vector().module(), 2.);
	double E2 = E / this->particle_COUNT;

	Vector3D<double>* dR = new Vector3D<double>(3);
	double dx = 0., dy = 0., dz = 0.;

	double SumRF = 0.;
	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];

			dR->set_Vector(dx, dy, dz);

			SumRF += dR->module() * this->ParticleMassive[i].get_Acceleration_Vector().module(); // m=1
		}
	}
	double P = (1. / (2. * pow(this->coord_MAX - this->coord_MIN, 2.))) * (E2 - SumRF);
	this->dump_FILE->print_p_Calculation(P);
	double lambda = sqrt(1. + (this->delta_T / _tau) * ((_P0 / P) - 1.));
	return lambda;
}
// метод расчета параметра лямбда для пропорционального термостата
double MD_VV::lambda_Proportional_Thermostat(double _T0, double _k)
{
	// k > 1
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	double T = E / (2. * this->k_B * this->particle_COUNT);
	double lambda = pow((_T0 / T), (1. / (2. * _k)));
	this->dump_FILE->print_t_Calculation(T);
	return lambda;
}
//
void MD_VV::load_LJ_forces()
{
	double rij = .5;
	double force = 0.;

	for (int i = 0; i < this->m_Size; i++)
	{
		rij += this->displ;
		force = 24. * (2. * pow((1. / rij), 13.) - pow((1. / rij), 7.));
		this->table_pot[i] = force;
	}
}
//
void MD_VV::load_forces()
{
	Vector3D<double>* dr = new Vector3D<double>(3);
	double dx = 0., dy = 0., dz = 0.;
	double particle_mass_i = 0., particle_mass_j = 0.;
	double d = this->coord_MAX - this->coord_MIN;
	int ri = 0;
	double dF = 0.;

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_Acceleration_Vector(0., 0., 0.);
	}

	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		particle_mass_i = this->ParticleMassive[i].get_Particle_Mass();

		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			particle_mass_j = this->ParticleMassive[j].get_Particle_Mass();

			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];

			dx = dx - (int)(dx / d) * d;
			dy = dy - (int)(dy / d) * d;
			double r = sqrt(dx * dx + dy * dy);
			dr->set_Vector(dx, dy, dz);

			if (r < this->r_Cut)
			{
				ri = (int)((r - .5) / this->displ);
				dF = this->table_pot[ri];

				this->ParticleMassive[i].set_Acceleration_Vector(this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] + dF / (particle_mass_i * r),			// aXi
																 this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] + dF / (particle_mass_i * r),			// aYi
																 0.);																													// aZi //this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] + this->force(*dR).get_Vector()[2] / this->ParticleMassive[i].get_Particle_Mass());
				this->ParticleMassive[j].set_Acceleration_Vector(this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[0] - dF / (particle_mass_j * r),			// aXj
																 this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[1] - dF / (particle_mass_j * r),			// aYj
																 0.);																													// aZj //this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[2] - this->force(*dR).get_Vector()[2] / this->ParticleMassive[j].get_Particle_Mass());
			}
		}
	}
}