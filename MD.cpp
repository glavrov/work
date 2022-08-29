#include "MD.h"
// конструкторы
MD::MD()
{
	this->particle_COUNT = 20;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->dump_FILE = new Dump();
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->coord_MIN = 0.;
	this->coord_MAX = 1.;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = 0.;
	this->particle_Temperature = 10.;
	this->coeff_Ax = 1;
	this->coeff_Ay = 1;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD::MD(int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->dump_FILE = new Dump();
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->coord_MIN = 0.;
	this->coord_MAX = 1.;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = 0.;
	this->particle_Temperature = 10.;
	this->coeff_Ax = 1;
	this->coeff_Ay = 1;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD::MD(int _particle_COUNT, double _delta_T)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->dump_FILE = new Dump();
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->coord_MIN = 0.;
	this->coord_MAX = 10.;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = _delta_T;
	this->particle_Temperature = 10.;
	this->coeff_Ax = 1;
	this->coeff_Ay = 1;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD::MD(int _particle_COUNT, double _coord_MIN, double _coord_MAX, double _velocity_MIN, double _velocity_MAX)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle[this->particle_COUNT];
	this->dump_FILE = new Dump();
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->coord_MIN = _coord_MIN;
	this->coord_MAX = _coord_MAX;
	this->velocity_MIN = _velocity_MIN;
	this->velocity_MAX = _velocity_MAX;
	this->delta_T = 0.;
	this->particle_Temperature = 10.;
	this->coeff_Ax = 1;
	this->coeff_Ay = 1;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// деструктор
MD::~MD()
{}
// метод инициализации массива частиц, сгенерированных случайным образом (равномерное распределение)
void MD::Initialization()
{
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_Particle_ID(i);
		this->ParticleMassive[i].set_Radius_Vector(this->coord_MIN, this->coord_MAX);
		this->ParticleMassive[i].set_Velocity_Vector(this->velocity_MIN, this->velocity_MAX);
	}
	// расчёт векторов средних скоростей и среднеквадратичных скоростей
	this->SumV->set_Vector(this->Average_velocity(0), this->Average_velocity(1), this->Average_velocity(2));
	this->SumV2->set_Vector(this->Average_velocity2(0), this->Average_velocity2(1), this->Average_velocity(2));
	// вывод радиус-веторов и скоростей частиц после генерации до вызова метода "нормализации"
	//this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
	//this->dump_FILE->print_vVector();
}
// метод нормализации векторов скоростей частиц
void MD::velocity_Normalization()
{
	// коэффициент нормализации
	this->velocity_Coefficient->set_Vector(sqrt((2. * this->particle_Temperature) / this->SumV2->get_Vector()[0]), sqrt((2. * this->particle_Temperature) / this->SumV2->get_Vector()[1]), 0.); // sqrt((2. * this->particle_Temperature) / this->SumV2->get_Vector()[2]));
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
	this->dump_FILE->print_vVector();
	// вычисление радиус-векторов на i0-1 шаге
	double ri0_1x = 0., ri0_1y = 0., ri0_1z = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		ri0_1x = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] * this->delta_T;
		ri0_1y = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] * this->delta_T;
		//ri0_1z = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] * this->delta_T;
		this->ParticleMassive[i].set_previous_rVector(ri0_1x, ri0_1y, ri0_1z);
	}
	// рассчет кинетической, потенциальной и полной энергий в момент времени t
	this->E_calculation();
}
// метод движения частиц - реализация классического метода Верле
void MD::Move_СV(int _i)
{
	Particle* temp = new Particle[this->particle_COUNT];	// временный массив указателей для копирования радиус-векторов частиц в момент времени t
	for (int i = 0; i < this->particle_COUNT; i++)			// записываем во временный массив указателей радиус-векторы в момент времени t
	{
		temp[i].set_Radius_Vector(this->ParticleMassive[i].get_Radius_Vector());
	}
	// рассчет ускорений
	this->compute_Acceleration();
	// рассчет радиус-векторов в момент времени t+Δt
	this->Next_Radius_Vector(_i);
	// рассчет скоростей частиц в момент времени t+Δt
	this->Next_Velocity_Vector();
	// рассчет кинетической, потенциальной и полной энергий в момент времени t+Δt
	this->E_calculation();
	// текущий радиус-вектор становится предыдущим
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_previous_rVector(temp[i].get_Radius_Vector());
	}
	// печать в файл
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
}
// метод для моделирования эволюции молекулярной динамики на основе классического метода Верле
void MD::Evolution_СV(int _step_COUNT, int _space_TYPE)
{
	cout << "[0%";
	for (int i = 0; i < _step_COUNT; i++)
	{
		if ((i + 1) % 100 == 0)
			this->v_Normalization();
		this->Move_СV(_space_TYPE);
		if (i % 1000 == 0)
			cout << ".";
	}
	cout << "100%]";
}
// рассчет радиус-векторов в момент времени t+Δt
void MD::Next_Radius_Vector(int _i)
{
	double x_NEW = 0., y_NEW = 0., z_NEW = 0.;
	if (_i == 1)
	{
		for (int i = 0; i < this->particle_COUNT; i++)
		{
			x_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[0] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T * this->delta_T;
			y_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[1] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T * this->delta_T;
			z_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[2] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T * this->delta_T;
			if (x_NEW > this->coord_MAX)
				x_NEW -= this->coord_MAX;
			else if (x_NEW < this->coord_MIN)
				x_NEW += this->coord_MAX - this->coord_MIN;
			if (y_NEW > this->coord_MAX)
				y_NEW -= this->coord_MAX;
			else if (y_NEW < this->coord_MIN)
				y_NEW += this->coord_MAX - this->coord_MIN;
			if (z_NEW > this->coord_MAX)
				z_NEW -= this->coord_MAX;
			else if (z_NEW < this->coord_MIN)
				z_NEW += this->coord_MAX - this->coord_MIN;
			this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
		}
	}
	if (_i == 2)
	{
		for (int i = 0; i < this->particle_COUNT; i++)
		{
			x_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[0] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T * this->delta_T * this->coeff_Ax;
			y_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[1] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T * this->delta_T * this->coeff_Ay;
			z_NEW = 2. * this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[i].get_previous_rVector().get_Vector()[2] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T * this->delta_T;
			if (x_NEW > this->coord_MAX)
			{
				x_NEW = this->coord_MAX;
				this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
				this->coeff_Ax = -this->coeff_Ax;
			}
			else if (x_NEW < this->coord_MIN)
			{
				x_NEW = this->coord_MIN;
				this->ParticleMassive[i].set_Velocity_Vector(-this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
				this->coeff_Ax = -this->coeff_Ax;
			}
			else
			{
				this->coeff_Ax = 1;
			}
			if (y_NEW > this->coord_MAX)
			{
				y_NEW = this->coord_MAX;
				this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], -this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
				this->coeff_Ay = -this->coeff_Ay;
			}
			else if (y_NEW < this->coord_MIN)
			{
				y_NEW = this->coord_MIN;
				this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], -this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
				this->coeff_Ay = -this->coeff_Ay;
			}
			else
			{
				this->coeff_Ay = 1;
			}
			/*if (z_NEW > this->coord_MAX)
			{
				z_NEW = this->coord_MAX;
				this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], -this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
			}
			else if (z_NEW < this->coord_MIN)
			{
				z_NEW = this->coord_MIN;
				this->ParticleMassive[i].set_Velocity_Vector(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0], -this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1], this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2]);
			}*/
			this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
		}
	}
}
// метод расчета скоростей частиц в момент t+Δt
void MD::Next_Velocity_Vector()
{
	double V_X = 0., V_Y = 0., V_Z = 0.;
	double lambda = 0.;

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		V_X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] * this->delta_T;
		V_Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] * this->delta_T;
		//V_Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] * this->delta_T);
		this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
	}
	// лямбда(t)
	//lambda = this->lambda_Proportional_Thermostat(.01, 2.);
	lambda = this->lambda_Thermostat(.01, .1);
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		V_X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] * lambda;
		V_Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] * lambda;
		//V_Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] * lambda;
		this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
	}
}
// метод нормализации скоростей частиц
void MD::v_Normalization()
{
	// координаты после нормализации
	double velocity_X_NEW = 0., velocity_Y_NEW = 0., velocity_Z_NEW = 0.;
	// нормализация векторов скоростей
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		velocity_X_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[0] - this->Average_velocity(0);
		velocity_Y_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[1] - this->Average_velocity(1);
		//velocity_Z_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[2] - this->Average_velocity(2)) * this->velocity_Coefficient->get_Vector()[2];
		this->ParticleMassive[i].set_Velocity_Vector(velocity_X_NEW, velocity_Y_NEW, velocity_Z_NEW);
	}
}
// расчет средней скорости
double MD::Average_velocity(int coordinate_Number)
{
	double average_Coordinate = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		average_Coordinate += this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[coordinate_Number];
	average_Coordinate /= this->particle_COUNT;
	return average_Coordinate;
}
// расчет среднеквадратической скорости
double MD::Average_velocity2(int coordiname_Number)
{
	double average_Coordinate2 = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		average_Coordinate2 += pow(this->ParticleMassive[i].get_Velocity_Vector().get_Vector()[coordiname_Number], 2.);
	}
	average_Coordinate2 /= this->particle_COUNT;
	return average_Coordinate2;
}
// рассчет ускорений для метода Верле
void MD::compute_Acceleration()
{
	Vector3D<double>* dR = new Vector3D<double>(3);
	double dx = 0., dy = 0., dz = 0.;
	double particle_RADIUS = 1. / sqrt(3.);

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_aVector(0., 0., 0.);
	}

	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];
			
			/*if (dx <= particle_RADIUS && dx >= 0)
				dx = particle_RADIUS;
			else if (dx >= -particle_RADIUS && dx < 0)
				dx = -particle_RADIUS;

			if (dy <= particle_RADIUS && dy >= 0)
				dy = particle_RADIUS;
			else if (dy >= -particle_RADIUS && dy < 0)
				dy = -particle_RADIUS;*/

			//if (dz <= particle_RADUIS && dz >= 0)
			//	dz = particle_RADUIS;
			//else if (dz >= -particle_RADUIS && dz < 0)
			//	dz = -particle_RADUIS;

			dR->set_Vector(dx, dy, dz);

			this->ParticleMassive[i].set_aVector(this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[0] + this->force(*dR).get_Vector()[0],			// aX
												 this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[1] + this->force(*dR).get_Vector()[1],			// aY
												 0.);//this->ParticleMassive[i].get_Acceleration_Vector().get_Vector()[2] + this->force(*dR).get_Vector()[2]);	// aZ
			this->ParticleMassive[j].set_aVector(this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[0] - this->force(*dR).get_Vector()[0],			// aX
												 this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[1] - this->force(*dR).get_Vector()[1],			// aY
												 0.);//this->ParticleMassive[j].get_Acceleration_Vector().get_Vector()[2] - this->force(*dR).get_Vector()[2]);	// aZ
		}
	}
}
// расчет кинетической, потенциальной и полной энергий
void MD::E_calculation()
{
	Vector3D<double>* dr = new Vector3D<double>(3);
	double dx = 0., dy = 0., dz = 0.;
	double particle_RADIUS = 1. / sqrt(3.);
	double at = 0.;
	// кинетическая энергия
	for (int i = 0; i < this->particle_COUNT; i++)
		this->E_kinetic += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	this->E_kinetic /= 2.;
	// потенциальная энергия
	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[0] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[0];
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[1] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[1];
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector()[2] - this->ParticleMassive[j].get_Radius_Vector().get_Vector()[2];

			if (dx <= particle_RADIUS && dx >= 0)
				dx = particle_RADIUS;
			else if (dx >= -particle_RADIUS && dx < 0)
				dx = -particle_RADIUS;

			if (dy <= particle_RADIUS && dy >= 0)
				dy = particle_RADIUS;
			else if (dy >= -particle_RADIUS && dy < 0)
				dy = -particle_RADIUS;

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
// рассчет силы, действующей на частицы, для метода Верле
Vector3D<double> MD::force(Vector3D<double>& _dR)
{
	double rl = _dR.module();
	double r6 = pow(rl, 6.);
	double fr = -((48. / r6) * (1. / r6 - .5)) / pow(rl, 2.);
	Vector3D<double>* force = new Vector3D<double>(3);
	force->set_Vector(fr * _dR.get_Vector()[0], fr * _dR.get_Vector()[1], fr * _dR.get_Vector()[2]);
	return *force;
}
// вычисление потенциала Леннард-Джонса
double MD::LJ_potential(Vector3D<double>& _dr)
{
	double rm = _dr.module();
	double rm6 = pow(rm, 6);
	double ur = (4. * ((1. / rm6) - 1)) / rm6;
	return ur;
}
// метод расчета параметра лямбда для термостата Берендсена
double MD::lambda_Thermostat(double _T0, double _tau)
{
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	double T = E / (2. * this->k_B * this->particle_COUNT);
	double lambda = sqrt(1. + (this->delta_T / _tau) * ((_T0 / T) - 1.));
	this->dump_FILE->print_t_Calculation(T);
	return lambda;
}
// метод расчета параметра лямбда для баростата Берендсена
double MD::lambda_Barostat(double _P0, double _tau)
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
	//this->dump_FILE->print_p_Calculation(P);
	double lambda = sqrt(1. + (this->delta_T / _tau) * ((_P0 / P) - 1.));
	return lambda;
}
// метод расчета параметра лямбда для пропорционального термостата
double MD::lambda_Proportional_Thermostat(double _T0, double _k)
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