#include "MD_Round_VV.h"
// контрукторы
MD_Round_VV::MD_Round_VV()
{
	this->particle_COUNT = 20;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->dump_FILE->set_particle_COUNT(this->particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->R = 8.;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = 0.;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_Round_VV::MD_Round_VV(int _particle_COUNT)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->dump_FILE->set_particle_COUNT(_particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->R = 8.;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = 0.;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_Round_VV::MD_Round_VV(int _particle_COUNT, double _R)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->dump_FILE->set_particle_COUNT(_particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->R = _R;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = 0.;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_Round_VV::MD_Round_VV(int _particle_COUNT, double _R, double _delta_T)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->dump_FILE->set_particle_COUNT(_particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->R = _R;
	this->velocity_MIN = -1.;
	this->velocity_MAX = 1.;
	this->delta_T = _delta_T;

	this->particle_RADIUS = 2. / (sqrt(this->particle_COUNT - 2));

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}

MD_Round_VV::MD_Round_VV(int _particle_COUNT, double _R, double _velocity_MIN, double _velocity_MAX)
{
	this->particle_COUNT = _particle_COUNT;
	this->ParticleMassive = new Particle_VV[this->particle_COUNT];
	this->dump_FILE = new Dump_VV();
	this->dump_FILE->set_particle_COUNT(_particle_COUNT);
	this->SumV = new Vector3D<double>(3);
	this->SumV2 = new Vector3D<double>(3);
	this->velocity_Coefficient = new Vector3D<double>(3);
	this->R = _R;
	this->velocity_MIN = _velocity_MIN;
	this->velocity_MAX = _velocity_MAX;
	this->delta_T = 0.;

	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// деструктор
MD_Round_VV::~MD_Round_VV()
{}
// метод инициализации массива частиц, сгенерированных случайным образом (равномерное распределение)
void MD_Round_VV::Initialization()
{
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		this->ParticleMassive[i].set_Particle_ID(i);
		this->ParticleMassive[i].set_Radius_Vector_Circle(0., 0., this->R);
		this->ParticleMassive[i].set_Velocity_Vector(this->velocity_MIN, this->velocity_MAX);
		this->ParticleMassive[i].set_Flag(0);
	}
	// расчет векторов средних скоростей и среднеквадратичных скоростей
	this->SumV->set_Vector(this->Average_velocity(0), this->Average_velocity(1), this->Average_velocity(2));
	this->SumV2->set_Vector(this->Average_velocity2(0), this->Average_velocity2(1), this->Average_velocity2(2));
	// вывод радиус-веторов и векторов скоростей частиц после генерации до вызова метода "нормализации"
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
	// печать координат границы области в файл
	this->dump_FILE->print_Borderline(this->R);
}
// метод нормализации векторов скоростей частиц
void MD_Round_VV::Velocity_Normalization()
{
	// коэффициент нормализации
	this->velocity_Coefficient->set_Vector(sqrt((2. * this->system_Temperature) / this->SumV2->get_Vector().at(0)), sqrt((2. * this->system_Temperature) / this->SumV2->get_Vector().at(1)), 0.); // sqrt((2. * this->particle_Temperature) / this->SumV2->get_Vector().at(2)));
	// рассчет средних скоростей по координатным осям
	double average_Velocity_X = this->Average_velocity(0), average_Velocity_Y = this->Average_velocity(1); // average_Velocity_Z = this->Average_velocity(2);
	// координаты после нормализации
	double velocity_X_NEW = 0., velocity_Y_NEW = 0., velocity_Z_NEW = 0.;
	// нормализация векторов скоростей
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		velocity_X_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(0) - average_Velocity_X) * this->velocity_Coefficient->get_Vector().at(0);
		velocity_Y_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(1) - average_Velocity_Y) * this->velocity_Coefficient->get_Vector().at(1);
		//velocity_Z_NEW = (this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(2) - average_Velocity_Z)* this->velocity_Coefficient->get_Vector().at(2);
		this->ParticleMassive[i].set_Velocity_Vector(velocity_X_NEW, velocity_Y_NEW, velocity_Z_NEW);
	}
	// печать векторов скоростей в файл
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_vVector();
	// расчет кинетической, потенциальной и полной энергий в момент времени t
	this->E_calculation_VV();
}
// метод движения частиц - реализация скоростного метода Верле
void MD_Round_VV::Move_VV()
{
	// расчет скоростей частиц в момент времени t+Δt/2
	this->Next_Velocity_Vector_05t();
	// расчет радиус-векторов в момент времени t+Δt
	this->Next_Radius_Vector();
	// расчет ускорений в момент времени t+Δt
	this->Compute_Acceleration();
	// расчет скоростей частиц в момент времени t+Δt
	this->Next_Velocity_Vector();
	// расчет кинетической, потенциальной и полной энергий в момент времени t+Δt
	this->E_calculation_VV();
	// печать в файл
	this->dump_FILE->set_ParticleMassive(this->ParticleMassive);
	this->dump_FILE->print_rVector();
	// обнуление параметра-флага у каждой частицы, у которой он отличен от 0 (нуля)
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		if (this->ParticleMassive[i].get_Flag() == 1) this->ParticleMassive[i].set_Flag(0);
	}
}
// метод моделирования эволюции молекулярной динамики на основе скоростного метода Верле
void MD_Round_VV::Evolution_VV(int _step_COUNT)
{
	cout << "[0%";
	for (int i = 0; i < _step_COUNT; i++)
	{
		//if ((i + 1) % 100 == 0)
		//	this->v_Normalization();
		this->Move_VV();
		if (i % 5000 == 0)
			cout << ".";
		//if ((i + 1) % 100000 == 0)
		//	this->T0 /= 10.;
		//if ((i + 1) % 100 == 0)
		//	this->dump_FILE->print_rVector();
	}
	cout << "100%]";
}
// расчет радиус-векторов в момент времени t+Δt
void MD_Round_VV::Next_Radius_Vector()
{
	double x = 0., y = 0.;							// текущие координаты радиус-вектора
	double x_NEW = 0., y_NEW = 0., z_NEW = 0.;		// новые координаты радиус-вектора
	double Vx = 0., Vy = 0., Vz = 0.;				// координаты вектора скорости
	double lambda = 0.;								// параметр для баростата Берендсена

	double k = 0., b = 0.;							// прямая, проходящая через точки Xi и Xi+1
	double x0 = 0., y0 = 0.;						// координаты точки пересечения касательной и окружности
	double D = 0.;									// дискриминант квадратного уравнения
	double x01 = 0., x02 = 0.;						// первый корень квадратного уравнения
	double y01 = 0., y02 = 0.;						// второй корень квадратного уравнения

	double Vx_ = 0., Vy_ = 0.;						// компоненты вектора скорости в новой системе координат
	double sinA = 0., cosA = 0.;					// sin, cos угла поворота координатных осей

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		x = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0);
		y = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1);

		x_NEW = x + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(0);
		y_NEW = y + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(1);
		//z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(2) + this->delta_T * this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(2);

		double r = sqrt(x_NEW * x_NEW + y_NEW * y_NEW + z_NEW * z_NEW);

		// если радиус-вектор частицы в момент времени t+Δt больше радиуса пространства
		if (r > this->R)
		{
			// у соответствующей частицы устанавливаем параметр флага равным 1, означающий начала рассчета процесса отражения от стенки
			this->ParticleMassive[i].set_Flag(1);

			Vx = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(0);
			Vy = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(1);

			k = (y - y_NEW) / (x - x_NEW);
			b = (x * y_NEW - x_NEW * y) / (x - x_NEW);

			// образуем систему уравнений из окружности и уравнения прямой, на которой лежит вектор,
			// при этом центр окружности совпадает с началом координат:
			// { y = kx + b
			// { x^2 + y^2 = R^2

			D = sqrt((pow(k, 2.) + 1) * pow(this->R, 2.) - pow(b, 2.));

			x01 = (-k * b + D) / (pow(k, 2.) + 1.);
			x02 = (-k * b - D) / (pow(k, 2.) + 1.);
			
			y01 = k * x01 + b;
			y02 = k * x02 + b;

			// выбираем координаты той точки пересечения прямой и окружности, расстояние от которой до частицы меньше
			if (sqrt(pow(x01 - x, 2.) + pow(y01 - y, 2.)) < sqrt(pow(x02 - x, 2.) + pow(y02 - y, 2.)))
			{
				x0 = x01;
				y0 = y01;
			}
			else
			{
				x0 = x02;
				y0 = y02;
			}

			x_NEW = (2. * x0 * ((x0 * x + y0 * y) / (pow(x0, 2.) + pow(y0, 2)))) - x;
			y_NEW = (2. * y0 * ((x0 * x + y0 * y) / (pow(x0, 2.) + pow(y0, 2)))) - y;

			sinA = y0 / this->R;
			cosA = x0 / this->R;

			Vx_ = Vx * cosA + Vy * sinA;	// отраженные компоненты скорости в новой системе координат
			Vy_ = -Vx * sinA + Vy * cosA;

			Vx_ *= -1.;

			Vx = Vx_ * cosA - Vy_ * sinA;	// отраженные компоненты скорости системе координат Oxy
			Vy = Vx_ * sinA + Vy_ * cosA;

			this->ParticleMassive[i].set_Velocity_Vector(Vx, Vy, Vz);
		}
		
		this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
	}
	// лямбда(t)
	//lambda = this->lambda_Barostat(.001, 1.2);
	//for (int i = 0; i < this->particle_COUNT; i++)
	//{
	//	x_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0) * lambda;
	//	y_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1) * lambda;
	//	z_NEW = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(2) * lambda;
	//	this->ParticleMassive[i].set_Radius_Vector(x_NEW, y_NEW, z_NEW);
	//S}
}
// расчет скоростей частиц в момент времени t+Δt/2
void MD_Round_VV::Next_Velocity_Vector_05t()
{
	double V_05X = 0., V_05Y = 0., V_05Z = 0.;
	double lambda = 0.;

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		V_05X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(0) + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(0) * this->delta_T / 2.;
		V_05Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(1) + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(1) * this->delta_T / 2.;
		//V_05Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(2) + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(2) * this->delta_T / 2.;

		this->ParticleMassive[i].set_Velocity_Vector_05t(V_05X, V_05Y, V_05Z);
	}
}
// расчет скоростей частиц в момент времени t+Δt
void MD_Round_VV::Next_Velocity_Vector()
{
	double V_X = 0., V_Y = 0., V_Z = 0.;
	double lambda = 1.;

	//lambda = this->lambda_Thermostat(.1);
	lambda = this->lambda_Proportional_Thermostat(2.);

	for (int i = 0; i < this->particle_COUNT; i++)
	{
		if (this->ParticleMassive[i].get_Flag() == 1)
			continue;
		else
		{
			V_X = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(0) * lambda + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(0) * this->delta_T / 2.; //
			V_Y = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(1) * lambda + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(1) * this->delta_T / 2.; // 
			//V_Z = this->ParticleMassive[i].get_Velocity_Vector_05t().get_Vector().at(2) * lambda + this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(2) * this->delta_T / 2.;
			this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
		}
	}
	// лямбда(t)
	//lambda = this->lambda_Thermostat(.1);
	//lambda = this->lambda_Proportional_Thermostat(2.);
	//for (int i = 0; i < this->particle_COUNT; i++)
	//{
	//	V_X = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(0) * lambda;
	//	V_Y = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(1) * lambda;
	//	//V_Z = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(2) * lambda;
	//	this->ParticleMassive[i].set_Velocity_Vector(V_X, V_Y, V_Z);
	//}
}
// метод нормализации скоростей частиц
void MD_Round_VV::v_Normalization()
{
	// средние координаты скорости
	double average_Velocity_X = this->Average_velocity(0), average_Velocity_Y = this->Average_velocity(1); // average_Velicity_Z = this->Average_velovity(2);
	// координаты после нормализации
	double velocity_X_NEW = 0., velocity_Y_NEW = 0., velocity_Z_NEW = 0.;
	// нормализация векторов скоростей
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		velocity_X_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(0) - average_Velocity_X;
		velocity_Y_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(1) - average_Velocity_Y;
		//velocity_Z_NEW = this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(2) - average_Velocity_Z;
		this->ParticleMassive[i].set_Velocity_Vector(velocity_X_NEW, velocity_Y_NEW, velocity_Z_NEW);
	}
}
// расчет средней скорости
double MD_Round_VV::Average_velocity(int coordinate_Number)
{
	double average_Coordinate = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		average_Coordinate += this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(coordinate_Number);
	}
	average_Coordinate /= this->particle_COUNT;
	return average_Coordinate;
}
// расчет среднеквадратической скорости
double MD_Round_VV::Average_velocity2(int coordiname_Number)
{
	double average_Coordinate2 = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		average_Coordinate2 += pow(this->ParticleMassive[i].get_Velocity_Vector().get_Vector().at(coordiname_Number), 2.);
	}
	average_Coordinate2 /= this->particle_COUNT;
	return average_Coordinate2;
}
// расчет ускорений для метода Верле
void MD_Round_VV::Compute_Acceleration()
{
	Vector3D<double>* dR = new Vector3D<double>(3);
	double aXi = 0., aYi = 0., aZi = 0.;
	double aXj = 0., aYj = 0., aZj = 0.;
	double particle_Mass_i = 0., particle_Mass_j = 0.;

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

//#pragma omp parallel for
	for (int i = 0; i < this->particle_COUNT; i++)
		this->ParticleMassive[i].set_Acceleration_Vector(0., 0., 0.);

	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		particle_Mass_i = this->ParticleMassive[i].get_Particle_Mass();

		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			particle_Mass_j = this->ParticleMassive[j].get_Particle_Mass();

			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(2) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(2);

			r = sqrt(dx * dx + dy * dy + dz * dz);

			if (r < this->particle_RADIUS)
			{
				x0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0);
				y0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1);
				xj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
				yj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);

				k = dy / dx;
				b = (x0 * yj - xj * y0) / dx;
				
				double D = /* 4. * */ pow(k * b - x0 - y0 * k, 2.) - /* 4. * */ (k * k + 1) * (x0 * x0 + b * b - 2. * y0 * b + y0 * y0 - this->particle_RADIUS * this->particle_RADIUS); // дискриминант кв. ур-ния (x-x0)^2 + (y-y0)^2 = pR^2
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
				
				if (sqrt(x_razn1 * x_razn1 + y_razn1 * y_razn1) <= sqrt(x_razn2 * x_razn2 + y_razn2 * y_razn2))
				{
					//this->ParticleMassive[j].set_Radius_Vector(u1, u2, 0.);

					//dx = x0 - x1;
					//dy = y0 - y1;

					if (x0 > xj)
						dx = x0 - x1;
					else
						dx = x1 - x0;

					if (y0 > yj)
						dy = y0 - y1;
					else
						dy = y1 - y0;
				}
				else
				{
					//this->ParticleMassive[j].set_Radius_Vector(u1, u2, 0.);

					//dx = x0 - x2;
					//dy = y0 - y2;

					if (x0 > xj)
						dx = x0 - x2;
					else
						dx = x2 - x0;

					if (y0 > yj)
						dy = y0 - y2;
					else
						dy = y2 - y0;
				}
			}

			dR->set_Vector(dx, dy, dz);

			/*Vector3D<double>* force = new Vector3D<double>(3);
			*force = this->force(*dR);

			this->ParticleMassive[i].set_Acceleration_Vector(this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(0) + force->get_Vector().at(0),
															 this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(1) + force->get_Vector().at(1), 
															 aZi);
			this->ParticleMassive[i].set_Acceleration_Vector(this->ParticleMassive[j].get_Acceleration_Vector().get_Vector().at(0) - force->get_Vector().at(0),
															 this->ParticleMassive[j].get_Acceleration_Vector().get_Vector().at(1) - force->get_Vector().at(1), 
															 aZi);*/

			this->ParticleMassive[i].set_Acceleration_Vector(this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(0) + this->force(*dR).get_Vector().at(0) / particle_Mass_i,	// aXi
															 this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(1) + this->force(*dR).get_Vector().at(1) / particle_Mass_i,	// aYi
															 aZi);	// aZi //this->ParticleMassive[i].get_Acceleration_Vector().get_Vector().at(2) + this->force(*dR).get_Vector().at(2) / this->ParticleMassive[i].get_Particle_Mass());
			this->ParticleMassive[j].set_Acceleration_Vector(this->ParticleMassive[j].get_Acceleration_Vector().get_Vector().at(0) - this->force(*dR).get_Vector().at(0) / particle_Mass_j,	// aXj
															 this->ParticleMassive[j].get_Acceleration_Vector().get_Vector().at(1) - this->force(*dR).get_Vector().at(1) / particle_Mass_j,	// aYj
															 aZj);	// aZj //this->ParticleMassive[j].get_Acceleration_Vector().get_Vector().at(2) - this->force(*dR).get_Vector().at(2) / this->ParticleMassive[j].get_Particle_Mass());

			this->E_potential += this->LJ_potential(*dR);
		}
	}

	/*
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
	*/
}
// метод торможения системы частиц
void MD_Round_VV::Deceleration()
{
	double sum_X = 0., sum_Y = 0., module = 0.;						// накопленная сумма по координатам и расстояние между частицами (i, j)
	double Xi = 0., Yi = 0., Zi = 0., Xj = 0., Yj = 0., Zj = 0.;	// значения координат частиц (i, j)
	
	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		Xi = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0);
		Yi = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1);

		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			Xj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
			Yj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);
			// нахождение расстояния между частицами i и j (два способа)
			// module = pow(sqrt(pow(Xi - Xj, 2.) + pow(Yi - Yj, 2.)), 3.);
			module = this->ParticleMassive[i].get_Radius_Vector().distance(this->ParticleMassive[j].get_Radius_Vector());
			
			sum_X += (Xi - Xj) / module;
			sum_Y += (Yi - Yj) / module;
		}
	}
	// расчет новых координат радиус-векторов системы частиц после торможения системы
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		Xi = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0);
		Yi = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1);

		Xi = -Xi + sum_X;
		Yi = -Yi + sum_Y;

		this->ParticleMassive[i].set_Radius_Vector(Xi, Yi, Zi);
	}
}
// расчет кинетической, потенциальной и полной энергий
void MD_Round_VV::E_calculation_VV()
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
		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(2) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(2);

			double r = sqrt(dx * dx + dy * dy + dz * dz);

			if (r < this->particle_RADIUS)
			{
				double x0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0);
				double y0 = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1);
				double xj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
				double yj = this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);

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

				double x_razn1 = x1 - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
				double y_razn1 = y1 - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);
				double x_razn2 = x2 - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
				double y_razn2 = y2 - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);

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
	// вывод расчитанных значений энергий в файл
	this->dump_FILE->print_e_Calculation(this->E_kinetic, this->E_potential, this->E_full);
	// обнуление полей, отвечающих за значения энергий
	this->E_potential = 0.;
	this->E_kinetic = 0.;
	this->E_full = 0.;
}
// расчет силы, действующей на частицы, для метода Верле
Vector3D<double> MD_Round_VV::force(Vector3D<double>& _dR)
{
	double rl = _dR.module();
	double r6 = pow(rl, 6.);
	double fr = ((48. / r6) * (1. / r6 - .5)) / pow(rl, 2.);			// полный потенциал Леннард-Джонса
	// double fr = ((48. / r6) * (1. / r6)) / pow(rl, 2.);				// только отталкивание
	// double fr = ((48. / (r6 * r6)) / pow(rl, 2.);					// только отталкивание
	// double fr = 1. / pow(rl, 2.);									// "чистый" Кулон
	Vector3D<double>* force = new Vector3D<double>(3);
	force->set_Vector(fr * _dR.get_Vector().at(0), fr * _dR.get_Vector().at(1), 0.); // fr* _dR.get_Vector().at(2));
	return *force;
}
// вычисление потенциала Леннард-Джонса
double MD_Round_VV::LJ_potential(Vector3D<double>& _dr)
{
	double rm = _dr.module();
	double rm6 = pow(rm, 6);
	double ur = (4. * ((1. / rm6) - 1.)) / rm6;		// полный потенциал Леннард-Джонса
	//double ur = (4. * ((1. / rm6))) / rm6;		// только отталкивание
	//double ur = 1. / rm;							// "чистый" Кулон
	return ur;
}
// метод расчета параметра лямбда для термостата Берендсена
double MD_Round_VV::lambda_Thermostat(double _tau)
{
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	double T = E / (2. * this->k_B * this->particle_COUNT);
	double lambda = sqrt(1. + ((this->delta_T / _tau) * ((this->T0 / T) - 1.)));
	this->dump_FILE->print_t_Calculation(T);
	return lambda;
}
// метод расчета параметра лямбда для баростата Берендсена
double MD_Round_VV::lambda_Barostat(double _P0, double _tau)
{
	double E = 0.;
	for (int i = 0; i < this->particle_COUNT; i++)
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Radius_Vector().module(), 2.);
	double E2 = E / this->particle_COUNT;

	Vector3D<double>* dR = new Vector3D<double>(3);
	double rXi = 0., rYi = 0., rZi = 0.;
	double dx = 0., dy = 0., dz = 0.;

	double SumRF = 0.;
	for (int i = 0; i < this->particle_COUNT - 1; i++)
	{
		for (int j = i + 1; j < this->particle_COUNT; j++)
		{
			dx = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(0) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(0);
			dy = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(1) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(1);
			//dz = this->ParticleMassive[i].get_Radius_Vector().get_Vector().at(2) - this->ParticleMassive[j].get_Radius_Vector().get_Vector().at(2);

			dR->set_Vector(dx, dy, dz);

			SumRF += dR->module() * this->ParticleMassive[i].get_Acceleration_Vector().module(); // m=1
		}
	}
	double P = (1. / (2. * M_PI * pow(this->R, 2.))) * (E2 - SumRF);
	this->dump_FILE->print_p_Calculation(P);
	double lambda = sqrt(1. + ((this->delta_T / _tau) * ((_P0 / P) - 1.)));
	return lambda;
}
// метод расчета параметра лямбда для пропорционального термостата
double MD_Round_VV::lambda_Proportional_Thermostat(double _k)
{
	// k > 1
	double E = 0.;
	double power = 1. / (2. * _k);
	for (int i = 0; i < this->particle_COUNT; i++)
	{
		E += this->ParticleMassive[i].get_Particle_Mass() * pow(this->ParticleMassive[i].get_Velocity_Vector().module(), 2.);
	}
	double T = E / (2. * this->k_B * this->particle_COUNT);
	double lambda = pow((this->T0 / T), power);
	this->dump_FILE->print_t_Calculation(T);
	return lambda;
}