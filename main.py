import materials_properties
import materials


a = materials.Material(materials_properties.Tungsten['melting_temperature'],
                       materials_properties.Tungsten['boiling_temperature'],
                       materials_properties.Tungsten['eps'])

print('mt:', a.melting_temperature)
"""

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <fstream>
#include "Matherial.h"

void MeshCreator(double const &, double const &, double const &, int const &, int const &, double*, double*, int &, int &);
void StateScreen(int const &, int const &);
void MeshCorrector(double const &, double const &, double const &, double const &, int &);
void OutputFunction(std::ofstream *, std::ofstream *, int const &, int const &, int const &, double const &, double *, double *, double **, double const &, double const, int &);
void LogOfSolution(int const &, double const &);
void BoilingFunction(double const &, double const &, double const &, double const &, int &);
void OutputMeltingSurface(std::ofstream *, int const &, int const &, double *, double *, double **, double const &, double const);

double DeltaOfPhaseTransition(int const &, int const &, double **, double const &);
double ThermalConductivityOfSteel40H(double const &);
double HeatCapacityOfSteel40H(double const &);
double ThermalConductivity(double const &, int &, int const &);
double HeatCapacity(double const &, int &, double const &, int const &);
double Density(double const &, int &, int const &);
double HeatFluxOfPlasmatron(double const &, double const &, double const &);

int main()
{
	int
		Indicator(-1),                              // Переменная-условие продолжения расчета.
		iteration(0),                               // Переменная хранящая кол-во итераций по времени.
		OutputIteration(0),                         // Переменная для хранения кол-во итераций вывода результатов в файл.
		SeparationSurface(0);                       // Переменная хранящая узел "разделяющий" слои материалов.
	const int
		K = 1000,                                    // Количество узлов по ширине.
		N = 500;                                    // Количество узлов по глубине.
	double
		PresentTime(0),                             // Переменная для хранения значения теперешнего времени.
		MeltingTemperature,                         // Переменная в которой хранится температура кипения материала.
		BoilingPoint,                               // Точка кипения.
		Delta,                                      // Температурная разность в области фазового перехода.
		a, b, c, d, l1, l2, r1, r2, h1, h2, y, eps, // Переменные используемые в разностной схеме.
		**temperature, *zn = new double[N + 1], *rk = new double[K + 1], *fz = new double[N + 1], *gz = new double[N + 1], *fr = new double[K + 1], *gr = new double[K + 1];
	const double
		BeginTemperature = 300,                     // Начальная температура образца (К).
//////////////////////////////////////////////////////
		H = 5e-4,                                   // Глубина образца (м).
		R = 11e-4,                                  // Ширина образца (м).
		ThicknessOfFirstLayer = 2e-4,               // Толщина внешнего слоя (м).
		EndTime = (1000)*1e-6;                      // Конечное время расчета (с) - а так-же удвоенное время воздействия плазменной струи.
	const int
		MatherialOfFirstLayer = 3,                  // Номер материала верхнего слоя.
//////////////////////////////////////////////////////
		NumberOfSteps = EndTime * 1e5;              // Количество итераций по времени.
	const double
		sigma = 5.67e-8,                            // Постоянная Стефана-Больцмана.
		TimeStep = EndTime / NumberOfSteps,         // Шаг по времени.
		OutputTime = EndTime / 50;                  // Вывод результатов в файл каждую 1/50 от общего времени.

	Matherial Steel40H(1600, 3000, 0.1);            // 1
	Matherial Tungsten(3695, 6000, 0.3);            // 2
	Matherial Molybdenum(2895, 5000, 0.1);          // 3

	std::ofstream TemperatureFields;
	TemperatureFields.open("TemperatureFields.csv");
	std::ofstream MaxTemperatureInTop;
	MaxTemperatureInTop.open("MaxTemperatureInTop.csv");
	std::ofstream TemperatureZR;
	TemperatureZR.open("TemperatureZR.csv");
	std::ofstream MeltingSurface;
	MeltingSurface.open("MeltingSurface.csv");

	if ((ThicknessOfFirstLayer == 0) || (MatherialOfFirstLayer == 1))
		Indicator = 9;
	switch (MatherialOfFirstLayer)
	{
	case 1:
		MeltingTemperature = Steel40H.GetTemperatureOfMelting();
		BoilingPoint = Steel40H.GetBoilingTemperature();
		eps = Steel40H.GetEps();
		break;
	case 2:
		MeltingTemperature = Tungsten.GetTemperatureOfMelting();
		BoilingPoint = Tungsten.GetBoilingTemperature();
		eps = Tungsten.GetEps();
		break;
	case 3:
		MeltingTemperature = Molybdenum.GetTemperatureOfMelting();
		BoilingPoint = Molybdenum.GetBoilingTemperature();
		eps = Molybdenum.GetEps();
		break;
	}
	temperature = new double *[N + 1];                    // Обьявление двумерного массива для хранения температур.
	for (int i(0); i < N + 1; i++)
		temperature[i] = new double[K + 1];
	for (int i(0); i <= N; i++)                           // Инициализация двумерного масива для хранения температур (НАЧАЛЬНОЕ УСЛОВИЕ).
		for (int j(0); j <= K; j++)
			temperature[i][j] = BeginTemperature;
	MeshCreator(H, R, ThicknessOfFirstLayer, K, N, zn, rk, SeparationSurface, Indicator);
	while ((PresentTime < EndTime) && (Indicator == -1))
	{
		MaxTemperatureInTop << "\n" << ";" << PresentTime  * 1e6 << ";" << temperature[0][0] << ";" << temperature[SeparationSurface + 1][0];
		iteration++;
		StateScreen(NumberOfSteps, iteration);
		PresentTime = iteration * TimeStep;
		Delta = DeltaOfPhaseTransition(N, K, temperature, BoilingPoint);
		for (int j(0); j <= K; j++)
		{
			h1 = zn[1] - zn[0];
			a = ThermalConductivity(temperature[0][j], Indicator, MatherialOfFirstLayer) / h1;
			b = -(a + eps * sigma * pow(temperature[0][j], 3));
			d = HeatFluxOfPlasmatron(rk[j], PresentTime, EndTime);
			fz[0] = -a / b;
			gz[0] = -d / b;
			for (int i(1); i <= SeparationSurface - 1; i++)
			{
				h1 = zn[i] - zn[i - 1];
				h2 = zn[i + 1] - zn[i];
				a = -2 * ThermalConductivity((temperature[i][j] + temperature[i + 1][j]) / 2.0, Indicator, MatherialOfFirstLayer) / ((h1 + h2) * h2);
				c = -2 * ThermalConductivity((temperature[i][j] + temperature[i - 1][j]) / 2.0, Indicator, MatherialOfFirstLayer) / ((h1 + h2) * h1);
				b = HeatCapacity(temperature[i][j], Indicator, Delta, MatherialOfFirstLayer) / TimeStep - a - c;
				d = -temperature[i][j] * HeatCapacity(temperature[i][j], Indicator, Delta, MatherialOfFirstLayer) / TimeStep;
				fz[i] = -a / (b + c * fz[i - 1]);
				gz[i] = -(d + c * gz[i - 1]) / (b + c * fz[i - 1]);
			}
			a = -ThermalConductivityOfSteel40H(temperature[SeparationSurface][j]) / (zn[SeparationSurface + 1] - zn[SeparationSurface]);
			c = -ThermalConductivity(temperature[SeparationSurface][j], Indicator, MatherialOfFirstLayer) / (zn[SeparationSurface] - zn[SeparationSurface - 1]);
			b = -c - a;
			d = 0;
			fz[SeparationSurface] = -a / (b + c * fz[SeparationSurface - 1]);
			gz[SeparationSurface] = -(d + c * gz[SeparationSurface - 1]) / (b + c * fz[SeparationSurface - 1]);
			for (int i(SeparationSurface + 1); i <= N - 1; i++)
			{
				h1 = zn[i] - zn[i - 1];
				h2 = zn[i + 1] - zn[i];
				a = -2 * ThermalConductivityOfSteel40H((temperature[i][j] + temperature[i + 1][j]) / 2.0) / ((h1 + h2) * h2);
				c = -2 * ThermalConductivityOfSteel40H((temperature[i][j] + temperature[i - 1][j]) / 2.0) / ((h1 + h2) * h1);
				b = HeatCapacityOfSteel40H(temperature[i][j]) / TimeStep - a - c;
				d = -temperature[i][j] * HeatCapacityOfSteel40H(temperature[i][j]) / TimeStep;
				fz[i] = -1.0*a / (b + c * fz[i - 1]);
				gz[i] = -1.0*(d + c * gz[i - 1]) / (b + c * fz[i - 1]);
			}
			b = 1;
			c = -1;
			d = 0;
			temperature[N][j] = -(d + c * gz[N - 1]) / (b + c * fz[N - 1]);
			for (int i(N - 1); i >= 0; i--)
				temperature[i][j] = fz[i] * temperature[i + 1][j] + gz[i];
		}
		for (int i(0); i <= SeparationSurface; i++)
		{
			a = 1;
			b = -1;
			d = 0;
			fr[0] = -a / b;
			gr[0] = -d / b;
			for (int j(1); j <= K - 1; j++)
			{
				l1 = rk[j] - rk[j - 1];
				l2 = rk[j + 1] - rk[j];
				r1 = rk[j] - l1 / 2;
				r2 = rk[j] + l2 / 2;
				y = (r2 * r2 - r1 * r1) / 2;
				a = -r2 * ThermalConductivity((temperature[i][j] + temperature[i][j + 1]) / 2, Indicator, MatherialOfFirstLayer) / (y * l2);
				c = -r1 * ThermalConductivity((temperature[i][j] + temperature[i][j - 1]) / 2, Indicator, MatherialOfFirstLayer) / (y * l1);
				b = HeatCapacity(temperature[i][j], Indicator, Delta, MatherialOfFirstLayer) / TimeStep - a - c;
				d = -temperature[i][j] * HeatCapacity(temperature[i][j], Indicator, Delta, MatherialOfFirstLayer) / TimeStep;
				fr[j] = -a / (b + c * fr[j - 1]);
				gr[j] = -(d + c *gr[j - 1]) / (b + c * fr[j - 1]);
			}
			b = 1;
			c = -1;
			d = 0;
			temperature[i][K] = -(d + c *gr[K - 1]) / (b + c *fr[K - 1]);
			for (int j(K - 1); j >= 0; j--)
				temperature[i][j] = fr[j] * temperature[i][j + 1] + gr[j];
		}
		for (int i(SeparationSurface + 1); i <= N; i++)
		{
			a = 1;
			b = -1;
			d = 0;
			fr[0] = -a / b;
			gr[0] = -d / b;
			for (int j(1); j <= K - 1; j++)
			{
				l1 = rk[j] - rk[j - 1];
				l2 = rk[j + 1] - rk[j];
				r1 = rk[j] - l1 / 2;
				r2 = rk[j] + l2 / 2;
				y = (r2 * r2 - r1 * r1) / 2;
				a = -r2 * ThermalConductivityOfSteel40H((temperature[i][j] + temperature[i][j + 1]) / 2) / (y * l2);
				c = -r1 * ThermalConductivityOfSteel40H((temperature[i][j] + temperature[i][j - 1]) / 2) / (y * l1);
				b = HeatCapacityOfSteel40H(temperature[i][j]) / TimeStep - a - c;
				d = -temperature[i][j] * HeatCapacityOfSteel40H(temperature[i][j]) / TimeStep;
				fr[j] = -a / (b + c * fr[j - 1]);
				gr[j] = -(d + c *gr[j - 1]) / (b + c * fr[j - 1]);
			}
			b = 1;
			c = -1;
			d = 0;
			temperature[i][K] = -(d + c *gr[K - 1]) / (b + c *fr[K - 1]);
			for (int j(K - 1); j >= 0; j--)
				temperature[i][j] = fr[j] * temperature[i][j + 1] + gr[j];
		}
		BoilingFunction(temperature[0][0], temperature[0][SeparationSurface + 1], BoilingPoint, Steel40H.GetBoilingTemperature(), Indicator);
		MeshCorrector(PresentTime, temperature[N][0], temperature[0][K], BeginTemperature, Indicator);
		if ((OutputIteration * OutputTime <= PresentTime) && (PresentTime != EndTime))
			OutputFunction(&TemperatureZR, &TemperatureFields, N, K, SeparationSurface, ThicknessOfFirstLayer, zn, rk, temperature, PresentTime, MeltingTemperature, OutputIteration);
		if (temperature[0][0] >= MeltingTemperature)
			OutputMeltingSurface(&MeltingSurface, N, K, zn, rk, temperature, PresentTime, MeltingTemperature);
		if (PresentTime == EndTime)
		{
			MaxTemperatureInTop << "\n" << ";" << PresentTime  * 1e6 << ";" << temperature[0][0] << ";" << temperature[SeparationSurface + 1][0];
			Indicator = 0;
		}
	}
	LogOfSolution(Indicator, PresentTime);
	delete[] zn, rk, fz, gz, fr, gr;
	for (int i(0); i < N + 1; i++)
		delete[] temperature[i];
	TemperatureFields.close();
	MaxTemperatureInTop.close();
	TemperatureZR.close();
	MeltingSurface.close();
	return 0;
}

// Функция создающая расчетную область.
void MeshCreator(double const &rH, double const &rR, double const &rThicknessOfFirstLayer, int const &rK, int const &rN, double *zn, double *rk, int &rSeparationSurface, int &rIndicator)
{
	int
		i, switcher(0);
	double
		qr, qz, hz, lr, delta, deltaPrev;
	const int
		qrK_1 = 2,                                  // Параметр прогрессии по ширине.
		qzN_1 = 20;                                 // Параметр прогрессии по глубине.
	if (rThicknessOfFirstLayer > rH)
		rIndicator = 1;
	qr = 1.0 * pow(qrK_1, 1.0 / (rK - 1));
	qz = 1.0 * pow(qzN_1, 1.0 / (rN - 1));
	hz = rH * (1 - qz) / (1 - pow(qz, rN));
	lr = rR * (1 - qr) / (1 - pow(qr, rK));
	zn[0] = 0;
	for (int i(1); i <= rN; i++)
		zn[i] = zn[i - 1] + hz * pow(qz, i - 1);
	rk[0] = 0;
	for (int i(1); i <= rK; i++)
		rk[i] = rk[i - 1] + lr * pow(qr, i - 1);
	if ((rThicknessOfFirstLayer != 0))                 // Блок который находит границу "раздела материалов".
	{
		deltaPrev = rThicknessOfFirstLayer;
		i = 1;
		switcher = 0;
		do
		{
			if (switcher == 0)
			{
				delta = rThicknessOfFirstLayer - zn[i];
				if (deltaPrev > abs(delta))
				{
					deltaPrev = delta;
					i++;
				}
				else
				{
					rSeparationSurface = i - 1;
					zn[i - 1] = rThicknessOfFirstLayer;
					switcher = 1;
				}
			}
			else
				break;
		} while (i < rK - 1);
	}
}
// Функция выводящая прогресс расчета на экран.
void StateScreen(int const &rNumberOfSteps, int const &riteration)
{
	system("cls");
	std::cout << "\n\n\t\t\t\t\t\t\t\t" << 100 * riteration / rNumberOfSteps << "%" << "(" << riteration << ")" << std::endl;
	std::cout << "\t";
	for (int i(0); i <= (1.0 * riteration / rNumberOfSteps) * 100; i++)
		std::cout << (char)219;
}
// Функция которая выдает сообщения которые помогают корректировать расчетную область.
void MeshCorrector(double const &rPresentTime, double const &rTemperatureAtTheBottomOfMesh, double const &rTemperatureAtTheRightOfMesh, double const &rBeginTemperature, int &rIndicator)
{
	double DifferenceOfTemperature = 50;           // Максимально допустимая разница между начальной температурой и минимальной температурой на границах расчетной области.
	if (rTemperatureAtTheBottomOfMesh - rBeginTemperature > DifferenceOfTemperature)
		rIndicator = 2;
	if (rTemperatureAtTheRightOfMesh - rBeginTemperature > DifferenceOfTemperature)
		rIndicator = 3;
}
// Функция вывода в файл результатов.
void OutputFunction(std::ofstream *TemperatureZR, std::ofstream *File, int const &rN, int const &rK, int const &rSeparationSurface, double const &rThicknessOfFirstLayer, double *zn, double *rk, double **temperature, double const &rPresentTime, double const rMeltingTemperature, int &rOutputIteration)
{
	int
		i(0), j(0);
	const int
		sz = 4,                                     // Вывод каждого ...-го распределения температуры в файл.
		sr = 4;
	rOutputIteration++;
	*File << "Time:" << ";" << rPresentTime * 1e6 << ";" << "mksec\n" << std::endl;
	*File << ";" << ";" << ";";
	while (j < rK)
	{
		*File << "r" << j << ";";
		j += sr;
	}
	*File << "r" << rK << std::endl;
	j = 0;
	*File << ";" << "mkm" << ";" << ";";
	while (j < rK)
	{
		*File << rk[j] * 1e6 << "; ";
		j = j + sr;
	}
	*File << rk[rK] * 1e6 << std::endl;
	*File << std::endl;
	i = 0;
	while (i < rSeparationSurface)
	{
		*File << "z" << i << ";" << zn[i] * 1e6 << ";";
		j = 0;
		while (j <= rK)
		{
			*File << ";" << temperature[i][j];
			j += sr;
		}
		*File << std::endl;
		i += sz;
	}
	*File << "z" << rSeparationSurface << ";" << zn[rSeparationSurface] * 1e6 << ";";
	j = 0;
	while (j <= rK)
	{
		*File << ";" << temperature[rSeparationSurface][j];
		j += sr;
	}
	*File << std::endl;
	i = rSeparationSurface + sr;
	while (i < rN)
	{
		*File << "z" << i << ";" << zn[i] * 1e6 << ";";
		j = 0;
		while (j < rK)
		{
			*File << ";" << temperature[i][j];
			j += sr;
		}
		*File << ";" << temperature[i][rK];
		*File << std::endl;
		i += sz;
	}
	*File << "z" << rN << ";" << zn[rN] * 1e6 << ";";
	j = 0;
	while (j < rK)
	{
		*File << ";" << temperature[rN][j];
		j += sr;
	}
	*File << ";" << temperature[rN][rK];
	*File << std::endl << std::endl << std::endl;
	////////////////////////////////////////////////////////////////////////////////
	*TemperatureZR << "\n" << std::endl;
	*TemperatureZR << "Time:" << ";" << rPresentTime * 1e6 << ";" << "mksec\n" << std::endl;
	*TemperatureZR << ";" << ";" << ";";
	j = 0;
	while (j < rK)
	{
		*TemperatureZR << "r" << j << ";";
		j += sr;
	}
	*TemperatureZR << "r" << rK << std::endl;
	j = 0;
	*TemperatureZR << "R (mkm):" << ";" << ";" << ";";
	while (j < rK)
	{
		*TemperatureZR << rk[j] * 1e6 << "; ";
		j = j + sr;
	}
	*TemperatureZR << rk[rK] * 1e6 << std::endl;
	*TemperatureZR << "Temperature (K):" << ";" << ";";
	j = 0;
	while (j <= rK)
	{
		*TemperatureZR << ";" << temperature[0][j];
		j += sr;
	}
	*TemperatureZR << "\n";
	*TemperatureZR << ";" << ";" << ";";
	i = 0;
	while (i < rN)
	{
		*TemperatureZR << "z" << i << ";";
		i += sr;
	}
	*TemperatureZR << "z" << rN << std::endl;
	i = 0;
	*TemperatureZR << "Z (mkm)" << ";" << ";" << ";";
	while (i < rN)
	{
		*TemperatureZR << zn[i] * 1e6 << "; ";
		i = i + sr;
	}
	*TemperatureZR << zn[rN] * 1e6 << std::endl;
	*TemperatureZR << "Temperature (K)" << ";" << ";";
	i = 0;
	while (i <= rN)
	{
		*TemperatureZR << ";" << temperature[i][0];
		i += sr;
	}
}
// Функция которая выводит информацию о статусе расчета (ошибках) по его завершеню.
void LogOfSolution(int const & rIndicator, double const & rPresentTime)
{
	std::cout << "\7";
	switch (rIndicator)
	{
	case 0: std::cout << "\n\nSuccess solution" << std::endl; break;
	case 1: std::cout << "\n\nFunction:\tMeshCreator" << "\n\tThikkness of first layer is bigger, than general" << std::endl; break;
	case 2: std::cout << "\n\nFunction:\tMeshCorrector" << "\ntime:\t\t " << rPresentTime * 1e6 << "\n\tIncrease H\n"; break;
	case 3: std::cout << "\n\nFunction:\tMeshCorrector" << "\ntime:\t\t " << rPresentTime * 1e6 << "\n\tIncrease R\n"; break;
	case 4: std::cout << "\n\nFunction:\tThermalConductivity" << "\ntime:\t\t " << rPresentTime * 1e6 << "\nERROR:\t\tinput temperature for tungsten is out of domain of definition of the function\n"; break;
	case 5: std::cout << "\n\nFunction:\tHeatCapacity" << "\ntime:\t\t " << rPresentTime * 1e6 << "\nERROR:\t\tinput temperature for tungsten is out of domain of definition of the function\n"; break;
	case 6: std::cout << "\n\nFunction:\tDensity" << "\ntime:\t\t" << rPresentTime * 1e6 << "\nERROR:\t\tinput is out of domain of definition of the function\n"; break;
	case 7: std::cout << "\n\nFunction:\tBoilingFunction" << "\ntime:\t\t" << rPresentTime * 1e6 << "\nERROR:\t\tFirst layer is boiling\n"; break;
	case 8: std::cout << "\n\nFunction:\tBoilingFunction" << "\ntime:\t\t" << rPresentTime * 1e6 << "\nERROR:\t\tSecond layer is boiling\n"; break;
	case 9: std::cout << "\n\n" << "\ntime:\t\t" << rPresentTime * 1e6 << "\nERROR:\t\tNeed to set Matherial of the first layer > 1 / Thickness of the first layer != 0\n"; break;
	case 10: std::cout << "\n\nFunction:\tDensity" << "\ntime:\t\t" << rPresentTime * 1e6 << "\nERROR:\t\tinput temperature for molybdenum is out of domain of definition of the function\n"; break;
	case 11: std::cout << "\n\nFunction:\tHeatCapacity" << "\ntime:\t\t " << rPresentTime * 1e6 << "\nERROR:\t\tinput temperature for molybdenum is out of domain of definition of the function\n"; break;
	case 12: std::cout << "\n\nFunction:\tThermalConductivity" << "\ntime:\t\t " << rPresentTime * 1e6 << "\nERROR:\t\tinput temperature for molybdenum is out of domain of definition of the function\n"; break;
	default: std::cout << "\n\nSomething else goes wrong" << std::endl;
	}
	std::cout << std::endl;
	system("pause");
}
// Функция предупреждающая о достижении температуры кипения.
void BoilingFunction(double const &rMaxTopTemperature, double const &rMaxDownTemperature, double const &rBoilingTemperature, double const &rBoilingTemperatureOfBottomMatherial, int &rIndicator)
{
	if (rMaxTopTemperature >= rBoilingTemperature)
		rIndicator = 7;
	if (rMaxDownTemperature >= rBoilingTemperatureOfBottomMatherial)
		rIndicator = 8;
}
// Функция выводящая поверхность фазового перехода.
void OutputMeltingSurface(std::ofstream *MeltingSurface, int const &rN, int const &rK, double *zn, double *rk, double **temperature, double const &rPresentTime, double const rMeltingTemperature)
{
	*MeltingSurface << ";" << "Time:" << ";" << rPresentTime  * 1e6 << std::endl;
	for (int i(0); i <= rN; i++)
	{
		if (temperature[i][0] >= rMeltingTemperature)
		{
			for (int j(0); j <= rK; j++)
			{
				if (temperature[i][j] < rMeltingTemperature)
				{
					*MeltingSurface << ";" << ";" << ";" << "z" << ";" << zn[i] * 1e6 << ";" << "r" << ";" << rk[j - 1] * 1e6 << std::endl;
					break;
				}
			}
		}
		else break;
	}
	*MeltingSurface << std::endl;
}

// Функция находящая толщину слоя фазового перехода.
double DeltaOfPhaseTransition(int const &rN, int const &rK, double **temperature, double const & rBoilingPoint)
{
	int sdel = 4;                   // Какой-то параметр.
	double  delta;
	for (int i(sdel); i < rN + 1 - sdel; i++)
	{
		if (temperature[i][0] < rBoilingPoint)
		{
			delta = temperature[i - 1][0] - temperature[i - 1 + sdel][0];
			if (temperature[i - sdel][0] - temperature[i][0] > delta)
			{
				delta = temperature[i - sdel][0] - temperature[i][0];
			}
			break;
		}
	}
	for (int j(sdel); j < rK + 1 - sdel; j++)
	{
		if (temperature[0][j] < rBoilingPoint)
		{
			if (temperature[0][j - 1] - temperature[0][j - 1 + sdel] > delta)
			{
				delta = temperature[0][j - 1] - temperature[0][j - 1 + sdel];
			}
			if (temperature[0][j - sdel] - temperature[0][j] > delta)
			{
				delta = temperature[0][j - sdel] - temperature[0][j];
			}
			break;
		}
	}
	return delta;
}
// Функция теплопроводности steel40H.
double ThermalConductivityOfSteel40H(double const &rtemperature)
{
	return 40.62 + 0.013 * rtemperature - 4.847 / 100000 * rtemperature * rtemperature + 2.405 / 100000000 * pow(rtemperature, 3);
}
// Функция теплоемкости steel40H.
double HeatCapacityOfSteel40H(double const &rtemperature)
{
	return (364.726 + 0.407 * rtemperature - 1.048 / 10000 * pow(rtemperature, 2)) * (7.918 * 1000 - 0.32 * rtemperature);
}
// Функция теплопроводности для поверхностного материала.
double ThermalConductivity(double const &rtemperature, int &rIndicator, int const & rNumberOfMatherial)
{
	switch (rNumberOfMatherial)
	{
	case 2:                                                         // Функция теплопроводности вольфрама.
		if ((rtemperature >= 300) && (rtemperature <= 3695))
			return 149.441 - 45.466 / 1000 * rtemperature + 13.193 / 1000000 * pow(rtemperature, 2) - 1.484 / 1000000000 * pow(rtemperature, 3) + 3.866 * 1000000 / pow(rtemperature, 2);
		if ((rtemperature > 3695) && (rtemperature <= 6000))
			return 66.6212 + 0.02086 * (rtemperature - 3695) - 3.7585 / 1000000 * pow((rtemperature - 3695), 2);
		else
			rIndicator = 4;
		break;
	case 3:                                                         // Функция теплопроводности Molybdenum.
		if ((rtemperature >= 293.15) && (rtemperature <= 2895))
			return 111.1349 - 0.01236 * rtemperature + 1.4122 / 1000000.0*rtemperature*rtemperature;
		if ((rtemperature > 2895) && (rtemperature <= 5000))
			return 11.4262 + 0.0208 * rtemperature;
		else
			rIndicator = 12;
		break;
	}
}
// Функция возвращающая произведение теплоемкости на плотность для поверхностного материала.
double HeatCapacity(double const &rtemperature, int &rIndicator, double const &rDelta, int const & rNumberOfMatherial)
{
	switch (rNumberOfMatherial)
	{
	case 2:                                                     // Функция возвращающая произведение теплоемкости вольфрама на его плотность.
		if ((rtemperature >= 293.15) && (rtemperature <= 3080))
			return ((21.868372 + 8.068661 / 1000 * rtemperature - 3.756196 / 1000000 * pow(rtemperature, 2) + 1.075862 / 1000000000 * pow(rtemperature, 3) + 1.406637 * 10000 / pow(rtemperature, 2)) / 0.18384) * Density(rtemperature, rIndicator, rNumberOfMatherial);
		if ((rtemperature > 3080) && (rtemperature < 3695 - rDelta))
			return (((2.022 + 1.315 / 100 * rtemperature) / 0.18384) * Density(rtemperature, rIndicator, rNumberOfMatherial));
		if ((rtemperature >= 3695 - rDelta) && (rtemperature <= 3695 + rDelta))
		{
			double
				HeatOfThePhaseTransition = 191000,                                               // (Дж/кг).
				AverageDensity = (((19.25 - 2.66207 / 10000 * (3695 - 293.15) - 3.0595 / 1000000000 * pow((3695 - 293.15), 2) - 9.5185 / 1000000000000.0 * pow((3695 - 293.15), 3)) * 1000) + ((16.267 - 7.679 / 10000 * (3695 - 3695) - 8.091 / 100000000.0 * pow((3695 - 3695), 2)) * 1000)) / 2.0,
				a = 19250,
				b = 2662.07 / 10000,
				c = 293.15,
				d = 3059.5 / 1000000000,
				f = 9518.5 / 1000000000000,
				g = 2.022 / 0.18384,
				h = 1.315 / 0.18384 / 100,
				a1 = 16267,
				b1 = 7679.0 / 10000,
				c1 = 3695,
				d1 = 8091.0 / 100000000,
				h1 = 51.03 / 0.18384;
			return (AverageDensity * HeatOfThePhaseTransition) / 2.0 / rDelta + 1.0 / 2 / rDelta * ((a*g + b*g*c - d*g*c*c + f*g*c*c*c) * (3695 - (3695 - rDelta)) + (-1.0 / 2 * b*g + d*g*c - 3.0 / 2 * f*g*c*c + 1.0 / 2 * a*h + 1.0 / 2 * b*h*c - 1.0 / 2 * d*h*c*c + 1.0 / 2 * f*h*c*c*c) * (3695.0 * 3695 - (3695 - rDelta)*(3695 - rDelta)) + (-1.0 / 3 * d*g + f*g*c - 1.0 / 3 * b*h + 2.0 / 3 * d*h*c - f*h*c*c) * (3695.0 * 3695 * 3695 - (3695 - rDelta)*(3695 - rDelta)*(3695 - rDelta)) + (-1.0 / 4 * f*g - 1.0 / 4 * d*h + 3.0 / 4 * f*h*c) * (3695.0 * 3695 * 3695 * 3695 - (3695 - rDelta)*(3695 - rDelta)*(3695 - rDelta)*(3695 - rDelta)) - 1.0 / 5 * f*h*(3695.0 * 3695 * 3695 * 3695 * 3695 - (3695 - rDelta)* (3695 - rDelta)*(3695 - rDelta)*(3695 - rDelta)*(3695 - rDelta))) + 1.0 / 2 / rDelta * ((a1*h1 + b1*c1*h1 - c1*c1*d1*h1)*((3695 + rDelta) - 3695) + (c1*d1*h1 - 1.0 / 2 * b1*h1)*((3695 + rDelta)*(3695 + rDelta) - 3695.0 * 3695) - 1.0 / 3 * d1*h1*((3695 + rDelta)*(3695 + rDelta)*(3695 + rDelta) - 3695.0 * 3695 * 3695));
		}
		if ((rtemperature > 3695 + rDelta) && (rtemperature <= 6000))
			return 51.3 / 0.18384 * Density(rtemperature, rIndicator, rNumberOfMatherial);
		else
			rIndicator = 5;
		break;
	case 3:                                            // Функция теплоемкости Molybdenum при постоянном давлении.
		if ((rtemperature >= 293.15) && (rtemperature < 2895 - rDelta))
			return (4.3 / 100000 * rtemperature*rtemperature - 0.029 * rtemperature + 274.252) * Density(rtemperature, rIndicator, rNumberOfMatherial);
		if ((rtemperature >= 2895 - rDelta) && (rtemperature <= 2895 + rDelta))
		{
			double
				HeatOfThePhaseTransition = 375234.5;
			return (Density(rtemperature, rIndicator, rNumberOfMatherial) * HeatOfThePhaseTransition) / 2.0 / rDelta + 1.0 / 2 * ((4.3 / 100000 * 2895.0*2895 - 0.029 * 2895 + 274.252) * Density(2895, rIndicator, rNumberOfMatherial)) + 1.0 / 2 * (420.4877 * Density(2895, rIndicator, rNumberOfMatherial));
		}
		if ((rtemperature > 2895 + rDelta) && (rtemperature < 5000))
			return 420.4877 * Density(rtemperature, rIndicator, rNumberOfMatherial);
		else
			rIndicator = 11;
		break;
	}

}
// Функция плотности.
double Density(double const &rtemperature, int &rIndicator, int const & rNumberOfMatherial)
{
	switch (rNumberOfMatherial)
	{
	case 2:                                            // Функция плотности вольфрама.
		if ((rtemperature >= 293.15) && (rtemperature < 3695))
			return (19.25 - 2.66207 / 10000 * (rtemperature - 293.15) - 3.0595 / 1000000000 * pow((rtemperature - 293.15), 2) - 9.5185 / 1000000000000 * pow((rtemperature - 293.15), 3)) * 1000;
		if ((rtemperature >= 3695) && (rtemperature <= 6000))
			return(16.267 - 7.679 / 10000 * (rtemperature - 3695) - 8.091 / 100000000 * pow((rtemperature - 3695), 2)) * 1000;
		else
			rIndicator = 6;
		break;
	case 3:                                            // Функция плотности молибдена.
		if ((rtemperature >= 293.15) && (rtemperature < 5000))
			return 1.02 * 10000 - 3.8 / 100 * rtemperature;
		else
			rIndicator = 10;
		break;
	}
}
// Функция вычисляющая воздействие плазмы на поверхность материала.
double HeatFluxOfPlasmatron(double const & rg, double const & rtg, double const & rEndTime)
{
	const double
		K(2.2),               // Отношение пятна нагрева к площади круга радиусом 1/гамма.
		L(10e-3),             // Растояние от среза сопла до поверхности материала.
		W(3e3),               // Энергия которая накапливается накопителем установки квазистационарного плазменного ускрителя.
		Eta0(0.8),            // Энергетический ККД разряда.
		Etat(0.3),            // Тепловой ККД импульса плазменного потока.
		r0(0.4e-3),           // Радиус сопла на срезе плазмотрона.
		alpha(1),             // Половина угла расхождения плазмового потока в его продольном разрезе.
		tg0(rEndTime / 2);    // Время действия импульса плазменной струи.
	double
		rl,                   // Радиус пятна нагрева.
		qm,                   // Максимальная плотность теплового потока которая воспринимается поверхностью материала от плазмы.
		ksi;                  // Коэфициент который определяет интенсивность изменения в радиальном направлении плотности этого теплового потока.

	rl = r0 + L * tan(alpha * M_PI / 180);
	ksi = pow(K, 0.5) / rl;
	qm = 4 * K * W * Eta0 * Etat * r0 * r0 / (M_PI * 1.0*pow(rl, 4));
	if (rtg <= tg0)
	{
		return qm * exp(-pow(rg * ksi, 2));
	}
	else return 0;
}
"""