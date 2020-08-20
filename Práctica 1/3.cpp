#include "util.h"
#include <cmath>
#include <time.h>
#include <stdlib.h>

using namespace std;

typedef vector<double> vd;

// Clase para el sistema de ecuaciones con parámetros A y B
class EqSystem
{
	double A, B; // parámetros del sistema

public:
	EqSystem(double a, double b) : A(a), B(b) {}
	~EqSystem() {}

	double fx(double x, double y)
	{
		return A - (B + 1) * x + pow(x, 2) * y;
	}

	double fy(double x, double y)
	{
		return B * x - pow(x, 2) * y;
	}
};

// Runge-Kutta órden 4
pair<vd, vd> RK4(EqSystem sys, double x0, double y0, double T, double h)
{
	vd x, y;									 // para guardar los resultados
	double k1x, k2x, k3x, k4x; // constantes para x
	double k1y, k2y, k3y, k4y; // constantes para y

	// valores iniciales
	x.push_back(x0);
	y.push_back(y0);

	for (int i = 0; i < (int)(T / h); ++i) // (T / h) pasos
	{
		// primera aproximación
		k1x = h * sys.fx(x[i], y[i]);
		k1y = h * sys.fy(x[i], y[i]);

		// RK2
		k2x = h * sys.fx(x[i] + 0.5 * k1x, y[i] + 0.5 * k1y);
		k2y = h * sys.fy(x[i] + 0.5 * k1x, y[i] + 0.5 * k1y);

		// RK2 con k2
		k3x = h * sys.fx(x[i] + 0.5 * k2x, y[i] + 0.5 * k2y);
		k3y = h * sys.fy(x[i] + 0.5 * k2x, y[i] + 0.5 * k2y);

		// k3 en el paso entero
		k4x = h * sys.fx(x[i] + k3x, y[i] + k3y);
		k4y = h * sys.fy(x[i] + k3x, y[i] + k3y);

		// agregar la solución
		x.push_back(x[i] + (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0);
		y.push_back(y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0);

		// lanzar excepción si el método diverge
		if (isinf(x[i + 1]) || isinf(y[i + 1]) || isnan(x[i + 1]) || isnan(y[i + 1]))
		{
			throw "Infinite or NAN value found with initial conditions (" + to_string(x0) + ", " + to_string(y0) + ") and h = " + to_string(h);
			return {x, y};
		}
	}

	return {x, y};
}

int main()
{

	double T = 10;								 // tiempo de evaluación
	vd vec_h = {1, 0.5, 0.2, 0.1}; // valores dados de h

	Plot<double> plt1("A = 1; B=1");

	EqSystem stable(1, 1); // sistema estable con punto fijo en (1,1)

	pair<vd, vd> v;
	for (double h : vec_h) // para cada paso
	{
		v = RK4(stable, 1.45, 1.45, T, h);

		plt1.add(v.first, "", true);
		plt1.add(v.second, "h = " + to_string(h).substr(0, 3), false);
	}

	plt1.setLabels({"x", "y"});
	plt1.plot("ej3_a_1");

	try
	{
		Plot<double> plt2("A = 1; B = 1");

		srand(time(0));
		for (int i = 0; i < 8; ++i)
		{
			// condiciones iniciales aleatorias
			double x0 = 2.0 * rand() / RAND_MAX;
			double y0 = 2.0 * rand() / RAND_MAX;

			for (int j = 0; j < (int)vec_h.size(); ++j) // y para cada valor de h
			{
				v = RK4(stable, x0, y0, T, vec_h[j]);

				plt2.add(v.first, "", true);
				plt2.add(v.second, (i == 0 ? "h = " + to_string(vec_h[j]).substr(0, 3) : ""), false, "w lp linetype " + to_string(j + 3) + " pointtype 7 pointsize 0.5");
			}
		}

		plt2.setLabels({"x", "y"});
		plt2.plot("ej3_a_2");
	}
	catch (string e)
	{
		cout << e << '\n';
	}

	EqSystem unstable(1, 3);
	Plot<double> plt3("A = 1; B = 3 (h = 0.1)"); // sistema inestable con ciclo límite
	plt3.setLabels({"x", "y"});

	double h = 0.1; // un paso conservador

	srand(time(0));
	for (int i = 0; i < 10; i++)
	{
		// condiciones iniciales aleatorias
		double x0 = 5.0 * rand() / RAND_MAX;
		double y0 = 4.0 * rand() / RAND_MAX;

		v = RK4(unstable, x0, y0, T, h);

		plt3.add(v.first, "", true);
		plt3.add(v.second, "(" + to_string(x0).substr(0, 3) + ", " + to_string(y0).substr(0, 3) + ") ", false, "w lp  pointsize 0.5 pointtype 7");
	}

	for (int i = 0; i < 5; i++)
	{
		// condiciones iniciales aleatorias pero cerca del punto fijo
		double x0 = 0.5 + 0.5 * rand() / RAND_MAX;
		double y0 = 2.5 + 0.5 * rand() / RAND_MAX;

		v = RK4(unstable, x0, y0, T, h);

		plt3.add(v.first, "", true);
		plt3.add(v.second, "(" + to_string(x0).substr(0, 3) + ", " + to_string(y0).substr(0, 3) + ") ", false, "w lp  pointsize 0.5 pointtype 7");
	}

	// con la condición inicial en el punto de equilibrio iniestable
	double x0 = 1;
	double y0 = 3;

	v = RK4(unstable, x0, y0, T, h);

	plt3.add(v.first, "", true);
	plt3.add(v.second, "(" + to_string(x0).substr(0, 3) + ", " + to_string(y0).substr(0, 3) + ") ", false, "w lp  pointsize 0.5 pointtype 7");

	plt3.plot("ej3_b");

	// condición inicial dentro del ciclo
	x0 = 1.5;
	y0 = 3.5;

	T = 10;
	// Aumentando h desde 0.1 a pasos de 0.01 hasta que diverja el método
	try
	{
		for (int j = 0; j < 50; ++j)
		{
			RK4(unstable, x0, y0, T, h + 0.01 * j);
		}
	}
	catch (const string e)
	{
		cout << e << '\n';
	}

	// ahora con condición inicial fuera del ciclo límite
	x0 = 3.1;
	y0 = 3;

	// Aumentando h desde 0.1 a pasos de 0.01 hasta que diverja el método
	try
	{
		for (int j = 0; j < 50; ++j)
		{
			RK4(unstable, x0, y0, T, h + 0.01 * j);
		}
	}
	catch (const string e)
	{
		cout << e << '\n';
	}

	return 0;
}