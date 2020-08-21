#include "util.h"
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <functional>

using namespace std;

// Habrá muchos vector<double> mejor ahorramos espacio y esfuerzo con una abreviación
typedef vector<double> vd;

// Clase para el sistema de ecuaciones
// solo con tiene una lista de las funciones f_i de
// d x_i / dt = f_i(x, t)
class EqSystem
{
	vector<function<double(vd, double)>> f;

public:
	// Constructor de la clase con los parámetros A y B
	EqSystem(double A, double B)
	{
		add([A, B](vd x, double t) -> double { return (A - (B + 1) * x[0] + pow(x[0], 2) * x[1]); }); // f_x
		add([A, B](vd x, double t) -> double { return (B * x[0] - pow(x[0], 2) * x[1]); });						// f_y
	}
	// Destructor
	~EqSystem()
	{
		f.clear(); // limpiar el arreglo de funciones
	}

	void add(function<double(vd, double)> _f)
	{
		f.push_back(_f);
	}

	// para llamar la función como eqSys[i](args)
	function<double(vd, double)> operator[](int i)
	{
		return f[i];
	}
};

// Para sumar vectores miembro a miembro
template <class Type>
vector<Type> operator+(vector<Type> a, vector<Type> b)
{
	vector<Type> v(a);
	for (int i = 0; i < (int)a.size(); ++i)
		v[i] += b[i];
	return v;
}

// Para restar vectores miembro a miembro
template <class Type>
vector<Type> operator-(vector<Type> a, vector<Type> b)
{
	for (int i = 0; i < (int)a.size(); ++i)
		a[i] -= b[i];
	return a;
}

// Multiplicación escalar * vector
template <class Type>
vector<Type> operator*(Type k, vector<Type> a)
{
	for (int i = 0; i < (int)a.size(); ++i)
		a[i] *= k;
	return a;
}

// Multiplicación vector * escalar
template <class Type>
vector<Type> operator*(vector<Type> a, Type k)
{
	return k * a;
}

// Traspuesta de una matriz de dimensión 2
template <class Type>
vector<vector<Type>> transpose(vector<vector<Type>> v)
{
	vector<vector<Type>> t;
	for (int i = 0; i < (int)v[0].size(); ++i)
	{
		vector<Type> col;
		for (int j = 0; j < (int)v.size(); ++j)
			col.push_back(v[j][i]);
		t.push_back(col);
	}

	return t;
}

// Norma 2, distancia cartesiana
template <class Type>
Type norm2(vector<Type> v)
{
	Type s = 0;
	for (auto x : v)
		s += pow(x, 2);
	return sqrt(s);
}

// Norma infinito (mayor componente)
template <class Type>
Type normInf(vector<Type> v)
{
	Type s = 0;
	for (auto x : v)
		s = max(abs(x), s);
	return s;
}

// Imprimir vector
template <class Type>
void print(vector<Type> v)
{
	for (int j = 0; j < (int)v.size(); ++j)
		cout << v[j] << " ";
}

// Imprimir matriz
template <class Type>
void print(vector<vector<Type>> v)
{
	for (int i = 0; i < (int)v.size(); ++i)
	{
		print(v[i]);
		cout << "\n";
	}
}

// Runge-Kutta órden 4 un paso
vd RK4(EqSystem sys, vd x0, double t, double h)
{
	int n = x0.size(); // tamaño del vector de las variables
	vd x;							 // para guardar los resultados
	vd k1, k2, k3, k4; // constantes de RK

	// primera aproximación
	for (int i = 0; i < n; ++i)
		k1.push_back(h * sys[i](x0, t));

	// RK2
	for (int i = 0; i < n; ++i)
		k2.push_back(h * sys[i](x0 + 0.5 * k1, t));

	// RK2 con k2
	for (int i = 0; i < n; ++i)
		k3.push_back(h * sys[i](x0 + 0.5 * k2, t));

	// k3 en el paso entero
	for (int i = 0; i < n; ++i)
		k4.push_back(h * sys[i](x0 + k3, t));

	// resultado final
	x = x0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (1 / 6.0);

	return x;
}

// Resolver el sistema con un paso fijo
vector<vd> Solve(vd method(EqSystem, vd, double, double), EqSystem sys, vd x0, const double t0, const double T, const double h)
{

	int steps = 0; // número de evaluaciones del método
	vector<vd> x;	 //vector de soluciones
	x.push_back(x0);

	for (int i = 0; i < (int)((T - t0) / h); ++i) // mientras t < T
	{
		x.push_back(method(sys, x[i], t0 + i * h, h)); // aplicar un paso del método
		++steps;																			 // una evaluación más del método

		//lanzar excepción si el método diverge
		for (double _x : x[i + 1])
		{
			if (isinf(_x) || isnan(_x))
			{
				string err = "Infinite or NAN value found with:\nx0 = (";
				for (double _x0 : x0)
					err += to_string(_x0) + ", ";
				err += ")\nh = " + to_string(h) + "\n";

				throw err;

				return x;
			}
		}
	}

	cout << "Con delta t = " << T - t0 << "\n# pasos: " << steps << "\n";
	return x;
}

// Resolver el sistema con paso adaptativo
// Devuelve un par:
// vector<vd>: un vector con la pisción en cada instante dada por en vector de doubles
// vd: un vector con cada uno de los instantes de tiempo
pair<vector<vd>, vd> adaptiveSolve(vd method(EqSystem, vd, double, double), EqSystem sys, vd x0, const double t0, const double T, double h, const double e = 1e-5, const double k = 1.5)
{

	double minh = h; // h minimo
	int steps = 0;	 // número de evaluaciones del método
	vector<vd> sol;	 // soluciones
	sol.push_back(x0);
	double t = t0; // insante
	vd ts = {t0};	 // vector de tiempo

	while (t < T)
	{
		minh = min(minh, h);
		// cálculo con paso h
		vd x1 = method(sys, sol.back(), t, h);
		++steps; // un paso más

		//lanzar excepción si el método diverge
		for (double _x : x1)
		{
			if (isinf(_x) || isnan(_x))
			{
				string err = "Infinite or NAN value found with:\nx0 = (";
				for (double _x0 : x0)
					err += to_string(_x0) + ", ";
				err += ")\nh = " + to_string(h) + "\n";

				throw err;

				return {sol, ts};
			}
		}

		// cálculo con dos pasos de h/2
		vd x2_1 = method(sys, sol.back(), t, h / 2);
		vd x2 = method(sys, x2_1, t + h / 2, h / 2);
		steps += 2; // dos pasos más

		//lanzar excepción si el método diverge
		for (double _x : x2)
		{
			if (isinf(_x) || isnan(_x))
			{
				string err = "Infinite or NAN value found with:\nx0 = (";
				for (double _x0 : x0)
					err += to_string(_x0) + ", ";
				err += ")\nh = " + to_string(h) + "\n";

				throw err;

				return {sol, ts};
			}
		}

		// diferencia entre soluciones
		double delta = normInf(x1 - x2);

		if (delta > e) // si es mayor que la tolerancia
			h = h / k;	 // disminuir el paso y calcular nuevamente
		else					 // si es menor que la tolerancia
		{
			t = t + h;					// actualizar el tiempo
			ts.push_back(t);		// guardar el instante t
			sol.push_back(x2);	// aceptar x2
			if (delta <= e / 2) // si es menos que e/2
				h = h * k;				//aumentar el paso
		}
	}

	cout << "Con e = " << e << " y delta t = " << T - t0 << "\n# pasos: " << steps << " h mínimo de " << minh << "\n";
	return {sol, ts};
}

int main()
{
	double t0 = 0, // instante inicial
			T = 10.0,	 // instante final
			h = 0.01,	 // paso
			e = 1e-5;	 // toletancia

	vd x0 = {3.0, 3.0}; // posición inicial
	EqSystem sys(1, 3); // sistema con A = 1; B = 3

	Plot<double> plt0("A =1; B=3");
	Plot<double> plt1("A =1; B=3 (e = " + to_string(e) + ")");
	Plot<double> plt2("A =1; B=3 (e = " + to_string(e) + ")");
	Plot<double> plt3("A =1; B=3");
	plt3.setLabels({"x(t)", "y(t)"});
	try
	{
		// Solución con paso fijo 0.01
		vector<vd> sol0 = transpose(Solve(RK4, sys, x0, t0, T, 0.01));

		// Plot en el espacio de fases
		plt0.setLabels({"x(t)", "y(t)"});
		plt0.add(sol0[0], "", true);
		plt0.add(sol0[1], "h = 0.01", false, "w lp pt 7 ps 0.5 lc \"black\" ");

		// Solución con paso fijo 0.1
		sol0 = transpose(Solve(RK4, sys, x0, t0, T, 0.1));

		// Plot en el espacio de fases
		plt0.setLabels({"x(t)", "y(t)"});
		plt0.add(sol0[0], "", true);
		plt0.add(sol0[1], "h = 0.1", false, "w lp pt 7 ps 0.5 lc rgb \"#e09c2d\" ");

		plt0.plot("ej4_phase_fixed");

		// Solución con paso adaptativo
		pair<vector<vd>, vd> sol = adaptiveSolve(RK4, sys, x0, t0, T, h, e);
		vector<vd> coords = transpose(sol.first);
		vd t = sol.second;

		// Plot en el espacio de fases
		plt1.setLabels({"x(t)", "y(t)"});
		plt1.add(coords[0], "", true);
		plt1.add(coords[1], "", false, "w lp pt 7 ps 0.5 lc rgb \"#e09c2d\" ");
		vd ticks(coords[0].size(), 0);
		plt1.add(ticks, "", false, "w p pt \"|\" lc rgb \"black\"");
		plt1.add(ticks, "", true);
		plt1.add(coords[1], "", false, "w p pt \"_\" lc rgb \"black\"");
		plt1.plot("ej4_phase_5");

		// Plot en el tiempo
		plt2.setLabels({"t", "coordenadas (t)"});
		plt2.add(t, "", true);
		plt2.add(coords[0], "x(t)", false, "w lp pt 7 ps 0.5  lc rgb \"#e09c2d\"");
		plt2.add(coords[1], "y(t)", false, "w lp pt 7 ps 0.5  lc rgb \"black\"");
		ticks = vd(t.size(), 0);
		plt2.add(ticks, "", false, "w p pt \"|\" lt -1");
		plt2.plot("ej4_time_5");

		plt3.add(coords[0], "", true);
		plt3.add(coords[1], "e = " + to_string(e), false, "w l lw 5 lc rgb \"#d5d5d5\"");

		// Solución con paso adaptativo con k = 2
		sol = adaptiveSolve(RK4, sys, x0, t0, T, h, e, 2);
		coords = transpose(sol.first);
		t = sol.second;

		// Plot en el tiempo
		plt2.add(t, "", true);
		plt2.add(coords[0], "k = 2.0", false, "w lp pt 4 ps 0.5  lc rgb \"#e09c2d\"");
		plt2.add(coords[1], "k = 2.0", false, "w lp pt 4 ps 0.5  lc rgb \"black\"");
		ticks = vd(t.size(), 5);
		plt2.add(ticks, "", false, "w p pt \"|\" lt -1");
		plt2.plot("ej4_time_5_2");

		plt3.add(coords[0], "", true);
		plt3.add(coords[1], "e = " + to_string(e), false, "w l lw 5 lc rgb \"#d5d5d5\"");
	}
	catch (string err)
	{
		cout << err;
	}

	for (int i = 0; i < 4; i++)
	{
		e *= 10;
		try
		{
			// // Solución con paso fijo
			// vector<vd> sol = transpose(Solve(RK4, sys, x0, T, h));

			// Solución con paso adaptativo
			pair<vector<vd>, vd> sol = adaptiveSolve(RK4, sys, x0, t0, T, h, e);
			vector<vd> coords = transpose(sol.first);
			vd t = sol.second;

			// Plot en el espacio de fases
			plt3.add(coords[0], "", true);
			plt3.add(coords[1], "e = " + to_string(e), false, "w p lc \"black\" ps 0.7");
		}
		catch (string err)
		{
			cout << err;
		}
	}

	plt3.plot("ej4_phase_cmp");
	return 0;
}