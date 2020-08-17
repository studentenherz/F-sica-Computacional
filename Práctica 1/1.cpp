#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include "util.h"

using namespace std;
const float fpi = acos(-1);	 // pi float
const double dpi = acos(-1); // pi double

// Derivative o(h) not centered
template <class Type>
Type d_h1(Type f(Type), Type x, Type h)
{
	return (f(x + h) - f(x)) / h;
}

// Derivative o(h^2) centered
template <class Type>
Type d_h2(Type f(Type), Type x, Type h)
{
	return (f(x + h) - f(x - h)) * 0.5 / h;
}

// Derivative o(h^4) centered
template <class Type>
Type d_h4(Type f(Type), Type x, Type h)
{
	return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
}

// Error calculation from method and actual derivative
template <class Type>
Type err(Type method(Type(Type), Type, Type), Type f(Type), Type df(Type), Type x, Type h)
{
	return abs(method(f, x, h) - df(x));
}

// Templated function
template <class Type>
Type f(Type x)
{
	return sin(x);
}

// Templated derivative
template <class Type>
Type df(Type x)
{
	return cos(x);
}

int main()
{
	float xf = fpi / 6;
	double xd = dpi / 6;

	vector<float> h, f1, f2, f4;
	vector<double> d1, d2, d4;

	// Los cálculos
	for (int i = 0; i < 20; ++i)
	{
		double h_double = pow(10, -i);
		h.push_back((float)h_double);

		// Con precisión simple
		f1.push_back(err(d_h1, f, df, xf, h[i]));
		f2.push_back(err(d_h2, f, df, xf, h[i]));
		f4.push_back(err(d_h4, f, df, xf, h[i]));

		// Con precisión doble
		d1.push_back(err(d_h1, f, df, xd, h_double));
		d2.push_back(err(d_h2, f, df, xd, h_double));
		d4.push_back(err(d_h4, f, df, xd, h_double));
	}

	// Convertir los valores calculados con presición simple a double
	// para poder ponerlos en el mismo plot todos
	vector<double> dh(h.begin(), h.end());
	vector<double> df1(f1.begin(), f1.end());
	vector<double> df2(f2.begin(), f2.end());
	vector<double> df4(f4.begin(), f4.end());

	// Esto es un poco de lo que me "inventé" para exportar en gnuplot
	Plot<double> plt;
	plt.add(dh);
	plt.add(df1, "float o(h)");
	plt.add(d1, "double o(h)", "w lp pointtype 5 pointsize 0.5");
	plt.add(df2, "float o(h^2)");
	plt.add(d2, "double o(h^2)", "w lp pointtype 5  pointsize 0.5");
	plt.add(df4, "float o(h^4)");
	plt.add(d4, "double o(h^4)", "w lp pointtype 5  pointsize 0.5");

	plt.setLabels({"h", "error"});
	plt.setLogScale({true, true});
	plt.setXRange(pow(10, -21), 10);
	plt.setYRange(pow(10, -14), 1000);

	plt.plot("ej1");
	return 0;
}
