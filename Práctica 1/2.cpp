#include "util.h"
#include <cmath>

using namespace std;

typedef vector<double> vd;

pair<vd, vd> euler(double x0, double p0, double T, double h)
{
	pair<vd, vd> v;
	v.first.push_back(x0);
	v.second.push_back(p0);
	for (int i = 0; i < (int)(T / h); ++i)
	{
		v.first.push_back(v.first[i] + v.second[i] * h);
		v.second.push_back(v.second[i] - v.first[i] * h);
	}

	return v;
}

vd energy(vd x, vd p)
{
	vd v;
	for (int i = 0; i < (int)x.size(); ++i)
	{
		v.push_back(0.5 * pow(x[i], 2) + 0.5 * pow(p[i], 2));
	}

	return v;
}

vd range(double start, double end, double step)
{
	vd v;
	for (int i = 0; i < (int)((end - start) / step); ++i)
		v.push_back(start + i * step);
	return v;
}

int main()
{

	Plot<double> plt, enPlt;

	plt.setXRange(-1.6, 1.6);
	plt.setYRange(-1.6, 1.6);

	enPlt.setYRange(0.95, 1.4);
	enPlt.setXRange(-0.5, 13.5);

	plt.addParametric("-sqrt(2) * sin(t - pi/4)", "sqrt(2) * cos(t - pi/4)", 0, 13, "exacta", "w l linewidth 2");
	enPlt.addParametric("t", "1", 0, 13, "energia real", "w l linewidth 2 lt rgb \"red\"");

	double h = 0.5;
	pair<vd, vd> v;
	for (int i = 1; i <= 5; ++i)
	{
		h /= 4;
		v = euler(1, 1, 13, h);
		vd vEnergy = energy(v.first, v.second), tRange = range(0, 13, h);

		enPlt.add(tRange, "", true);
		enPlt.add(vEnergy, "h = " + to_string(h), false, "w l linewidth 0.5");

		int nPoints = 100;
		v.first = takeNth(v.first, v.first.size() / nPoints);
		v.second = takeNth(v.second, v.second.size() / nPoints);

		plt.add(v.first, "", true);
		plt.add(v.second, "h = " + to_string(h), false, "w lp pointsize 0.5 linewidth 0.5");
	}

	plt.setLabels({"x", "p"});
	enPlt.setLabels({"t", "Energía Mecánica"});

	plt.plot("prueba");
	enPlt.plot("energy");
	return 0;
}