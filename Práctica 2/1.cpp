#include "util.h"
#include <functional>
#include <time.h>

// segunda derivada centrada y preriódica
vector<vd> d2_CP(vector<vd> v, double dx)
{
	int n = v.size();
	vector<vd> dv;
	for (int i = 0; i < n; ++i)
	{
		dv.push_back((1 / pow(dx, 2)) * (v[(n + i + 1) % n] - 2.0 * v[(n + i) % n] + v[(n + i - 1) % n]));
	}
	return dv;
}

// primera derivada centrada y periódica o(h^4)
vector<double> d_h4(vector<double> v, double dx)
{
	vector<double> d1;
	int n = v.size();
	for (int i = 0; i < n; ++i)
		d1.push_back((-v[(n + i + 2) % n] + 8 * v[(n + i + 1) % n] - 8 * v[(n + i - 1) % n] + v[(n + i - 2) % n]) / (12 * dx));
	return d1;
}

template <class Type>
vector<Type> operator*(vector<Type> a, vector<Type> b)
{
	vector<Type> p;
	if (a.size() != b.size())
		throw "Not matching size for scalar product";
	for (int i = 0; i < (int)a.size(); ++i)
		p.push_back(a[i] * b[i]);
	return p;
}

class EqSystem
{
	vector<function<double(vd)>> f;

public:
	EqSystem() {}
	~EqSystem()
	{
		f.clear();
	}

	void add(double _f(vd))
	{
		f.push_back(_f);
	}

	vd operator()(vd v)
	{
		vd ans;
		for (int i = 0; i < (int)f.size(); ++i)
		{
			ans.push_back(f[i](v));
		}
		return ans;
	}

	vector<vd> operator()(vector<vd> v)
	{
		vector<vd> ans;
		for (int i = 0; i < (int)v.size(); ++i)
		{
			ans.push_back(this->operator()(v[i]));
		}
		return ans;
	}
};

pair<double, vector<vd>> FTCS(vector<vd> v0, vd D, EqSystem f, double t0, double T, double dt, double dx)
{
	vector<vd> sol = v0;
	double t = t0;

	while (t < T)
	{
		vector<vd> d2_v = d2_CP(sol, dx); // segunda derivada espacial

		// trow exceptions if inf or nan
		for (int i = 0; i < (int)d2_v.size(); ++i)
			for (int j = 0; j < (int)d2_v[0].size(); ++j)
				if (isnan(d2_v[i][j]) || isinf(d2_v[i][j]))
				{
					string s = "Found nan or inf at t = " + to_string(t);
					throw s;
				}

		sol = sol + dt * (D * d2_v + f(sol));
		t += dt; // un paso de tiempo
	}
	return {t, sol};
}

double rand(double a, double b)
{
	return a + (b - a) * rand() / RAND_MAX;
}

double period(vector<double> v0, double dx)
{
	vd v = v0;
	v.insert(v.end(), v0.begin(), v0.end()); // repetir el vector por si hay un solo máximo

	vd d1 = d_h4(v, dx);
	vd d2 = d_h4(d1, dx);

	vector<int> extremas;
	for (int i = 0; i < (int)v.size(); ++i)
	{
		if (d1[i] * d1[i + 1] < 0)
		{
			extremas.push_back(i);
		}
	}

	vector<int> maximas;
	for (int i : extremas)
	{
		if (d2[i] > 0)
			maximas.push_back(i);
	}

	double x = 0;
	for (int i = 1; i < (int)maximas.size(); ++i)
	{
		x += dx * (maximas[i] - maximas[i - 1]);
	}

	return x / (maximas.size() - 1);
}

int main()
{
	EqSystem g;
	g.add([](vd v) -> double { return (-7 * pow(v[0], 2) - 50 * v[0] * v[1] + 57) / 32.0; });
	g.add([](vd v) -> double { return (7 * pow(v[0], 2) + 50 * v[0] * v[1] - 2 * v[1] - 55) / 32.0; });

	double L = 1, dt = 0.001, dx = 0.01;
	vector<double> x;
	vd v_ast = {1, 1};
	vector<vd> v0;
	double a = 0.15;

	srand(time(0));

	for (int i = 0; i < (int)(L / dx); ++i)
	{
		v0.push_back(v_ast + a * vd({rand(-1, 1), rand(-1, 1)}));
		x.push_back(dx * i);
	}

	vector<vd> v;
	vector<double> T = {0, 5, 20, 50, 200, 300};
	vd Ds = {0.01, 0.0025, 0.0015, 0.00075, 0.0005};

	for (double D : Ds)
	{
		v = v0;
		double t = 0;
		for (double tau : T)
		{
			try
			{
				pair<double, vector<vd>> ans = FTCS(v, {D, D / 2}, g, t, tau, dt, dx);
				v = ans.second;
				t = ans.first;

				vector<vector<double>> v1 = transpose(v);
				Plot<double> plt("D =" + to_string(D) + ", t = " + to_string(tau).substr(0, 3));
				plt.setLabels({"x", "u, v"});
				plt.add(x, "", true);
				double l1 = period(v1[0], dx);
				double l2 = period(v1[1], dx);
				plt.add(v1[0], "u(x)", false, "w l lc 1");
				plt.add(v1[1], "v(x)", false, "w l lc 2");
				plt.plot("D=" + to_string(D) + "t=" + to_string(tau).substr(0, 3));

				cout << "lambda 1 = " << period(v1[0], dx) << "\nlambda 2 = " << period(v1[1], dx) << "\n============================\n";
			}
			catch (string e)
			{
				cout << e << '\n';
			}
		}
	}

	cout << "ok";
	return 0;
}
