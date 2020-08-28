#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

// Habrá muchos vector<double> mejor ahorramos espacio y esfuerzo con una abreviación
typedef vector<double> vd;

// Para sumar vectores miembro a miembro
template <class Type>
vector<Type> operator+(vector<Type> a, vector<Type> b)
{
	vector<Type> v(a);
	for (int i = 0; i < (int)a.size(); ++i)
		v[i] = v[i] + b[i];
	return v;
}

// Para restar vectores miembro a miembro
template <class Type>
vector<Type> operator-(vector<Type> a, vector<Type> b)
{
	for (int i = 0; i < (int)a.size(); ++i)
		a[i] = a[i] - b[i];
	return a;
}

// Multiplicación escalar * vector
template <class Type>
vector<Type> operator*(Type k, vector<Type> a)
{
	for (int i = 0; i < (int)a.size(); ++i)
		a[i] = k * a[i];
	return a;
}

// Multiplicación vector * escalar
template <class Type>
vector<Type> operator*(vector<Type> a, Type k)
{
	return k * a;
}

// // Multiplicación vector * vector = vector, término a término
// template <class Type>
// vector<Type> operator*(vector<Type> a, vector<Type> b)
// {
// 	for (int i = 0; i < (int)a.size(); ++i)
// 	{
// 		a[i] = a[i] * b[i];
// 	}
// 	return a;
// }

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

template <class Type>
vector<vector<Type>> operator*(vector<vector<Type>> v, Type a)
{
	vector<vector<Type>> ans;
	for (auto x : v)
	{
		ans.push_back(a * x);
	}
	return ans;
}

template <class Type>
vector<vector<Type>> operator*(Type a, vector<vector<Type>> v)
{
	vector<vector<Type>> ans;
	for (auto x : v)
	{
		ans.push_back(a * x);
	}
	return ans;
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

// Take every n-th element from vector
template <class Type>
vector<Type> takeNth(vector<Type> v, int n, int start = 0)
{
	vector<Type> u;
	for (int i = start; i < (int)v.size(); i += n)
	{
		u.push_back(v[i]);
	}
	return u;
}

template <class Type>
class Plot
{
	vector<vector<Type>> data;
	vector<string> titles;
	vector<string> styles;
	string title = "";
	string xLabel = "", yLabel = "";
	bool logScaleX = false, logScaleY = false;
	Type minX, maxX, minY, maxY;
	bool autorange = true;
	vector<bool> x;
	vector<string> param;
	vector<string> paramRange;
	vector<string> paramTitles;
	vector<string> paramStyle;

public:
	Plot() {}
	Plot(string _title) : title(_title) {}
	~Plot() {}

	void addParametric(string fx, string fy, double tmin, double tmax, string title = "", string style = "w l")
	{
		param.push_back("(" + fx + "):(" + fy + ")");
		paramRange.push_back("[t=" + to_string(tmin) + ":" + to_string(tmax) + "]");
		paramTitles.push_back(title);
		paramStyle.push_back(style);
	}

	void add(vector<Type> &col, string title = "", bool isX = false, string style = "")
	{
		if (data.empty())
		{
			x.push_back(true);
		}
		else
		{
			x.push_back(isX);
		}
		data.push_back(col);

		titles.push_back(title);
		if (style == "")
			styles.push_back("w lp pointtype 7 pointsize 0.5");
		else
			styles.push_back(style);
	}

	void plot(string fname)
	{
		_make_files(fname);
		string _temp = "gnuplot " + fname + ".gnuplot-script";
		const char *s = _temp.c_str();
		system(s);

		cout << "[*] Plotted " << fname << ".svg\n";
	}

	void setLabels(pair<string, string> labels)
	{
		xLabel = labels.first;
		yLabel = labels.second;
	}

	void setLogScale(pair<bool, bool> logscale)
	{
		logScaleX = logscale.first;
		logScaleY = logscale.second;
	}

	void setXRange(Type xmin, Type xmax)
	{
		autorange = false;
		minX = xmin;
		maxX = xmax;
	}

	void setYRange(Type ymin, Type ymax)
	{
		autorange = false;
		minY = ymin;
		maxY = ymax;
	}

private:
	void _make_files(string fname)
	{
		ofstream fData, fGnuScript;
		fData.open(fname + ".dat", ofstream::out);
		fGnuScript.open(fname + ".gnuplot-script", ofstream::out);

		if (!data.empty())
		{
			int max_size = 0;
			for (auto col : data)
			{
				max_size = max(max_size, (int)col.size());
			}

			for (int i = 0; i < max_size; ++i)
			{
				for (auto col : data)
				{
					if (i < col.size())
						fData << col[i] << ' ';
					else
						fData << ". ";
				}
				fData << '\n';
			}
		}

		fGnuScript << "set encoding utf8\n";

		fGnuScript << "set terminal svg\n"
							 << "set size square\n"
							 << "set output \'" << fname + ".svg\'\n";

		if (!(logScaleX || logScaleY) && autorange)
		{
			bool cmpredx = false, cmpredy = false;
			for (int i = 0; i < (int)data.size(); ++i)
			{
				if (x[i])
				{
					if (!cmpredx)
					{
						minX = maxX = data[i][0];
						cmpredx = true;
					}

					for (Type val : data[i])
					{
						minX = min(minX, val);
						maxX = max(maxX, val);
					}
				}
				else
				{
					if (!cmpredy)
					{
						minY = maxY = data[i][0];
						cmpredy = true;
					}

					for (Type val : data[i])
					{
						minY = min(minY, val);
						maxY = max(maxY, val);
					}
				}
			}

			// make plot range 10% bigger than actual range
			Type lx = (maxX - minX);
			Type ly = (maxY - minY);

			minX -= 0.05 * lx;
			maxX += 0.05 * lx;
			minY -= 0.05 * ly;
			maxY += 0.05 * ly;
		}

		if (!(logScaleX || logScaleY))
		{
			fGnuScript << "set xrange [" << minX << ":" << maxX << "]\n";
			fGnuScript << "set yrange [" << minY << ":" << maxY << "]\n";
		}

		if (title != "")
			fGnuScript << "set title \'" << title << "\'\n";

		if (xLabel != "")
			fGnuScript << "set xlabel \'" << xLabel << "\'\n";

		if (xLabel != "")
			fGnuScript << "set ylabel \'" << yLabel << "\'\n";

		if (logScaleX)
			fGnuScript << "set logscale x\n";

		if (logScaleY)
			fGnuScript << "set logscale y\n";

		fGnuScript << "set key outside right center\n";

		if (!data.empty())
		{
			fGnuScript << "plot ";
			int lastX = 0;
			for (int i = 1; i < data.size(); ++i)
			{
				if (x[i])
				{
					lastX = i;
					continue;
				}
				fGnuScript << "\'" << fname + ".dat\' using " << lastX + 1 << ":" << i + 1 << " " << styles[i];
				if (titles[i] != "")
					fGnuScript << " title \'" << titles[i] << "\'";
				else
					fGnuScript << " notitle";
				fGnuScript << ", ";
			}
		}

		if (!param.empty())
		{
			for (int i = 0; i < (int)param.size(); ++i)
			{
				fGnuScript << paramRange[i] << " \'+\' using " << param[i] << " title \'" << paramTitles[i] << "\' " << paramStyle[i] << " , ";
			}
		}
	}
};

#endif