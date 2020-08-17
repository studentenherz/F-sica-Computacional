#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <vector>
#include <string>
// #include <pair>

using namespace std;

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

public:
	Plot() {}
	Plot(string _title) : title(_title) {}
	~Plot() {}

	void add(vector<Type> &col, string title = "", string style = "")
	{
		data.push_back(col);
		if (title == "")
			titles.push_back(to_string(data.size()));
		else
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

	void setAutorange()
	{
		autorange = true;
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

		for (int i = 0; i < data[0].size(); ++i)
		{
			for (auto col : data)
			{
				fData << col[i] << ' ';
			}
			fData << '\n';
		}

		fGnuScript << "set terminal svg\n"
							 << "set output \'" << fname + ".svg\'\n";

		if (autorange)
			fGnuScript << "set autoscale\n";
		else
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

		fGnuScript << "set key above\n";

		fGnuScript << "plot ";
		for (int i = 1; i < data.size(); ++i)
		{
			fGnuScript << "\'" << fname + ".dat\' using 1:" << i + 1 << " " << styles[i] << " title \'" << titles[i] << "\', ";
		}
	}
};

#endif