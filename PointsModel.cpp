/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "PointsModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace DNest3;

PointsModel::PointsModel()
{

}

void PointsModel::fromPrior()
{
	xc = -10. + 20.*randomU();
	yc = -10. + 20.*randomU();
	sig = exp(log(1E-2) + log(1E3)*randomU());
}

double PointsModel::perturb()
{
	double logH = 0.;

	int which = randInt(3);
	if(which == 0)
	{
		xc += 20.*pow(10., 1.5 - 6.*randomU())*randn();
		xc = mod(xc + 10., 20.) - 10.;
	}
	else if(which == 1)
	{
		yc += 20.*pow(10., 1.5 - 6.*randomU())*randn();
		yc = mod(xc + 10., 20.) - 10.;
	}
	else
	{
		sig = log(sig);
		sig += log(1E3)*pow(10., 1.5 - 6.*randomU())*randn();
		sig = mod(sig - log(1E-2), log(1E3)) + log(1E-2);	
		sig = exp(sig);
	}

	return logH;
}


double PointsModel::logLikelihood() const
{
	// Load data
	static bool loaded = false;
	static vector<double> x, y;
	if(!loaded)
	{
		cout<<"# Loading data..."<<endl;
		fstream fin("point_data.txt", ios::in);
		if(!fin)
			cerr<<"# Can't load data."<<endl;
		if(fin.peek() == '#')
			fin.ignore(100000, '\n');
		double temp1, temp2;
		while(fin>>temp1 && fin>>temp2)
		{
			x.push_back(temp1);
			y.push_back(temp2);
		}
		cout<<"# Loaded "<<x.size()<<" data points."<<endl;
		fin.close();
		loaded = true;
	}

	double logL = -(x.size()*log(2.*M_PI*pow(sig, 2)));
	for(size_t i=0; i<x.size(); i++)
	{
		logL += -0.5*pow((x[i] - xc)/sig, 2) - 0.5*pow((y[i] - yc)/sig, 2);
	}
	return logL;
}

void PointsModel::print(std::ostream& out) const
{
	out<<xc<<' '<<yc<<' '<<sig;
}

string PointsModel::description() const
{
	return string("xc, yc, sig");
}

