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

#include "GalaxyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest3;

GalaxyModel::GalaxyModel()
:image(200, vector<double>(200))
{

}

void GalaxyModel::fromPrior()
{
	// Generates initial values from the prior
	rho = exp(log(1E-3) + log(1E6)*randomU());
	rc = randomU();
	gamma = 6*randomU();

	L1 = exp(log(1E-3) + log(1E6)*randomU());
	L2 = exp(log(1E-3) + log(1E6)*randomU());
	w = randomU();

	computeImage();
}

double GalaxyModel::perturb()
{
	int which = randInt(6);

	if(which == 0)
	{
		rho = log(rho);
		rho += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		rho = mod(rho - log(1E-3), log(1E6)) + log(1E-3);
		rho = exp(rho);
		computeImage();
	}
	else if(which == 1)
	{
		rc += pow(10., 1.5 - 6.*randomU())*randn();
		rc = mod(rc, 1.);
		computeImage();
	}
	else if(which == 2)
	{
		gamma += 6.*pow(10., 1.5 - 6.*randomU())*randn();
		gamma = mod(gamma, 6.);
		computeImage();
	}
	else if(which == 3)
	{
		L1 = log(L1);
		L1 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L1 = mod(L1 - log(1E-3), log(1E6)) + log(1E-3);
		L1 = exp(L1);
	}
	else if(which == 4)
	{
		L2 = log(L2);
		L2 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L2 = mod(L2 - log(1E-3), log(1E6)) + log(1E-3);
		L2 = exp(L2);
	}
	else
	{
		w += pow(10., 1.5 - 6.*randomU())*randn();
		w = mod(w, 1.);
	}

	return 0.;
}


double GalaxyModel::logLikelihood() const
{

	double logL = 0.;

	double diff;
	for(int i=0; i<200; i++)
	{
		for(int j=0; j<200; j++)
		{
			diff = Data::get_data()(i, j) - image[i][j];
			if(diff > 0)
				logL += log(w) - log(L2) - diff/L2;
			else
				logL += log(1.-w) - log(L1) + diff/L1; 
		}
	}

	return logL;
}

void GalaxyModel::computeImage()
{
	double dx = 2./200;
	double dy = 2./200;

	double x, y, r;

	for(int i=0; i<200; i++)
	{
		y = 1. - (i + 0.5)*dy;
		for(int j=0; j<200; j++)
		{
			x = -1. + (j + 0.5)*dx;

			r = sqrt(x*x + y*y);
			image[i][j] = rho/pow(1 + r/rc, gamma);

		}
	}
}

void GalaxyModel::print(std::ostream& out) const
{
	out<<rho<<' '<<rc<<' '<<gamma<<' '<<L1<<' '<<L2<<' '<<w;
}

string GalaxyModel::description() const
{
	return string("rho, rc, gamma");
}

