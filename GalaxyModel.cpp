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
#include <gsl/gsl_sf_gamma.h>

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
	xc = Data::get_data().get_xMin() + Data::get_data().get_xRange();
	yc = Data::get_data().get_yMin() + Data::get_data().get_yRange();

	L1 = exp(log(1E-3) + log(1E6)*randomU());
	L2 = exp(log(1E-3) + log(1E6)*randomU());
	nu1 = exp(log(0.1) + log(30./0.1)*randomU());
	nu2 = exp(log(0.1) + log(30./0.1)*randomU());
	w = randomU();

	computeImage();
}

double GalaxyModel::perturb()
{
	int which = randInt(9);

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
		double scale = pow(10., 1.5 - 6.*randomU())*randn();
		xc += Data::get_data().get_xRange()*scale*randn();
		yc += Data::get_data().get_yRange()*scale*randn();
		xc = mod(xc - Data::get_data().get_xMin(),
			Data::get_data().get_xRange()) + Data::get_data().get_xMin();
		yc = mod(yc - Data::get_data().get_yMin(),
			Data::get_data().get_yRange()) + Data::get_data().get_yMin();

	}
	else if(which == 4)
	{
		L1 = log(L1);
		L1 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L1 = mod(L1 - log(1E-3), log(1E6)) + log(1E-3);
		L1 = exp(L1);
	}
	else if(which == 5)
	{
		L2 = log(L2);
		L2 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L2 = mod(L2 - log(1E-3), log(1E6)) + log(1E-3);
		L2 = exp(L2);
	}
	else if(which == 6)
	{
		nu1 = log(nu1);
		nu1 += log(30./0.1)*pow(10., 1.5 - 6.*randomU())*randn();
		nu1 = mod(nu1 - log(0.1), log(30./0.1)) + log(0.1);
		nu1 = exp(nu1);
	}
	else if(which == 7)
	{
		nu2 = log(nu2);
		nu2 += log(30./0.1)*pow(10., 1.5 - 6.*randomU())*randn();
		nu2 = mod(nu2 - log(0.1), log(30./0.1)) + log(0.1);
		nu2 = exp(nu2);
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

	double terms1 = gsl_sf_lngamma(0.5*(nu1 + 1.))
					- 0.5*log(nu1*M_PI) - gsl_sf_lngamma(0.5*nu1)
					- log(L1);
	double terms2 = gsl_sf_lngamma(0.5*(nu2 + 1.))
					- 0.5*log(nu2*M_PI) - gsl_sf_lngamma(0.5*nu2)
					- log(L2);

	double diff;
	for(int i=0; i<200; i++)
	{
		for(int j=0; j<200; j++)
		{
			diff = Data::get_data()(i, j) - image[i][j];
			if(diff > 0)
			{
				logL += log(w) + terms1 - 0.5*(nu1 + 1.)*
					(1. + pow(diff/L1, 2)/nu1);
			}
			else
			{
				diff *= -1.;
				logL += log(1. - w) + terms2 - 0.5*(nu2 + 1.)*
					(1. + pow(diff/L2, 2)/nu2); 
			}
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

			r = sqrt(pow(x - xc, 2) + pow(y - yc, 2));
			image[i][j] = rho/pow(1 + r/rc, gamma);
		}
	}
}

void GalaxyModel::print(std::ostream& out) const
{
	out<<rho<<' '<<rc<<' '<<gamma<<' '<<xc<<' '<<yc<<' '<<L1<<' '<<L2<<' '<<nu1<<' '<<nu2<<' '<<w;
}

string GalaxyModel::description() const
{
	return string("rho, rc, gamma, L1, L2, nu1, nu2, w");
}

