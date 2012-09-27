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
:image(Data::get_data().get_ni(), vector<double>(Data::get_data().get_nj()))
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
	q = exp(randn());
	theta = 2.*M_PI*randomU(); cosTheta = cos(theta); sinTheta = sin(theta);

	L1 = exp(log(1E-3) + log(1E6)*randomU());
	L2 = exp(log(1E-3) + log(1E6)*randomU());
	nu1 = exp(log(0.1) + log(30./0.1)*randomU());
	nu2 = exp(log(0.1) + log(30./0.1)*randomU());
	w = randomU();

	computeImage();
}

double GalaxyModel::perturb()
{
	int which = randInt(11);
	double logH = 0.;

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
		computeImage();

	}
	else if(which == 4)
	{
		q = log(q);
		logH -= -0.5*q*q;
		q += pow(10., 1.5 - 6.*randomU())*randn();
		logH += -0.5*q*q;
		q = exp(q);
		computeImage();
	}
	else if(which == 5)
	{
		theta += 2.*M_PI*pow(10., 1.5 - 6.*randomU())*randn();
		theta = mod(theta, 2.*M_PI); cosTheta = cos(theta); sinTheta = sin(theta);
		computeImage();
	}
	else if(which == 6)
	{
		L1 = log(L1);
		L1 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L1 = mod(L1 - log(1E-3), log(1E6)) + log(1E-3);
		L1 = exp(L1);
	}
	else if(which == 7)
	{
		L2 = log(L2);
		L2 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L2 = mod(L2 - log(1E-3), log(1E6)) + log(1E-3);
		L2 = exp(L2);
	}
	else if(which == 8)
	{
		nu1 = log(nu1);
		nu1 += log(30./0.1)*pow(10., 1.5 - 6.*randomU())*randn();
		nu1 = mod(nu1 - log(0.1), log(30./0.1)) + log(0.1);
		nu1 = exp(nu1);
	}
	else if(which == 9)
	{
		nu2 = log(nu2);
		nu2 += log(30./0.1)*pow(10., 1.5 - 6.*randomU())*randn();
		nu2 = mod(nu2 - log(0.1), log(30./0.1)) + log(0.1);
		nu2 = exp(nu2);
	}
	else if(which == 10)
	{
		w += pow(10., 1.5 - 6.*randomU())*randn();
		w = mod(w, 1.);
	}

	return logH;
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
	for(size_t i=0; i<image.size(); i++)
	{
		for(size_t j=0; j<image[i].size(); j++)
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
	double dx = Data::get_data().get_dx();
	double dy = Data::get_data().get_dy();

	double x, y, r, xx, yy;

	for(size_t i=0; i<image.size(); i++)
	{
		y = 1. - (i + 0.5)*dy;
		for(size_t j=0; j<image[i].size(); j++)
		{
			x = -1. + (j + 0.5)*dx;

			xx =  cosTheta*(x - xc) + sinTheta*(y - yc);
			yy = -sinTheta*(x - xc) + cosTheta*(y - yc);
			r = sqrt(q*pow(xx, 2) + pow(yy, 2)/q);
			image[i][j] = rho/pow(1 + r/rc, gamma);
		}
	}
}

void GalaxyModel::print(std::ostream& out) const
{
	out<<rho<<' '<<rc<<' '<<gamma<<' '<<xc<<' '<<yc<<' '<<q<<' '<<theta<<
			' '<<L1<<' '<<L2<<' '<<nu1<<' '<<nu2<<' '<<w;
}

string GalaxyModel::description() const
{
	return string("rho, rc, gamma, xc, yc, q, theta, L1, L2, nu1, nu2, w");
}

