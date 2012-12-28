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
	rc = exp(log(1E-2) + log(10./1E-2)*randomU());
	gamma = 6*randomU();
	xc = Data::get_data().get_xMin() + Data::get_data().get_xRange();
	yc = Data::get_data().get_yMin() + Data::get_data().get_yRange();
	q = exp(randn());
	theta = 2.*M_PI*randomU(); cosTheta = cos(theta); sinTheta = sin(theta);

	sig0 = exp(log(1E-3) + log(1E6)*randomU());
	sig1 = exp(log(1E-3) + log(1E6)*randomU());

	beta = 0.1 + 1.9*randomU();
	L = exp(log(1E-3) + log(1E6)*randomU());
	w = randomU();

	computeImage();

}

double GalaxyModel::perturb()
{
	int which = randInt(10);
	double logH = 0.;

	if(which == 0)
	{
		double old = rho;

		rho = log(rho);
		rho += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		rho = mod(rho - log(1E-3), log(1E6)) + log(1E-3);
		rho = exp(rho);

		double ratio = rho/old;
		for(size_t i=0; i<image.size(); i++)
			for(size_t j=0; j<image[i].size(); j++)
				image[i][j] *= ratio;
		
	}
	else if(which == 1)
	{
		rc = log(rc);
		rc += log(10./1E-2)*pow(10., 1.5 - 6.*randomU())*randn();
		rc = mod(rc + log(1E-2), log(10./1E-2)) + log(1E-2);
		rc = exp(rc);
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
		sig0 = log(sig0);
		sig0 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		sig0 = mod(sig0 - log(1E-3), log(1E6)) + log(1E-3);
		sig0 = exp(sig0);

		sig1 = log(sig1);
		sig1 += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		sig1 = mod(sig1 - log(1E-3), log(1E6)) + log(1E-3);
		sig1 = exp(sig1);
	}
	else if(which == 7)
	{
		beta += 1.9*pow(10., 1.5 - 6.*randomU())*randn();
		beta = mod(beta - 0.1, 1.9) + 0.1;
	}
	else if(which == 8)
	{
		L = log(L);
		L += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		L = mod(L - log(1E-3), log(1E6)) + log(1E-3);
		L = exp(L);
	}
	else if(which == 9)
	{
		w += pow(10., 1.5 - 6.*randomU())*randn();
		w = mod(w, 1.);
	}

	return logH;
}


double GalaxyModel::logLikelihood() const
{

	double logL = 0.;

	double term1, term2, diff, var;
	for(size_t i=0; i<image.size(); i++)
	{
		for(size_t j=0; j<image[i].size(); j++)
		{
			if(Data::get_data()(i, j) > -1E250)
			{
				diff = Data::get_data()(i, j) - image[i][j];

				var = sig0*sig0 + sig1*image[i][j];
				term1 = log(w) - 0.5*log(2.*M_PI*var)
					-0.5*pow(diff, 2)/var;

				if(diff < 0)
				{
					logL += term1;
				}
				else
				{
					term2 = log(1. - w) - log(L) - diff/L;
					logL += logsumexp(term1, term2);
				}

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
			image[i][j] = rho*exp(-pow(r/rc, gamma));
		}
	}
}

void GalaxyModel::print(std::ostream& out) const
{
	out<<rho<<' '<<rc<<' '<<gamma<<' '<<xc<<' '<<yc<<' '<<q<<' '<<theta<<
			' '<<sig0<<' '<<sig1<<' '<<beta<<' '<<L<<' '<<w;
}

string GalaxyModel::description() const
{
	return string("rho, rc, gamma, xc, yc, q, theta, sig, beta, L, w");
}

