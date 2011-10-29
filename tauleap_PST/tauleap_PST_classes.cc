/*
 * tauleap_PST_classes.cc
 *
 *  Created on: 28/10/2011
 *      Author: Tim Vaughan
 */

#include <iostream>
#include <string>
#include <vector>

#include <cstdlib>
#include <cmath>

#include "poissonian.h"
#include "tauleap_PST_classes.h"


// StateVec member functions:

StateVec::StateVec (int length) {
	L = length;
	neighbourNum = 3*L;

	X = 0.0;
	Y.resize(L+1,0.0);
	V.resize(L+1,0.0);
}

StateVec::StateVec (const StateVec & src) {
	L = src.L;
	neighbourNum = src.neighbourNum;

	X = src.X;
	Y = src.Y;
	V = src.V;
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;
	neighbourNum = src.neighbourNum;

	X = src.X;
	Y = src.Y;
	V = src.V;

	return *this;
}


// Reaction member functions:

Reaction::Reaction(int p_inX, int p_inY, int p_inV,
		int p_outX, int p_outY, int p_outV,
		bool p_mutY, bool p_mutV,
		double p_rate, double p_mutrate) {

	inX = p_inX;
	inY = p_inY;
	inV = p_inV;
	outX = p_outX;
	outY = p_outY;
	outV = p_outV;

	mutY = p_mutY;
	mutV = p_mutV;

	rate = p_rate;
	mutrate = p_mutrate;

	if (mutY || mutV)
		mutation = true;
	else
		mutation = false;

	if ((inY == 0) && (inV == 0) && (outY == 0) && (outV == 0))
		onlyX = true;
	else
		onlyX = false;
}

Reaction::Reaction() { };

/**
 * Calculates the number of distinct mutations which begin at
 * a given sequence located in h1 and end with any sequence
 * located in h2.
 */
const double Reaction::get_gcond(int h2, int h1, int L)
{
	switch (h2-h1) {
		case -1:
			return h1;
		case 0:
			return 2*h1;
		case 1:
			return 3*(L-h1);
		default:
			return 0.0;
	}
}

/**
 * Calculate reaction propensities and return safe leap distance.
 *
 * tau:		desired leap distance
 * sv:		StateVec used to calculate propensities
 * buf:		RNG buffer
 *
 * returns:	safe leap distance
 */
double Reaction::getLeapDistance (double tau, double Ncrit, const StateVec & sv, unsigned short *buf) {

	// Calculate X portion of propensity:
	aX = rate;
	for (int m=0; m<inX; m++)
		aX *= sv.X-m;

	if (onlyX) {

		// Check for critical X reaction:
		double dX = Ncrit*(outX - inX);
		if ((dX<0) && (sv.X + dX < 0)) {

			critX = true;

//			std::cout << "M"; // DEBUG

			double newtaucrit = -log(erand48(buf))/aX;
			if (newtaucrit<tau)
				return newtaucrit;
		} else
			critX = false;

		return tau;
	}

	double taucrit = tau;

	// Ensure propensity vectors are allocated:
	if (mutation) {
		amut.resize(3*(sv.L+1), 0.0);
		critmut.resize(3*(sv.L+1), false);
	} else {
		a.resize(sv.L+1, 0.0);
		crit.resize(sv.L+1, false);
	}

	for (int h=0; h<=sv.L; h++) {

		double atmp = aX;
		for (int m=0; m<inY; m++)
			atmp *= sv.Y[h]-m;
		for (int m=0; m<inV; m++)
			atmp *= sv.V[h]-m;

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv.L ? h+1 : sv.L); hp++) {

				int idx = 3*h + hp - h + 1;

				amut[idx] = atmp*mutrate*get_gcond(hp,h,sv.L);
				if (hp==h)
					amut[idx] += atmp*(1-sv.neighbourNum*mutrate);

				// Check for critical mutation reaction:
				double dY = -Ncrit*inY;
				double dV = -Ncrit*inV;
				if (((dY<0) && (sv.Y[h] + dY < 0))
						|| ((dV<0) && (sv.V[h] + dV < 0))) {

					critmut[idx] = true;

					double newtaucrit = -log(erand48(buf))/amut[idx];
					if (newtaucrit < taucrit) {
						taucrit = newtaucrit;
						critreact = idx;
					}
				} else
					critmut[idx] = false;

			}

		} else {

			a[h] = atmp;

			// Check for critical non-mutation reaction:
			double dY = Ncrit*(outY-inY);
			double dV = Ncrit*(outV-inV);
			if (((dY<0) && (sv.Y[h] + dY < 0))
					|| ((dV<0) && (sv.V[h] + dV < 0))) {

				crit[h] = true;

//				std::cout << "M"; // DEBUG

				double newtaucrit = -log(erand48(buf))/a[h];
				if (newtaucrit < taucrit) {
					taucrit = newtaucrit;
					critreact = h;
				}
			} else
				crit[h] = false;
		}
	}

	return taucrit;
}

/**
 * Uses tau-leaping to implement non-critical reactions.
 *
 * tau:		length of leap
 * new_sv:	StateVec to modify
 * buf:		RNG buffer
 *
 * returns:	false if leap resulted in negative pop, true otherwise.
 */
bool Reaction::tauleap(double dt, StateVec & sv_new, unsigned short int *buf) {

	if (onlyX) {

		// Skip if reaction critical:
		if (critX)
			return true;

		// Implement reaction:
		double q = poissonian(dt*aX, buf);
		sv_new.X += q*(outX-inX);

		// Check for negative population:
		if (sv_new.X<0)
			return false;

		return true;
	}

	for (int h=0; h<=sv_new.L; h++) {

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv_new.L ? h+1 : sv_new.L); hp++) {
				int idx = 3*h + hp - h + 1;

				// Skip if reaction critical:
				if (critmut[idx])
					continue;

				// Implement reaction:
				double q = poissonian(dt*amut[idx], buf);
				sv_new.X += q*(outX - inX);
				if (mutY) {
					sv_new.Y[h] -= q*inY;
					sv_new.Y[hp] += q*outY;
					sv_new.V[h] += q*(outV - inV);
				}
				if (mutV) {
					sv_new.Y[h] += q*(outY - inY);
					sv_new.V[h] -= q*inV;
					sv_new.V[hp] += q*outV;
				}

				// Check for negative populations:
				if (sv_new.X < 0 || sv_new.Y[h] < 0 || sv_new.V[h] < 0)
					return false;
			}

		} else {

			// Skip if reaction critical:
			if (crit[h])
				continue;

			// Implement reaction:
			double q = poissonian(dt*a[h], buf);
			sv_new.X += q*(outX - inX);
			sv_new.Y[h] += q*(outY - inY);
			sv_new.V[h] += q*(outV - inV);

			// Check for negative populations:
			if (sv_new.X < 0 || sv_new.Y[h] < 0 || sv_new.V[h] < 0)
				return false;
		}

	}

	return true;
}

/**
 * Implement "next" critical reaction.
 *
 * sv_new:		StateVec to update
 * buf:			RNG buffer
 */
void Reaction::doCritical(StateVec & sv_new, unsigned short int *buf) {

	// Implement X component
	sv_new.X += outX - inX;

	// If X is all there is, we're done:
	if (critX)
		return;

	if (mutation) {

		// Implement Y and V components for mutation

		int h = critreact/3;
		int hp = h + (critreact%3) - 1;

		if (mutY) {
			sv_new.Y[h] -= inY;
			sv_new.Y[hp] += outY;
			sv_new.V[h] += outV-inV;
		}
		if (mutV) {
			sv_new.Y[h] += outY-inY;
			sv_new.V[h] -= inV;
			sv_new.V[hp] += outV;
		}

	} else {

		// Implement Y and V components for non-mutation

		int h = critreact;
		sv_new.Y[h] += outY-inY;
		sv_new.V[h] += outV-inV;
	}

	return;

}


// MomentScalar member functions

MomentScalar::MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name) {
	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;

	mean.resize(Nsamples);
	var.resize(Nsamples);
	sem.resize(Nsamples);
}
MomentScalar::MomentScalar() { };

/**
 * Sample moment
 */
void MomentScalar::sample(const StateVec & sv, int samp) {
	double tmp = (*samplefunc)(sv);
	mean[samp] += tmp;
	var[samp] += tmp*tmp;
}

/**
 * Perform post-processing of moment data
 */
void MomentScalar::normalise (int Npaths) {
	for (int s=0; s<Nsamples; s++) {
		mean[s] /= Npaths;
		var[s] = var[s]/Npaths - mean[s]*mean[s];
		sem[s] = sqrt(var[s]/Npaths);
	}
}


// MomentVector member functions:

MomentVector::MomentVector(int p_Nsamples, std::vector<double> (*p_samplefunc)(const StateVec &),
		std::string p_name, int p_sequenceL) {

	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;
	sequenceL = p_sequenceL;

	mean.resize(Nsamples*(sequenceL+1));
	var.resize(Nsamples*(sequenceL+1));
	sem.resize(Nsamples*(sequenceL+1));
}

MomentVector::MomentVector() { };

/**
 * Sample moment
 */
void MomentVector::sample(const StateVec & sv, int samp) {

	std::vector<double> tmp = (*samplefunc)(sv);

	for (int h=0; h<=sequenceL; h++) {
		int idx = (sequenceL+1)*samp + h;

		mean[idx] += tmp[h];
		var[idx] += tmp[h]*tmp[h];
	}
}

/**
 * Perform post-processing of moment data
 */
void MomentVector::normalise (int Npaths) {

	for (int samp=0; samp<Nsamples; samp++) {
		for (int h=0; h<=sequenceL; h++) {
			int idx = (sequenceL+1)*samp + h;

			mean[idx] /= Npaths;
			var[idx] = var[idx]/Npaths - mean[idx]*mean[idx];
			sem[idx] = sqrt(var[idx]/Npaths);
		}
	}

}
