/*
 * tauleap_PST_classes.cc
 *
 *  Created on: 28/10/2011
 *      Author: Tim Vaughan
 */

#include <string>
#include <vector>

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

	Y.resize(L+1,0.0);
	V.resize(L+1,0.0);

	for (int h=0; h<=L; h++) {
		Y[h] = src.Y[h];
		V[h] = src.V[h];
	}
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;
	neighbourNum = src.neighbourNum;

	X = src.X;

	for (int h=0; h<=L; h++) {
		Y[h] = src.Y[h];
		V[h] = src.V[h];
	}

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
		scalar = true;
	else
		scalar = false;
}

Reaction::Reaction() { };

/**
 * Calculates the number of distinct mutations which begin at
 * a given sequence located in h1 and end with any sequence
 * located in h2.
 */
double Reaction::get_gcond(int h2, int h1, int L)
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
 * Implements the reaction as a tau-leap.
 *
 * tau: 	size of leap
 * sv:		StateVec used to calculate propensities
 * sv_new:	StateVec modified as a result of the reaction(s)
 * buf:		RNG buffer
 *
 * Returns:	true or false depending on the
 */
bool Reaction::leap (double tau, const StateVec & sv, StateVec & sv_new,
		unsigned short *buf) {

	double aX = rate;
	for (int m=0; m<inX; m++)
		aX *= sv.X-m;

	if (scalar) {
		int q = poissonian(aX*tau, buf);
		sv_new.X += q*(outX - inX);

		if (sv_new.X<0)
			return false;

		return true;
	}

	for (int h=0; h<=sv.L; h++) {

		double a = aX;
		for (int m=0; m<inY; m++)
			a *= sv.Y[h]-m;
		for (int m=0; m<inV; m++)
			a *= sv.V[h]-m;

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv.L ? h+1 : sv.L); hp++) {

				double ap = a*mutrate*get_gcond(hp,h,sv.L);
				if (hp==h)
					ap += a*(1-sv.neighbourNum*mutrate);

				int q = poissonian(ap*tau, buf);

				sv_new.Y[h] -= q*inY;
				if (mutY)
					sv_new.Y[hp] += q*outY;
				else
					sv_new.Y[h] += q*outY;

				sv_new.V[h] -= q*inV;
				if (mutV)
					sv_new.V[hp] += q*outV;
				else
					sv_new.V[h] += q*outV;
			}

		} else {

			int q = poissonian(a*tau, buf);

			sv_new.X += q*(outX-inX);
			sv_new.Y[h] += q*(outY-inY);
			sv_new.V[h] += q*(outV-inV);

			if (sv_new.X<0 || sv_new.Y[h]<0 || sv_new.V[h]<0)
				return false;
		}

	}

	return true;
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
