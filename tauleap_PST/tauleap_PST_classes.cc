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

StateVec::StateVec (const StateVec & sv) {
	L = sv.L;
	neighbourNum = sv.neighbourNum;

	X = sv.X;
	Y = sv.Y;
	V = sv.V;
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
bool Reaction::leap (double tau, StateVec & sv, StateVec & sv_new, unsigned short *buf) {

	for (int h=0; h<=sv.L; h++) {

		double a = rate;
		for (int m=0; m<inX; m++)
			a *= sv.X-m;
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

				sv_new.X += q*(outX-inX);
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
