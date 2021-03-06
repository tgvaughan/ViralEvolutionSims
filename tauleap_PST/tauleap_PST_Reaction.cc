/*
 * tauleap_PST_Reaction.cc
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#include <vector>

#include <cstdlib>
#include <cmath>

#include "poissonian.h"
#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_Reaction.h"

/**
 * Define reaction in terms of number of reactants and products
 * together with the reaction rate and probability of mutation.
 */
Reaction::Reaction(int p_inX, int p_inY, int p_inYL, int p_inV,
		int p_outX, int p_outY, int p_outYL, int p_outV,
		bool p_mutY, bool p_mutYL, bool p_mutV,
		double p_rate, double p_mutrate) {

	inX = p_inX;
	inY = p_inY;
	inYL = p_inYL;
	inV = p_inV;
	outX = p_outX;
	outY = p_outY;
	outYL = p_outYL;
	outV = p_outV;

	mutY = p_mutY;
	mutYL = p_mutYL;
	mutV = p_mutV;

	rate = p_rate;
	mutrate = p_mutrate;

	dX = outX-inX;
	dY = (mutY ? -inY : outY-inY);
	dYL = (mutYL ? -inYL : outYL-inYL);
	dV = (mutV ? -inV : outV-inV);

	if (mutY || mutYL || mutV) {
		mutation = true;
	} else
		mutation = false;

	if ((inY == 0) && (inYL == 0) && (inV == 0) && (outY == 0) && (outYL == 0) && (outV == 0))
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
const double Reaction::get_gcond(int h1, int h2, int L)
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
double Reaction::getLeapDistance (double tau, double Ncrit, int critcond, const StateVec & sv, unsigned short *buf) {

	// Calculate X portion of propensity:
	aX = rate;
	for (int m=0; m<inX; m++)
		aX *= sv.X-m;

	if (onlyX) {

		// Check for critical X reaction:
		critX = false;

		switch(critcond) {
		case 0:
			break;

		case 1:
			if (sv.X + Ncrit*dX < 0)
				critX = true;
			break;

		case 2:
			double deltaX = tau*aX*dX;
			if ((deltaX<0) && (sv.X + deltaX - Ncrit*sqrt(-deltaX) < 0))
				critX = true;
			break;
		}

		// If critical, determine reaction time:
		if (critX) {
			double newtaucrit = -log(erand48(buf))/aX;
			if (newtaucrit<tau)
				return newtaucrit;
		}

		return tau;
	}

	double taucrit = tau;

	// Ensure propensity vectors are allocated:
	if (mutation) {
		amut.resize(3*(sv.maxHD+1), 0.0);
		critmut.resize(3*(sv.maxHD+1), false);
	} else {
		a.resize(sv.maxHD+1, 0.0);
		crit.resize(sv.maxHD+1, false);
	}

	for (int h=0; h<=sv.maxHD; h++) {

		double atmp = aX;
		for (int m=0; m<inY; m++)
			atmp *= sv.Y[h]-m;
		for (int m=0; m<inYL; m++)
			atmp *= sv.YL[h]-m;
		for (int m=0; m<inV; m++)
			atmp *= sv.V[h]-m;

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv.maxHD ? h+1 : sv.maxHD); hp++) {

				int idx = 3*h + hp - h + 1;

				amut[idx] = atmp*mutrate/(3.0*sv.L)*get_gcond(h,hp,sv.L);
				if (hp==h)
					amut[idx] += atmp*(1.0-mutrate);

				// Check for critical mutation reaction:
				critmut[idx] = false;

				switch(critcond) {
				case 0:
					break;

				case 1:
					if ((sv.Y[h] + Ncrit*dY < 0)
							|| (sv.YL[h] + Ncrit*dYL < 0)
							|| (sv.V[h] + Ncrit*dV < 0))
						critmut[idx] = true;
					break;

				case 2:
					double deltaY = tau*amut[idx]*dY;
					double deltaYL = tau*amut[idx]*dYL;
					double deltaV = tau*amut[idx]*dV;
					if (((deltaY <0) && (sv.Y[h] + deltaY - Ncrit*sqrt(-deltaY) < 0))
							|| ((deltaYL<0) && (sv.YL[h] + deltaYL - Ncrit*sqrt(-deltaYL) < 0))
							|| ((deltaV<0) && (sv.V[h] + deltaV - Ncrit*sqrt(-deltaV) < 0)))
						critmut[idx] = true;
					break;
				}

				// If critical, determine reaction time:
				if (critmut[idx]) {
					double newtaucrit = -log(erand48(buf))/amut[idx];
					if (newtaucrit < taucrit) {
						taucrit = newtaucrit;
						critreact = idx;
					}
				}
			}

		} else {

			a[h] = atmp;

			// Check for critical non-mutation reaction:
			crit[h] = false;
			switch(critcond) {
			case 0:
				break;
			case 1:
				if ((sv.Y[h] + Ncrit*dY < 0)
						|| (sv.YL[h] + Ncrit*dYL < 0)
						|| (sv.V[h] + Ncrit*dV < 0))
					crit[h] = true;
				break;

			case 2:
				double deltaY = tau*a[h]*dY;
				double deltaYL = tau*a[h]*dYL;
				double deltaV = tau*a[h]*dV;
				if (((deltaY<0) && (sv.Y[h] + deltaY - Ncrit*sqrt(-deltaY) < 0))
						|| ((deltaYL<0) && (sv.YL[h] + deltaYL - Ncrit*sqrt(-deltaYL) < 0))
						|| ((deltaV<0) && (sv.V[h] + deltaV - Ncrit*sqrt(-deltaV) < 0)))
					crit[h] = true;
				break;
			}

			// If critical, determine reaction time:
			if (crit[h]) {

				double newtaucrit = -log(erand48(buf))/a[h];
				if (newtaucrit < taucrit) {
					taucrit = newtaucrit;
					critreact = h;
				}
			}
		}
	}

	return taucrit;
}

/**
 * Implement "next" critical reaction.
 *
 * sv_new:		StateVec to update
 */
void Reaction::doCritical(StateVec & sv_new) {

	// Implement X component
	sv_new.X += outX - inX;

	// If X is all there is, we're done:
	if (onlyX)
		return;

	if (mutation) {

		// Implement Y and V components for mutation

		int h = critreact/3;
		int hp = h + (critreact%3) - 1;

		if (mutY) {
			sv_new.Y[h] += dY;
			sv_new.Y[hp] += outY;
			sv_new.YL[h] += dYL;
			sv_new.V[h] += dV;
		}
		if (mutYL) {
			sv_new.Y[h] += dY;
			sv_new.YL[h] += dYL;
			sv_new.YL[hp] += outYL;
			sv_new.V[h] += dV;
		}
		if (mutV) {
			sv_new.Y[h] += dY;
			sv_new.YL[h] += dYL;
			sv_new.V[h] += dV;
			sv_new.V[hp] += outV;
		}

	} else {

		// Implement Y and V components for non-mutation

		int h = critreact;
		sv_new.Y[h] += dY;
		sv_new.YL[h] += dYL;
		sv_new.V[h] += dV;

	}

	return;

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
			return false;

		// Implement reaction:
		double q = poissonian(dt*aX, buf);
		sv_new.X += q*(outX-inX);

		// Check for negative population:
		if (sv_new.X<0) {
			sv_new.X = 0;
			return true;
		}

		return false;
	}

	bool negativePop = false;

	for (int h=0; h<=sv_new.maxHD; h++) {

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv_new.maxHD ? h+1 : sv_new.maxHD); hp++) {
				int idx = 3*h + hp - h + 1;

				// Skip if reaction critical:
				if (critmut[idx])
					continue;

				// Implement reaction:
				double q = poissonian(dt*amut[idx], buf);

				sv_new.X += q*dX;
				sv_new.Y[h] += q*dY;
				sv_new.YL[h] += q*dYL;
				sv_new.V[h] += q*dV;

				if (mutY)
					sv_new.Y[hp] += q*outY;
				if (mutYL)
					sv_new.YL[hp] += q*outYL;
				if (mutV)
					sv_new.V[hp] += q*outV;

				// Check for negative populations:
				if (sv_new.X < 0) {
					sv_new.X = 0;
					negativePop = true;
				}
				if (sv_new.Y[h] < 0) {
					sv_new.Y[h] = 0;
					negativePop = true;
				}
				if (sv_new.YL[h] < 0) {
					sv_new.YL[h] = 0;
					negativePop = true;
				}
				if (sv_new.V[h] < 0) {
					sv_new.V[h] = 0;
					negativePop = true;
				}
			}

		} else {

			// Skip if reaction critical:
			if (crit[h])
				continue;

			// Implement reaction:
			double q = poissonian(dt*a[h], buf);

			sv_new.X += q*dX;
			sv_new.Y[h] += q*dY;
			sv_new.YL[h] += q*dYL;
			sv_new.V[h] += q*dV;

			// Check for negative populations:
			if (sv_new.X < 0) {
				sv_new.X = 0;
				negativePop = true;
			}
			if (sv_new.Y[h] < 0) {
				sv_new.Y[h] = 0;
				negativePop = true;
			}
			if (sv_new.YL[h] < 0) {
				sv_new.YL[h] = 0;
				negativePop = true;
			}
			if (sv_new.V[h] < 0) {
				sv_new.V[h] = 0;
				negativePop = true;
			}
		}

	}

	return negativePop;
}
