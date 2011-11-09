/*
 * tauleap_PST_SemiImplicit.cc
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#include <iostream>
#include <vector>

#include <boost/math/distributions.hpp>

#include "poissonian.h"
#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_Reaction.h"
#include "tauleap_PST_SemiImplicit.h"

/**
 * Transform uniform random variate into poissonian random variate using
 * the inverse CDF approach.
 *
 * lambda: mean of Poissonian distribution
 * u:      sample from uniform distribution on [0,1]
 *
 * returns: corresponding sample from Poissonian distribution
 */
double SemiImplicit::unif2poisson(double lambda, double u)
{
	return boost::math::quantile(boost::math::poisson(lambda), u);
}

/**
 * Increment state vector using Peter's semi-implicit algorithm.
 *
 * sv:         StateVec to increment
 * reactions:  array of reactions
 * Nreactions: number of reactions in array
 * tau:        size of time step
 * buf:        RNG buffer
 *
 * returns: true if negative populations were generated, false otherwise.
 */
bool SemiImplicit::step(StateVec & sv, Reaction *reactions, int Nreactions, double tau, unsigned short  *buf)
{
	bool negativePop = false;

	// TODO: Implement semi-implicit algorithm.


	return negativePop;
}
