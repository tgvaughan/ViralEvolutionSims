/*
 * tauleap_PST_SemiImplicit.h
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_SEMIIMPLICIT_H_
#define TAULEAP_PST_SEMIIMPLICIT_H_

class SemiImplicit {
public:
	static double unif2poisson(double lambda, double u);
	static bool step(StateVec & sv, std::vector<Reaction> & reactions, double tau, unsigned short *buf);
};

#endif /* TAULEAP_PST_SEMIIMPLICIT_H_ */
