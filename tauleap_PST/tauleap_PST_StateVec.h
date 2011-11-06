/*
 * tauleap_PST_StateVec.h
 *
 *  Created on: 06/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_STATEVEC_H_
#define TAULEAP_PST_STATEVEC_H_

class StateVec {
public:

	// Sequence length:
	int L;

	// Target cell population:
	double X;

	// Infected cell and virion populations:
	std::vector<double> Y, V;

	StateVec(int length);
	StateVec(const StateVec & sv);
	StateVec operator=(const StateVec & sv);

};

#endif /* TAULEAP_PST_STATEVEC_H_ */
