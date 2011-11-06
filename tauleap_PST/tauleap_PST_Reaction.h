/*
 * tauleap_PST_Reaction.h
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_REACTION_H_
#define TAULEAP_PST_REACTION_H_

class Reaction {
public:

	// Variables describing reaction stoichiometry:
	int inX, inY, inV;
	int outX, outY, outV;

	// Whether or not new virions or infected cells have mutated genome:
	bool mutY, mutV;

	// Set to true for reactions which involve mutation:
	bool mutation;

	// Set to true for reactions which only involve X:
	bool onlyX; // this is dumb

	// Base reaction rate and mutation rate:
	double rate, mutrate;

	// Propensities:
	double aX;
	std::vector<double> a, amut;

	// Critical reaction flags:
	bool critX;
	std::vector<bool> crit, critmut;

	// Critical reaction identifier
	int critreact;

	Reaction(int inX, int inY, int inV, int outX, int outY, int outV, bool mutY,
			bool mutV, double rate, double mutrate);
	Reaction();

	const double get_gcond(int h2, int h1, int sequenceL);
	double getLeapDistance(double tau, double alpha, bool newcritcond,
			const StateVec & sv, unsigned short *buf);
	bool tauleap(double dt, StateVec & sv_new, unsigned short *buf);
	void doCritical(StateVec & sv_new);

};

#endif /* TAULEAP_PST_REACTION_H_ */
