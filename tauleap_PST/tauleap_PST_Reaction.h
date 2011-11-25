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
	int inX, inY, inYL, inV;
	int outX, outY, outYL, outV;

	// Whether or not new virions or infected cells have mutated genome:
	bool mutY, mutYL, mutV;

	// Reactant changes taking into account mutation:
	int dX, dY, dYL, dV;

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

	Reaction(int inX, int inY, int inYL, int inV,
			int outX, int outY, int outYL, int outV,
			bool mutY, bool mutYL, bool mutV, double rate, double mutrate);
	Reaction();

	const double get_gcond(int h2, int h1, int sequenceL);
	double getLeapDistance(double tau, double alpha, int critcond,
			const StateVec & sv, unsigned short *buf);
	void doCritical(StateVec & sv_new);
	bool tauleap(double dt, StateVec & sv_new, unsigned short *buf);
};

#endif /* TAULEAP_PST_REACTION_H_ */
