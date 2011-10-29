/*
 * tauleap_PST_classes.h
 *
 *  Created on: 28/10/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_CLASSES_H_
#define TAULEAP_PST_CLASSES_H_


class StateVec {
    public:

        // Sequence length:
        int L;
		int neighbourNum;

        // Target cell population:
        double X;

        // Infected cell and virion populations:
        std::vector<double> Y, V;

        StateVec (int length);
        StateVec (const StateVec & sv);
        StateVec operator= (const StateVec & sv);

};

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
		std::vector<double> crit, critmut;

		// Critical reaction identifier
		int critreact;

		Reaction(int inX, int inY, int inV, int outX, int outY, int outV,
				bool mutY, bool mutV, double rate, double mutrate);
		Reaction();

		const double get_gcond(int h2, int h1, int sequenceL);
		double getLeapDistance(double tau, double alpha, const StateVec & sv, unsigned short *buf);
		bool tauleap(double dt, StateVec & sv_new, unsigned short *buf);
		void doCritical(StateVec & sv_new, unsigned short *buf);

};


class MomentVector {
	public:

		int Nsamples;
		std::vector<double> (*samplefunc)(const StateVec &);
		std::string name;

		std::vector<double> mean;
		std::vector<double> var;
		std::vector<double> sem;

		int sequenceL;

		MomentVector(int p_Nsamples, std::vector<double> (*p_samplefunc)(const StateVec &),
				std::string p_name, int p_sequenceL);
		MomentVector();

		void sample(const StateVec & sv, int samp);
		void normalise (int Npaths);

};

class MomentScalar {
	public:

		int Nsamples;
		double (*samplefunc)(const StateVec &);
		std::string name;

		std::vector<double> mean;
		std::vector<double> var;
		std::vector<double> sem;

		// Constructor:
		MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name);
		MomentScalar();

		void sample(const StateVec & sv, int samp);
		void normalise (int Npaths);

};

#endif /* TAULEAP_PST_CLASSES_H_ */
