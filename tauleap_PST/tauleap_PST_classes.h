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
        StateVec (const StateVec &);
        StateVec operator= (const StateVec &);

};

class Reaction {
	public:

		int inX, inY, inV;
		int outX, outY, outV;

		bool mutY, mutV;

		bool mutation;

		double rate, mutrate;

		Reaction(int,int,int, int,int,int, bool,bool, double,double);
		Reaction();

		double get_gcond(int, int, int);
		bool leap(double, StateVec &, StateVec &, unsigned short *);

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

		// Constructor:
		MomentVector(int p_Nsamples, std::vector<double> (*p_samplefunc)(const StateVec &),
				std::string p_name, int p_sequenceL) {

			Nsamples = p_Nsamples;
			samplefunc = p_samplefunc;
			name = p_name;
			sequenceL = p_sequenceL;

			mean.resize(Nsamples*(sequenceL+1));
			var.resize(Nsamples*(sequenceL+1));
			sem.resize(Nsamples*(sequenceL+1));
		}
		MomentVector() { };

		// Sample moment:
		void sample(const StateVec & sv, int samp) {

			std::vector<double> tmp = (*samplefunc)(sv);

			for (int h=0; h<=sequenceL; h++) {
				int idx = (sequenceL+1)*samp + h;

				mean[idx] += tmp[h];
				var[idx] += tmp[h]*tmp[h];
			}
		}

		// Post processing of moment data:
		void normalise (int Npaths) {

			for (int samp=0; samp<Nsamples; samp++) {
				for (int h=0; h<=sequenceL; h++) {
					int idx = (sequenceL+1)*samp + h;

					mean[idx] /= Npaths;
					var[idx] = var[idx]/Npaths - mean[idx]*mean[idx];
					sem[idx] = sqrt(var[idx]/Npaths);
				}
			}

		}
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
		MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name) {
			Nsamples = p_Nsamples;
			samplefunc = p_samplefunc;
			name = p_name;

			mean.resize(Nsamples);
			var.resize(Nsamples);
			sem.resize(Nsamples);
		}
		MomentScalar() { };

		// Sample moment:
		void sample(const StateVec & sv, int samp) {
			double tmp = (*samplefunc)(sv);
			mean[samp] += tmp;
			var[samp] += tmp*tmp;
		}

		// Post processing of moment data:
		void normalise (int Npaths) {
			for (int s=0; s<Nsamples; s++) {
				mean[s] /= Npaths;
				var[s] = var[s]/Npaths - mean[s]*mean[s];
				sem[s] = sqrt(var[s]/Npaths);
			}
		}
};

#endif /* TAULEAP_PST_CLASSES_H_ */
