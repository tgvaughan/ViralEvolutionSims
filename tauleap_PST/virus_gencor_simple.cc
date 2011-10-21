// Within-host viral evolution simulation.  Makes use of a projected
// sequence space of dramatically reduced dimension.

#include <iostream>
#include <vector>

#include <cmath>
#include <cstdlib>

#include "poissonian.h"

class StateVec {
    public:

        // Sequence length:
        int L;
		int neighbourNum;

        // Target cell population:
        double X;

        // Infected cell and virion populations:
        std::vector<double> Y, V;

        // Constructor:
        StateVec (int length) {
            L = length;
			neighbourNum = 3*L;

            X = 0.0;
            Y.resize(L+1,0.0);
            V.resize(L+1,0.0);
        }

        // Copy constructor:
        StateVec (const StateVec & sv) {
            L = sv.L;
            Y = sv.Y;
            V = sv.V;
        }

};

// Number of neighbouring sequences of a sequence on
// h1 which lie on h2:
double gcond(int h2, int h1, int L)
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

// Reaction class
class Reaction {
	public:

		int inX, inY, inV;
		int outX, outY, outV;

		bool mutY, mutV;

		bool mutation;

		double rate, mutrate;

		// Constructor:
		Reaction(int p_inX, int p_inY, int p_inV,
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
		Reaction () { };

		// Implement tau-leap:
		void leap (double tau, StateVec & sv0, StateVec & sv_res, unsigned short *buf) {

			for (int h=0; h<sv0.L; h++) {

				double a = rate;
				for (int m=0; m<inX; m++)
					a *= sv0.X-m;
				for (int m=0; m<inY; m++)
					a *= sv0.Y[h]-m;
				for (int m=0; m<inV; m++)
					a *= sv0.V[h]-m;

				if (mutation) {

					for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv0.L ? h+1 : sv0.L); hp++) {

						double ap = a*mutrate*gcond(hp,h,sv0.L);
						if (hp==h)
							ap += a*(1-sv0.neighbourNum*mutrate);

						int q = poissonian(ap*tau, buf);

						sv_res.X += q*(outX-inX);
						sv_res.Y[h] -= q*inY;

						if (mutY)
							sv_res.Y[hp] += q*outY;
						else
							sv_res.Y[h] += q*outY;

						sv_res.V[h] -= q*inV;
						if (mutV)
							sv_res.V[hp] += q*outV;
						else
							sv_res.V[h] += q*outV;
					}

				} else {

					int q = poissonian(a*tau, buf);

					sv_res.X += q*(outX-inX);
					sv_res.Y[h] += q*(outY-inY);
					sv_res.V[h] += q*(outV-inV);

				}

			}
			
		}

};

class MomentScalar {
	public:

		int Nsamples;
		double (*samplefunc)(StateVec &);
		std::string name;

		std::vector<double> mean;
		std::vector<double> var;
		std::vector<double> sem;

		// Constructor:
		MomentScalar(int p_Nsamples, double (*p_samplefunc)(StateVec &), std::string p_name) {
			Nsamples = p_Nsamples;
			samplefunc = p_samplefunc;
			name = p_name;

			mean.resize(Nsamples);
			var.resize(Nsamples);
			sem.resize(Nsamples);
		}
		MomentScalar() { };

		// Sample moment:
		void sample(StateVec & sv, int samp) {
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

class MomentVector {
	public:

		int Nsamples;
		std::vector<double> (*samplefunc)(StateVec &);
		std::string name;

		std::vector<std::vector<double> > mean;
		std::vector<std::vector<double> > var;
		std::vector<std::vector<double> > sem;

		int sequenceL;

		// Constructor:
		MomentVector(int p_Nsamples, std::vector<double> (*p_samplefunc)(StateVec &),
				std::string p_name, int p_sequenceL) {

			Nsamples = p_Nsamples;
			samplefunc = p_samplefunc;
			name = p_name;
			sequenceL = p_sequenceL;

			mean.resize(Nsamples);
			var.resize(Nsamples);
			sem.resize(Nsamples);

			for (int s=0; s<Nsamples; s++) {
				mean[s].resize(sequenceL);
				var[s].resize(sequenceL);
				sem[s].resize(sequenceL);
			}
		}
		MomentVector() { };

		// Sample moment:
		void sample(StateVec & sv, int samp) {
			std::vector<double> tmp = (*samplefunc)(sv);
			for (int h=0; h<=sequenceL; h++) {
				mean[samp][h] += tmp[h];
				var[samp][h] += tmp[h]*tmp[h];
			}
		}

		// Post processing of moment data:
		void normalise (int Npaths) {
			for (int s=0; s<Nsamples; s++) {
				for (int h=0; h<=sequenceL; h++) {
					mean[s][h] /= Npaths;
					var[s][h] = var[s][h]/Npaths - mean[s][h]*mean[s][h];
					sem[s][h] = sqrt(var[s][h]/Npaths);
				}
			}
		}
};

double samplefunc_mean_X (StateVec & sv) {
	return sv.X;
}
std::vector<double> samplefunc_mean_Y (StateVec & sv) {
	return sv.Y;
}
std::vector<double> samplefunc_mean_V (StateVec & sv) {
	return sv.V;
}

int main (int argc, char **argv)
{
    using namespace std;

    // Parse command line parameters:
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " outfile" << endl;
        exit(0);
    }
    char *ofname = argv[1];

    // Simulation parameters:
	double T = 5.0;
	int Nt = 1001;
	int Nsamples = 1001;
	int Npaths = 1;

	// Derived simulation parameters:
	int steps_per_sample = (Nt-1)/(Nsamples-1);
	double dt = T/(Nt-1);

	// Genetic parameters:
	int sequenceL = 35*3; // DNA sequence length coresponding to V3
	double mu = 2e-5/3.0; // Mutation rate per character given outcome

    // Demographic parameters:
	double param_d = 1e-3;
	double param_a = 1.0;
	double param_u = 3.0;
	double param_lambda = 2.5e8;
	double param_beta = 5e-13;
	double param_k = 1e3;

	// Set up reactions:
	int Nreactions = 6;
	Reaction reactions[6];

	reactions[0] = Reaction(0,0,0, 1,0,0,
			false,false,
			param_lambda,0.0);

	reactions[1] = Reaction(1,0,1, 0,1,0,
			true,false,
			param_beta,mu);

	reactions[2] = Reaction(0,1,0, 0,1,1,
			false,false,
			param_k,0.0);

	reactions[3] = Reaction(1,0,0, 0,0,0,
			false,false,
			param_d,0.0);

	reactions[4] = Reaction(0,1,0, 0,1,0,
			false,false,
			param_a,0.0);

	reactions[5] = Reaction(0,0,1, 0,0,0,
			false,false,
			param_u,0.0);

	// Set up moments:
	int NScalarMoments = 1;
	MomentScalar scalarMoments[1];
	scalarMoments[0] = MomentScalar (Nsamples, samplefunc_mean_X, "X");

	int NVectorMoments = 2;
	MomentVector vectorMoments[2];
	vectorMoments[0] = MomentVector (Nsamples, samplefunc_mean_Y, "Y", sequenceL);
	vectorMoments[1] = MomentVector (Nsamples, samplefunc_mean_V, "V", sequenceL);

	// Set up initial condition:
	StateVec sv0(sequenceL);
	sv0.X = param_lambda/param_d;
	sv0.V[0] = 100;

    // Initialise RNG:
	unsigned short buf[3] = {42, 53, 1234};
	
	// Loop over paths:
	for (int path=0; path<Npaths; path++) {

		// Initialise state vector:
		StateVec sv = sv0;

		// Perform initial sample:
		for (int i=0; i<NScalarMoments; i++)
			scalarMoments[i].sample(sv, 0);
		for (int i=0; i<NVectorMoments; i++)
			vectorMoments[i].sample(sv, 0);
		

		// Simulation loop:
		for (int t_idx=0; t_idx < Nt; t_idx++) {

			// Sample if necessary:
			if (t_idx % steps_per_sample == 0) {
				for (int i=0; i<NScalarMoments; i++)
					scalarMoments[i].sample(sv, t_idx/steps_per_sample);
				for (int i=0; i<NVectorMoments; i++)
					vectorMoments[i].sample(sv, t_idx/steps_per_sample);
			}

			// Implement reactions:
			StateVec sv_new = sv;
			for (int r=0; r<Nreactions; r++)
				reactions[r].leap(dt, sv, sv_new, buf);

			sv = sv_new;
		}
	}

	// Perform moment post-processing:
	for (int i=0; i<NScalarMoments; i++)
		scalarMoments[i].normalise(Npaths);
	for (int i=0; i<NVectorMoments; i++)
		vectorMoments[i].normalise(Npaths);


	// TODO: Write moments to disk.

		
	// Done!
	exit(0);
}
