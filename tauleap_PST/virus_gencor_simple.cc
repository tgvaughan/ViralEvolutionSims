// Within-host viral evolution simulation.  Makes use of a projected
// sequence space of dramatically reduced dimension.

#include <iostream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "poissonian.h"


// Class for system state vectors:
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
            neighbourNum = sv.neighbourNum;

            X = sv.X;
            Y = sv.Y;
            V = sv.V;
        }

        // Assignment operator:
        StateVec operator= (const StateVec & src) {

        	L = src.L;
        	neighbourNum = src.neighbourNum;

        	X = src.X;
        	Y = src.Y;
        	V = src.V;

        	return *this;
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
		// (Returns true on success, false if negative populations were generated.)
		bool leap (double tau, StateVec & sv, StateVec & sv_new, unsigned short *buf) {

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

						double ap = a*mutrate*gcond(hp,h,sv.L);
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

};

// Class for scalar moments:
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

// Class for vector moments:
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

// Sampling functions for use in moment calculations:
double samplefunc_X (const StateVec & sv) {
	return sv.X;
}
std::vector<double> samplefunc_Y (const StateVec & sv) {
	return sv.Y;
}
std::vector<double> samplefunc_V (const StateVec & sv) {
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

    // Ensure output filename ends with ".h5" and record basename:
    string ofname = argv[1];
    string ofbasename;
    if (ofname.length()>3 && (ofname.substr(ofname.length()-4)==".h5"))
    	ofbasename = ofname.substr(0,ofname.length()-4);
    else
    	ofbasename = ofname;
    ofname = ofbasename + ".h5";

    // Attempt to create output file:
    // (Better to fail before the calculation begins.)
    hid_t file_id = H5Fcreate(ofname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id<0) {
    	cout << "Error: Cannot open file '" << ofname << "' for writing."
    			<< endl;
    	exit(1);
    }

    // Simulation parameters:
	double T = 30.0;
	int Nt = 10001;
	int Nsamples = 1001;
	int Npaths = 1;

	// Derived simulation parameters:
	int steps_per_sample = (Nt-1)/(Nsamples-1);
	double dt = T/(Nt-1);
	double sample_dt = T/(Nsamples-1);

	// Genetic parameters:
	int sequenceL = 35*3; // DNA sequence length coresponding to V3
	double mu = 2e-5/3.0; // Mutation rate per character given outcome

    // Demographic parameters:
	/*double param_d = 1e-3;
	double param_a = 1.0;
	double param_u = 3.0;
	double param_lambda = 2.5e8;
	double param_beta = 5e-13;
	double param_k = 1e3;*/

	double param_d = 0.1;
	double param_a = 0.5;
	double param_u = 5.0;
	double param_lambda = 1e5;
	double param_beta = 2e-7;
	double param_k = 100;

	// Set up reactions:
	int Nreactions = 6;
	Reaction reactions[6];

	reactions[0] = Reaction(0,0,0, 1,0,0,
			false,false,
			param_lambda,0.0);

	reactions[1] = Reaction(1,0,1, 0,1,0,
			false,false,
			param_beta,mu);

	reactions[2] = Reaction(0,1,0, 0,1,1,
			false,false,
			param_k,0.0);

	reactions[3] = Reaction(1,0,0, 0,0,0,
			false,false,
			param_d,0.0);

	reactions[4] = Reaction(0,1,0, 0,0,0,
			false,false,
			param_a,0.0);

	reactions[5] = Reaction(0,0,1, 0,0,0,
			false,false,
			param_u,0.0);

	// Set up moments:
	int NScalarMoments = 1;
	MomentScalar scalarMoments[1];
	scalarMoments[0] = MomentScalar (Nsamples, samplefunc_X, "X");

	int NVectorMoments = 2;
	MomentVector vectorMoments[2];
	vectorMoments[0] = MomentVector (Nsamples, samplefunc_Y, "Y", sequenceL);
	vectorMoments[1] = MomentVector (Nsamples, samplefunc_V, "V", sequenceL);

	// Set up initial condition:
	StateVec sv0(sequenceL);
	sv0.X = param_lambda/param_d;
	sv0.V[0] = 100;

    // Initialise RNG:
	unsigned short buf[3] = {42, 53, 1534};
	
	// Loop over paths:
	for (int path=0; path<Npaths; path++) {

		// Initialize state vector:
		StateVec sv = sv0;
		StateVec sv_new = sv0;

		// Perform initial sample:
		for (int i=0; i<NScalarMoments; i++)
			scalarMoments[i].sample(sv, 0);
		for (int i=0; i<NVectorMoments; i++)
			vectorMoments[i].sample(sv, 0);
		

		// Simulation loop:
		for (int t_idx=1; t_idx < Nt; t_idx++) {

			// Sample if necessary:
			if (t_idx % steps_per_sample == 0) {
				for (int i=0; i<NScalarMoments; i++)
					scalarMoments[i].sample(sv, t_idx/steps_per_sample);
				for (int i=0; i<NVectorMoments; i++)
					vectorMoments[i].sample(sv, t_idx/steps_per_sample);

				std::cout << "Recording sample " << t_idx/steps_per_sample + 1 << " of " << Nsamples << std::endl;
			}

			// Implement reactions:
			for (int r=0; r<Nreactions; r++) {
				if (!reactions[r].leap(dt, sv, sv_new, buf)) {
					cout << "Error: negative population generated. Exiting..." << endl;
					H5close();
					exit(1);
				}
			}

			sv = sv_new;
		}
	}

	// Perform moment post-processing:
	for (int i=0; i<NScalarMoments; i++)
		scalarMoments[i].normalise(Npaths);
	for (int i=0; i<NVectorMoments; i++)
		vectorMoments[i].normalise(Npaths);

	// Create group with same name as output file
	hid_t group_id = H5Gcreate(file_id, ("/"+ofbasename).c_str(),
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write scalar moments to file:
	int scalar_rank = 1;
	hsize_t scalar_dims[1] = {Nsamples};

	for (int m=0; m<NScalarMoments; m++) {

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_mean").c_str(), scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].mean[0]));

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_var").c_str(), scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].var[0]));

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_sem").c_str(), scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].sem[0]));

	}

	// Write vector moments to file:
	int vector_rank = 2;
	hsize_t vector_dims[2] = {Nsamples, sequenceL+1};

	for (int m=0; m<NVectorMoments; m++) {

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_mean").c_str(), vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].mean[0]));

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_var").c_str(), vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].var[0]));

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_sem").c_str(), vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].sem[0]));

	}

	// Write sample times to file:
	vector<double> t_vec;
	t_vec.resize(Nsamples);
	for (int s=0; s<Nsamples; s++)
		t_vec[s] = sample_dt*s;
	H5LTmake_dataset(group_id, "t", scalar_rank, scalar_dims,
			H5T_NATIVE_DOUBLE, &(t_vec[0]));

	// Close HDF file:
	H5Fclose(file_id);

	// Done!
	exit(0);
}
