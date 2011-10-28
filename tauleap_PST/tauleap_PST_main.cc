/*
 * tauleap_PST_main.cc
 *
 * Within-host viral evolution simulation.  Makes use of a projected
 * sequence space of dramatically reduced dimension.
 *
 *  Created on: 28/10/2011
 *      Author: Tim Vaughan
 *
 */

#include <iostream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <ctime>

// Ensure macro points to correct function:
#define H5Gcreate_vers 2

#include <hdf5.h>
#include <hdf5_hl.h>

#include "poissonian.h"
#include "tauleap_PST_classes.h"


// Functions which calculate the moments to be sampled:
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
        cout << "Usage: " << argv[0] << " outfile[.h5]" << endl;
        exit(0);
    }

    // Ensure output filename ends with ".h5" and record basename:
    string ofname = argv[1];
    string ofbasename;
    if (ofname.length()>3 && (ofname.compare(ofname.length()-4,3,".h5")))
    	ofbasename = ofname.substr(0,ofname.length()-3);
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
	int Nt = 3001;
	int Nsamples = 1001;
	int Npaths = 1;

	// Derived simulation parameters:
	int steps_per_sample = (Nt-1)/(Nsamples-1);
	double dt = T/(Nt-1);
	double sample_dt = T/(Nsamples-1);

	// Genetic parameters:
	int sequenceL = 35*3; // DNA sequence length corresponding to V3
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
	unsigned short buf[3] = {42, 53, time(NULL)};
	
	// Loop over paths:
	for (int path=0; path<Npaths; path++) {

		cout << "Path " << path+1 << " of " << Npaths << "...";
		cout.flush();

		// Initialise state vector:
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

		cout << "done." << endl;
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
