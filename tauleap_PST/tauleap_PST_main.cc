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
#include <map>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <boost/program_options.hpp>

// Ensure macro points to correct function:
#define H5Gcreate_vers 2

#include <hdf5.h>
#include <hdf5_hl.h>

#include <mpi.h>

#include "poissonian.h"
#include "tauleap_PST_classes.h"
#include "tauleap_PST_OptionParser.h"

// Functions which calculate the moments to be sampled:
double samplefunc_X (const StateVec & sv) {
	return sv.X;
}
double samplefunc_Vtot (const StateVec & sv) {
	double res = 0;
	for (int h=0; h<=sv.L; h++)
		res += sv.V[h];

	return res;
}
std::vector<double> samplefunc_Y (const StateVec & sv) {
	return sv.Y;
}
std::vector<double> samplefunc_V (const StateVec & sv) {
	return sv.V;
}
std::vector<double> samplefunc_V0Vh (const StateVec & sv) {

	std::vector<double> V0Vh(sv.L + 1, 0.0);

	for (int h=0; h<=sv.L; h++)
		V0Vh[h] = sv.V[0]*sv.V[h];

	return V0Vh;
}
std::vector<double> samplefunc_VVh (const StateVec & sv) {

	std::vector<double> res(sv.L+1, 0.0);

	for (int d=0; d<=sv.L; d++) {
		for (int h=0; h<=sv.L; h++) {
			int hp=(h+d)%(sv.L+1);
			res[h] += sv.V[h]*sv.V[hp];
		}
		res[d] /= sv.L+1;
	}

	return res;

}

int main (int argc, char **argv)
{
    using namespace std;

    // Parse command line parameters:
    boost::program_options::variables_map vm = OptionParser::parse(argc, argv);

    // Initialize MPI:
    MPI::Init(argc, argv);
    int mpi_size = MPI::COMM_WORLD.Get_size();
    int mpi_rank = MPI::COMM_WORLD.Get_rank();
    atexit(MPI::Finalize);

    // Ensure output filename ends with ".h5" and record basename:
    string ofbasename, ofname = vm["outfile"].as<string>();
    if (ofname.length()>3 && (ofname.compare(ofname.length()-3,3,".h5")==0))
    	ofbasename = ofname.substr(0,ofname.length()-3);
    else
    	ofbasename = ofname;
    ofname = ofbasename + ".h5";

    // Attempt to create output file:
    // (Better to fail before the calculation begins.)
    if (mpi_rank==0) {
    	hid_t file_id = H5Fcreate(ofname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id<0) {
    	   	cout << "Error: Cannot open file '" << ofname << "' for writing."
    	    		<< endl;
    	   	exit(1);
        }
        H5Fclose(file_id);
    }

    // Simulation parameters:
	double T = vm["simulation.T"].as<double>();
	int Nt = vm["algorithm.Nt"].as<int>();
	int Nsamples = vm["algorithm.Nsamples"].as<int>();
	int Npaths = vm["algorithm.Npaths"].as<int>();

	int Ncrit = vm["algorithm.Ncrit"].as<int>(); // Reaction criticality parameter
	bool newcritcond = false; // Use new criticality condition

	// Derived simulation parameters:
	int steps_per_sample = (Nt-1)/(Nsamples-1);
	double dt = T/(Nt-1);
	double sample_dt = T/(Nsamples-1);

	// Viral genome length
	int sequenceL = (int)(vm["model.sequenceL"].as<double>());

	// Set up initial condition:
	StateVec sv0(sequenceL);
	sv0.X = vm["simulation.X0"].as<double>();
	sv0.Y[0] = vm["simulation.Y0"].as<double>();
	sv0.V[0] = vm["simulation.V0"].as<double>();

	// Set up reactions:
	int Nreactions = 6;
	Reaction reactions[6];

	reactions[0] = Reaction(0,0,0, 1,0,0,
			false,false,
			vm["model.lambda"].as<double>(),0.0);

	reactions[1] = Reaction(1,0,1, 0,1,0,
			true,false,
			vm["model.beta"].as<double>(),vm["model.mu"].as<double>());

	reactions[2] = Reaction(0,1,0, 0,1,1,
			false,false,
			vm["model.k"].as<double>(),0.0);

	reactions[3] = Reaction(1,0,0, 0,0,0,
			false,false,
			vm["model.d"].as<double>(),0.0);

	reactions[4] = Reaction(0,1,0, 0,0,0,
			false,false,
			vm["model.a"].as<double>(),0.0);

	reactions[5] = Reaction(0,0,1, 0,0,0,
			false,false,
			vm["model.u"].as<double>(),0.0);

	// Set up moments:
	int NScalarMoments = 2;
	MomentScalar scalarMoments[2];
	scalarMoments[0] = MomentScalar (Nsamples, samplefunc_X, "X");
	scalarMoments[1] = MomentScalar (Nsamples, samplefunc_Vtot, "Vtot");

	int NVectorMoments = 4;
	MomentVector vectorMoments[4];
	vectorMoments[0] = MomentVector (Nsamples, samplefunc_Y, "Y", sequenceL);
	vectorMoments[1] = MomentVector (Nsamples, samplefunc_V, "V", sequenceL);
	vectorMoments[2] = MomentVector (Nsamples, samplefunc_V0Vh, "V0Vh", sequenceL);
	vectorMoments[3] = MomentVector (Nsamples, samplefunc_VVh, "VVh", sequenceL);

    // Initialise RNG:
	unsigned short buf[3] = {53, time(NULL), mpi_rank};
	
	// Determine number of paths to integrate on this node:
	int chunk_size = Npaths/mpi_size;
	if (Npaths%mpi_size>0)
		chunk_size++;
	Npaths = chunk_size*mpi_size;

	// Loop over paths:
	for (int path=0; path<chunk_size; path++) {

		cout << "Rank " << mpi_rank << ": Integrating path " << path+1 << " of " << chunk_size << "..." << endl;;

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

			// Loop to get us through a single interval of length dt
			double t=0.0;
			do {

				// Determine leap time to use
				double tau = dt-t;
				int crit_react = -1;
				for (int r=0; r<Nreactions; r++) {
					double thistau = reactions[r].getLeapDistance(tau, Ncrit, newcritcond, sv, buf);
					if (thistau<tau) {
						tau = thistau;
						crit_react = r;
					}
				}

				if (crit_react>=0) {
					// Implement critical reaction:
					reactions[crit_react].doCritical(sv_new);
				}

				// Implement reactions:
				for (int r=0; r<Nreactions; r++) {
					if (reactions[r].tauleap(tau, sv_new, buf)) {
						cout << "Rank " << mpi_rank << ": Warning: negative population generated at time t="
								<< dt*t_idx << " by reaction " << r << "." << endl;
					}
				}

				t += tau;

				sv = sv_new;

			} while (t<dt);
		}

	}

	// Non-root nodes send sampled data to root and exit.
	if (mpi_rank != 0) {

		cout << "Rank " << mpi_rank << ": Sending data to root node." << endl;

		int tag = 0;
		for (int i=0; i<NScalarMoments; i++)
			scalarMoments[i].mpi_send(mpi_rank, tag);
		for (int i=0; i<NVectorMoments; i++)
			vectorMoments[i].mpi_send(mpi_rank, tag);

		exit(0);
	}

	// Collect sampled data from non-root nodes
	// and perform moment post-processing:
	for (int recv_rank = 1; recv_rank<mpi_size; recv_rank++) {

		cout << "Rank 0: Receiving data from rank " << recv_rank << "." << endl;

		int tag = 0;
		for (int i=0; i<NScalarMoments; i++)
			scalarMoments[i].mpi_recv(recv_rank, tag);
		for (int i=0; i<NVectorMoments; i++)
			vectorMoments[i].mpi_recv(recv_rank, tag);
	}

	// Normalise results:
	for (int i=0; i<NScalarMoments; i++)
		scalarMoments[i].normalise(Npaths);
	for (int i=0; i<NVectorMoments; i++)
		vectorMoments[i].normalise(Npaths);

	cout << "Rank 0: Writing results to disk...";

	// Open database file:
	hid_t file_id = H5Fcreate(ofname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Create group with same name as output file
	hid_t group_id = H5Gcreate(file_id, ("/"+ofbasename).c_str(),
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write scalar moments to file:
	int scalar_rank = 1;
	hsize_t scalar_dims[1] = {Nsamples};

	for (int m=0; m<NScalarMoments; m++) {

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_mean").c_str(),
				scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].mean[0]));

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_var").c_str(),
				scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].var[0]));

		H5LTmake_dataset(group_id, (scalarMoments[m].name + "_sem").c_str(),
				scalar_rank, scalar_dims,
				H5T_NATIVE_DOUBLE, &(scalarMoments[m].sem[0]));

	}

	// Write vector moments to file:
	int vector_rank = 2;
	hsize_t vector_dims[2] = {Nsamples, sequenceL+1};

	for (int m=0; m<NVectorMoments; m++) {

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_mean").c_str(),
				vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].mean[0]));

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_var").c_str(),
				vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].var[0]));

		H5LTmake_dataset(group_id, (vectorMoments[m].name + "_sem").c_str(),
				vector_rank, vector_dims,
				H5T_NATIVE_DOUBLE, &(vectorMoments[m].sem[0]));

	}

	// Write sample times to file:
	vector<double> t_vec;
	t_vec.resize(Nsamples);
	for (int s=0; s<Nsamples; s++)
		t_vec[s] = sample_dt*s;
	H5LTmake_dataset(group_id, "t", scalar_rank, scalar_dims,
			H5T_NATIVE_DOUBLE, &(t_vec[0]));

	// Write simulation parameters to file:
	boost::program_options::variables_map::iterator it;
	for (it = vm.begin(); it != vm.end(); it++) {

		if (it->first.find("model.") == 0)
			H5LTset_attribute_double(file_id, ("/"+ofbasename).c_str(),
					it->first.c_str(), &(it->second.as<double>()), 1);

		if (it->first.find("simulation.") == 0)
			H5LTset_attribute_double(file_id, ("/"+ofbasename).c_str(),
					it->first.c_str(), &(it->second.as<double>()), 1);

		if (it->first.find("algorithm.") == 0)
			H5LTset_attribute_int(file_id, ("/"+ofbasename).c_str(),
					it->first.c_str(), &(it->second.as<int>()), 1);
	}

	// Close HDF file:
	H5Fclose(file_id);

	cout << "done." << endl;

	// Done!
	exit(0);
}
