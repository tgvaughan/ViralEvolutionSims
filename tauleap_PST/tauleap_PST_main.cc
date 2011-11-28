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
#include <boost/math/distributions.hpp>

// Ensure HDF group creation macro points to correct function:
#define H5Gcreate_vers 2

#include <hdf5.h>
#include <hdf5_hl.h>

#include <mpi.h>

#include "poissonian.h"
#include "tauleap_PST_OptionParser.h"
#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_Reaction.h"
#include "tauleap_PST_MomentScalar.h"
#include "tauleap_PST_MomentVector.h"
#include "tauleap_PST_SemiImplicit.h"

// Functions which calculate the moments to be sampled:
double samplefunc_X (const StateVec & sv) {
	return sv.X;
}

double samplefunc_Ytot (const StateVec & sv) {
	double res = 0;
	for (int h=0; h<=sv.L; h++)
		res += sv.Y[h];

	return res;
}

double samplefunc_YLtot (const StateVec & sv) {
	double res = 0;
	for (int h=0; h<=sv.L; h++)
		res += sv.YL[h];

	return res;
}

double samplefunc_Vtot (const StateVec & sv) {
	double res = 0;
	for (int h=0; h<=sv.L; h++)
		res += sv.V[h];

	return res;
}

double samplefunc_Vdiv (const StateVec & sv) {

	double N2 = 0.0;
	double N = 0.0;

	for (int h=0; h<=sv.L; h++) {
		N2 += sv.V[h]*sv.V[h];
		N += sv.V[h];
	}

	if (N>0)
		return exp(2*log(N)-log(N2));
	else
		return 1;
}

double samplefunc_Vext (const StateVec & sv) {

	int h;
	for (h=sv.L; h>=0; h--) {
		if (sv.V[h]>0)
			break;
	}

	if (h>=0)
		return h;
	else
		return 0;
}

void samplefunc_Y (const StateVec & sv, std::vector<double> & res) {
	res = sv.Y;
}

void samplefunc_YL (const StateVec & sv, std::vector<double> & res) {
	res = sv.YL;
}

void samplefunc_V (const StateVec & sv, std::vector<double> & res) {
	res = sv.V;
}

void samplefunc_V0Vh (const StateVec & sv, std::vector<double> & res) {

	for (int h=0; h<=sv.L; h++)
		res[h] = sv.V[0]*sv.V[h];
}

void samplefunc_VhVhp (const StateVec & sv, std::vector<double> & res) {

	for (int h1=0; h1<=sv.L; h1++) {
		for (int h2=0; h2<=sv.L; h2++)
			res[h1*(sv.L+1)+h2] = sv.V[h1]*sv.V[h2];
	}
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

	// Select criticality condition to use in hybrid scheme:
	// 0: none, tau-leaping only
	// 1: standard
	// 2: modified (variance-dependent)
	int critcond = 1;

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
	sv0.YL[0] = vm["simulation.YL0"].as<double>();
	sv0.V[0] = vm["simulation.V0"].as<double>();

	// Set up reactions:
	std::vector<Reaction> reactions;

	reactions.push_back(Reaction(0,0,0,0, 1,0,0,0,
			false,false,false,
			vm["model.lambda"].as<double>(),0.0));

	reactions.push_back(Reaction(1,0,0,1, 0,1,0,0,
			true,false,false,
			vm["model.beta"].as<double>()*(1.0-vm["model.lat_p"].as<double>()),
			vm["model.mu"].as<double>()));

	reactions.push_back(Reaction(1,0,0,1, 0,0,1,0,
			false,true,false,
			vm["model.beta"].as<double>()*vm["model.lat_p"].as<double>(),
			vm["model.mu"].as<double>()));

	reactions.push_back(Reaction(0,0,1,0, 0,1,0,0,
			false,false,false,
			vm["model.lat_a"].as<double>(), 0.0));

	reactions.push_back(Reaction(0,1,0,0, 0,1,0,1,
			false,false,true,
			vm["model.k"].as<double>(),vm["model.mu_RNAP"].as<double>()));

	reactions.push_back(Reaction(1,0,0,0, 0,0,0,0,
			false,false,false,
			vm["model.d"].as<double>(),0.0));

	reactions.push_back(Reaction(0,1,0,0, 0,0,0,0,
			false,false,false,
			vm["model.a"].as<double>(),0.0));

	reactions.push_back(Reaction(0,0,1,0, 0,0,0,0,
			false,false,false,
			vm["model.lat_d"].as<double>(),0.0));

	reactions.push_back(Reaction(0,0,0,1, 0,0,0,0,
			false,false,false,
			vm["model.u"].as<double>(),0.0));

	// Set up moments:
	std::vector<MomentScalar> scalarMoments;
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_X, "X"));
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_Ytot, "Ytot"));
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_YLtot, "YLtot"));
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_Vtot, "Vtot"));
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_Vdiv, "Vdiv"));
	scalarMoments.push_back(MomentScalar (Nsamples, samplefunc_Vext, "Vext"));

	std::vector<MomentVector> vectorMoments;
	vectorMoments.push_back(MomentVector (Nsamples, samplefunc_Y, "Y", sequenceL+1));
	vectorMoments.push_back(MomentVector (Nsamples, samplefunc_YL, "YL", sequenceL+1));
	vectorMoments.push_back(MomentVector (Nsamples, samplefunc_V, "V", sequenceL+1));
	vectorMoments.push_back(MomentVector (Nsamples, samplefunc_V0Vh, "V0Vh", sequenceL+1));
	//vectorMoments.push_back(MomentVector (Nsamples, samplefunc_VhVhp, "VhVhp", (sequenceL+1)*(sequenceL+1)));

    // Initialise RNG:
	//unsigned short buf[3] = {53, time(NULL), mpi_rank};
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
		for (unsigned int i=0; i<scalarMoments.size(); i++)
			scalarMoments[i].sample(sv, 0);
		for (unsigned int i=0; i<vectorMoments.size(); i++)
			vectorMoments[i].sample(sv, 0);

		// Simulation loop:
		for (int t_idx=1; t_idx < Nt; t_idx++) {

			// Sample if necessary:
			if (t_idx % steps_per_sample == 0) {
				for (unsigned int i=0; i<scalarMoments.size(); i++)
					scalarMoments[i].sample(sv, t_idx/steps_per_sample);
				for (unsigned int i=0; i<vectorMoments.size(); i++)
					vectorMoments[i].sample(sv, t_idx/steps_per_sample);
			}

			// Loop to get us through a single interval of length dt
			double t=0.0;
			do {

				// Determine leap time to use
				double tau = dt-t;
				int crit_react = -1;
				for (unsigned int r=0; r<reactions.size(); r++) {
					double thistau = reactions[r].getLeapDistance(tau, Ncrit, critcond, sv, buf);
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
				for (unsigned int r=0; r<reactions.size(); r++) {
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
		for (unsigned int i=0; i<scalarMoments.size(); i++)
			scalarMoments[i].mpi_send(mpi_rank, tag);
		for (unsigned int i=0; i<vectorMoments.size(); i++)
			vectorMoments[i].mpi_send(mpi_rank, tag);

		exit(0);
	}

	// Collect sampled data from non-root nodes
	// and perform moment post-processing:
	for (int recv_rank = 1; recv_rank<mpi_size; recv_rank++) {

		cout << "Rank 0: Receiving data from rank " << recv_rank << "." << endl;

		int tag = 0;
		for (unsigned int i=0; i<scalarMoments.size(); i++)
			scalarMoments[i].mpi_recv(recv_rank, tag);
		for (unsigned int i=0; i<vectorMoments.size(); i++)
			vectorMoments[i].mpi_recv(recv_rank, tag);
	}

	// Normalise results:
	for (unsigned int i=0; i<scalarMoments.size(); i++)
		scalarMoments[i].normalise(Npaths);
	for (unsigned int i=0; i<vectorMoments.size(); i++)
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

	for (unsigned int m=0; m<scalarMoments.size(); m++) {

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

	for (unsigned int m=0; m<vectorMoments.size(); m++) {

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
