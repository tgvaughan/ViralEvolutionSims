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
#include "tauleap_PST_Moment.h"
#include "tauleap_PST_SemiImplicit.h"

// Functions which calculate the moments to be sampled:
void samplefunc_X (const StateVec & sv, std::vector<double> & res) {
	res[0] = sv.X;
}

void samplefunc_Ytot (const StateVec & sv, std::vector<double> & res) {
	res[0] = 0;
	for (int h=0; h<=sv.L; h++)
		res[0] += sv.Y[h];
}

void samplefunc_YLtot (const StateVec & sv, std::vector<double> & res) {
	res[0] = 0;
	for (int h=0; h<=sv.L; h++)
		res[0] += sv.YL[h];
}

void samplefunc_Vtot (const StateVec & sv, std::vector<double> & res) {
	res[0] = 0;
	for (int h=0; h<=sv.L; h++)
		res[0] += sv.V[h];
}

void samplefunc_Vdiv (const StateVec & sv, std::vector<double> & res) {

	double N2 = 0.0;
	double N = 0.0;

	for (int h=0; h<=sv.L; h++) {
		N2 += sv.V[h]*sv.V[h];
		N += sv.V[h];
	}

	if (N>0)
		res[0] = exp(2*log(N)-log(N2));
	else
		res[0] = 1;
}

void samplefunc_Vext (const StateVec & sv, std::vector<double> & res) {

	int h;
	for (h=sv.L; h>=0; h--) {
		if (sv.V[h]>0)
			break;
	}

	if (h>=0)
		res[0] = h;
	else
		res[0] = 0;
}

void samplefunc_Y (const StateVec & sv, std::vector<double> & res) {
	for (int h=0; h<=20; h++)
		res[h] = sv.Y[h];
}

void samplefunc_YL (const StateVec & sv, std::vector<double> & res) {
	for (int h=0; h<=20; h++)
		res[h] = sv.YL[h];
}

void samplefunc_V (const StateVec & sv, std::vector<double> & res) {
	for (int h=0; h<=20; h++)
		res[h] = sv.V[h];
}

void samplefunc_YhYhp (const StateVec & sv, std::vector<double> & res) {

	for (int h1=0; h1<20; h1++) {
		for (int h2=0; h2<=20; h2++)
			res[h1*21+h2] = sv.Y[h1]*sv.Y[h2];
	}
}

void samplefunc_YLhYLhp (const StateVec & sv, std::vector<double> & res) {

	for (int h1=0; h1<20; h1++) {
		for (int h2=0; h2<=20; h2++)
			res[h1*21+h2] = sv.Y[h1]*sv.Y[h2];
	}
}

void samplefunc_VhVhp (const StateVec & sv, std::vector<double> & res) {

	for (int h1=0; h1<20; h1++) {
		for (int h2=0; h2<=20; h2++)
			res[h1*21+h2] = sv.V[h1]*sv.V[h2];
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
			vm["model.mu_RT"].as<double>()));

	reactions.push_back(Reaction(1,0,0,1, 0,0,1,0,
			false,true,false,
			vm["model.beta"].as<double>()*vm["model.lat_p"].as<double>(),
			vm["model.mu_RT"].as<double>()));

	reactions.push_back(Reaction(0,0,1,0, 0,1,0,0,
			false,false,false,
			vm["model.lat_a"].as<double>(), 0.0));

	reactions.push_back(Reaction(0,1,0,0, 0,1,0,1,
			false,false,true,
			vm["model.k"].as<double>(),
			vm["model.mu_RNAP"].as<double>()));

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

	std::vector<Moment> moments;
	std::vector<int> dim;

	moments.push_back(Moment (Nsamples, samplefunc_X, "X", dim));
	moments.push_back(Moment (Nsamples, samplefunc_Ytot, "Ytot", dim));
	moments.push_back(Moment (Nsamples, samplefunc_YLtot, "YLtot", dim));
	moments.push_back(Moment (Nsamples, samplefunc_Vtot, "Vtot", dim));
	moments.push_back(Moment (Nsamples, samplefunc_Vdiv, "Vdiv", dim));
	moments.push_back(Moment (Nsamples, samplefunc_Vext, "Vext", dim));

	//dim.push_back(sequenceL+1);
	dim.push_back(21);
	moments.push_back(Moment (Nsamples, samplefunc_Y, "Y", dim));
	moments.push_back(Moment (Nsamples, samplefunc_YL, "YL", dim));
	moments.push_back(Moment (Nsamples, samplefunc_V, "V", dim));

	//dim.push_back(sequenceL+1);
	dim.push_back(21);
	moments.push_back(Moment (Nsamples, samplefunc_YhYhp, "YhYhp", dim));
	moments.push_back(Moment (Nsamples, samplefunc_YhYhp, "YLhYLhp", dim));
	moments.push_back(Moment (Nsamples, samplefunc_VhVhp, "VhVhp", dim));

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
		for (unsigned int i=0; i<moments.size(); i++)
			moments[i].sample(sv, 0);

		// Simulation loop:
		for (int t_idx=1; t_idx < Nt; t_idx++) {

			// Sample if necessary:
			if (t_idx % steps_per_sample == 0) {
				for (unsigned int i=0; i<moments.size(); i++)
					moments[i].sample(sv, t_idx/steps_per_sample);
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
		for (unsigned int i=0; i<moments.size(); i++)
			moments[i].mpi_send(mpi_rank, tag);

		exit(0);
	}

	// Collect sampled data from non-root nodes
	// and perform moment post-processing:
	for (int recv_rank = 1; recv_rank<mpi_size; recv_rank++) {

		cout << "Rank 0: Receiving data from rank " << recv_rank << "." << endl;

		int tag = 0;
		for (unsigned int i=0; i<moments.size(); i++)
			moments[i].mpi_recv(recv_rank, tag);
	}

	// Normalise results:
	for (unsigned int i=0; i<moments.size(); i++)
		moments[i].normalise(Npaths);

	cout << "Rank 0: Writing results to disk...";

	// Open database file:
	hid_t file_id = H5Fcreate(ofname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Create group with same name as output file
	hid_t group_id = H5Gcreate(file_id, ("/"+ofbasename).c_str(),
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write moments moments to file:
	for (unsigned int m=0; m<moments.size(); m++) {

		std::vector<hsize_t> dims(moments[m].dim.size()+1);
		dims[0] = Nsamples;
		for (unsigned int i=0; i<moments[m].dim.size(); i++)
			dims[i+1] = (hsize_t)moments[m].dim[i];

		H5LTmake_dataset(group_id, (moments[m].name + "_mean").c_str(),
				dims.size(), &(dims[0]),
				H5T_NATIVE_DOUBLE, &(moments[m].mean[0]));

		H5LTmake_dataset(group_id, (moments[m].name + "_var").c_str(),
				dims.size(), &(dims[0]),
				H5T_NATIVE_DOUBLE, &(moments[m].var[0]));

		H5LTmake_dataset(group_id, (moments[m].name + "_sem").c_str(),
				dims.size(), &(dims[0]),
				H5T_NATIVE_DOUBLE, &(moments[m].sem[0]));
	}

	// Write sample times to file:
	vector<double> t_vec;
	t_vec.resize(Nsamples);
	for (int s=0; s<Nsamples; s++)
		t_vec[s] = sample_dt*s;
	vector<hsize_t> t_dims(1,Nsamples);
	H5LTmake_dataset(group_id, "t", t_dims.size(), &(t_dims[0]),
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
