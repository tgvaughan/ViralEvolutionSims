// virus_gencor_tauleap
//
// Uses the hybrid tau-leaping algorithm of Cao et al. to
// simulate the dynamics of a genetically diverse viral
// infection.
//

#include <iostream>
#include <fstream>

#include <vector>
#include <set>

#include <cstdlib>
#include <cmath>
#include <ctime>

#include <mpi.h>

#include "poissonian.h"
#include "virus_gencor_tauleap.h"

// Moment calculation functions:
double samplefunc_x(StateVec x) {
	return x[0];
}
double samplefunc_y(StateVec x) {
	return x[1];
}
double samplefunc_v(StateVec x) {
	return x[2];
}
double samplefunc_xy(StateVec x) {
	return x[0]*x[1];
}
double samplefunc_xv(StateVec x) {
	return x[0]*x[2];
}
double samplefunc_yv(StateVec x) {
	return x[1]*x[2];
}
double samplefunc_clear(StateVec x) {
	return (double)(x[1]<0.5 && x[2]<0.5);
}

int main (int argc, char **argv)
{
	using namespace std;

	// Parse command line parameters:
	
	if (argc<4) {
		cout << "Usage: " << argv[0] << " Nv0 seed outfile" << endl;
		exit(0);
	}
	double Nv0 = (int)strtod(argv[1], NULL);
	unsigned short seed = (unsigned short)strtol(argv[2], NULL, 10);
	char *ofname = argv[3];

	// MPI Initialisation:
	MPI::Init(argc, argv);
	int mpi_rank = MPI::COMM_WORLD.Get_rank();
	int mpi_size = MPI::COMM_WORLD.Get_size();

	// Simulation parameters:
	double T = 200.0;		// Total simulation time
	int Ntraj = 2048;		// Number of trajectories to generate
	int Nt_full = 20001;	// Number of full-sized tau-leaps
	int Nc = 100;			// Critical reaction number
	int Nsamples = 10001;	// Number of samples to record

	// Calculation conditional on no extinction?
	bool conditional = false;

	// Model parameters:
	double d = 1e-3;
	double a = 1.0;
	double u = 3.0;
	double lambda = 2.5e8;
	double beta = 5e-13;
	double k = 1e3;

	// Derived simulation parameters:
	int Nt[2] = {Nt_full, (Nt_full-1)*2 + 1};
	double dt[2] = {T/(Nt[0]-1), T/(Nt[1]-1)};
	int steps_per_sample[2] = {(Nt[0]-1)/(Nsamples-1), (Nt[1]-1)/(Nsamples-1)};

	// Set up reactions:
	int Nreactions = 6;
	Reaction reactions[6];
	StateVec in(3), out(3);

	// T cell production
	in[0] = 0; in[1] = 0; in[2] = 0;
	out[0] = 1; out[1] = 0; out[2] = 0;
	reactions[0] = Reaction(in, out, lambda);

	// T cell infection
	in[0] = 1; in[1] = 0; in[2] = 1;
	out[0] = 0; out[1] = 1; out[2] = 0;
	reactions[1] = Reaction(in, out, beta);

	// Virus production
	in[0] = 0; in[1] = 1; in[2] = 0;
	out[0] = 0; out[1] = 1; out[2] = 1;
	reactions[2] = Reaction(in, out, k);

	// T cell death
	in[0] = 1; in[1] = 0; in[2] = 0;
	out[0] = 0; out[1] = 0; out[2] = 0;
	reactions[3] = Reaction(in, out, d);

	// Infected T cell death
	in[0] = 0; in[1] = 1; in[2] = 0;
	out[0] = 0; out[1] = 0; out[2] = 0;
	reactions[4] = Reaction(in, out, a);

	// Virion clearance
	in[0] = 0; in[1] = 0; in[2] = 1;
	out[0] = 0; out[1] = 0; out[2] = 0;
	reactions[5] = Reaction(in, out, u);

	// Initial state:
	StateVec x0(3);
	x0[0] = 2.5e11; // Uninfected
	x0[1] = 0; 		// Infected
	x0[2] = Nv0;	// Virus

	// Allocate memory for recording moments:
	int Nmoments = 7;
	Moment moments[7] = {
		Moment(Nsamples, &samplefunc_x, "x"),
		Moment(Nsamples, &samplefunc_y, "y"),
		Moment(Nsamples, &samplefunc_v, "v"),
		Moment(Nsamples, &samplefunc_xy, "xy"),
		Moment(Nsamples, &samplefunc_xv, "xv"),
		Moment(Nsamples, &samplefunc_yv, "yv"),
		Moment(Nsamples, &samplefunc_clear, "clear")
	};

	// Initialise PRNG:
	unsigned short buf[3];
	buf[0] = 42;
	buf[1] = seed;
	buf[2] = (unsigned short)mpi_rank;

	// Divvy out trajectories:
	int chunk_size = Ntraj/mpi_size;
	if (Ntraj%mpi_size > 0)
		chunk_size++;
	Ntraj = mpi_size*chunk_size;
	// (Makes sure no node is sitting idle and
	// _at least_ Ntraj trajectories are calculated.)

	// Start timer:
	time_t sim_time = time(NULL);


	/////////////////////////////////
	// Trajectory loop:
	//

	for (int traj=0; traj<chunk_size; traj++) {

		for (int phase=0; phase<=1; phase++) {

			// Report simulation phase:
			switch(phase) {
				case 0:
					cout << "Rank " << mpi_rank << ": Integrating FULL step trajectory "
						<< traj+1 << " of " << chunk_size << "..." << endl;
					break;
				case 1:
					cout << "Rank " << mpi_rank << ": Integrating HALF step trajectory "
						<< traj+1 << " of " << chunk_size << "..." << endl;
					break;
			}

			// Initialise system state:
			StateVec x(3);
			x[0] = poissonian(x0[0], buf);
			x[1] = poissonian(x0[1], buf);
			//x[2] = poissonian(x0[2], buf);
			x[2] = x0[2];

			// Initialise sampling index:
			int sidx = 0;


			/////////////////////////////////
			// Integration loop:
			//

			for (int tidx=0; tidx<Nt[phase]; tidx++) {

				// Check for extinction in conditional calculation:
				if (conditional && x[1]<0.5 && x[2]<0.5)
					break;

				// Sample if necessary:
				if (tidx % steps_per_sample[phase] == 0) {

					// Sample moments:
					for (int m=0; m<Nmoments; m++) {
						if (phase==0)
							moments[m].record(x, sidx);
						else
							moments[m].record_half(x, sidx);
					}

					// Increment sampling index:
					sidx++;

				}

				// Hybrid integration step:
				double t = 0.0;
				while (true) {

					int critReaction = -1;
					double delta = dt - t;

					// Calculate propensities and determine critical reactions:
					for (int r=0; r<Nreactions; r++) {

						reactions[r].calcPropensity(x, buf);

						if (reactions[r].isCritical() && reactions.criticalDelta() < delta) {
							critReaction = r;
							delta = reactions.criticalDelta();
						}

					}

					// Perform tau-leap:
					for (int r=0; r<Nreactions; r++)
						x = reactions[r].tauLeap(x, delta, buf);

					// Will a critical reaction occur before end of current interval?
					if (critReaction<0)
						break;

					// Implement chosen critical reaction:
					x = reactions[critReaction].implementCritical(x);

					// Increment within-interval time:
					t += delta;

				}

				// Check for negative values:
				if (x.isnegative()) {
					cout << "FATAL ERROR: Negative population size generated!  Aborting..." << endl;
					MPI::COMM_WORLD.Abort(1);
				}

			}

			// Check for extinction event:
			if (conditional && x[1]<0.5 && x[2]<0.5) {
				cout << "Rank " << mpi_rank << " extinction." << endl;
				phase--;
				continue;
			}

			// Record last sample:
			for (int m=0; m<Nmoments; m++) {
				if (phase==0)
					moments[m].record(x, sidx);
				else
					moments[m].record_half(x, sidx);
			}

		}

		// Accept recorded samples:
		for (int m=0; m<Nmoments; m++)
			moments[m].add();

	}

	cout << "Rank " << mpi_rank << ": finished trajectory integrations." << endl;

	// Stop timer:
	sim_time = time(NULL) - sim_time;


	/////////////////////////////////
	// Collate integration results:
	//

	if (mpi_rank>0) {

		// Send moments to root node:
		
		cout << "Rank " << mpi_rank << ": Sending results to rank 0...";

		int tag=0;

		for (int m=0; m<Nmoments; m++)
			moments[m].mpi_send(&tag);

		cout << " sent." << endl;
	
	} else {

		// Collate moments from other nodes:

		for (int from_rank=1; from_rank<mpi_size; from_rank++) {

			cout << "Rank 0: Receiving results from rank " << from_rank << "...";

			int tag=0;

			for (int m=0; m<Nmoments; m++)
				moments[m].mpi_receive(&tag, from_rank);

			cout << " received." << endl;

		}


		// Normalise moments:

		for (int m=0; m<Nmoments; m++)
			moments[m].normalise(Ntraj);

		// Write results to file:
	
		ofstream ofile;
		ofile.open(ofname);

		ofile << "# Tau-leaping simulation of stochastic viral infection dynamics" << endl;
		ofile << "# " << endl;
		ofile << "# Ntraj = " << Ntraj << endl;
		ofile << "# T = " << T << endl;
		ofile << "# Nt = " << Nt << endl;
		ofile << "# Nsamples = " << Nsamples << endl;
		ofile << "# " << endl;
		ofile << "# Model Parameters: " << endl;
		ofile << "# d = " << d << endl;
		ofile << "# a = " << a << endl;
		ofile << "# u = " << u << endl;
		ofile << "# lambda = " << lambda << endl;
		ofile << "# k = " << k << endl;
		ofile << "# beta = " << beta << endl;
		ofile << "# " << endl;
		ofile << "# Simulation took " << sim_time << " seconds to complete." << endl << endl;

		ofile << "t ";
		for (int m=0; m<Nmoments; m++) {
			ofile << moments[m].name << "_mean ";
			ofile << moments[m].name << "_half_mean ";
			ofile << moments[m].name << "_var ";
			ofile << moments[m].name << "_sem ";
			ofile << moments[m].name << "_ftse ";
		}
		ofile << endl;

		int sidx = 0;
		for (int s=0; s<Nsamples; s++) {

			ofile << dt[0]*(double)(s*steps_per_sample[0]) << " ";

			for (int m=0; m<Nmoments; m++) {
				ofile << moments[m].mean[s] << " ";
				ofile << moments[m].half_mean[s] << " ";
				ofile << moments[m].var[s] << " ";
				ofile << moments[m].sem[s] << " ";
				ofile << moments[m].ftse[s] << " ";
			}
			ofile << endl;
		}

		ofile.close();

	}

	// MPI clean-up:
	MPI::Finalize();

	// Done!
	exit(0);
}
