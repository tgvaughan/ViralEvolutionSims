/* virus_gencor: perform simulations based on the standard virus
   model with mutation in order to probe the qualities of the
   resulting genetic correlations. */

#include <iostream>
#include <fstream>

#include <map>
#include <vector>

#include <cstdlib>
#include <cmath>
#include <ctime>

#include <mpi.h>

#include "virus_gencor.h"

void sample(unsigned int genome_L,
		std::vector<std::vector<unsigned long> > & jointpcount,
		int samp,
		StateVec & sv, std::vector<std::vector<double> > & momentvec)
{
	if (!sv.totals_valid)
		sv.calc_totals();

	std::vector<int> seq0(genome_L, 0);
	int m = 0; // moment vector index

	// SCALAR: mean_NX
	momentvec[m++][samp] += (double)sv.X;

	// SCALAR: mean_NY
	momentvec[m++][samp] += (double)sv.NY;

	// SCALAR: mean_NV
	momentvec[m++][samp] += (double)sv.NV;

	// SCALAR: mean_NVsq
	momentvec[m++][samp] += (double)sv.NV*(double)sv.NV;

	// VECTOR: averages of HD-unique virus subpopulations
	for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
			it != sv.V.end(); ++it)
	{
		int HD = 0;
		for (unsigned int i=0; i<genome_L; i++)
			if ((it->first)[i] != seq0[i])
				HD++;

		momentvec[m+HD][samp] += (double)(it->second);
	}
	m += genome_L+1;

	// VECTOR: averages of products between HD-unique virus subpopulations
	std::vector<double> thisV(genome_L+1, 0);
	for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
			it != sv.V.end(); ++it)
	{
		int HD = 0;
		for (unsigned int i=0; i<genome_L; i++)
			if ((it->first)[i] != seq0[i])
				HD++;

		thisV[HD] += it->second;
	}
	for (unsigned int HD=0; HD<=genome_L; HD++)
		momentvec[m+HD][samp] += thisV[0]*thisV[HD];
	m += genome_L+1;

	// VECTOR: cross-correlations between virus subpopulations
	for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
			it != sv.V.end(); ++it)
	{
		for (std::map<std::vector<int>,int>::iterator itp=sv.V.begin();
				itp != sv.V.end(); ++itp)
		{
			int HD = 0;

			for (unsigned int i=0; i<genome_L; i++)
				if ((it->first)[i] != (itp->first)[i])
					HD++;

			momentvec[m+HD][samp] += (double)(it->second)*(double)(itp->second);
		}
	}
	m += genome_L+1;


	///// Probabilistic quantities:
	int p = 0; // probability sub-vector index

	// VECTOR: Sequence drawing probabilities of HD-unique virus subpopulations
	if (sv.NV>=1) {
		for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
				it != sv.V.end(); ++it)
		{
			int HD = 0;
			for (unsigned int i=0; i<genome_L; i++)
				if ((it->first)[i] != seq0[i])
					HD++;

			momentvec[m+HD][samp] += (double)(it->second)/(double)sv.NV;
		}
		for (unsigned int HD=0; HD<=genome_L; HD++)
			jointpcount[p+HD][samp]++;
	}
	m += genome_L+1;
	p += genome_L+1;

	// VECTOR: Sequence drawing probabilities of virus sequence pair separations
	if (sv.NV>=1) {
		for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
				it != sv.V.end(); ++it)
		{
			for (std::map<std::vector<int>,int>::iterator itp=sv.V.begin();
					itp != sv.V.end(); ++itp)
			{
				int HD = 0;

				for (unsigned int i=0; i<genome_L; i++)
					if ((it->first)[i] != (itp->first)[i])
						HD++;

				momentvec[m+HD][samp] += (double)(it->second)*(double)(itp->second)
					/pow((double)sv.NV,2.0);
			}
		}
		for (unsigned int HD=0; HD <= genome_L; HD++)
			jointpcount[p+HD][samp]++;
	}
	m += genome_L+1;
	p += genome_L+1;

	// SCALAR: diversity_V
	if (sv.NV>=1) {
		double thisdiv = 0.0;
		
		for (std::map<std::vector<int>,int>::iterator it=sv.V.begin();
				it != sv.V.end(); ++it)
			thisdiv += (double)it->second*(double)it->second;
		thisdiv = (double)sv.NV*(double)sv.NV/thisdiv;

		momentvec[m++][samp] += thisdiv;
		jointpcount[p++][samp]++;
	} 

}


// Write averaged results to file:
void dump_results(unsigned int Ntraj, double samp_dt, unsigned int genome_L,
		const std::vector<std::vector<double> > & momentvec,
		const std::vector<std::vector<unsigned long> > & jointpcount,
		char *ofilename,
		time_t start_time, time_t end_time)
{
	unsigned int Nmoments = momentvec.size();
	unsigned int Nsamples = momentvec[0].size();
	unsigned int Nprobs = jointpcount.size();

	// Calculate total simulation time:
	int time_h = (end_time-start_time)/3600;
	int time_m = ((end_time-start_time)%3600)/60;
	int time_s = (end_time-start_time)%60;

	// Open output file:
	std::ofstream of;
	of.open(ofilename);

	// Write header:
	of << "# Gillepie simulation of HIV infection." << std::endl
		<< "#" << std::endl
		<< "# Simulation required "
		<< time_h << " hours, "
		<< time_m << " minutes and "
		<< time_s << " seconds to complete." << std::endl
		<< "#" << std::endl
		<< "t NX NY NV NVsq";
	for (unsigned int i=0; i<=genome_L; i++)
		of << " V" << i;
	for (unsigned int i=0; i<=genome_L; i++)
		of << " V0V" << i;
	for (unsigned int i=0; i<=genome_L; i++)
		of << " VV" << i;
	for (unsigned int i=0; i<=genome_L; i++)
		of << " P" << i;
	for (unsigned int i=0; i<=genome_L; i++)
		of << " PP" << i;
	of << " div_V" << std::endl;

	// Switch to scientific notation:
	//of << std::scientific;

	// Write normalized moments:
	for (unsigned int s=0; s<Nsamples; s++) {
		of << samp_dt*(double)s << " ";
		for (unsigned int m=0; m<Nmoments-Nprobs; m++) {
			of << momentvec[m][s] / (double)Ntraj << " ";
		}
		for (unsigned int p=0; p<Nprobs; p++)
			of << momentvec[Nmoments-Nprobs+p][s] / (double)jointpcount[p][s] << " ";
		of << std::endl;
	}

	// Close file:
	of.close();
}


int main (int argc, char **argv)
{

	// Parse command line arguments:
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " outfile" << std::endl;
		exit(0);
	}
	char *ofilename = argv[1];

	// Simulation parameters:

	unsigned int Nsamples = 1001;
	unsigned int Ntraj = 1024;

	// Model parameters:
	unsigned int genome_L = 5;
	unsigned int genome_base = 4;
	double T = 30.0;

	double r_xprod = 1e2;
	double r_inf = 0.01;
	double r_vprod = 2.0;
	double r_xdeath = 0.1;
	double r_ydeath = 0.5;
	double r_vdeath = 5.0;
	double mutprob_rt = 0.2;
	double mutprob_rnap = 0.0;

	std::vector<int> seq0(genome_L, 0);
	StateVec n0;
	n0.X = 1000;
	n0.V[seq0] = 10;

	// Derived simulation parameters:
	double samp_dt = T/((double)(Nsamples-1));

	// Define sampling vectors:

	unsigned int Nmoments = 5 + 5*(1+genome_L);
	unsigned int Nprobs = 1 + 2*(1+genome_L); // Probabilistic moments

	std::vector<std::vector<double> > momentvec;
	momentvec.resize(Nmoments);
	for (unsigned int i=0; i<Nmoments; i++)
		momentvec[i].assign(Nsamples, 0.0);

	// Fire up MPI:
	MPI::Init(argc, argv);
	int mpi_size = MPI::COMM_WORLD.Get_size();
	int mpi_rank = MPI::COMM_WORLD.Get_rank();
	atexit(MPI::Finalize);

	// Determine chunk size:
	unsigned int chunk = Ntraj/mpi_size;
	if (mpi_size*chunk < Ntraj)
		chunk++;
	Ntraj = mpi_size*chunk;

	// Initialise PRNG:
	unsigned short *seedbuff = new unsigned short [3];
	seedbuff[0] = 42;
	seedbuff[1] = 53;
	seedbuff[2] = mpi_rank;

	// Define variable to track number of trajectories used
	// in joint probability calculations:
	std::vector<std::vector<unsigned long> > jointpcounts;
	jointpcounts.resize(Nprobs);
	for (unsigned int p=0; p<Nprobs; p++)
		jointpcounts[p].assign(Nsamples, 0);

	// Start timer:
	time_t start_time = std::time(NULL);

	// Loop over different trajectories:
	for (unsigned int traj=0; traj<chunk; traj++) {

		// Initialise trajectory state and time:
		StateVec n = n0;
		double t = 0.0;

		// Perform initial sample:
		unsigned int samp = 0;
		sample(genome_L, jointpcounts, samp++, n, momentvec);

		////////////////////////////////////////////
		// BEGIN SIMULATION LOOP

		while (true) {

			// Calculate rates:

			std::vector<double> R(6,0.0);
			R[0] = r_xprod;
			R[1] = r_inf*n.get_inf_rate();
			R[2] = r_vprod*n.get_vprod_rate();
			R[3] = r_xdeath*n.get_xdeath_rate();
			R[4] = r_ydeath*n.get_ydeath_rate();
			R[5] = r_vdeath*n.get_vdeath_rate();

			double Rtot = R[0]+R[1]+R[2]+R[3]+R[4]+R[5];


			// Determine time of next reaction:

			t += -log(erand48(seedbuff))/Rtot; 

			//std::cout << "t=" << t << " ";


			// Sample and/or break if necessary:

			while (samp_dt*(double)samp <= fmin(t,T)) {
				std::cout << "Rank " << mpi_rank
					<< " Traj " << traj+1 << " of " << chunk
					<< " t=" << samp_dt*(double)samp << std::endl;
				sample(genome_L, jointpcounts, samp++, n, momentvec);
			}

			if (t>=T) {
				// Dump final configuration:
				n.calc_totals();
				std::cout << "traj " << traj
					<< " final conf: X=" << n.X << " Y=" << n.NY << " V=" << n.NV
					<< std::endl;
				break;
			}


			// Determine reaction to implement:

			double u = Rtot*erand48(seedbuff);
			int react=0;
			double R_cumul = 0.0;
			while ((R_cumul+=R[react]) < u)
				react++;


			// Implement chosen reaction:

			std::vector<int> seq;
			switch(react) {
				case 0:
					// X production
					//std::cout << "X production" << std::endl;

					n.X++;

					break;

				case 1:
					// Infection
					//std::cout << "Infection ";

					// Remove uninfected cell
					n.X--;

					// Remove random virion
					seq = n.choose_V(seedbuff);
					n.dec_V(seq);

					if (erand48(seedbuff)<mutprob_rt) {
						// RT mutation
						unsigned int rnd = nrand48(seedbuff) % (genome_L*(genome_base-1));
						unsigned int loc = rnd/(genome_base-1);
						unsigned int rot = rnd%(genome_base-1) + 1;
						seq[loc] = (seq[loc]+rot)%genome_base;

						//std::cout << "with mutation";
					}

					//std::cout << std::endl;

					// Increment corresponding infected cell count
					n.inc_Y(seq);

					break;

				case 2:
					// V production
					//std::cout << "V production ";

					// Choose random infected cell
					seq = n.choose_Y(seedbuff);

					if (erand48(seedbuff)<mutprob_rnap) {
						// RNAP mutation
						unsigned int rnd = nrand48(seedbuff) % (genome_L*(genome_base-1));
						unsigned int loc = rnd/(genome_base-1);
						unsigned int rot = rnd%(genome_base-1) + 1;
						seq[loc] = (seq[loc]+rot)%genome_base;

						//std::cout << "with mutation";
					}

					//std::cout << std::endl;

					// Increment corresponding virion count
					n.inc_V(seq);

					break;

				case 3:
					// X death
					//std::cout << "X death" << std::endl;

					n.X--;

					break;

				case 4:
					// Y death
					//std::cout << "Y death" << std::endl;

					// Remove random infected cell
					n.dec_Y(n.choose_Y(seedbuff));

					break;

				case 5:
					// V death
					//std::cout << "V death" << std::endl;

					// Remove random virion
					n.dec_V(n.choose_V(seedbuff));

					break;
			}

		}

		// END SIMULATION LOOP
		////////////////////////////////////////////

		//n.dump_vseqs();
	}

	// End timer:
	time_t end_time = time(NULL);

	// Collate results from MPI processes:

	if (mpi_rank == 0) {

		// This is the root node. Aggregate results:

		double *recv_moment = new double [Nsamples];
		unsigned long *recv_jointpcounts = new unsigned long [Nsamples];
		for (int recv_rank = 1; recv_rank<mpi_size; recv_rank++) {

			for (unsigned int m=0; m<Nmoments; m++) {
				MPI::COMM_WORLD.Recv(recv_moment, Nsamples, MPI::DOUBLE, recv_rank, m);

				for (unsigned int i=0; i<Nsamples; i++)
					momentvec[m][i] += recv_moment[i];
			}
			for (unsigned int p=0; p<Nprobs; p++) {
				MPI::COMM_WORLD.Recv(recv_jointpcounts, Nsamples, MPI::UNSIGNED_LONG,
						recv_rank, Nmoments+p);

				for (unsigned int i=0; i<Nsamples; i++)
					jointpcounts[p][i] += recv_jointpcounts[i];
			}
		}

		delete [] recv_moment;
		delete [] recv_jointpcounts;


		// Write normalized results to file:
		dump_results(Ntraj, samp_dt, genome_L,
				momentvec, jointpcounts, ofilename, start_time, end_time);


	} else {
		// Send results to root node:

		for (unsigned int m=0; m<Nmoments; m++)
			MPI::COMM_WORLD.Send(&(momentvec[m][0]), Nsamples, MPI::DOUBLE, 0, m);

		for (unsigned int p=0; p<Nprobs; p++)
			MPI::COMM_WORLD.Send(&(jointpcounts[p][0]), Nsamples, MPI::UNSIGNED_LONG, 0, Nmoments+p);
	}


	// Success in Unix-ese:
	exit(0);
}
