/*
 * tauleap_PST_MomentVector.cc
 *
 *  Created on: 28/11/2011
 *      Author: Tim Vaughan
 */

#include <vector>
#include <cmath>

#include <mpi.h>

#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_MomentVector.h"

MomentVector::MomentVector(int p_Nsamples,
		void (*p_samplefunc)(const StateVec &, std::vector<double> &),
		std::string p_name, int p_length) {

	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;
	length = p_length;

	res.resize(length);

	mean.resize(Nsamples*length,0);
	var.resize(Nsamples*length,0);
	sem.resize(Nsamples*length,0);
}

MomentVector::MomentVector() { };

/**
 * Sample moment
 */
void MomentVector::sample(const StateVec & sv, int samp) {

	// Collect sample in res:
	(*samplefunc)(sv, res);

	for (int h=0; h<length; h++) {
		int idx = length*samp + h;

		mean[idx] += res[h];
		var[idx] += res[h]*res[h];
	}
}

/**
 * Send sampled data to root node
 */
void MomentVector::mpi_send(int mpi_rank, int & tag) {

	MPI::COMM_WORLD.Send(&(mean[0]), mean.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(var[0]), var.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(sem[0]), sem.size(), MPI_DOUBLE, 0, tag++);

}

/**
 * Obtain sampled data from non-root nodes
 */
void MomentVector::mpi_recv(int recv_rank, int & tag) {

	std::vector<double> recv_mean(mean.size());
	std::vector<double> recv_var(var.size());
	std::vector<double> recv_sem(sem.size());

	MPI::COMM_WORLD.Recv(&(recv_mean[0]), mean.size(), MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(&(recv_var[0]), var.size(), MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(&(recv_sem[0]), sem.size(), MPI_DOUBLE, recv_rank, tag++);

	for (unsigned int i=0; i<mean.size(); i++) {
		mean[i] += recv_mean[i];
		var[i] += recv_var[i];
		sem[i] += recv_sem[i];
	}
}

/**
 * Perform post-processing of moment data
 */
void MomentVector::normalise (int Npaths) {

	for (unsigned int i=0; i<mean.size(); i++) {
		mean[i] /= Npaths;
		var[i] = var[i]/Npaths - mean[i]*mean[i];
		sem[i] = sqrt(var[i]/Npaths);
	}
}
