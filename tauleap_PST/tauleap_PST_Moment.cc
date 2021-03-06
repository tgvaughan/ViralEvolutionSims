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
#include "tauleap_PST_Moment.h"

Moment::Moment(int p_Nsamples,
		void (*p_samplefunc)(const StateVec &, std::vector<double> &),
		std::string p_name, std::vector<int> & p_dim) {

	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;
	dim = p_dim;

	length = 1;
	for (unsigned int i=0; i<dim.size(); i++)
		length *= dim[i];

	res.resize(length);

	mean.resize(Nsamples*length,0);
	var.resize(Nsamples*length,0);
	sem.resize(Nsamples*length,0);
}

Moment::Moment() { };

/**
 * Sample moment
 */
void Moment::sample(const StateVec & sv, int samp) {

	// Collect sample in res:
	(*samplefunc)(sv, res);

	for (int i=0; i<length; i++) {
		int idx = length*samp + i;

		mean[idx] += res[i];
		var[idx] += res[i]*res[i];
	}
}

/**
 * Send sampled data to root node
 */
void Moment::mpi_send(int mpi_rank, int & tag) {

	MPI::COMM_WORLD.Send(&(mean[0]), mean.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(var[0]), var.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(sem[0]), sem.size(), MPI_DOUBLE, 0, tag++);

}

/**
 * Obtain sampled data from non-root nodes
 */
void Moment::mpi_recv(int recv_rank, int & tag) {

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
void Moment::normalise (int Npaths) {

	for (unsigned int i=0; i<mean.size(); i++) {
		mean[i] /= Npaths;
		var[i] = var[i]/Npaths - mean[i]*mean[i];
		sem[i] = sqrt(var[i]/Npaths);
	}
}
