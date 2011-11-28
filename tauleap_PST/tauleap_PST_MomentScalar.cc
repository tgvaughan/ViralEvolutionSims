/*
 * tauleap_PST_MomentScalar.cc
 *
 *  Created on: 28/11/2011
 *      Author: Tim Vaughan
 */

#include <vector>
#include <cmath>

#include <mpi.h>

#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_MomentScalar.h"

MomentScalar::MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name) {
	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;

	mean.resize(Nsamples,0);
	var.resize(Nsamples,0);
	sem.resize(Nsamples,0);
}
MomentScalar::MomentScalar() { };

/**
 * Sample moment
 */
void MomentScalar::sample(const StateVec & sv, int samp) {
	double tmp = (*samplefunc)(sv);
	mean[samp] += tmp;
	var[samp] += tmp*tmp;
}

/**
 * Send sampled data to root node
 */
void MomentScalar::mpi_send(int mpi_rank, int & tag) {
	MPI::COMM_WORLD.Send(&(mean[0]), mean.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(var[0]), var.size(), MPI_DOUBLE, 0, tag++);
	MPI::COMM_WORLD.Send(&(sem[0]), sem.size(), MPI_DOUBLE, 0, tag++);
}

/**
 * Obtain sampled data from non-root nodes
 */
void MomentScalar::mpi_recv(int recv_rank, int & tag) {

	double *recv_mean = new double[mean.size()];
	double *recv_var = new double[var.size()];
	double *recv_sem = new double[sem.size()];

	MPI::COMM_WORLD.Recv(recv_mean, Nsamples, MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(recv_var, Nsamples, MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(recv_sem, Nsamples, MPI_DOUBLE, recv_rank, tag++);

	for (int s=0; s<Nsamples; s++) {
		mean[s] += recv_mean[s];
		var[s] += recv_var[s];
		sem[s] += recv_sem[s];
	}

	delete [] recv_mean;
	delete [] recv_var;
	delete [] recv_sem;

}

/**
 * Perform post-processing of moment data
 */
void MomentScalar::normalise (int Npaths) {
	for (int s=0; s<Nsamples; s++) {
		mean[s] /= Npaths;
		var[s] = var[s]/Npaths - mean[s]*mean[s];
		sem[s] = sqrt(var[s]/Npaths);
	}
}
