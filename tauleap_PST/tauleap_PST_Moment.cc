/*
 * tauleap_PST_Moment.cc
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#include <vector>
#include <cmath>

#include <mpi.h>

#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_Moment.h"

// MomentScalar member functions

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


// MomentVector member functions:

MomentVector::MomentVector(int p_Nsamples, std::vector<double> (*p_samplefunc)(const StateVec &),
		std::string p_name, int p_sequenceL) {

	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;
	sequenceL = p_sequenceL;

	mean.resize(Nsamples*(sequenceL+1),0);
	var.resize(Nsamples*(sequenceL+1),0);
	sem.resize(Nsamples*(sequenceL+1),0);
}

MomentVector::MomentVector() { };

/**
 * Sample moment
 */
void MomentVector::sample(const StateVec & sv, int samp) {

	std::vector<double> tmp = (*samplefunc)(sv);

	for (int h=0; h<=sequenceL; h++) {
		int idx = (sequenceL+1)*samp + h;

		mean[idx] += tmp[h];
		var[idx] += tmp[h]*tmp[h];
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
