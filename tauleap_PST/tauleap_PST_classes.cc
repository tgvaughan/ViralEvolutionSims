/*
 * tauleap_PST_classes.cc
 *
 *  Created on: 28/10/2011
 *      Author: Tim Vaughan
 */

#include <iostream>
#include <string>
#include <vector>

#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "poissonian.h"
#include "tauleap_PST_classes.h"


// StateVec member functions:

StateVec::StateVec (int length) {
	L = length;

	X = 0.0;
	Y.resize(L+1,0.0);
	V.resize(L+1,0.0);
}

StateVec::StateVec (const StateVec & src) {
	L = src.L;

	X = src.X;
	Y = src.Y;
	V = src.V;
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;

	X = src.X;
	Y = src.Y;
	V = src.V;

	return *this;
}


// Reaction member functions:

Reaction::Reaction(int p_inX, int p_inY, int p_inV,
		int p_outX, int p_outY, int p_outV,
		bool p_mutY, bool p_mutV,
		double p_rate, double p_mutrate) {

	inX = p_inX;
	inY = p_inY;
	inV = p_inV;
	outX = p_outX;
	outY = p_outY;
	outV = p_outV;

	mutY = p_mutY;
	mutV = p_mutV;

	rate = p_rate;
	mutrate = p_mutrate;

	if (mutY || mutV)
		mutation = true;
	else
		mutation = false;

	if ((inY == 0) && (inV == 0) && (outY == 0) && (outV == 0))
		onlyX = true;
	else
		onlyX = false;
}

Reaction::Reaction() { };

/**
 * Calculates the number of distinct mutations which begin at
 * a given sequence located in h1 and end with any sequence
 * located in h2.
 */
const double Reaction::get_gcond(int h1, int h2, int L)
{
	switch (h2-h1) {
		case -1:
			return h1;
		case 0:
			return 2*h1;
		case 1:
			return 3*(L-h1);
		default:
			return 0.0;
	}
}

/**
 * Calculate reaction propensities and return safe leap distance.
 *
 * tau:		desired leap distance
 * sv:		StateVec used to calculate propensities
 * buf:		RNG buffer
 *
 * returns:	safe leap distance
 */
double Reaction::getLeapDistance (double tau, double alpha, bool newcritcond, const StateVec & sv, unsigned short *buf) {

	// Calculate X portion of propensity:
	aX = rate;
	for (int m=0; m<inX; m++)
		aX *= sv.X-m;

	if (onlyX) {

		// Check for critical X reaction:
		critX = false;
		if (newcritcond) {
			double dX = tau*aX*(outX - inX);
			if ((dX<0) && (sv.X + dX - alpha*sqrt(-dX)< 0))
				critX = true;
		} else {
			if (sv.X + alpha*(outX - inX) < 0)
				critX = true;
		}

		// If critical, determine reaction time:
		if (critX) {
			double newtaucrit = -log(erand48(buf))/aX;
			if (newtaucrit<tau)
				return newtaucrit;
		}

		return tau;
	}

	double taucrit = tau;

	// Ensure propensity vectors are allocated:
	if (mutation) {
		amut.resize(3*(sv.L+1), 0.0);
		critmut.resize(3*(sv.L+1), false);
	} else {
		a.resize(sv.L+1, 0.0);
		crit.resize(sv.L+1, false);
	}

	for (int h=0; h<=sv.L; h++) {

		double atmp = aX;
		for (int m=0; m<inY; m++)
			atmp *= sv.Y[h]-m;
		for (int m=0; m<inV; m++)
			atmp *= sv.V[h]-m;

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv.L ? h+1 : sv.L); hp++) {

				int idx = 3*h + hp - h + 1;

				amut[idx] = atmp*mutrate/(3.0*sv.L)*get_gcond(h,hp,sv.L);
				if (hp==h)
					amut[idx] += atmp*(1.0-mutrate);

				// Check for critical mutation reaction:
				critmut[idx] = false;
				if (newcritcond) {
					double dY = -tau*amut[idx]*inY;
					double dV = -tau*amut[idx]*inV;
					if ((sv.Y[h] + dY - alpha*sqrt(-dY) < 0)
							|| (sv.V[h] + dV -alpha*sqrt(-dV) < 0))
						critmut[idx] = true;
				} else {
					if ((sv.Y[h] - alpha*inY < 0) || (sv.V[h] - alpha*inV < 0))
						critmut[idx] = true;
				}

				// If critical, determine reaction time:
				if (critmut[idx]) {
					double newtaucrit = -log(erand48(buf))/amut[idx];
					if (newtaucrit < taucrit) {
						taucrit = newtaucrit;
						critreact = idx;
					}
				}
			}

		} else {

			a[h] = atmp;

			// Check for critical non-mutation reaction:
			crit[h] = false;
			if (newcritcond) {
				double dY = tau*a[h]*(outY-inY);
				double dV = tau*a[h]*(outV-inV);
				if (((dY<0) && (sv.Y[h] + dY - alpha*sqrt(-dY) < 0))
						|| ((dV<0) && (sv.V[h] + dV - alpha*sqrt(-dV) < 0)))
					crit[h] = true;
			} else {
				if ((sv.Y[h] + alpha*(outY-inY) < 0) || (sv.V[h] + alpha*(outV-inV) < 0))
					crit[h] = true;
			}

			// If critical, determine reaction time:
			if (crit[h]) {

				double newtaucrit = -log(erand48(buf))/a[h];
				if (newtaucrit < taucrit) {
					taucrit = newtaucrit;
					critreact = h;
				}
			}
		}
	}

	return taucrit;
}

/**
 * Implement "next" critical reaction.
 *
 * sv_new:		StateVec to update
 */
void Reaction::doCritical(StateVec & sv_new) {

	// Implement X component
	sv_new.X += outX - inX;

	// If X is all there is, we're done:
	if (onlyX)
		return;

	if (mutation) {

		// Implement Y and V components for mutation

		int h = critreact/3;
		int hp = h + (critreact%3) - 1;

		if (mutY) {
			sv_new.Y[h] -= inY;
			sv_new.Y[hp] += outY;
			sv_new.V[h] += outV-inV;
		}
		if (mutV) {
			sv_new.Y[h] += outY-inY;
			sv_new.V[h] -= inV;
			sv_new.V[hp] += outV;
		}

	} else {

		// Implement Y and V components for non-mutation

		int h = critreact;
		sv_new.Y[h] += outY-inY;
		sv_new.V[h] += outV-inV;

	}

	return;

}

/**
 * Uses tau-leaping to implement non-critical reactions.
 *
 * tau:		length of leap
 * new_sv:	StateVec to modify
 * buf:		RNG buffer
 *
 * returns:	false if leap resulted in negative pop, true otherwise.
 */
bool Reaction::tauleap(double dt, StateVec & sv_new, unsigned short int *buf) {

	if (onlyX) {

		// Skip if reaction critical:
		if (critX)
			return true;

		// Implement reaction:
		double q = poissonian(dt*aX, buf);
		sv_new.X += q*(outX-inX);

		// Check for negative population:
		if (sv_new.X<0) {
			sv_new.X = 0;
			return true;
		}

		return false;
	}

	bool negativePop = false;

	for (int h=0; h<=sv_new.L; h++) {

		if (mutation) {

			for (int hp=(h>0 ? h-1 : 0); hp<=(h<sv_new.L ? h+1 : sv_new.L); hp++) {
				int idx = 3*h + hp - h + 1;

				// Skip if reaction critical:
				if (critmut[idx])
					continue;

				// Implement reaction:
				double q = poissonian(dt*amut[idx], buf);
				sv_new.X += q*(outX - inX);
				if (mutY) {
					sv_new.Y[h] -= q*inY;
					sv_new.Y[hp] += q*outY;
					sv_new.V[h] += q*(outV - inV);
				}
				if (mutV) {
					sv_new.Y[h] += q*(outY - inY);
					sv_new.V[h] -= q*inV;
					sv_new.V[hp] += q*outV;
				}

				// Check for negative populations:
				if (sv_new.X < 0) {
					sv_new.X = 0;
					negativePop = true;
				}
				if (sv_new.Y[h] < 0) {
					sv_new.Y[h] = 0;
					negativePop = true;
				}
				if (sv_new.V[h] < 0) {
					sv_new.V[h] = 0;
					negativePop = true;
				}
			}

		} else {

			// Skip if reaction critical:
			if (crit[h])
				continue;

			// Implement reaction:
			double q = poissonian(dt*a[h], buf);
			sv_new.X += q*(outX - inX);
			sv_new.Y[h] += q*(outY - inY);
			sv_new.V[h] += q*(outV - inV);

			// Check for negative populations:
			if (sv_new.X < 0) {
				sv_new.X = 0;
				negativePop = true;
			}
			if (sv_new.Y[h] < 0) {
				sv_new.Y[h] = 0;
				negativePop = true;
			}
			if (sv_new.V[h] < 0) {
				sv_new.V[h] = 0;
				negativePop = true;
			}
		}

	}

	return negativePop;
}

// MomentScalar member functions

MomentScalar::MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name) {
	Nsamples = p_Nsamples;
	samplefunc = p_samplefunc;
	name = p_name;

	mean.resize(Nsamples);
	var.resize(Nsamples);
	sem.resize(Nsamples);
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

	mean.resize(Nsamples*(sequenceL+1));
	var.resize(Nsamples*(sequenceL+1));
	sem.resize(Nsamples*(sequenceL+1));
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

	double *recv_mean = new double[mean.size()];
	double *recv_var = new double[var.size()];
	double *recv_sem = new double[sem.size()];

	MPI::COMM_WORLD.Recv(recv_mean, mean.size(), MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(recv_var, var.size(), MPI_DOUBLE, recv_rank, tag++);
	MPI::COMM_WORLD.Recv(recv_sem, sem.size(), MPI_DOUBLE, recv_rank, tag++);

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
void MomentVector::normalise (int Npaths) {

	for (int samp=0; samp<Nsamples; samp++) {
		for (int h=0; h<=sequenceL; h++) {
			int idx = (sequenceL+1)*samp + h;

			mean[idx] /= Npaths;
			var[idx] = var[idx]/Npaths - mean[idx]*mean[idx];
			sem[idx] = sqrt(var[idx]/Npaths);
		}
	}

}
