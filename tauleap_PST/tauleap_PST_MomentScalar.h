/*
 * tauleap_PST_MomentScalar.h
 *
 *  Created on: 28/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_MOMENTSCALAR_H_
#define TAULEAP_PST_MOMENTSCALAR_H_

class MomentScalar {
public:

	int Nsamples;
	double (*samplefunc)(const StateVec &);
	std::string name;

	std::vector<double> mean;
	std::vector<double> var;
	std::vector<double> sem;

	MomentScalar(int p_Nsamples, double (*p_samplefunc)(const StateVec &), std::string p_name);
	MomentScalar();

	void sample(const StateVec & sv, int samp);

	void mpi_send(int mpi_rank, int & tag);
	void mpi_recv(int recv_rank, int & tag);

	void normalise(int Npaths);
};

#endif /* TAULEAP_PST_MOMENTSCALAR_H_ */
