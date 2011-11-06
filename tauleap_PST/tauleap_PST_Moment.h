/*
 * tauleap_PST_Moment.h
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_MOMENT_H_
#define TAULEAP_PST_MOMENT_H_

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


class MomentVector {
public:

	int Nsamples;
	std::vector<double> (*samplefunc)(const StateVec &);
	std::string name;

	std::vector<double> mean;
	std::vector<double> var;
	std::vector<double> sem;

	int sequenceL;

	MomentVector(int p_Nsamples,
			std::vector<double>(*p_samplefunc)(const StateVec &),
			std::string p_name, int p_sequenceL);
	MomentVector();

	void sample(const StateVec & sv, int samp);

	void mpi_send(int mpi_rank, int & tag);
	void mpi_recv(int recv_rank, int & tag);

	void normalise(int Npaths);

};

#endif /* TAULEAP_PST_MOMENT_H_ */
