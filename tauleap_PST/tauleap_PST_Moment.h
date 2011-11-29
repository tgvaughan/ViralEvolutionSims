/*
 * tauleap_PST_Moment.h
 *
 *  Created on: 28/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_MOMENT_H_
#define TAULEAP_PST_MOMENT_H_

class Moment {
public:

	int Nsamples;
	void (*samplefunc)(const StateVec &, std::vector<double> &);
	std::string name;

	std::vector<double> res;

	std::vector<double> mean;
	std::vector<double> var;
	std::vector<double> sem;

	std::vector<int> dim;
	int length;

	Moment(int p_Nsamples,
			void (*p_samplefunc)(const StateVec &, std::vector<double> &),
			std::string p_name, std::vector<int> & p_dim);
	Moment();

	void sample(const StateVec & sv, int samp);

	void mpi_send(int mpi_rank, int & tag);
	void mpi_recv(int recv_rank, int & tag);

	void normalise(int Npaths);

};

#endif /* TAULEAP_PST_MOMENT_H_ */
