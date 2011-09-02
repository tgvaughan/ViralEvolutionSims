// Definition of various classes used by virus_gencor_tauleap.cc.

// Class for genetic sequences:
class Sequence : public std::vector<int> {
	public:

		// Number of characters in sequence:
		int nChar;

		// Vector of neighbouring sequences:
		std::vector<Sequence> neighbours;

		// Flag to mark whether neighbours has been populated:
		bool neighboursPopulated;

		// Constructors:
		Sequence (int length, int p_nChar) : std::vector<int> (length, 0) {
			nChar = p_nChar;
			neighboursPopulated = false;
		}
		Sequence () {
			neighboursPopulated = false;
		};
		Sequence (const Sequence & s) : std::vector<int> (s) {

			nChar = s.nChar;
			/*
			resize(s.size());
			for (int i=0; i<s.size(); i++)
				operator[](i) = s[i];
			*/

			neighboursPopulated = false;
		}

		// Populate vector of neighbouring sequences:
		void genNeighbours() {

			neighbours.resize(size()*(nChar-1));

			int idx=0;
			for (int i=0; i<size(); i++) {
				for (int j=1; j<nChar; j++) {
					Sequence a = *this;
					a[i] = (a[i]+j) % nChar;
					neighbours[idx++] = a;
				}
			}
		}

		// Return vector containing neighbouring sequences:
		std::vector<Sequence> getNeighbours() {

			if (!neighboursPopulated)
				genNeighbours();

			return neighbours;
		}

		// Return the number of neighbouring sequences:
		int getNumNeighbours () {
			return size()*(nChar - 1);
		}

		// Return a particular neighbouring sequence:
		Sequence getNeighbour (int n) {

			// Copy current sequence:
			Sequence a = *this;

			// Extract locus and character shift (offset) from n:
			int locus = n/(nChar - 1);
			int offset = n%(nChar - 1) + 1;

			// Modify character at locus:
			a[locus] = (a[locus]+offset) % nChar;

			// Return neighbour sequence:
			return a;
		}

		// Select and return a random neighbouring sequence:
		Sequence chooseNeighbour (unsigned short *buf) {

			return getNeighbour(nrand48(buf)%getNumNeighbours());
		}
};

// Function to handle stream output of sequence objects:
std::ostream & operator<< (std::ostream & output, const Sequence & seq)
{
	output << "[";

	for (int i=0; i<seq.size(); i++) {

//		if (i>0)
//			output << ",";

		if (seq[i]>0)
			output << seq[i];
		else
			output << '.';
	}

	output << "]";

	return output;
}


// Class for genetically diverse populations
class GenPopulation {
	public:

		std::map<Sequence, double> pop;

		// Constructors:
		GenPopulation() {};
		GenPopulation(const GenPopulation & p) {
			pop = p.pop;
		}

		// Assignment:
		GenPopulation operator= (GenPopulation arg) {
			pop = arg.pop;

			return *this;
		}

		// Return iterators for pop:
		std::map<Sequence, double>::iterator begin() {
			return pop.begin();
		}
		std::map<Sequence, double>::iterator end() {
			return pop.end();
		}

		// Element indexing operator:
		double operator[] (Sequence s) {
			std::map<Sequence, double>::iterator it = pop.find(s);
			if (it == pop.end())
				return 0;

			return it->second;
		}

		// Negativity check:
		bool isnegative () {

			std::map<Sequence, double>::iterator it;

			for (it = pop.begin(); it != pop.end(); it++)
				if (it->second < 0)
					return true;

			return false;
		}

		// Get total size:
		double popSize() {

			std::map<Sequence, double>::iterator it;

			double size = 0.0;

			for (it = pop.begin(); it != pop.end(); it++)
				size += it->second;

			return size;
		}

		// Add particles to population with zero-removal:
		double popAdd (double n, Sequence s) {
			double newval = pop[s] + n;
			if (newval != 0)
				pop[s] = newval;
			else
				pop.erase(s);

			return newval;
		}
};


// Class for genetically homogeneous populations
class NonGenPopulation {
	public:

		double n;

		// Constructors:
		NonGenPopulation (double p_n) {
			n = p_n;
		}
		NonGenPopulation () {
		}
		NonGenPopulation (const NonGenPopulation & p) {
			n = p.n;
		}

		// Assignment:
		NonGenPopulation operator= (const NonGenPopulation & arg) {
			n = arg.n;

			return *this;
		}

		// Negativity check:
		bool isnegative () {
			return n<0;
		}

		// Get total population size:
		double popSize() {
			return n;
		}
};


// Class for system state vectors
class StateVec {
	public:

		std::vector<NonGenPopulation> nonGenetic;
		std::vector<GenPopulation> genetic;

		// Constructors:
		StateVec(int nNonGeneticSpecies, int nGeneticSpecies)
		{
			nonGenetic.resize(nNonGeneticSpecies);
			genetic.resize(nGeneticSpecies);
		}
		StateVec() {};
		StateVec(const StateVec & arg) {
			nonGenetic = arg.nonGenetic;
			genetic = arg.genetic;
		}

		// Check for negative populations:
		bool isnegative() {

			for (int i=0; i<nonGenetic.size(); i++)
				if (nonGenetic[i].isnegative())
					return true;

			for (int i=0; i<genetic.size(); i++)
				if (genetic[i].isnegative())
					return true;

			return false;
		}

};


// Class describing individual reactions:
class Reaction {
	public:

		double rate;
		std::vector<int> inNonGen, inGen;
		std::vector<int> outNonGen, outGen;
		std::vector<bool> mutate;

		// Critical reaction parameter:
		double alpha;

		// Does reaction involve genetically diverse populations?
		bool isGen;
		bool isMutation;

		// Genetic reaction propensities:
		std::map<Sequence, double> propensities;
		Sequence criticalSeq;

		// Nongenetic reaction propensity:
		double propensity;

		// Time of next critical reaction:
		double criticalDelta;

		// Constructors:
		Reaction(std::vector<int> p_inNonGen,
				std::vector<int> p_inGen,
				std::vector<int> p_outNonGen,
				std::vector<int> p_outGen,
				std::vector<bool> p_mutate,
				double p_alpha, double p_rate)
		{
			inNonGen = p_inNonGen;
			inGen = p_inGen;
			outNonGen = p_outNonGen;
			outGen = p_outGen;
			mutate = p_mutate;
			rate = p_rate;

			isGen = false;
			for (int i=0; i<inGen.size(); i++) {
				if (inGen[i]>0 || outGen[i]>0)
					isGen = true;
			}

			isMutation = false;
			for (int i=0; i<mutate.size(); i++)
				if (mutate[i])
					isMutation = true;

			alpha = p_alpha;

		}
		Reaction() {};

		// Calculate reaction propensities for given state,
		// determine critical reactions and time and sequence
		// corresponding to next critical reaction:
		void calcPropensity(StateVec & x, double delta, unsigned short *buf)
		{
			// Reset critical state:
			criticalDelta = -1.0;

			if (isGen) {

				// Empty propensity and critical sequence containers:
				propensities.clear();
				criticalSeq.clear();

				// Determine reactant population with smallest diversity:
				// (Rationale: the subpopulation of each genotype must be
				// non-zero for each of the reactant species, so it makes sense
				// to consider only those genotypes extant in the reactant
				// population with the smallest diversity.)
				int mindiversity = 0;
				int imin = -1;
				for (int i=0; i<x.genetic.size(); i++) {
					if (inGen[i]>0) {
						if (imin<0 || x.genetic[i].pop.size()<mindiversity) {
							mindiversity = x.genetic[i].pop.size();
							imin = i;
						}
					}
				}

				// Consider each extant genotype separately:
				std::map<Sequence, double>::iterator it;
				for (it = x.genetic[imin].begin(); it != x.genetic[imin].end(); it++) {

					Sequence thisSeq = it->first;

					bool crit = false;

					// Determine reaction propensity:
					double a = rate;
					/*
					if (isMutation)
						a *= thisSeq.size()*3;
					*/

					for (int i=0; i<x.genetic.size(); i++) {
						for (int m=0; m<inGen[i]; m++)
							a *= x.genetic[i][thisSeq] - m;
					}
					for (int i=0; i<x.nonGenetic.size(); i++) {
						for (int m=0; m<inNonGen[i]; m++)
							a *= x.nonGenetic[i].n - m;
					}

					//if (isMutation)
					//	std::cout << "Mutation process rate = " << a << std::endl;

					if (a>0) {

						// Determine critical reaction number:
						//double Nc = delta*a + alpha*sqrt(delta*a);
						double Nc = 1000.0;

						// Check for criticality:
						for (int i=0; i<x.genetic.size(); i++) {
							if (!mutate[i]) {
								if (x.genetic[i][thisSeq] < Nc*(inGen[i]-outGen[i]))
									crit = true;
							} else {
								if (x.genetic[i][thisSeq] < Nc*inGen[i])
									crit = true;
							}
						}
						for (int i=0; i<x.nonGenetic.size(); i++) {
							if (x.nonGenetic[i].n < Nc*(inNonGen[i] - outNonGen[i]))
									crit = true;
						}

						if (crit) {

							// Critical reaction: choose reaction time

							if (isMutation) {

								for (int i=0; i<thisSeq.size()*(thisSeq.nChar-1); i++) {
									
									double thisCriticalDelta = -log(erand48(buf))/a;
									if (criticalDelta < 0 || thisCriticalDelta < criticalDelta) {
										criticalDelta = thisCriticalDelta;
										criticalSeq = thisSeq;
									}
								}

							} else {

								double thisCriticalDelta = -log(erand48(buf))/a;
								if (criticalDelta < 0 || thisCriticalDelta < criticalDelta) {
									criticalDelta = thisCriticalDelta;
									criticalSeq = thisSeq;
								}
							}

						} else {

							// Non-critical reaction: record propensity
							propensities[it->first] = a;
						}
					}
				}

			} else {

				bool crit = false;

				// Determine reaction propensity:
				double a = rate;
				for (int i=0; i<x.nonGenetic.size(); i++) {
					for (int m=0; m<inNonGen[i]; m++) {
						a *= x.nonGenetic[i].n - m;
					}
				}

				propensity = 0.0;

				if (a>0) {

					// Determine critical reaction number:
					//double Nc = a*delta + alpha*sqrt(a*delta);
					double Nc = 1000.0;

					// Check for criticality
					for (int i=0; i<x.nonGenetic.size(); i++) {
						if (x.nonGenetic[i].n < Nc*(inNonGen[i] - outNonGen[i]))
							crit = true;
					}

					if (crit) {

						// Critical reaction: choose reaction time
						criticalDelta = -log(erand48(buf))/a;

					} else {

						// Non-critical reaction: record propensity
						propensity = a;
					}
				}
			}
		}

		// Check whether state is within Nc reactions of bottoming out:
		bool isCritical ()
		{
			return criticalDelta > 0.0;
		}

		// Perform tau-leaping integration step:
		void tauLeap(StateVec & x, double dt, unsigned short *buf)
		{

			if (isGen) {

				std::map<Sequence, double>::iterator it;
				for (it = propensities.begin(); it != propensities.end(); it++) {

					Sequence thisSeq = it->first;
					double thisProp = it->second;

					if (isMutation) {

						// Adjust propensity to account for number of neighbouring sequences:
						thisProp /= thisSeq.getNumNeighbours();

						for (int neighbourNum = 0; neighbourNum < thisSeq.getNumNeighbours(); neighbourNum++) {

							double nreacts = poissonian(thisProp*dt, buf);

							if (nreacts>0) {

								Sequence mutantSeq = thisSeq.getNeighbour(neighbourNum);

								for (int i=0; i<x.genetic.size(); i++) {
									if (!mutate[i]) {
										if (outGen[i]-inGen[i] != 0)
											x.genetic[i].popAdd(nreacts*(outGen[i]-inGen[i]), thisSeq);
									} else {
										x.genetic[i].popAdd(-nreacts*inGen[i], thisSeq);
										x.genetic[i].popAdd(nreacts*outGen[i], mutantSeq);
									}
								}
								for (int i=0; i<x.nonGenetic.size(); i++) {
									if (outNonGen[i]-inNonGen[i] != 0)
										x.nonGenetic[i].n += nreacts*(outNonGen[i]-inNonGen[i]);
								}
							}
						}

					} else {

						double nreacts = poissonian(thisProp*dt, buf);

						if (nreacts>0) {
							for (int i=0; i<x.genetic.size(); i++) {
								if ((outGen[i] - inGen[i]) != 0)
									x.genetic[i].popAdd(nreacts*(outGen[i]-inGen[i]), thisSeq);
							}
							for (int i=0; i<x.nonGenetic.size(); i++) {
								if ((outNonGen[i] - inNonGen[i]) != 0)
									x.nonGenetic[i].n += nreacts*(outNonGen[i]-inNonGen[i]);
							}
						}

					}

				}

			} else {

				double nreacts = poissonian(propensity*dt, buf);

				for (int i=0; i<x.nonGenetic.size(); i++) {
					if (outNonGen[i]-inNonGen[i] != 0) {
						x.nonGenetic[i].n += nreacts*(outNonGen[i]-inNonGen[i]);
					}
				}

			}
		}

		// Implement critical reaction on given state:
		void implementCritical(StateVec & x, unsigned short *buf)
		{

			if (isGen) {

				if (isMutation) {

					Sequence mutantSeq = criticalSeq.chooseNeighbour(buf);

					for (int i=0; i<x.genetic.size(); i++) {
						if (mutate[i]) {
							if (inGen[i] > 0)
								x.genetic[i].pop[criticalSeq] -= inGen[i];
							if (outGen[i] > 0)
								x.genetic[i].pop[mutantSeq] += outGen[i];
						} else {
							if (outGen[i]-inGen[i] != 0)
								x.genetic[i].pop[criticalSeq] += outGen[i]-inGen[i];
						}
					}
					for (int i=0; i<x.nonGenetic.size(); i++)
						x.nonGenetic[i].n += outNonGen[i]-inNonGen[i];

				} else {

					for (int i=0; i<x.genetic.size(); i++) {
						if (outGen[i]-inGen[i] != 0)
							x.genetic[i].pop[criticalSeq] += outGen[i]-inGen[i];
					}
					for (int i=0; i<x.nonGenetic.size(); i++)
						x.nonGenetic[i].n += outNonGen[i]-inNonGen[i];
				}

			} else {

				for (int i=0; i<x.nonGenetic.size(); i++)
					x.nonGenetic[i].n += outNonGen[i]-inNonGen[i];
			}

		}
};


// Class for sampling moments.
class Moment {
	public:

		std::string name;
		int Nsamples;
		double *mean;
		double *half_mean;
		double *var;

		double *tmp_mean;
		double *tmp_half_mean;
		double *tmp_var;

		double *ftse;
		double *sem;

		double (*samplefunc)(StateVec);

		bool initialised;

		// Constructor:
		Moment(int p_Nsamples, double (*p_samplefunc)(StateVec), std::string p_name)
		{
			Nsamples = p_Nsamples;
			samplefunc = p_samplefunc;
			name = p_name;

			mean = new double[Nsamples];
			half_mean = new double[Nsamples];
			var = new double[Nsamples];

			tmp_mean = new double[Nsamples];
			tmp_half_mean = new double[Nsamples];
			tmp_var = new double[Nsamples];

			ftse = new double[Nsamples];
			sem = new double[Nsamples];

			for (int s=0; s<Nsamples; s++) {
				mean[s] = 0.0;
				half_mean[s] = 0.0;
				var[s] = 0.0;
			}

			initialised = true;

		}

		// Default constructor:
		Moment()
		{
			initialised = false;
		}

		// Destructor:
		~Moment() {
			if (initialised) {
				delete [] mean;
				delete [] half_mean;
				delete [] var;

				delete [] tmp_mean;
				delete [] tmp_half_mean;
				delete [] tmp_var;

				delete [] ftse;
				delete [] sem;
			}
		}

		// Record moment sample to temporary array:
		void record(StateVec x, int samp) {
			tmp_mean[samp] = (*samplefunc)(x);
			tmp_var[samp] = tmp_mean[samp]*tmp_mean[samp];
		}

		// Record half-step moment to temporary array:
		void record_half(StateVec x, int samp) {
			tmp_half_mean[samp] = (*samplefunc)(x);
		}

		// Incorporate temporary array contents into moment calculation:
		void add() {
			for (int s=0; s<Nsamples; s++) {
				mean[s] += tmp_mean[s];
				half_mean[s] += tmp_half_mean[s];
				var[s] += tmp_var[s];
			}
		}

		// Transmit moment calculation results to root mpi node:
		void mpi_send(int *tag) {
			MPI::COMM_WORLD.Send(mean, Nsamples, MPI::DOUBLE, 0, (*tag)++);
			MPI::COMM_WORLD.Send(half_mean, Nsamples, MPI::DOUBLE, 0, (*tag)++);
			MPI::COMM_WORLD.Send(var, Nsamples, MPI::DOUBLE, 0, (*tag)++);
		}

		// Receive moment calculation results from non-root mpi node:
		void mpi_receive(int *tag, int from_rank) {
			MPI::COMM_WORLD.Recv(tmp_mean, Nsamples, MPI::DOUBLE, from_rank, (*tag)++);
			MPI::COMM_WORLD.Recv(tmp_half_mean, Nsamples, MPI::DOUBLE, from_rank, (*tag)++);
			MPI::COMM_WORLD.Recv(tmp_var, Nsamples, MPI::DOUBLE, from_rank, (*tag)++);

			for (int s=0; s<Nsamples; s++) {
				mean[s] += tmp_mean[s];
				half_mean[s] += tmp_half_mean[s];
				var[s] += tmp_var[s];
			}
		}

		// Post-processing of moment calculations:
		void normalise(int Ntraj) {

			for (int s=0; s<Nsamples; s++) {
				mean[s] /= Ntraj;
				half_mean[s] /= Ntraj;
				var[s] = var[s]/Ntraj - mean[s]*mean[s];
				sem[s] = sqrt(var[s]/Ntraj);
				ftse[s] = fabs(mean[s]-half_mean[s]);
			}

		}

};
