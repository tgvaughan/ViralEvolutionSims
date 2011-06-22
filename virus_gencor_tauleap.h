// Definition of various classes used by virus_gencor_tauleap.cc.

// Class for genetic sequences:
class Sequence : public std::vector<int> {
	public:

		int nChar;

		// Constructors:
		Sequence (int length, int p_nChar) {
			nChar = p_nChar;
			resize(length, 0);  // Initial sequence filled with zeros.
		}
		Sequence () {};
		Sequence (Sequence & s) {
			nChar = s.nChar;
			for (int i=0; i<s.size(); i++)
				operator[](i) = s[i];
		}

		// Return set of neighbouring sequences:
		std::set<Sequence> getNeighbours() {

			std::set<Sequence> neighbours;

			for (int i=0; i<size(); i++) {
				for (int j=0; j<nChar-1; j++) {
					Sequence a = *this;
					a[i] = (a[i]+j) % nChar;
					neighbours.insert(a);
				}
			}

			return neighbours;
		}
};


// Abstract base class for populations
class Population {
	public:

		bool isGen;

		virtual bool isnegative() =0;
		virtual int operator[] (Sequence s) =0;

		virtual Population operator= (Population p) =0;
		virtual Population operator+ (Population p) =0;
		virtual Population operator- (Population p) =0;
		virtual Population operator* (double p) =0;
};


// Class for genetically diverse populations
class GenPopulation : public Population {
	public:

		std::map<Sequence, double> pop;
		int seqLen;

		// Constructors:
		GenPopulation(int p_seqLen) {
			isGen = true;
			seqLen = p_seqLen;
		}
		GenPopulation() {
			isGen = true;
		}
		GenPopulation(GenPopulation & p) {
			isGen = true;
			seqLen = p.seqLen;
			pop = p.pop;
		}

		// Element indexing operator:
		double operator[] (Sequence s) {
			std::map<Sequence, double>::iterator it = pop.find(s);
			if (it == pop.end())
				return 0;

			return it->second;
		}

		// Return iterators:
		map<Sequence, double>::iterator begin() {
			return pop.begin();
		}
		map<Sequence, double>::iterator end() {
			return pop.end();
		}

		// Assignment operator:
		GenPopulation operator= (GenPopulation p) {
			seqLen = p.seqLen;
			pop = p.pop;

			return *this;
		}

		//Arithmetic operators:
		GenPopulation operator+ (GenPopulation p) {

			GenPopulation res = *this;

			std::map<Sequence, double>::iterator it;

			for (it = p.pop.begin(); it != p.pop.end(); it++) {
				res.pop[it->first] += it->second;

				if (res[it->first] == 0)
					res.pop.erase(it->first);
			}

			return res;
		}

		GenPopulation operator- (GenPopulation p) {

			GenPopulation res = *this;

			std::map<Sequence, double>::iterator it;

			for (it = p.pop.begin(); it != p.pop.end(); it++) {
				res.pop[it->first] -= it->second;

				if (res.pop[it->first] == 0)
					res.pop.erase(it->first);
			}

			return res;
		}

		// Multiplication by scalar:
		GenPopulation operator* (double p) {

			GenPopulation res;

			if (p > 0) {
				std::map<Sequence, double>::iterator it;

				for (it = res.pop.begin(); it != res.pop.end(); it++)
					it->second *= p;
			}

			return res;
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
};


// Class for genetically homogeneous populations
class NongenPopulation : public Population {
	public:

		double n;

		// Constructors:
		NongenPopulation (int p_n) {
			n = p_n;
			isGen = false;
		}
		NongenPopulation () {
			isGen = false;
		}
		NongenPopulation (NongenPopulation & p) {
			isGen = false;
			n = p.n;
		}

		// Assignment operator:
		NongenPopulation operator= (NongenPopulation p) {
			n = p.n;

			return *this;
		}

		// Arithmetic operators:
		NongenPopulation operator+ (NongenPopulation p) {
			NongenPopulation res = *this;

			res.n += p.n;

			return res;
		}

		NongenPopulation operator- (NongenPopulation p) {
			NongenPopulation res = *this;

			res.n -= p.n;

			return res;
		}

		// Multiplication by scalar:
		NongenPopulation operator* (double p) {

			NongenPopulation res = *this;

			res.n *= p;

			return res;
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
class StateVec : public std::vector<Population> {
	public:

		// Constructors:
		StateVec(int nSpecies)
		{
			resize(nSpecies);
		}
		StateVec() {};

		// State vector algebra:
		StateVec operator+ (StateVec p) {

			StateVec res = p;
			for (int i=0; i<size(); i++)
				res[i] = res[i] + operator[](i);

			return res;
		}

		StateVec operator- (StateVec p) {

			StateVec res = *this;
			for (int i=0; i<size(); i++)
				res[i] = res[i] - p[i];

			return res;
		}

		// Multiplication by scalar:
		StateVec operator* (double p) {

			StateVec res = *this;
			for (int i=0; i<size(); i++)
				res[i] = res[i] * p;

			return res;
		}

		// Check for negative populations:
		bool isnegative() {

			for (int i=0; i<size(); i++)
				if (operator[](i).isnegative())
					return true;

			return false;
		}

};


// Class describing individual reactions:
class Reaction {
	public:

		double rate;
		std::vector<int> in, out;
		std::vector<bool> mutate;

		// Critical reaction number:
		int Nc;

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
		Reaction(std::vector<int> p_in, std::vector<int> p_out,
				std::vector<bool> p_mutate,
				double p_rate, bool p_isGen, int p_Nc)
		{
			rate = p_rate;
			in = p_in;
			out = p_out;
			mutate = p_mutate;

			isGen = p_isGen;
			isMutation = false;
			for (int i=0; i<mutate.size(); i++)
				if (mutate[i])
					isMutation = true;

			Nc = p_Nc;
		}
		Reaction() {};

		// Calculate reaction propensities for given state,
		// determine critical reactions and time and sequence
		// corresponding to next critical reaction:
		void calcPropensity(StateVec x, unsigned short *buf)
		{
			// Reset critical state:
			criticalDelta = -1.0;

			if (isGen) {

				// Empty propensity and critical sequence containers:
				propensities.clear();
				criticalSeq.clear();

				// Determine reactant population with smallest diversity:
				int mindiversity = 0;
				int imin = 0;
				for (int i=0; i<x.size(); i++) {
					if (in[i]>0 && x[i].isGen) {
						if (imin>=0 || x[i].pop.size()<mindiversity) {
							mindiversity = x[i].pop.size();
							imin = i;
						}
					}
				}

				// Determine all reaction propensities:
				std::map<Sequence, double>::iterator it;
				for (it = x[imin].begin(); it != x[imin].end(); it++) {

					double a = 1.0;
					bool crit = false;
					for (int i=0; i<x.size(); i++) {
						if (x[i].isGen) {

							// Check for criticality
							if (!mutate[i]) {
							   	if (x[i][it->first] < Nc*(in[i]-out[i]))
									crit = true;
							} else {
							   	if (x[i][it->first] < Nc*in[i])
									crit = true;
							}

							// Calculate propensity contribution
							for (int m=0; m<in[i]; m++)
								a *= x[i][it->first] - m;

						} else {

							// Check for criticality
							if (x[i] < Nc*(in[i] - out[i]))
									crit = true;

							// Calculate propensity contribution
							for (int m=0; m<in[i]; m++)
								a *= x[i].popSize() - m;
						}
					}

					if (a>0) {
						if (crit) {

							// Critical reaction: choose reaction time

							if (isMutation) {

								for (int i=0; i<it->first.size()*(it->first.nChar-1); i++) {
									
									double delta = -log(erand48(buf))/a;
									if (criticalDelta < 0 || delta < criticalDelta) {
										criticalDelta = delta;
										criticalSeq = it->first;
									}
								}

							} else {

								double delta = -log(erand48(buf))/a;
								if (criticalDelta < 0 || delta < criticalDelta) {
									criticalDelta = delta;
									criticalSeq = it->first;
								}
							}

						} else {

							// Non-critical reaction: record propensity
							propensities[it->first] = a;
						}
					}
				}

			} else {

				double a = 1.0;
				bool crit = false;
				for (int i=0; i<x.size(); i++) {

					if (in[i] == 0)
						continue;

					// Check for criticality
					if (x[i] < Nc*(in[i] - out[i]))
						crit = true;

					// Calculate propensity contribution
					for (int m=0; m<in[i]; m++) {
						a *= x[i].popSize() - m;
					}
				}

				propensity = 0.0;

				if (a>0) {
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
			return criticalDelta < 0.0;
		}

		// Perform tau-leaping integration step:
		StateVec tauLeap(StateVec x, double dt, unsigned short *buf)
		{
			if (isGen) {

				std::map<Sequence, double>::iterator it;
				for (it = propensities.begin(); it != propensities.end(); it++) {

					Sequence thisSeq = it->first;
					double thisProp = it->second;

					if (isMutation) {
						std::set<Sequence> mutSequences = thisSeq.getNeighbours();

						std::set<Sequence>::iterator itm;
						for (itm = mutSequences.begin(); itm != mutSequences.end(); itm++) {
							Sequence mutant = itm->first;

							double nreacts = poissonian(thisProp*dt, buf);

							for (int i=0; i<x.size(); i++) {
								if (x[i].isGen) {
									if (!mutate[i]) {
										if (out[i]-in[i] != 0)
											x[i].pop[thisSeq] += nreacts*(out[i]-in[i]);
									} else {
										x[i].pop[thisSeq] -= nreacts*in[i];
										x[i].pop[mutantSeq] += nreacts*out[i];
									}
								} else {
									if (out[i]-in[i] != 0)
										x[i].n += nreacts*(out[i]-in[i]);
								}
							}
						}

					} else {

						double nreacts = poissonian(thisProp*dt, buf);

						for (int i=0; i<x.size(); i++) {
							if (out[i] - in[i] != 0) {
								if (x[i].isGen)
									x[i].pop[thisSeq] += nreacts*(out[i]-in[i]);
								else
									x[i].n += nreacts*(out[i]-in[i]);
							}
						}

					}

				}

			} else {

				double nreacts = poissonian(propensity*dt, buf);

				for (int i=0; i<x.size(); i++) {
					if (in[i]-out[i] != 0) {
						x[i].n += nreacts*(out[i]-in[i]);
					}
				}

			}
		}

		// Implement critical reaction on given state:
		StateVec implementCritical(StateVec x, std::vector<double> )
		{

			if (isGen) {

				if (isMutation) {

					std::set<Sequence> mutSequences = criticalSeq.getNeighbours();
					int nNeighbours = criticalSeq.size()*(criticalSeq.nChar-1);
					Sequence mutantSeq = criticalSeq.getNaighbours()[drand48(buf)%nNeighbours];

					for (int i=0; i<x.size(); i++) {
						if (x[i].isGen) {
							if (mutate[i]) {
								if (in[i] > 0)
									x[i].pop[criticalSeq] -= in[i];
								if (out[i] > 0)
									x[i].pop[mutantSeq] += out[i];
							} else {
								if (out[i]-in[i] != 0)
									x[i].pop[criticalSeq] += out[i]-in[i];
							}
						} else {
							x[i].n += out[i]-in[i];
						}
					}

				} else {

					for (int i=0; i<x.size(); i++) {
						if (x[i].isGen) {
							if (out[i]-in[i] != 0)
								x[i].pop[criticalSeq] += out[i]-in[i];
						} else
							x[i].n += out[i]-in[i];
					}

				}

			} else {

				for (int i=0; i<x.size(); i++)
					x[i].n += out[i]-in[i];

			}

			return x;
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
