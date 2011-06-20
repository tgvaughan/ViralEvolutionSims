// Definition of various classes used by virus_gencor_tauleap.cc.

// Abstract base class for populations
class Population {
	public:
		bool isGen;

		virtual bool isnegative() =0;
		virtual int operator[] (std::vector<double>) =0;

		virtual Population operator= (Population p) =0;
		virtual Population operator+ (Population p) =0;
		virtual Population operator- (Population p) =0;
		virtual Population operator* (double p) =0;
};

// Class for genetically diverse populations
class GenPopulation : public Population {
	public:
		std::map<vector<double>, double> pop;
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
		double operator[] (std::vector<double> p) {
			std::map<vector<double>, double>::iterator it = pop.find(p);
			if (it == pop.end())
				return 0;

			return it->second;
		}

		// Return iterators:
		map<vector<double>, double>::iterator begin() {
			return pop.begin();
		}
		map<vector<double>, double>::iterator end() {
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

			std::map<vector<double>, double>::iterator it;

			for (it = p.pop.begin(); it != p.pop.end(); it++) {
				res.pop[it->first] += it->second;

				if (res[it->first] == 0)
					res.pop.erase(it->first);
			}

			return res;
		}

		GenPopulation operator- (GenPopulation p) {

			GenPopulation res = *this;

			std::map<vector<double>, double>::iterator it;

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
				std::map<vector<double>. double>::iterator it;

				for (it = res.pop.begin(); it != res.pop.end(); it++)
					it->second *= p;
			}

			return res;
		}

		// Negativity check:
		bool isnegative () {

			std::map<vector<double>, double>::iterator it;

			for (it = pop.begin(); it != pop.end(); it++)
				if (it->second < 0)
					return true;

			return false;
		}

		// Get total size:
		double size() {

			std::map<vector<double>, double>::iterator it;

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
		double size() {
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

		// Constructors:
		Reaction(std::vector<int> p_in, std::vector<int> p_out,
				std::vector<bool> p_mutate,
				double p_rate)
		{
			rate = p_rate;
			in = p_in;
			out = p_out;
			mutate = p_mutate;
		}
		Reaction() {};

		// Calculate reaction propensity for given state:
		double propensity(StateVec x)
		{
			double a = rate;

			for (int i = 0; i<in.size(); i++) {
				if (!x[i].isgen) {
					for (int m=0; m<in[i]; m++)
						a *= x[i].size() - m;
				} else {

					std::map<std::vector<double>, double>::iterator it;
					for (it = x[i].begin(); it != x[i].end(); it++) {

					}

				}
			}

			return a;
		}

		// Implement reaction on given state:
		StateVec implement(StateVec x, unsigned short *buf)
		{
			// TODO
		}

		// Check whether state is within Nc reactions of bottoming out:
		bool iscritical (StateVec x, int Nc)
		{
			// TODO
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

