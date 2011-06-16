
// Class for system state vectors
class StateVec : public std::vector<double> {
	public:

		// Constructors:
		StateVec(int nSpecies)
		{
			resize(nSpecies);
		}
		StateVec() {};

		// State vector algebra:
		StateVec operator+ (StateVec p)
		{

			StateVec res = p;
			for (int i=0; i<size(); i++)
				res[i] += operator[](i);

			return res;
		}

		StateVec operator- (StateVec p)
		{

			StateVec res = *this;
			for (int i=0; i<size(); i++)
				res[i] -= p[i];

			return res;
		}

		// Multiplication by scalar:
		StateVec operator* (double p)
		{

			StateVec res = *this;
			for (int i=0; i<size(); i++)
				res[i] *= p;

			return res;
		}

		// Check for negative populations:
		bool isnegative()
		{

			for (int i=0; i<size(); i++)
				if (operator[](i) < 0)
					return true;

			return false;
		}

};


// Class describing individual reactions:
class Reaction {
	public:
		double rate;
		StateVec in, out, delta;

		// Constructors:
		Reaction(StateVec p_in, StateVec p_out, double p_rate)
		{
			rate = p_rate;
			in = p_in;
			out = p_out;
			delta = out - in;
		}
		Reaction() {};

		// Calculate reaction propensity for given state:
		double propensity(StateVec x)
		{
			double a = rate;

			for (int i=0; i<x.size(); i++) {
				for (int m=0; m<in[i]; m++)
					a *= x[i]-(double)m;
			}

			if (a<0)
				return 0;

			return a;
		}

		// Implement reaction on given state:
		StateVec implement(StateVec x)
		{
			return x + delta;
		}

		// Check whether state is within Nc reactions of bottoming out:
		bool iscritical (StateVec x, int Nc)
		{
			return (x + delta*(double)Nc).isnegative();
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

