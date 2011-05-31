// virus_gencor.h: Class definitions


void print_seq(const std::vector<int> & seq)
{
	for (unsigned int i=0; i<seq.size(); i++)
		std::cout << seq[i];
}


class StateVec {
	public:

		int X;
		std::map<std::vector<int>,int> Y;
		std::map<std::vector<int>,int> V;

		int NY, NV;
		bool totals_valid;

		// Constructor:
		StateVec () {};

		// Destructor:
		~StateVec () {};

		// Copy Constructor:
		StateVec (const StateVec & sv)
		{
			X = sv.X;
			Y = sv.Y;
			V = sv.V;

			totals_valid = false;
		};

		// Assignment operator:
		StateVec operator= (const StateVec & sv)
		{
			X = sv.X;
			Y = sv.Y;
			V = sv.V;

			totals_valid = false;

			return *this;
		}

		// Calculate total population sizes:
		void calc_totals()
		{
			std::map<std::vector<int>,int>::iterator it;

			NY = 0, NV = 0;
			for (it = Y.begin(); it != Y.end(); ++it)
				NY += it->second;
			for (it = V.begin(); it != V.end(); ++it)
				NV += it->second;

			totals_valid = true;
		}

		// Return process rate multipliers:
		double get_inf_rate()
		{
			if (!totals_valid)
				calc_totals();

			return (double)X*(double)NV;
		}
		double get_vprod_rate()
		{
			if (!totals_valid)
				calc_totals();

			return (double)NY;
		}
		double get_xdeath_rate()
		{
			return (double)X;
		}
		double get_ydeath_rate()
		{
			if (!totals_valid)
				calc_totals();

			return (double)NY;
		}
		double get_vdeath_rate()
		{
			if (!totals_valid)
				calc_totals();

			return (double)NV;
		}

		// Select random infected cell:
		std::vector<int> choose_Y(unsigned short *seedbuff)
		{
			if (!totals_valid)
				calc_totals();

			std::map<std::vector<int>,int>::iterator it;

			int u = nrand48(seedbuff)%NY;
			int n_cumul = 0;
			for (it = Y.begin(); it != Y.end(); ++it) {
				n_cumul += it->second;
				if (u<n_cumul)
					break;
			}

			return it->first;
		}

		// Select random virion:
		std::vector<int> choose_V(unsigned short *seedbuff)
		{
			if (!totals_valid)
				calc_totals();

			std::map<std::vector<int>,int>::iterator it;

			int u = nrand48(seedbuff)%NV;
			int n_cumul = 0;
			for (it = V.begin(); it != V.end(); ++it) {
				n_cumul += it->second;
				if (u<n_cumul)
					break;
			}

			return it->first;
		}

		// Increment or decrement Y and V populations:
		void inc_Y(std::vector<int> & seq)
		{
			++Y[seq];

			totals_valid = false;
		}
		void inc_V(std::vector<int> & seq)
		{
			++V[seq];

			totals_valid = false;
		}
		void dec_Y(std::vector<int> seq)
		{
			if (--Y[seq] <= 0)
				Y.erase(seq);

			totals_valid = false;
		}
		void dec_V(std::vector<int> seq)
		{
			if (--V[seq] <= 0)
				V.erase(seq);

			totals_valid = false;
		}

		// Arithmetic operations:
		StateVec operator+ (const StateVec & arg)
		{
			StateVec result;

			result.X = X + arg.X;

			std::map<std::vector<int>,int>::iterator it;

			result.Y = arg.Y;
			for (it = Y.begin(); it != Y.end(); ++it)
				result.Y[it->first] += it->second;

			result.V = arg.V;
			for (it = V.begin(); it != V.end(); ++it)
				result.V[it->first] += it->second;

			return result;
		}

		// Debugging methods:
		void dump_vseqs()
		{
			for (std::map<std::vector<int>,int>::iterator it = V.begin();
					it != V.end(); ++it) {

				print_seq(it->first);
				std::cout << ": " << it->second << std::endl;
				
			}
		}
};
