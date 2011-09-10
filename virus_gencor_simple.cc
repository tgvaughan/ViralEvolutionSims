// Within-host viral evolution simulation.  Makes use of a projected
// sequence space of dramatically reduced dimension.

#include <iostream>

#include <cmath>

#include "poissonian.h"

class StateVec {
    public:

        // Sequence length:
        int L;

        // Target cell population:
        double X;

        // Infected cell and virion populations:
        std::vector<double> Y, V;

        // Constructor:
        StateVec (int length) {
            L = length;

            X = 0.0;
            Y.resize(L+1,0.0);
            V.resize(L+1,0.0);
        }

        // Copy constructor:
        StateVec (const StateVec & sv) {
            L = sv.L;
            Y = sv.Y;
            V = sv.V;
        }

};

// Number of neighbouring sequences of a sequence on
// h1 which lie on h2:
double gcond(int h2, int h1)
{
	switch (h2) {
		case h1:
		case h1+1:
		case h1-1:
		default:
			return 0.0;
	}

}

// Reaction class
class Reaction {
	public:

		int inX, inY, inV;
		int outX, outY, outV;
		bool mutY, mutV;

		// Constructor:
		Reaction(int p_inX, int p_inY, int p_inV,
				int p_outX, int p_outY, int p_outV,
				bool p_mutY, bool p_mutV) {

			inX = p_inX;
			inY = p_inY;
			inV = p_inV;
			outX = p_outX;
			outY = p_outY;
			outV = p_outV;

			mutY = p_mutY;
			mutV = p_mutV;
		}

		// Get reaction propensity:
		void getPropensity(StateVec & sv) {
		}

};

int main (int argc, char **argv)
{
    using namespace std;

    // Parse command line parameters:
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " outfile" << endl;
        exit(0);
    }
    char *ofname = argv[1];

    // Simulation parameters:

    // Model parameters:

    // Initialise RNG:


}
