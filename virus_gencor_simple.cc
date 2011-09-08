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
            Y.resize(L,0.0);
            V.resize(L,0.0);
        }

        // Copy constructor:
        StateVec (const StateVec & sv) {
            L = sv.L;
            Y = sv.Y;
            V = sv.V;
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
