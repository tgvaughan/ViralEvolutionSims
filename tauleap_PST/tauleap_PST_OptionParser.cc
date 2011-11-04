/*
 * tauleap_PST_OptionParser.cc
 *
 *  Created on: 04/11/2011
 *      Author: Tim Vaughan
 */

#include <string>

#include <boost/program_options.hpp>
#include "tauleap_PST_OptionParser.h"

boost::program_options::variables_map OptionParser::parse(int argc, char **argv)
{
	namespace po = boost::program_options;

	po::options_description desc_io("I/O options");
	desc_io.add_options()
			("outfile", po::value<std::string>(), "location of HDF5 file to write results to")
			("config", po::value<std::string>(), "location of configuration file")
			;

	po::options_description desc_model("Model parameters");
	desc_model.add_options()
			("lambda", po::value<double>(), "target cell production rate")
			("beta", po::value<double>(), "cellular infection rate")
			("k", po::value<double>(), "virion production rate")
			("d", po::value<double>(), "uninfected cell death rate")
			("a", po::value<double>(), "infected cell death rate")
			("u", po::value<double>(), "virion clearance rate")
			;

	po::options_description desc_sim("Simulation parameters");
	desc_sim.add_options()
			("T", po::value<double>(), "length of simulation (days)")
			("X0", po::value<double>(), "initial number of uninfected cells")
			("Y0", po::value<double>(), "initial number of infected cells")
			("V0", po::value<double>(), "initial number of virions")
			;

	po::options_description desc_alg("Algorithm parameters");
	desc_alg.add_options()
			("Nt", po::value<int>(), "number tau-leaping time steps")
			("alpha", po::value<double>(), "critical reaction parameter")
			("Npaths", po::value<int>(), "number of trajectories to generate")
			("Nsamples", po::value<int>(), "number of samples to record")
			;

	po::positional_options_description posdesc;
	posdesc.add("config", 1).add("outfile",1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc_io).positional(posdesc).run(), vm);

	return vm;
}

