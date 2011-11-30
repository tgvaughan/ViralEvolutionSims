/*
 * tauleap_PST_OptionParser.cc
 *
 *  Created on: 04/11/2011
 *      Author: Tim Vaughan
 */

#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>
#include "tauleap_PST_OptionParser.h"

/**
 * Obtains simulation parameters from command line and/or config file.
 *
 * argc: argument count
 * argv: argument value array
 *
 * returns: variable_map populated with parameters to use
 */
boost::program_options::variables_map OptionParser::parse(int argc, char **argv)
{
	namespace po = boost::program_options;

	po::options_description desc_model("Model parameters");
	desc_model.add_options()
			("model.lambda", po::value<double>()->default_value(2.5e8), "target cell production rate")
			("model.beta", po::value<double>()->default_value(5e-13), "cellular infection rate")
			("model.k", po::value<double>()->default_value(1e3), "virion production rate")
			("model.d", po::value<double>()->default_value(1e-3), "uninfected cell death rate")
			("model.a", po::value<double>()->default_value(1.0), "infected cell death rate")
			("model.u", po::value<double>()->default_value(3.0), "virion clearance rate")

			("model.lat_p", po::value<double>()->default_value(0.0), "probability of producing latently infected cell")
			("model.lat_d", po::value<double>()->default_value(0.0), "latently infected cellular death rate")
			("model.lat_a", po::value<double>()->default_value(0.0), "activation rate of latently infected cell")

			("model.sequenceL", po::value<double>()->default_value(105), "virus genome length")
			("model.mu_RT", po::value<double>()->default_value(105*2e-5), "mutation probability per replication")
			("model.mu_RNAP", po::value<double>()->default_value(0), "mutation probability per replication (RNAP)")
			;

	po::options_description desc_sim("Simulation parameters");
	desc_sim.add_options()
			("simulation.T", po::value<double>()->default_value(100), "length of simulation (days)")
			("simulation.X0", po::value<double>()->default_value(2.5e11), "initial number of uninfected cells")
			("simulation.Y0", po::value<double>()->default_value(0), "initial number of infected cells")
			("simulation.YL0", po::value<double>()->default_value(0), "initial number of latently infected cells")
			("simulation.V0", po::value<double>()->default_value(100), "initial number of virions")
			;

	po::options_description desc_alg("Algorithm parameters");
	desc_alg.add_options()
			("algorithm.Nt", po::value<int>()->default_value(10001), "number tau-leaping time steps")
			("algorithm.Ncrit", po::value<int>()->default_value(500), "critical reaction parameter")
			("algorithm.Npaths", po::value<int>()->default_value(1), "number of trajectories to generate")
			("algorithm.Nsamples", po::value<int>()->default_value(1001), "number of samples to record")
			;

	po::options_description desc_clhidden;
	desc_clhidden.add_options()
			("outfile", po::value<std::string>(), "location of output file")
			;

	po::options_description desc_clvisible;
	desc_clvisible.add_options()
			("help", "display usage information and exit")
			("config", po::value<std::string>(), "location of configuration file")
			;
	desc_clvisible.add(desc_model).add(desc_sim).add(desc_alg);

	po::positional_options_description posdesc;
	posdesc.add("outfile", 1);

	po::options_description desc_commandline;
	desc_commandline.add(desc_clvisible).add(desc_clhidden);

	// Process command line options:
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc_commandline).positional(posdesc).run(), vm);
	po::notify(vm);

	// Print usage info if outfile was not provided:
	if (vm.count("help") || vm.count("outfile")==0) {
		std::cout << "Usage: " << argv[0] << " [options] outfile[.h5]" << std::endl
				<< desc_clvisible << std::endl;
		exit(0);
	}

	// Process configuration file if filename provided:
	if (vm.count("config")) {
		std::ifstream cfile;
		cfile.open(vm["config"].as<std::string>().c_str());
		if (!cfile.is_open()) {
			std::cout << "Error opening config file '" << vm["config"].as<std::string>()
					<< "'." << std::endl;
			exit(1);
		}

		po::options_description desc_config;
		desc_config.add(desc_model).add(desc_sim).add(desc_alg);
		po::store(po::parse_config_file(cfile, desc_config, false), vm);
		po::notify(vm);
	}

	return vm;
}
