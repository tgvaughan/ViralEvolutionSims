/*
 * tauleap_PST_OptionParser.h
 *
 *  Created on: 04/11/2011
 *      Author: Tim Vaughan
 */

#ifndef TAULEAP_PST_OPTIONPARSER_H_
#define TAULEAP_PST_OPTIONPARSER_H_

class OptionParser {
public:
	static boost::program_options::variables_map parse(int argc, char **argv);
};

#endif /* TAULEAP_PST_OPTIONPARSER_H_ */
