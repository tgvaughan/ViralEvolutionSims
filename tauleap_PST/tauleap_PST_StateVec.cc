/*
 * tauleap_PST_StateVec.cc
 *
 *  Created on: 06/11/2011
 *      Author: Tim Vaughan
 */

#include <vector>

#include "tauleap_PST_StateVec.h"

StateVec::StateVec (int length) {
	L = length;

	X = 0.0;
	Y.resize(L+1,0.0);
	V.resize(L+1,0.0);
}

StateVec::StateVec (const StateVec & src) {
	L = src.L;

	X = src.X;
	Y = src.Y;
	V = src.V;
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;

	X = src.X;
	Y = src.Y;
	V = src.V;

	return *this;
}
