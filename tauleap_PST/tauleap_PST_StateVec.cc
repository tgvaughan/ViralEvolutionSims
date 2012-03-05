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
	vecL = length+1;

	X = 0.0;
	Y.resize(vecL,0.0);
	YL.resize(vecL,0.0);
	V.resize(vecL,0.0);
}

StateVec::StateVec(int length, int vecLength) {
	L = length;
	vecL = vecLength;

	X = 0.0;
	Y.resize(vecL);
	YL.resize(vecL);
	V.resize(vecL);
}

StateVec::StateVec (const StateVec & src) {
	L = src.L;
	vecL = src.vecL;

	X = src.X;
	Y = src.Y;
	YL = src.YL;
	V = src.V;
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;
	vecL = src.vecL;

	X = src.X;
	Y = src.Y;
	YL = src.YL;
	V = src.V;

	return *this;
}
