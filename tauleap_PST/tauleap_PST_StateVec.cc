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
	maxHD = length;

	X = 0.0;
	Y.resize(maxHD+1,0.0);
	YL.resize(maxHD+1,0.0);
	V.resize(maxHD+1,0.0);
}

StateVec::StateVec(int length, int maxHammingDistance) {
	L = length;
	if (maxHammingDistance>0)
		maxHD = maxHammingDistance;
	else
		maxHD = length;

	X = 0.0;
	Y.resize(maxHD+1);
	YL.resize(maxHD+1);
	V.resize(maxHD+1);
}

StateVec::StateVec (const StateVec & src) {
	L = src.L;
	maxHD = src.maxHD;

	X = src.X;
	Y = src.Y;
	YL = src.YL;
	V = src.V;
}

StateVec StateVec::operator= (const StateVec & src) {

	L = src.L;
	maxHD = src.maxHD;

	X = src.X;
	Y = src.Y;
	YL = src.YL;
	V = src.V;

	return *this;
}
