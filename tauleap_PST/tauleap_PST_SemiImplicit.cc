/*
 * tauleap_PST_SemiImplicit.cc
 *
 *  Created on: 07/11/2011
 *      Author: Tim Vaughan
 */

#include "poissonian.h"
#include "tauleap_PST_StateVec.h"
#include "tauleap_PST_Reaction.h"
#include "tauleap_PST_SemiImplicit.h"

bool SemiImplicit::step(StateVec & sv, Reaction *reactions, int Nreactions, double tau, unsigned short  *buf)
{
	bool negativePop = false;

	return negativePop;
}

